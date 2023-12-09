/* This file is part of NBI, a library for Nystrom Boundary Integral solvers
 *
 * Copyright (C) 2021, 2023 Michael Carley
 *
 * NBI is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version. NBI is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with NBI.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <glib.h>

#include <nbi.h>

#ifdef HAVE_SQT
#include <sqt.h>
#endif /*HAVE_SQT*/

#ifdef HAVE_BLASWRAP
#include <blaswrap.h>
#endif /*HAVE_BLASWRAP*/

#include "nbi-private.h"

NBI_REAL **_Ku = NULL, **_Ki = NULL ;
gint _Nk[16] ;

nbi_surface_t *NBI_FUNCTION_NAME(nbi_surface_alloc)(gint nnmax, gint npmax)

{
  nbi_surface_t *s ;

  s = (nbi_surface_t *)g_malloc(sizeof(nbi_surface_t)) ;

  nbi_surface_node_number(s) = 0 ;
  nbi_surface_node_number_max(s) = nnmax ;
  nbi_surface_patch_number(s) = 0 ;
  nbi_surface_patch_number_max(s) = npmax ;

  s->xc  = g_malloc0(nnmax*NBI_SURFACE_NODE_LENGTH*sizeof(NBI_REAL)) ;
  s->pcr = g_malloc0(npmax*NBI_SURFACE_PATCH_DATA_LENGTH*sizeof(NBI_REAL)) ;
  s->ip  = (gint *)g_malloc0(NBI_SURFACE_PATCH_LENGTH*npmax*sizeof(gint)) ;

  s->fpsize = sizeof(NBI_REAL) ;
  
  return s ;
}


gint NBI_FUNCTION_NAME(nbi_surface_write)(nbi_surface_t *s, FILE *f)

{
  gint i, j ;
  NBI_REAL *e ;
  gchar header[NBI_HEADER_LENGTH], buf[40] ;

  nbi_header_init(header, "NBI", "1.0", "GEO", "A") ;

  sprintf(buf, "%d %d",
	  nbi_surface_node_number(s),
	  nbi_surface_patch_number(s)) ;

  nbi_header_insert_string(header, NBI_HEADER_DATA, 40, buf) ;

  nbi_header_write(f, header) ;
  
  for ( i = 0 ; i < nbi_surface_patch_number(s) ; i ++ ) {
    fprintf(f, "%d %d\n",
	    nbi_surface_patch_node(s, i),
	    nbi_surface_patch_node_number(s, i)) ;
  }

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    e = (NBI_REAL *)nbi_surface_node(s,i) ;
    for ( j = 0 ; j < NBI_SURFACE_NODE_LENGTH ; j ++ )
      fprintf(f, " %1.16e", e[j]) ;
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

nbi_surface_t *NBI_FUNCTION_NAME(nbi_surface_read)(FILE *f)

{
  nbi_surface_t *s ;
  gint nn, np, i, j ;
  NBI_REAL *e ;
  gchar header[NBI_HEADER_LENGTH] ;

  nbi_header_read(f, header) ;

  sscanf(&(header[NBI_HEADER_DATA]), "%d %d", &nn, &np) ;

  s = nbi_surface_alloc(nn, np) ;
  nbi_surface_node_number(s) = nn ;
  nbi_surface_patch_number(s) = np ;

  for ( i = 0 ; i < np ; i ++ ) {
    fscanf(f, "%d", &(nbi_surface_patch_node(s, i))) ;
    fscanf(f, "%d", &(nbi_surface_patch_node_number(s, i))) ;
  }

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    e = (NBI_REAL *)nbi_surface_node(s,i) ;    
    for ( j = 0 ; j < NBI_SURFACE_NODE_LENGTH ; j ++ )
      fscanf(f, "%lg", &(e[j])) ;
  }
  
  return s ;
}

gint NBI_FUNCTION_NAME(nbi_surface_patch_centroid)(NBI_REAL *x, gint xstr,
						   NBI_REAL *w, gint wstr,
						   gint nx,
						   NBI_REAL *c)

{
  gint i ;
  NBI_REAL A ;
  
  c[0] = c[1] = c[2] = A = 0.0 ;

  for ( i = 0 ; i < nx ; i ++ ) {
    c[0] += x[i*xstr+0]*w[i*wstr] ;
    c[1] += x[i*xstr+1]*w[i*wstr] ;
    c[2] += x[i*xstr+2]*w[i*wstr] ;
    A    +=             w[i*wstr] ;
  }

  c[0] /= A ; c[1] /= A ; c[2] /= A ; 
  
  return 0 ;
}

NBI_REAL NBI_FUNCTION_NAME(nbi_surface_patch_radius)(NBI_REAL *x, gint xstr,
						     gint nx, NBI_REAL *c)

{
  NBI_REAL r, ri ;
  gint i ;

  r = 0.0 ;
  for ( i = 0 ; i < nx ; i ++ ) {
    ri = nbi_vector_distance(&(x[i*xstr]), c) ;    
    r = MAX(r, ri) ;
  }
  
  return r ;
}

gint NBI_FUNCTION_NAME(nbi_patch_neighbours)(NBI_REAL *c, NBI_REAL r,
					     NBI_REAL *x, gint xstr, gint nx,
					     gint n0, gint n1,
					     gint *nbrs, gint *nnbrs,
					     gint nnmax)

{
  NBI_REAL r2 = r*r, rx2 ;
  gint i ;
  
  for ( i = n0 ; (i < n1) && (*nnbrs < nnmax) ; i ++ ) {
    if ( (rx2 = nbi_vector_distance2(c, &(x[i*xstr]))) < r2 ) {
      nbrs[(*nnbrs)] = i ; (*nnbrs) ++ ;
    }
  }
  
  return 0 ;
}

static gboolean patch_boxes_encroach(NBI_REAL *c, NBI_REAL r,
				     NBI_REAL *ct, NBI_REAL rt)

{
  if ( fabs(c[0] - ct[0]) > r+rt ) return FALSE ;
  if ( fabs(c[1] - ct[1]) > r+rt ) return FALSE ;
  if ( fabs(c[2] - ct[2]) > r+rt ) return FALSE ;
  
  return TRUE ;
}

gint NBI_FUNCTION_NAME(nbi_surface_patch_neighbours)(nbi_surface_t *s,
						     gint p, NBI_REAL r,
						     gint *nbrs, gint *nnbrs,
						     gint nnmax)

{
  gint i, n0, n1 ;
  NBI_REAL *cp, rp, *ci, ri, *x ;

  cp = (NBI_REAL *)nbi_surface_patch_centre(s, p) ;
  rp = *((NBI_REAL *)nbi_surface_patch_sphere_radius(s, p)) ;
  
  x = (NBI_REAL *)nbi_surface_node(s, 0) ;
  for ( i = 0 ; i < p ; i ++ ) {
    ci = (NBI_REAL *)nbi_surface_patch_centre(s, i) ;
    ri = *((NBI_REAL *)nbi_surface_patch_sphere_radius(s, i)) ;
    if ( patch_boxes_encroach(cp, rp+r, ci, ri) ) {
      n0 = nbi_surface_patch_node(s, i) ;
      n1 = n0 + nbi_surface_patch_node_number(s, i) ;
      nbi_patch_neighbours(cp, r,
			   x, NBI_SURFACE_NODE_LENGTH,
			   nbi_surface_node_number(s),
			   n0, n1, 
			   nbrs, nnbrs, nnmax) ;
    }
  }

  for ( i = p+1 ; i < nbi_surface_patch_number(s) ; i ++ ) {
    ci = (NBI_REAL *)nbi_surface_patch_centre(s, i) ;
    ri = *((NBI_REAL *)nbi_surface_patch_sphere_radius(s, i)) ;
    if ( patch_boxes_encroach(cp, rp+r, ci, ri) ) {
      n0 = nbi_surface_patch_node(s, i) ;
      n1 = n0 + nbi_surface_patch_node_number(s, i) ;
      nbi_patch_neighbours(cp, r,
			   x, NBI_SURFACE_NODE_LENGTH,
			   nbi_surface_node_number(s),
			   n0, n1, 
			   nbrs, nnbrs, nnmax) ;
    }
  }
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_element_interp_matrix)(gint ns, NBI_REAL **K,
						  gint *Nk)

{
  gint nqi[] = {7, 25, 54, 85, 126, 175, 453} ;
  gint i, nq = 7, order ;
  NBI_REAL *st ;
  
  if ( _Ki == NULL ) _Ki = (NBI_REAL **)g_malloc0(nq*   sizeof(NBI_REAL *)) ;
  
  for ( i = 0 ; i < nq ; i ++ ) if ( nqi[i] == ns ) break ;

  if ( _Ki[i] == NULL ) {
    NBI_FUNCTION_NAME(sqt_quadrature_select)(ns, &st, &order) ;
    _Ki[i] = (NBI_REAL *)g_malloc0(ns*ns*sizeof(NBI_REAL)) ;
    _Nk[i] = NBI_FUNCTION_NAME(sqt_koornwinder_interp_matrix)(&(st[0]), 3,
							      &(st[1]), 3,
							      &(st[2]), 3,
							      ns, _Ki[i]) ;
  }

  *K = _Ki[i] ; *Nk = _Nk[i] ;
  
  return 0 ;
}

NBI_REAL *NBI_FUNCTION_NAME(nbi_patch_upsample_matrix)(gint ns, gint nu)

{
  gint nqi[] = {7, 25, 54, 85, 126, 175, 453} ;
  gint i, j, nq = 7, order, Nk ;
  NBI_REAL *st, work[453*3], *Ki ;
  
  if ( _Ku == NULL ) _Ku = (NBI_REAL **)g_malloc0(nq*nq*sizeof(NBI_REAL *)) ;
  
  for ( i = 0 ; i < nq ; i ++ ) if ( nqi[i] == ns ) break ;
  for ( j = 0 ; j < nq ; j ++ ) if ( nqi[j] == nu ) break ;

  NBI_FUNCTION_NAME(nbi_element_interp_matrix)(ns, &Ki, &Nk) ;
  
  if ( _Ku[i*nq+j] == NULL ) {
    /*generate the matrix*/
    NBI_FUNCTION_NAME(sqt_quadrature_select)(nu, &st, &order) ;
    _Ku[i*nq+j] = (NBI_REAL *)g_malloc0((ns+8)*nu*sizeof(NBI_REAL)) ;
    NBI_FUNCTION_NAME(sqt_interp_matrix)(Ki, ns, Nk,
					 &(st[0]), 3, &(st[1]), 3, nu,
					 _Ku[i*nq+j], work) ;
  }

  return _Ku[i*nq+j] ;
}

gint NBI_FUNCTION_NAME(nbi_surface_set_patch_data)(nbi_surface_t *s)

{
  gint i, ip ;

  for ( i = 0 ; i < nbi_surface_patch_number(s) ; i ++ ) {
    NBI_REAL *w, *n, *c, *r ;
    gint nnodes ;
    
    ip = nbi_surface_patch_node(s, i) ;
    nnodes = nbi_surface_patch_node_number(s, i),
    n = (NBI_REAL *) nbi_surface_node(s,ip) ;
    w = (NBI_REAL *) nbi_surface_node_weight(s, ip) ;
    c = (NBI_REAL *) nbi_surface_patch_centre(s, i) ;
    r = (NBI_REAL *)nbi_surface_patch_sphere_radius(s, i) ;

    NBI_FUNCTION_NAME(nbi_surface_patch_centroid)(n, NBI_SURFACE_NODE_LENGTH,
						  w, NBI_SURFACE_NODE_LENGTH,
						  nnodes,
						  c) ;
    *r = NBI_FUNCTION_NAME(nbi_surface_patch_radius)(n,
						     NBI_SURFACE_NODE_LENGTH, 
						     nnodes, c) ;
  }    
  
  return 0 ;
}

nbi_matrix_t *NBI_FUNCTION_NAME(nbi_matrix_new)(nbi_surface_t *s)

{
  nbi_matrix_t *m ;

  g_assert(s->fpsize == sizeof(NBI_REAL)) ;

  m = (nbi_matrix_t *)g_malloc(sizeof(nbi_matrix_t)) ;

  m->diag = 0.0 ;
  
  m->problem = NBI_PROBLEM_UNDEFINED ;
  m->potential = NBI_POTENTIAL_UNDEFINED ;
  m->s = s ;
  m->fpsize = sizeof(NBI_REAL) ;

  m->ustr = m->pstr = m->nstr = 0 ;
  m->idxu = m->idx = NULL ;
  m->xu = m->p = m->pn = NULL ;

  m->tree    = NULL ;
  m->targets = NULL ;
  m->shifts  = NULL ;

  m->test = NULL ;
  
  return m ;
}

static void read_upsampled_patches(FILE *f,
				   gint np, gint ustr,
				   nbi_matrix_t *m)

{
  gint i, j ;
  NBI_REAL *xu, *bc ;

  m->ustr = ustr ;
  m->idxu = (gint *)g_malloc((np+1)*sizeof(gint)) ;
  for ( i = 0 ; i < np + 1 ; i ++ ) {
    fscanf(f, "%d", &j) ;
    g_assert(j == i) ;
    fscanf(f, "%d", &(m->idxu[i])) ;
  }
  
  m->xu = g_malloc(m->idxu[np]*(m->ustr)*sizeof(NBI_REAL)) ;
  xu = (NBI_REAL *)(m->xu) ;
  for ( i = 0 ; i < (m->idxu)[np]*(m->ustr) ; i ++ )
    fscanf(f, "%lg", &(xu[i])) ;

  switch ( m->problem ) {
  default: g_assert_not_reached() ; break ;
  case NBI_PROBLEM_LAPLACE  : m->pstr = m->nstr = 2 ; break ;
  case NBI_PROBLEM_HELMHOLTZ: m->pstr = m->nstr = 4 ; break ;
  }
  m->bc = g_malloc(m->idxu[np]*(m->pstr)*sizeof(NBI_REAL)) ;
  bc = (NBI_REAL *)(m->bc) ;
  
  m->p  = (gchar *)(&(bc[0])) ;
  m->pn = (gchar *)(&(bc[(m->pstr)/2])) ;
  
  return ;
}

static gint read_correction_matrices(FILE *f,
				     gint np, gint nst,
				     nbi_matrix_t *m)

{
  gint i, j, lda ;
  NBI_REAL *Ast ;

  m->idxp = (gint *)g_malloc0((np+1)*sizeof(gint)) ;

  for ( i = 0 ; i < np+1 ; i ++ ) {
    fscanf(f, "%d", &j) ;
    g_assert(j == i) ;
    fscanf(f, "%d", &(m->idxp[i])) ;
  }
  
  m->idx = (gint *)g_malloc0((m->idxp[np])*sizeof(gint)) ;
  for ( i = 0 ; i < m->idxp[np] ; i ++ ) {
    fscanf(f, "%d", &(m->idx[i])) ;
  }

  switch ( m->problem ) {
  default:
    g_error("%s: unrecognized problem %d", __FUNCTION__, m->problem) ;
    break ;
  case NBI_PROBLEM_LAPLACE:   lda = 2*nst ; break ;
  case NBI_PROBLEM_HELMHOLTZ: lda = 4*nst ; break ;
  }

  Ast = (NBI_REAL *)g_malloc0(lda*(m->idxp[np])*sizeof(NBI_REAL)) ;
  m->Ast = (gchar *)Ast ;
  for ( i = 0 ; i < lda*(m->idxp[np]) ; i ++ ) {
    fscanf(f, "%lg", &(Ast[i])) ;
  }

  return np ;
}

gint NBI_FUNCTION_NAME(nbi_matrix_read)(FILE *f, nbi_matrix_t *m)

{
  gchar header[NBI_HEADER_LENGTH] ;
  gint np, str, nst ;
  
  nbi_header_read(f, header) ;

  sscanf(&(header[NBI_HEADER_DATA]), "%d %d %d", &np, &str, &nst) ;
  
  read_upsampled_patches(f, np, str, m) ;
  read_correction_matrices(f, np, nst, m) ;

  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_matrix_write)(FILE *f, nbi_matrix_t *m)

{
  gchar header[80], buf[40] ;
  NBI_REAL *Ast, *xu ;
  gint i, j, xstr, lda, nst, nntot ;
  
  nst = nbi_surface_patch_node_number(m->s,0) ;
  xstr = NBI_SURFACE_NODE_LENGTH ;

  switch ( m->problem ) {
  default:
    g_error("%s: unrecognized problem %d", __FUNCTION__, m->problem) ;
    break ;
  case NBI_PROBLEM_LAPLACE:   lda = 2*nst ; break ;
  case NBI_PROBLEM_HELMHOLTZ: lda = 4*nst ; break ;
  }
  
  nntot = m->idxp[nbi_surface_patch_number(m->s)] ;
  
  nbi_header_init(header, "NBI", "1.0", "MAT", "A") ;
  sprintf(buf, "%d %d %d %s",
	  nbi_surface_patch_number(m->s), xstr, nst,
	  nbi_problem_type_string(m->problem)) ;

  nbi_header_insert_string(header, NBI_HEADER_DATA, 40, buf) ;

  nbi_header_write(f, header) ;
  
  for ( i = 0 ; i < nbi_surface_patch_number(m->s)+1 ; i ++ ) {
    fprintf(f, "%d %d\n", i, m->idxu[i]) ;
  }

  xu = (NBI_REAL *)(m->xu) ;  
  for ( i = 0 ; i < m->idxu[nbi_surface_patch_number(m->s)] ; i ++ ) {
    for ( j = 0 ; j < xstr ; j ++ )
      fprintf(f, " %1.16e", xu[i*xstr+j]) ;
    fprintf(f, "\n") ;
  }
  
  for ( i = 0 ; i < nbi_surface_patch_number(m->s)+1 ; i ++ ) {
    fprintf(f, "%d %d\n", i, m->idxp[i]) ;
  }

  for ( i = 0 ; i < m->idxp[nbi_surface_patch_number(m->s)] ; i ++ ) {
    fprintf(f, "%d\n", m->idx[i]) ;
  }

  Ast = (NBI_REAL *)(m->Ast) ;
  for ( i = 0 ; i < nntot ; i ++ ) {
    for ( j = 0 ; j < lda ; j ++ )
      fprintf(f, " %1.16e", Ast[i*lda+j]) ;
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_data_write)(FILE *f, NBI_REAL *dat,
				       gint dstr, gint ne, gint nd)

{
  gchar header[NBI_HEADER_LENGTH], buf[40] ;
  gint i, j ;
  
  nbi_header_init(header, "NBI", "1.0", "DAT", "A") ;
  sprintf(buf, "%d %d", ne, nd) ;
  nbi_header_insert_string(header, NBI_HEADER_DATA, 40, buf) ;
  nbi_header_write(f, header) ;
  
  for ( i = 0 ; i < nd ; i ++ ) {
    for ( j = 0 ; j < ne ; j ++ ) {
      fprintf(f, " %1.16e", dat[i*dstr+j]) ;
    }
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

NBI_REAL *NBI_FUNCTION_NAME(nbi_data_read)(FILE *f, gint *nd, gint *ne)

{
  NBI_REAL *dat ;
  gchar header[NBI_HEADER_LENGTH] ;
  gint i ;
  
  nbi_header_read(f, header) ;

  sscanf(&(header[NBI_HEADER_DATA]), "%d %d", ne, nd) ;
  
  dat = (NBI_REAL *)g_malloc0((*ne)*(*nd)*sizeof(NBI_REAL)) ;

  for ( i = 0 ; i < (*nd)*(*ne) ; i ++ ) 
    fscanf(f, "%lg", &(dat[i])) ;
  
  return dat ;
}

gint NBI_FUNCTION_NAME(nbi_matrix_fmm_init)(nbi_matrix_t *m,
					    nbi_problem_t problem,
					    wbfmm_shift_operators_t *shifts,
					    guint *order_s, gint sstr,
					    guint *order_r, gint rstr,
					    gint depth,
					    NBI_REAL dtree,
					    gboolean shift_bw,
					    gboolean precompute_local,
					    NBI_REAL *work)

{
  gint fmmpstr, ustr, nsrc, i, nqfmm ;
  wbfmm_source_t source ;
  NBI_REAL xtree[3], xtmax[3], *xu, D, k ;
  guint order_max, field ;
  nbi_surface_t *s = m->s ;

  nqfmm = 1 ;
  order_max = 0 ;
  for ( i = 0 ; i <= depth ; i ++ ) {
    order_max = MAX(order_s[i*sstr], order_max) ;
    order_max = MAX(order_r[i*rstr], order_max) ;
  }

  xu = (NBI_REAL *)(m->xu) ;
  ustr = m->ustr ;
  nsrc = m->idxu[nbi_surface_patch_number(s)] ;
  
  source = WBFMM_SOURCE_MONOPOLE | WBFMM_SOURCE_DIPOLE ;
  field = WBFMM_FIELD_SCALAR ;
  xtree[0] = xtree[1] = xtree[2] = 0.0 ;
  wbfmm_points_origin_width(xu, ustr, nsrc, xtree, xtmax, &D, TRUE) ;
  wbfmm_points_origin_width((NBI_REAL *)nbi_surface_node(s,0),
			    NBI_SURFACE_NODE_LENGTH,
			    nbi_surface_node_number(s),
			    xtree, xtmax, &D, FALSE) ;

  xtree[0] -= dtree ; xtree[1] -= dtree ; xtree[2] -= dtree ;
  D += 2.0*dtree ;

  fmmpstr = ustr*sizeof(NBI_REAL) ;
  m->tree = wbfmm_tree_new(xtree, D, 2*nsrc) ;

  if ( shifts != NULL )
    m->shifts = shifts ;
  else
    m->shifts = wbfmm_shift_operators_new(order_max, shift_bw, work) ;

  wbfmm_tree_add_points(m->tree,
			(gpointer)xu, fmmpstr,
			(gpointer)(&(xu[3])), fmmpstr, nsrc, FALSE) ;
  
  for ( i = 0 ; i < depth ; i ++ ) wbfmm_tree_refine(m->tree) ;

  switch ( problem ) {
  default: g_assert_not_reached() ; break ;
  case NBI_PROBLEM_LAPLACE:
    wbfmm_tree_problem(m->tree) = WBFMM_PROBLEM_LAPLACE ;
    wbfmm_tree_source_size(m->tree) = nqfmm ;
    for ( i = 1 ; i <= depth ; i ++ ) {
      wbfmm_tree_laplace_coefficient_init(m->tree, i,
					  order_r[i*rstr], order_s[i*sstr]) ;
    }

    m->targets = NULL ;
    if ( !precompute_local ) return 0 ;

    m->targets = wbfmm_target_list_new(m->tree, nbi_surface_node_number(s)) ;
    wbfmm_target_list_coefficients_init(m->targets, field) ;
    wbfmm_target_list_add_points(m->targets,
				 (gpointer)nbi_surface_node(s,0),
				 NBI_SURFACE_NODE_LENGTH*sizeof(gdouble),
				 nbi_surface_node_number(s)) ;
  
    wbfmm_laplace_target_list_local_coefficients(m->targets, source, work) ;

    return 0 ;
  case NBI_PROBLEM_HELMHOLTZ:
    wbfmm_tree_problem(m->tree) = WBFMM_PROBLEM_HELMHOLTZ ;
    wbfmm_tree_source_size(m->tree) = nqfmm ;
    k = nbi_matrix_wavenumber(m) ;
    g_assert(k != 0.0) ;
    for ( i = 1 ; i <= depth ; i ++ ) {
      wbfmm_tree_coefficient_init(m->tree, i,
				  order_r[i*rstr], order_s[i*sstr]) ;
    }
    for ( i = 1 ; i <= depth ; i ++ ) {
      wbfmm_shift_operators_coaxial_SR_init(m->shifts, D, i,
					    order_s[i*sstr],
					    k, work) ;
    }
    for ( i = 2 ; i <= depth ; i ++ ) {
      wbfmm_shift_operators_coaxial_SS_init(m->shifts, D, i, 
					    order_s[(i-1)*sstr+0], 
					    k, work) ;
    }

    m->targets = NULL ;
    if ( !precompute_local ) return 0 ;

    g_assert_not_reached() ;

    m->targets = wbfmm_target_list_new(m->tree, nbi_surface_node_number(s)) ;
    wbfmm_target_list_coefficients_init(m->targets, field) ;
    wbfmm_target_list_add_points(m->targets,
				 (gpointer)nbi_surface_node(s,0),
				 NBI_SURFACE_NODE_LENGTH*sizeof(gdouble),
				 nbi_surface_node_number(s)) ;
  
    wbfmm_laplace_target_list_local_coefficients(m->targets, source, work) ;

    return 0 ;
    break ;
  }

  return 0 ;
}
						    
gint NBI_FUNCTION_NAME(nbi_matrix_multiply)(nbi_matrix_t *A,
					    NBI_REAL *x, gint xstr, NBI_REAL al,
					    NBI_REAL *y, gint ystr, NBI_REAL bt,
					    gint nthreads,
					    NBI_REAL *work)

/*
  y := al*A*x + bt*y
*/
  
{
  switch ( A-> problem ) {
  default: g_assert_not_reached() ; break ;
  case NBI_PROBLEM_LAPLACE:
    NBI_FUNCTION_NAME(nbi_matrix_multiply_laplace)(A,
						   x, xstr, al,
						   y, ystr, bt,
						   nthreads,
						   work) ;
    break ;
  case NBI_PROBLEM_HELMHOLTZ:
    NBI_FUNCTION_NAME(nbi_matrix_multiply_helmholtz)(A,
						     x, xstr, al,
						     y, ystr, bt,
						     nthreads,
						     work) ;
    break ;
  }  

  return 0 ;
}

gint nbi_matrix_neighbour_number_max(nbi_matrix_t *m)

{
  gint i, nnmax ;

  nnmax = 0 ;
  for ( i = 0 ; i < nbi_surface_patch_number(m->s) ; i ++ ) {
    nnmax = MAX(nnmax, m->idxp[i+1] - m->idxp[i]) ;
  }  
  
  return nnmax ;
}

static gint boundary_condition_laplace(nbi_surface_t *s,
				       NBI_REAL *p , gint pstr,
				       NBI_REAL *pn, gint nstr,
				       nbi_boundary_condition_t *b)

{
  gdouble *x, *n ;
  gint i ;
  gchar *terms[] = {"rp", "rdp", NULL} ;

  if ( pstr < 1 && p != NULL ) {
    g_error("%s: stride for Laplace problem must be at least 1 (pstr == %d)",
	    __FUNCTION__, pstr) ;
  }
  if ( nstr < 1 && pn != NULL ) {
    g_error("%s: stride for Laplace problem must be at least 1 (nstr == %d)",
	    __FUNCTION__, nstr) ;
  }

  /*check that required boundary condition terms are defined*/
  for ( i = 0 ; terms[i] != NULL ; i ++ ) {
    if ( !nbi_boundary_condition_defined(b, terms[i]) ) {
      fprintf(stderr, "%s: boundary condition \"%s\" not defined\n",
	      __FUNCTION__, terms[i]) ;
      exit(1) ;
    }
  }

  if ( p == NULL && pn == NULL ) return 0 ;
  
  if ( p == NULL && pn != NULL ) {
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      x = (NBI_REAL *)nbi_surface_node(s, i) ;
      n = (NBI_REAL *)nbi_surface_normal(s, i) ;
      
      nbi_boundary_condition_eval(b, x, n) ;
      
      /* p [i*pstr] += nbi_boundary_condition_p (b) ; */
      pn[i*nstr] += nbi_boundary_condition_dp(b) ;
    }

    return 0 ;
  }

  if ( p != NULL && pn == NULL ) {
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      x = (NBI_REAL *)nbi_surface_node(s, i) ;
      n = (NBI_REAL *)nbi_surface_normal(s, i) ;
      
      nbi_boundary_condition_eval(b, x, n) ;
      
      p [i*pstr] += nbi_boundary_condition_p (b) ;
      /* pn[i*nstr] += nbi_boundary_condition_dp(b) ; */
    }

    return 0 ;
  }

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    x = (NBI_REAL *)nbi_surface_node(s, i) ;
    n = (NBI_REAL *)nbi_surface_normal(s, i) ;

    nbi_boundary_condition_eval(b, x, n) ;

    p [i*pstr] += nbi_boundary_condition_p (b) ;
    pn[i*nstr] += nbi_boundary_condition_dp(b) ;
  }
  
  return 0 ;
}

static gint boundary_condition_helmholtz(nbi_surface_t *s,
					 NBI_REAL *p , gint pstr,
					 NBI_REAL *pn, gint nstr,
					 nbi_boundary_condition_t *b)

{
  gdouble *x, *n ;
  gint i ;
  gchar *terms[] = {"rp", "ip", "rdp", "idp", NULL} ;
  
  if ( pstr < 2 && p != NULL ) {
    g_error("%s: stride for Helmholtz problem must be at least 2 (pstr == %d)",
	    __FUNCTION__, pstr) ;
  }
  if ( nstr < 2 && pn != NULL ) {
    g_error("%s: stride for Helmholtz problem must be at least 2 (nstr == %d)",
	    __FUNCTION__, nstr) ;
  }

  /*check that required boundary condition terms are defined*/
  for ( i = 0 ; terms[i] != NULL ; i ++ ) {
    if ( !nbi_boundary_condition_defined(b, terms[i]) ) {
      fprintf(stderr, "%s: boundary condition \"%s\" not defined\n",
	      __FUNCTION__, terms[i]) ;
      exit(1) ;
    }
  }

  if ( p == NULL && pn == NULL ) return 0 ;
  
  if ( p == NULL && pn != NULL ) {
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      x = (NBI_REAL *)nbi_surface_node(s, i) ;
      n = (NBI_REAL *)nbi_surface_normal(s, i) ;
      
      nbi_boundary_condition_eval(b, x, n) ;
      
      /* p [i*pstr+0] += nbi_boundary_condition_p_real(b) ; */
      /* p [i*pstr+1] += nbi_boundary_condition_p_imag(b) ; */
      pn[i*nstr+0] += nbi_boundary_condition_dp_real(b) ;
      pn[i*nstr+1] += nbi_boundary_condition_dp_imag(b) ;
    }

    return 0 ;
  }

  if ( p != NULL && pn == NULL ) {
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      x = (NBI_REAL *)nbi_surface_node(s, i) ;
      n = (NBI_REAL *)nbi_surface_normal(s, i) ;
      
      nbi_boundary_condition_eval(b, x, n) ;
      
      p [i*pstr+0] += nbi_boundary_condition_p_real(b) ;
      p [i*pstr+1] += nbi_boundary_condition_p_imag(b) ;
      /* pn[i*nstr+0] += nbi_boundary_condition_dp_real(b) ; */
      /* pn[i*nstr+1] += nbi_boundary_condition_dp_imag(b) ; */
    }

    return 0 ;
  }

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    x = (NBI_REAL *)nbi_surface_node(s, i) ;
    n = (NBI_REAL *)nbi_surface_normal(s, i) ;

    nbi_boundary_condition_eval(b, x, n) ;

    p [i*pstr+0] += nbi_boundary_condition_p_real(b) ;
    p [i*pstr+1] += nbi_boundary_condition_p_imag(b) ;
    pn[i*nstr+0] += nbi_boundary_condition_dp_real(b) ;
    pn[i*nstr+1] += nbi_boundary_condition_dp_imag(b) ;
  }
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_boundary_condition_set)(nbi_surface_t *s,
						   NBI_REAL *p , gint pstr,
						   NBI_REAL *pn, gint nstr,
						   nbi_boundary_condition_t *b)

{
  if ( nbi_boundary_condition_problem(b) == NBI_PROBLEM_LAPLACE ) {
    return boundary_condition_laplace(s, p, pstr, pn, nstr, b) ;
  }

  if ( nbi_boundary_condition_problem(b) == NBI_PROBLEM_HELMHOLTZ ) {
    return boundary_condition_helmholtz(s, p, pstr, pn, nstr, b) ;
  }

  g_error("%s: unrecognized problem (%d) for boundary condition",
	  __FUNCTION__, nbi_boundary_condition_problem(b)) ;

  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_surface_set_weights)(nbi_surface_t *s)

{
  NBI_REAL K[453*453], ce[453*3], al, bt, *xi, xx[3], work[3*453], J, *st ;
  gint Nk, i3 = 3, xstr, i, j, k, order, nq ;

  nq = nbi_surface_patch_node_number(s, 0) ;
  NBI_FUNCTION_NAME(sqt_quadrature_select)(nq, &st, &order) ;
  Nk = NBI_FUNCTION_NAME(sqt_koornwinder_interp_matrix)(&(st[0]), 3,
							&(st[1]), 3,
							&(st[2]), 3,
							nq, K) ;
  al = 1.0 ; bt = 0.0 ; xstr = NBI_SURFACE_NODE_LENGTH ;
  for ( i = 0 ; i < nbi_surface_patch_number(s) ; i ++ ) {
    /* NBI_FUNCTION_NAME(sqt_quadrature_select)(nq, &st, &order) ; */
    /* sqt_quadrature_select(nq, &st, &order) ; */
    g_assert(nbi_surface_patch_node_number(s, i) == nq) ;
    j = nbi_surface_patch_node(s, i) ;
    xi = (NBI_REAL *)nbi_surface_node(s, j) ;
    blaswrap_dgemm(FALSE, FALSE, nq, i3, nq, al, K, nq, xi, xstr,
		   bt, ce, i3) ;
    for ( k = 0 ; k < nq ; k ++ ) {
      NBI_FUNCTION_NAME(sqt_element_interp)(ce, nq, Nk,
					    st[3*k+0], st[3*k+1], xx,
					    (NBI_REAL *)
					    nbi_surface_normal(s, j+k),
					    &J, NULL, work) ;
      *((NBI_REAL *)nbi_surface_node_weight(s, j+k)) = st[3*k+2]*J ;
    }
  }

  return 0 ;
}
