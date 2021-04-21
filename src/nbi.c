/* This file is part of NBI, a library for Nystrom Boundary Integral solvers
 *
 * Copyright (C) 2021 Michael Carley
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <glib.h>

#include <sqt.h>

#include <blaswrap.h>

#include <nbi.h>

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
  
  fprintf(f, "%d %d\n",
	  nbi_surface_node_number(s),
	  nbi_surface_patch_number(s)) ;

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
  
  fscanf(f, "%d", &nn) ;
  fscanf(f, "%d", &np) ;

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

  m->problem = 0 ;
  m->s = s ;
  m->fpsize = sizeof(NBI_REAL) ;

  m->ustr = m->pstr = m->nstr = 0 ;
  m->idxu = NULL ;
  m->xu = m->p = m->pn = NULL ;

  m->tree    = NULL ;
  m->targets = NULL ;
  m->shifts  = NULL ;
  
  return m ;
}

static gint read_upsampled_patches(FILE *f, nbi_matrix_t *m)

{
  gint np, i, j ;
  NBI_REAL *xu, *bc ;
  
  fscanf(f, "%d", &np) ;
  fscanf(f, "%d", &(m->ustr)) ;

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
  case NBI_PROBLEM_LAPLACE: m->pstr = m->nstr = 2 ; break ;
  }
  m->bc = g_malloc(m->idxu[np]*(m->pstr)*sizeof(NBI_REAL)) ;
  bc = (NBI_REAL *)(m->bc) ;
  
  m->p  = (gchar *)(&(bc[0])) ;
  m->pn = (gchar *)(&(bc[(m->pstr)/2])) ;
  
  return 0 ;
}

static gint read_correction_matrices(FILE *f, nbi_matrix_t *m)

{
  gint i, j, np, nst ;
  NBI_REAL *Ast ;
  
  fscanf(f, "%d", &np) ;
  fscanf(f, "%d", &nst) ;

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

  Ast = (NBI_REAL *)g_malloc0(nst*2*(m->idxp[np])*sizeof(NBI_REAL)) ;
  m->Ast = (gchar *)Ast ;
  for ( i = 0 ; i < 2*nst*(m->idxp[np]) ; i ++ ) {
    fscanf(f, "%lg", &(Ast[i])) ;
  }  

  return np ;
}

gint NBI_FUNCTION_NAME(nbi_matrix_read)(FILE *input, nbi_matrix_t *m)

{
  read_upsampled_patches(input, m) ;
  read_correction_matrices(input, m) ;

  return 0 ;
}

static gint upsample_sources(nbi_surface_t *s,
			     NBI_REAL *p, gint pstr, NBI_REAL *pn, gint nstr,
			     NBI_REAL *wt, gint wstr, gint *idxu,
			     NBI_REAL *pu, gint pustr, NBI_REAL pwt,
			     NBI_REAL *pnu, gint nustr, NBI_REAL nwt)

{
  gint i, j, pt, ns, nu ;
  gdouble *K, al, bt ;

  al =  1.0 ; bt = 0.0 ;
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    ns = nbi_surface_patch_node_number(s, pt) ;
    nu = idxu[pt+1] - idxu[pt] ;
    K = NBI_FUNCTION_NAME(nbi_patch_upsample_matrix)(ns, nu) ;

    i = nbi_surface_patch_node(s, pt) ;
    j = idxu[pt] ;
    /* if ( pstr != 0 )  */
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(p[i*pstr]), pstr, bt, &(pu[j*pustr]), pustr) ;
    /* if ( nstr != 0 )  */
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(pn[i*pstr]), nstr, bt, &(pnu[j*pustr]), nustr) ;
    for ( i = 0 ; i < nu ; i ++ ) {
      pu [(j+i)*pustr] *= pwt*wt[(j+i)*wstr] ;
      pnu[(j+i)*nustr] *= nwt*wt[(j+i)*wstr] ;
    }
  }
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_matrix_upsample_laplace)(nbi_matrix_t *m,
						    NBI_REAL *p, gint pstr,
						    NBI_REAL *pn, gint nstr,
						    NBI_REAL pwt, NBI_REAL nwt)

{
  NBI_REAL *xu, *su, *sn ;

  su = (NBI_REAL *)(m->p ) ;
  sn = (NBI_REAL *)(m->pn) ;
  xu = (NBI_REAL *)(m->xu) ;
  upsample_sources(m->s, p, pstr, pn, nstr, 
		   &(xu[6]), m->ustr, m->idxu,
		   su, m->pstr, pwt, sn, m->nstr, nwt) ;
  
  return 0 ;
}
						    
