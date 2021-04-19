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

  s->xc  = (NBI_REAL *)
    g_malloc0(nnmax*NBI_SURFACE_NODE_LENGTH*sizeof(NBI_REAL)) ;
  s->pcr = (NBI_REAL *)
    g_malloc0(npmax*NBI_SURFACE_PATCH_DATA_LENGTH*sizeof(NBI_REAL)) ;
  s->ip  = (gint *)g_malloc0(NBI_SURFACE_PATCH_LENGTH*npmax*sizeof(gint)) ;

  s->fpsize = sizeof(NBI_REAL) ;
  
  return s ;
}

gint NBI_FUNCTION_NAME(nbi_surface_write)(nbi_surface_t *s, FILE *f)

{
  gint i, j ;
  
  fprintf(f, "%d %d\n",
	  nbi_surface_node_number(s),
	  nbi_surface_patch_number(s)) ;

  for ( i = 0 ; i < nbi_surface_patch_number(s) ; i ++ ) {
    fprintf(f, "%d %d\n",
	    nbi_surface_patch_node(s, i),
	    nbi_surface_patch_node_number(s, i)) ;
  }

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    for ( j = 0 ; j < NBI_SURFACE_NODE_LENGTH ; j ++ )
      fprintf(f, " %1.16e", nbi_surface_node_element(s,i,j)) ;
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

nbi_surface_t *NBI_FUNCTION_NAME(nbi_surface_read)(FILE *f)

{
  nbi_surface_t *s ;
  gint nn, np, i, j ;

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
    for ( j = 0 ; j < NBI_SURFACE_NODE_LENGTH ; j ++ )
      fscanf(f, "%lg", &(nbi_surface_node_element(s,i,j))) ;
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

  cp = nbi_surface_patch_centre(s, p) ;
  rp = nbi_surface_patch_sphere_radius(s, p) ;
  
  x = nbi_surface_node(s, 0) ;
  for ( i = 0 ; i < p ; i ++ ) {
    ci = nbi_surface_patch_centre(s, i) ;
    ri = nbi_surface_patch_sphere_radius(s, i) ;
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
    ci = nbi_surface_patch_centre(s, i) ;
    ri = nbi_surface_patch_sphere_radius(s, i) ;
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
    ip = nbi_surface_patch_node(s, i) ;
    NBI_FUNCTION_NAME(nbi_surface_patch_centroid)(nbi_surface_node(s,ip),
						 NBI_SURFACE_NODE_LENGTH,
						 &(nbi_surface_node_weight(s, ip)),
						 NBI_SURFACE_NODE_LENGTH,
						 nbi_surface_patch_node_number(s, i),
						 nbi_surface_patch_centre(s, i)) ;
    nbi_surface_patch_sphere_radius(s, i) =
      NBI_FUNCTION_NAME(nbi_surface_patch_radius)(nbi_surface_node(s,ip),
						 NBI_SURFACE_NODE_LENGTH, 
						 nbi_surface_patch_node_number(s, i),
						 nbi_surface_patch_centre(s, i)) ;
			       
  }    
  
  return 0 ;
}
