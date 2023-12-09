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

/**
 * @ingroup surfaces
 *
 * @{
 * 
 */

static gint sphere_patch_pair(nbi_surface_t *s, NBI_REAL r,
			      NBI_REAL th0, NBI_REAL th1,
			      NBI_REAL ph0, NBI_REAL ph1,
			      NBI_REAL *st, gint nq)
  
{
  gint np = nbi_surface_patch_number(s) ;
  gint nn = nbi_surface_node_number(s) ;
  
  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nq ;
  g_assert(nn + nq <= nbi_surface_node_number_max(s)) ;

  sqt_patch_nodes_sphere(r, th0, ph0, th0, ph1, th1, ph1,
			 &(st[0]), 3, &(st[1]), 3, &(st[2]), 3, nq,
			 (NBI_REAL *)nbi_surface_node(s,nn),
			 NBI_SURFACE_NODE_LENGTH,
			 (NBI_REAL *)nbi_surface_normal(s,nn),
			 NBI_SURFACE_NODE_LENGTH,
			 (NBI_REAL *)nbi_surface_node_weight(s,nn),
			 NBI_SURFACE_NODE_LENGTH) ;
  nbi_surface_patch_number(s) ++ ;
  np = nbi_surface_patch_number(s) ;
  nn = (nbi_surface_node_number(s) += nq) ;  

  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nq ;
  g_assert(nn + nq <= nbi_surface_node_number_max(s)) ;

  sqt_patch_nodes_sphere(r, th0, ph0, th1, ph1, th1, ph0,
			 &(st[0]), 3, &(st[1]), 3, &(st[2]), 3, nq,
			 (NBI_REAL *)nbi_surface_node(s,nn),
			 NBI_SURFACE_NODE_LENGTH,
			 (NBI_REAL *)nbi_surface_normal(s,nn),
			 NBI_SURFACE_NODE_LENGTH,
			 (NBI_REAL *)
			 nbi_surface_node_weight(s,nn),
			 NBI_SURFACE_NODE_LENGTH) ;			     
  nbi_surface_patch_number(s) ++ ;
  nbi_surface_node_number(s) += nq ;
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_geometry_sphere)(nbi_surface_t *s,
					    NBI_REAL r, gint nth, gint nph,
					    gint nq)

{
  NBI_REAL th0, th1, ph0, ph1, *st ;
  gint i, j, order ;
  
  sqt_quadrature_select(nq, &st, &order) ;

  nbi_surface_node_number(s) = 0 ;
  nbi_surface_patch_number(s) = 0 ;
  
  for ( j = 0 ; j < nph ; j ++ ) {
    ph0 = M_PI*j/nph ; ph1 = M_PI*(j+1)/nph ;
    for ( i = 0 ; i < nth ; i ++ ) {
      th0 = 2.0*M_PI*i/nth ; th1 = 2.0*M_PI*(i+1)/nth ;

      sphere_patch_pair(s, r, th0, th1, ph0, ph1, st, nq) ;  
    }
  }

  nbi_surface_set_weights(s) ;
  
  return 0 ;
}

static gint ellipsoid_patch_pair(nbi_surface_t *s,
				 NBI_REAL a, NBI_REAL b, NBI_REAL c,
				 NBI_REAL th0, NBI_REAL th1,
				 NBI_REAL ph0, NBI_REAL ph1,
				 NBI_REAL *st, gint nq)
  
{
  gint np = nbi_surface_patch_number(s) ;
  gint nn = nbi_surface_node_number(s) ;
  
  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nq ;
  g_assert(nn + nq <= nbi_surface_node_number_max(s)) ;

  NBI_FUNCTION_NAME(sqt_patch_nodes_ellipsoid)(a, b, c,
					       th0, ph0, th0, ph1, th1, ph1,
					       &(st[0]), 3, &(st[1]), 3,
					       &(st[2]), 3, nq,
					       (NBI_REAL *)
					       nbi_surface_node(s,nn),
					       NBI_SURFACE_NODE_LENGTH,
					       (NBI_REAL *)
					       nbi_surface_normal(s,nn),
					       NBI_SURFACE_NODE_LENGTH,
					       (NBI_REAL *)
					       nbi_surface_node_weight(s,nn),
					       NBI_SURFACE_NODE_LENGTH) ;
  nbi_surface_patch_number(s) ++ ;
  np = nbi_surface_patch_number(s) ;
  nn = (nbi_surface_node_number(s) += nq) ;  

  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nq ;
  g_assert(nn + nq <= nbi_surface_node_number_max(s)) ;

  NBI_FUNCTION_NAME(sqt_patch_nodes_ellipsoid)(a, b, c,
					       th0, ph0, th1, ph1, th1, ph0,
					       &(st[0]), 3, &(st[1]), 3,
					       &(st[2]), 3, nq,
					       (NBI_REAL *)
					       nbi_surface_node(s,nn),
					       NBI_SURFACE_NODE_LENGTH,
					       (NBI_REAL *)
					       nbi_surface_normal(s,nn),
					       NBI_SURFACE_NODE_LENGTH,
					       (NBI_REAL *)
					       nbi_surface_node_weight(s,nn),
					       NBI_SURFACE_NODE_LENGTH) ;
  nbi_surface_patch_number(s) ++ ;
  nbi_surface_node_number(s) += nq ;
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_geometry_ellipsoid)(nbi_surface_t *s,
					       NBI_REAL a, NBI_REAL b,
					       NBI_REAL c,
					       gint nth, gint nph,
					       gint nq)

{
  NBI_REAL th0, th1, ph0, ph1, *st ;
  gint i, j, k, order ;

  NBI_FUNCTION_NAME(sqt_quadrature_select)(nq, &st, &order) ;

  nbi_surface_node_number(s) = 0 ;
  nbi_surface_patch_number(s) = 0 ;
  
  for ( j = 0 ; j < nph ; j ++ ) {
    ph0 = acos(1.0 - 2.0*(j+0)/nph) ;
    ph1 = acos(1.0 - 2.0*(j+1)/nph) ;

    for ( i = 0 ; i < nth ; i ++ ) {
      th0 = 2.0*M_PI*i/nth ; th1 = 2.0*M_PI*(i+1)/nth ;

      ellipsoid_patch_pair(s, a, b, c, th0, th1, ph0, ph1, st, nq) ;  
    }
  }
  
  {
    NBI_REAL K[453*453], ce[453*3], al, bt, *xi, xx[3], work[3*453], J ;
    gint Nk, i3 = 3, xstr ;

    Nk = NBI_FUNCTION_NAME(sqt_koornwinder_interp_matrix)(&(st[0]), 3,
							  &(st[1]), 3,
							  &(st[2]), 3,
							  nq, K) ;
    al = 1.0 ; bt = 0.0 ; xstr = NBI_SURFACE_NODE_LENGTH ;
    for ( i = 0 ; i < nbi_surface_patch_number(s) ; i ++ ) {
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
      }
    }
  }
  
  return 0 ;
}

static gint cart2sph(NBI_REAL *x, NBI_REAL *r, NBI_REAL *th, NBI_REAL *ph)

{
  *r = nbi_vector_length(x) ;

  *ph = acos(x[2]/(*r)) ;
  *th = atan2(x[1], x[0]) ;

  if ( *th < 0.0 ) *th += 2.0*M_PI ;
  
  return 0 ;
}

static gint sphere_ico_add_patches(nbi_surface_t *s,
				   NBI_REAL *x0, NBI_REAL *x1, NBI_REAL *x2,
				   gint nr, NBI_REAL *st, gint nst)

{
  if ( nr == 0 ) {
    NBI_REAL xe[9] ;
    gint np = nbi_surface_patch_number(s) ;
    gint nn = nbi_surface_node_number(s) ;

    memcpy(&(xe[3*0]), x0, 3*sizeof(NBI_REAL)) ;
    memcpy(&(xe[3*1]), x1, 3*sizeof(NBI_REAL)) ;
    memcpy(&(xe[3*2]), x2, 3*sizeof(NBI_REAL)) ;

    nbi_surface_patch_node(s, np) = nn ;
    nbi_surface_patch_node_number(s, np) = nst ;
    g_assert(nn + nst <= nbi_surface_node_number_max(s)) ;

    NBI_FUNCTION_NAME(sqt_patch_nodes_tri)(xe, 3, 3, &(st[0]), 3, &(st[1]), 3,
					   &(st[2]), 3,
					   nst, 
					   (NBI_REAL *)nbi_surface_node(s,nn),
					   NBI_SURFACE_NODE_LENGTH,
					   (NBI_REAL *)nbi_surface_normal(s,nn),
					   NBI_SURFACE_NODE_LENGTH,
					   (NBI_REAL *)
					   nbi_surface_node_weight(s,nn),
					   NBI_SURFACE_NODE_LENGTH) ;
    nbi_surface_patch_number(s) ++ ;
    np = nbi_surface_patch_number(s) ;
    nn = (nbi_surface_node_number(s) += nst) ;  

    return 0 ;
  }
 
  /*split and recurse*/
  NBI_REAL x01[3], x12[3], x20[3] ;

  x01[0] = 0.5*(x0[0] + x1[0]) ;
  x01[1] = 0.5*(x0[1] + x1[1]) ;
  x01[2] = 0.5*(x0[2] + x1[2]) ;

  x12[0] = 0.5*(x1[0] + x2[0]) ;
  x12[1] = 0.5*(x1[1] + x2[1]) ;
  x12[2] = 0.5*(x1[2] + x2[2]) ;

  x20[0] = 0.5*(x2[0] + x0[0]) ;
  x20[1] = 0.5*(x2[1] + x0[1]) ;
  x20[2] = 0.5*(x2[2] + x0[2]) ;

  sphere_ico_add_patches(s,  x0, x01, x20, nr-1, st, nst) ;
  sphere_ico_add_patches(s,  x1, x12, x01, nr-1, st, nst) ;
  sphere_ico_add_patches(s,  x2, x20, x12, nr-1, st, nst) ;
  sphere_ico_add_patches(s, x01, x12, x20, nr-1, st, nst) ;
  
  return 0 ;
}

static gint basic_sphere_ico(nbi_surface_t *s, NBI_REAL r, gint nr, gint nq)

/*http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html*/  

{
  NBI_REAL *st, rt, xs[60] ;
  gint i, order ;
  gint faces[] = {0,  11,  5,
		  0 ,  5,  1,
		  0 ,  1,  7,
		  0 ,  7, 10,
		  0 , 10, 11,
		  1 ,  5,  9,
		  5 , 11,  4,
		  11, 10,  2,
		  10,  7,  6,
		  7 ,  1,  8,
		  3 ,  9,  4,
		  3 ,  4,  2,
		  3 ,  2,  6,
		  3 ,  6,  8,
		  3 ,  8,  9,
		  4 ,  9,  5,
		  2 ,  4, 11,
		  6 ,  2, 10,
		  8 ,  6,  7,
		  9 ,  8,  1} ;
  
  sqt_quadrature_select(nq, &st, &order) ;

  nbi_surface_node_number(s) = 0 ;
  nbi_surface_patch_number(s) = 0 ;

  rt = 0.5*(1.0 + sqrt(5.0)) ;

  xs[ 0*3 + 0] = -1.0 ; xs[ 0*3 + 1] =  rt  ; xs[ 0*3 + 2] =  0.0 ; 
  xs[ 1*3 + 0] =  1.0 ; xs[ 1*3 + 1] =  rt  ; xs[ 1*3 + 2] =  0.0 ; 
  xs[ 2*3 + 0] = -1.0 ; xs[ 2*3 + 1] = -rt  ; xs[ 2*3 + 2] =  0.0 ; 
  xs[ 3*3 + 0] =  1.0 ; xs[ 3*3 + 1] = -rt  ; xs[ 3*3 + 2] =  0.0 ; 
  
  xs[ 4*3 + 0] =  0.0 ; xs[ 4*3 + 1] = -1.0 ; xs[ 4*3 + 2] =   rt ; 
  xs[ 5*3 + 0] =  0.0 ; xs[ 5*3 + 1] =  1.0 ; xs[ 5*3 + 2] =   rt ; 
  xs[ 6*3 + 0] =  0.0 ; xs[ 6*3 + 1] = -1.0 ; xs[ 6*3 + 2] =  -rt ; 
  xs[ 7*3 + 0] =  0.0 ; xs[ 7*3 + 1] =  1.0 ; xs[ 7*3 + 2] =  -rt ; 

  xs[ 8*3 + 0] =  rt  ; xs[ 8*3 + 1] =  0.0 ; xs[ 8*3 + 2] =  -1.0 ; 
  xs[ 9*3 + 0] =  rt  ; xs[ 9*3 + 1] =  0.0 ; xs[ 9*3 + 2] =   1.0 ; 
  xs[10*3 + 0] = -rt  ; xs[10*3 + 1] =  0.0 ; xs[10*3 + 2] =  -1.0 ; 
  xs[11*3 + 0] = -rt  ; xs[11*3 + 1] =  0.0 ; xs[11*3 + 2] =   1.0 ; 

  for ( i = 0 ; i < 20 ; i ++ ) {
    sphere_ico_add_patches(s,
			   &(xs[3*faces[3*i+0]]),
			   &(xs[3*faces[3*i+1]]),
			   &(xs[3*faces[3*i+2]]),
			   nr, st, nq) ;
  }

  return 0 ;
}
  
gint NBI_FUNCTION_NAME(nbi_geometry_sphere_ico)(nbi_surface_t *s,
						NBI_REAL r, gint nr, gint nq)
  
{
  gint i ;

  basic_sphere_ico(s, r, nr, nq) ;
  
  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    NBI_REAL *n, *x, th, ph, rho ;
    x = (NBI_REAL *)nbi_surface_node(s, i) ;
    n = (NBI_REAL *)nbi_surface_normal(s, i) ;
    cart2sph(x, &rho, &th, &ph) ;
    n[0] = x[0]/rho ;
    n[1] = x[1]/rho ;
    n[2] = x[2]/rho ;
    x[0] = r*n[0] ; x[1] = r*n[1] ; x[2] = r*n[2] ; 
  }
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_geometry_ellipsoid_ico)(nbi_surface_t *s,
						   NBI_REAL a, NBI_REAL b,
						   NBI_REAL c,
						   gint nr, gint nq)

{
  gint i ;

  basic_sphere_ico(s, 1.0, nr, nq) ;
  
  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    NBI_REAL *n, *x, th, ph, rho, Cp, Sp, Ct, St ;
    x = (NBI_REAL *)nbi_surface_node(s, i) ;
    n = (NBI_REAL *)nbi_surface_normal(s, i) ;
    cart2sph(x, &rho, &th, &ph) ;
    Cp = cos(ph) ; Sp = sin(ph) ;
    Ct = cos(th) ; St = sin(th) ;
    n[0] = b*c*Ct*Sp*Sp ;
    n[1] = a*c*St*Sp*Sp ;
    n[2] = a*b*Sp*Cp ;
    rho = nbi_vector_length(n) ;
    n[0] /= rho ; n[1] /= rho ; n[2] /= rho ;
    x[0] = a*Ct*Sp ;
    x[1] = b*St*Sp ;
    x[2] = c*   Cp ;
  }

  nbi_surface_set_weights(s) ;
  
  return 0 ;
}

static void add_grid_patch(nbi_surface_t *s,
			   NBI_REAL s0, NBI_REAL t0,
			   NBI_REAL s1, NBI_REAL t1,
			   NBI_REAL s2, NBI_REAL t2,
			   NBI_REAL *q, gint nq)

{
  gint i, np, nn, xstr ;
  NBI_REAL *x, *st ;
  
  np = nbi_surface_patch_number(s) ;
  nn = nbi_surface_node_number(s) ;

  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nq ;
  g_assert(nn + nq <= nbi_surface_node_number_max(s)) ;

  x = (NBI_REAL *)nbi_surface_node(s, nn) ;
  xstr = NBI_SURFACE_NODE_LENGTH ;
  
  for ( i = 0 ; i < nq ; i ++ ) {
    st = &(q[3*i]) ;
    x[i*xstr+0] = s0*st[0] + s1*st[1] + s2*(1.0 - st[0] - st[1]) ;
    x[i*xstr+1] = t0*st[0] + t1*st[1] + t2*(1.0 - st[0] - st[1]) ;
    x[i*xstr+6] = st[2] ;
  }
  
  nbi_surface_patch_number(s) ++ ;
  nbi_surface_node_number(s) += nq ;

  return ;
}

gint NBI_FUNCTION_NAME(nbi_geometry_grid)(nbi_surface_t *s,
					  NBI_REAL smin, NBI_REAL smax,
					  gint ns,
					  NBI_REAL tmin, NBI_REAL tmax,
					  gint nt,
					  gint nq)

{
  NBI_REAL *st, s0, t0, s1, t1 ;
  gint i, j, order ;
  
  sqt_quadrature_select(nq, &st, &order) ;

  nbi_surface_node_number(s) = 0 ;
  nbi_surface_patch_number(s) = 0 ;

  for ( i = 0 ; i < ns-1 ; i ++ ) {
    s0 = smin + (smax - smin)*(i+0)/(ns-1) ;
    s1 = smin + (smax - smin)*(i+1)/(ns-1) ;
    for ( j = 0 ; j < nt-1 ; j ++ ) {
      t0 = tmin + (tmax - tmin)*(j+0)/(nt-1) ;
      t1 = tmin + (tmax - tmin)*(j+1)/(nt-1) ;
      g_assert(nbi_surface_patch_number(s)+1 <=
	       nbi_surface_patch_number_max(s)) ;
      g_assert(nbi_surface_node_number(s)+nq <=
	       nbi_surface_node_number_max(s)) ;

      add_grid_patch(s, s0, t0, s1, t0, s1, t1, st, nq) ;

      g_assert(nbi_surface_patch_number(s)+1 <=
	       nbi_surface_patch_number_max(s)) ;
      g_assert(nbi_surface_node_number(s)+nq <=
	       nbi_surface_node_number_max(s)) ;

      add_grid_patch(s, s0, t0, s1, t1, s0, t1, st, nq) ;      
    }
  }
  
  return 0 ;
}
					 

/**
 *
 * @}
 *
 */
