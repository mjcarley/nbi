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

#include <blaswrap.h>

#include <sqt.h>

#include <nbi.h>

static gint sphere_patch_pair(nbi_surface_t *s, gdouble r,
			      gdouble th0, gdouble th1,
			      gdouble ph0, gdouble ph1,
			      gdouble *st, gint nq)
  
{
  gint np = nbi_surface_patch_number(s) ;
  gint nn = nbi_surface_node_number(s) ;
  
  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nq ;
  g_assert(nn + nq <= nbi_surface_node_number_max(s)) ;

  sqt_patch_nodes_sphere(r, th0, ph0, th0, ph1, th1, ph1,
			 &(st[0]), 3, &(st[1]), 3, &(st[2]), 3, nq,
			 nbi_surface_node(s,nn),
			 NBI_SURFACE_NODE_LENGTH,
			 nbi_surface_normal(s,nn),
			 NBI_SURFACE_NODE_LENGTH,
			 &(nbi_surface_node_weight(s,nn)),
			 NBI_SURFACE_NODE_LENGTH) ;
  nbi_surface_patch_number(s) ++ ;
  np = nbi_surface_patch_number(s) ;
  nn = (nbi_surface_node_number(s) += nq) ;  

  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nq ;
  g_assert(nn + nq <= nbi_surface_node_number_max(s)) ;

  sqt_patch_nodes_sphere(r, th0, ph0, th1, ph1, th1, ph0,
			 &(st[0]), 3, &(st[1]), 3, &(st[2]), 3, nq,
			 nbi_surface_node(s,nn),
			 NBI_SURFACE_NODE_LENGTH,
			 nbi_surface_normal(s,nn),
			 NBI_SURFACE_NODE_LENGTH,
			 &(nbi_surface_node_weight(s,nn)),
			 NBI_SURFACE_NODE_LENGTH) ;			     
  nbi_surface_patch_number(s) ++ ;
  nbi_surface_node_number(s) += nq ;
  
  return 0 ;
}

gint nbi_geometry_sphere(nbi_surface_t *s, gdouble r, gint nth, gint nph,
			 gint nq)

{
  gdouble th0, th1, ph0, ph1, *st ;
  gint i, j, k, order ;
  
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

  {
    gdouble K[453*453], ce[453*3], al, bt, *xi, xx[3], work[3*453], J ;
    gint Nk, i3 = 3, xstr ;

    Nk = sqt_koornwinder_interp_matrix(&(st[0]), 3, &(st[1]), 3, &(st[2]), 3,
				     nq, K) ;
    al = 1.0 ; bt = 0.0 ; xstr = NBI_SURFACE_NODE_LENGTH ;
    for ( i = 0 ; i < nbi_surface_patch_number(s) ; i ++ ) {
      j = nbi_surface_patch_node(s, i) ;
      xi = nbi_surface_node(s, j) ;
      blaswrap_dgemm(FALSE, FALSE, nq, i3, nq, al, K, nq, xi, xstr,
		     bt, ce, i3) ;
      for ( k = 0 ; k < nq ; k ++ ) {
      sqt_element_interp(ce, nq, Nk, st[3*k+0], st[3*k+1], xx,
			 nbi_surface_normal(s, j+k), &J, NULL, work) ;
      }
    }
  }
  
  return 0 ;
}

static gint ellipsoid_patch_pair(nbi_surface_t *s,
				 gdouble a, gdouble b, gdouble c,
				 gdouble th0, gdouble th1,
				 gdouble ph0, gdouble ph1,
				 gdouble *st, gint nq)
  
{
  gint np = nbi_surface_patch_number(s) ;
  gint nn = nbi_surface_node_number(s) ;
  
  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nq ;
  g_assert(nn + nq <= nbi_surface_node_number_max(s)) ;

  sqt_patch_nodes_ellipsoid(a, b, c, th0, ph0, th0, ph1, th1, ph1,
			    &(st[0]), 3, &(st[1]), 3, &(st[2]), 3, nq,
			    nbi_surface_node(s,nn),
			    NBI_SURFACE_NODE_LENGTH,
			    nbi_surface_normal(s,nn),
			    NBI_SURFACE_NODE_LENGTH,
			    &(nbi_surface_node_weight(s,nn)),
			    NBI_SURFACE_NODE_LENGTH) ;
  nbi_surface_patch_number(s) ++ ;
  np = nbi_surface_patch_number(s) ;
  nn = (nbi_surface_node_number(s) += nq) ;  

  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nq ;
  g_assert(nn + nq <= nbi_surface_node_number_max(s)) ;

  sqt_patch_nodes_ellipsoid(a, b, c, th0, ph0, th1, ph1, th1, ph0,
			    &(st[0]), 3, &(st[1]), 3, &(st[2]), 3, nq,
			    nbi_surface_node(s,nn),
			    NBI_SURFACE_NODE_LENGTH,
			    nbi_surface_normal(s,nn),
			    NBI_SURFACE_NODE_LENGTH,
			    &(nbi_surface_node_weight(s,nn)),
			    NBI_SURFACE_NODE_LENGTH) ;			     
  nbi_surface_patch_number(s) ++ ;
  nbi_surface_node_number(s) += nq ;
  
  return 0 ;
}

gint nbi_geometry_ellipsoid(nbi_surface_t *s,
			    gdouble a, gdouble b, gdouble c,
			    gint nth, gint nph,
			    gint nq)

{
  gdouble th0, th1, ph0, ph1, *st ;
  gint i, j, k, order ;
  
  sqt_quadrature_select(nq, &st, &order) ;

  nbi_surface_node_number(s) = 0 ;
  nbi_surface_patch_number(s) = 0 ;
  
  for ( j = 0 ; j < nph ; j ++ ) {
    /* ph0 = M_PI*j/nph ; ph1 = M_PI*(j+1)/nph ; */
    ph0 = acos(1.0 - 2.0*(j+0)/nph) ;
    ph1 = acos(1.0 - 2.0*(j+1)/nph) ;

    for ( i = 0 ; i < nth ; i ++ ) {
      th0 = 2.0*M_PI*i/nth ; th1 = 2.0*M_PI*(i+1)/nth ;

      ellipsoid_patch_pair(s, a, b, c, th0, th1, ph0, ph1, st, nq) ;  
    }
  }
  
  {
    gdouble K[453*453], ce[453*3], al, bt, *xi, xx[3], work[3*453], J ;
    gint Nk, i3 = 3, xstr ;

    Nk = sqt_koornwinder_interp_matrix(&(st[0]), 3, &(st[1]), 3, &(st[2]), 3,
				     nq, K) ;
    al = 1.0 ; bt = 0.0 ; xstr = NBI_SURFACE_NODE_LENGTH ;
    for ( i = 0 ; i < nbi_surface_patch_number(s) ; i ++ ) {
      j = nbi_surface_patch_node(s, i) ;
      xi = nbi_surface_node(s, j) ;
      blaswrap_dgemm(FALSE, FALSE, nq, i3, nq, al, K, nq, xi, xstr,
		     bt, ce, i3) ;
      for ( k = 0 ; k < nq ; k ++ ) {
      sqt_element_interp(ce, nq, Nk, st[3*k+0], st[3*k+1], xx,
			 nbi_surface_normal(s, j+k), &J, NULL, work) ;
      }
    }
  }
  
  return 0 ;
}
