/* This file is part of NBI, a library for Nystrom Boundary Integral solvers
 *
 * Copyright (C) 2022 Michael Carley
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

#include <wbfmm.h>

#include <blaswrap.h>

#include <nbi.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include "nbi-private.h"

#ifdef HAVE_AGG

#include <agg.h>

static void add_agg_grid_patch(agg_body_t *b,
			       agg_grid_t *g, gint t,
			       agg_workspace_t *w,
			       agg_parser_t *p,
			       gdouble *st, gint nst,
			       nbi_surface_t *s)

{
  gint i, np, nn, i3 = 3, xstr ;
  gdouble *x, work[3*453], K[454*454], ci[453*453], N ;
  gdouble al, bt, u, v ;
  
  np = nbi_surface_patch_number(s) ;
  nn = nbi_surface_node_number(s) ;
  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nst ;

  g_assert(nn + nst <= nbi_surface_node_number_max(s)) ;

  N = sqt_koornwinder_interp_matrix(&(st[0]), 3, &(st[1]), 3, &(st[2]), 3,
				    nst, K) ;

  for ( i = 0 ; i < nst ; i ++ ) {
    x = (gdouble *)nbi_surface_node(s, nn+i) ;
    agg_grid_element_interpolate(g, t, st[3*i+0], st[3*i+1], &u, &v) ;
    agg_body_point_eval(b, u, v, x, w) ;
  }
  al = 1.0 ; bt = 0.0 ; xstr = NBI_SURFACE_NODE_LENGTH ;
  x = (gdouble *)nbi_surface_node(s, nn) ;
  blaswrap_dgemm(FALSE, FALSE, nst, i3, nst, al, K, nst, x, xstr, bt, ci, i3) ;

  for ( i = 0 ; i < nst ; i ++ ) {
    x = (gdouble *)nbi_surface_node(s, nn+i) ;
    sqt_element_interp(ci, nst, N, st[3*i+0], st[3*i+1],
		       &(x[0]), &(x[3]), &(x[6]), NULL, work) ;
    x[6] *= st[3*i+2] ;
  }

  nbi_surface_node_number(s) += nst ;
  nbi_surface_patch_number(s) ++ ;
  
  return ;
}

nbi_surface_t *nbi_agg_mesh(gint fid, gint nq)

{
  nbi_surface_t *s = NULL ;
  gint i, j, order, npts, ntri ;
  agg_parser_t *p ;
  agg_body_t *b ;
  agg_crowd_t c ;
  agg_distribution_t *d ;
  agg_workspace_t *w ;
  GScanner *scanner ;
  gdouble *st ;
  
  g_assert(sizeof(NBI_REAL) == sizeof(gdouble)) ;
  
  p = agg_parser_alloc() ;
  b = agg_body_alloc() ;
  scanner = agg_scanner_alloc() ;
  g_scanner_input_file(scanner, fid) ;

  c.p = p ;

  agg_parser_crowd_read(scanner, &c) ;

  fprintf(stderr, "%d bod%s in crowd\n",
	  agg_crowd_body_number(&c),
	  (agg_crowd_body_number(&c) == 1 ? "y" : "ies")) ;

  npts = ntri = 0 ;
  for ( i = 0 ; i < agg_crowd_body_number(&c) ; i ++ ) {
    if ( c.b[i]->g != NULL ) {
      fprintf(stderr, "  body %d: %d grid points; %d triangles\n", i,
	      agg_grid_point_number_max(c.b[i]->g),
	      agg_grid_triangle_number_max(c.b[i]->g)) ;
      npts += agg_grid_point_number_max(c.b[i]->g) ;
      ntri += agg_grid_triangle_number_max(c.b[i]->g) ;
    }
  }

  fprintf(stderr, "   total: %d grid points; %d triangles\n", npts, ntri) ;

  agg_parser_expressions_evaluate(c.p) ;

  w = agg_workspace_alloc(32) ;
  
  s = nbi_surface_alloc(nq*ntri, ntri) ;
  sqt_quadrature_select(nq, &st, &order) ;

  nbi_surface_node_number(s) = 0 ;
  nbi_surface_patch_number(s) = 0 ;

  for ( i = 0 ; i < agg_crowd_body_number(&c) ; i ++ ) {
    fprintf(stderr, "meshing body %d\n", i) ;
    b = agg_crowd_body(&c, i) ;
    for ( j = 0 ; j < agg_body_distribution_number(b) ; j ++ ) {
      d = agg_body_distribution(b,j) ;
      agg_distribution_interpolation_weights(d) ;
    }
    if ( b->g != NULL ) {
      for ( j = 0 ; j < agg_grid_triangle_number(b->g) ; j ++ ) {
        add_agg_grid_patch(b, b->g, j, w, p, st, nq, s) ;
      }
    }
  }
  
  return s ;
}

#endif /*HAVE_AGG*/
