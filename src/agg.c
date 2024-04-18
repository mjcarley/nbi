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

#ifdef HAVE_AGG
#include <agg.h>

static void shapefunc(gdouble s, gdouble t, gdouble L[])

{
  L[0] = 1.0 - s - t ;
  L[1] =       s     ;
  L[2] =           t ;
    
  return ;
}

static void element_interp(agg_mesh_t *m, gint surf, gint *nodes,
			   gdouble *s, gint sstr, gdouble *t, gint tstr,
			   gint nst, gdouble *x, gint xstr,
			   agg_surface_workspace_t *w)

{
  gdouble L[3], si, ti, st[6], u, v, *p, xc[3] ;
  agg_surface_t *S ;
  agg_patch_t *P ;
  gint i, tag ;

  S = agg_mesh_surface(m, surf) ;
  P = agg_mesh_patch(m, surf) ;
  tag = agg_mesh_point_tag(m, nodes[0]) ;

  if ( tag >= 0 ) {
    for ( i = 0 ; i < 3 ; i ++ ) {
      g_assert(agg_mesh_point_tag(m, nodes[i]) == tag) ;
      st[2*i+0] = agg_mesh_point_s(m, nodes[i]) ;
      st[2*i+1] = agg_mesh_point_t(m, nodes[i]) ;
      agg_patch_map(P, st[2*i+0], st[2*i+1], &u, &v) ;
      agg_surface_point_eval(S, u, v, xc, w) ;

      p = agg_mesh_point(m, nodes[i]) ;
      if ( fabs(p[0]-xc[0]) > 1e-3 ||
	   fabs(p[1]-xc[1]) > 1e-3 ||
	   fabs(p[2]-xc[2]) > 1e-3 )
	g_error("%s: incorrect point (%d) found at (s,t)=(%lg, %lg) "
		"(%lg, %lg, %lg)=/=(%lg, %lg, %lg)",
		__FUNCTION__, nodes[i], st[2*i+0], st[2*i+1],
		xc[0], xc[1], xc[2], p[0], p[1], p[2]) ;
    }

    for ( i = 0 ; i < nst ; i ++ ) {
      /* shapefunc(s[i*sstr], t[i*tstr], L) ; */
      /* si = st[2*0+0]*L[0] + st[2*1+0]*L[1] + st[2*2+0]*L[2] ; */
      /* ti = st[2*0+1]*L[0] + st[2*1+1]*L[1] + st[2*2+1]*L[2] ; */
      agg_patch_triangle_interp(P,
				st[2*0+0], st[2*0+1],
				st[2*1+0], st[2*1+1],
				st[2*2+0], st[2*2+1],
				s[i*sstr], t[i*tstr], &si, &ti) ;
      agg_patch_map(P, si, ti, &u, &v) ;
      agg_surface_point_eval(S, u, v, &(x[i*xstr]), w) ;    
    }

    return ;
  }

  /*dealing with a surface blend*/
  tag = -1 - tag ;
  for ( i = 0 ; i < 3 ; i ++ ) {
    g_assert(agg_mesh_point_tag(m, nodes[i]) == -1-tag) ;
    st[2*i+0] = agg_mesh_point_s(m, nodes[i]) ;
    st[2*i+1] = agg_mesh_point_t(m, nodes[i]) ;
    agg_surface_blend_evaluate(agg_mesh_surface_blend(m,tag),
			       st[2*i+0], st[2*i+1], xc, w) ;
    p = agg_mesh_point(m, nodes[i]) ;
    if ( fabs(p[0]-xc[0]) > 1e-3 ||
	 fabs(p[1]-xc[1]) > 1e-3 ||
	 fabs(p[2]-xc[2]) > 1e-3 )
      g_error("%s: incorrect point (%d) found at (s,t)=(%lg, %lg) "
	      "(%lg, %lg, %lg)=/=(%lg, %lg, %lg)",
	      __FUNCTION__, nodes[i], st[2*i+0], st[2*i+1],
	      xc[0], xc[1], xc[2], p[0], p[1], p[2]) ;
  }

  for ( i = 0 ; i < nst ; i ++ ) {
    shapefunc(s[i*sstr], t[i*tstr], L) ;
    si = st[2*0+0]*L[0] + st[2*1+0]*L[1] + st[2*2+0]*L[2] ;
    ti = st[2*0+1]*L[0] + st[2*1+1]*L[1] + st[2*2+1]*L[2] ;
    agg_surface_blend_evaluate(agg_mesh_surface_blend(m,tag), si, ti,
			       &(x[i*xstr]), w) ;
  }
  
  return ;
}

static void add_nbi_triangle(nbi_surface_t *s, agg_mesh_t *m, gint *nodes,
			     gint surf, agg_surface_workspace_t *w,
			     gdouble *st, gint nst)

{
  gint i, np, nn, i3 = 3, xstr ;
  gdouble *x, work[3*453], K[454*454], ci[453*453], N ;
  gdouble al, bt ;
  
  np = nbi_surface_patch_number(s) ;
  nn = nbi_surface_node_number(s) ;
  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nst ;
  
  N = sqt_koornwinder_interp_matrix(&(st[0]), 3, &(st[1]), 3, &(st[2]), 3,
				    nst, K) ;

  /*generate the mesh nodes*/
  al = 1.0 ; bt = 0.0 ; xstr = NBI_SURFACE_NODE_LENGTH ;
  x = (gdouble *)nbi_surface_node(s, nn) ;
  element_interp(m, surf, nodes, &(st[0]), 3, &(st[1]), 3, nst, x, xstr, w) ;
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

static void add_nbi_patch(nbi_surface_t *s, agg_mesh_t *m,
			  gint *nodes, gint nnodes,
			  gint surf, agg_surface_workspace_t *w,
			  gdouble *st, gint nst)

{
  gint ntri[3] ;
  
  if ( nnodes == 3 ) {
    add_nbi_triangle(s, m, nodes, surf, w, st, nst) ;
    
    return ;
  }

  if ( nnodes == 4 ) {
    ntri[0] = nodes[0] ; ntri[1] = nodes[1] ; ntri[2] = nodes[2] ; 
    add_nbi_triangle(s, m, ntri, surf, w, st, nst) ;
    
    ntri[0] = nodes[0] ; ntri[1] = nodes[2] ; ntri[2] = nodes[3] ;
    add_nbi_triangle(s, m, ntri, surf, w, st, nst) ;

    return ;
  }

  g_assert_not_reached() ;
  
  return ;
}

nbi_surface_t *nbi_agg_mesh(gchar *file, gint nq)

{
  nbi_surface_t *s = NULL ;
  gint i, pps, nodes[4], nnodes ;
  gint surf, order ;
  gdouble *st ;
  agg_body_t *b ;
  agg_surface_workspace_t *w ;
  agg_mesh_t *m ;
  
  g_assert(sizeof(NBI_REAL) == sizeof(gdouble)) ;
  pps = 2 ;
  
  b = agg_body_new(64, 64) ;
  w = agg_surface_workspace_new() ;
  g_debug("%s: reading AGG body", __FUNCTION__) ;
  agg_body_read(b, file, FALSE) ;
  g_debug("%s: compiling AGG globals", __FUNCTION__) ;
  agg_body_globals_compile(b) ;
  g_debug("%s: evaluating AGG globals", __FUNCTION__) ;
  agg_body_globals_eval(b) ;

  g_debug("%s: making new mesh", __FUNCTION__) ;
  m = agg_mesh_new(65536, 65536, 65536) ;
  g_assert(m != NULL) ;
  g_debug("%s: generating AGG mesh", __FUNCTION__) ;
  agg_mesh_body(m, b, pps, w) ;
  g_debug("%s: AGG mesh generated", __FUNCTION__) ;
  s = nbi_surface_alloc(nq*2*agg_mesh_element_number(m),
			2*agg_mesh_element_number(m)) ;
  sqt_quadrature_select(nq, &st, &order) ;
  
  g_debug("%s: generating NBI surface", __FUNCTION__) ;
  for ( i = 0 ; i < agg_mesh_element_number(m) ; i ++ ) {
    agg_mesh_element_nodes(m, i, nodes, &nnodes, &surf) ;
    add_nbi_patch(s, m, nodes, nnodes, surf, w, st, nq) ;
  }
  g_debug("%s: NBI surface generated", __FUNCTION__) ;
  
  return s ;
}

#endif /*HAVE_AGG*/
