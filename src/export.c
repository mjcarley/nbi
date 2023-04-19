/* This file is part of NBI, a library for Nystrom Boundary Integral solvers
 *
 * Copyright (C) 2023 Michael Carley
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

#include <nbi.h>

#include <sqt.h>
#include <blaswrap.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include "nbi-private.h"

static gint patch_split_recursion(NBI_REAL *ce, gint ne, gint Nk,
				  NBI_REAL *cf, gint nf,
				  gint d, gint dmax, gint order, NBI_REAL tol,
				  NBI_REAL *st,
				  NBI_REAL *xi, gint xistr,
				  NBI_REAL *fi, gint fistr,
				  gint *tri, gint tstr,
				  gint *np, gint *nel,
				  gint npmax, gint ntmax)

{
  NBI_REAL st0[6], st1[6], st2[6], st3[6], s, t, K[453], al, bt ;
  gint i3 = 3, i1 = 1, i ;

  if ( d == dmax ) {
    /*generate elements*/
    al = 1.0 ; bt = 0.0 ;
    for ( i = 0 ; i < 3 ; i ++ ) {
      s = st[2*i+0] ; t = st[2*i+1] ;
      sqt_koornwinder_nm(Nk, s, t, K, 1, ne) ;
      blaswrap_dgemm(FALSE, FALSE, i1, i3, ne, al, K, ne, ce, i3, bt,
		     &(xi[(*np)*xistr]), xistr) ;
      if ( fi != NULL ) {
	blaswrap_dgemm(FALSE, FALSE, i1, nf, ne, al, K, ne, cf, nf, bt,
		       &(fi[(*np)*fistr]), fistr) ;
      }
      tri[3*(*nel)+i] = (*np) ;
      (*np) ++ ;
      if ( (*np) >= npmax ) {
	g_error("%s: not enough space allocated (%d) for number "
		"of points (%d)", __FUNCTION__, npmax, *np) ;
      }
    }
    
    (*nel) ++ ;
    if ( (*nel) >= ntmax ) {
      g_error("%s: not enough space allocated (%d) for number "
		"of triangles (%d)", __FUNCTION__, ntmax, *nel) ;
    }
    
    return 0 ;
  }

  sqt_triangle_divide_loop30(st, st0) ;
  patch_split_recursion(ce, ne, Nk, cf, nf, d+1, dmax, order, tol, st0,
			xi, xistr, fi, fistr, tri, tstr, np, nel,
			npmax, ntmax) ;
  sqt_triangle_divide_loop31(st, st1) ;
  patch_split_recursion(ce, ne, Nk, cf, nf, d+1, dmax, order, tol, st1,
			xi, xistr, fi, fistr, tri, tstr, np, nel,
			npmax, ntmax) ;
  sqt_triangle_divide_loop32(st, st2) ;
  patch_split_recursion(ce, ne, Nk, cf, nf, d+1, dmax, order, tol, st2,
			xi, xistr, fi, fistr, tri, tstr, np, nel,
			npmax, ntmax) ;
  sqt_triangle_divide_loop33(st, st3) ;
  patch_split_recursion(ce, ne, Nk, cf, nf, d+1, dmax, order, tol, st3,
			xi, xistr, fi, fistr, tri, tstr, np, nel,
			npmax, ntmax) ;
  
  return 0 ;
}

static gint split_patch_elements(NBI_REAL *xp, gint xstr, gint np,
				 NBI_REAL *f, gint fstr, gint nf,
				 NBI_REAL *K, gint Nk,
				 gint dmax, gint order, NBI_REAL tol,
				 NBI_REAL *xi, gint xistr,
				 NBI_REAL *fi, gint fistr,
				 gint *tri, gint tstr,
				 gint npmax, gint ntmax,
				 gint *ni, gint *ne)

{
  NBI_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0} ;
  NBI_REAL ce[3*453], fe[8*543], al, bt ;
  gint i3 = 3 ;

  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, np, i3, np, al, K, np, xp, xstr, bt, ce, i3) ;
  if ( fi != NULL ) {
    blaswrap_dgemm(FALSE, FALSE, np, nf, np, al, K, np, f, fstr, bt, fe, nf) ;
  }
  
  patch_split_recursion(ce, np, Nk, fe, nf, 0, dmax, order, tol, st,
			xi, xistr, fi, fistr, tri, tstr, ni, ne,
			npmax, ntmax) ;
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_mesh_triangulate)(nbi_surface_t *s,
					     gint dmax, gint order,
					     NBI_REAL tol,
					     NBI_REAL *K, gint Nk,
					     NBI_REAL *x, gint xstr,
					     gint npmax,
					     gint *tri, gint tstr,
					     gint ntmax,
					     NBI_REAL *f, gint fstr,
					     NBI_REAL *fi, gint fistr,
					     gint *np, gint *nt)

{
  gint pt, str, ip, ppp ;
  NBI_REAL *xp ;
  
  *np = *nt = 0 ;
  /*points per patch*/
  ppp = nbi_surface_patch_node(s, 1) - nbi_surface_patch_node(s, 0) ;

  str = NBI_SURFACE_NODE_LENGTH ;

  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    ip = nbi_surface_patch_node(s, pt) ;
    xp = (NBI_REAL *)nbi_surface_node(s,ip) ;

    split_patch_elements(xp, str, ppp, &(f[ip*fstr]), fstr, fstr, K, Nk,
			 dmax, order, tol, x, xstr, fi, fistr,
			 tri, tstr, npmax, ntmax, np, nt) ;
  }
  
  return 0 ;
}  

gint NBI_FUNCTION_NAME(nbi_mesh_export_gmsh)(FILE *f, gchar *view,
					     NBI_REAL *x, gint xstr, gint np,
					     gint offp,
					     gint *tri, gint tstr, gint nt,
					     gint offt,
					     NBI_REAL *data, gint dstr)

{
  gint i, j ;
  gchar *vn ;

  if ( view == NULL ) {
    vn = g_strdup("NBI view") ;
  } else {
    vn = view ;
  }
  
  fprintf(f,
	  "$MeshFormat\n"
	  "4.1 0 8\n"
	  "$EndMeshFormat\n"
	  "$Nodes\n") ;
  fprintf(f, "1 %d 1 %d\n", np, np) ;
  fprintf(f, "2 1 0 %d\n", np) ;
  for ( i = 0 ; i < np ; i ++ )
    fprintf(f, "%d\n", offp + i + 1) ;
  for ( i = 0 ; i < np ; i ++ ) {
    fprintf(f, "%lg %lg %lg\n", x[i*xstr+0], x[i*xstr+1], x[i*xstr+2]) ;
  }
  fprintf(f,
	  "$EndNodes\n"
	  "$Elements\n") ;
  fprintf(f, "1 %d 1 %d\n", nt, nt) ;
  fprintf(f, "2 1 2 %d\n", nt) ;
  for ( i = 0 ; i < nt ; i ++ ) {
    fprintf(f, "%d", offt + i + 1) ;
    for ( j = 0 ; j < 3 ; j ++ ) {
      fprintf(f, " %d", offp + tri[i*tstr + j]+1) ;
    }
    fprintf(f, "\n") ;
  }
  fprintf(f, "$EndElements\n") ;

  if ( data == NULL ) return 0 ;

  fprintf(f, "$NodeData\n") ;
  fprintf(f,
	  "1\n"
	  "\"%s\"\n"
	  "1\n"
	  "0.0\n"
	  "3\n"
	  "0\n"
	  "1\n"
	  "%d\n", vn, np) ;
  for ( i = 0 ; i < np ; i ++ )
    fprintf(stdout, "%d %lg\n", offp+i+1, data[i*dstr]) ;

  fprintf(f, "$EndNodeData\n") ;

  return 0 ;
}
