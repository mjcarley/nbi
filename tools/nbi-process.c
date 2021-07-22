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

#include <nbi.h>

#include <wbfmm.h>
#include <sqt.h>

#include <blaswrap.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include "nbi-private.h"

GTimer *timer ;
gchar *progname ;

static gint patch_split_recursion(NBI_REAL *ce, gint ne, gint Nk,
				  NBI_REAL *cf, gint nf,
				  gint d, gint dmax, gint order, NBI_REAL tol,
				  NBI_REAL *st,
				  NBI_REAL *xi, gint xistr,
				  NBI_REAL *fi, gint fistr,
				  gint *ni, gint *nel)

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
		     &(xi[(*ni)*xistr]), xistr) ;
      if ( fi != NULL ) {
	blaswrap_dgemm(FALSE, FALSE, i1, nf, ne, al, K, ne, cf, nf, bt,
		       &(fi[(*ni)*fistr]), fistr) ;
      }
      (*ni) ++ ;
    }

    (*nel) ++ ;
    
    return 0 ;
  }

  sqt_triangle_divide_loop30(st, st0) ;
  patch_split_recursion(ce, ne, Nk, cf, nf, d+1, dmax, order, tol, st0,
			xi, xistr, fi, fistr, ni, nel) ;
  sqt_triangle_divide_loop31(st, st1) ;
  patch_split_recursion(ce, ne, Nk, cf, nf, d+1, dmax, order, tol, st1,
			xi, xistr, fi, fistr, ni, nel) ;
  sqt_triangle_divide_loop32(st, st2) ;
  patch_split_recursion(ce, ne, Nk, cf, nf, d+1, dmax, order, tol, st2,
			xi, xistr, fi, fistr, ni, nel) ;
  sqt_triangle_divide_loop33(st, st3) ;
  patch_split_recursion(ce, ne, Nk, cf, nf, d+1, dmax, order, tol, st3,
			xi, xistr, fi, fistr, ni, nel) ;
  
  return 0 ;
}

static gint split_patch_elements(NBI_REAL *xp, gint xstr, gint np,
				 NBI_REAL *f, gint fstr, gint nf,
				 NBI_REAL *K, gint Nk,
				 gint dmax, gint order, NBI_REAL tol,
				 NBI_REAL *xi, gint xistr,
				 NBI_REAL *fi, gint fistr,
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
			xi, xistr, fi, fistr, ni, ne) ;
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  gchar *gfile, *dfile, ch ;
  nbi_surface_t *s ;
  FILE *input, *output ;
  gint np, pt, xstr, fstr, ip, dmax, order, Nk, xistr, ni, fistr ;
  gint ntri, i, j, ne, nd, ppe, field ;
  NBI_REAL *xp, *K, *xi, *fi, tol, *f ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  dmax = 2 ; xistr = 3 ; tol = 0.0 ;
  fi = NULL ; fistr = 0 ; fstr = 2 ;

  /*linear output elements, 3 nodes per element*/
  order = 1 ; ppe = 3 ;
  
  input = stdin ; output = stdout ;
  gfile = NULL ; dfile = NULL ;

  field = 0 ;

  while ( (ch = getopt(argc, argv, "d:g:r:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'd': dfile = g_strdup(optarg) ; break ;
    case 'f': field = atoi(optarg) ; break ;
    case 'g': gfile = g_strdup(optarg) ; break ;
    case 'r': dmax = atoi(optarg) ; break ;
    }
  }

  if ( gfile != NULL ) {
    input = fopen(gfile, "r") ;
    if ( input == NULL ) {
      fprintf(stderr, "%s: cannot open geometry file %s\n",
	      progname, gfile) ;
      exit(1) ;
    }
  }

  s = nbi_surface_read(input) ;

  if ( gfile != NULL ) fclose(input) ;

  if ( dfile == NULL ) {
    fprintf(stderr, "%s: no data file specified\n", progname) ;
    exit(1) ;
  }
  
  input = fopen(dfile, "r") ;
  if ( input == NULL ) {
    fprintf(stderr, "%s: cannot open data file %s\n",
	    progname, dfile) ;
    exit(1) ;
  }

  f = nbi_data_read(input, &nd, &fstr) ;
  if ( nd < nbi_surface_node_number(s) ) {
    fprintf(stderr,
	    "%s: not enough data points (%d) for surface (%d nodes)\n",
	    progname, nd, nbi_surface_node_number(s)) ;
    exit(1) ;
  }
  
  if ( dfile != NULL ) fclose(input) ;

  if ( field >= fstr ) {
    fprintf(stderr,
	    "%s: data file has %d fields, cannot extract field %d\n",
	    progname, fstr, field) ;
  }
  
  /*size workspaces: maximum number of triangles per patch*/
  ntri = 2 << (2*dmax) ;

  fistr = fstr ;
  xi = (NBI_REAL *)g_malloc0(nbi_surface_patch_number(s)*ntri*xistr*ppe*
			     sizeof(NBI_REAL)) ;
  fi = (NBI_REAL *)g_malloc0(nbi_surface_patch_number(s)*ntri*fistr*ppe*
			     sizeof(NBI_REAL)) ;
  ni = 0 ; ne = 0 ;
  np = nbi_surface_patch_node(s, 1) - nbi_surface_patch_node(s, 0) ;
  nbi_element_interp_matrix(np, &K, &Nk) ;
  xstr = NBI_SURFACE_NODE_LENGTH ;

  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    ip = nbi_surface_patch_node(s, pt) ;
    xp = (NBI_REAL *)nbi_surface_node(s,ip) ;

    split_patch_elements(xp, xstr, np, &(f[ip*fstr]), fstr, fstr, K, Nk,
			 dmax, order, tol, xi, xistr, fi, fistr, &ni, &ne) ;
  }
  
  fprintf(stderr, "%d nodes, %d elements generated\n", ni, ne) ;

  fprintf(output,
	  "$MeshFormat\n"
	  "4.1 0 8\n"
	  "$EndMeshFormat\n"
	  "$Nodes\n") ;
  fprintf(output, "1 %d 1 %d\n", ni, ni) ;
  fprintf(output, "2 1 0 %d\n", ni) ;
  for ( i = 0 ; i < ni ; i ++ )
    fprintf(output, "%d\n", i+1) ;
  for ( i = 0 ; i < ni ; i ++ ) {
    fprintf(output, "%lg %lg %lg\n",
	    xi[i*xistr+0], xi[i*xistr+1], xi[i*xistr+2]) ;
  }
  fprintf(output,
	  "$EndNodes\n"
	  "$Elements\n") ;
  fprintf(output, "1 %d 1 %d\n", ne, ne) ;
  fprintf(output, "2 1 2 %d\n", ne) ;
  for ( i = 0 ; i < ne ; i ++ ) {
    fprintf(output, "%d", i+1) ;
    for ( j = 0 ; j < ppe ; j ++ ) {
      fprintf(output, " %d", i*ppe+j+1) ;
    }
    fprintf(output, "\n") ;
  }
  fprintf(output, "$EndElements\n") ;

  if ( f == NULL ) return 0 ;

  fprintf(output, "$NodeData\n") ;
  fprintf(output,
	  "1\n"
	  "\"NBI view\"\n"
	  "1\n"
	  "0.0\n"
	  "3\n"
	  "0\n"
	  "1\n"
	  "%d\n", ni) ;
  for ( i = 0 ; i < ni ; i ++ )
    fprintf(stdout, "%d %lg\n", i+1, fi[i*fistr]) ;

  fprintf(output, "$EndNodeData\n") ;  
  
  return 0 ;
}
