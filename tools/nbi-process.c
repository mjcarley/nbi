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

static void print_help_text(FILE *f, gint offt, gint field, gint offp,
			    gint dmax)

{
  fprintf(f, 
	  "Usage:\n\n"
	  "  %s <options>\n\n",
	  progname) ;

  fprintf(f,
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -d # input data file\n"
	  "  -e # element index offset (%d)\n"
	  "  -f # field element index (%d)\n"
	  "  -g # geometry file\n"
	  "  -p # point index offset (%d)\n"
	  "  -r # recursion depth for triangle generation (%d)\n"
	  "  -v # name of view\n",
	  offt, field, offp, dmax) ;
	  
  return ;
}

gint main(gint argc, gchar **argv)

{
  gchar *gfile, *dfile, ch, *view ;
  nbi_surface_t *s ;
  FILE *input, *output ;
  gint np, fstr, dmax, order, Nk, xistr, ni, fistr, npmax, ntmax ;
  gint ntri, ne, nd, ppe, field, *tri, tstr, offp, offt ;
  NBI_REAL *K, *xi, *fi, tol, *f ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  dmax = 2 ; xistr = 3 ; tol = 0.0 ;
  fi = NULL ; fistr = 0 ; fstr = 2 ;
  offp = 0 ; offt = 0 ;

  /* offp = 100000 ; */

  /*linear output elements, 3 nodes per element*/
  order = 1 ; ppe = 3 ;
  
  input = stdin ; output = stdout ;
  gfile = NULL ; dfile = NULL ;
  view = NULL ;
  
  field = 0 ;

  while ( (ch = getopt(argc, argv, "hd:e:f:g:p:r:v:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'h':
      print_help_text(stderr, offt, field, offp, dmax) ;
      return 0 ;
      break ;
    case 'd': dfile = g_strdup(optarg) ; break ;
    case 'e': offt = atoi(optarg) ; break ;
    case 'f': field = atoi(optarg) ; break ;
    case 'g': gfile = g_strdup(optarg) ; break ;
    case 'r': dmax = atoi(optarg) ; break ;
    case 'p': offp = atoi(optarg) ; break ;
    case 'v': view = g_strdup(optarg) ; break ;
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
    f  = (NBI_REAL *)g_malloc0(sizeof(NBI_REAL)) ; fstr = 0 ;
  } else {  
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
    fclose(input) ;    

    if ( field >= fstr ) {
      fprintf(stderr,
	      "%s: data file has %d fields, cannot extract field %d\n",
	      progname, fstr, field) ;
      exit(1) ;
    }
  }
  
  /*size workspaces: maximum number of triangles per patch*/
  ntri = 2 << (2*dmax) ;

  fistr = fstr ;
  npmax = nbi_surface_patch_number(s)*ntri*ppe ;
  xi = (NBI_REAL *)g_malloc0(npmax*xistr*sizeof(NBI_REAL)) ;
  tstr = 3 ;
  ntmax = nbi_surface_patch_number(s)*ntri ;
  tri = (gint *)g_malloc0(ntmax*tstr*sizeof(gint)) ;
  if ( dfile != NULL ) {
    fi = (NBI_REAL *)g_malloc0(nbi_surface_patch_number(s)*ntri*fistr*ppe*
			       sizeof(NBI_REAL)) ;
  } else {
    fi = (NBI_REAL *)g_malloc0(sizeof(NBI_REAL)) ;
  }

  np = nbi_surface_patch_node(s, 1) - nbi_surface_patch_node(s, 0) ;
  nbi_element_interp_matrix(np, &K, &Nk) ;

  nbi_mesh_triangulate(s, dmax, order, tol, K, Nk,
		       xi, xistr, npmax, tri, tstr, ntmax,
		       &(f[field]), fstr, fi, fistr,
		       &ni, &ne) ;
  
  fprintf(stderr, "%d nodes, %d elements generated\n", ni, ne) ;

  fprintf(stderr, "final point: %d\n", offp+ni) ;
  fprintf(stderr, "final element: %d\n", offt+ne) ;
  
  nbi_mesh_export_gmsh(output, view,
		       xi, xistr, ni, offp,
		       tri, tstr, ne, offt, fi, fistr) ;

  return 0 ;
}
