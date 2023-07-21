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

#include <wbfmm.h>
#include <sqt.h>

#include <blaswrap.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include "nbi-private.h"

GTimer *timer ;
gchar *progname ;

static void print_help_text(FILE *f, gint field)

{
  fprintf(f, 
	  "%s: compute Helmholtz potential field around a surface\n\n"
	  "Usage:\n\n"
	  "  %s <options>\n\n",
	  progname, progname) ;

  fprintf(f,
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -d # input data file\n"
	  "  -F # field geometry file\n"
	  "  -f # field element index (%d)\n"
	  "  -g # geometry file\n"
	  "  -k # wavenumber\n"
	  "  -S use surface geometry and data, and do not compute field\n",
	  /* "  -r # recursion depth for triangle generation (%d)\n", */
	  field) ;
	  
  return ;
}

gint main(gint argc, gchar **argv)

{
  gchar *gfile, *dfile, ch, *ffile, *bfile ;
  nbi_surface_t *s, *sf ;
  nbi_boundary_condition_t *bc ;
  gboolean surface_data ;
  FILE *input, *output ;
  gint xstr, i, field, nd, fstr ;
  NBI_REAL *f, k, x[3], p[2], al[] = {1,0}, bt[] = {0,0} ;
  NBI_REAL *pf, *xg ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  fstr = 2 ;

  input = stdin ; output = stdout ;
  bfile = NULL ; dfile = NULL ; ffile = NULL ; gfile = NULL ;
  
  field = 0 ; k = 0.0 ; bc = NULL ;
  surface_data = FALSE ;
  
  while ( (ch = getopt(argc, argv, "hb:d:F:f:g:k:S")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'h':
      print_help_text(stderr, field) ;
      return 0 ;
      break ;
    case 'b': bfile = g_strdup(optarg) ; break ;
    case 'd': dfile = g_strdup(optarg) ; break ;
    case 'F': ffile = g_strdup(optarg) ; break ;
    case 'f': field = atoi(optarg) ; break ;
    case 'g': gfile = g_strdup(optarg) ; break ;
    case 'k': k = atof(optarg) ; break ;
    case 'S': surface_data = TRUE ; break ;
    /* case 'r': dmax = atoi(optarg) ; break ; */
    }
  }

  if ( gfile == NULL ) {
    fprintf(stderr, "%s: no geometry file specified (use -g option)\n",
	    progname) ;
    exit(1) ;
  }

  if ( dfile == NULL ) {
    fprintf(stderr, "%s: no surface data file specified (use -d option)\n",
	    progname) ;
    exit(1) ;
  }
  
  input = fopen(gfile, "r") ;
  if ( input == NULL ) {
    fprintf(stderr, "%s: cannot open geometry file %s\n",
	    progname, gfile) ;
    exit(1) ;
  }

  fprintf(stderr, "%s: geometry %s\n", progname, gfile) ;
  
  s = nbi_surface_read(input) ;

  fclose(input) ;
  
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

  if ( fstr < 4 ) {
    fprintf(stderr,
	    "%s: not enough surface data per node (%d); at least four "
	    "required\n", progname, fstr) ;
    exit(1) ;
  }
  
  if ( field >= fstr ) {
    fprintf(stderr,
	    "%s: data file has %d fields, cannot extract field %d\n",
	    progname, fstr, field) ;
    exit(1) ;
  }

  if ( bfile != NULL ) {
    fprintf(stderr, "%s: reading boundary conditions from %s\n",
	    progname, bfile) ;
    if ( (input = fopen(bfile, "r")) == NULL ) {
      fprintf(stderr, "%s: cannot open boundary condition file %s\n",
	      progname, bfile) ;
      return 1 ;
    }
    
    bc = nbi_boundary_condition_new(NBI_PROBLEM_HELMHOLTZ) ;
    
    nbi_boundary_condition_read(input, bc) ;

    /*so that the field calculation increments the incident field*/
    bt[0] = 1.0 ; bt[1] = 0.0 ;
    
    fclose(input) ;
  }
  
  /*read points from stdin and stream output*/
  if ( ffile == NULL && !surface_data ) {
    input = stdin ;
    while ( fscanf(input, "%lg %lg %lg", &(x[0]), &(x[1]), &(x[2])) != EOF ) {

      if ( bc != NULL ) {
	g_assert_not_reached() ; /*need to write a function to do this*/
      }
      
      nbi_surface_field_helmholtz(s, k, &(f[field]), fstr, al, bt, x, p) ;
      
      fprintf(stdout, "%e %e %e %e %e\n", x[0], x[1], x[2], p[0], p[1]) ;
    }

    return 0 ;
  }

  if ( !surface_data ) {
    input = fopen(ffile, "r") ;
    if ( input == NULL ) {
      fprintf(stderr, "%s: cannot open field geometry file %s\n",
	      progname, ffile) ;
      exit(1) ;
    }
    
    sf = nbi_surface_read(input) ;
    
    fclose(input) ;
  } else {
    sf = s ;
  }

  pf = (NBI_REAL *)g_malloc0(nbi_surface_node_number(sf)*sizeof(NBI_REAL)*2) ;
  xg = (NBI_REAL *)nbi_surface_node(sf,0) ;
  xstr = NBI_SURFACE_NODE_LENGTH ;

  if ( bc != NULL ) {
    nbi_boundary_condition_set(sf, &(pf[0]), 2, NULL, 0, bc) ;
  }

  if ( !surface_data ) {
    for ( i = 0 ; i < nbi_surface_node_number(sf) ; i ++ ) {
      nbi_surface_field_helmholtz(s, k, &(f[field]), fstr, al, bt,
				  &(xg[i*xstr]), &(pf[i*2])) ;    
    }
  } else {
    for ( i = 0 ; i < nbi_surface_node_number(sf) ; i ++ ) {
      pf[2*i+0] += f[field+i*fstr+0] ;
      pf[2*i+1] += f[field+i*fstr+1] ;
    }    
  }

  output = stdout ;

  nbi_data_write(output, pf, 2, 2, nbi_surface_node_number(sf)) ;
  
  return 0 ;
}
