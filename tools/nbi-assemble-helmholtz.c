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

#include <sqt.h>

#include <blaswrap.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include "nbi-private.h"

GTimer *timer ;
gchar *progname ;

static void print_help_text(FILE *f, gint nqa, gint dmax, gdouble tol,
			    gint N, gdouble eta, gint nthreads, gint nqu)


{
  fprintf(f, 
	  "Usage:\n\n"
	  "  %s <options>\n\n",
	  progname) ;

  fprintf(f,
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -a # number of nodes in adaptive quadrature rule (%d)\n"
	  "  -d # maximum depth of adaptive quadrature (%d)\n"
	  "  -e # tolerance for adaptive quadrature (%lg)\n"
	  "  -g # input geometry file\n"
	  "  -k # wavenumber\n"
	  "  -m # output matrix file\n"
	  "  -N # order of singular quadrature (%d)\n"
	  "  -n # separation parameter for selection of quadratures (%lg)\n"
	  "  -T # number of threads (%d)\n"
	  "  -u # number of nodes in upsample quadrature rule (%d)\n",
	  nqa, dmax, tol, N, eta, nthreads, nqu) ;
  return ;
}

gint main(gint argc, gchar **argv)

{
  nbi_surface_t *s ;
  nbi_matrix_t *m ;
  gint nqa, dmax, N, nnmax, nqu ;
  gdouble eta, tol, t, k ;
  FILE *output, *input ;
  gchar ch, *mfile, *gfile ;
  gint nthreads, nproc ;
  
  nthreads = 1 ;

#ifdef _OPENMP
  nproc = g_get_num_processors() ;
  nthreads = -1 ;
#else  /*_OPENMP*/
  nproc = 1 ;
#endif /*_OPENMP*/

  output = stdout ; input = stdin ;

  nqu = 54 ;
  eta = 1.25 ; dmax = 8 ; tol = 1e-12 ; N = 8 ; nqa = 54 ;
  nnmax = 8192 ; 
  gfile = NULL ; mfile = NULL ;

  k = 0.0 ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  while ( (ch = getopt(argc, argv, "ha:d:e:g:k:m:N:n:T:u:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'h':
      print_help_text(stderr, nqa, dmax, tol, N, eta, nthreads, nqu) ;
      return 0 ;
      break ;
    case 'a': nqa  = atoi(optarg) ; break ;
    case 'd': dmax = atoi(optarg) ; break ;
    case 'e': tol  = atof(optarg) ; break ;
    case 'g': gfile = g_strdup(optarg) ; break ;
    case 'k': k = atof(optarg) ; break ;
    case 'm': mfile = g_strdup(optarg) ; break ;
    case 'N': N    = atoi(optarg) ; break ;
    case 'n': eta  = atof(optarg) ; break ;      
    case 'T': nthreads = atoi(optarg) ; break ;
    case 'u': nqu  = atoi(optarg) ; break ;      
    }
  }

  if ( k == 0.0 ) {
    fprintf(stderr, "%s: wavenumber must be set for Helmholtz problem"
	    " (use option -k)", progname) ;
    return 1 ;
  }
  
  if ( nthreads < 0 ) nthreads = nproc ;
  
  if ( mfile == NULL ) mfile = g_strdup("matrix.dat") ;
  
  timer = g_timer_new() ;
  
  fprintf(stderr, "%s: reading geometry; t=%lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
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

  fprintf(stderr, "%s: starting matrix assembly; t=%lg\n",
	  progname, t = g_timer_elapsed(timer, NULL)) ;

  m = nbi_surface_assemble_matrix_helmholtz(s, k, eta, nqa, dmax, tol,
					    N, nqu, nnmax,
					    nthreads) ;
  fprintf(stderr, "%s: matrix assembly complete; t=%lg\n",
	  progname, t = g_timer_elapsed(timer, NULL)) ;

  
  output = fopen(mfile, "w") ;
  nbi_matrix_write(output, m) ;
  fclose(output) ;
  
  return 0 ;
}
