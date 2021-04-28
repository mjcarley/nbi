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

#include "nbi-private.h"

GTimer *timer ;
gchar *progname ;

gint nbi_gmres_real(nbi_matrix_t *A, 
		    gdouble *x, gint xstr, gdouble *b, gint bstr,
		    gint m, gint max_it, gdouble tol, gdouble *error,
		    gdouble *work) ;

static gint make_sources(nbi_surface_t *s,
			 gdouble *xs, gdouble q,
			 gdouble *p , gint pstr,
			 gdouble *pn, gint nstr)

{
  gdouble R, Rn, r[3], *x, *n, G ;
  gint i ;

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    x = (NBI_REAL *)nbi_surface_node(s, i) ;
    n = (NBI_REAL *)nbi_surface_normal(s, i) ;

    nbi_vector_diff(r, x, xs) ;
    R = nbi_vector_length(r) ;
    Rn = nbi_vector_scalar(r,n)/R ;
    G = q*0.25*M_1_PI/R ;
    p [i*pstr] += G ;
    pn[i*nstr] -= G*Rn/R ;    
  }
  
  return 0 ;
}

static gdouble greens_function_laplace(gdouble *x, gdouble *y)

{
  gdouble G, R ;

  R = nbi_vector_distance(x, y) ;
  G = 0.25*M_1_PI/R ;
  
  return G ;
}


gint main(gint argc, gchar **argv)

{
  nbi_surface_t *s ;
  nbi_matrix_t *matrix ;
  gdouble xs[512], *f, *xp, *src, t ;
  gdouble emax, fmax, G ;
  gdouble dtree, pwt, nwt ;
  FILE *output, *input ;
  gchar ch, *gfile, *mfile ;
  gdouble *work, tol ;
  gint fmm_work_size, nqfmm, order_fmm, order_inc, i, fstr, solver_work_size ;
  gint gmres_max_iter, gmres_restart ;
  guint depth, order[48] = {0}, order_s, order_r, order_max ;
  gboolean fmm, shift_bw, greens_id, layer_potentials ;
  
  output = stdout ;
  mfile = NULL ; gfile = NULL ;
  
  dtree = 1e-2 ; fmm = FALSE ; nqfmm = 1 ; shift_bw = TRUE ;
  greens_id = FALSE ; layer_potentials = FALSE ;
  
  pwt = 2.0 ; nwt = 2.0 ;
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  order_fmm = 12 ; order_inc = 2 ; depth = 4 ;

  solver_work_size = 0 ;
  gmres_max_iter = 128 ; gmres_restart = 40 ; tol = 1e-9 ;
    
  fstr = 3 ;
  while ( (ch = getopt(argc, argv, "d:fGg:Lm:o:T:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'd': order_inc = atoi(optarg) ; break ;
    case 'f': fmm = TRUE ; break ;
    case 'G': greens_id = TRUE ; break ;
    case 'g': gfile = g_strdup(optarg) ; break ;
    case 'L': layer_potentials = TRUE ; break ;
    case 'm': mfile = g_strdup(optarg) ; break ;
    case 'o': order_fmm = atoi(optarg) ; break ;
    case 'T': depth = atoi(optarg) ; break ;
    }
  }

  if ( gfile == NULL ) gfile = g_strdup("geometry.dat") ;
  if ( mfile == NULL ) mfile = g_strdup("matrix.dat") ;
  
  fprintf(stderr, "%s: reading geometry from %s\n", progname, gfile) ;
  input = fopen(gfile, "r") ;

  s = nbi_surface_read(input) ;

  fclose(input) ;

  matrix = nbi_matrix_new(s) ;
  matrix->problem = NBI_PROBLEM_LAPLACE ;
  
  fprintf(stderr, "%s: reading matrix from %s\n", progname, mfile) ;
  input = fopen(mfile, "r") ;

  nbi_matrix_read(input, matrix) ;
  
  fclose(input) ;

  timer = g_timer_new() ;

  fprintf(stderr,
	  "%s: geometry initialized, %d nodes, %d patches; t=%lg\n",
	  progname, nbi_surface_node_number(s), nbi_surface_patch_number(s),
	  g_timer_elapsed(timer, NULL)) ;

  /*boundary point sources*/
  fprintf(stderr, "%s: making surface sources; t=%lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  src = (gdouble *)g_malloc0(nbi_surface_node_number(s)*2*sizeof(gdouble)) ;
  f   = (gdouble *)g_malloc0(nbi_surface_node_number(s)*fstr*sizeof(gdouble)) ;

  xs[0] = 0.3 ; xs[1] = -0.4 ; xs[2] = 0.2 ;

  make_sources(s, xs, 1.0, &(src[0]), 2, &(src[1]), 2) ;

  if ( !greens_id && !layer_potentials ) {
    /*solver settings for GMRES*/
    i = nbi_surface_node_number(s) ;
        
    solver_work_size = 2*i + 4*gmres_restart +
      (gmres_restart+1)*(i+gmres_restart) + 2 ;
  }
  
  if ( fmm ) {
    order_s = order_fmm ; order_r = order_fmm ;
    order[2*depth+0] = order_s ; 
    order[2*depth+1] = order_r ; 
    order_max = MAX(order_s, order_r) ;
    for ( i = depth-1 ; i > 0 ; i -- ) {
      order[2*i+0] = order[2*(i+1)+0] ;
      order[2*i+0] = order[2*(i+1)+0] + order_inc ;
      order[2*i+1] = order[2*(i+1)+1] + order_inc ;
      order_max = MAX(order_max, order[2*i+0]) ;
      order_max = MAX(order_max, order[2*i+1]) ;
    }
    fmm_work_size = wbfmm_element_number_rotation(2*(order_max+2)) ;
    fmm_work_size = MAX(fmm_work_size, (order_max+1)*(order_max+1)*nqfmm*16) ;

    fmm_work_size += solver_work_size ;
    
    work = (gdouble *)g_malloc0(fmm_work_size*sizeof(gdouble)) ;
    wbfmm_laplace_coaxial_translate_init(order_max+1) ;    

    fprintf(stderr, "%s: building tree; t=%lg\n",
	    progname, t = g_timer_elapsed(timer, NULL)) ;
    wbfmm_shift_angle_table_init() ;

    nbi_matrix_fmm_init(matrix, NBI_PROBLEM_LAPLACE,
			NULL, &(order[0]), 2, &(order[1]), 2,
			depth, dtree, shift_bw, work) ;    
  } else {
    fmm_work_size = 16384 ;
    work = (gdouble *)g_malloc0(fmm_work_size*sizeof(gdouble)) ;    
  }
  
  if ( greens_id ) {
    fprintf(stderr, "%s: evaluating Green's identity; t=%lg\n",
	    progname, t = g_timer_elapsed(timer, NULL)) ;
    nbi_surface_greens_identity_laplace(matrix,
					&(src[0]), 2, pwt,
					&(src[1]), 2, nwt,
					f, fstr, work) ;
    fprintf(stderr, "%s: surface integration complete; t=%lg (%lg)\n",
	    progname,
	    g_timer_elapsed(timer, NULL), g_timer_elapsed(timer, NULL) - t) ;
    
    emax = 0.0 ; fmax = 0.0 ;
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      xp = (NBI_REAL *)nbi_surface_node(s,i) ;
      G = greens_function_laplace(xp, xs) ;
      fmax = MAX(fmax, G) ;
      emax = MAX(emax, fabs(G - f[i*fstr])) ;
      fprintf(output, "%lg %lg %lg %lg %lg\n",
	      xp[0], xp[1], xp[2], f[i*fstr], fabs(G - f[i*fstr])) ;
    }
    
    fprintf(stderr, "L_inf norm: %lg\n", emax/fmax) ;

    return 0 ;
  }

  if ( layer_potentials ) {
    fprintf(stderr, "%s: evaluating double-layer potential; t=%lg\n",
    	    progname, t = g_timer_elapsed(timer, NULL)) ;
    matrix->potential = NBI_POTENTIAL_DOUBLE ;
    nbi_matrix_multiply(matrix, &(src[0]), 2, 1.0, f, fstr, 0.0, work) ;
    fprintf(stderr, "%s: evaluating single-layer potential; t=%lg\n",
    	    progname, t = g_timer_elapsed(timer, NULL)) ;
    matrix->potential = NBI_POTENTIAL_SINGLE ;
    nbi_matrix_multiply(matrix, &(src[1]), 2, -2.0, f, fstr, 2.0, work) ;
    fprintf(stderr, "%s: surface integration complete; t=%lg (%lg)\n",
	    progname,
	    g_timer_elapsed(timer, NULL), g_timer_elapsed(timer, NULL) - t) ;
    
    emax = 0.0 ; fmax = 0.0 ;
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      xp = (NBI_REAL *)nbi_surface_node(s,i) ;
      G = greens_function_laplace(xp, xs) ;
      fmax = MAX(fmax, G) ;
      emax = MAX(emax, fabs(G - f[i*fstr])) ;
      fprintf(output, "%lg %lg %lg %lg %lg\n",
	      xp[0], xp[1], xp[2], f[i*fstr], fabs(G - f[i*fstr])) ;
    }
    
    fprintf(stderr, "L_inf norm: %lg\n", emax/fmax) ;

    return 0 ;
  }

  /*if we get to here, we're doing a solve*/
  gdouble *rhs, *p, error = 0.0 ;

  rhs = (gdouble *)g_malloc0(nbi_surface_node_number(s)*sizeof(gdouble)) ;
  p   = (gdouble *)g_malloc0(nbi_surface_node_number(s)*sizeof(gdouble)) ;
  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) p[i] = 1.0 ;
  
  /*form right hand side*/
  matrix->diag = 0.0 ;
  matrix->potential = NBI_POTENTIAL_SINGLE ;
  fprintf(stderr, "%s: forming right hand side, t=%lg\n", progname,
	  t=g_timer_elapsed(timer, NULL)) ;

  nbi_matrix_multiply(matrix, &(src[1]), 2, 1.0, rhs, 1, 0.0, work) ;
  
  matrix->diag = -0.5 ;
  matrix->potential = NBI_POTENTIAL_DOUBLE ;
  fprintf(stderr, "%s: starting solver\n", progname) ;

  i = nbi_gmres_real(matrix, p, 1, rhs, 1, gmres_restart, gmres_max_iter, tol,
		     &error, work) ;

  fprintf(stderr, "%s: %d iterations; error = %lg, t=%lg (%lg)\n",
	  progname, i, error,
	  g_timer_elapsed(timer, NULL),
	  g_timer_elapsed(timer, NULL) - t) ;

  emax = fmax = 0.0 ;
  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    emax = MAX(emax,fabs(p[i]-src[2*i])) ;
    fmax = MAX(fmax, fabs(src[2*i])) ;
  }

  fprintf(stderr, "%s: L_inf norm %lg\n", progname, emax/fmax) ;
  
  return 0 ;
}
