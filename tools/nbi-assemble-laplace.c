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

#define NBI_THREAD_NUMBER_MAX     8

GTimer *timer ;
gchar *progname ;


static gint upsample_patch(gdouble *xs, gint sstr, gint ns,
			   gint nu, gdouble *xu, gint ustr, gdouble *work)

/*
 * workspace minimum size: 6*ns+2
 */
  
{
  gdouble *K, *qu, s, t, w, al, bt, *ci ;
  gint Nk, order, i, i3 = 3 ;

  ci = &(work[3*ns+1]) ;
  /*check that workspaces have been sized correctly*/
  ci[3*ns] = G_MAXDOUBLE ;
  work[3*ns] = G_MAXDOUBLE ;
  
  nbi_element_interp_matrix(ns, &K, &Nk) ;
  sqt_quadrature_select(nu, &qu, &order) ;

  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemm(FALSE, FALSE, ns, i3, ns, al, K, ns, xs, sstr, bt, ci, i3) ;
  
  for ( i = 0 ; i < nu ; i ++ ) {
    s = qu[3*i+0] ; t = qu[3*i+1] ; w = qu[3*i+2] ;
    sqt_element_interp(ci, ns, Nk, s, t,
		       &(xu[i*ustr+0]),
		       &(xu[i*ustr+3]),
		       &(xu[i*ustr+6]),
		       NULL, work) ;
    xu[i*ustr+6] *= w ;
  }

  /*check final elements are untouched*/
  g_assert(ci[3*ns] == G_MAXDOUBLE) ;
  g_assert(work[3*ns] == G_MAXDOUBLE) ;
  
  return nu ;
}

static gint correct_matrix(gdouble *xp, gint pstr, gint np,
			   gdouble *xu, gint ustr, gint nu,
			   gdouble *x , gint xstr,
			   gint *nbrs, gint nnbrs,
			   gdouble *Ku, gdouble *Ast)

{
  gint lda = 2*np, i, j, one = 1 ;
  gdouble R2, *xf, *xup, r[3], al ;

  for ( i = 0 ; i < nnbrs ; i ++ ) {
    /*field point for correction*/
    xf = &(x[nbrs[i]*xstr]) ;
    for ( j = 0 ; j < nu ; j ++ ) {
      xup = &(xu[j*ustr]) ;
      nbi_vector_diff(r, xf, xup) ; 
      R2 = nbi_vector_length2(r) ;
      if ( R2 > NBI_LOCAL_CUTOFF_RADIUS*NBI_LOCAL_CUTOFF_RADIUS ) {
	al = -0.25*M_1_PI/sqrt(R2)*xup[6] ;
	blaswrap_daxpy(np, al, &(Ku[j*np]), one, &(Ast[i*lda+0*np]), one) ;
	/* al *= nbi_vector_scalar(r,&(xup[3]))/R/R ; */
	al *= nbi_vector_scalar(r,&(xup[3]))/R2 ;
	blaswrap_daxpy(np, al, &(Ku[j*np]), one, &(Ast[i*lda+1*np]), one) ;
      }
    }
  }
  
  return 0 ;
}

static gint local_correction_matrices(nbi_matrix_t *m,
				      gdouble *st, gint nst,
				      gdouble *K0, gint NK0,
				      gint nqu,
				      gdouble *q, gint nq,
				      gdouble tol, gint dmax,
				      gint N,
				      gint pt, 
				      gdouble *work, gint wsize)

{
  gint *nbrs, ip, xstr, nnbrs, lda, *idxu, *idx, *idxp ;
  gdouble *xs, *Ku, *Ast, *xu ;
  nbi_surface_t *s = m->s ;

  Ast = (gdouble *)(m->Ast) ;
  xu  = (gdouble *)(m->xu) ;
  idx = m->idx ;
  idxp = m->idxp ;
  idxu = m->idxu ;
  
  lda = 2*nst ;
  nnbrs = idxp[pt+1] - idxp[pt] ;
  nbrs = &(idx[idxp[pt]]) ;
  ip = nbi_surface_patch_node(s, pt) ;
  xs = (NBI_REAL *)nbi_surface_node(s,ip) ;
  xstr = NBI_SURFACE_NODE_LENGTH ;
  /* fprintf(stderr, "patch %d/%d; %d neighbours (%lg)\n", */
  /* 	  pt, nbi_surface_patch_number(s), nnbrs, */
  /* 	  g_timer_elapsed(timer, NULL)) ; */
  g_assert(4*dmax*2*nst*(nnbrs-nst) + 12*nst + 3*nst +
	   2*(nnbrs-nst) <= wsize) ;
  sqt_laplace_source_indexed_kw_adaptive(xs, xstr, nst, q, nq, K0, NK0,
					 tol, dmax,
					 (gdouble *)nbi_surface_node(s,0),
					 xstr,
					 nbrs, nnbrs-nst,
					 &(Ast[idxp[pt]*lda]), work) ;
  sqt_laplace_source_target_kw_self(xs, xstr, nst,
				    K0, NK0, N,
				    &(st[0]), 3, &(st[1]), 3,
				    &(Ast[idxp[pt]*lda + (nnbrs-nst)*lda]),
				    work) ;
  /*correct for the upsampled point sources*/
  Ku = nbi_patch_upsample_matrix(nst, nqu) ;
  correct_matrix(xs, xstr, nst,
		 &(xu[idxu[pt]*xstr]), xstr, nqu,
		 (gdouble *)nbi_surface_node(s,0), xstr,
		 nbrs, nnbrs, Ku, &(Ast[idxp[pt]*lda])) ;

  return 0 ;
}

static gpointer local_correction_thread(gpointer tdata)

{
  gpointer *mdata = tdata ;
  gint th  = GPOINTER_TO_INT(mdata[NBI_THREAD_MAIN_DATA_THREAD]) ;
  gint nth = GPOINTER_TO_INT(mdata[NBI_THREAD_MAIN_DATA_NTHREAD]) ;
  gpointer *data = mdata[NBI_THREAD_MAIN_DATA_DATA] ;
  NBI_REAL *work = mdata[NBI_THREAD_MAIN_DATA_WORK] ;

  nbi_matrix_t *m = data[NBI_THREAD_DATA_MATRIX] ;
  gint *idata = data[NBI_THREAD_DATA_INT] ;
  NBI_REAL *ddata = data[NBI_THREAD_DATA_REAL] ;
  NBI_REAL **dpdata = data[NBI_THREAD_DATA_REAL_POINTER] ;
  gint NK0, nqu, nq, dmax, N, nst, wsize, np, pt0, pt1, pt ;
  gdouble *st, *K0, *q, tol ;

  nst   = idata[0] ;
  NK0   = idata[1] ;
  nqu   = idata[2] ;
  nq    = idata[3] ;
  dmax  = idata[4] ;
  N     = idata[5] ;
  wsize = idata[6] ;
  
  st = dpdata[0] ;
  K0 = dpdata[1] ;
  q  = dpdata[2] ;
  tol = ddata[0] ;
  
  np = nbi_surface_patch_number(m->s) ;

  pt0 = th*(np/nth) ;
  pt1 = (th+1)*(np/nth) ;

  if ( th == nth - 1 ) {
    if ( pt1 < np ) pt1 = np ;
  }
  
  /* fprintf(stderr, "thread %d/%d: %d-%d (%d)\n", */
  /* 	  th, nth, pt0, pt1, np) ; */
  
  for ( pt = pt0 ; pt < pt1 ; pt ++ ) {
    local_correction_matrices(m, st, nst, K0, NK0, nqu, q, nq, tol, dmax, N,
			      pt, work, wsize)  ;
  }

  return NULL ;
}


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
  gdouble eta, tol, t ;
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
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  while ( (ch = getopt(argc, argv, "ha:d:e:g:m:N:n:T:u:")) != EOF ) {
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
    case 'm': mfile = g_strdup(optarg) ; break ;
    case 'N': N    = atoi(optarg) ; break ;
    case 'n': eta  = atof(optarg) ; break ;      
    case 'T': nthreads = atoi(optarg) ; break ;
    case 'u': nqu  = atoi(optarg) ; break ;      
    }
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

  m = nbi_matrix_assemble_laplace(s, eta, nqa, dmax, tol, N, nqu, nnmax,
					  nthreads) ;
  fprintf(stderr, "%s: matrix assembly complete; t=%lg\n",
	  progname, t = g_timer_elapsed(timer, NULL)) ;

  
  output = fopen(mfile, "w") ;
  nbi_matrix_write(output, m) ;
  fclose(output) ;
  
  return 0 ;
}
