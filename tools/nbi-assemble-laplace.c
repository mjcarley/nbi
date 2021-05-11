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

nbi_matrix_t *nbi_surface_assemble_matrix(nbi_surface_t *s, gdouble eta,
					  gint nq, gint dmax, gdouble tol,
					  gint N, gint nu,
					  gint nnmax, gint nthreads) ;

static gint upsample_patch(gdouble *xs, gint sstr, gint ns,
			   gint nu, gdouble *xu, gint ustr)

{
  gdouble *K, *qu, ci[453*3], s, t, w, al, bt, work[453*3] ;
  gint Nk, order, i, i3 = 3 ;
  
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
  
  return nu ;
}

static gint correct_matrix(gdouble *xp, gint pstr, gint np,
			   gdouble *xu, gint ustr, gint nu,
			   gdouble *x , gint xstr,
			   gint *nbrs, gint nnbrs,
			   gdouble *Ku, gdouble *Ast)

{
  gint lda = 2*np, i, j, one = 1 ;
  gdouble R, *xf, *xup, r[3], al ;

  for ( i = 0 ; i < nnbrs ; i ++ ) {
    /*field point for correction*/
    xf = &(x[nbrs[i]*xstr]) ;
    for ( j = 0 ; j < nu ; j ++ ) {
      xup = &(xu[j*ustr]) ;
      nbi_vector_diff(r, xf, xup) ; 
      R = nbi_vector_length(r) ;
      if ( R > NBI_LOCAL_CUTOFF_RADIUS ) {
	al = -0.25*M_1_PI/R*xup[6] ;
	blaswrap_daxpy(np, al, &(Ku[j*np]), one, &(Ast[i*lda+0*np]), one) ;
	al *= nbi_vector_scalar(r,&(xup[3]))/R/R ;
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
  gint th = GPOINTER_TO_INT(mdata[0]) ;
  gint nth = GPOINTER_TO_INT(mdata[3]) ;
  gpointer *data = mdata[1] ;
  NBI_REAL *work = mdata[2] ;
  gint *idata = data[0] ;
  gdouble **dpdata = data[1] ;
  gdouble *ddata = data[2] ;
  nbi_matrix_t *m = data[3] ;
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
  
  fprintf(stderr, "thread %d/%d: %d-%d (%d)\n",
  	  th, nth, pt0, pt1, np) ;
  
  for ( pt = pt0 ; pt < pt1 ; pt ++ ) {
    local_correction_matrices(m, st, nst, K0, NK0, nqu, q, nq, tol, dmax, N,
			      pt, work, wsize) ;
  }

  return NULL ;
}

nbi_matrix_t *nbi_surface_assemble_matrix(nbi_surface_t *s, gdouble eta,
					  gint nq, gint dmax, gdouble tol,
					  gint N, gint nu,
					  gint nnmax, gint nthreads)

{
  nbi_matrix_t *m ;
  gint xstr, nst, wstr, wsize, usize, ksize, NK0, nntot, pt ;
  gint ip, order, nnbrs, i, asize, lda ;
  gdouble *st, *K0, r, *q, *work, *xu ;
  
  m = nbi_matrix_new(s) ;
  m->problem = NBI_PROBLEM_LAPLACE ;

  m->idx  = (gint *)g_malloc0(nbi_surface_patch_number(s)*nnmax*sizeof(gint)) ;
  m->idxp = (gint *)g_malloc0((nbi_surface_patch_number(s)+1)*sizeof(gint)) ;
  m->idxu = (gint *)g_malloc0((nbi_surface_patch_number(s)+1)*sizeof(gint)) ;

  xstr = NBI_SURFACE_NODE_LENGTH ;

  sqt_quadrature_select(nq, &q, &order) ;

  nst = nbi_surface_patch_node_number(s,0) ;
  sqt_quadrature_select(nst, &st, &order) ;

  /*size workspace for adaptive quadrature*/
  wstr = 4*dmax*2*nst*nnmax + 12*nst + 3*nst ;    
  wsize = MAX(1,nthreads)*wstr ;
  usize = nu*nbi_surface_patch_number(s)*NBI_SURFACE_NODE_LENGTH ;
  ksize = nst*nst ;
  work = (gdouble *)g_malloc((wsize+ksize)*sizeof(gdouble)) ;
  K0  = &(work[wsize]) ;

  m->xu = g_malloc(usize*sizeof(gdouble)) ;
  m->ustr = NBI_SURFACE_NODE_LENGTH ;
  
  NK0 = sqt_koornwinder_interp_matrix(&(st[0]), 3, &(st[1]), 3, &(st[2]), 3,
				      nst, K0) ;
			   
  fprintf(stderr, "%s: starting local corrections; t=%lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL)) ;

  /*find neighbours and set up sparse matrix skeleton*/
  memset(m->idx, 0, nbi_surface_patch_number(s)*nnmax*sizeof(gint)) ;

  nntot = 0 ; m->idxu[0] = 0 ;

  nbi_surface_set_patch_data(s) ;

  xu = (NBI_REAL *)(m->xu) ;
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    /*loop on patches treated as sources*/
    g_assert(nbi_surface_patch_node_number(s, pt) == nst) ;
    ip = nbi_surface_patch_node(s, pt) ;
    xstr = NBI_SURFACE_NODE_LENGTH ;
    r = *((NBI_REAL *)nbi_surface_patch_sphere_radius(s, pt)) ;
    
    /*find the neighbours*/
    nnbrs = 0 ;
    nbi_surface_patch_neighbours(s, pt, eta*r, 
				 &(m->idx[nntot]), &nnbrs, nnmax-nst) ;
    if ( nnbrs >= nnmax - nst )
      g_error("%s: too many neighbours (%d) for limit (%d)",
	      __FUNCTION__, nnmax-nst, nnbrs) ;
    for ( i = 0 ; i < nbi_surface_patch_node_number(s, pt) ; i ++ )
      m->idx[nntot+nnbrs+i] = ip + i ;
    nnbrs += nbi_surface_patch_node_number(s, pt) ;
    nntot += nnbrs ;
    m->idxp[pt+1] = nntot ;

    /*upsample into xu*/
    m->idxu[pt+1] = m->idxu[pt] +
      upsample_patch((NBI_REAL *)nbi_surface_node(s,ip), xstr, nst,
		     nu, &(xu[m->idxu[pt]*xstr]), xstr) ;
  }

  lda = 2*nst ;
  asize = nntot*lda ;
  
  m->Ast = g_malloc(asize*sizeof(gdouble)) ;

#ifdef _OPENMP
  if ( nthreads == 0 ) {
    for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
      local_correction_matrices(m, st, nst, K0, NK0, nu, q, nq, tol, dmax, N,
				pt, work, wsize) ;
    }
  } else {
    GThread *threads[NBI_THREAD_NUMBER_MAX] ;
    gpointer data[8], main_data[NBI_THREAD_NUMBER_MAX*4] ;
    gint *ipdata[8], idata[16] ;
    gdouble *dpdata[8], ddata[8] ;

    data[0] = idata ;
    data[1] = dpdata ; data[2] = ddata ;
    data[3] = m ;

    idata[0] = nst ;
    idata[1] = NK0 ;
    idata[2] = nu ;
    idata[3] = nq ;
    idata[4] = dmax ;
    idata[5] = N ;
    idata[6] = wstr ;

    dpdata[0] = st ;
    dpdata[1] = K0 ;
    dpdata[2] = q ;
    ddata[0] = tol ;
    
    for ( i = 0 ; i < nthreads ; i ++ ) {
      main_data[4*i+0] = GINT_TO_POINTER(i) ;
      main_data[4*i+1] = data ;
      main_data[4*i+2] = &(work[i*wstr]) ;
      main_data[4*i+3] = GINT_TO_POINTER(nthreads) ;
      
      threads[i] = g_thread_new(NULL, local_correction_thread,
				&(main_data[4*i])) ;
    }
    /*make sure all threads complete before we move on*/
    for ( i = 0 ; i < nthreads ; i ++ ) g_thread_join(threads[i]) ;
    fprintf(stderr, "threads joined\n") ;
  }
#else /*_OPENMP*/
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    local_correction_matrices(m, st, nst, K0, NK0, nu, q, nq,
			      tol, dmax, N, pt, work, wsize) ;
  }
#endif /*_OPENMP*/

  return m ;
}

gint main(gint argc, gchar **argv)

{
  nbi_surface_t *s ;
  nbi_matrix_t *m ;
  gint nqa, dmax, N, *idx, *idxp, nnmax, nqu, *idxu ;
  gdouble r, eta, tol, t ;
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
  
  r = 1.0 ; nqu = 54 ;
  eta = 1.25 ; dmax = 8 ; tol = 1e-12 ; N = 8 ; nqa = 54 ;
  nnmax = 2048 ; 
  gfile = NULL ; mfile = NULL ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  while ( (ch = getopt(argc, argv, "a:d:e:g:m:N:n:r:T:u:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'a': nqa  = atoi(optarg) ; break ;      
    case 'd': dmax = atoi(optarg) ; break ;
    case 'e': tol  = atof(optarg) ; break ;      
    case 'g': gfile = g_strdup(optarg) ; break ;
    case 'm': mfile = g_strdup(optarg) ; break ;
    case 'N': N    = atoi(optarg) ; break ;
    case 'n': eta  = atof(optarg) ; break ;      
    case 'r': r    = atof(optarg) ; break ;
    case 'T': nthreads = atoi(optarg) ; break ;
    case 'u': nqu  = atoi(optarg) ; break ;      
    }
  }

  if ( nthreads < 0 ) nthreads = nproc ;
  
  if ( mfile == NULL ) mfile = g_strdup("matrix.dat") ;
  
  timer = g_timer_new() ;
  
  fprintf(stderr, "%s: initializing geometry r=%lg; t=%lg\n",
	  progname, r, g_timer_elapsed(timer, NULL)) ;
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

  idx = (gint *)g_malloc0(nbi_surface_patch_number(s)*nnmax*sizeof(gint)) ;
  idxp = (gint *)g_malloc0((nbi_surface_patch_number(s)+1)*sizeof(gint)) ;
  idxu = (gint *)g_malloc0((nbi_surface_patch_number(s)+1)*sizeof(gint)) ;
  fprintf(stderr, "%s: starting surface assembly; t=%lg\n",
	  progname, t = g_timer_elapsed(timer, NULL)) ;

  m = nbi_surface_assemble_matrix(s, eta, nqa, dmax, tol, N, nqu, nnmax,
  				  nthreads) ;
  
  output = fopen(mfile, "w") ;
  nbi_matrix_write(output, m) ;
  fclose(output) ;
  
  return 0 ;
}
