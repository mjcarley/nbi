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

#include <sqt.h>

#include <blaswrap.h>

#include <nbi.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include "nbi-private.h"


static int upsample_patch(gdouble *xs, gint sstr, gint ns,
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

static void correct_matrix_helmholtz(gdouble k,
				     gdouble *xp, gint pstr, gint np,
				     gdouble *xu, gint ustr, gint nu,
				     gdouble *x , gint xstr,
				     gint *nbrs, gint nnbrs,
				     gdouble *Ku, gdouble *Ast)

{
  gint lda = 4*np, i, j, one = 1, two = 2 ;
  gdouble R, R2, *xf, *xup, r[3], al, C, S, w, G[2], dG[2], rn ;

  for ( i = 0 ; i < nnbrs ; i ++ ) {
    /*field point for correction*/
    xf = &(x[nbrs[i]*xstr]) ;
    for ( j = 0 ; j < nu ; j ++ ) {
      xup = &(xu[j*ustr]) ;
      nbi_vector_diff(r, xf, xup) ; 
      R2 = nbi_vector_length2(r) ;
      /* if ( R2 > NBI_LOCAL_CUTOFF_RADIUS*NBI_LOCAL_CUTOFF_RADIUS ) { */
      /* 	R = sqrt(R2) ; */
      R = sqrt(R2) ;
      if ( R > NBI_LOCAL_CUTOFF_RADIUS ) {
	C = cos(k*R) ; S = sin(k*R) ;
	w = xup[6] ;
	rn = -nbi_vector_scalar(r,&(xup[3])) ;
	G[0] = 0.25*M_1_PI*C/R ;
	G[1] = 0.25*M_1_PI*S/R ;
	dG[0] = -0.25*M_1_PI*(C + k*R*S)*rn/R/R/R ;
	dG[1] =  0.25*M_1_PI*(k*R*C - S)*rn/R/R/R ;
	al = -G[0]*w ;
	blaswrap_daxpy(np, al, &(Ku[j*np]), one, &(Ast[i*lda+0*2*np+0]), two) ;
	al = -G[1]*w ;
	blaswrap_daxpy(np, al, &(Ku[j*np]), one, &(Ast[i*lda+0*2*np+1]), two) ;
	al = -dG[0]*w ;
	blaswrap_daxpy(np, al, &(Ku[j*np]), one, &(Ast[i*lda+2*np+0]), two) ;
	/* al = (k*R*C - S)*w ; */
	al = -dG[1]*w ;
	blaswrap_daxpy(np, al, &(Ku[j*np]), one, &(Ast[i*lda+2*np+1]), two) ;
      }
    }
  }
  
  return ;
}

static void local_correction_matrices_helmholtz(nbi_matrix_t *m,
						gdouble k,
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
  
  lda = 4*nst ;
  nnbrs = idxp[pt+1] - idxp[pt] ;
  nbrs = &(idx[idxp[pt]]) ;
  ip = nbi_surface_patch_node(s, pt) ;
  xs = (NBI_REAL *)nbi_surface_node(s,ip) ;
  xstr = NBI_SURFACE_NODE_LENGTH ;
  /* fprintf(stderr, "patch %d/%d; %d neighbours (%lg)\n", */
  /* 	  pt, nbi_surface_patch_number(s), nnbrs, */
  /* 	  g_timer_elapsed(timer, NULL)) ; */
  g_assert(8*dmax*2*nst*(nnbrs-nst) + 12*nst + 3*nst +
	   2*(nnbrs-nst) <= wsize) ;
  sqt_helmholtz_source_indexed_kw_adaptive(k, xs, xstr, nst, q, nq, K0, NK0,
					   tol, dmax,
					   (gdouble *)nbi_surface_node(s,0),
					   xstr,
					   nbrs, nnbrs-nst,
					   &(Ast[idxp[pt]*lda]), work) ;
  sqt_helmholtz_source_target_kw_self(k, xs, xstr, nst,
				      K0, NK0, N,
				      &(st[0]), 3, &(st[1]), 3,
				      &(Ast[idxp[pt]*lda + (nnbrs-nst)*lda]),
				      work) ;
  /*correct for the upsampled point sources*/
  Ku = nbi_patch_upsample_matrix(nst, nqu) ;
  correct_matrix_helmholtz(k, xs, xstr, nst,
			   &(xu[idxu[pt]*xstr]), xstr, nqu,
			   (gdouble *)nbi_surface_node(s,0), xstr,
			   nbrs, nnbrs, Ku, &(Ast[idxp[pt]*lda])) ;

  return ;
}

static gpointer local_correction_thread_helmholtz(gpointer tdata)

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
  gdouble *st, *K0, *q, tol, k ;

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
  k   = ddata[1] ;
  
  np = nbi_surface_patch_number(m->s) ;

  pt0 = th*(np/nth) ;
  pt1 = (th+1)*(np/nth) ;

  if ( th == nth - 1 ) {
    if ( pt1 < np ) pt1 = np ;
  }
  
  /* fprintf(stderr, "thread %d/%d: %d-%d (%d)\n", */
  /* 	  th, nth, pt0, pt1, np) ; */
  
  for ( pt = pt0 ; pt < pt1 ; pt ++ ) {
    local_correction_matrices_helmholtz(m, k, st, nst, K0, NK0, nqu, q, nq, tol,
					dmax, N, pt, work, wsize)  ;
  }

  return NULL ;
}

nbi_matrix_t *nbi_matrix_assemble_helmholtz(nbi_surface_t *s,
					    gdouble k,
					    gdouble eta,
					    gint nqa, gint dmax,
					    gdouble tol,
					    gint N, gint nu,
					    gint nnmax, gint nthreads)

{
  nbi_matrix_t *m ;
  gint xstr, nst, wstr, wsize, usize, ksize, NK0, nntot, pt ;
  gint ip, order, nnbrs, i, asize, lda ;
  gdouble *st, *K0, r, *qa, *work, *xu ;
  
  m = nbi_matrix_new(s) ;
  m->problem = NBI_PROBLEM_HELMHOLTZ ;

  m->idx  = (gint *)g_malloc0(nbi_surface_patch_number(s)*nnmax*sizeof(gint)) ;
  m->idxp = (gint *)g_malloc0((nbi_surface_patch_number(s)+1)*sizeof(gint)) ;
  m->idxu = (gint *)g_malloc0((nbi_surface_patch_number(s)+1)*sizeof(gint)) ;

  nbi_matrix_wavenumber(m) = k ;
  
  xstr = NBI_SURFACE_NODE_LENGTH ;

  sqt_quadrature_select(nqa, &qa, &order) ;

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
			   
  /* fprintf(stderr, "%s: starting local corrections; t=%lg\n", */
  /* 	  __FUNCTION__, g_timer_elapsed(timer, NULL)) ; */

  /*find neighbours and set up sparse matrix skeleton*/
  memset(m->idx, 0, nbi_surface_patch_number(s)*nnmax*sizeof(gint)) ;

  nntot = 0 ; m->idxu[0] = 0 ;

  nbi_surface_set_patch_data(s) ;

  /* fprintf(stderr, "%s: generating upsampled sources; t=%lg\n", */
  /* 	  __FUNCTION__, g_timer_elapsed(timer, NULL)) ;   */

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
		     nu, &(xu[m->idxu[pt]*xstr]), xstr, work) ;
  }

  lda = 4*nst ;
  asize = nntot*lda ;
  
  m->Ast = g_malloc(asize*sizeof(gdouble)) ;

  /* fprintf(stderr, "%s: generating local corrections; t=%lg\n", */
  /* 	  __FUNCTION__, g_timer_elapsed(timer, NULL)) ;   */
#ifdef _OPENMP
  if ( nthreads == 0 ) {
    for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
      local_correction_matrices_helmholtz(m, k, st, nst, K0, NK0, nu, qa, nqa,
					  tol, dmax, N,	pt, work, wsize) ;
    }
  } else {
    /* g_assert_not_reached() ; /\*get the basic version working first before  */
    /* 			       threading*\/ */
    GThread *threads[NBI_THREAD_NUMBER_MAX] ;
    gpointer data[NBI_THREAD_DATA_SIZE],
      main_data[NBI_THREAD_NUMBER_MAX*NBI_THREAD_MAIN_DATA_SIZE] ;
    gint idata[NBI_THREAD_DATA_INT_SIZE] ;
    gdouble
      ddata[NBI_THREAD_DATA_REAL_SIZE],
      *dpdata[NBI_THREAD_DATA_REAL_PTR_SIZE] ;

    data[NBI_THREAD_DATA_MATRIX] = m ;
    data[NBI_THREAD_DATA_INT] = idata ;
    data[NBI_THREAD_DATA_REAL] = ddata ;
    data[NBI_THREAD_DATA_REAL_POINTER] = dpdata ;

    idata[0] = nst ;
    idata[1] = NK0 ;
    idata[2] = nu ;
    idata[3] = nqa ;
    idata[4] = dmax ;
    idata[5] = N ;
    idata[6] = wstr ;

    dpdata[0] = st ;
    dpdata[1] = K0 ;
    dpdata[2] = qa ;
    ddata[0] = tol ;
    ddata[1] = k ;
    
    for ( i = 0 ; i < nthreads ; i ++ ) {
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_THREAD] =
	GINT_TO_POINTER(i) ;
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_DATA] =
	data ;
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_WORK] =
	&(work[i*wstr]) ;
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_NTHREAD] =
	GINT_TO_POINTER(nthreads) ;
      
      threads[i] = g_thread_new(NULL, local_correction_thread_helmholtz,
				&(main_data[NBI_THREAD_MAIN_DATA_SIZE*i])) ;
    }
    /*make sure all threads complete before we move on*/
    for ( i = 0 ; i < nthreads ; i ++ ) g_thread_join(threads[i]) ;
    /* fprintf(stderr, "threads joined\n") ; */
  }
#else /*_OPENMP*/
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    local_correction_matrices_helmholtz(m, k, st, nst, K0, NK0, nu, qa, nqa,
					tol, dmax, N, pt, work, wsize) ;
  }
#endif /*_OPENMP*/

  return m ;
}

static void correct_matrix_laplace(gdouble *xp, gint pstr, gint np,
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
	al *= nbi_vector_scalar(r,&(xup[3]))/R2 ;
	blaswrap_daxpy(np, al, &(Ku[j*np]), one, &(Ast[i*lda+1*np]), one) ;
      }
    }
  }
  
  return ;
}

static void local_correction_matrices_laplace(nbi_matrix_t *m,
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
  correct_matrix_laplace(xs, xstr, nst,
			 &(xu[idxu[pt]*xstr]), xstr, nqu,
			 (gdouble *)nbi_surface_node(s,0), xstr,
			 nbrs, nnbrs, Ku, &(Ast[idxp[pt]*lda])) ;

  return ;
}

static gpointer local_correction_thread_laplace(gpointer tdata)

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
    local_correction_matrices_laplace(m, st, nst, K0, NK0, nqu, q, nq, tol,
				      dmax, N, pt, work, wsize)  ;
  }

  return NULL ;
}

nbi_matrix_t *nbi_matrix_assemble_laplace(nbi_surface_t *s,
					  gdouble eta,
					  gint nqa, gint dmax,
					  gdouble tol,
					  gint N, gint nu,
					  gint nnmax, gint nthreads)
{
  nbi_matrix_t *m ;
  gint xstr, nst, wstr, wsize, usize, ksize, NK0, nntot, pt ;
  gint ip, order, nnbrs, i, asize, lda ;
  gdouble *st, *K0, r, *qa, *work, *xu ;
  
  m = nbi_matrix_new(s) ;
  m->problem = NBI_PROBLEM_LAPLACE ;

  m->idx  = (gint *)g_malloc0(nbi_surface_patch_number(s)*nnmax*sizeof(gint)) ;
  m->idxp = (gint *)g_malloc0((nbi_surface_patch_number(s)+1)*sizeof(gint)) ;
  m->idxu = (gint *)g_malloc0((nbi_surface_patch_number(s)+1)*sizeof(gint)) ;

  xstr = NBI_SURFACE_NODE_LENGTH ;

  sqt_quadrature_select(nqa, &qa, &order) ;

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
			   
  /* fprintf(stderr, "%s: starting local corrections; t=%lg\n", */
  /* 	  __FUNCTION__, g_timer_elapsed(timer, NULL)) ; */

  /*find neighbours and set up sparse matrix skeleton*/
  memset(m->idx, 0, nbi_surface_patch_number(s)*nnmax*sizeof(gint)) ;

  nntot = 0 ; m->idxu[0] = 0 ;

  nbi_surface_set_patch_data(s) ;

  /* fprintf(stderr, "%s: generating upsampled sources; t=%lg\n", */
  /* 	  __FUNCTION__, g_timer_elapsed(timer, NULL)) ;   */

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
		     nu, &(xu[m->idxu[pt]*xstr]), xstr, work) ;
  }

  lda = 2*nst ;
  asize = nntot*lda ;
  
  m->Ast = g_malloc(asize*sizeof(gdouble)) ;

  /* fprintf(stderr, "%s: generating local corrections; t=%lg\n", */
  /* 	  __FUNCTION__, g_timer_elapsed(timer, NULL)) ;   */
#ifdef _OPENMP
  if ( nthreads == 0 ) {
    for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
      local_correction_matrices_laplace(m, st, nst, K0, NK0, nu, qa, nqa, tol,
					dmax, N, pt, work, wsize) ;
    }
  } else {
    GThread *threads[NBI_THREAD_NUMBER_MAX] ;
    gpointer data[NBI_THREAD_DATA_SIZE],
      main_data[NBI_THREAD_NUMBER_MAX*NBI_THREAD_MAIN_DATA_SIZE] ;
    gint idata[NBI_THREAD_DATA_INT_SIZE] ;
    gdouble
      ddata[NBI_THREAD_DATA_REAL_SIZE],
      *dpdata[NBI_THREAD_DATA_REAL_PTR_SIZE] ;

    data[NBI_THREAD_DATA_MATRIX] = m ;
    data[NBI_THREAD_DATA_INT] = idata ;
    data[NBI_THREAD_DATA_REAL] = ddata ;
    data[NBI_THREAD_DATA_REAL_POINTER] = dpdata ;

    idata[0] = nst ;
    idata[1] = NK0 ;
    idata[2] = nu ;
    idata[3] = nqa ;
    idata[4] = dmax ;
    idata[5] = N ;
    idata[6] = wstr ;

    dpdata[0] = st ;
    dpdata[1] = K0 ;
    dpdata[2] = qa ;
    ddata[0] = tol ;
    
    for ( i = 0 ; i < nthreads ; i ++ ) {
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_THREAD] =
	GINT_TO_POINTER(i) ;
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_DATA] =
	data ;
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_WORK] =
	&(work[i*wstr]) ;
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_NTHREAD] =
	GINT_TO_POINTER(nthreads) ;
      
      threads[i] = g_thread_new(NULL, local_correction_thread_laplace,
				&(main_data[NBI_THREAD_MAIN_DATA_SIZE*i])) ;
    }
    /*make sure all threads complete before we move on*/
    for ( i = 0 ; i < nthreads ; i ++ ) g_thread_join(threads[i]) ;
    /* fprintf(stderr, "threads joined\n") ; */
  }
#else /*_OPENMP*/
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    local_correction_matrices_laplace(m, st, nst, K0, NK0, nu, qa, nqa,
				      tol, dmax, N, pt, work, wsize) ;
  }
#endif /*_OPENMP*/

  return m ;
}
