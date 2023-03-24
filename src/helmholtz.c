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

#include <wbfmm.h>

#include <blaswrap.h>

#include <nbi.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include "nbi-private.h"
#define wbfmm_tree_point_index(_t,_i)			\
  ((NBI_REAL *)(&((_t)->points[(_i)*((_t)->pstr)])))

static void point_source_field_helmholtz(NBI_REAL k,
					 NBI_REAL *xs, gint xstr, gint ns,
					 NBI_REAL *p , gint pstr, NBI_REAL pwt,
					 NBI_REAL *pn, gint nstr, NBI_REAL nwt,
					 NBI_REAL *x, NBI_REAL wt, NBI_REAL *f,
					 NBI_REAL *al, NBI_REAL *bt)

/*
 * f := bt*f + al*SUM(pwt*dG*\phi - nwt*G*d\phi)
 *
 * pwt: weighting for potential source (double layer potential)
 * nwt: weighting for normal derivative source (single layer potential)
 */
  
{
  gint i ;  
  NBI_REAL R, r[3], C, S, rn, G[2], dG[2], df[2] ;

  /*multiply by \beta*/
  df[0] = bt[0]*f[0] - bt[1]*f[1] ;
  f [1] = bt[1]*f[0] + bt[0]*f[1] ;
  f [0] = df[0] ;
   
  for ( i = 0 ; i < ns ; i ++ ) {
    nbi_vector_diff(r, x, &(xs[i*xstr])) ;
    R = nbi_vector_length(r) ;
    if ( R > NBI_LOCAL_CUTOFF_RADIUS ) {
      C = cos(k*R) ; S = sin(k*R) ;
      /*negative because r = x - x_{1} and rn is dR/dn_{1}*/
      rn = -nbi_vector_scalar(r,&(xs[i*xstr+3])) ;
      G[0] = 0.25*M_1_PI*C/R ;
      G[1] = 0.25*M_1_PI*S/R ;
      dG[0] = -0.25*M_1_PI*(C + k*R*S)*rn/R/R/R ;
      dG[1] =  0.25*M_1_PI*(k*R*C - S)*rn/R/R/R ;
      df[0] = wt*(pwt*(dG[0]*p [i*pstr + 0] - dG[1]*p [i*pstr + 1]) -
		  nwt*( G[0]*pn[i*nstr + 0] -  G[1]*pn[i*nstr + 1])) ;
      df[1] = wt*(pwt*(dG[1]*p [i*pstr + 0] + dG[0]*p [i*pstr + 1]) -
		  nwt*( G[1]*pn[i*nstr + 0] +  G[0]*pn[i*nstr + 1])) ;
      f[0] += al[0]*df[0] - al[1]*df[1] ;
      f[1] += al[1]*df[0] + al[0]*df[1] ;
    }
  }
  
  return ;
}

static void upsample_sources(nbi_surface_t *s,
			     NBI_REAL *p, gint pstr, NBI_REAL *pn, gint nstr,
			     NBI_REAL *wt, gint wstr, gint *idxu,
			     NBI_REAL *pu, gint pustr, NBI_REAL pwt,
			     NBI_REAL *pnu, gint nustr, NBI_REAL nwt)

/*
 * interpolate from source nodes to upsampled nodes:
 *
 * pwt: weight for potential source (double layer potential)
 * nwt: weight for normal derivative source (single layer potential)
 */
  
{
  gint i, j, pt, ns, nu ;
  NBI_REAL *K, al, bt ;

  al = 1.0 ; bt = 0.0 ;
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    ns = nbi_surface_patch_node_number(s, pt) ;
    nu = idxu[pt+1] - idxu[pt] ;
    K = NBI_FUNCTION_NAME(nbi_patch_upsample_matrix)(ns, nu) ;

    i = nbi_surface_patch_node(s, pt) ;
    j = idxu[pt] ;
#ifdef NBI_SINGLE_PRECISION
    g_assert_not_reached() ; /*unchecked code*/
    blaswrap_sgemv(FALSE, nu, ns, al, K, ns,
		   &(p[i*pstr]), pstr, bt, &(pu[j*pustr]), pustr) ;
    blaswrap_sgemv(FALSE, nu, ns, al, K, ns,
		   &(pn[i*pstr]), nstr, bt, &(pnu[j*pustr]), nustr) ;
#else  /*NBI_SINGLE_PRECISION*/
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(p[i*pstr+0]), pstr, bt, &(pu[j*pustr+0]), pustr) ;
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(p[i*pstr+1]), pstr, bt, &(pu[j*pustr+1]), pustr) ;
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(pn[i*nstr+0]), nstr, bt, &(pnu[j*nustr+0]), nustr) ;
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(pn[i*nstr+1]), nstr, bt, &(pnu[j*nustr+1]), nustr) ;
#endif /*NBI_SINGLE_PRECISION*/
  }

  nu = idxu[nbi_surface_patch_number(s)] ;
  for ( i = 0 ; i < nu ; i ++ ) {
    pnu[i*nustr+0] *= nwt*wt[i*wstr] ;
    pnu[i*nustr+1] *= nwt*wt[i*wstr] ;
    pu [i*pustr+0] *= pwt*wt[i*wstr] ;
    pu [i*pustr+1] *= pwt*wt[i*wstr] ;
  }
  
  return ;
}

/* static void upsample_source_write(FILE *f, */
/* 				  nbi_surface_t *s, gint *idxu, */
/* 				  NBI_REAL *xu, gint ustr, */
/* 				  NBI_REAL *pu, gint pustr, */
/* 				  NBI_REAL *pnu, gint nustr) */

/* { */
/*   gint i, j, pt, ns, nu ; */
/*   NBI_REAL *x ; */

/*   for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) { */
/*     ns = nbi_surface_patch_node_number(s, pt) ; */
/*     nu = idxu[pt+1] - idxu[pt] ; */

/*     j = idxu[pt] ; */
/*     for ( i = 0 ; i < nu ; i ++ ) { */
/*       x = &(xu[(j+i)*ustr]) ; */
/*       fprintf(f, "%e %e %e %e %e %e %e %e %e %e %e\n", */
/* 	      x[0], x[1], x[2], x[3], x[4], x[5], x[6], */
/* 	      pu [(j+i)*pustr+0], pu[(j+i)*pustr+1], */
/* 	      pnu[(j+i)*nustr+0], pnu[(j+i)*nustr+1]) ; */
/*     } */
/*   } */
  
/*   return ; */
/* } */

gint NBI_FUNCTION_NAME(nbi_matrix_upsample_helmholtz)(nbi_matrix_t *m,
						      NBI_REAL *p, gint pstr,
						      NBI_REAL pwt,
						      NBI_REAL *pn, gint nstr,
						      NBI_REAL nwt)

{
  NBI_REAL *xu, *su, *sn ;

  su = (NBI_REAL *)(m->p ) ;
  sn = (NBI_REAL *)(m->pn) ;
  xu = (NBI_REAL *)(m->xu) ;
  upsample_sources(m->s, p, pstr, pn, nstr, 
		   &(xu[6]), m->ustr, m->idxu,
		   su, m->pstr, pwt, sn, m->nstr, nwt) ;
  
  return 0 ;
}

static void point_source_summation(nbi_matrix_t *m,
				   NBI_REAL *f, gint fstr,
				   NBI_REAL *al, NBI_REAL *bt,
				   NBI_REAL *work, gint nthreads)
/*
 * f := bt*f + al*SUM(dG\phi - Gd\phi)
 */
  
{
  gint i, nu, np, ustr, *idxu, pustr, nustr ;
  nbi_surface_t *s ;
  NBI_REAL *xu, *pu, *pnu, k ;

  k = nbi_matrix_wavenumber(m) ;
  ustr = m->ustr ;
  idxu = m->idxu ; pustr = m->pstr ; nustr = m->nstr ;
  xu = (NBI_REAL *)(m->xu) ;
  pu  = (NBI_REAL *)(m->p) ;
  pnu = (NBI_REAL *)(m->pn) ;
  s = m->s ;

  np = nbi_surface_patch_number(s) ; nu = idxu[np] ;

  if ( m->tree == NULL ) {
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      point_source_field_helmholtz(k,
				   xu, ustr, nu,
				   pu, pustr, 1.0,
				   pnu, nustr, 1.0,
				   (NBI_REAL *)nbi_surface_node(s, i),
				   1.0, &(f[i*fstr]), al, bt) ;
    }
  } else {
    /*FMM summation*/
    gint level, depth ;
    NBI_REAL kc[] = {0, -k} ;
    wbfmm_tree_t *tree ;
    wbfmm_target_list_t *targets ;
    wbfmm_shift_operators_t *shifts ;
    
    tree = m->tree ; targets = m->targets ; shifts = m->shifts ;

    depth = wbfmm_tree_depth(tree) ;
    for ( level = 2 ; level <= depth ; level ++ ) {
      wbfmm_tree_coefficients_zero(tree, level) ;
    }
    wbfmm_tree_leaf_expansions(tree, k,
			       pnu, nustr,
			       &(xu[3]), ustr,
			       pu, pustr,
			       TRUE, work) ;
    for ( level = depth ; level >= 3 ; level -- ) {
      wbfmm_upward_pass(tree, shifts, level, work) ;
    }
    for ( level = 2 ; level <= depth ; level ++ ) {      
      wbfmm_downward_pass(tree, shifts, level, work, nthreads) ;
    }
    if ( targets != NULL ) {
      g_assert_not_reached() ;
      wbfmm_target_list_local_field(targets, pnu, nustr, pu, pustr, f, fstr) ;
    } else {
      for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
	guint64 box ;
	NBI_REAL df[2] ;
	box = wbfmm_point_box(tree, tree->depth,
			      (NBI_REAL *)nbi_surface_node(s, i)) ;
	df[0] = bt[0]*f[i*fstr+0] - bt[1]*f[i*fstr+1] ;
	f [i*fstr+1] = bt[1]*f[i*fstr+0] + bt[0]*f[i*fstr+1] ;
	f [i*fstr+0] = df[0] ;
	df[0] = df[1] = 0.0 ;
	wbfmm_tree_box_local_field(tree, tree->depth, box,
				   k,
				   (NBI_REAL *)nbi_surface_node(s,i),
				   df, 2,
				   pnu, nustr,
				   &(xu[3]), ustr,
				   pu, pustr,
				   TRUE, WBFMM_FIELD_SCALAR, work) ;
	/*wbfmm calculates h_{0}(kR)/4\pi rather than the Green's
	  function, so we need to multiply by -j k*/
	nbi_scale_complex(df, kc) ;
	f[i*fstr+0] += al[0]*df[0] - al[1]*df[1] ;
	f[i*fstr+1] += al[1]*df[0] + al[0]*df[1] ;
      }
    }
  }
  
  return ;
}

static void local_matrix_correction(nbi_matrix_t *m,
				    NBI_REAL *p, gint pstr, NBI_REAL pwt,
				    NBI_REAL *pn, gint nstr, NBI_REAL nwt,
				    gint pt0, gint pt1,
				    NBI_REAL *al, NBI_REAL *bt,
				    NBI_REAL *buf, NBI_REAL *work)

{
  gint nsts, ip, *nbrs, nnbrs, lda, i, one = 1, pt, *idx, *idxp, str ;
  NBI_REAL *Ast, Al[2], Bt[] = {0, 0}, *A ;
  nbi_surface_t *s = m->s ;

  A = (NBI_REAL *)(m->Ast) ;
  idx = m->idx ; idxp = m->idxp ;

  for ( pt = pt0 ; pt < pt1 ; pt ++ ) {
    nsts = nbi_surface_patch_node_number(s, pt) ;

    ip = nbi_surface_patch_node(s, pt) ;
    nnbrs = idxp[pt+1] - idxp[pt] ;
    nbrs = &(idx[idxp[pt]]) ;
    Ast = &(A[4*nsts*idxp[pt]]) ;

    lda = 2*nsts ;
#ifdef NBI_SINGLE_PRECISION
    g_assert_not_reached() ; /*unchecked code*/
    /* blaswrap_sgemv(FALSE, nnbrs, nsts, al, &(Ast[1*nsts]), lda, */
    /* 		    &(p[ip*pstr]) , pstr, bt, work, one) ; */
    /* al = -nwt*sgn ; bt = 1.0 ; */
    /* blaswrap_sgemv(FALSE, nnbrs, nsts, al, &(Ast[0*nsts]), lda, */
    /* 		   &(pn[ip*nstr]), nstr, bt, work, one) ; */
#else  /*NBI_SINGLE_PRECISION*/
    str = pstr/2 ;
    Al[0] =  pwt*al[0] ; Al[1] = pwt*al[1] ; Bt[0] = 0 ; Bt[1] = 0 ;
    blaswrap_zgemv(FALSE, nnbrs, nsts, Al, &(Ast[2*nsts]), lda,
		   &(p[ip*pstr]), str, Bt, work, one) ;
    Al[0] =  -nwt*al[0] ; Al[1] = -nwt*al[1] ; Bt[0] = 1 ; Bt[1] = 0 ;
    str = nstr/2 ;
    blaswrap_zgemv(FALSE, nnbrs, nsts, Al, &(Ast[0*nsts]), lda,
		   &(pn[ip*nstr]), str, Bt, work, one) ;
#endif /*NBI_SINGLE_PRECISION*/    
    for ( i = 0 ; i < nnbrs ; i ++ ) {
      buf[2*nbrs[i]+0] += work[2*i+0] ;
      buf[2*nbrs[i]+1] += work[2*i+1] ;
    }
  }

  return ;
}

static gpointer local_correction_thread(gpointer tdata)

{
  gpointer *mdata = tdata ;
  gint th  = GPOINTER_TO_INT(mdata[NBI_THREAD_MAIN_DATA_THREAD]) ;
  gint nth = GPOINTER_TO_INT(mdata[NBI_THREAD_MAIN_DATA_NTHREAD]) ;
  gpointer *data = mdata[NBI_THREAD_MAIN_DATA_DATA] ;
  NBI_REAL *work = mdata[NBI_THREAD_MAIN_DATA_WORK] ;
  NBI_REAL *buf = mdata[NBI_THREAD_MAIN_DATA_BUFFER] ;
  nbi_matrix_t *m = data[NBI_THREAD_DATA_MATRIX  ] ;
  gint *idata = data[NBI_THREAD_DATA_INT] ;
  NBI_REAL *ddata = data[NBI_THREAD_DATA_REAL] ;
  NBI_REAL **dpdata = data[NBI_THREAD_DATA_REAL_POINTER] ;
  gint np, pt0, pt1, pstr, nstr ;
  NBI_REAL *p, *pn, pwt, nwt, *al, *bt ;

  np = nbi_surface_patch_number(m->s) ;

  pt0 = th*(np/nth) ;
  pt1 = (th+1)*(np/nth) ;

  if ( th == nth - 1 ) {
    if ( pt1 < np ) pt1 = np ;
  }

  p = dpdata[NBI_THREAD_DATA_REAL_PTR_P ]  ;
  pn = dpdata[NBI_THREAD_DATA_REAL_PTR_PN] ;
  /* f = dpdata[NBI_THREAD_DATA_REAL_PTR_F ]  ; */

  pstr = idata[NBI_THREAD_DATA_INT_PSTR] ;
  nstr = idata[NBI_THREAD_DATA_INT_NSTR] ;
  /* fstr = idata[NBI_THREAD_DATA_INT_FSTR] ; */
  
  pwt = ddata[NBI_THREAD_DATA_REAL_WT1 ] ;
  nwt = ddata[NBI_THREAD_DATA_REAL_WT2 ] ;

  al = dpdata[NBI_THREAD_DATA_REAL_PTR_WT1] ;
  bt = dpdata[NBI_THREAD_DATA_REAL_PTR_WT2] ;

  local_matrix_correction(m, p, pstr, pwt, pn, nstr, nwt,
			  pt0, pt1, al, bt, buf, work) ;  
  
  return NULL ;
}

gint NBI_FUNCTION_NAME(nbi_surface_greens_identity_helmholtz)(nbi_matrix_t *m,
							      NBI_REAL *p ,
							      gint pstr,
							      NBI_REAL pwt,
							      NBI_REAL *pn,
							      gint nstr,
							      NBI_REAL nwt,
							      NBI_REAL *f,
							      gint fstr,
							      gint nthreads,
							      NBI_REAL *work)

{
  gint i, nnmax, nth, npts, i2 = 2 ;
  NBI_REAL al[] = {1, 0}, bt[] = {0, 0}, *buffer, d1 = 1 ;

  /*f is complex*/
  g_assert(fstr >= 2) ;
  
  if ( nthreads < 0 ) nth = g_get_num_processors() ; else nth = nthreads ;
  
  NBI_FUNCTION_NAME(nbi_matrix_upsample_helmholtz)(m,
						   p , pstr, pwt,
						   pn, nstr, nwt) ;

  /*point source approximation (FMM handled internally)*/
  point_source_summation(m, f, fstr, al, bt, work, nth) ;
  
  nnmax = nbi_matrix_neighbour_number_max(m) ;
  npts = nbi_surface_node_number(m->s) ;

  /*local corrections*/
#ifdef _OPENMP
  if ( nth > 0 ) {
    gint idata[NBI_THREAD_DATA_INT_SIZE] ;
    NBI_REAL
      ddata[NBI_THREAD_DATA_REAL_SIZE],
      *dpdata[NBI_THREAD_DATA_REAL_PTR_SIZE] ;
    GThread *threads[NBI_THREAD_NUMBER_MAX] ;
    gpointer data[NBI_THREAD_DATA_SIZE],
      main_data[NBI_THREAD_NUMBER_MAX*NBI_THREAD_MAIN_DATA_SIZE] ;

    buffer = &(work[nth*2*nnmax]) ;
    memset(buffer, 0, nth*2*npts*sizeof(NBI_REAL)) ;
    
    dpdata[NBI_THREAD_DATA_REAL_PTR_P ] = p  ;
    dpdata[NBI_THREAD_DATA_REAL_PTR_PN] = pn ;
    dpdata[NBI_THREAD_DATA_REAL_PTR_F ] = f  ;
    dpdata[NBI_THREAD_DATA_REAL_PTR_WT1] = al ;
    dpdata[NBI_THREAD_DATA_REAL_PTR_WT2] = bt ;

    idata[NBI_THREAD_DATA_INT_PSTR] = pstr ;
    idata[NBI_THREAD_DATA_INT_NSTR] = nstr ;
    idata[NBI_THREAD_DATA_INT_FSTR] = fstr ;
  
    ddata[NBI_THREAD_DATA_REAL_WT1 ] = pwt ;
    ddata[NBI_THREAD_DATA_REAL_WT2 ] = nwt ;

    data[NBI_THREAD_DATA_MATRIX] = m ;
    data[NBI_THREAD_DATA_INT] = idata ;
    data[NBI_THREAD_DATA_REAL] = ddata ;
    data[NBI_THREAD_DATA_REAL_POINTER] = dpdata ;
    
    for ( i = 0 ; i < nth ; i ++ ) {
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_THREAD] =
	GINT_TO_POINTER(i) ;
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_DATA] =
	data ;
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_WORK] =
	&(work[i*2*nnmax]) ;
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_NTHREAD] =
	GINT_TO_POINTER(nth) ;
      
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_BUFFER] =
	&(buffer[i*2*npts]) ;

      threads[i] = g_thread_new(NULL, local_correction_thread,
				&(main_data[NBI_THREAD_MAIN_DATA_SIZE*i])) ;
    }

    /*make sure all threads complete before we move on*/
    for ( i = 0 ; i < nth ; i ++ ) g_thread_join(threads[i]) ;
  } else {
    buffer = &(work[2*nnmax]) ;
    memset(buffer, 0, 2*npts*sizeof(NBI_REAL)) ;
    local_matrix_correction(m, p, pstr, pwt, pn, nstr, nwt,
			    0, nbi_surface_patch_number(m->s),
			    al, bt, buffer, work) ;
  }
  
#else /*_OPENMP*/
  buffer = &(work[2*nnmax]) ;
  memset(buffer, 0, 2*npts*sizeof(NBI_REAL)) ;
  local_matrix_correction(m, p, pstr, pwt, pn, nstr, nwt,
			  0, nbi_surface_patch_number(m->s),
			  al, bt, buffer, work) ;
#endif /*_OPENMP*/
  i = 0 ;
  blaswrap_daxpy(npts, d1, &(buffer[i*2*npts+0]), i2, &(f[0]), fstr) ;
  blaswrap_daxpy(npts, d1, &(buffer[i*2*npts+1]), i2, &(f[1]), fstr) ;
  for ( i = 1 ; i < nth ; i ++ ) {
    blaswrap_daxpy(npts, d1, &(buffer[i*2*npts+0]), i2, &(f[0]), fstr) ;
    blaswrap_daxpy(npts, d1, &(buffer[i*2*npts+1]), i2, &(f[1]), fstr) ;
  }
  
  return 0 ;
}

static void upsample_sources_single(nbi_surface_t *s,
				    NBI_REAL *pn, gint nstr,
				    NBI_REAL *wt, gint wstr, gint *idxu,
				    NBI_REAL *pu, gint pustr, NBI_REAL pwt,
				    NBI_REAL *pnu, gint nustr, NBI_REAL nwt)

/*
 * upsample sources for single layer potential (normal derivative source)
 */
  
{
  gint i, j, pt, ns, nu ;
  NBI_REAL *K, al, bt ;

  al =  1.0 ; bt = 0.0 ;
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    ns = nbi_surface_patch_node_number(s, pt) ;
    nu = idxu[pt+1] - idxu[pt] ;
    K = NBI_FUNCTION_NAME(nbi_patch_upsample_matrix)(ns, nu) ;

    i = nbi_surface_patch_node(s, pt) ;
    j = idxu[pt] ;
#ifdef NBI_SINGLE_PRECISION
    g_assert_not_reached() ; /*unchecked code*/
    blaswrap_sgemv(FALSE, nu, ns, al, K, ns,
		   &(pn[i*pstr]), nstr, bt, &(pnu[j*pustr]), nustr) ;
#else  /*NBI_SINGLE_PRECISION*/
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(pn[i*nstr+0]), nstr, bt, &(pnu[j*nustr+0]), nustr) ;
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(pn[i*nstr+1]), nstr, bt, &(pnu[j*nustr+1]), nustr) ;
#endif /*NBI_SINGLE_PRECISION*/
  }

  nu = idxu[nbi_surface_patch_number(s)] ;
  for ( i = 0 ; i < nu ; i ++ ) {
    pu [i*pustr+0] = pu [i*pustr+1] = 0.0 ;
    pnu[i*nustr+0] *= nwt*wt[i*wstr] ;
    pnu[i*nustr+1] *= nwt*wt[i*wstr] ;
  }
  
  return ;
}

static void upsample_sources_double(nbi_surface_t *s,
				    NBI_REAL *p, gint pstr,
				    NBI_REAL *wt, gint wstr, gint *idxu,
				    NBI_REAL *pu, gint pustr, NBI_REAL pwt,
				    NBI_REAL *pnu, gint nustr, NBI_REAL nwt)

/*
 * upsample sources for double layer potential (potential source)
 */

{
  gint i, j, pt, ns, nu ;
  NBI_REAL *K, al, bt ;

  al =  1.0 ; bt = 0.0 ;
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    ns = nbi_surface_patch_node_number(s, pt) ;
    nu = idxu[pt+1] - idxu[pt] ;
    K = NBI_FUNCTION_NAME(nbi_patch_upsample_matrix)(ns, nu) ;

    i = nbi_surface_patch_node(s, pt) ;
    j = idxu[pt] ;
#ifdef NBI_SINGLE_PRECISION
    g_assert_not_reached() ; /*unchecked code*/
    blaswrap_sgemv(FALSE, nu, ns, al, K, ns,
    		   &(x[i*xstr]), xstr, bt, &(pu[j*pustr]), pustr) ;
#else  /*NBI_SINGLE_PRECISION*/
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(p[i*pstr+0]), pstr, bt, &(pu[j*pustr+0]), pustr) ;
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(p[i*pstr+1]), pstr, bt, &(pu[j*pustr+1]), pustr) ;
#endif /*NBI_SINGLE_PRECISION*/
  }

  nu = idxu[nbi_surface_patch_number(s)] ;
  for ( i = 0 ; i < nu ; i ++ ) {
    pnu[i*nustr+0] = pnu[i*nustr+1] = 0.0 ;
    pu [i*pustr+0] *= pwt*wt[i*wstr] ;
    pu [i*pustr+1] *= pwt*wt[i*wstr] ;
  }
  
  return ;
}

static void local_matrix_multiply(nbi_matrix_t *m, NBI_REAL *p, gint pstr,
				  gint off, NBI_REAL *al,
				  NBI_REAL *buf, NBI_REAL *work,
				  /* NBI_REAL *f, gint fstr, */
				  gint pt0, gint pt1,
				  gint nthreads)

{
  gint i, ip, nnbrs, *nbrs, one = 1, nsts, lda, pt, str ;
  NBI_REAL *Ast, *A, c0[] = {0, 0} ;
  
  A = (NBI_REAL *)(m->Ast) ;
  A = &(A[off]) ;
  /*loop on patches treated as sources*/
  nsts = nbi_surface_patch_node_number(m->s, 0) ;
  lda = 2*nsts ;
  
  for ( pt = pt0 ; pt < pt1 ; pt ++ ) {
    ip = nbi_surface_patch_node(m->s, pt) ;
    nnbrs = m->idxp[pt+1] - m->idxp[pt] ;
    g_assert(nnbrs != 0) ;
    nbrs = &(m->idx[m->idxp[pt]]) ;
    Ast = &(A[4*nsts*(m->idxp[pt])]) ;
    
#ifdef NBI_SINGLE_PRECISION
    g_assert_not_reached() ; /*unchecked code*/
    blaswrap_sgemv(FALSE, nnbrs, nsts, al, Ast, lda,
		   &(p[ip*pstr]) , pstr, c0, work, one) ;
#else  /*NBI_SINGLE_PRECISION*/
    str = pstr/2 ;
    blaswrap_zgemv(FALSE, nnbrs, nsts, al, Ast, lda,
		   &(p[ip*pstr]), str, c0, work, one) ;
#endif /*NBI_SINGLE_PRECISION*/    
    
    /* for ( i = 0 ; i < nnbrs ; i ++ ) { */
    /*   f[fstr*nbrs[i]+0] += work[2*i+0] ; */
    /*   f[fstr*nbrs[i]+1] += work[2*i+1] ; */
    /* } */
    for ( i = 0 ; i < nnbrs ; i ++ ) {
      buf[2*nbrs[i]+0] += work[2*i+0] ;
      buf[2*nbrs[i]+1] += work[2*i+1] ;
    }
  }

  return ;
}

#ifdef _OPENMP

static gpointer matrix_multiply_thread(gpointer tdata)

{
  gpointer *mdata = tdata ;
  gint th  = GPOINTER_TO_INT(mdata[NBI_THREAD_MAIN_DATA_THREAD]) ;
  gint nth = GPOINTER_TO_INT(mdata[NBI_THREAD_MAIN_DATA_NTHREAD]) ;
  gpointer *data = mdata[NBI_THREAD_MAIN_DATA_DATA] ;
  NBI_REAL *work = mdata[NBI_THREAD_MAIN_DATA_WORK] ;
  NBI_REAL *buf = mdata[NBI_THREAD_MAIN_DATA_BUFFER] ;
  nbi_matrix_t *m = data[NBI_THREAD_DATA_MATRIX] ;
  gint *idata = data[NBI_THREAD_DATA_INT] ;
  NBI_REAL **dpdata = data[NBI_THREAD_DATA_REAL_POINTER] ;
  NBI_REAL *p, *al ;
  gint pt0, pt1, np, pstr, off ;

  np = nbi_surface_patch_number(m->s) ;
  pt0 = th*(np/nth) ;
  pt1 = (th+1)*(np/nth) ;

  if ( th == nth - 1 ) { pt1 = np ; }

  p = dpdata[NBI_THREAD_DATA_REAL_PTR_P] ;
  al = dpdata[NBI_THREAD_DATA_REAL_PTR_WT1] ;

  pstr = idata[NBI_THREAD_DATA_INT_PSTR] ;
  off  = idata[NBI_THREAD_DATA_INT_OFFSET] ;
  
  local_matrix_multiply(m, p, pstr, off, al, buf, work, pt0, pt1, nth) ;

  return NULL ;
}

static void local_matrix_multiply_thread(nbi_matrix_t *m,
					 NBI_REAL *p, gint pstr,
					 gint off, NBI_REAL *al, NBI_REAL *bt,
					 NBI_REAL *work, gint nthreads)
{
  gint nth, i, nnmax, idata[NBI_THREAD_DATA_INT_SIZE], npts ;
  NBI_REAL
    *buffer,
    ddata[NBI_THREAD_DATA_REAL_SIZE],
    *dpdata[NBI_THREAD_DATA_REAL_PTR_SIZE] ;
  GThread *threads[NBI_THREAD_NUMBER_MAX] ;
  gpointer data[NBI_THREAD_DATA_SIZE],
    main_data[NBI_THREAD_NUMBER_MAX*NBI_THREAD_MAIN_DATA_SIZE] ;

  nnmax = nbi_matrix_neighbour_number_max(m) ;
  npts = nbi_surface_node_number(m->s) ;
  
  if ( nthreads < 0 ) nth = g_get_num_processors() ; else nth = nthreads ;

  dpdata[NBI_THREAD_DATA_REAL_PTR_P] = p ;
  dpdata[NBI_THREAD_DATA_REAL_PTR_WT1] = al ;
  dpdata[NBI_THREAD_DATA_REAL_PTR_WT2] = bt ;

  idata[NBI_THREAD_DATA_INT_PSTR]   = pstr ;
  idata[NBI_THREAD_DATA_INT_OFFSET] = off ;

  data[NBI_THREAD_DATA_MATRIX] = m ;
  data[NBI_THREAD_DATA_INT] = idata ;
  data[NBI_THREAD_DATA_REAL] = ddata ;
  data[NBI_THREAD_DATA_REAL_POINTER] = dpdata ;

  buffer = &(work[nth*2*nnmax]) ;
  for ( i = 0 ; i < nth ; i ++ ) {
    main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_THREAD] =
      GINT_TO_POINTER(i) ;
    main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_DATA] =
      data ;
    main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_WORK] =
      &(work[i*2*nnmax]) ;
    main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_NTHREAD] =
      GINT_TO_POINTER(nth) ;
    main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_BUFFER] =
      &(buffer[i*2*npts]) ;
  }

  /* memset(buffer, 0, nth*2*npts*sizeof(NBI_REAL)) ; */
  
  for ( i = 0 ; i < nth ; i ++ ) {
    threads[i] = g_thread_new(NULL, matrix_multiply_thread,
			      &(main_data[NBI_THREAD_MAIN_DATA_SIZE*i])) ;    
  }

  for ( i = 0 ; i < nth ; i ++ ) g_thread_join(threads[i]) ;

  /* for ( i = 0 ; i < nth ; i ++ ) { */
  /*   blaswrap_daxpy(npts, d1, &(buffer[i*2*npts+0]), i2, &(f[0]), fstr) ; */
  /*   blaswrap_daxpy(npts, d1, &(buffer[i*2*npts+1]), i2, &(f[1]), fstr) ; */
  /* } */
  
  return ;
}
#endif /*_OPENMP*/

static void matrix_multiply_single(nbi_matrix_t *m,
				   NBI_REAL *pn, gint nstr, NBI_REAL *Al,
				   NBI_REAL *f, gint fstr,
				   gint nthreads,
				   NBI_REAL *work) 

{
  gint nth, nnmax, npts, ustr, *idxu, pustr, nustr, i2 = 2, i ;
  NBI_REAL al[2], bt[2], *xu, *pu, *pnu, *buf, d1 = 1 ;
  nbi_surface_t *s ;

  ustr = m->ustr ;
  idxu = m->idxu ;
  pustr = m->pstr ; nustr = m->nstr ;
  xu = (NBI_REAL *)(m->xu) ;
  pu  = (NBI_REAL *)(m->p) ;
  pnu = (NBI_REAL *)(m->pn) ;
  s = m->s ;
  nnmax = nbi_matrix_neighbour_number_max(m) ;
  npts = nbi_surface_node_number(m->s) ;

  if ( nthreads < 0 ) nth = g_get_num_processors() ; else nth = nthreads ;

  upsample_sources_single(s, pn, nstr, &(xu[6]), ustr, idxu,
			  pu, pustr, 0.0, pnu, nustr, 1.0) ;

  al[0] = -Al[0] ; al[1] = -Al[1] ;
  bt[0] = 1 ; bt[1] = 0 ;
  point_source_summation(m, f, fstr, al, bt, work, nth) ;

  /*local corrections*/
  al[0] = Al[0] ; al[1] = Al[1] ;
  bt[0] = 0.0 ; bt[1] = 0.0 ;
#ifdef _OPENMP
  if ( nth > 0 ) {
    buf = &(work[2*nth*nnmax]) ;
    memset(buf, 0, 2*nth*npts*sizeof(NBI_REAL)) ;
    local_matrix_multiply_thread(m, pn, nstr, 0, al, bt, work, nth) ;
  } else {
    buf = &(work[2*nnmax]) ;
    memset(buf, 0, 2*npts*sizeof(NBI_REAL)) ;
    local_matrix_multiply(m, pn, nstr, 0, al, buf, work,
			  0, nbi_surface_patch_number(s), nth) ;
  }
#else /*_OPENMP*/
  g_assert_not_reached() ; /*untested code*/
  buf = &(work[2*nnmax]) ;
  memset(buf, 0, 2*npts*sizeof(NBI_REAL)) ;
  local_matrix_multiply(m, pn, nstr, 0, al, buf, work,
			0, nbi_surface_patch_number(s), nth) ;
#endif /*_OPENMP*/
  
  i = 0 ;
  blaswrap_daxpy(npts, d1, &(buf[i*2*npts+0]), i2, &(f[0]), fstr) ;
  blaswrap_daxpy(npts, d1, &(buf[i*2*npts+1]), i2, &(f[1]), fstr) ;
  for ( i = 1 ; i < nth ; i ++ ) {
    blaswrap_daxpy(npts, d1, &(buf[i*2*npts+0]), i2, &(f[0]), fstr) ;
    blaswrap_daxpy(npts, d1, &(buf[i*2*npts+1]), i2, &(f[1]), fstr) ;
  }

  return ;
}

static void matrix_multiply_double(nbi_matrix_t *m,
				   NBI_REAL *p, gint pstr, NBI_REAL *al,
				   NBI_REAL *f, gint fstr,
				   gint nthreads,
				   NBI_REAL *work) 

{
  gint nsts, nth, nnmax, npts, ustr, *idxu, pustr, nustr, i2 = 2, i ;
  NBI_REAL bt[2], *xu, *pu, *pnu, *buf, d1 = 1 ;
  nbi_surface_t *s ;

  ustr = m->ustr ;
  idxu = m->idxu ; pustr = m->pstr ; nustr = m->nstr ;
  xu = (NBI_REAL *)(m->xu) ;
  pu  = (NBI_REAL *)(m->p) ;
  pnu = (NBI_REAL *)(m->pn) ;
  s = m->s ;
  nnmax = nbi_matrix_neighbour_number_max(m) ;
  npts = nbi_surface_node_number(m->s) ;

  if ( nthreads < 0 ) nth = g_get_num_processors() ; else nth = nthreads ;

  upsample_sources_double(s, p, pstr, &(xu[6]), ustr, idxu,
			  pu , pustr, 1.0, pnu, nustr, 0.0) ;

  /*point source approximation (FMM handled internally)*/
  bt[0] = 1 ; bt[1] = 0 ;
  point_source_summation(m, f, fstr, al, bt, work, nth) ;
  
  /*local corrections*/
  nsts = nbi_surface_patch_node_number(s, 0) ;
  bt[0] = 0.0 ; bt[1] = 0.0 ;
  
#ifdef _OPENMP
  if ( nth > 0 ) {
    buf = &(work[2*nth*nnmax]) ;
    memset(buf, 0, 2*nth*npts*sizeof(NBI_REAL)) ;
    local_matrix_multiply_thread(m, p, pstr, 2*nsts, al, bt, work,
				 nth) ;
  } else {
    buf = &(work[2*nnmax]) ;
    memset(buf, 0, 2*npts*sizeof(NBI_REAL)) ;
    local_matrix_multiply(m, p, pstr, 2*nsts, al, buf, work,
			  0, nbi_surface_patch_number(s), nth) ;
  }
#else /*_OPENMP*/
  g_assert_not_reached() ; /*untested code*/
  buf = &(work[2*nnmax]) ;
  memset(buf, 0, 2*npts*sizeof(NBI_REAL)) ;
  local_matrix_multiply(m, p, pstr, 2*nsts, al, buf, work,
			0, nbi_surface_patch_number(s), nth) ;
#endif /*_OPENMP*/
  
  i = 0 ;
  blaswrap_daxpy(npts, d1, &(buf[i*2*npts+0]), i2, &(f[0]), fstr) ;
  blaswrap_daxpy(npts, d1, &(buf[i*2*npts+1]), i2, &(f[1]), fstr) ;
  for ( i = 1 ; i < nth ; i ++ ) {
    blaswrap_daxpy(npts, d1, &(buf[i*2*npts+0]), i2, &(f[0]), fstr) ;
    blaswrap_daxpy(npts, d1, &(buf[i*2*npts+1]), i2, &(f[1]), fstr) ;
  }

  return ;
}

gint NBI_FUNCTION_NAME(nbi_matrix_multiply_helmholtz)(nbi_matrix_t *A,
						      NBI_REAL *x, gint xstr,
						      NBI_REAL al,
						      NBI_REAL *y, gint ystr,
						      NBI_REAL bt,
						      gint nthreads,
						      NBI_REAL *work)

/*
  y := al*A*x + bt*y
*/
  
{
  NBI_REAL diag, Al[2] ;
  gint n ;

  g_assert(A->problem == NBI_PROBLEM_HELMHOLTZ) ;

  if ( xstr < 2 )
    g_error("%s: xstr (%d) too small for complex variables",
	    __FUNCTION__, xstr) ;
  if ( ystr < 2 )
    g_error("%s: ystr (%d) too small for complex variables",
	    __FUNCTION__, ystr) ;
  
  diag = ((NBI_REAL)(A->diag))*al ;
  
#ifdef NBI_SINGLE_PRECISION
  g_assert_not_reached() ; /*untested code*/
  blaswrap_sscal(nbi_surface_node_number(A->s), bt, y, ystr) ;
  blaswrap_saxpy(nbi_surface_node_number(A->s), diag, x, xstr, y, ystr) ;
#else /*NBI_SINGLE_PRECISION*/
  n = nbi_surface_node_number(A->s) ;
  blaswrap_dscal(n, bt, &(y[0]), ystr) ;
  blaswrap_dscal(n, bt, &(y[1]), ystr) ;
  blaswrap_daxpy(n, diag, &(x[0]), xstr, &(y[0]), ystr) ;
  blaswrap_daxpy(n, diag, &(x[1]), xstr, &(y[1]), ystr) ;
#endif /*NBI_SINGLE_PRECISION*/

  Al[0] = al ; Al[1] = 0.0 ;

  switch ( A->potential ) {
  default:
  case NBI_POTENTIAL_UNDEFINED:
    g_error("%s: required potential must be set for multiplication",
  	    __FUNCTION__) ;
    break ;
  case NBI_POTENTIAL_SINGLE:
    matrix_multiply_single(A, x, xstr, Al, y, ystr, nthreads, work) ;
    return 0 ;
    break ;
  case NBI_POTENTIAL_DOUBLE:
    matrix_multiply_double(A, x, xstr, Al, y, ystr, nthreads, work) ;
    return 0 ;
    break ;
  }
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_calc_field_helmholtz)(nbi_matrix_t *m,
						 NBI_REAL *p ,
						 gint pstr,
						 NBI_REAL pwt,
						 NBI_REAL *pn,
						 gint nstr,
						 NBI_REAL nwt,
						 NBI_REAL *x,
						 NBI_REAL *f,
						 gint nthreads,
						 NBI_REAL *work)
  
{
  gint nu, ustr, *idxu, pustr, nustr ;
  nbi_surface_t *s ;
  NBI_REAL *xu, *pu, *pnu, k, al[2] = {1, 0}, bt[] = {0, 0} ;
  
  NBI_FUNCTION_NAME(nbi_matrix_upsample_helmholtz)(m,
						   p , pstr, pwt,
						   pn, nstr, nwt) ;

  s = m->s ;
  k = nbi_matrix_wavenumber(m) ;
  ustr = m->ustr ;
  idxu = m->idxu ;
  pustr = m->pstr ; nustr = m->nstr ;
  xu = (NBI_REAL *)(m->xu) ;
  pu  = (NBI_REAL *)(m->p) ;
  pnu = (NBI_REAL *)(m->pn) ;
  nu = idxu[nbi_surface_patch_number(s)] ;

  point_source_field_helmholtz(k,
			       xu, ustr, nu,
			       pu, pustr, 1.0,
			       pnu, nustr, 1.0,
			       x, 1.0, f, al, bt) ;
  
  return 0 ;
}

