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

#include <sqt.h>

#include <wbfmm.h>

#include <blaswrap.h>

#include <nbi.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#ifdef __GNUC__
#if __GNUC__ == 11
#warning "duff optimizer"
#pragma GCC optimize("O0")
/* vect-cost-model=very-cheap") */
/* push_options ("string"...) */
#endif /*__GNUC__ == 11*/
#endif /*__GNUC__*/

#include "nbi-private.h"
#define wbfmm_tree_point_index(_t,_i)			\
  ((NBI_REAL *)(&((_t)->points[(_i)*((_t)->pstr)])))

#if 0
static gint point_source_field_laplace(NBI_REAL *xs, gint xstr, gint ns,
				       NBI_REAL *p , gint pstr, NBI_REAL pwt,
				       NBI_REAL *pn, gint nstr, NBI_REAL nwt,
				       NBI_REAL *x, NBI_REAL wt, NBI_REAL *f)
  
{
  gint i ;  
  NBI_REAL R, r[3] ;

  for ( i = 0 ; i < ns ; i ++ ) {
    nbi_vector_diff(r, x, &(xs[i*xstr])) ;
    R = nbi_vector_length(r) ;
    if ( R > NBI_LOCAL_CUTOFF_RADIUS ) {
      *f += wt*(pwt*p[i*pstr]*nbi_vector_scalar(r,&(xs[i*xstr+3]))/R/R -
		nwt*pn[i*nstr])*0.25*M_1_PI/R ;
    }
  }
  
  return 0 ;
}

static gint upsample_sources(nbi_surface_t *s,
			     NBI_REAL *p, gint pstr, NBI_REAL *pn, gint nstr,
			     NBI_REAL *wt, gint wstr, gint *idxu,
			     NBI_REAL *pu, gint pustr, NBI_REAL pwt,
			     NBI_REAL *pnu, gint nustr, NBI_REAL nwt)

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
		   &(p[i*pstr]), pstr, bt, &(pu[j*pustr]), pustr) ;
    blaswrap_sgemv(FALSE, nu, ns, al, K, ns,
		   &(pn[i*pstr]), nstr, bt, &(pnu[j*pustr]), nustr) ;
#else  /*NBI_SINGLE_PRECISION*/
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(p[i*pstr]), pstr, bt, &(pu[j*pustr]), pustr) ;
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(pn[i*pstr]), nstr, bt, &(pnu[j*pustr]), nustr) ;
#endif /*NBI_SINGLE_PRECISION*/
  }

  nu = idxu[nbi_surface_patch_number(s)] ;
  for ( i = 0 ; i < nu ; i ++ ) {
    pnu[i*nustr] *= nwt*wt[i*wstr] ;
    pu [i*pustr] *= pwt*wt[i*wstr] ;
  }
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_matrix_upsample_laplace)(nbi_matrix_t *m,
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

static void local_source_field(wbfmm_tree_t *t,
			       guint level,
			       guint64 b,
			       NBI_REAL *x,
			       NBI_REAL *f,
			       NBI_REAL *src, gint sstr,
			       NBI_REAL *normals, gint nstr,
			       NBI_REAL *d, gint dstr)

{
  NBI_REAL xb[3], wb, *C, *xs, r, rr[3], nr, g ;
  wbfmm_box_t *boxes, box ;
  guint64 neighbours[27] ;
  gint nnbr, i, k, idx, nq ;
  guint j ;
  
  nnbr = wbfmm_box_neighbours(level, b, neighbours) ;
  g_assert(nnbr >= 0 && nnbr < 28) ;

  nq = wbfmm_tree_source_size(t) ;

  boxes = t->boxes[level] ;
  for ( i = 0 ; i < nnbr ; i ++ ) {
    box = boxes[neighbours[i]] ;
    for ( j = 0 ; j < box.n ; j ++ ) {
      idx = t->ip[box.i+j] ;
      xs = wbfmm_tree_point_index(t, idx) ;
      rr[0] = x[0] - xs[0] ; 
      rr[1] = x[1] - xs[1] ; 
      rr[2] = x[2] - xs[2] ;
      r = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2] ;
      if ( r > WBFMM_LOCAL_CUTOFF_RADIUS*WBFMM_LOCAL_CUTOFF_RADIUS ) {
	r = SQRT(r) ;
	nr =
	  rr[0]*normals[idx*nstr+0] + 
	  rr[1]*normals[idx*nstr+1] + 
	  rr[2]*normals[idx*nstr+2] ;
	g = 0.25*M_1_PI/r ;
	for ( k = 0 ; k < nq ; k ++ ) {
	  f[k] += (d[idx*dstr+k]*nr/r/r + src[idx*sstr+k])*g ;
	}
      }
    }
    
  }   
}

static gint point_source_summation(nbi_matrix_t *m,
				   NBI_REAL *f, gint fstr,
				   NBI_REAL *work, gint nthreads)
{
  gint i, nu, np, ustr, *idxu, pustr, nustr ;
  nbi_surface_t *s ;
  NBI_REAL *xu, *pu, *pnu ;
  
  ustr = m->ustr ;
  idxu = m->idxu ; pustr = m->pstr ; nustr = m->nstr ;
  xu = (NBI_REAL *)(m->xu) ;
  pu  = (NBI_REAL *)(m->p) ;
  pnu = (NBI_REAL *)(m->pn) ;
  s = m->s ;

  np = nbi_surface_patch_number(s) ; nu = idxu[np] ;

  if ( m->tree == NULL ) {
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      point_source_field_laplace(xu, ustr, nu,
				 pu, pustr, 1.0,
				 pnu, nustr, 1.0,
				 (NBI_REAL *)nbi_surface_node(s, i),
				 1.0, &(f[i*fstr])) ;
    }
  } else {
    /*FMM summation*/
    gint level, depth ;
    wbfmm_tree_t *tree ;
    wbfmm_target_list_t *targets ;
    wbfmm_shift_operators_t *shifts ;
    
    tree = m->tree ; targets = m->targets ; shifts = m->shifts ;

    depth = wbfmm_tree_depth(tree) ;
    for ( level = 2 ; level <= depth ; level ++ ) {
      wbfmm_tree_coefficients_zero(tree, level) ;
    }
    wbfmm_tree_laplace_leaf_expansions(tree,
				       pnu, nustr,
				       &(xu[3]), ustr,
				       pu, pustr,
				       TRUE, work) ;
    /* fprintf(stderr, "%s: upward pass\n", __FUNCTION__) ; */
    for ( level = depth ; level >= 3 ; level -- ) {
      wbfmm_laplace_upward_pass(tree, shifts, level, work) ;
    }  
    /* fprintf(stderr, "%s: downward pass\n", __FUNCTION__) ; */
    for ( level = 2 ; level <= depth ; level ++ ) {      
      wbfmm_laplace_downward_pass(tree, shifts, level, work, nthreads) ;
    }
    if ( targets != NULL ) {
      wbfmm_target_list_local_field(targets, pnu, nustr, pu, pustr, f, fstr) ;
    } else {
      for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
	/* NBI_REAL fl[32] ; */
	guint64 box ;
	/* fprintf(stderr, "%s: surface node %d\n", __FUNCTION__, i) ; */
	box = wbfmm_point_box(tree, tree->depth,
			      (NBI_REAL *)nbi_surface_node(s, i)) ;
	wbfmm_tree_laplace_box_local_field(tree, tree->depth, box,
					   (NBI_REAL *)nbi_surface_node(s,i),
					   &(f[i*fstr]),
					   pnu, nustr,
					   &(xu[3]), ustr,
					   pu, pustr,
					   TRUE, work) ;
	/* local_source_field(tree, tree->depth, box, */
	/* 		   (NBI_REAL *)nbi_surface_node(s,i), */
	/* 		   &(f[i*fstr]), */
	/* 		   pnu, nustr, */
	/* 		   &(xu[3]), ustr, */
	/* 		   pu, pustr) ; */
      }
    }
  }

  return 0 ;
}

static gint local_matrix_correction(nbi_matrix_t *m,
				    gdouble *p, gint pstr, gdouble pwt,
				    gdouble *pn, gint nstr, gdouble nwt,
				    gdouble sgn,
				    gint pt0, gint pt1,
				    gdouble *f, gint fstr,
				    gdouble *work)

{
  gint nsts, ip, *nbrs, nnbrs, lda, i, one = 1, pt, *idx, *idxp ;
  gdouble *Ast, al, bt, *A ;
  nbi_surface_t *s = m->s ;
  
  A = (NBI_REAL *)(m->Ast) ;
  idx = m->idx ; idxp = m->idxp ;

  for ( pt = pt0 ; pt < pt1 ; pt ++ ) {
    nsts = nbi_surface_patch_node_number(s, pt) ;

    ip = nbi_surface_patch_node(s, pt) ;
    nnbrs = idxp[pt+1] - idxp[pt] ;
    nbrs = &(idx[idxp[pt]]) ;
    Ast = &(A[2*nsts*idxp[pt]]) ;

    lda = 2*nsts ;
    al =  pwt ; bt = 0.0 ;
#ifdef NBI_SINGLE_PRECISION
    g_assert_not_reached() ; /*unchecked code*/
    blaswrap_sgemv(FALSE, nnbrs, nsts, al, &(Ast[1*nsts]), lda,
		   &(p[ip*pstr]) , pstr, bt, work, one) ;
    al = -nwt*sgn ; bt = 1.0 ;
    blaswrap_sgemv(FALSE, nnbrs, nsts, al, &(Ast[0*nsts]), lda,
		   &(pn[ip*nstr]), nstr, bt, work, one) ;
#else  /*NBI_SINGLE_PRECISION*/
    blaswrap_dgemv(FALSE, nnbrs, nsts, al, &(Ast[1*nsts]), lda,
		   &(p[ip*pstr]) , pstr, bt, work, one) ;
    al = -nwt*sgn ; bt = 1.0 ;
    blaswrap_dgemv(FALSE, nnbrs, nsts, al, &(Ast[0*nsts]), lda,
		   &(pn[ip*nstr]), nstr, bt, work, one) ;
#endif /*NBI_SINGLE_PRECISION*/    
    
    for ( i = 0 ; i < nnbrs ; i ++ ) {
      f[fstr*nbrs[i]] += work[i] ;
    }
  }

  return 0 ;
}

static gpointer local_correction_thread(gpointer tdata)

{
  gpointer *mdata = tdata ;
  gint th  = GPOINTER_TO_INT(mdata[NBI_THREAD_MAIN_DATA_THREAD]) ;
  gint nth = GPOINTER_TO_INT(mdata[NBI_THREAD_MAIN_DATA_NTHREAD]) ;
  gpointer *data = mdata[NBI_THREAD_MAIN_DATA_DATA] ;
  NBI_REAL *work = mdata[NBI_THREAD_MAIN_DATA_WORK] ;
  nbi_matrix_t *m = data[NBI_THREAD_DATA_MATRIX  ] ;
  gint *idata = data[NBI_THREAD_DATA_INT] ;
  NBI_REAL *ddata = data[NBI_THREAD_DATA_REAL] ;
  NBI_REAL **dpdata = data[NBI_THREAD_DATA_REAL_POINTER] ;
  gint np, pt0, pt1, pstr, nstr, fstr ;
  NBI_REAL *p, *pn, *f, pwt, nwt, sgn ;

  np = nbi_surface_patch_number(m->s) ;

  pt0 = th*(np/nth) ;
  pt1 = (th+1)*(np/nth) ;

  if ( th == nth - 1 ) {
    if ( pt1 < np ) pt1 = np ;
  }

  p = dpdata[NBI_THREAD_DATA_REAL_PTR_P ]  ;
  pn = dpdata[NBI_THREAD_DATA_REAL_PTR_PN] ;
  f = dpdata[NBI_THREAD_DATA_REAL_PTR_F ]  ;

  pstr = idata[NBI_THREAD_DATA_INT_PSTR] ;
  nstr = idata[NBI_THREAD_DATA_INT_NSTR] ;
  fstr = idata[NBI_THREAD_DATA_INT_FSTR] ;
  
  pwt = ddata[NBI_THREAD_DATA_REAL_WT1 ] ;
  nwt = ddata[NBI_THREAD_DATA_REAL_WT2 ] ;
  sgn = ddata[NBI_THREAD_DATA_REAL_SIGN] ;
  
  local_matrix_correction(m, p, pstr, pwt, pn, nstr, nwt, sgn,
			  pt0, pt1, f, fstr, work) ;  
  
  return NULL ;
}

gint NBI_FUNCTION_NAME(nbi_surface_greens_identity_laplace)(nbi_matrix_t *m,
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
  gint i, nnmax, nth ;
  NBI_REAL sgn ;

  sgn = 1.0 ;
  if ( m->tree != NULL ) { sgn = -1.0 ; }
  nwt *= sgn ;

  if ( nthreads < 0 ) nth = g_get_num_processors() ; else nth = nthreads ;
  
  NBI_FUNCTION_NAME(nbi_matrix_upsample_laplace)(m,
						 p , pstr, pwt,
						 pn, nstr, nwt) ;
  /*point source approximation (FMM handled internally)*/
  /* fprintf(stderr, "%s: calling point source summation\n", __FUNCTION__) ; */
  point_source_summation(m, f, fstr, work, nth) ;
  /* fprintf(stderr, "%s: point source summation complete\n", __FUNCTION__) ; */

  nnmax = nbi_matrix_neighbour_number_max(m) ;
  
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

    dpdata[NBI_THREAD_DATA_REAL_PTR_P ] = p  ;
    dpdata[NBI_THREAD_DATA_REAL_PTR_PN] = pn ;
    dpdata[NBI_THREAD_DATA_REAL_PTR_F ] = f  ;

    idata[NBI_THREAD_DATA_INT_PSTR]   = pstr ;
    idata[NBI_THREAD_DATA_INT_NSTR]   = nstr ;
    idata[NBI_THREAD_DATA_INT_FSTR]   = fstr ;
  
    ddata[NBI_THREAD_DATA_REAL_WT1 ] = pwt ;
    ddata[NBI_THREAD_DATA_REAL_WT2 ] = nwt ;
    ddata[NBI_THREAD_DATA_REAL_SIGN] = sgn ;

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
	&(work[i*nnmax]) ;
      main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_NTHREAD] =
	GINT_TO_POINTER(nth) ;
      
      threads[i] = g_thread_new(NULL, local_correction_thread,
				&(main_data[NBI_THREAD_MAIN_DATA_SIZE*i])) ;
    }

    /*make sure all threads complete before we move on*/
    for ( i = 0 ; i < nth ; i ++ ) g_thread_join(threads[i]) ;
  } else {
    local_matrix_correction(m, p, pstr, pwt, pn, nstr, nwt, sgn,
			    0, nbi_surface_patch_number(m->s), f, fstr,
			    work) ;  
  }
  
#else /*_OPENMP*/
  local_matrix_correction(m, p, pstr, pwt, pn, nstr, nwt, sgn,
			  0, nbi_surface_patch_number(m->s), f, fstr, work) ;  
#endif /*_OPENMP*/
  
  return 0 ;
}

static gint upsample_sources_single(nbi_surface_t *s,
				    NBI_REAL *x, gint xstr,
				    NBI_REAL *wt, gint wstr, gint *idxu,
				    NBI_REAL *pu, gint pustr, NBI_REAL pwt,
				    NBI_REAL *pnu, gint nustr, NBI_REAL nwt)

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
		   &(x[i*xstr]), xstr, bt, &(pnu[j*pustr]), nustr) ;
#endif /*NBI_SINGLE_PRECISION*/
  }

  nu = idxu[nbi_surface_patch_number(s)] ;
  for ( i = 0 ; i < nu ; i ++ ) {
    pu [i*pustr] = 0.0 ;
    pnu[i*nustr] *= nwt*wt[i*wstr] ;
  }
  
  return 0 ;
}

static gint upsample_sources_double(nbi_surface_t *s,
				    NBI_REAL *x, gint xstr,
				    NBI_REAL *wt, gint wstr, gint *idxu,
				    NBI_REAL *pu, gint pustr, NBI_REAL pwt,
				    NBI_REAL *pnu, gint nustr, NBI_REAL nwt)

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
    		   &(x[i*xstr]), xstr, bt, &(pu[j*pustr]), pustr) ;
#endif /*NBI_SINGLE_PRECISION*/
  }

  nu = idxu[nbi_surface_patch_number(s)] ;
  for ( i = 0 ; i < nu ; i ++ ) {
    pnu[i*nustr] = 0.0 ;
    pu [i*pustr] *= pwt*wt[i*wstr] ;
  }
  
  return 0 ;
}

static gint local_matrix_multiply(nbi_matrix_t *m, NBI_REAL *p, gint pstr,
				  gint off, NBI_REAL al, NBI_REAL bt,
				  NBI_REAL *work,
				  NBI_REAL *f, gint fstr,
				  gint pt0, gint pt1,
				  gint nthreads)

{
  gint i, ip, nnbrs, *nbrs, one = 1, nsts, lda, pt ;
  NBI_REAL *Ast, *A ;

  A = (NBI_REAL *)(m->Ast) ;
  A = &(A[off]) ;
  /*loop on patches treated as sources*/
  nsts = nbi_surface_patch_node_number(m->s, 0) ;
  lda = 2*nsts ;

  /* for ( pt = pt0 ; pt < nbi_surface_patch_number(m->s) ; pt ++ ) { */
  for ( pt = pt0 ; pt < pt1 ; pt ++ ) {
    ip = nbi_surface_patch_node(m->s, pt) ;
    nnbrs = m->idxp[pt+1] - m->idxp[pt] ;
    nbrs = &(m->idx[m->idxp[pt]]) ;
    Ast = &(A[2*nsts*(m->idxp[pt])]) ;
    
#ifdef NBI_SINGLE_PRECISION
    g_assert_not_reached() ; /*unchecked code*/
    blaswrap_sgemv(FALSE, nnbrs, nsts, al, Ast, lda,
		   &(p[ip*pstr]) , pstr, bt, work, one) ;
#else  /*NBI_SINGLE_PRECISION*/
    blaswrap_dgemv(FALSE, nnbrs, nsts, al, Ast, lda,
		   &(p[ip*pstr]) , pstr, bt, work, one) ;
#endif /*NBI_SINGLE_PRECISION*/    
    
    for ( i = 0 ; i < nnbrs ; i ++ ) {
      f[fstr*nbrs[i]] += work[i] ;
    }
  }

  return 0 ;
}

#ifdef _OPENMP

static gpointer matrix_multiply_thread(gpointer tdata)

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
  NBI_REAL *p, *f, al, bt ;
  gint pt0, pt1, np, pstr, off, fstr ;
  
  np = nbi_surface_patch_number(m->s) ;
  pt0 = th*(np/nth) ;
  pt1 = (th+1)*(np/nth) ;

  if ( th == nth - 1 ) {
    if ( pt1 < np ) pt1 = np ;
  }

  p = dpdata[NBI_THREAD_DATA_REAL_PTR_P] ;
  f = dpdata[NBI_THREAD_DATA_REAL_PTR_F] ;

  pstr = idata[NBI_THREAD_DATA_INT_PSTR] ;
  off  = idata[NBI_THREAD_DATA_INT_OFFSET] ;
  fstr = idata[NBI_THREAD_DATA_INT_FSTR] ;
  
  al = ddata[NBI_THREAD_DATA_REAL_WT1] ;
  bt = ddata[NBI_THREAD_DATA_REAL_WT2] ;
  
  local_matrix_multiply(m, p, pstr, off, al, bt, work, f, fstr,
			pt0, pt1, nth) ;

  return NULL ;
}

static gint local_matrix_multiply_thread(nbi_matrix_t *m,
					 NBI_REAL *p, gint pstr,
					 gint off, NBI_REAL al, NBI_REAL bt,
					 NBI_REAL *work,
					 NBI_REAL *f, gint fstr,
					 gint nthreads)
{
  gint nth, i, nnmax, idata[NBI_THREAD_DATA_INT_SIZE] ;
  NBI_REAL
    ddata[NBI_THREAD_DATA_REAL_SIZE],
    *dpdata[NBI_THREAD_DATA_REAL_PTR_SIZE] ;
  GThread *threads[NBI_THREAD_NUMBER_MAX] ;
  gpointer data[NBI_THREAD_DATA_SIZE],
    main_data[NBI_THREAD_NUMBER_MAX*NBI_THREAD_MAIN_DATA_SIZE] ;

  nnmax = nbi_matrix_neighbour_number_max(m) ;
  
  if ( nthreads < 0 ) nth = g_get_num_processors() ; else nth = nthreads ;

  dpdata[NBI_THREAD_DATA_REAL_PTR_P] = p ;
  dpdata[NBI_THREAD_DATA_REAL_PTR_F] = f ;

  idata[NBI_THREAD_DATA_INT_PSTR]   = pstr ;
  idata[NBI_THREAD_DATA_INT_OFFSET] = off ;
  idata[NBI_THREAD_DATA_INT_FSTR]   = fstr ;
  
  ddata[NBI_THREAD_DATA_REAL_WT1] = al ;
  ddata[NBI_THREAD_DATA_REAL_WT2] = bt ;

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
      &(work[i*nnmax]) ;
    main_data[NBI_THREAD_MAIN_DATA_SIZE*i+NBI_THREAD_MAIN_DATA_NTHREAD] =
      GINT_TO_POINTER(nth) ;

    threads[i] = g_thread_new(NULL, matrix_multiply_thread,
			      &(main_data[NBI_THREAD_MAIN_DATA_SIZE*i])) ;    
  }

  for ( i = 0 ; i < nth ; i ++ ) g_thread_join(threads[i]) ;
  
  return 0 ;
}
#endif /*_OPENMP*/

static gint matrix_multiply_single(nbi_matrix_t *m,
				   NBI_REAL *pn, gint nstr, NBI_REAL nwt,
				   NBI_REAL *f, gint fstr,
				   gint nthreads,
				   NBI_REAL *work) 

{
  gint nth, ustr, *idxu, pustr, nustr ;
  NBI_REAL al, bt, *xu, *pu, *pnu, sgn ;
  nbi_surface_t *s ;
  
  ustr = m->ustr ;
  idxu = m->idxu ;
  pustr = m->pstr ; nustr = m->nstr ;
  xu = (NBI_REAL *)(m->xu) ;
  pu  = (NBI_REAL *)(m->p) ;
  pnu = (NBI_REAL *)(m->pn) ;
  s = m->s ;

  if ( nthreads < 0 ) nth = g_get_num_processors() ; else nth = nthreads ;

  sgn = 1.0 ; nwt = -nwt ;
  if ( m->tree != NULL ) { sgn = -1.0 ; }
  nwt *= sgn ;

  upsample_sources_single(s, pn, nstr, &(xu[6]), ustr, idxu,
			  pu, pustr, 0.0, pnu, nustr, nwt) ;
			  
  /*point source approximation (FMM handled internally)*/
  /* fprintf(stderr, "%s: calling point source summation\n", __FUNCTION__) ; */
  point_source_summation(m, f, fstr, work, nth) ;
  /* fprintf(stderr, "%s: point source summation complete\n", __FUNCTION__) ; */
  /* point_source_summation(m, f, fstr, work, nth) ; */

  /*local corrections*/
  al = -nwt*sgn ; bt = 0.0 ;

#ifdef _OPENMP
  local_matrix_multiply_thread(m, pn, nstr, 0, al, bt, work, f, fstr,
			       nthreads) ;
#else /*_OPENMP*/
  local_matrix_multiply(m, pn, nstr, 0, al, bt, work, f, fstr,
			0, nbi_surface_patch_number(s), nthreads) ;
#endif /*_OPENMP*/
  
  return 0 ;
}

static gint matrix_multiply_double(nbi_matrix_t *m,
				   NBI_REAL *p, gint pstr, NBI_REAL pwt,
				   NBI_REAL *f, gint fstr,
				   gint nthreads,
				   NBI_REAL *work) 

{
  gint nsts, nth ;
  gint ustr, *idxu, pustr, nustr ;
  NBI_REAL al, bt, *xu, *pu, *pnu ;
  nbi_surface_t *s ;
  
  ustr = m->ustr ;
  idxu = m->idxu ; pustr = m->pstr ; nustr = m->nstr ;
  xu = (NBI_REAL *)(m->xu) ;
  pu  = (NBI_REAL *)(m->p) ;
  pnu = (NBI_REAL *)(m->pn) ;
  s = m->s ;
  
  if ( nthreads < 0 ) nth = g_get_num_processors() ; else nth = nthreads ;

  upsample_sources_double(s, p, pstr, &(xu[6]), ustr, idxu,
			  pu , pustr, pwt, pnu, nustr, 0.0) ;
  /*point source approximation (FMM handled internally)*/
  /* fprintf(stderr, "%s: calling point source summation\n", __FUNCTION__) ; */
  point_source_summation(m, f, fstr, work, nth) ;
  /* fprintf(stderr, "%s: point source summation complete\n", __FUNCTION__) ; */
  /* point_source_summation(m, f, fstr, work, nth) ; */

  /*local corrections*/
  nsts = nbi_surface_patch_node_number(s, 0) ;
  al =  pwt ; bt = 0.0 ;

#ifdef _OPENMP
  local_matrix_multiply_thread(m, p, pstr, nsts, al, bt, work, f, fstr,
			       nthreads) ;
#else /*_OPENMP*/
  local_matrix_multiply(m, p, pstr, nsts, al, bt, work, f, fstr,
			0, nbi_surface_patch_number(s), nthreads) ;
#endif /*_OPENMP*/
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_matrix_multiply_laplace)(nbi_matrix_t *A,
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
  NBI_REAL diag ;
  
  g_assert(A->problem == NBI_PROBLEM_LAPLACE) ;

  diag = (NBI_REAL)(A->diag)*al ;

#ifdef NBI_SINGLE_PRECISION
  g_assert_not_reached() ; /*untested code*/
  blaswrap_sscal(nbi_surface_node_number(A->s), bt, y, ystr) ;
  blaswrap_saxpy(nbi_surface_node_number(A->s), diag, x, xstr, y, ystr) ;
#else /*NBI_SINGLE_PRECISION*/
  blaswrap_dscal(nbi_surface_node_number(A->s), bt, y, ystr) ;
  blaswrap_daxpy(nbi_surface_node_number(A->s), diag, x, xstr, y, ystr) ;
#endif /*NBI_SINGLE_PRECISION*/
  
  switch ( A->potential ) {
  case NBI_POTENTIAL_UNDEFINED:
    g_error("%s: required potential must be set for multiplication",
  	    __FUNCTION__) ;
    break ;
  case NBI_POTENTIAL_SINGLE:
    matrix_multiply_single(A, x, xstr, al, y, ystr, nthreads, work) ;
    return 0 ;
    break ;
  case NBI_POTENTIAL_DOUBLE:
    matrix_multiply_double(A, x, xstr, al, y, ystr, nthreads, work) ;
    return 0 ;
    break ;
  }
  
  return 0 ;
}
#endif
