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

#include "nbi-private.h"

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

static gint point_source_summation(nbi_matrix_t *m,
				   NBI_REAL *f, gint fstr,
				   NBI_REAL *work)
							    
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
    for ( level = depth ; level >= 3 ; level -- ) {
      wbfmm_laplace_upward_pass(tree, shifts, level, work) ;
    }  
    for ( level = 2 ; level <= depth ; level ++ ) {      
      wbfmm_laplace_downward_pass(tree, shifts, level, work) ;
    }
    if ( targets != NULL ) {
      wbfmm_target_list_local_field(targets, pnu, nustr, pu, pustr, f, fstr) ;
    } else {
      for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
	guint64 box ;
	box = wbfmm_point_box(tree, tree->depth,
			      (NBI_REAL *)nbi_surface_node(s, i)) ;
	wbfmm_tree_laplace_box_local_field(tree, tree->depth, box,
					   (NBI_REAL *)nbi_surface_node(s,i),
					   &(f[i*fstr]),
					   pnu, nustr,
					   &(xu[3]), ustr,
					   pu, pustr,
					   TRUE, work) ;
      }
    }
  }

  return 0 ;
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
							    NBI_REAL *work) 

{
  gint i, ip, nsts, nu, lda, one = 1, pt, *nbrs, nnbrs ;
  gint *idx, *idxp, ustr, *idxu, pustr, nustr ;
  NBI_REAL *A, *xu, *pu, *pnu, *Ast, al, bt, sgn ;
  nbi_surface_t *s ;
  
  idx = m->idx ; idxp = m->idxp ; ustr = m->ustr ;
  idxu = m->idxu ; pustr = m->pstr ; nustr = m->nstr ;
  A = (NBI_REAL *)(m->Ast) ;
  xu = (NBI_REAL *)(m->xu) ;
  pu  = (NBI_REAL *)(m->p) ;
  pnu = (NBI_REAL *)(m->pn) ;
  s = m->s ;

  sgn = 1.0 ;
  if ( m->tree != NULL ) { sgn = -1.0 ; }
  nwt *= sgn ;
  
  NBI_FUNCTION_NAME(nbi_matrix_upsample_laplace)(m,
						 p , pstr, pwt,
						 pn, nstr, nwt) ;
  /*point source approximation (FMM handled internally)*/
  point_source_summation(m, f, fstr, work) ;

  /*local corrections*/
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    /*loop on patches treated as sources*/
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

    ip = idxu[pt] ; nu = idxu[pt+1] - idxu[pt] ;
    
    for ( i = 0 ; i < nnbrs ; i ++ ) {
      f[fstr*nbrs[i]] += work[i] ;
      point_source_field_laplace(&(xu[ip*ustr]), ustr, nu,
      				 &(pu [ip*pustr]), pustr, 1.0,
      				 &(pnu[ip*nustr]), nustr, sgn,
      				 (NBI_REAL *)nbi_surface_node(s, nbrs[i]),
      				 -1.0, &(f[fstr*nbrs[i]])) ;
    }
  }
  
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

static gint matrix_multiply_single(nbi_matrix_t *m,
				   NBI_REAL *pn, gint nstr, NBI_REAL nwt,
				   NBI_REAL *f, gint fstr,
				   NBI_REAL *work) 

{
  gint i, ip, nsts, nu, lda, one = 1, pt ;
  gint *nbrs, nnbrs, *idx, *idxp, ustr, *idxu, pustr, nustr ;
  NBI_REAL *Ast, al, bt, *A, *xu, *pu, *pnu, sgn ;
  nbi_surface_t *s ;
  
  idx = m->idx ; idxp = m->idxp ; ustr = m->ustr ;
  idxu = m->idxu ; pustr = m->pstr ; nustr = m->nstr ;
  A = (NBI_REAL *)(m->Ast) ;
  xu = (NBI_REAL *)(m->xu) ;
  pu  = (NBI_REAL *)(m->p) ;
  pnu = (NBI_REAL *)(m->pn) ;
  s = m->s ;

  sgn = 1.0 ; nwt = -nwt ;
  if ( m->tree != NULL ) { sgn = -1.0 ; }
  nwt *= sgn ;

  upsample_sources_single(s, pn, nstr, &(xu[6]), ustr, idxu,
			  pu, pustr, 0.0, pnu, nustr, nwt) ;
			  
  /*point source approximation (FMM handled internally)*/
  point_source_summation(m, f, fstr, work) ;

  /*local corrections*/
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    /*loop on patches treated as sources*/
    nsts = nbi_surface_patch_node_number(s, pt) ;

    ip = nbi_surface_patch_node(s, pt) ;
    nnbrs = idxp[pt+1] - idxp[pt] ;
    nbrs = &(idx[idxp[pt]]) ;
    Ast = &(A[2*nsts*idxp[pt]]) ;

    lda = 2*nsts ;
    al = -nwt*sgn ; bt = 0.0 ;
#ifdef NBI_SINGLE_PRECISION
    g_assert_not_reached() ; /*unchecked code*/
    blaswrap_sgemv(FALSE, nnbrs, nsts, al, &(Ast[0*nsts]), lda,
		   &(pn[ip*nstr]), nstr, bt, work, one) ;
#else  /*NBI_SINGLE_PRECISION*/
    blaswrap_dgemv(FALSE, nnbrs, nsts, al, &(Ast[0*nsts]), lda,
		   &(pn[ip*nstr]), nstr, bt, work, one) ;
#endif /*NBI_SINGLE_PRECISION*/    

    ip = idxu[pt] ; nu = idxu[pt+1] - idxu[pt] ;
    
    for ( i = 0 ; i < nnbrs ; i ++ ) {
      f[fstr*nbrs[i]] += work[i] ;
      point_source_field_laplace(&(xu[ip*ustr]), ustr, nu,
      				 &(pu [ip*pustr]), pustr, 1.0,
      				 &(pnu[ip*nustr]), nustr, sgn,
      				 (NBI_REAL *)nbi_surface_node(s, nbrs[i]),
      				 -1.0, &(f[fstr*nbrs[i]])) ;
    }
  }
  
  return 0 ;
}

static gint matrix_multiply_double(nbi_matrix_t *m,
				   NBI_REAL *p, gint pstr, NBI_REAL pwt,
				   NBI_REAL *f, gint fstr,
				   NBI_REAL *work) 

{
  gint i, ip, nsts, nu, lda, one = 1, pt ;
  gint *nbrs, nnbrs, *idx, *idxp, ustr, *idxu, pustr, nustr ;
  NBI_REAL *Ast, al, bt, *A, *xu, *pu, *pnu ;
  nbi_surface_t *s ;
  
  idx = m->idx ; idxp = m->idxp ; ustr = m->ustr ;
  idxu = m->idxu ; pustr = m->pstr ; nustr = m->nstr ;
  A = (NBI_REAL *)(m->Ast) ;
  xu = (NBI_REAL *)(m->xu) ;
  pu  = (NBI_REAL *)(m->p) ;
  pnu = (NBI_REAL *)(m->pn) ;
  s = m->s ;
  
  upsample_sources_double(s, p, pstr, &(xu[6]), ustr, idxu,
			  pu , pustr, pwt, pnu, nustr, 0.0) ;
  /*point source approximation (FMM handled internally)*/
  point_source_summation(m, f, fstr, work) ;

  /*local corrections*/
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    /*loop on patches treated as sources*/
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
#else  /*NBI_SINGLE_PRECISION*/
    blaswrap_dgemv(FALSE, nnbrs, nsts, al, &(Ast[1*nsts]), lda,
		   &(p[ip*pstr]) , pstr, bt, work, one) ;
#endif /*NBI_SINGLE_PRECISION*/    

    ip = idxu[pt] ; nu = idxu[pt+1] - idxu[pt] ;
    
    for ( i = 0 ; i < nnbrs ; i ++ ) {
      f[fstr*nbrs[i]] += work[i] ;
      point_source_field_laplace(&(xu[ip*ustr]), ustr, nu,
      				 &(pu [ip*pustr]), pustr, 1.0,
      				 &(pnu[ip*nustr]), nustr, 1.0,
      				 (NBI_REAL *)nbi_surface_node(s, nbrs[i]),
      				 -1.0, &(f[fstr*nbrs[i]])) ;
    }
  }
  
  return 0 ;
}

gint NBI_FUNCTION_NAME(nbi_matrix_multiply_laplace)(nbi_matrix_t *A,
						    NBI_REAL *x, gint xstr,
						    NBI_REAL al,
						    NBI_REAL *y, gint ystr,
						    NBI_REAL bt,
						    NBI_REAL *work)

/*
  y := al*A*x + bt*y
*/
  
{
  NBI_REAL diag ;
  
  g_assert(A->problem == NBI_PROBLEM_LAPLACE) ;

  diag = (NBI_REAL)(A->diag)*al ;

#ifdef NBI_SINGLE_PRECISION
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
    matrix_multiply_single(A, x, xstr, al, y, ystr, work) ;
    return 0 ;
    break ;
  case NBI_POTENTIAL_DOUBLE:
    matrix_multiply_double(A, x, xstr, al, y, ystr, work) ;
    return 0 ;
    break ;
  }
  
  return 0 ;
}
