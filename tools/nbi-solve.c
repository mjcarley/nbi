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

#define LOCAL_CUTOFF_RADIUS 1e-6

GTimer *timer ;
gchar *progname ;

gint nbi_surface_greens_identity_laplace(nbi_surface_t *s, gint *idx,
					 gint *idxp,
					 gdouble *A,
					 gdouble *p , gint pstr,
					 gdouble *pn, gint nstr,
					 gdouble *xu, gint ustr, gint *idxu,
					 gdouble *pu, gint pustr, gdouble pwt,
					 gdouble *pnu, gint nustr, gdouble nwt,
					 wbfmm_tree_t *tree,
					 wbfmm_target_list_t *targets,
					 wbfmm_shift_operators_t *shifts,
					 gdouble *work,
					 gdouble *f) ;

static gint make_sources(nbi_surface_t *s,
			 gdouble *xs, gdouble q,
			 gdouble *p , gint pstr,
			 gdouble *pn, gint nstr)

{
  gdouble R, Rn, r[3], *x, *n, G ;
  gint i ;

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    x = nbi_surface_node(s, i) ;
    n = nbi_surface_normal(s, i) ;

    nbi_vector_diff(r, x, xs) ;
    R = nbi_vector_length(r) ;
    Rn = nbi_vector_scalar(r,n)/R ;
    /* G = q*0.25*M_1_PI/R*nbi_surface_node_weight(s,i) ; */
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

static gint point_source_field_laplace(/* nbi_surface_t *s, */
				       gdouble *xs, gint xstr, gint ns,
				       gdouble *p , gint pstr, gdouble pwt,
				       gdouble *pn, gint nstr, gdouble nwt,
				       gdouble *x, gdouble wt, gdouble *f)
  
{
  gint i ;  
  gdouble R, r[3] ;

  for ( i = 0 ; i < ns ; i ++ ) {
    nbi_vector_diff(r, x, &(xs[i*xstr])) ;
    R = nbi_vector_length(r) ;
    if ( R > LOCAL_CUTOFF_RADIUS ) {
      *f += wt*(pwt*p[i*pstr]*nbi_vector_scalar(r,&(xs[i*xstr+3]))/R/R -
		nwt*pn[i*nstr])*0.25*M_1_PI/R ;
    }
  }
  
  return 0 ;
}

gint nbi_surface_greens_identity_laplace(nbi_surface_t *s, gint *idx,
					 gint *idxp,
					 gdouble *A,
					 gdouble *p , gint pstr,
					 gdouble *pn, gint nstr,
					 gdouble *xu, gint ustr, gint *idxu,
					 gdouble *pu, gint pustr, gdouble pwt,
					 gdouble *pnu, gint nustr, gdouble nwt,
					 wbfmm_tree_t *tree,
					 wbfmm_target_list_t *targets,
					 wbfmm_shift_operators_t *shifts,
					 gdouble *work,
					 gdouble *f)

{
  gint i, ip, nsts, nu, np ;
  gint lda, one = 1, pt = 0 ;
  gdouble *Ast,  al, bt ;
  gint *nbrs, nnbrs ;
  
  memset(f, 0, nbi_surface_node_number(s)*sizeof(gdouble)) ;

  /* work = (gdouble *)g_malloc(16384*sizeof(gdouble)) ; */

  if ( tree == NULL ) {
    fprintf(stderr, "%s: starting point source summation; t=%lg\n",
	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ;
    np = nbi_surface_patch_number(s) ;
    nu = idxu[np] ;
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      point_source_field_laplace(xu, ustr, nu,
				 pu, pustr, pwt,
				 pnu, nustr, nwt,
				 nbi_surface_node(s, i),
				 1.0, &(f[i])) ;
    }
    fprintf(stderr, "%s: point source summation complete; t=%lg\n",
	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ;  
  } else {
    /*FMM summation*/
    gint level, depth ;
    /* nsrc = idxu[nbi_surface_patch_number(s)] ; */

    depth = wbfmm_tree_depth(tree) ;
    fprintf(stderr, "%s: initializing leaf expansions; %lg\n",
	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ;
    wbfmm_tree_laplace_leaf_expansions(tree,
				       pnu, nustr,
				       &(xu[3]), ustr,
				       pu, pustr,
				       TRUE, work) ;  
    fprintf(stderr, "%s: leaf expansions initialized; %lg\n",
	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ;
    
    fprintf(stderr, "%s: upward pass; %lg\n",
	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ;
    for ( level = depth ; level >= 3 ; level -- ) {
      wbfmm_laplace_upward_pass(tree, shifts, level, work) ;
    }  
    fprintf(stderr, "%s: upward pass completed; %lg\n",
	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ;
    
    fprintf(stderr, "%s: downward pass; %lg\n",
	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ;
    for ( level = 2 ; level <= depth ; level ++ ) {
      wbfmm_laplace_downward_pass(tree, shifts, level, work) ;
    }
    fprintf(stderr, "%s: downward pass completed; %lg\n",
	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ;

    if ( targets != NULL ) {
      wbfmm_target_list_local_field(targets, pnu, nustr, pu, pustr, f, 1) ;
    } else {
      for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
	guint64 box ;
	box = wbfmm_point_box(tree, tree->depth, nbi_surface_node(s, i)) ;
	wbfmm_tree_laplace_box_local_field(tree, tree->depth, box,
					   nbi_surface_node(s,i),
					   &(f[i]),
					   pnu, nustr,
					   &(xu[3]), ustr,
					   pu, pustr,
					   TRUE, work) ;
      }
    }
  }

  /* return 0 ; */
  
  fprintf(stderr, "%s: starting local corrections; t=%lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL)) ;

  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    /*loop on patches treated as sources*/
    nsts = nbi_surface_patch_node_number(s, pt) ;

    ip = nbi_surface_patch_node(s, pt) ;
    nnbrs = idxp[pt+1] - idxp[pt] ;
    nbrs = &(idx[idxp[pt]]) ;
    Ast = &(A[2*nsts*idxp[pt]]) ;

    lda = 2*nsts ;
    al =  1.0 ; bt = 0.0 ;
    blaswrap_dgemv(FALSE, nnbrs, nsts, al, &(Ast[1*nsts]), lda,
		   &(p[ip*pstr]) , pstr, bt, work, one) ;
    al = -1.0 ; bt = 1.0 ;
    blaswrap_dgemv(FALSE, nnbrs, nsts, al, &(Ast[0*nsts]), lda,
		   &(pn[ip*nstr]), nstr, bt, work, one) ;

    ip = idxu[pt] ; nu = idxu[pt+1] - idxu[pt] ;
    for ( i = 0 ; i < nnbrs ; i ++ ) {
      f[nbrs[i]] += work[i] ;
      point_source_field_laplace(&(xu[ip*ustr]), ustr, nu,
      				 &(pu [ip*pustr]), pustr, pwt,
      				 &(pnu[ip*nustr]), nustr, nwt,
      				 nbi_surface_node(s, nbrs[i]),
      				 -1.0, &(f[nbrs[i]])) ;
    }
  }
  
  return 0 ;
}

static gint read_correction_matrices(FILE *f, gint **idx, gint **idxp,
				     gdouble **Ast)

{
  gint i, j, np, nst ;

  fscanf(f, "%d", &np) ;
  fscanf(f, "%d", &nst) ;

  *idxp = (gint *)g_malloc0((np+1)*sizeof(gint)) ;

  for ( i = 0 ; i < np+1 ; i ++ ) {
    fscanf(f, "%d", &j) ;
    g_assert(j == i) ;
    fscanf(f, "%d", &((*idxp)[i])) ;
  }
  
  *idx = (gint *)g_malloc0(((*idxp)[np])*sizeof(gint)) ;
  for ( i = 0 ; i < (*idxp)[np] ; i ++ ) {
    fscanf(f, "%d", &((*idx)[i])) ;
  }

  *Ast = (gdouble *)g_malloc0(nst*2*((*idxp)[np])*sizeof(gdouble)) ;
  for ( i = 0 ; i < 2*nst*((*idxp)[np]) ; i ++ ) {
    fscanf(f, "%lg", &((*Ast)[i])) ;
  }  

  return np ;
}

static gint read_upsampled_patches(FILE *f, gint **idxu,
				   gdouble **xu, gint *ustr)

{
  gint np, i, j ;

  fscanf(f, "%d", &np) ;
  fscanf(f, "%d",  ustr) ;

  *idxu = (gint *)g_malloc((np+1)*sizeof(gint)) ;
  for ( i = 0 ; i < np + 1 ; i ++ ) {
    fscanf(f, "%d", &j) ;
    g_assert(j == i) ;
    fscanf(f, "%d", &((*idxu)[i])) ;
  }
  
  *xu   = (gdouble *)g_malloc((*idxu)[np]*(*ustr)*sizeof(gdouble)) ;
  for ( i = 0 ; i < (*idxu)[np]*(*ustr) ; i ++ )
    fscanf(f, "%lg", &((*xu)[i])) ;
  
  return 0 ;
}

static gint upsample_sources(nbi_surface_t *s,
			     gdouble *p, gint pstr, gdouble *pn, gint nstr,
			     gdouble *wt, gint wstr, gint *idxu,
			     gdouble *pu, gint pustr, gdouble pwt,
			     gdouble *pnu, gint nustr, gdouble nwt)

{
  gint i, j, pt, ns, nu ;
  gdouble *K, al, bt ;

  al =  1.0 ; bt = 0.0 ;
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    ns = nbi_surface_patch_node_number(s, pt) ;
    nu = idxu[pt+1] - idxu[pt] ;
    K = nbi_patch_upsample_matrix(ns, nu) ;

    i = nbi_surface_patch_node(s, pt) ;
    j = idxu[pt] ;
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(p[i*pstr]), pstr, bt, &(pu[j*pustr]), pustr) ;
    blaswrap_dgemv(FALSE, nu, ns, al, K, ns,
		   &(pn[i*pstr]), nstr, bt, &(pnu[j*pustr]), nustr) ;
    for ( i = 0 ; i < nu ; i ++ ) {
      pu [(j+i)*pustr] *= pwt*wt[(j+i)*wstr] ;
      pnu[(j+i)*nustr] *= nwt*wt[(j+i)*wstr] ;
    }
  }
  
  return 0 ;
}
		      
gint main(gint argc, gchar **argv)

{
  nbi_surface_t *s ;
  gint i, *idx, *idxp, ustr, *idxu, nsrc ;
  gdouble xs[512], *f, *xp, *src, t, *xu, *su ;
  gdouble emax, fmax, G, *Ast ;
  gdouble xtree[3], xtmax[3], D, dtree, pwt, nwt ;
  FILE *output, *input ;
  gchar ch, *gfile, *mfile ;
  wbfmm_tree_t *tree ;
  wbfmm_target_list_t *targets ;
  wbfmm_shift_operators_t *shifts ;
  gdouble *work ;
  gint fmm_work_size, nqfmm, order_fmm, order_inc ;
  guint depth, order[48] = {0}, order_s, order_r, order_max, field ;
  gboolean fmm, shift_bw, target_list, greens_id ;
  
  output = stdout ;
  mfile = NULL ; gfile = NULL ;
  
  dtree = 1e-2 ; fmm = FALSE ; nqfmm = 1 ; shift_bw = TRUE ;
  tree = NULL ; targets = NULL ; shifts = NULL ;
  field = WBFMM_FIELD_SCALAR ;
  target_list = FALSE ;
  greens_id = TRUE ;
  
  pwt = 1.0 ; nwt = 1.0 ;
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  order_fmm = 12 ; order_inc = 2 ; depth = 4 ;
  
  while ( (ch = getopt(argc, argv, "d:fGg:lm:o:T:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'd': order_inc = atoi(optarg) ; break ;
    case 'f': fmm = TRUE ; break ;
    case 'G': greens_id = TRUE ; break ;
    case 'g': gfile = g_strdup(optarg) ; break ;
    case 'l': target_list = TRUE ; break ;
    case 'm': mfile = g_strdup(optarg) ; break ;
    case 'o': order_fmm = atoi(optarg) ; break ;
    case 'T': depth = atoi(optarg) ; break ;
    }
  }

  if ( gfile == NULL ) gfile = g_strdup("geometry.dat") ;
  if ( mfile == NULL ) mfile = g_strdup("matrix.dat") ;

  if ( fmm == TRUE ) {
    /*FMM solver sums contributions from sources directly, so single
      layer potential sources need to be negative to give correct
      contribution to Green's identity*/
    pwt = 1.0 ; nwt = -1.0 ;
  }
  
  fprintf(stderr, "%s: reading geometry from %s\n", progname, gfile) ;
  input = fopen(gfile, "r") ;

  s = nbi_surface_read(input) ;

  fclose(input) ;

  fprintf(stderr, "%s: reading matrix from %s\n", progname, mfile) ;
  input = fopen(mfile, "r") ;

  read_upsampled_patches(input, &idxu, &xu, &ustr) ;
  read_correction_matrices(input, &idx, &idxp, &Ast) ;

  nsrc = idxu[nbi_surface_patch_number(s)] ;

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
  f   = (gdouble *)g_malloc0(nbi_surface_node_number(s)*  sizeof(gdouble)) ;

  xs[0] = 0.3 ; xs[1] = -0.4 ; xs[2] = 0.2 ;

  make_sources(s, xs, 1.0, &(src[0]), 2, &(src[1]), 2) ;

  su = (gdouble *)g_malloc0(2*idxu[nbi_surface_patch_number(s)]*
			    sizeof(gdouble)) ;
  upsample_sources(s, &(src[0]), 2, &(src[1]), 2,
		   &(xu[6]), ustr, idxu,
		   &(su[0]), 2, pwt, &(su[1]), 2, nwt) ;

  if ( fmm ) {
    gint fmmpstr ;
    wbfmm_source_t source ;

    source = WBFMM_SOURCE_MONOPOLE | WBFMM_SOURCE_DIPOLE ;
    /* order_s = 20 ; order_r = 20 ; */
    order_s = order_fmm ; order_r = order_fmm ;
    order[2*depth+0] = order_s ; 
    order[2*depth+1] = order_r ; 
    order_max = MAX(order_s, order_r) ;
    for ( i = depth-1 ; i > 0 ; i -- ) {
      order[2*i+0] = order[2*(i+1)+0] ;
      order[2*i+0] = order[2*(i+1)+0] + order_inc ;
      order[2*i+1] = order[2*(i+1)+1] + order_inc ;
      /* order[2*i+1] = order[2*(i+1)+1] ; */
      order_max = MAX(order_max, order[2*i+0]) ;
      order_max = MAX(order_max, order[2*i+1]) ;
    }
    
    fprintf(stderr, "%s: building tree; t=%lg\n",
	    progname, t = g_timer_elapsed(timer, NULL)) ;
    xtree[0] = xtree[1] = xtree[2] = 0.0 ;
    wbfmm_points_origin_width(xu, ustr, nsrc, xtree, xtmax, &D, TRUE) ;
    wbfmm_points_origin_width(nbi_surface_node(s,0), NBI_SURFACE_NODE_LENGTH,
			      nbi_surface_node_number(s),
			      xtree, xtmax, &D, FALSE) ;

    xtree[0] -= dtree ; xtree[1] -= dtree ; xtree[2] -= dtree ;
    D += 2.0*dtree ;
    fprintf(stderr, "%s: tree built; t=%lg\n",
	    progname, t = g_timer_elapsed(timer, NULL)) ;

    fmmpstr = ustr*sizeof(gdouble) ;
    /* fmmfstr = strf*sizeof(gdouble) ; */
    tree = wbfmm_tree_new(xtree, D, 2*nsrc) ;
    fmm_work_size = wbfmm_element_number_rotation(2*order_max) ;
    fmm_work_size = MAX(fmm_work_size, (order_max+1)*(order_max+1)*nqfmm*16) ;
    work = (gdouble *)g_malloc0(2*fmm_work_size*sizeof(gdouble)) ;
    fprintf(stderr, "%s: initializing shift rotation operators; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
    wbfmm_shift_angle_table_init() ;
    shifts = wbfmm_shift_operators_new(order_max, shift_bw, work) ;
    fprintf(stderr, "%s: shift rotation operators initialized; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
    
    fprintf(stderr, "%s: initializing coaxial translation coefficients; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
    wbfmm_laplace_coaxial_translate_init(order_max+2) ;
    fprintf(stderr, "%s: coaxial translation coefficients initialized; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;

    wbfmm_tree_add_points(tree,
			  (gpointer)xu, fmmpstr,
			  (gpointer)(&(xu[3])), fmmpstr, nsrc) ;
    for ( i = 0 ; i < depth ; i ++ ) wbfmm_tree_refine(tree) ;
    wbfmm_tree_problem(tree) = WBFMM_PROBLEM_LAPLACE ;
    wbfmm_tree_source_size(tree) = nqfmm ;
    for ( i = 1 ; i <= depth ; i ++ ) {
      wbfmm_tree_laplace_coefficient_init(tree, i,
					  order[2*i+1], order[2*i+0]) ;
    }

    if ( target_list ) {
      fprintf(stderr, "%s: initializing target point list; %lg\n",
	      progname, g_timer_elapsed(timer, NULL)) ;
      targets = wbfmm_target_list_new(tree, nbi_surface_node_number(s)) ;
      wbfmm_target_list_coefficients_init(targets, field) ;
      wbfmm_target_list_add_points(targets,
				   nbi_surface_node(s,0),
				   NBI_SURFACE_NODE_LENGTH*sizeof(gdouble),
				   nbi_surface_node_number(s)) ;
      wbfmm_laplace_target_list_local_coefficients(targets, source, work) ;
      fprintf(stderr, "%s: target point list initialized; %lg\n",
	      progname, g_timer_elapsed(timer, NULL)) ;
    }
  } else {
      work = (gdouble *)g_malloc0(16384*sizeof(gdouble)) ;    
  }

  if ( greens_id ) {
    fprintf(stderr, "%s: evaluating Green's identity; t=%lg\n",
	    progname, t = g_timer_elapsed(timer, NULL)) ;
    nbi_surface_greens_identity_laplace(s, idx, idxp, Ast,
					&(src[0]), 2, &(src[1]), 2,
					xu, ustr, idxu,
					&(su[0]), 2, pwt, &(su[1]), 2, nwt,
					tree, targets, shifts, work,
					f) ;
    fprintf(stderr, "%s: surface integration complete; t=%lg (%lg)\n",
	    progname,
	    g_timer_elapsed(timer, NULL), g_timer_elapsed(timer, NULL) - t) ;
    
    emax = 0.0 ; fmax = 0.0 ;
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      xp = nbi_surface_node(s,i) ;
      G = greens_function_laplace(xp, xs)*0.5 ;
      fmax = MAX(fmax, G) ;
      emax = MAX(emax, fabs(G - f[i])) ;
      fprintf(output, "%lg %lg %lg %lg %lg\n",
	      xp[0], xp[1], xp[2], f[i], fabs(G - f[i])) ;
    }
    
    fprintf(stderr, "L_inf norm: %lg\n", emax/fmax) ;

    return 0 ;
  }

  
  return 0 ;
}
