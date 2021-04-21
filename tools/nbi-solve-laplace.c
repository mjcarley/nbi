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

gint nbi_surface_greens_identity_laplace(nbi_matrix_t *m,
					 gdouble *p , gint pstr,
					 gdouble *pn, gint nstr,
					 gdouble pwt, gdouble nwt,
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

static gint point_source_field_laplace(gdouble *xs, gint xstr, gint ns,
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

gint nbi_surface_greens_identity_laplace(nbi_matrix_t *m,
					 gdouble *p , gint pstr,
					 gdouble *pn, gint nstr,
					 gdouble pwt, gdouble nwt,
					 gdouble *f,
					 gdouble *work) 

{
  gint i, ip, nsts, nu, np ;
  gint lda, one = 1, pt = 0 ;
  gdouble *Ast,  al, bt ;
  gint *nbrs, nnbrs ;
  nbi_surface_t *s ;
  gint *idx, *idxp, ustr, *idxu, pustr, nustr ;
  gdouble *A, *xu, *pu, *pnu ;
  wbfmm_tree_t *tree ;
  wbfmm_target_list_t *targets ;
  wbfmm_shift_operators_t *shifts ;
  gdouble sgn ;
  
  tree = m->tree ; targets = m->targets ; shifts = m->shifts ;

  idx = m->idx ; idxp = m->idxp ; ustr = m->ustr ;
  idxu = m->idxu ; pustr = m->pstr ; nustr = m->nstr ;
  A = (gdouble *)(m->Ast) ;
  xu = (gdouble *)(m->xu) ;
  pu  = (gdouble *)(m->p) ;
  pnu = (gdouble *)(m->pn) ;
  s = m->s ;

  np = nbi_surface_patch_number(s) ; nu = idxu[np] ;

  sgn = 1.0 ;
  if ( m->tree != NULL ) { nwt = -nwt ; sgn = -1.0 ; }
  nbi_matrix_upsample_laplace(m, p, pstr, pn, nstr, pwt, nwt) ;
  
  if ( tree == NULL ) {
    /* fprintf(stderr, "%s: starting point source summation; t=%lg\n", */
    /* 	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ; */
    np = nbi_surface_patch_number(s) ;
    nu = idxu[np] ;
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      point_source_field_laplace(xu, ustr, nu,
				 pu, pustr, 1.0,
				 pnu, nustr, sgn,
				 (NBI_REAL *)nbi_surface_node(s, i),
				 1.0, &(f[i])) ;
    }
    /* fprintf(stderr, "%s: point source summation complete; t=%lg\n", */
    /* 	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ;   */
  } else {
    /*FMM summation*/
    gint level, depth ;

    depth = wbfmm_tree_depth(tree) ;
    for ( level = 2 ; level <= depth ; level ++ ) {
      wbfmm_tree_coefficients_zero(tree, level) ;
    }
    /* fprintf(stderr, "%s: initializing leaf expansions; %lg\n", */
    /* 	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ; */
    wbfmm_tree_laplace_leaf_expansions(tree,
				       pnu, nustr,
				       &(xu[3]), ustr,
				       pu, pustr,
				       TRUE, work) ;  
    /* fprintf(stderr, "%s: leaf expansions initialized; %lg\n", */
    /* 	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ; */
    
    /* fprintf(stderr, "%s: upward pass; %lg\n", */
    /* 	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ; */
    for ( level = depth ; level >= 3 ; level -- ) {
      wbfmm_laplace_upward_pass(tree, shifts, level, work) ;
    }  
    /* fprintf(stderr, "%s: upward pass completed; %lg\n", */
    /* 	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ; */
    
    /* fprintf(stderr, "%s: downward pass; %lg\n", */
    /* 	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ; */
    for ( level = 2 ; level <= depth ; level ++ ) {      
      wbfmm_laplace_downward_pass(tree, shifts, level, work) ;
    }
    /* fprintf(stderr, "%s: downward pass completed; %lg\n", */
    /* 	    __FUNCTION__, g_timer_elapsed(timer, NULL)) ; */

    if ( targets != NULL ) {
      wbfmm_target_list_local_field(targets, pnu, nustr, pu, pustr, f, 1) ;
    } else {
      for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
	guint64 box ;
	box = wbfmm_point_box(tree, tree->depth,
			      (NBI_REAL *)nbi_surface_node(s, i)) ;
	wbfmm_tree_laplace_box_local_field(tree, tree->depth, box,
					   (NBI_REAL *)nbi_surface_node(s,i),
					   &(f[i]),
					   pnu, nustr,
					   &(xu[3]), ustr,
					   pu, pustr,
					   TRUE, work) ;
      }
    }
  }

  /* fprintf(stderr, "%s: starting local corrections; t=%lg\n", */
  /* 	  __FUNCTION__, g_timer_elapsed(timer, NULL)) ; */

  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    /*loop on patches treated as sources*/
    nsts = nbi_surface_patch_node_number(s, pt) ;

    ip = nbi_surface_patch_node(s, pt) ;
    nnbrs = idxp[pt+1] - idxp[pt] ;
    nbrs = &(idx[idxp[pt]]) ;
    Ast = &(A[2*nsts*idxp[pt]]) ;

    lda = 2*nsts ;
    al =  pwt ; bt = 0.0 ;
    blaswrap_dgemv(FALSE, nnbrs, nsts, al, &(Ast[1*nsts]), lda,
		   &(p[ip*pstr]) , pstr, bt, work, one) ;
    al = -nwt*sgn ; bt = 1.0 ;
    blaswrap_dgemv(FALSE, nnbrs, nsts, al, &(Ast[0*nsts]), lda,
		   &(pn[ip*nstr]), nstr, bt, work, one) ;
    
    ip = idxu[pt] ; nu = idxu[pt+1] - idxu[pt] ;
    for ( i = 0 ; i < nnbrs ; i ++ ) {
      f[nbrs[i]] += work[i] ;
      point_source_field_laplace(&(xu[ip*ustr]), ustr, nu,
      				 &(pu [ip*pustr]), pustr, 1.0,
      				 &(pnu[ip*nustr]), nustr, sgn,
      				 (NBI_REAL *)nbi_surface_node(s, nbrs[i]),
      				 -1.0, &(f[nbrs[i]])) ;
    }
  }
  
  return 0 ;
}

gint nbi_matrix_fmm_init(nbi_matrix_t *m,
			 nbi_problem_t problem,
			 wbfmm_shift_operators_t *shifts,
			 guint *order_s, gint sstr,
			 guint *order_r, gint rstr,
			 gint depth,
			 NBI_REAL dtree,
			 gboolean shift_bw,
			 NBI_REAL *work)

{
  gint fmmpstr, ustr, nsrc, i, nqfmm ;
  wbfmm_source_t source ;
  NBI_REAL xtree[3], xtmax[3], *xu, D ;
  guint order_max, field ;
  nbi_surface_t *s = m->s ;

  nqfmm = 1 ;
  order_max = 0 ;
  for ( i = 0 ; i <= depth ; i ++ ) {
    order_max = MAX(order_s[i*sstr], order_max) ;
    order_max = MAX(order_r[i*rstr], order_max) ;
  }

  xu = (NBI_REAL *)(m->xu) ;
  ustr = m->ustr ;
  nsrc = m->idxu[nbi_surface_patch_number(s)] ;
  
  source = WBFMM_SOURCE_MONOPOLE | WBFMM_SOURCE_DIPOLE ;
  field = WBFMM_FIELD_SCALAR ;
  xtree[0] = xtree[1] = xtree[2] = 0.0 ;
  wbfmm_points_origin_width(xu, ustr, nsrc, xtree, xtmax, &D, TRUE) ;
  wbfmm_points_origin_width((NBI_REAL *)nbi_surface_node(s,0),
			    NBI_SURFACE_NODE_LENGTH,
			    nbi_surface_node_number(s),
			    xtree, xtmax, &D, FALSE) ;

  xtree[0] -= dtree ; xtree[1] -= dtree ; xtree[2] -= dtree ;
  D += 2.0*dtree ;

  fmmpstr = ustr*sizeof(NBI_REAL) ;
  m->tree = wbfmm_tree_new(xtree, D, 2*nsrc) ;

  if ( shifts != NULL )
    m->shifts = shifts ;
  else
    m->shifts = wbfmm_shift_operators_new(order_max, shift_bw, work) ;

  wbfmm_tree_add_points(m->tree,
			(gpointer)xu, fmmpstr,
			(gpointer)(&(xu[3])), fmmpstr, nsrc) ;
  
  for ( i = 0 ; i < depth ; i ++ ) wbfmm_tree_refine(m->tree) ;

  switch ( problem ) {
  default: g_assert_not_reached() ; break ;
  case NBI_PROBLEM_LAPLACE:
    wbfmm_tree_problem(m->tree) = WBFMM_PROBLEM_LAPLACE ;
    wbfmm_tree_source_size(m->tree) = nqfmm ;
    for ( i = 1 ; i <= depth ; i ++ ) {
      wbfmm_tree_laplace_coefficient_init(m->tree, i,
					  order_r[i*rstr], order_s[i*sstr]) ;
    }
    
    m->targets = wbfmm_target_list_new(m->tree, nbi_surface_node_number(s)) ;
    wbfmm_target_list_coefficients_init(m->targets, field) ;
    wbfmm_target_list_add_points(m->targets,
				 (gpointer)nbi_surface_node(s,0),
				 NBI_SURFACE_NODE_LENGTH*sizeof(gdouble),
				 nbi_surface_node_number(s)) ;
  
    wbfmm_laplace_target_list_local_coefficients(m->targets, source, work) ;
    return 0 ;
    break ;
  }

  return 0 ;
}

gint nbi_matrix_multiply(nbi_matrix_t *A, nbi_potential_t p,
			 gdouble *x, gint xstr, gdouble al,
			 gdouble *y, gint ystr, gdouble bt,
			 gdouble *work)

/*
  y := al*A*x + bt*y
*/
  
{
  nbi_surface_t *s ;
  gint i ;

  g_assert(A->problem == NBI_PROBLEM_LAPLACE) ;

  s = A->s ;
  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) y[i*ystr] *= bt ;
  
  switch ( p ) {
  case NBI_POTENTIAL_SINGLE:
    /*single layer potential, source term is normal derivative and
      surface potential is zeroed*/
    nbi_surface_greens_identity_laplace(A, x, xstr, x, xstr,
					0.0, -al, y, work) ;
    break ;
  case NBI_POTENTIAL_DOUBLE:
    /*double layer potential, source term is surface potential and
      normal derivative is zeroed*/
    nbi_surface_greens_identity_laplace(A, x, xstr, x, xstr,
					al, 0.0, y, work) ;
    break ;
  }
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  nbi_surface_t *s ;
  nbi_matrix_t *matrix ;
  gint i ;
  gdouble xs[512], *f, *xp, *src, t ;
  gdouble emax, fmax, G ;
  gdouble dtree, pwt, nwt ;
  FILE *output, *input ;
  gchar ch, *gfile, *mfile ;
  gdouble *work ;
  gint fmm_work_size, nqfmm, order_fmm, order_inc ;
  guint depth, order[48] = {0}, order_s, order_r, order_max ;
  gboolean fmm, shift_bw, greens_id, layer_potentials ;
  
  output = stdout ;
  mfile = NULL ; gfile = NULL ;
  
  dtree = 1e-2 ; fmm = FALSE ; nqfmm = 1 ; shift_bw = TRUE ;
  greens_id = FALSE ; layer_potentials = FALSE ;
  
  pwt = 2.0 ; nwt = 2.0 ;
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  order_fmm = 12 ; order_inc = 2 ; depth = 4 ;
  
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
  f   = (gdouble *)g_malloc0(nbi_surface_node_number(s)*  sizeof(gdouble)) ;

  xs[0] = 0.3 ; xs[1] = -0.4 ; xs[2] = 0.2 ;

  make_sources(s, xs, 1.0, &(src[0]), 2, &(src[1]), 2) ;

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
					&(src[0]), 2, &(src[1]), 2,
					pwt, nwt,
					f, work) ;
    fprintf(stderr, "%s: surface integration complete; t=%lg (%lg)\n",
	    progname,
	    g_timer_elapsed(timer, NULL), g_timer_elapsed(timer, NULL) - t) ;
    
    emax = 0.0 ; fmax = 0.0 ;
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      xp = (NBI_REAL *)nbi_surface_node(s,i) ;
      G = greens_function_laplace(xp, xs) ;
      fmax = MAX(fmax, G) ;
      emax = MAX(emax, fabs(G - f[i])) ;
      fprintf(output, "%lg %lg %lg %lg %lg\n",
	      xp[0], xp[1], xp[2], f[i], fabs(G - f[i])) ;
    }
    
    fprintf(stderr, "L_inf norm: %lg\n", emax/fmax) ;

    return 0 ;
  }

  if ( layer_potentials ) {
    fprintf(stderr, "%s: evaluating double-layer potential; t=%lg\n",
    	    progname, t = g_timer_elapsed(timer, NULL)) ;
    nbi_matrix_multiply(matrix, NBI_POTENTIAL_DOUBLE,
    			&(src[0]), 2, 1.0, f, 1, 0.0, work) ;
    fprintf(stderr, "%s: evaluating single-layer potential; t=%lg\n",
    	    progname, t = g_timer_elapsed(timer, NULL)) ;
    nbi_matrix_multiply(matrix, NBI_POTENTIAL_SINGLE,
     			&(src[1]), 2, -2.0, f, 1, 2.0, work) ;
    fprintf(stderr, "%s: surface integration complete; t=%lg (%lg)\n",
	    progname,
	    g_timer_elapsed(timer, NULL), g_timer_elapsed(timer, NULL) - t) ;
    
    emax = 0.0 ; fmax = 0.0 ;
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      xp = (NBI_REAL *)nbi_surface_node(s,i) ;
      G = greens_function_laplace(xp, xs) ;
      fmax = MAX(fmax, G) ;
      emax = MAX(emax, fabs(G - f[i])) ;
      fprintf(output, "%lg %lg %lg %lg %lg\n",
	      xp[0], xp[1], xp[2], f[i], G) ; /* fabs(G - f[i])) ; */
    }
    
    fprintf(stderr, "L_inf norm: %lg\n", emax/fmax) ;

    return 0 ;
  }
  
  return 0 ;
}
