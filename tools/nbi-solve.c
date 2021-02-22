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

#include "nbi-private.h"


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
    G = q*0.25*M_1_PI/R*nbi_surface_node_weight(s,i) ;
    p [i*pstr] += G ;
    pn[i*nstr] -= G*Rn/R ;    
  }
  
  return 0 ;
}


static gint source_field_laplace(nbi_surface_t *s,
				 gdouble *p , gint pstr,
				 gdouble *pn, gint nstr,
				 gdouble *x, gdouble *f)

{
  gint i ;
  gdouble G, dG, *y, *n, r[3], R, Rn ;
  
  *f = 0 ;

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    y = nbi_surface_node(s, i) ;
    n = nbi_surface_normal(s, i) ;

    nbi_vector_diff(r, x, y) ;
    R = nbi_vector_length(r) ;
    Rn = nbi_vector_scalar(r,n)/R ;
    G  = 0.25*M_1_PI/R ;
    dG = G*Rn/R ;

    *f += p[i*pstr]*dG - pn[i*nstr]*G ;
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

static gboolean patch_encroaches(gdouble *c, gdouble r,
				 gdouble *xp, gint xstr, gint np)
{
  gint i ;

  for ( i = 0 ; i < np ; i ++ ) {
    if ( nbi_vector_distance2(c, &(xp[i*xstr])) < r*r )
      return TRUE ;
  }
  
  return FALSE ;
}

static gint matrix_scale_weights(gdouble *A, gint nr, gint nc,
				 gdouble *w, gint wstr)

{
  gint i, j ;

  for ( i = 0 ; i < nr ; i ++ ) {
    for ( j = 0 ; j < nc ; j ++ ) {
      A[i*2*nc   +j] /= w[j*wstr] ;
      A[i*2*nc+nc+j] /= w[j*wstr] ; 
    }
  }
  
  return 0 ;
}

static gint source_target_correction(gdouble *xt, gint xstr, gint nt,
				     gdouble *xs, gint sstr, gint ns,
				     gdouble *A)

{
  gint i, j ;
  gdouble G, dG, R, Rn, r[3] ;
  
  for ( i = 0 ; i < nt ; i ++ ) {
    for ( j = 0 ; j < ns ; j ++ ) {
      nbi_vector_diff(r, &(xt[i*xstr]), &(xs[j*sstr])) ;
      R = nbi_vector_length(r) ;
      if ( R > 1e-12 ) {
	Rn = nbi_vector_scalar(r,&(xs[j*sstr+3]))/R ;
	G  = 0.25*M_1_PI/R ;
	dG = G*Rn/R ;
	A[i*2*ns +      j] -=  G ;
	A[i*2*ns + ns + j] -= dG ;
      }
    }
  }
  
  return 0 ;
}

static gint point_source_field_laplace(nbi_surface_t *s,
				       gdouble *p , gint pstr,
				       gdouble *pn, gint nstr,
				       gdouble *x, gdouble *f)
  
{
  gint i ;  
  gdouble R, r[3] ;

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    nbi_vector_diff(r, x, nbi_surface_node(s, i)) ;
    R = nbi_vector_length(r) ;
    if ( R > 1e-12 ) {
      *f += (p[i*pstr]*nbi_vector_scalar(r,nbi_surface_normal(s,i))/R/R -
	     pn[i*nstr])*0.25*M_1_PI/R ;
    }
  }
  
  return 0 ;
}

gint nbi_surface_integrate_self(nbi_surface_t *s, gdouble eta,
				gint nq, gint dmax, gdouble tol, gint N,
				gdouble *p , gint pstr,
				gdouble *pn, gint nstr,
				gint pt, gdouble *f)

{
  gint i, j, order, nst0, nst, xstr ;
  gint lda, one = 1 ;
  gdouble c[3], r, NK0, K0[453*453], *q, *st, Ast[2*453*453] ;
  gdouble *xs, *xt, al, bt, *work ;

  xstr = NBI_SURFACE_NODE_LENGTH ;
  
  sqt_quadrature_select(nq, &q, &order) ;

  g_assert(pt < nbi_surface_patch_node_number(s, pt)) ;
  
  nst0 = nbi_surface_patch_node_number(s, pt) ;

  tol /= nst0*nst0 ;
  
  memset(f, 0, nst0*sizeof(gdouble)) ;

  work = (gdouble *)g_malloc(4*dmax*2*nst0*nst0*sizeof(gdouble)) ;
  
  sqt_quadrature_select(nst0, &st, &order) ;

  NK0 = sqt_koornwinder_interp_matrix(&(st[0]), 3, &(st[1]), 3, &(st[2]), 3,
				      nst0, K0) ;
  
  i = nbi_surface_patch_node(s, pt) ;
  xt = nbi_surface_node(s,i) ; xstr = NBI_SURFACE_NODE_LENGTH ;
  nbi_surface_patch_centroid(nbi_surface_node(s,i),
			     NBI_SURFACE_NODE_LENGTH,
			     &(nbi_surface_node_weight(s, i)),
			     NBI_SURFACE_NODE_LENGTH,
			     nst0, c) ;
  r = nbi_surface_patch_radius(nbi_surface_node(s,i),
			       NBI_SURFACE_NODE_LENGTH, nst0, c) ;

  for ( i = 0 ; i < nst0 ; i ++ ) {
    point_source_field_laplace(s, p, pstr, pn, nstr, &(xt[i*xstr]), &(f[i])) ;
  }

  for ( i = 0 ; i < nbi_surface_patch_number(s) ; i ++ ) {
    /*near-field and self terms with point source terms subtracted*/
    nst = nbi_surface_patch_node_number(s, i) ;
    lda = 2*nst ;
    j = nbi_surface_patch_node(s, i) ;
    g_assert(nst == nst0) ;
    xs = nbi_surface_node(s, j) ;    
    if ( patch_encroaches(c, eta*r, xs, xstr, nst) ) {
      if ( i != pt ) {
	sqt_laplace_source_target_kw_adaptive(xs, xstr, nst,
					      q, nq, K0, NK0,
					      tol, dmax,
					      xt, xstr, nst0,
					      Ast, work) ;
	matrix_scale_weights(Ast, nst0, nst,
			     &(nbi_surface_node_weight(s, j)),
			     NBI_SURFACE_NODE_LENGTH) ;
      } else {
	sqt_laplace_source_target_kw_self(xt, xstr, nst0,
					  K0, NK0, N,
					  &(st[0]), 3, &(st[1]), 3,
					  Ast) ;
	matrix_scale_weights(Ast, nst0, nst,
			     &(nbi_surface_node_weight(s, j)),
			     NBI_SURFACE_NODE_LENGTH) ;
      }
      source_target_correction(xt, xstr, nst0, xs, xstr, nst, Ast) ;
      al =  1.0 ; bt = 1.0 ;
      blaswrap_dgemv(FALSE, nst0, nst, al, &(Ast[1*nst]), lda,
		     &(p[j*pstr]) , pstr, bt, f, one) ;
      al = -1.0 ; bt = 1.0 ;
      blaswrap_dgemv(FALSE, nst0, nst, al, &(Ast[0*nst]), lda,
		     &(pn[j*nstr]), nstr, bt, f, one) ;      
    }
  }
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  nbi_surface_t *s ;
  gint nth, nph, nq, i, np, pt, nqa, dmax, N ;
  gdouble r, S, x[3], xs[512], f[453], eta, *xp, *src, tol ;
  FILE *output ;
  gchar ch ;
  
  output = stdout ;
  
  r = 1.0 ; nth = 16 ; nph = 8 ; nq = 25 ;
  eta = 1.25 ; dmax = 8 ; tol = 1e-12 ; N = 8 ; nqa = 54 ;
  pt = 0 ;
  
  while ( (ch = getopt(argc, argv, "a:d:e:l:N:n:p:q:r:t:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'a': nqa  = atoi(optarg) ; break ;      
    case 'd': dmax = atoi(optarg) ; break ;
    case 'e': tol  = atof(optarg) ; break ;      
    case 'l': pt   = atoi(optarg) ; break ;      
    case 'N': N    = atoi(optarg) ; break ;
    case 'n': eta  = atof(optarg) ; break ;      
    case 'p': nph  = atoi(optarg) ; break ;
    case 'q': nq   = atoi(optarg) ; break ;
    case 'r': r    = atof(optarg) ; break ;
    case 't': nth  = atoi(optarg) ; break ;      
    }
  }
  
  np = 2*nth*nph ;
  fprintf(stderr, "nth = %d; nph = %d; nq = %d\n", nth, nph, nq) ;
  fprintf(stderr, "r = %lg\n", r) ;
  fprintf(stderr, "allocating surface with %d nodes, %d patches\n", np*nq, np) ;
  
  s = nbi_surface_alloc(np*nq, np) ;
  
  /* nbi_geometry_sphere(s, r, nth, nph, nq) ;   */
  nbi_geometry_ellipsoid(s, r, 1.0, 1.0, nth, nph, nq) ;  

  fprintf(stderr, "spherical surface with %d nodes, %d patches generated\n",
	  nbi_surface_node_number(s), nbi_surface_patch_number(s)) ;

  S = 0.0 ;
  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ )
    S += nbi_surface_node_weight(s, i) ;
  
  fprintf(stderr, "surface area: %lg (%lg, %lg)\n",
	  S, 4.0*M_PI*r*r, fabs(4.0*M_PI*r*r - S)) ;
  
  /* nbi_surface_write(s, output) ; */

  /*do a radiation calculation to test surface generation and
    integration*/

  /*boundary point sources*/
  src = (gdouble *)g_malloc0(nbi_surface_node_number(s)*2*sizeof(gdouble)) ;

  xs[0] = 0.3 ; xs[1] = -0.4 ; xs[2] = 0.2 ;

  make_sources(s, xs, 1.0, &(src[0]), 2, &(src[1]), 2) ;
  
  nbi_surface_integrate_self(s, eta, nqa, dmax, tol, N,
			     &(src[0]), 2, &(src[1]), 2, pt, f) ;

  for ( i = 0 ; i < nbi_surface_patch_node_number(s, pt) ; i ++ ) {
    xp = nbi_surface_patch_local_node(s,pt,i) ;
    fprintf(stdout, "%lg %lg %lg %lg\n", xp[0], xp[1], xp[2], f[i]) ;
  }

  return 0 ;
  
  x[0] = 1.5 ; x[1] = 0.5 ; x[2] = -3.9 ;
  /* x[0] = 0.0 ; x[1] = 0.0 ; x[2] = 0.0 ; */
  x[0] = 3.0 ; x[1] = -1.5 ; x[2] = 2.0 ;

  x[0] = r*1.05 ; x[1] = r*1.3 ; x[2] = -r*0.7 ;
  
  source_field_laplace(s, &(src[0]), 2, &(src[1]), 2, x, f) ;

  fprintf(stdout, "%lg %lg %lg %lg %lg (%lg)\n",
  	  x[0], x[1], x[2], f[0], greens_function_laplace(x, xs),
  	  fabs(f[0] - greens_function_laplace(x, xs))) ;
  
  return 0 ;
}
