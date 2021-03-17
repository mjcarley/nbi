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

GTimer *timer ;
gchar *progname ;


static gint matrix_scale_weights(gdouble *A, gint nr, gint nc,
				 gdouble *w, gint wstr)

{
  gint i, j ;
  /* gdouble al ; */

  /*this seems to be the quickest approach (divide by weight, don't
    use BLAS)*/
  for ( j = 0 ; j < nc ; j ++ ) {
    /* al = 1.0/w[wstr*j] ; */
    for ( i = 0 ; i < nr ; i ++ ) {
      /* A[i*2*nc   +j] *= al ; A[i*2*nc+nc+j] *= al ; */
      A[i*2*nc   +j] /= w[wstr*j] ; A[i*2*nc+nc+j] /= w[wstr*j] ;
    }
  }
  
  /* gint n = 2*nc ; */
  /* for ( j = 0 ; j < nc ; j ++ ) { */
  /*   al = 1.0/w[j*wstr] ;   */
  /*   blaswrap_dscal(nr, al, &(A[0*nc+j]), n) ; */
  /*   blaswrap_dscal(nr, al, &(A[1*nc+j]), n) ; */
  /* } */
  
  return 0 ;
}

static gint source_target_correction_indexed(gdouble *xt, gint xstr,
					     gint *idx, gint nt,
					     gdouble *xs, gint sstr, gint ns,
					     gdouble *A)

{
  gint i, j, k ;
  gdouble G, dG, R, Rn, r[3] ;
  
  for ( k = 0 ; k < nt ; k ++ ) {
    i = idx[k] ;
    for ( j = 0 ; j < ns ; j ++ ) {
      nbi_vector_diff(r, &(xt[i*xstr]), &(xs[j*sstr])) ;
      R = nbi_vector_length(r) ;
      if ( R > 1e-12 ) {
	Rn = nbi_vector_scalar(r,&(xs[j*sstr+3]))/R ;
	G  = 0.25*M_1_PI/R ;
	dG = G*Rn/R ;
	A[k*2*ns +      j] -=  G ;
	A[k*2*ns + ns + j] -= dG ;
      }
    }
  }
  
  return 0 ;
}

gint upsample_patch(gdouble *xs, gint sstr, gint ns,
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

gint nbi_surface_assemble_write(nbi_surface_t *s, gdouble eta,
				gint nq, gint dmax, gdouble tol, gint N,
				gint nu,				
				gint *idx, gint *idxp, gint nnmax,
				gint *idxu, gint nqu,
				FILE *output)

{
  gint i, j, ip, order, nsts, nst, xstr, nntot, lda, pt, *nbrs, nnbrs ;
  gdouble c[3], r, NK0, K0[453*453], Ku[453*543], *q, *st, *Ast, *xs, *work ;
  gdouble *xu ;
  
  xstr = NBI_SURFACE_NODE_LENGTH ;

  sqt_quadrature_select(nq, &q, &order) ;

  nst = nbi_surface_patch_node_number(s,0) ;
  sqt_quadrature_select(nst, &st, &order) ;

  NK0 = sqt_koornwinder_interp_matrix(&(st[0]), 3, &(st[1]), 3, &(st[2]), 3,
				      nst, K0) ;

  work = (gdouble *)g_malloc(nnmax*4*dmax*2*nst*sizeof(gdouble)) ;
  Ast  = (gdouble *)g_malloc((nnmax+nst)*2*nst*sizeof(gdouble)) ;
  xu = (gdouble *)g_malloc(nu*nbi_surface_patch_number(s)*
			   NBI_SURFACE_NODE_LENGTH*sizeof(gdouble)) ;
			   
  fprintf(stderr, "%s: starting local corrections; t=%lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL)) ;

  /*find neighbours and set up sparse matrix skeleton*/
  memset(idx, 0, nbi_surface_node_number(s)*nnmax*sizeof(gint)) ;

  nntot = 0 ; idxu[0] = 0 ;
  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    /*loop on patches treated as sources*/
    nsts = nbi_surface_patch_node_number(s, pt) ;
  
    ip = nbi_surface_patch_node(s, pt) ;
    xs = nbi_surface_node(s,ip) ; xstr = NBI_SURFACE_NODE_LENGTH ;
    nbi_surface_patch_centroid(nbi_surface_node(s,ip),
			       NBI_SURFACE_NODE_LENGTH,
			       &(nbi_surface_node_weight(s, ip)),
			       NBI_SURFACE_NODE_LENGTH,
			       nsts, c) ;
    r = nbi_surface_patch_radius(nbi_surface_node(s,ip),
				 NBI_SURFACE_NODE_LENGTH, nsts, c) ;

    /*find the neighbours*/
    nnbrs = 0 ;
    nbi_patch_neighbours(c, eta*r,
			 nbi_surface_node(s, 0), xstr,
			 nbi_surface_node_number(s),
			 0, ip,
			 &(idx[nntot]), &nnbrs, nnmax-nsts) ;
    nbi_patch_neighbours(c, eta*r,
			 nbi_surface_node(s, 0), xstr,
			 nbi_surface_node_number(s),
			 ip+nbi_surface_patch_node_number(s,pt),
			 nbi_surface_node_number(s),
			 &(idx[nntot]), &nnbrs, nnmax-nsts) ;
    g_assert(nnbrs < nnmax-nsts) ;
    for ( i = 0 ; i < nbi_surface_patch_node_number(s, pt) ; i ++ )
      idx[nntot+nnbrs+i] = ip + i ;
    nnbrs += nbi_surface_patch_node_number(s, pt) ;
    nntot += nnbrs ;
    idxp[pt+1] = nntot ;

    /*upsample into xu*/
    idxu[pt+1] = idxu[pt] +
      upsample_patch(nbi_surface_node(s,ip), xstr, nst,
		     nqu, &(xu[idxu[pt]*xstr]), xstr) ;
  }

  fprintf(output, "%d %d\n", nbi_surface_patch_number(s), xstr) ;
  for ( i = 0 ; i < nbi_surface_patch_number(s)+1 ; i ++ ) {
    fprintf(output, "%d %d\n", i, idxu[i]) ;
  }
  for ( i = 0 ; i < idxu[nbi_surface_patch_number(s)] ; i ++ ) {
    for ( j = 0 ; j < xstr ; j ++ )
      fprintf(output, " %1.16e", xu[i*xstr+j]) ;
    fprintf(output, "\n") ;
  }
  
  fprintf(output, "%d %d\n", nbi_surface_patch_number(s), nst) ;
  for ( i = 0 ; i < nbi_surface_patch_number(s)+1 ; i ++ ) {
    fprintf(output, "%d %d\n", i, idxp[i]) ;
  }

  for ( i = 0 ; i < idxp[nbi_surface_patch_number(s)] ; i ++ ) {
    fprintf(output, "%d\n", idx[i]) ;
  }

  for ( pt = 0 ; pt < nbi_surface_patch_number(s) ; pt ++ ) {
    /*local correction matrices*/
    nnbrs = idxp[pt+1] - idxp[pt] ;
    nbrs = &(idx[idxp[pt]]) ;
    nsts = nbi_surface_patch_node_number(s, pt) ;
    ip = nbi_surface_patch_node(s, pt) ;
    xs = nbi_surface_node(s,ip) ; xstr = NBI_SURFACE_NODE_LENGTH ;
    fprintf(stderr, "patch %d/%d; %d neighbours (%lg)\n",
	    pt, nbi_surface_patch_number(s), nnbrs,
	    g_timer_elapsed(timer, NULL)) ;
    sqt_laplace_source_indexed_kw_adaptive(xs, xstr, nsts, q, nq, K0, NK0,
    					   tol, dmax,
    					   nbi_surface_node(s,0), xstr,
    					   nbrs, nnbrs-nsts,
    					   Ast, work) ;
    sqt_laplace_source_target_kw_self(xs, xstr, nsts,
				      K0, NK0, N,
				      &(st[0]), 3, &(st[1]), 3,
				      &(Ast[2*(nnbrs-nsts)*nsts])) ;
    /* matrix_scale_weights(Ast, nnbrs, nsts, */
    /* 			 &(nbi_surface_node_weight(s, ip)), */
    /* 			 NBI_SURFACE_NODE_LENGTH) ; */
    /* source_target_correction_indexed(nbi_surface_node(s,0), xstr, */
    /* 				     &(idx[idxp[pt]]), nnbrs, */
    /* 				     xs, xstr, nsts, Ast) ; */
    lda = 2*nsts ;

    for ( i = 0 ; i < nnbrs ; i ++ ) {
      for ( j = 0 ; j < lda ; j ++ )
	fprintf(output, " %1.16e", Ast[i*lda+j]) ;
      fprintf(output, "\n") ;
    }
  }

  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  nbi_surface_t *s ;
  gint nth, nph, nq, np, nqa, dmax, N, *idx, *idxp, nnmax, nqu, *idxu ;
  gdouble r, eta, tol, t ;
  FILE *output ;
  gchar ch, *mfile, *gfile ;
  
  output = stdout ;
  
  r = 1.0 ; nth = 16 ; nph = 8 ; nq = 25 ; nqu = 54 ;
  eta = 1.25 ; dmax = 8 ; tol = 1e-12 ; N = 8 ; nqa = 54 ;
  nnmax = 10000 ; 
  gfile = NULL ; mfile = NULL ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  while ( (ch = getopt(argc, argv, "a:d:e:g:m:N:n:p:q:r:t:u:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'a': nqa  = atoi(optarg) ; break ;      
    case 'd': dmax = atoi(optarg) ; break ;
    case 'e': tol  = atof(optarg) ; break ;      
    case 'g': gfile = g_strdup(optarg) ; break ;
    case 'm': mfile = g_strdup(optarg) ; break ;
    case 'N': N    = atoi(optarg) ; break ;
    case 'n': eta  = atof(optarg) ; break ;      
    case 'p': nph  = atoi(optarg) ; break ;
    case 'q': nq   = atoi(optarg) ; break ;
    case 'r': r    = atof(optarg) ; break ;
    case 't': nth  = atoi(optarg) ; break ;      
    case 'u': nqu  = atoi(optarg) ; break ;      
    }
  }
  
  if ( gfile == NULL ) gfile = g_strdup("geometry.dat") ;
  if ( mfile == NULL ) mfile = g_strdup("matrix.dat") ;
  
  np = 2*nth*nph ;
  fprintf(stderr, "nth = %d; nph = %d; nq = %d\n", nth, nph, nq) ;
  fprintf(stderr, "r = %lg\n", r) ;
  fprintf(stderr, "allocating surface with %d nodes, %d patches\n", np*nq, np) ;
  
  s = nbi_surface_alloc(np*nq, np) ;

  timer = g_timer_new() ;
  
  /* nbi_geometry_sphere(s, r, nth, nph, nq) ;   */
  fprintf(stderr, "%s: initializing geometry r=%lg; t=%lg\n",
	  progname, r, g_timer_elapsed(timer, NULL)) ;
  output = fopen(gfile, "w") ;
  nbi_geometry_ellipsoid(s, r, 1.0, 1.0, nth, nph, nq) ;  
  fprintf(stderr,
	  "%s: geometry initialized, %d nodes, %d patches generated; t=%lg\n",
	  progname, nbi_surface_node_number(s), nbi_surface_patch_number(s),
	  g_timer_elapsed(timer, NULL)) ;
  nbi_surface_write(s, output) ;
  fclose(output) ;
  
  idx = (gint *)g_malloc0(nbi_surface_node_number(s)*nnmax*sizeof(gint)) ;
  idxp = (gint *)g_malloc0((nbi_surface_patch_number(s)+1)*sizeof(gint)) ;
  idxu = (gint *)g_malloc0((nbi_surface_patch_number(s)+1)*sizeof(gint)) ;
  fprintf(stderr, "%s: starting surface assembly; t=%lg\n",
	  progname, t = g_timer_elapsed(timer, NULL)) ;

  output = fopen(mfile, "w") ;  
  nbi_surface_assemble_write(s, eta, nqa, dmax, tol, N, nqu, idx, idxp, nnmax,
  			     idxu, nqu, output) ;
  fclose(output) ;
  
  return 0 ;
}
