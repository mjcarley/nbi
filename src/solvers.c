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

#include <blaswrap.h>

#include <nbi.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#ifdef HAVE_PETSC
#include <petscksp.h>
#endif

#include "nbi-private.h"

static void gmres_update(gdouble *x, gint xstr, gint k, gdouble *H,
			 gint n, gint m,
			 gdouble *s, gdouble *V, gdouble *y)

{
  gint one = 1, ldv = m+1, info ;
  gdouble al=1.0, bt=1.0 ;
  
  blaswrap_dcopy(m, s, one, y, one) ;
  blaswrap_dtrtrs(TRUE, FALSE, FALSE, k, one, H, m, y, ldv, &info) ;
  blaswrap_dgemv(FALSE, n, k, al, V, ldv, y, one, bt, x, one) ;
  
  return ;
}

gint NBI_FUNCTION_NAME(nbi_gmres_real)(nbi_matrix_t *A, 
				       NBI_REAL *x, gint xstr,
				       NBI_REAL *b, gint bstr,
				       gint m, gint max_it,
				       NBI_REAL tol, NBI_REAL *error,
				       gint nthreads,
				       NBI_REAL *work)

/*
  A is n x n matrix, vectors are appropriate sizes

  workspace size: 2n + 4m + (m+1)*(n+m) + 2

  m is restart interval
*/
  
{
  gint iter, one = 1, i, k, str, n ;
  gdouble bnrm2, al, bt, norm, tmp ;
  gdouble *r, *V, *H, *cs, *sn, *s, *w, *y, *mwork ;

  n = nbi_surface_node_number(A->s) ;
  
  iter = 0 ;
  r = &(work[0]) ;
  V = &(r[n]) ;
  memset(V, 0, (m+1)*n*sizeof(gdouble)) ;
  H = &(V[n*(m+1)]) ;
  memset(H, 0, (m+1)*m*sizeof(gdouble)) ;
  cs = &(H[(m+1)*m]) ;
  sn = &(cs[m]) ;
  w = &(sn[m]) ;
  s = &(w[n]) ;
  y = &(s[m+1]) ;
  mwork = &(y[m]) ;
  
  bnrm2 = blaswrap_dnrm2(n,b,bstr) ;
  if ( bnrm2 == 0.0 ) bnrm2 = 1.0 ;

  al = -1.0 ; bt = 1.0 ;
  blaswrap_dcopy(n, b, bstr, r, one) ;
  nbi_matrix_multiply(A, x, xstr, al, r, one, bt, nthreads, mwork) ;
  
  *error = blaswrap_dnrm2(n, r, one)/bnrm2 ;

  for ( iter = 0 ; iter < max_it ; iter ++ ) {
    fprintf(stderr, "%s: iteration %d\n", __FUNCTION__, iter) ;
    al = -1.0 ; bt = 1.0 ;
    blaswrap_dcopy(n, b, bstr, r, one) ;
    nbi_matrix_multiply(A, x, xstr, al, r, one, bt, nthreads, mwork) ;

    str = m+1 ;
    blaswrap_dcopy(n, r, one, &(V[0]), str) ;
    norm = blaswrap_dnrm2(n,r,one) ;

    memset(s, 0, m*sizeof(gdouble)) ;
    s[0] = norm ;
    
    norm = 1.0/norm ;
    blaswrap_dscal(n, norm, &(V[0]), str) ;

    for ( i = 0 ; i < m ; i ++ ) {
      al = 1.0 ; bt = 0.0 ; str = m+1 ;
      nbi_matrix_multiply(A, &(V[i]), str, al, w, one, bt, nthreads, mwork) ;
      for ( k = 0 ; k <= i ; k ++ ) {
	H[k*m+i] = blaswrap_ddot(n, w, one, &(V[k]), str) ;
	bt = -H[k*m+i] ; str = m+1 ;
	blaswrap_daxpy(n, bt, &(V[k]), str, w, one) ;
      }

      H[(i+1)*m+i] = blaswrap_dnrm2(n, w, one) ;
      str = m + 1 ;
      blaswrap_dcopy(n, w, one, &(V[i+1]), str) ;
      norm = 1.0/H[(i+1)*m+i] ;
      blaswrap_dscal(n, norm, &(V[i+1]), str) ;

      for ( k = 0 ; k <= i-1 ; k ++ ) {
	tmp = cs[k]*H[k*m+i] + sn[k]*H[(k+1)*m+i] ;
	H[(k+1)*m+i] = -sn[k]*H[k*m+i] + cs[k]*H[(k+1)*m+i] ;
	H[k*m+i] = tmp ;
      }
      blaswrap_drotg(H[i*m+i],H[(i+1)*m+i],&(cs[i]),&(sn[i])) ;
      tmp = cs[i]*s[i] ;
      s[i+1] = -sn[i]*s[i] ;
      s[i] = tmp ;
      H[i*m+i] = cs[i]*H[i*m+i] + sn[i]*H[(i+1)*m+i] ;
      H[(i+1)*m+i] = 0.0 ;
      *error = fabs(s[i+1])/bnrm2 ;
      fprintf(stderr, "%s: substep %d; error: %lg\n",
	      __FUNCTION__, i, *error) ;
      if ( *error <= tol ) {
	gmres_update(x, xstr, i+1, H, n, m, s, V, y) ;
	/* fprintf(stderr, "break i=%d/%d; iter=%d\n", i, m, iter) ; */
	break ;
      }
    }

    if ( *error <= tol ) break ;

    gmres_update(x, xstr, m, H, n, m, s, V, y) ;
    al = -1.0 ; bt = 1.0 ;
    blaswrap_dcopy(n, b, bstr, r, one) ;
    nbi_matrix_multiply(A, x, xstr, al, r, one, bt, nthreads, mwork) ;

    s[i+1] = blaswrap_dnrm2(n, r, one) ;
    *error = s[i+1]/bnrm2 ;
    if ( *error <= tol ) break ;
  }

  if ( *error > tol ) return -iter ;

  return iter ;
}

static void gmres_update_c(gdouble *x, gint xstr, gint k, gdouble *H,
			   gint n, gint m,
			   gdouble *s, gdouble *V, gdouble *y)

{
  gint one = 1, ldv = m+1, info ;
  gdouble al[]={1.0, 0}, bt[]={1.0,0} ;
  
  blaswrap_zcopy(m, s, one, y, one) ;
  blaswrap_ztrtrs(TRUE, FALSE, FALSE, k, one, H, m, y, ldv, &info) ;
  blaswrap_zgemv(FALSE, n, k, al, V, ldv, y, one, bt, x, one) ;
  
  return ;
}

static void local_zrotg(NBI_REAL *H, NBI_REAL *Hp1,
			NBI_REAL *c, NBI_REAL *s)

/*
 * From Saad, Iterative Methods for Sparse Linear Systems, eq. 6.81
 */
  
{
  NBI_REAL rt ;

  rt = sqrt(H[0]*H[0] + H[1]*H[1] + Hp1[0]*Hp1[0]) ;

  s[0] = Hp1[0]/rt ; s[1] = 0.0 ;
  c[0] = H[0]/rt ; c[1] = H[1]/rt ;
  
  return ;
}

gint NBI_FUNCTION_NAME(nbi_gmres_complex)(nbi_matrix_t *Ai, 
					  NBI_REAL *x, gint xstr,
					  NBI_REAL *b, gint bstr,
					  gint m, gint max_it,
					  NBI_REAL tol, NBI_REAL *error,
					  gint nthreads,
					  NBI_REAL *work)

/*
  A is n x n matrix, vectors are appropriate sizes

  workspace size: 2*(2n + 4m + (m+1)*(n+m) + 2)

  m is restart interval
*/
  
{
  gint iter, one = 1, two = 2, i, k, str, n ;
  gdouble bnrm2, al, bt, norm, tmp, ac[2], tc[2] ;
  gdouble *r, *V, *H, *cs, *sn, *s, *w, *y, *mwork ;
  gint lda, ldb ;
  gdouble *A ;

  A = (gdouble *)Ai ;

  n = 2 ;
  
  g_assert(xstr == 2) ;
  g_assert(bstr == 2) ;
  
  /* n  = nbi_surface_node_number(A->s) ; */
  memset(work, 0, 2*(2*n + 4*m + (m+1)*(n+m) + 2)*sizeof(gdouble)) ;
  
  iter = 0 ;
  r = &(work[0]) ;
  V = &(r[2*n]) ;
  memset(V, 0, 2*(m+1)*n*sizeof(gdouble)) ;
  H = &(V[2*n*(m+1)]) ;
  memset(H, 0, 2*(m+1)*m*sizeof(gdouble)) ;
  cs = &(H[2*(m+1)*m]) ;
  sn = &(cs[2*(m+1)]) ;
  w = &(sn[2*(m+1)]) ;
  s = &(w[2*(n+1)]) ;
  y = &(s[2*(m+1)]) ;
  mwork = &(y[2*m]) ;
  
  /* bnrm2 = blaswrap_dnrm2(n2,b,bstr) ; */
  /* if ( bnrm2 == 0.0 ) bnrm2 = 1.0 ; */

  blaswrap_zcopy(n, b, one, r, one) ;
  bnrm2 = blaswrap_dznrm2(n, b, one) ;
  if ( bnrm2 == 0.0 ) bnrm2 = 1.0 ;

  al = -1.0 ; bt = 1.0 ;
  lda = n ;
  ac[0] = al ; ac[1] = 0 ;
  tc[0] = bt ; tc[1] = 0 ;
  blaswrap_zgemv(FALSE, n, n, ac, A, lda, x, one, tc, r, one) ;
  /* nbi_matrix_multiply(A, x, xstr, al, r, 2, bt, nthreads, mwork) ; */
  
  *error = blaswrap_dznrm2(n, r, one)/bnrm2 ;

  for ( iter = 0 ; iter < max_it ; iter ++ ) {
    fprintf(stderr, "%s: iteration %d\n", __FUNCTION__, iter) ;
    al = -1.0 ; bt = 1.0 ;
    blaswrap_zcopy(n, b, one, r, one) ;
    ac[0] = al ; ac[1] = 0 ;
    tc[0] = bt ; tc[1] = 0 ;
    blaswrap_zgemv(FALSE, n, n, ac, A, lda, x, one, tc, r, one) ;
    /* nbi_matrix_multiply(A, x, xstr, al, r, 2, bt, nthreads, mwork) ; */

    str = m+1 ;
    blaswrap_zcopy(n, r, one, &(V[0]), str) ;
    norm = blaswrap_dznrm2(n,r,one) ;

    memset(s, 0, 2*m*sizeof(gdouble)) ;
    s[2*0+0] = norm ;
    
    /* norm = 1.0/norm ; */
    ac[0] = 1.0/norm ; ac[1] = 0.0 ;
    blaswrap_zscal(n, ac, &(V[0]), str) ;

    for ( i = 0 ; i < m ; i ++ ) {
      al = 1.0 ; bt = 0.0 ; str = m+1 ;
      ac[0] = al ; ac[1] = 0 ;
      tc[0] = bt ; tc[1] = 0 ;
      blaswrap_zgemv(FALSE, n, n, ac, A, lda, &(V[2*i]), str, tc, w, one) ;
      /* nbi_matrix_multiply(A, &(V[2*i]), 2*str, al, w, 2, bt, nthreads, mwork) ; */
      for ( k = 0 ; k <= i ; k ++ ) {
	/* H[k*m+i] = blaswrap_ddot(n, w, one, &(V[k]), str) ; */
	blaswrap_zdotu(&(H[2*(k*m+i)]), n, w, one, &(V[2*k]), str) ;
	bt = -H[k*m+i] ; str = m+1 ;
	/* blaswrap_daxpy(n2, bt, &(V[k]), str, w, one) ; */
	ac[0] = -H[2*(k*m+i)+0] ; ac[1] = -H[2*(k*m+i)+1] ;
	blaswrap_zaxpy(n, ac, &(V[2*k]), str, w, one) ;
      }

      H[2*((i+1)*m+i)+0] = blaswrap_dznrm2(n, w, one) ;
      H[2*((i+1)*m+i)+1] = 0 ;
      str = m + 1 ;
      /* blaswrap_dcopy(n2, w, one, &(V[i+1]), str) ; */
      blaswrap_zcopy(n, w, one, &(V[2*(i+1)]), str) ;
      ac[0] = 1.0/H[2*((i+1)*m+i)+0] ; ac[1] = 0.0 ;
      blaswrap_zscal(n, ac, &(V[2*(i+1)]), str) ;

      for ( k = 0 ; k <= i-1 ; k ++ ) {
	/* tmp = cs[k]*H[k*m+i] + sn[k]*H[(k+1)*m+i] ; */
	/* H[(k+1)*m+i] = -sn[k]*H[k*m+i] + cs[k]*H[(k+1)*m+i] ; */
	/* H[k*m+i] = tmp ; */
	tc[0] =
	  cs[2*k+0]*H[2*(k*m+i)+0] - cs[2*k+1]*H[2*(k*m+i)+1] +
	  sn[2*k+0]*H[2*((k+1)*m+i)+0] - sn[2*k+1]*H[2*((k+1)*m+i)+1] ;
	tc[1] =
	  cs[2*k+1]*H[2*(k*m+i)+0]     + cs[2*k+0]*H[2*(k*m+i)+1] +
	  sn[2*k+1]*H[2*((k+1)*m+i)+0] + sn[2*k+0]*H[2*((k+1)*m+i)+1] ;

	ac[0] =
	  cs[2*k+0]*H[2*((k+1)*m+i)+0] - cs[2*k+1]*H[2*((k+1)*m+i)+1] -
	  sn[2*k+0]*H[2*(k*m+i)+0] + sn[2*k+1]*H[2*(k*m+i)+1] ;
	ac[1] =
	  cs[2*k+1]*H[2*((k+1)*m+i)+0] + cs[2*k+0]*H[2*((k+1)*m+i)+1] -
	  sn[2*k+1]*H[2*(k*m+i)+0] - sn[2*k+0]*H[2*(k*m+i)+1] ;
	
	H[2*(k*m+i)+0] = tc[0] ; H[2*(k*m+i)+1] = tc[1] ; 
	H[2*((k+1)*m+i)+0] = ac[0] ; H[2*((k+1)*m+i)+1] = ac[1] ; 
      }
      /* blaswrap_drotg(H[i*m+i],H[(i+1)*m+i],&(cs[i]),&(sn[i])) ; */
      local_zrotg(&(H[2*(i*m+i)]),&(H[2*((i+1)*m+i)]),
		  &(cs[2*i]),&(sn[2*i])) ;
      /* tmp = cs[i]*s[i] ; */
      tc[0] = cs[2*i+0]*s[2*i+0] - cs[2*i+1]*s[2*i+1] ;
      tc[1] = cs[2*i+1]*s[2*i+0] + cs[2*i+0]*s[2*i+1] ;
      /* s[i+1] = -sn[i]*s[i] ; */
      s[2*(i+1)+0] = -(sn[2*i+0]*s[2*i+0] - sn[2*i+1]*s[2*i+1]) ;
      s[2*(i+1)+1] = -(sn[2*i+1]*s[2*i+0] + sn[2*i+0]*s[2*i+1]) ;
      s[2*i+0] = tc[0] ; s[2*i+1] = tc[1] ;
      /* H[i*m+i] = cs[i]*H[i*m+i] + sn[i]*H[(i+1)*m+i] ; */
      /* H[(i+1)*m+i] = 0.0 ; */
      tc[0] = 
	cs[2*i+0]*H[2*(i*m+i)+0]     - cs[2*i+1]*H[2*(i*m+i)+1] +
	sn[2*i+0]*H[2*((i+1)*m+i)+0] - sn[2*i+1]*H[2*((i+1)*m+i)+1] ;
      tc[1] = 
	cs[2*i+1]*H[2*(i*m+i)+0]     + cs[2*i+0]*H[2*(i*m+i)+1] +
	sn[2*i+1]*H[2*((i+1)*m+i)+0] + sn[2*i+0]*H[2*((i+1)*m+i)+1] ;
      H[2*(i*m+i)+0] = tc[0] ; 
      H[2*(i*m+i)+1] = tc[1] ; 

      H[2*((i+1)*m+i)+0] = 0.0 ;
      H[2*((i+1)*m+i)+1] = 0.0 ;
      /* *error = fabs(s[i+1])/bnrm2 ; */
      *error = sqrt(s[2*(i+1)+0]*s[2*(i+1)+0] + s[2*(i+1)+1]*s[2*(i+1)+1]) ;
      *error /= bnrm2 ;
      fprintf(stderr, "%s: substep %d; error: %lg\n",
	      __FUNCTION__, i, *error) ;
      if ( *error <= tol ) {
	gmres_update_c(x, xstr, i+1, H, n, m, s, V, y) ;
	/* fprintf(stderr, "break i=%d/%d; iter=%d\n", i, m, iter) ; */
	break ;
      }
    }

    /* if ( *error <= tol ) break ; */

    gmres_update_c(x, xstr, m, H, n, m, s, V, y) ;
    al = -1.0 ; bt = 1.0 ;
    blaswrap_zcopy(n, b, one, r, one) ;
    /* blaswrap_dcopy(n, b, bstr, r, one) ; */
    /* blaswrap_dgemv(FALSE, n, n, al, A, lda, x, xstr, bt, r, one) ; */
    /* nbi_matrix_multiply(A, x, xstr, al, r, two, bt, nthreads, mwork) ; */
    ac[0] = al ; ac[1] = 0 ;
    tc[0] = bt ; tc[1] = 0 ;
    blaswrap_zgemv(FALSE, n, n, ac, A, lda, x, one, tc, r, one) ;

    s[2*(i+1)+0] = blaswrap_dznrm2(n, r, one) ;
    /* *error = s[i+1]/bnrm2 ; */
    *error = sqrt(s[2*(i+1)+0]*s[2*(i+1)+0] + s[2*(i+1)+1]*s[2*(i+1)+1]) ;
    *error /= bnrm2 ;
    if ( *error <= tol ) break ;
  }

  if ( *error > tol ) return -iter ;

  return iter ;
}
