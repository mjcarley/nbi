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

gint NBI_FUNCTION_NAME(nbi_gmres_real)(nbi_matrix_t *A, gint n,
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
  gint iter, one = 1, i, k, str ;
  gdouble bnrm2, al, bt, norm, tmp ;
  gdouble *r, *V, *H, *cs, *sn, *s, *w, *y, *mwork ;

  /* n = nbi_surface_node_number(A->s) ; */
  
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

gint NBI_FUNCTION_NAME(nbi_gmres_complex)(nbi_matrix_t *A,
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

  this is strictly a real solver, with modifications to treat the
  complex matrix and vectors as generating real outputs

  rhs and solution are assumed densely packed
*/
  
{
  gint iter, one = 1, two = 2, i, k, str, n, j ;
  gdouble bnrm2, al, bt, norm, tmp ;
  gdouble *r, *V, *H, *cs, *sn, *s, *w, *y, *mwork, *buf ;

  /*size of problem as number of entries in vectors and matrices*/
  n = 2*nbi_surface_node_number(A->s) ;
  
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
  buf = &(y[m]) ;
  /* mwork = &(y[m]) ; */
  mwork = &(buf[n]) ;
  
  /* bnrm2 = blaswrap_dnrm2(n,b,bstr) ; */
  bnrm2 = blaswrap_dnrm2(n,b,one) ;
  if ( bnrm2 == 0.0 ) bnrm2 = 1.0 ;

  al = -1.0 ; bt = 1.0 ;
  /* blaswrap_dcopy(n, b, bstr, r, one) ; */
  blaswrap_dcopy(n, b, one, r, one) ;
  
  *error = blaswrap_dnrm2(n, r, one)/bnrm2 ;

  for ( iter = 0 ; iter < max_it ; iter ++ ) {
    fprintf(stderr, "%s: iteration %d\n", __FUNCTION__, iter) ;
    al = -1.0 ; bt = 1.0 ;
    /* blaswrap_dcopy(n, b, bstr, r, one) ; */
    blaswrap_dcopy(n, b, one, r, one) ;
    /* nbi_matrix_multiply(A, x, xstr, al, r, one, bt, nthreads, mwork) ; */
    nbi_matrix_multiply(A, x, xstr, al, r, two, bt, nthreads, mwork) ;

    str = m+1 ;
    blaswrap_dcopy(n, r, one, &(V[0]), str) ;
    norm = blaswrap_dnrm2(n,r,one) ;

    memset(s, 0, m*sizeof(gdouble)) ;
    s[0] = norm ;
    
    norm = 1.0/norm ;
    blaswrap_dscal(n, norm, &(V[0]), str) ;

    for ( i = 0 ; i < m ; i ++ ) {
      al = 1.0 ; bt = 0.0 ; str = m+1 ;
      /* nbi_matrix_multiply(A, &(V[i]), str, al, w, one, bt, nthreads, mwork) ; */
      for ( j = 0 ; j < n ; j ++ ) {
	buf[j] = V[i+str*j] ;
      }
      /* nbi_matrix_multiply(A, &(V[i]), str, al, w, two, bt, nthreads, mwork) ; */
      nbi_matrix_multiply(A, buf, two, al, w, two, bt, nthreads, mwork) ;
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
	/* gmres_update(x, xstr, i+1, H, n, m, s, V, y) ; */
	gmres_update(x, one, i+1, H, n, m, s, V, y) ;
	/* fprintf(stderr, "break i=%d/%d; iter=%d\n", i, m, iter) ; */
	break ;
      }
    }

    if ( *error <= tol ) break ;

    /* gmres_update(x, xstr, m, H, n, m, s, V, y) ; */
    gmres_update(x, one, m, H, n, m, s, V, y) ;
    al = -1.0 ; bt = 1.0 ;
    /* blaswrap_dcopy(n, b, bstr, r, one) ; */
    blaswrap_dcopy(n, b, one, r, one) ;
    /* nbi_matrix_multiply(A, x, xstr, al, r, one, bt, nthreads, mwork) ; */
    nbi_matrix_multiply(A, x, xstr, al, r, two, bt, nthreads, mwork) ;

    s[i+1] = blaswrap_dnrm2(n, r, one) ;
    *error = s[i+1]/bnrm2 ;
    if ( *error <= tol ) break ;
  }

  if ( *error > tol ) return -iter ;

  return iter ;
}

gint nbi_gmres_workspace_size_complex(gint n, gint m)

/*
 * n: problem size (number of nodes)
 * m: GMRES restart interval
 */
  
{
  gint size ;

  g_assert(n > 0) ; g_assert(m > 0) ;

  size =
    2*n +
    2*n*(m+1) +
    (m+1)*m +
    m +
    m +
    2*n +
    m + 1 +
    2*n ;
  
  return size ;
}

gint nbi_gmres_workspace_size_real(gint n, gint m)

/*
 * n: problem size (number of nodes)
 * m: GMRES restart interval
 */
  
{
  gint size ;

  g_assert(n > 0) ; g_assert(m > 0) ;

  size =
    n +
    n*(m+1) +
    (m+1)*m +
    m +
    m +
    n +
    m + 1 ;
  
  return size ;
}
