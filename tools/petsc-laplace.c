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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

/* #undef HAVE_PETSC */

#ifdef HAVE_PETSC
#include <petscksp.h>
#include <petscoptions.h>

PetscErrorCode nbi_petsc_MatMult_real(Mat mat, Vec x, Vec y) ;
PetscErrorCode nbi_petsc_MatGetDiagonal(Mat mat, Vec v) ;
PetscErrorCode nbi_petsc_MatSetValues(Mat mat, gint ni, gint *i,
				      gint nj, gint *j, gdouble *value,
				      guint insert) ;

/* PetscErrorCode my_monitor(KSP ksp, PetscInt i, PetscReal ee, void *) */

/* { */
/*   fprintf(stdout, "iter %d: err: %lg\n", i, ee) ; */
  
/*   return 0 ; */
/* } */

#endif



/*operations for user-defined matrices in PETSc*/
#ifdef HAVE_PETSC
PetscErrorCode nbi_petsc_MatMult_real(Mat mat, Vec x, Vec y)

{
  nbi_matrix_t *m ;
  gint nthreads ;
  gdouble *work ;
  const PetscScalar *xx ;
  PetscScalar *yy ;
  gpointer *petsc_ctx ;
  
  PetscCall(MatShellGetContext(mat, &petsc_ctx)) ;  

  m    = petsc_ctx[NBI_SOLVER_DATA_MATRIX] ;
  work = petsc_ctx[NBI_SOLVER_DATA_WORK] ;
  nthreads = *((gint *)(petsc_ctx[NBI_SOLVER_DATA_NTHREADS])) ;

  PetscCall(VecGetArrayRead(x, &xx)) ;
  PetscCall(VecGetArrayWrite(y, &yy)) ;
  nbi_matrix_multiply(m, (gdouble *)xx, 1, 1.0, yy, 1, 0.0, nthreads, work) ;
  PetscCall(VecRestoreArrayRead(x, (&xx))) ;
  PetscCall(VecRestoreArrayWrite(y, (&yy))) ;

  return 0 ;
}

PetscErrorCode nbi_petsc_MatGetDiagonal(Mat mat, Vec v)

{
  nbi_matrix_t *m ;

  PetscCall(MatShellGetContext(mat, &m)) ;  
  
  if ( m->test != NULL ) {
    /*test code using an explicitly defined matrix*/
    gdouble *A, *vv, t ;
    gint i, n ;
    
    A = m->test ; n = m->nstr ;
    PetscCall(VecGetArrayRead(v, (const gdouble **)(&vv))) ;

    PetscCall(VecGetSize(v, &i)) ;
    g_assert(i == n) ;
    
    for ( i = 0 ; i < n ; i ++ ) {
      vv[i] = A[i*n+i] ;
      PetscCall(VecGetValues(v, 1, &i, &t)) ;
      g_assert(vv[i] == t) ;
    }
    
    return 0 ;
  }

  g_assert_not_reached() ;

  return 0 ;
}

PetscErrorCode nbi_petsc_MatSetValues(Mat mat, gint ni, gint *i,
				      gint nj, gint *j, gdouble *value,
				      guint insert)

{
  nbi_matrix_t *m ;
  gint n ;

  PetscCall(MatShellGetContext(mat, &m)) ;  
  
  if ( m->test != NULL ) {
    /*test code using an explicitly defined matrix*/
    gdouble *A ;

    A = m->test ; n = m->nstr ;
    A[(*i)*n+(*j)] = value[0] ;
  }  
  
  return 0 ;
}

#endif /*HAVE_PETSC*/
