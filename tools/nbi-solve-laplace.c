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

PetscErrorCode nbi_petsc_MatMult(Mat mat, Vec x, Vec y) ;
PetscErrorCode nbi_petsc_MatGetDiagonal(Mat mat, Vec v) ;
PetscErrorCode nbi_petsc_MatSetValues(Mat mat, gint ni, gint *i,
				      gint nj, gint *j, gdouble *value,
				      guint insert) ;

PetscErrorCode my_monitor(KSP ksp, PetscInt i, PetscReal ee, void *)

{
  fprintf(stdout, "iter %d: err: %lg\n", i, ee) ;
  
  return 0 ;
}


#endif

#include "nbi-private.h"

GTimer *timer ;
gchar *progname ;


static void print_help_text(FILE *f, gint depth,
			    gint order_inc, gint order_fmm,
			    gint nthreads, gdouble tol)

{
  fprintf(f, 
	  "Usage:\n\n"
	  "  %s <options>\n\n",
	  progname) ;

  fprintf(f,
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -b # boundary condition file\n"
	  "  -D # FMM tree depth (%d)\n"
	  "  -d # FMM tree increment in expansion order between levels (%d)\n"
	  "  -f use FMM\n"
	  "  -G evaluate Green's identity\n"
	  "  -g # geometry file name\n"
	  "  -L evaluate single and double layer potentials\n"
	  "  -m # matrix file name\n"
	  "  -o # FMM order (%d)\n"
	  "  -p precompute local interactions in FMM\n"
	  "  -T # number of threads (%d)\n"
	  "  -t # GMRES solution tolerance (%lg)\n",
	  depth, order_inc, order_fmm, nthreads, tol) ;
  return ;
}

gint main(gint argc, gchar **argv)

{
  nbi_surface_t *s ;
  nbi_matrix_t *matrix ;
  nbi_boundary_condition_t *bc ;
  gdouble *f, *xp, *src, t, emax, fmax, e2, f2, G, dtree, pwt, nwt ;
  FILE *output, *input ;
  gchar ch, *gfile, *mfile, *bfile ;
  gdouble *work, tol ;
  gint fmm_work_size, nqfmm, order_fmm, order_inc, i, fstr, solver_work_size ;
  gint gmres_max_iter, gmres_restart ;
  gint nthreads, nproc ;
  guint depth, order[48] = {0}, order_s, order_r, order_max ;
  gboolean fmm, shift_bw, greens_id, layer_potentials, precompute_local ;

  nthreads = 1 ;

#ifdef _OPENMP
  nproc = g_get_num_processors() ;
  nthreads = -1 ;
#else  /*_OPENMP*/
  nproc = 1 ;
#endif /*_OPENMP*/
  
  output = stdout ;
  mfile = NULL ; gfile = NULL ; bfile = NULL ;
  
  dtree = 1e-2 ; fmm = FALSE ; nqfmm = 1 ; shift_bw = TRUE ;
  greens_id = FALSE ; layer_potentials = FALSE ;
  precompute_local = FALSE ;
  
  pwt = 2.0 ; nwt = 2.0 ;
  /* pwt = 1.0 ; nwt = 1.0 ; */
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  order_fmm = 12 ; order_inc = 2 ; depth = 4 ;

  solver_work_size = 0 ;
  gmres_max_iter = 128 ; gmres_restart = 10 ; tol = 1e-9 ;
    
  fstr = 3 ;
  while ( (ch = getopt(argc, argv, "hBb:D:d:fGg:Lm:o:pr:T:t:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'h':
      print_help_text(stderr, depth, order_inc, order_fmm, nthreads, tol) ;
      return 0 ;
      break ;
    case 'B':
      fprintf(stderr, "%s: built-in boundary condition functions\n\n",
	      progname) ;
      nbi_functions_list(stderr, TRUE) ;
      return 0 ;
      break ;
    case 'b': bfile = g_strdup(optarg) ; break ;
    case 'D': depth = atoi(optarg) ; break ;
    case 'd': order_inc = atoi(optarg) ; break ;
    case 'f': fmm = TRUE ; break ;
    case 'G': greens_id = TRUE ; break ;
    case 'g': gfile = g_strdup(optarg) ; break ;
    case 'L': layer_potentials = TRUE ; break ;
    case 'm': mfile = g_strdup(optarg) ; break ;
    case 'o': order_fmm = atoi(optarg) ; break ;
    case 'p': precompute_local = TRUE ; break ;
    case 'r': gmres_restart = atoi(optarg) ; break ;
    case 'T': nthreads = atoi(optarg) ; break ;
    case 't': tol = atof(optarg) ; break ;
    }
  }

  if ( gfile == NULL ) gfile = g_strdup("geometry.dat") ;
  if ( mfile == NULL ) mfile = g_strdup("matrix.dat") ;

  fprintf(stderr, "%s: %d threads (%d processors)\n",
	  progname, nthreads, nproc) ;
  
  fprintf(stderr, "%s: reading geometry from %s\n", progname, gfile) ;
  if ( (input = fopen(gfile, "r")) == NULL ) {
    fprintf(stderr, "%s: cannot open geometry file %s\n",
	    progname, gfile) ;
    return -1 ;
  }

  s = nbi_surface_read(input) ;

  fclose(input) ;

  matrix = nbi_matrix_new(s) ;
  matrix->problem = NBI_PROBLEM_LAPLACE ;
  
  fprintf(stderr, "%s: reading matrix from %s\n", progname, mfile) ;
  if ( (input = fopen(mfile, "r")) == NULL ) {
    fprintf(stderr, "%s: cannot open matrix file %s\n",
	    progname, mfile) ;
    return -1 ;
  }

  nbi_matrix_read(input, matrix) ;
  
  fclose(input) ;

  timer = g_timer_new() ;

  fprintf(stderr,
	  "%s: geometry initialized, %d nodes, %d patches; t=%lg\n",
	  progname, nbi_surface_node_number(s), nbi_surface_patch_number(s),
	  g_timer_elapsed(timer, NULL)) ;

  /*boundary point sources*/
  fprintf(stderr, "%s: setting boundary conditions; t=%lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  src = (gdouble *)g_malloc0(nbi_surface_node_number(s)*2*sizeof(gdouble)) ;

  if ( bfile == NULL ) {
    fprintf(stderr, "%s: no boundary condition specified\n", progname) ;
    exit(1) ;
  }

  fprintf(stderr, "%s: reading boundary conditions from %s\n",
	  progname, bfile) ;
  if ( (input = fopen(bfile, "r")) == NULL ) {
    fprintf(stderr, "%s: cannot open boundary condition file %s\n",
	    progname, bfile) ;
    return 1 ;
  }

  bc = nbi_boundary_condition_new(NBI_PROBLEM_LAPLACE) ;
  nbi_boundary_condition_read(input, bc) ;
  
  fclose(input) ;
  
  nbi_boundary_condition_set(s, &(src[0]), 2, &(src[1]), 2, bc) ;
  
  if ( !greens_id && !layer_potentials ) {
    /*solver settings for GMRES*/
    i = nbi_surface_node_number(s) ;
        
    solver_work_size = 2*i + 4*gmres_restart +
      (gmres_restart+1)*(i+gmres_restart) + 2 ;
    solver_work_size *= 2 ;
  }
  
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
    fmm_work_size = MAX((guint)fmm_work_size,
			(order_max+1)*(order_max+1)*nqfmm*16) ;
    
    fmm_work_size += solver_work_size ;

    fprintf(stderr, "%s: allocating %d elements to FMM work space\n",
	    progname, fmm_work_size) ;
    
    work = (gdouble *)g_malloc0(fmm_work_size*sizeof(gdouble)) ;
    wbfmm_laplace_coaxial_translate_init(order_max+1) ;    

    fprintf(stderr, "%s: building tree; t=%lg\n",
	    progname, t = g_timer_elapsed(timer, NULL)) ;
    wbfmm_shift_angle_table_init() ;

    nbi_matrix_fmm_init(matrix, NBI_PROBLEM_LAPLACE,
			NULL, &(order[0]), 2, &(order[1]), 2,
			depth, dtree, shift_bw, precompute_local,
			work) ;    
    fprintf(stderr, "%s: FMM matrix initialized; t=%lg\n",
	    progname, t = g_timer_elapsed(timer, NULL)) ;
  } else {
    fmm_work_size = 16384 ;
    fmm_work_size += solver_work_size ;
    
    work = (gdouble *)g_malloc0(fmm_work_size*sizeof(gdouble)) ;    
  }
  
  if ( greens_id ) {
    fprintf(stderr, "%s: evaluating Green's identity; t=%lg\n",
	    progname, t = g_timer_elapsed(timer, NULL)) ;
    f = (gdouble *)g_malloc0(nbi_surface_node_number(s)*fstr*sizeof(gdouble)) ;
    nbi_surface_greens_identity_laplace(matrix,
					&(src[0]), 2, pwt,
					&(src[1]), 2, nwt,
					f, fstr, nthreads, work) ;
    fprintf(stderr, "%s: surface integration complete; t=%lg (%lg)\n",
	    progname,
	    g_timer_elapsed(timer, NULL), g_timer_elapsed(timer, NULL) - t) ;
    
    emax = fmax = e2 = f2 = 0.0 ;
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      xp = (NBI_REAL *)nbi_surface_node(s,i) ;
      /* G = greens_function_laplace(xp, xs) ; */
      G = src[2*i+0] ;
      fmax = MAX(fmax, G) ;
      emax = MAX(emax, fabs(G - f[i*fstr])) ;
      e2 += (G - f[i*fstr])*(G - f[i*fstr]) ;
      f2 += G*G ;
      fprintf(output, "%lg %lg %lg %lg %lg\n",
	      xp[0], xp[1], xp[2], f[i*fstr], fabs(G - f[i*fstr])) ;
    }
    
    fprintf(stderr, "L_inf norm: %lg; L_2 norm: %lg\n",
	    emax/fmax, sqrt(e2/f2)) ;

    return 0 ;
  }

  if ( layer_potentials ) {
    fprintf(stderr, "%s: evaluating double-layer potential; t=%lg\n",
    	    progname, t = g_timer_elapsed(timer, NULL)) ;
    matrix->potential = NBI_POTENTIAL_DOUBLE ;
    f = (gdouble *)g_malloc0(nbi_surface_node_number(s)*fstr*sizeof(gdouble)) ;
    nbi_matrix_multiply(matrix, &(src[0]), 2, 1.0, f, fstr, 0.0, nthreads,
			work) ;
    fprintf(stderr, "%s: evaluating single-layer potential; t=%lg\n",
    	    progname, t = g_timer_elapsed(timer, NULL)) ;
    matrix->potential = NBI_POTENTIAL_SINGLE ;
    nbi_matrix_multiply(matrix, &(src[1]), 2, -2.0, f, fstr, 2.0, nthreads,
			work) ;
    fprintf(stderr, "%s: surface integration complete; t=%lg (%lg)\n",
	    progname,
	    g_timer_elapsed(timer, NULL), g_timer_elapsed(timer, NULL) - t) ;
    
    emax = fmax = e2 = f2 = 0.0 ;
    for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
      xp = (NBI_REAL *)nbi_surface_node(s,i) ;
      /* G = greens_function_laplace(xp, xs) ; */
      G = src[2*i+0] ;
      fmax = MAX(fmax, G) ;
      emax = MAX(emax, fabs(G - f[i*fstr])) ;
      e2 += (G - f[i*fstr])*(G - f[i*fstr]) ;
      f2 += G*G ;
      fprintf(output, "%lg %lg %lg %lg %lg\n",
	      xp[0], xp[1], xp[2], f[i*fstr], fabs(G - f[i*fstr])) ;
    }
    
    fprintf(stderr, "L_inf norm: %lg; L_2 norm: %lg\n",
	    emax/fmax, sqrt(e2/f2)) ;

    return 0 ;
  }

#ifdef HAVE_PETSC
  /*if we get to here, we're doing a solve: set up PETSc*/
  Vec         x, b, u, rhs, sol ; /* approx solution, RHS, exact solution */
  Mat         A ;       /* linear system matrix */
  KSP         ksp ;     /* linear solver context */
  PC          pc ;      /* preconditioner context */
  PetscReal   norm ;    /* norm of solution error */
  PetscInt    n = 10 ;
  PetscMPIInt size ;
  PetscScalar value[3] ;
  gpointer petsc_ctx[NBI_SOLVER_DATA_SIZE] ;
  gdouble *buf ;
  gchar *help = "" ;
  FILE *fid ;
  
  PetscFunctionBeginUser;
  i = 0 ;
  PetscCall(PetscInitialize(&i, &argv, "options", help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCheck(size == 1, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE,
	     "This is a uniprocessor example only!");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create vectors.  Note that we form 1 vector from scratch and
     then duplicate as needed.
  */
  n = nbi_surface_node_number(s) ;
  PetscCall(VecCreate(PETSC_COMM_SELF, &rhs));
  PetscCall(PetscObjectSetName((PetscObject)rhs, "rhs"));
  PetscCall(VecSetSizes(rhs, PETSC_DECIDE, n));
  PetscCall(VecSetFromOptions(rhs));
  PetscCall(VecDuplicate(rhs, &sol));
  PetscCall(PetscObjectSetName((PetscObject)sol, "solution"));

  petsc_ctx[NBI_SOLVER_DATA_MATRIX]   =  matrix ;
  petsc_ctx[NBI_SOLVER_DATA_WORK]     =  work ;
  petsc_ctx[NBI_SOLVER_DATA_NTHREADS] = &nthreads ;
  
  PetscCall(MatCreateShell(PETSC_COMM_SELF, n, n, n, n, &petsc_ctx, &A)) ;
  PetscCall(MatShellSetOperation(A, MATOP_MULT, (void *)nbi_petsc_MatMult)) ;
  PetscCall(MatShellSetOperation(A, MATOP_SET_VALUES,
				 (void *)nbi_petsc_MatSetValues)) ;

  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));

  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

  /*set right hand side*/
  matrix->diag = 0.0 ;
  matrix->potential = NBI_POTENTIAL_SINGLE ;
  fprintf(stderr, "%s: forming right hand side, t=%lg\n", progname,
	  t=g_timer_elapsed(timer, NULL)) ;
  PetscCall(VecGetArrayRead(rhs, &buf)) ;

  nbi_matrix_multiply(matrix, &(src[1]), 2, 1.0, buf, 1, 0.0, nthreads, work) ;

  matrix->diag = -0.5 ;
  matrix->potential = NBI_POTENTIAL_DOUBLE ;
  fprintf(stderr, "%s: starting PETSc solver\n", progname) ;

  PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
  PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED)) ;
  PetscCall(KSPSetType(ksp, KSPGMRES));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPGetPC(ksp, &pc));
  PetscCall(PCSetType(pc, PCNONE)) ;
  PetscCall(KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT,
			     PETSC_DEFAULT));

  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSolve(ksp, rhs, sol));

  PetscCall(VecGetArrayRead(sol, &buf)) ;

  PetscCall(KSPGetIterationNumber(ksp, &i)) ;
  fprintf(stderr, "%s: %d iterations; t=%lg (%lg)\n",
	  progname, i,
	  g_timer_elapsed(timer, NULL),
	  g_timer_elapsed(timer, NULL) - t) ;
  PetscCall(VecDestroy(&rhs));
  PetscCall(VecDestroy(&sol));
  PetscCall(MatDestroy(&A));
  PetscCall(KSPDestroy(&ksp));

  PetscCall(PetscFinalize());
#else
  gdouble *rhs, *p, error = 0.0 ;

  fprintf(stderr,
	  "%s: compiled without PETSc, using internal GMRES solver\n\n"
	  "  If you want solver support, recompile NBI with PETSc\n"
	  "  available from https://petsc.org/\n\n",
	  progname) ;
  
  rhs = (gdouble *)g_malloc0(nbi_surface_node_number(s)*sizeof(gdouble)) ;
  p   = (gdouble *)g_malloc0(nbi_surface_node_number(s)*sizeof(gdouble)) ;
  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) p[i] = 1.0 ;
  
  /*form right hand side*/
  matrix->diag = 0.0 ;
  matrix->potential = NBI_POTENTIAL_SINGLE ;
  fprintf(stderr, "%s: forming right hand side, t=%lg\n", progname,
	  t=g_timer_elapsed(timer, NULL)) ;

  nbi_matrix_multiply(matrix, &(src[1]), 2, 1.0, rhs, 1, 0.0, nthreads, work) ;
  
  matrix->diag = -0.5 ;
  matrix->potential = NBI_POTENTIAL_DOUBLE ;
  fprintf(stderr, "%s: starting solver\n", progname) ;

  i = nbi_gmres_real(matrix, p, 1, rhs, 1, gmres_restart, gmres_max_iter, tol,
		     &error, nthreads, work) ;

  fprintf(stderr, "%s: %d iterations; error = %lg, t=%lg (%lg)\n",
	  progname, i, error,
	  g_timer_elapsed(timer, NULL),
	  g_timer_elapsed(timer, NULL) - t) ;
  buf = p ;
#endif
  
  emax = fmax = e2 = f2 = 0.0 ;
  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    emax = MAX(emax,fabs(buf[i]-src[2*i])) ;
    fmax = MAX(fmax, fabs(src[2*i])) ;
    e2 += (buf[i]-src[2*i])*(buf[i]-src[2*i]) ;
    f2 += src[2*i]*src[2*i] ;
  }

  fprintf(stderr, "%s: emax = %lg; fmax = %lg;\n", progname, emax, fmax) ;
  fprintf(stderr, "%s: L_inf norm = %lg; L_2 norm = %lg\n",
	  progname, emax/fmax, sqrt(e2/f2)) ;

  nbi_data_write(stdout, buf, 1, 1, nbi_surface_node_number(s)) ;

  return 0 ;
}
