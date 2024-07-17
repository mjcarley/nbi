/* This file is part of NBI, a library for Nystrom Boundary Integral solvers
 *
 * Copyright (C) 2021, 2023 Michael Carley
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <glib.h>

#include <nbi.h>

#include "tinyexpr.h"

#include "nbi-private.h"

const char *variables[] = {"x", "y", "z", "nx", "ny", "nz", "rp", "ip",
  "rdp", "idp", NULL} ;

/**
 * @ingroup boundary
 *
 * @{
 * 
 */

const te_variable functions[] = {
  {"laplace_G" , nbi_function_gfunc_laplace_G , TE_FUNCTION3},
  {"laplace_dG", nbi_function_gfunc_laplace_dG, TE_FUNCTION6},
  {"helmholtz_Gr" , nbi_function_gfunc_helmholtz_G_real , TE_FUNCTION4},
  {"helmholtz_Gi" , nbi_function_gfunc_helmholtz_G_imag , TE_FUNCTION4},
  {"helmholtz_dGr" , nbi_function_gfunc_helmholtz_dG_real, TE_FUNCTION7},
  {"helmholtz_dGi" , nbi_function_gfunc_helmholtz_dG_imag, TE_FUNCTION7},
  {"helmholtz_ring_r", nbi_function_gfunc_helmholtz_ring_real, TE_FUNCTION6},
  {"helmholtz_ring_i", nbi_function_gfunc_helmholtz_ring_imag, TE_FUNCTION6},
  {"helmholtz_ring_x_r", nbi_function_gfunc_helmholtz_ring_real_x,
   TE_FUNCTION6},
  {"helmholtz_ring_x_i", nbi_function_gfunc_helmholtz_ring_imag_x,
   TE_FUNCTION6},
  {"helmholtz_ring_y_r", nbi_function_gfunc_helmholtz_ring_real_y,
   TE_FUNCTION6},
  {"helmholtz_ring_y_i", nbi_function_gfunc_helmholtz_ring_imag_y,
   TE_FUNCTION6},
  {"helmholtz_ring_z_r", nbi_function_gfunc_helmholtz_ring_real_z,
   TE_FUNCTION6},
  {"helmholtz_ring_z_i", nbi_function_gfunc_helmholtz_ring_imag_z,
   TE_FUNCTION6},
  {"sphere_scattered_r", nbi_function_sphere_scattered_r,
   TE_FUNCTION5},
  {"sphere_scattered_i", nbi_function_sphere_scattered_i,
   TE_FUNCTION5},
  {NULL , NULL, TE_FUNCTION0}
} ;

const char *function_help[] = {
  "laplace_G(x, y, z)",
  " = 1/4 PI R; R^2 = x^2 + y^2 + z^2;\n"
  "  Green's function for Laplace equation",
  "laplace_dG(x, y, z, nx, ny, nz)",
  " = (nx*dG/dx + ny*dG/dy + nz*dG/dz)\n"
  "  normal derivative of Green's function for Laplace equation",
  "helmholtz_Gr(k, x, y, z)",
  " = cos(kR)/4 PI R;\n"
  "  real part of Green's function for Helmholtz equation",
  "helmholtz_Gi(k, x, y, z)",
  " = sin(kR)/4 PI R;\n"
  "  imaginary part of Green's function for Helmholtz equation",
  "helmholtz_dGr(k, x, y, z, nx, ny, nx)",
  " = (nx*dGr/dx + ny*dGr/dy + nz*dGr/dz);\n"
  "  real part of normal derivative of Green's function for Helmholtz equation",
  "helmholtz_dGi(k, x, y, z, nx, ny, nx)",
  " = (nx*dGi/dx + ny*dGi/dy + nz*dGi/dz);\n"
  "  imaginary part of normal derivative of Green's function for "
  "Helmholtz equation",
  "helmholtz_ring_r(a, n, k, x, y, z)",
  "\n  real part of field from ring source of radius a, azimuthal order n,"
  "\n  and wavenumber k",
  "helmholtz_ring_i(a, n, k, x, y, z)",
  "\n  imaginary part of field from ring source of radius a, azimuthal "
  "order n,\n  and wavenumber k",
  "helmholtz_ring_x_r(a, n, k, x, y, z)",
  "\n  real part of x derivative of field from ring source",
  "helmholtz_ring_x_i(a, n, k, x, y, z)",
  "\n  imaginary part of x derivative of field from ring source",
  "helmholtz_ring_y_r(a, n, k, x, y, z)",
  "\n  real part of y derivative of field from ring source",
  "helmholtz_ring_y_i(a, n, k, x, y, z)",
  "\n  imaginary part of y derivative of field from ring source",
  "helmholtz_ring_z_r(a, n, k, x, y, z)",
  "\n  real part of z derivative of field from ring source",
  "helmholtz_ring_z_i(a, n, k, x, y, z)",
  "\n  imaginary part of z derivative of field from ring source",
  "sphere_scattered_r(a, k, x, y, z)",
  "\n  real part of potential scattered from sphere of radius a by incident"
  "\n  plane wave exp[j k z]",
  "sphere_scattered_i(a, k, x, y, z)",
  "\n  imaginary part of potential scattered from sphere of radius a by"
  "\n  incident plane wave exp[j k z]",
  NULL, NULL
} ;

const char *nbi_function_help(char *f)

{
  gint i ;

  for ( i = 0 ; function_help[2*i+0] != NULL ; i ++ ) {
    if ( strncmp(function_help[2*i+0], f, strlen(f)) == 0 )
      return function_help[2*i+1] ;
  }
  
  return NULL ;
}

/** 
 * Write a list of built-in functions to file
 * 
 * @param f output file stream;
 * @param help if TRUE, also write help text for functions.
 * 
 * @return 0 on success.
 */

gint nbi_functions_list(FILE *f, gboolean help)

{

  gint i ;

  if ( help ) {
    for ( i = 0 ; function_help[2*i+0] != NULL ; i ++ ) {
      fprintf(f, "%s", function_help[2*i+0]) ;
      fprintf(f, "%s\n\n", function_help[2*i+1]) ;
    }

    return 0 ;
  }

  for ( i = 0 ; function_help[2*i+0] != NULL ; i ++ ) {
    fprintf(f, "%s\n\n", function_help[2*i+0]) ;
  }

  return 0 ;
}

/** 
 * Allocate a new boundary condition evaluator for a specified problem.
 * 
 * @param problem problem to which boundary condition applies.
 * 
 * @return pointer to newly allocated ::nbi_boundary_condition_t
 */

nbi_boundary_condition_t *nbi_boundary_condition_new(nbi_problem_t problem)

{
  nbi_boundary_condition_t *b ;
  te_variable *vars ;
  gint i ;
  
  b = (nbi_boundary_condition_t *)g_malloc0(sizeof(nbi_boundary_condition_t)) ;

  nbi_boundary_condition_problem(b) = problem ;
  
  b->vars = g_malloc0(NBI_BOUNDARY_CONDITION_VARIABLE_NUMBER*
		      sizeof(te_variable)) ;
  vars = (te_variable *)(b->vars) ;

  /*initial variables*/
  for ( b->nvars = 0 ; variables[b->nvars] != NULL ; b->nvars ++ ) {
    vars[b->nvars].name = g_strdup(variables[b->nvars]) ;
    vars[b->nvars].address = &(b->x[b->nvars]) ;
    vars[b->nvars].type = TE_VARIABLE ;
    vars[b->nvars].context = NULL ;
    b->expr[b->nvars] = NULL ;
  }

  for ( i = 0 ; functions[i].name != NULL ; i ++ ) {
    vars[b->nvars] = functions[i] ;
    vars[b->nvars].context = NULL ;
    b->nvars ++ ;
  }  

  /*check to make sure nothing has gone wrong in macros*/
  g_assert(strcmp(vars[NBI_BOUNDARY_CONDITION_POINT].name, "x") == 0) ;
  g_assert(strcmp(vars[NBI_BOUNDARY_CONDITION_NORMAL].name, "nx") == 0) ;
  g_assert(strcmp(vars[NBI_BOUNDARY_CONDITION_P_REAL].name, "rp") == 0) ;
  g_assert(strcmp(vars[NBI_BOUNDARY_CONDITION_P_IMAG].name, "ip") == 0) ;
  g_assert(strcmp(vars[NBI_BOUNDARY_CONDITION_DP_REAL].name, "rdp") == 0) ;
  g_assert(strcmp(vars[NBI_BOUNDARY_CONDITION_DP_IMAG].name, "idp") == 0) ;
  
  return b ;
}

static gint getlongline(FILE *f, char *buf, gint nbmax, char *cont)

{
  char *line ;
  gint i ;
  gint len ;
  gsize n ;
  ssize_t nc ;
  
  buf[0] = '\0' ;
  n = 0 ;
  nc = getline(&line, &n, f) ;
  if ( nc == -1 ) {
    return -1 ;    
  }

  /*check for comment or empty line*/
  line = g_strchug(line) ;
  if ( line[0] == '#' || strlen(line) == 0 ) {
    return getlongline(f, buf, nbmax, cont) ;
  }
  
  if ( n + strlen(buf) > nbmax )
    g_error("%s: string buffer overrun (%lu>%d)",
	    __FUNCTION__, n + strlen(buf), nbmax) ;

  strcpy(buf, line) ;
  len = strlen(buf) ;
  
  /*check for continuation character*/
  for ( i = 0 ; i < strlen(cont) ; i ++ ) {
    if ( buf[len-2] == cont[i] ) {
      return getlongline(f, &(buf[len-2]), nbmax-len, cont) ;
    }
  }

  return 0 ;
}

/** 
 * Read the definition of a boundary condition from file
 * 
 * @param f input file stream;
 * @param b boundary condition evaluator, allocated with 
 * ::nbi_boundary_condition_new.
 * 
 * @return 0 on success.
 */

gint nbi_boundary_condition_read(FILE *f, nbi_boundary_condition_t *b)

{
  char **tokens, buf[16384], *cont = "\\" ;
  te_variable *vars ;
  ssize_t nc ;
  gint nt, i, err, nbmax = 16384 ;
  
  for ( i = 0 ; i < NBI_BOUNDARY_CONDITION_VARIABLE_NUMBER ; i ++ ) {
    b->compiled[i] = NULL ;
  }
  nc = 0 ;
  vars = b->vars ;
  nc = getlongline(f, buf, nbmax, cont) ;
  while ( nc != -1 ) {
    tokens = g_strsplit(buf, "=", 0) ;

    for ( nt = 0 ; tokens[nt] != NULL ; nt ++ ) ;
    g_assert(nt == 2) ;
    for ( i = 0 ; i < nt ; i ++ ) {
      tokens[i] = g_strchug(g_strchomp(tokens[i])) ;
    }

    i = nbi_boundary_condition_has_variable(b, tokens[0]) ;

    if ( i != -1 ) {
      b->expr[i] = g_strdup(tokens[1]) ;
    } else {
      g_assert(b->nvars < NBI_BOUNDARY_CONDITION_VARIABLE_NUMBER) ;
      i = b->nvars ;
      vars[i].name = g_strdup(tokens[0]) ;
      vars[i].address = &(b->x[i]) ;
      vars[i].type = TE_VARIABLE ;

      b->expr[i] = g_strdup(tokens[1]) ;
      
      b->nvars ++ ;
    }
    nc = getlongline(f, buf, nbmax, cont) ;
  }

  /*compile expressions*/
  for ( i = 6 ; i < b->nvars ; i ++ ) {
    if ( vars[i].type == TE_VARIABLE && b->expr[i] != NULL ) 
      b->compiled[i] = te_compile(b->expr[i], vars, b->nvars, &err) ;

    if ( err != 0 ) {
      fprintf(stderr,
	      "%s: cannot parse expression\n"
	      "    %s\n"
	      "error at character %d\n", __FUNCTION__, b->expr[i], err) ;
      exit(-1) ;
    } 
  }
  
  return 0 ;
}

/** 
 * Check if a boundary condition evaluator includes a particular variable
 * 
 * @param b boundary condition evaluator to check;
 * @param v name of variable to find.
 * 
 * @return index of \a v in the variable list of \a b, or -1 if not found.
 */

gint nbi_boundary_condition_has_variable(nbi_boundary_condition_t *b,
					 char *v)

{
  gint i ;
  te_variable *vars = (te_variable *)(b->vars) ;
  
  for ( i = 0 ; i < b->nvars ; i ++ ) {
    if ( strcmp(v, vars[i].name) == 0 ) return i ;
  }
  
  return -1 ;
}

/** 
 * Write a boundary condition evaluator to file
 * 
 * @param f output file stream;
 * @param b boundary condition evaluator.
 * 
 * @return 0 on success.
 */

gint nbi_boundary_condition_write(FILE *f, nbi_boundary_condition_t *b)

{
  gint i ;
  te_variable *vars = b->vars ;
  
  for ( i = 6 ; i < b->nvars ; i ++ ) {
    if ( vars[i].type == TE_VARIABLE ) 
      fprintf(f, "%s = %s\n", vars[i].name, b->expr[i]) ;
  }
  
  return 0 ;
}

/** 
 * Evaluate a boundary condition at some point
 * 
 * @param b boundary condition evaluator;
 * @param x point for evaluation;
 * @param n surface normal at \a x.
 * 
 * @return 0 on success.
 */

gint nbi_boundary_condition_eval(nbi_boundary_condition_t *b, gdouble *x,
				 gdouble *n)

{
  te_variable *vars = b->vars ;
  gint i ;
  
  b->x[NBI_BOUNDARY_CONDITION_POINT+0] = x[0] ;
  b->x[NBI_BOUNDARY_CONDITION_POINT+1] = x[1] ;
  b->x[NBI_BOUNDARY_CONDITION_POINT+2] = x[2] ;
  
  b->x[NBI_BOUNDARY_CONDITION_NORMAL+0] = n[0] ;
  b->x[NBI_BOUNDARY_CONDITION_NORMAL+1] = n[1] ;
  b->x[NBI_BOUNDARY_CONDITION_NORMAL+2] = n[2] ;

  /*update the user-supplied variables*/
  for ( i = NBI_BOUNDARY_CONDITION_P_REAL + 4 ; i < b->nvars ; i ++ ) {
    if ( vars[i].type == TE_VARIABLE ) {
      b->x[i] = te_eval((te_expr *)(b->compiled[i])) ;
    }
  }

  /*now evaluate any expressions for the boundary conditions proper*/
  for ( i = NBI_BOUNDARY_CONDITION_P_REAL ;
	i < NBI_BOUNDARY_CONDITION_P_REAL + 4 ; i ++ ) {
    if ( b->expr[i] != NULL ) {
      b->x[i] = te_eval((te_expr *)(b->compiled[i])) ;
    }
  }

  return 0 ;
}

/** 
 * Check if variable is defined in a boundary condition (simple
 * version of ::nbi_boundary_condition_has_variable)
 * 
 * @param b boundary condition evaluator to check;
 * @param v name of variable to find.
 * 
 * @return TRUE if \a v is defined in \a b, FALSE otherwise.
 */

gboolean nbi_boundary_condition_defined(nbi_boundary_condition_t *b,
					char *v)

{
  te_variable *vars = b->vars ;
  gint i ;

  for ( i = 0 ; i < b->nvars ; i ++ ) {
    if ( (strcmp(vars[i].name, v) == 0) &&
	 b->expr[i] != NULL ) return TRUE ;
  }
  
  return FALSE ;
}

/**
 *
 * @}
 *
 */
