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

#include "tinyexpr.h"

#include "nbi-private.h"

const gchar *variables[] = {"x", "y", "z", "nx", "ny", "nz", NULL} ;

nbi_expression_t *nbi_expression_new(gchar *expression)

{
  gint i, err ;  
  nbi_expression_t *e ;
  te_variable *vars ;

  e = (nbi_expression_t *)g_malloc0(sizeof(nbi_expression_t)) ;

  e->vars = g_malloc(NBI_EXPRESSION_VARIABLE_NUMBER*sizeof(te_variable)) ;
  vars = e->vars ;

  for ( i = 0 ; i < NBI_EXPRESSION_VARIABLE_NUMBER ; i ++ ) {
    vars[i].name = g_strdup(variables[i]) ;
    vars[i].address = &(e->x[i]) ;
    vars[i].type = TE_VARIABLE ;
    vars[i].context = NULL ;
  }

  e->compiled = te_compile(expression, vars, NBI_EXPRESSION_VARIABLE_NUMBER,
			    &err) ;

  if ( err != 0 ) {
    fprintf(stderr,
	    "%s: cannot parse expression\n"
	    "    %s\n"
	    "error at character %d\n", __FUNCTION__, expression, err) ;
    exit(-1) ;
  }

  e->expression = g_strdup(expression) ;
  
  return e ;
}

gdouble nbi_expression_eval(nbi_expression_t *e, gdouble *x, gdouble *n)

{
  gdouble f ;
  /* te_variable *vars = e->vars ; */

  e->x[0] = x[0] ; e->x[1] = x[1] ; e->x[2] = x[2] ; 
  e->x[3] = n[0] ; e->x[4] = n[1] ; e->x[5] = n[2] ; 
  
  f = te_eval((te_expr *)(e->compiled)) ;
  
  return f ;
}  

nbi_boundary_condition_t *nbi_boundary_condition_new(nbi_problem_t problem)

{
  nbi_boundary_condition_t *b ;

  b = (nbi_boundary_condition_t *)
    g_malloc0(sizeof(nbi_boundary_condition_t)) ;

  b->problem = problem ;
  b->e[0] = b->e[1] = b->e[2] = b->e[3] = NULL ;
  
  return b ;
}

gint nbi_boundary_condition_add(nbi_boundary_condition_t *b, gchar *e)

{
  switch ( b->problem ) {
  default: g_assert_not_reached() ; break ;
  case NBI_PROBLEM_LAPLACE:
    if ( b->e[0] == NULL ) {
      /*surface potential*/
      b->e[0] = nbi_expression_new(e) ;
      return 0 ;
    }
    if ( b->e[1] == NULL ) {
      /*surface potential derivative*/
      b->e[1] = nbi_expression_new(e) ;
      return 0 ;
    }
    fprintf(stderr,
	    "%s: only two boundary condition terms permitted for Laplace "
	    "problem\n", __FUNCTION__) ;
    exit(-1) ;
    break ;
  }
  
  return 0 ;
}
