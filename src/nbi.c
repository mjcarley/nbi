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

#include "nbi-private.h"

nbi_surface_t *nbi_surface_alloc(gint nnmax, gint npmax)

{
  nbi_surface_t *s ;

  s = (nbi_surface_t *)g_malloc(sizeof(nbi_surface_t)) ;

  nbi_surface_node_number(s) = 0 ;
  nbi_surface_node_number_max(s) = nnmax ;
  nbi_surface_patch_number(s) = 0 ;
  nbi_surface_patch_number_max(s) = npmax ;

  s->x = (gdouble *)g_malloc0(nnmax*NBI_SURFACE_NODE_LENGTH*sizeof(gdouble)) ;
  s->ip = (gint *)g_malloc0(2*npmax*sizeof(gint)) ;
  
  return s ;
}

gint nbi_surface_write(nbi_surface_t *s, FILE *f)

{
  gint i, j ;
  
  fprintf(f, "%d %d\n",
	  nbi_surface_node_number(s), nbi_surface_patch_number(s)) ;

  for ( i = 0 ; i < nbi_surface_patch_number(s) ; i ++ ) {
    fprintf(f, "%d %d\n",
	    nbi_surface_patch_node(s, i),
	    nbi_surface_patch_node_number(s, i)) ;
  }

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    for ( j = 0 ; j < NBI_SURFACE_NODE_LENGTH ; j ++ )
      fprintf(f, " %lg", nbi_surface_node_element(s,i,j)) ;
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

gint nbi_surface_patch_centroid(gdouble *x, gint xstr,
				gdouble *w, gint wstr,
				gint nx,
				gdouble *c)

{
  gint i ;
  gdouble A ;
  
  c[0] = c[1] = c[2] = A = 0.0 ;

  for ( i = 0 ; i < nx ; i ++ ) {
    c[0] += x[i*xstr+0]*w[i*wstr] ;
    c[1] += x[i*xstr+1]*w[i*wstr] ;
    c[2] += x[i*xstr+2]*w[i*wstr] ;
    A    +=             w[i*wstr] ;
  }

  c[0] /= A ; c[1] /= A ; c[2] /= A ; 
  
  return 0 ;
}

gdouble nbi_surface_patch_radius(gdouble *x, gint xstr, gint nx, gdouble *c)

{
  gdouble r, ri ;
  gint i ;

  r = 0.0 ;
  for ( i = 0 ; i < nx ; i ++ ) {
    ri = nbi_vector_distance(&(x[i*xstr]), c) ;    
    r = MAX(r, ri) ;
  }
  
  return r ;
}
