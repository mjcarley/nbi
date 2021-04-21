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

#ifndef NBI_H_INCLUDED
#define NBI_H_INCLUDED

#define NBI_SURFACE_NODE_LENGTH   7
#define NBI_UPSAMPLE_NODE_LENGTH  7
#define NBI_SURFACE_PATCH_LENGTH  2
#define NBI_SURFACE_PATCH_DATA_LENGTH     4

#include <stdio.h>

#include <wbfmm.h>

typedef enum
  {
   NBI_PROBLEM_LAPLACE = 1, /**< Laplace equation */
   NBI_PROBLEM_HELMHOLTZ = 2 /**< Helmholtz equation */	      
  } nbi_problem_t ;	      

typedef enum
  {
   NBI_POTENTIAL_SINGLE = 1, /*< single-layer potential 
			       \f$\int G\sigma dS\f$*/
   NBI_POTENTIAL_DOUBLE = 2  /*< double-layer potential 
			       \f$\int \partial G/\partial n\sigma dS\f$*/
  } nbi_potential_t ; 

typedef struct _nbi_surface_t nbi_surface_t ;

struct _nbi_surface_t {
  gint
  nn,         /*< number of nodes */
    nnmax,    /*< maximum number of nodes (allocated) */
    np,       /*< number of patches */
    npmax,    /*< maximum number of patches */
    *ip ;     /*< index of first node and number of nodes on each patch*/
  gchar
  *pcr,       /*< patch data (centroids and radii)*/
    *xc ;       /*< collocation nodes, normals, and quadrature weights*/
  gsize
  fpsize ; /*< size of floating point data (single or double precision)*/
} ;

typedef struct _nbi_matrix_t nbi_matrix_t ;


#define nbi_surface_node_number(_s)  ((_s)->nn)
#define nbi_surface_node_number_max(_s)  ((_s)->nnmax)

#define nbi_surface_patch_number(_s)  ((_s)->np)
#define nbi_surface_patch_number_max(_s)  ((_s)->npmax)

#define nbi_surface_patch_node(_s,_i)		\
  ((_s)->ip[NBI_SURFACE_PATCH_LENGTH*(_i)+0])
#define nbi_surface_patch_node_number(_s,_i)	\
  ((_s)->ip[NBI_SURFACE_PATCH_LENGTH*(_i)+1])

#define nbi_surface_patch_centre(_s,_i)			\
  (&((_s)->pcr[(NBI_SURFACE_PATCH_DATA_LENGTH*(_i)+0)*((_s)->fpsize)]))
#define nbi_surface_patch_sphere_radius(_s,_i)		\
  (&((_s)->pcr[(NBI_SURFACE_PATCH_DATA_LENGTH*(_i)+3)*((_s)->fpsize)]))

#define nbi_surface_node(_s,_i)			\
  &((_s)->xc[(_i)*NBI_SURFACE_NODE_LENGTH*((_s)->fpsize)])
#define nbi_surface_normal(_s,_i)		\
  &((_s)->xc[((_i)*NBI_SURFACE_NODE_LENGTH+3)*((_s)->fpsize)])

#define nbi_surface_node_weight(_s,_i)		\
  &((_s)->xc[((_i)*NBI_SURFACE_NODE_LENGTH+6)*((_s)->fpsize)])

struct _nbi_matrix_t {
  nbi_surface_t *s ;
  gsize fpsize ; /*< floating point data size */
  nbi_problem_t problem ;
  gint
  ustr,     /*< upsampled patch node stride */
    pstr,   /*< upsampled node potential stride */
    nstr,   /*< upsampled node normal derivative stride */
    *idxu, /*< upsampled node indices */
    *idxp,
    *idx ;
  gchar
  *Ast,    /*< correction matrices */
    *xu,   /*< upsampled nodes */
    *bc,    /*< upsampled boundary condition buffer */
    *p, *pn ;
  wbfmm_tree_t 
  *tree ;  /*< FMM tree if applicable */
  wbfmm_target_list_t
  *targets ;  /*< FMM target data if applicable */
  wbfmm_shift_operators_t
  *shifts ;   /*< FMM shift operators if applicable */
} ;

nbi_surface_t *nbi_surface_alloc(gint nnmax, gint npmax) ;
gint nbi_surface_write(nbi_surface_t *s, FILE *f) ;
nbi_surface_t *nbi_surface_read(FILE *f) ;
gint nbi_surface_patch_centroid(gdouble *x, gint xstr,
				gdouble *w, gint wstr,
				gint nx,
				gdouble *c) ;
gdouble nbi_surface_patch_radius(gdouble *x, gint xstr, gint nx, gdouble *c) ;
gint nbi_surface_set_patch_data(nbi_surface_t *s) ;

gint nbi_geometry_sphere(nbi_surface_t *s, gdouble r, gint nth, gint nph,
			 gint nq) ;
gint nbi_geometry_ellipsoid(nbi_surface_t *s,
			    gdouble a, gdouble b, gdouble c,
			    gint nth, gint nph,
			    gint nq) ;
gint nbi_geometry_sphere_ico(nbi_surface_t *s, gdouble r, gint nr, gint nq) ;
gint nbi_geometry_ellipsoid_ico(nbi_surface_t *s,
				gdouble a, gdouble b, gdouble c,
				gint nr, gint nq) ;

gint nbi_patch_neighbours(gdouble *c, gdouble r,
			  gdouble *x, gint xstr, gint nx,
			  gint n0, gint n1,
			  gint *nbrs, gint *nnbrs, gint nnmax) ;
gint nbi_surface_patch_neighbours(nbi_surface_t *s,
				  gint p, gdouble r,
				  gint *nbrs, gint *nnbrs, gint nnmax) ;

gdouble *nbi_patch_upsample_matrix(gint ns, gint nu) ;
gint nbi_element_interp_matrix(gint ns, gdouble **K, gint *Nk) ;


nbi_matrix_t *nbi_matrix_new(nbi_surface_t *s) ;
gint nbi_matrix_read(FILE *input, nbi_matrix_t *m) ;
gint nbi_matrix_upsample_laplace(nbi_matrix_t *m,
				 gdouble *p, gint pstr,
				 gdouble *pn, gint nstr,
				 gdouble pwt, gdouble nwt) ;
  
#endif /*NBI_H_INCLUDED*/
