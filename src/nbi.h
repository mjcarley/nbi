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

#ifndef NBI_H_INCLUDED
#define NBI_H_INCLUDED

#define NBI_SURFACE_NODE_LENGTH   7
#define NBI_UPSAMPLE_NODE_LENGTH  7
#define NBI_SURFACE_PATCH_LENGTH  2
#define NBI_SURFACE_PATCH_DATA_LENGTH     4

#define NBI_HEADER_LENGTH  80
#define NBI_HEADER_ID       0
#define NBI_HEADER_VERSION  4
#define NBI_HEADER_TYPE     8
#define NBI_HEADER_FORMAT  12

#define NBI_HEADER_DATA    40

#include <stdio.h>

#include <wbfmm.h>

#define NBI_SOLVER_DATA_SIZE     8
#define NBI_SOLVER_DATA_MATRIX   0
#define NBI_SOLVER_DATA_WORK     1
#define NBI_SOLVER_DATA_NTHREADS 2

#ifdef DOXYGEN
/**
 * @file   nbi.h
 * @author  <michael@michael.paraffinalia.co.uk>
 * @date   Thu Aug 24 12:48:32 2023
 * 
 * @brief  
 * 
 * 
 */
#endif /*DOXYGEN*/


#ifdef DOXYGEN
/**
 * @typedef nbi_problem_t
 * 
 * Enumerated data type to specify problem to be solved
 */
#endif /*DOXYGEN*/
typedef enum
  {
   NBI_PROBLEM_UNDEFINED = 0, /**< problem not defined */
   NBI_PROBLEM_LAPLACE = 1, /**< Laplace equation */
   NBI_PROBLEM_HELMHOLTZ = 2 /**< Helmholtz equation */	      
  } nbi_problem_t ;	      

#ifdef DOXYGEN
/**
 * @typedef nbi_potential_t
 * 
 * Enumerated data type to select single or double layer potential
 */
#endif /*DOXYGEN*/
typedef enum
  {
   NBI_POTENTIAL_UNDEFINED = 0,
   NBI_POTENTIAL_SINGLE = 1, /**< single-layer potential 
				\f$\int G\sigma dS\f$*/
   NBI_POTENTIAL_DOUBLE = 2  /**< double-layer potential 
				\f$\int \partial G/\partial n\sigma dS\f$*/
  } nbi_potential_t ; 

#ifdef DOXYGEN
/**
 * @{
 * 
 * @ingroup surfaces
 */

/**
 * @typedef nbi_surface_t
 * 
 * Basic data type for NBI surfaces, including element and quadrature
 * information
 * 
 */
typedef nbi_surface_t ;

/**
 * @brief number of nodes on an ::nbi_surface_t
 */
#define nbi_surface_node_number(s)  
/**
 * @brief maximum number of nodes on an ::nbi_surface_t \a s
 */
#define nbi_surface_node_number_max(s)
/**
 * @brief number of surface patches on an ::nbi_surface_t \a s
 */
#define nbi_surface_patch_number(s)  
/**
 * @brief maximum number of surface patches on an ::nbi_surface_t \a s
 */
#define nbi_surface_patch_number_max(s)
/**
 * @brief index of first node on patch \a i of surface \a s
 */
#define nbi_surface_patch_node(s,i)
/**
 * @brief number of nodes on patch \a i of surface \a s
 */
#define nbi_surface_patch_node_number(s,i)
/**
 * @brief coordinates of estimated centroid of patch \a i of surface \a s
 */
#define nbi_surface_patch_centre(s,i)			
/**
 * @brief estimated radius of circumsphere of patch \a i of surface \a s
 */
#define nbi_surface_patch_sphere_radius(s,i)		
/**
 * @brief pointer to node \a i of surface \a s
 */
#define nbi_surface_node(s,i)			
/**
 * @brief pointer to normal \a i of surface \a s
 */
#define nbi_surface_normal(s,i)		
/**
 * @brief pointer to quadrature weight of node \a i of surface \a s
 */
#define nbi_surface_node_weight(s,i)		
/**
 * @}
 */

#else  /*DOXYGEN*/
typedef struct _nbi_surface_t nbi_surface_t ;

struct _nbi_surface_t {
  gint
  nn,         /*< number of nodes */
    nnmax,    /*< maximum number of nodes (allocated) */
    np,       /*< number of patches */
    npmax,    /*< maximum number of patches */
    *ip ;     /*< index of first node and number of nodes on each patch*/
  gchar  *pcr,       /*< patch data (centroids and radii)*/
    *xc ;       /*< collocation nodes, normals, and quadrature weights*/
  gsize
  fpsize ; /*< size of floating point data (single or double precision)*/
} ;

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
#endif /*DOXYGEN*/

#ifdef DOXYGEN
/**
 * @{
 * 
 * @ingroup matrix
 */
/**
 * @typedef nbi_matrix_t
 * 
 * Basic data type for NBI matrices, containing surface information,
 * local correction matrices and Fast Multipole Method data
 * 
 */
typedef nbi_matrix_t ;

#else  /*DOXYGEN*/
typedef struct _nbi_matrix_t nbi_matrix_t ;

#define NBI_PROBLEM_DATA_SIZE       8
#define NBI_PROBLEM_DATA_WAVENUMBER 0

struct _nbi_matrix_t {
  nbi_surface_t *s ;
  gsize fpsize ; /*< floating point data size */
  gdouble diag, /*< in multiplications add diag*I */
    pdata[NBI_PROBLEM_DATA_SIZE] ; /*< problem specific data */
  nbi_problem_t problem ;
  nbi_potential_t potential ;
  gint
  ustr,     /*< upsampled patch node stride */
    pstr,   /*< upsampled node potential stride */
    nstr,   /*< upsampled node normal derivative stride */
    *idxu,  /*< upsampled node indices */
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
  gdouble *test ; /*used for testing solvers (and only for that)*/
} ;

#define nbi_matrix_wavenumber(_m) ((_m)->pdata[NBI_PROBLEM_DATA_WAVENUMBER])
#endif /*DOXYGEN*/

#ifdef DOXYGEN
/**
 * @typedef nbi_boundary_condition_t
 *
 * Specification of boundary conditions using analytical expressions
 * 
 */

typedef nbi_boundary_condition_t ;
/**
 * @brief problem type (::nbi_problem_t) for which boundary condition
 * is evaluated
 */
#define nbi_boundary_condition_problem(b)
/**
 * @brief boundary potential in real problems (synonym for
 * ::nbi_boundary_condition_p_real)
 */
#define nbi_boundary_condition_p(b)		
/**
 * @brief real part of boundary potential in complex problems
 */
#define nbi_boundary_condition_p_real(b)	
/**
 * @brief imaginary part of boundary potential in complex problems
 */
#define nbi_boundary_condition_p_imag(b)	
/**
 * @brief 
 */
/**
 * @brief boundary potential normal derivative in real problems (synonym for
 * ::nbi_boundary_condition_dp_real)
 */
#define nbi_boundary_condition_dp(b)	
/**
 * @brief real part of boundary potential gradient in complex problems
 */
#define nbi_boundary_condition_dp_real(b)	
/**
 * @brief imaginary part of boundary potential gradient in complex problems
 */
#define nbi_boundary_condition_dp_imag(b)	

#else /*DOXYGEN*/

typedef struct _nbi_boundary_condition_t nbi_boundary_condition_t ;

#define NBI_BOUNDARY_CONDITION_VARIABLE_NUMBER 64
struct _nbi_boundary_condition_t {
  nbi_problem_t problem ;
  gint nvars ;
  gdouble x[NBI_BOUNDARY_CONDITION_VARIABLE_NUMBER] ;
  gpointer compiled[NBI_BOUNDARY_CONDITION_VARIABLE_NUMBER] ;
  /*< compiled version for evaluation by tinyexpr */
  gpointer vars ;
  gchar *expr[NBI_BOUNDARY_CONDITION_VARIABLE_NUMBER] ;
} ;

#define NBI_BOUNDARY_CONDITION_POINT   0
#define NBI_BOUNDARY_CONDITION_NORMAL  3
#define NBI_BOUNDARY_CONDITION_P_REAL  6
#define NBI_BOUNDARY_CONDITION_P_IMAG  7
#define NBI_BOUNDARY_CONDITION_DP_REAL 8
#define NBI_BOUNDARY_CONDITION_DP_IMAG 9

#define nbi_boundary_condition_problem(_b) ((_b)->problem) 

#define nbi_boundary_condition_p(_b)		\
  ((_b)->x[NBI_BOUNDARY_CONDITION_P_REAL])
#define nbi_boundary_condition_p_real(_b)	\
  ((_b)->x[NBI_BOUNDARY_CONDITION_P_REAL])
#define nbi_boundary_condition_p_imag(_b)	\
  ((_b)->x[NBI_BOUNDARY_CONDITION_P_IMAG])
#define nbi_boundary_condition_dp(_b)	\
  ((_b)->x[NBI_BOUNDARY_CONDITION_DP_REAL])
#define nbi_boundary_condition_dp_real(_b)	\
  ((_b)->x[NBI_BOUNDARY_CONDITION_DP_REAL])
#define nbi_boundary_condition_dp_imag(_b)	\
  ((_b)->x[NBI_BOUNDARY_CONDITION_DP_IMAG])
#endif /*DOXYGEN*/

#ifdef DOXYGEN
/**
 * @}
 */
#endif /*DOXYGEN*/

typedef gint (*nbi_mesh_export_func_t)(FILE *,
				       gdouble *, gint, gint,
				       gint *, gint, gint,
				       gdouble *) ;

nbi_surface_t *nbi_surface_alloc(gint nnmax, gint npmax) ;
gint nbi_surface_write(nbi_surface_t *s, FILE *f) ;
nbi_surface_t *nbi_surface_read(FILE *f) ;
gint nbi_surface_set_weights(nbi_surface_t *s) ;

gint nbi_header_read(FILE *f, gchar header[]) ;
gint nbi_header_write(FILE *f, gchar header[]) ;
gint nbi_header_init(gchar *header,
		     gchar *id,
		     gchar *version,
		     gchar *type,
		     gchar *format) ;
gint nbi_header_insert_string(gchar *header, gint i, gint len, gchar *str) ;
gchar *nbi_problem_type_string(nbi_problem_t p) ;
nbi_problem_t nbi_problem_from_string(gchar *p) ;

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
gint nbi_geometry_grid(nbi_surface_t *s,
		       gdouble smin, gdouble smax, gint ns,
		       gdouble tmin, gdouble tmax, gint nt,
		       gint nq) ;

gint nbi_patch_neighbours(gdouble *c, gdouble r,
			  gdouble *x, gint xstr, gint nx,
			  gint n0, gint n1,
			  gint *nbrs, gint *nnbrs, gint nnmax) ;
gint nbi_surface_patch_neighbours(nbi_surface_t *s,
				  gint p, gdouble r,
				  gint *nbrs, gint *nnbrs, gint nnmax) ;

gint nbi_matrix_neighbour_number_max(nbi_matrix_t *m) ;

gdouble *nbi_patch_upsample_matrix(gint ns, gint nu) ;
gint nbi_element_interp_matrix(gint ns, gdouble **K, gint *Nk) ;


nbi_matrix_t *nbi_matrix_new(nbi_surface_t *s) ;
gint nbi_matrix_read(FILE *input, nbi_matrix_t *m) ;
gint nbi_matrix_write(FILE *f, nbi_matrix_t *m) ;

gint nbi_matrix_upsample_laplace(nbi_matrix_t *m,
				 gdouble *p, gint pstr, gdouble pwt, 
				 gdouble *pn, gint nstr, gdouble nwt) ;
gint nbi_matrix_upsample_helmholtz(nbi_matrix_t *m,
				   gdouble *p, gint pstr,
				   gdouble pwt,
				   gdouble *pn, gint nstr,
				   gdouble nwt) ;

gint nbi_matrix_fmm_init(nbi_matrix_t *m,
			 nbi_problem_t problem,
			 wbfmm_shift_operators_t *shifts,
			 guint *order_s, gint sstr,
			 guint *order_r, gint rstr,
			 gint depth,
			 gdouble dtree,
			 gboolean shift_bw,
			 gboolean precompute_local,
			 gdouble *work) ;
gint nbi_matrix_multiply(nbi_matrix_t *A,
			 gdouble *x, gint xstr, gdouble al,
			 gdouble *y, gint ystr, gdouble bt,
			 gint nthreads,
			 gdouble *work) ;
gint nbi_matrix_multiply_helmholtz(nbi_matrix_t *A,
				   gdouble *x, gint xstr,
				   gdouble al,
				   gdouble *y, gint ystr,
				   gdouble bt,
				   gint nthreads,
				   gdouble *work) ;

gint nbi_calc_field_helmholtz(nbi_matrix_t *m,
			      gdouble *p ,
			      gint pstr,
			      gdouble pwt,
			      gdouble *pn,
			      gint nstr,
			      gdouble nwt,
			      gdouble *x,
			      gdouble *f,
			      gint nthreads,
			      gdouble *work) ;
gint nbi_surface_field_helmholtz(nbi_surface_t *s, gdouble k,
				 gdouble *ps, gint pstr,
				 gdouble *al, gdouble *bt,
				 gdouble *x, gdouble *p) ;
gint nbi_surface_field_laplace(nbi_surface_t *s,
			       gdouble *ps,
			       gint pstr,
			       gdouble al, gdouble bt,
			       gdouble *x, gdouble *p) ;

gint nbi_surface_greens_identity_helmholtz(nbi_matrix_t *m,
					   gdouble *p ,
					   gint pstr,
					   gdouble pwt,
					   gdouble *pn,
					   gint nstr,
					   gdouble nwt,
					   gdouble *f,
					   gint fstr,
					   gint nthreads,
					   gdouble *work) ;
gint nbi_surface_greens_identity_laplace(nbi_matrix_t *m,
					 gdouble *p , gint pstr, gdouble pwt,
					 gdouble *pn, gint nstr, gdouble nwt,
					 gdouble *f , gint fstr,
					 gint nthreads,
					 gdouble *work) ;
gint nbi_matrix_multiply_laplace(nbi_matrix_t *A,
				 gdouble *x, gint xstr, gdouble al,
				 gdouble *y, gint ystr, gdouble bt,
				 gint nthreads,
				 gdouble *work) ;

nbi_boundary_condition_t *nbi_boundary_condition_new(nbi_problem_t problem) ;
gint nbi_boundary_condition_add(nbi_boundary_condition_t *b, gchar *e) ;
gint nbi_boundary_condition_set(nbi_surface_t *s,
				gdouble *p , gint pstr,
				gdouble *pn, gint nstr,
				nbi_boundary_condition_t *b) ;
gint nbi_boundary_condition_has_variable(nbi_boundary_condition_t *b,
					 gchar *v) ;
gint nbi_boundary_condition_write(FILE *f, nbi_boundary_condition_t *b) ;
gint nbi_boundary_condition_eval(nbi_boundary_condition_t *b, gdouble *x,
				 gdouble *n) ;
gboolean nbi_boundary_condition_defined(nbi_boundary_condition_t *b,
					gchar *v) ;

gint nbi_gmres_real(nbi_matrix_t *A, gint n,
		    gdouble *x, gint xstr,
		    gdouble *b, gint bstr,
		    gint m, gint max_it,
		    gdouble tol, gdouble *error,
		    gint nthreads,
		    gdouble *work) ;
gint nbi_gmres_complex(nbi_matrix_t *A, 
		       gdouble *x, gint xstr,
		       gdouble *b, gint bstr,
		       gint m, gint max_it,
		       gdouble tol, gdouble *error,
		       gint nthreads,
		       gdouble *work) ;
gint nbi_gmres_workspace_size_complex(gint n, gint m) ;
gint nbi_gmres_workspace_size_real(gint n, gint m) ;

gdouble *nbi_data_read(FILE *f, gint *nd, gint *ne) ;
gint nbi_data_write(FILE *f, gdouble *dat, gint dstr, gint ne, gint nd) ;
gint nbi_boundary_condition_read(FILE *f, nbi_boundary_condition_t *bc) ;

nbi_surface_t *nbi_agg_mesh(gchar *file, gint nq) ;

nbi_surface_t *nbi_gmsh_mesh(gchar *file, gint nq) ;

gint nbi_mesh_triangulate(nbi_surface_t *s,
			  gint dmax,
			  gdouble *K, gint Nk,
			  gdouble *x, gint xstr, gint npmax,
			  gint *tri, gint tstr, gint ntmax,
			  gdouble *f, gint fstr,
			  gdouble *fi, gint fistr,
			  gint *np, gint *nt) ;
gint nbi_mesh_export_gmsh(FILE *f, gchar *view,
			  gdouble *x, gint xstr, gint np, gint offp,
			  gint *tri, gint tstr, gint nt, gint offt,
			  gdouble *data, gint dstr) ;

const gchar *nbi_function_help(gchar *f) ;
gint nbi_functions_list(FILE *f, gboolean help) ;

nbi_matrix_t *nbi_matrix_assemble_helmholtz(nbi_surface_t *s,
					    gdouble k,
					    gdouble eta,
					    gint nqa, gint dmax,
					    gdouble tol,
					    gint N, gint nu,
					    gint nnmax, gint nthreads) ;
nbi_matrix_t *nbi_matrix_assemble_laplace(nbi_surface_t *s,
					  gdouble eta,
					  gint nqa, gint dmax,
					  gdouble tol,
					  gint N, gint nu,
					  gint nnmax, gint nthreads) ;

#endif /*NBI_H_INCLUDED*/
