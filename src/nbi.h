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

#define NBI_SURFACE_NODE_LENGTH 7

typedef struct _nbi_surface_t nbi_surface_t ;

struct _nbi_surface_t {
  gint
  nn,         /*< number of nodes */
    nnmax,    /*< maximum number of nodes (allocated) */
    np,       /*< number of patches */
    npmax,    /*< maximum number of patches */
    *ip ;     /*< index of first node and number of nodes on each patch*/
  gdouble
  *x ;        /*< nodes, normals, and quadrature weights*/
} ;

#define nbi_surface_node_number(_s)  ((_s)->nn)
#define nbi_surface_node_number_max(_s)  ((_s)->nnmax)
#define nbi_surface_patch_number(_s)  ((_s)->np)
#define nbi_surface_patch_number_max(_s)  ((_s)->npmax)
#define nbi_surface_patch_node(_s,_i)  ((_s)->ip[2*(_i)+0])
#define nbi_surface_patch_node_number(_s,_i)  ((_s)->ip[2*(_i)+1])
#define nbi_surface_node(_s,_i) &((_s)->x[(_i)*NBI_SURFACE_NODE_LENGTH])
#define nbi_surface_normal(_s,_i) &((_s)->x[(_i)*NBI_SURFACE_NODE_LENGTH+3])
#define nbi_surface_patch_local_node(_s,_i,_j)				\
  (nbi_surface_node((_s),(nbi_surface_patch_node((_s),(_i))+(_j))))
#define nbi_surface_node_element(_s,_i,_j)	\
  ((_s)->x[(_i)*NBI_SURFACE_NODE_LENGTH+(_j)])
#define nbi_surface_node_weight(_s,_i)		\
  ((_s)->x[(_i)*NBI_SURFACE_NODE_LENGTH+6])
  
nbi_surface_t *nbi_surface_alloc(gint nnmax, gint npmax) ;
gint nbi_surface_write(nbi_surface_t *s, FILE *f) ;
gint nbi_surface_patch_centroid(gdouble *x, gint xstr,
				gdouble *w, gint wstr,
				gint nx,
				gdouble *c) ;
gdouble nbi_surface_patch_radius(gdouble *x, gint xstr, gint nx, gdouble *c) ;

gint nbi_geometry_sphere(nbi_surface_t *s, gdouble r, gint nth, gint nph,
			 gint nq) ;
gint nbi_geometry_ellipsoid(nbi_surface_t *s,
			    gdouble a, gdouble b, gdouble c,
			    gint nth, gint nph,
			    gint nq) ;

#endif /*NBI_H_INCLUDED*/
