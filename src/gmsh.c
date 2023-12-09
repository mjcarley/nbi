/* This file is part of NBI, a library for Nystrom Boundary Integral solvers
 *
 * Copyright (C) 2023 Michael Carley
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


#ifdef HAVE_SQT
#include <sqt.h>
#endif /*HAVE_SQT*/

#ifdef HAVE_BLASWRAP
#include <blaswrap.h>
#endif /*HAVE_BLASWRAP*/

#ifdef HAVE_GMSHC_H
#include <gmshc.h>
#endif /*HAVE_GMSHC_H*/

#include "nbi-private.h"

#define GMSH_NODE_SIZE      5
#define GMSH_NODE_DATA_SIZE 3
#define GMSH_ELEMENT_SIZE   4

static gint compare_tags(gconstpointer e1, gconstpointer e2)

{
  const gint *t1 = e1 ;
  const gint *t2 = e2 ;

  if ( t1[0] < t2[0] ) return -1 ;
  if ( t1[0] > t2[0] ) return  1 ;
  
  return 0 ;
}

static gint find_tag(gsize *tags, gint ntags, gsize tag)

{
  gint i ;

  for ( i = 0 ; i < ntags ; i ++ ) if ( tags[i] == tag ) return i ;
  
  return -1 ;
}

static void gmesh_get(gdouble **p, gint **pt, gint *np,
		      gint **t, gdouble **uvt, gint *nt)

/*
 * on exit:
 * x:  mesh points [x, y, z, u, v]
 * pt: point tags [gmsh index,  entity tag]
 * np: number of points
 * t:  triangles [p1 p2 p3]
 * uvt: parametric coordinates of triangle nodes on each entity
 * nt: number of triangles
 */
  
{
  gint ierr, i, j, k, s, ip, it, **etypes, *node, *tags, offp ;
  gsize *ntags, nsurf, *etypes_n, ***etags, **etags_n, etags_nn,
    ***entags, **entags_n, entags_nn, *nx, nuv, **ptags ;
  gdouble **uv, **x ;

  /*number of points and surfaces in mesh*/
  gmshModelGetEntities(&tags, &nsurf, 2, &ierr) ;
  nsurf /= 2 ;
  x  = (gdouble **)g_malloc0(nsurf*sizeof(gdouble *)) ;
  uv = (gdouble **)g_malloc0(nsurf*sizeof(gdouble *)) ;
  nx = (gsize *)g_malloc0(nsurf*sizeof(gsize)) ;
  ptags = (gsize **)g_malloc0(nsurf*sizeof(gsize *)) ;
  ntags = (gsize *)g_malloc0(nsurf*sizeof(gsize)) ;

  etypes = (gint **)g_malloc0(nsurf*sizeof(gint **)) ;
  etypes_n = (gsize *)g_malloc0(nsurf*sizeof(gsize)) ;
  etags = (gsize ***)g_malloc0(nsurf*sizeof(gsize **)) ;
  entags = (gsize ***)g_malloc0(nsurf*sizeof(gsize **)) ;
  entags_n = (gsize **)g_malloc0(nsurf*sizeof(gsize *)) ;
  etags_n = (gsize **)g_malloc0(nsurf*sizeof(gsize *)) ;
  
  *np = *nt = 0 ;
  for ( s = 0 ; s < nsurf ; s ++ ) {
    gmshModelMeshGetNodes(&(ptags[s]), &(ntags[s]), &(x[s]), &(nx[s]),
			  &(uv[s]), &nuv, 2, tags[2*s+1], TRUE, TRUE, &ierr) ;
    if ( nx[s] != 0 ) {
      nx[s] /= 3 ;
      g_assert(nuv/2 == nx[s]) ;
      *np += nx[s] ;
    }
    gmshModelMeshGetElements(&(etypes[s]), &(etypes_n[s]),
			     &(etags[s]), &(etags_n[s]), &etags_nn,
			     &(entags[s]), &(entags_n[s]), &entags_nn,
			     2, tags[2*s+1], &ierr) ;
    /*make sure we only have triangles*/
    g_assert(etypes_n[s] == 1) ;
    g_assert(etypes[s][0] == 2) ;
    
    *nt += etags_n[s][0] ;
  }

  /*read nodes*/
  *p  = (gdouble *)g_malloc0(GMSH_NODE_SIZE*(*np)*sizeof(gdouble)) ;
  *pt = (gint    *)g_malloc0(GMSH_NODE_DATA_SIZE*(*np)*sizeof(gint)) ;
  *t   = (gint *)g_malloc0(GMSH_ELEMENT_SIZE*(*nt)*sizeof(gint)) ;
  *uvt = (gdouble *)g_malloc0(6*(*nt)*sizeof(gdouble)) ;

  ip = 0 ;
  it = 0 ;
  offp = 0 ;
  for ( s = 0 ; s < nsurf ; s ++ ) {
    for ( j = 0 ; j < nx[s] ; j ++ ) {
      (*p)[ip*GMSH_NODE_SIZE+0] =  x[s][3*j+0] ;
      (*p)[ip*GMSH_NODE_SIZE+1] =  x[s][3*j+1] ;
      (*p)[ip*GMSH_NODE_SIZE+2] =  x[s][3*j+2] ;
      (*p)[ip*GMSH_NODE_SIZE+3] = uv[s][2*j+0] ;
      (*p)[ip*GMSH_NODE_SIZE+4] = uv[s][2*j+1] ;
      
      (*pt)[ip*GMSH_NODE_DATA_SIZE+0] = ptags[s][j] ;
      (*pt)[ip*GMSH_NODE_DATA_SIZE+1] = ip ;
      (*pt)[ip*GMSH_NODE_DATA_SIZE+2] = s ;
      ip ++ ;
    }
    for ( j = 0 ; j < etags_n[s][0] ; j ++ ) {
      (*t)[it*GMSH_ELEMENT_SIZE+3] = s ; /* tags[2*s+1] ; */
      for ( k = 0 ; k < 3 ; k ++ ) {
	(*t)[it*GMSH_ELEMENT_SIZE+k] = entags[s][0][3*j+k] ;
	i = offp + find_tag(ptags[s], ntags[s], entags[s][0][3*j+k]) ;
	(*uvt)[6*it+2*k+0] = (*p)[GMSH_NODE_SIZE*i+3] ;
	(*uvt)[6*it+2*k+1] = (*p)[GMSH_NODE_SIZE*i+4] ;	
      }
      it ++ ;
    }
    gmshFree(etags_n[s]) ;
    gmshFree(etypes[s]) ;
    gmshFree(entags[s]) ;
    offp += nx[s] ;
  }
  /*sort the data list to use as a lookup table*/
  qsort((*pt), (*np), GMSH_NODE_DATA_SIZE*sizeof(gint), compare_tags) ;

  /*tie the element node indices to the point list entries*/
  for ( i = 0 ; i < (*nt) ; i ++ ) {
    for ( j = 0 ; j < 3 ; j ++ ) {
      node = bsearch(&((*t)[i*GMSH_ELEMENT_SIZE+j]), (*pt), (*np),
		     GMSH_NODE_DATA_SIZE*sizeof(gint), compare_tags) ;
      g_assert(node[0] == (*t)[i*GMSH_ELEMENT_SIZE+j]) ;
      (*t)[i*GMSH_ELEMENT_SIZE+j] = node[1] ;
    }
  }

  gmshFree(ptags) ; gmshFree(x) ; gmshFree(uv) ;
  g_free(etypes) ; g_free(entags) ;

  g_free(etags_n) ;
  
  return ;
}

static void shapefunc(gdouble s, gdouble t, gdouble L[])

{
  L[0] = 1.0 - s - t ;
  L[1] =       s     ;
  L[2] =           t ;
    
  return ;
}

static void element_interp(gint *tags, gint *tri, gdouble *uvtri,
			   gdouble *s, gint sstr,
			   gdouble *t, gint tstr, gint nst,
			   gdouble *xi, gint xstr)

{
  gdouble *buf, uv[1024], L[3] ;
  gint i, surf, ierr ;
  gsize n ;
  
  g_assert(nst <= 512) ;
  surf = tri[3] ;

  for ( i = 0 ; i < nst ; i ++ ) {
    shapefunc(s[i*sstr], t[i*tstr], L) ;
    uv[2*i+0] = uvtri[2*0+0]*L[0] + uvtri[2*1+0]*L[1] + uvtri[2*2+0]*L[2] ;
    uv[2*i+1] = uvtri[2*0+1]*L[0] + uvtri[2*1+1]*L[1] + uvtri[2*2+1]*L[2] ;
  }
  
  gmshModelGetValue(tags[2*surf+0], tags[2*surf+1], uv, 2*nst, &buf, &n,
		    &ierr) ;
  g_assert(ierr == 0) ;
  
  /* memcpy(xi, buf, 3*nst*sizeof(gdouble)) ; */

  for ( i = 0 ; i < nst ; i ++ ) {
    xi[i*xstr+0] = buf[3*i+0] ;
    xi[i*xstr+1] = buf[3*i+1] ;
    xi[i*xstr+2] = buf[3*i+2] ;
  }
  
  gmshFree(buf) ;
  
  return ;
}

static void add_gmsh_grid_patch(gint *tags, gint *tri, gdouble *uvtri,
				gdouble *st, gint nst,
				nbi_surface_t *s)
{
  gint i, np, nn, i3 = 3, xstr ;
  gdouble *x, work[3*453], K[454*454], ci[453*453], N ;
  gdouble al, bt ;

  np = nbi_surface_patch_number(s) ;
  nn = nbi_surface_node_number(s) ;
  nbi_surface_patch_node(s, np) = nn ;
  nbi_surface_patch_node_number(s, np) = nst ;

  N = sqt_koornwinder_interp_matrix(&(st[0]), 3, &(st[1]), 3, &(st[2]), 3,
				    nst, K) ;

  /*generate the mesh nodes*/
  al = 1.0 ; bt = 0.0 ; xstr = NBI_SURFACE_NODE_LENGTH ;
  x = (gdouble *)nbi_surface_node(s, nn) ;
  element_interp(tags, tri, uvtri, &(st[0]), 3, &(st[1]), 3, nst, x, xstr) ;
  blaswrap_dgemm(FALSE, FALSE, nst, i3, nst, al, K, nst, x, xstr, bt, ci, i3) ;
  
  for ( i = 0 ; i < nst ; i ++ ) {
    x = (gdouble *)nbi_surface_node(s, nn+i) ;
    sqt_element_interp(ci, nst, N, st[3*i+0], st[3*i+1],
		       &(x[0]), &(x[3]), &(x[6]), NULL, work) ;
    x[6] *= st[3*i+2] ;
  }

  nbi_surface_node_number(s) += nst ;
  nbi_surface_patch_number(s) ++ ;

  return ;
}

/** 
 * @ingroup surfaces
 *
 * Generate an ::nbi_surface_t from a GMSH .geo file (experimental,
 * but has been reliable up to now)
 * 
 * @param file name of GMSH file;
 * @param nq number of quadrature points per patch.
 * 
 * @return pointer to newly allocated ::nbi_surface_t containing
 * surface generated from \a file.
 */

nbi_surface_t *nbi_gmsh_mesh(gchar *file, gint nq)

{
  gint ierr, *tags, npts, *tri, ntri, *pdata, t ;
  gsize ntags ;
  gdouble *pts, *uvtri ;
  gdouble *st  ;
  gint order ;
  nbi_surface_t *s ;

  gmshOpen(file, &ierr) ;
  gmshModelGeoSynchronize(&ierr);
  gmshModelGetEntities(&tags, &ntags, 2, &ierr) ;
  gmshModelMeshGenerate(2, &ierr) ;

  gmesh_get(&pts, &pdata, &npts, &tri, &uvtri, &ntri) ;

  s = nbi_surface_alloc(nq*ntri, ntri) ;
  sqt_quadrature_select(nq, &st, &order) ;

  for ( t = 0 ; t < ntri ; t ++ ) {
    add_gmsh_grid_patch(tags,
			&(tri[t*GMSH_ELEMENT_SIZE]), &(uvtri[6*t]),
			st, nq, s) ;
  }
  
  return s ;
}
