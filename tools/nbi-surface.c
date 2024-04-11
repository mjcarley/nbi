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

#ifdef HAVE_WBFMM
#include <wbfmm.h>
#endif /*HAVE_WBFMM*/

#ifdef HAVE_SQT
#include <sqt.h>
#endif /*HAVE_SQT*/

#ifdef HAVE_BLASWRAP
#include <blaswrap.h>
#endif /*HAVE_BLASWRAP*/

#ifdef HAVE_AGG
#include <agg.h>
#endif /*HAVE_AGG*/

#include "nbi-private.h"

#ifdef HAVE_GMSHC_H
#include <gmshc.h>
#endif /*HAVE_GMSHC_H*/

GTimer *timer ;
gchar *progname ;

/**
 * @page nbisurface Generating surfaces for NBI calculations
 *
 * @c nbi-surface generates and discretizes surface geometries and
 * outputs them in the NBI file format. 
 * 
 * @verbatim
 nbi-surface -h
 @endverbatim
 * gives a current list of command line options. The most important of 
 * these are:
 * - @c -q number of quadrature points per surface patch;
 * - @c -o output file name;
 * - @c -g built-in surface geometry;
 * - @c -G list available built-in geometries;
 * - @c -m GMSH geometry file name;
 * - @c -a AGG geometry file name;
 * - @c -i add an integer argument to pass to the built-in geometry generator;
 * - @c -f add a floating point argument to pass to the built-in geometry 
 * generator.
 *
 * Option @c -f, @c -m, or @c -a must be specified. 
 * 
 */

typedef nbi_surface_t *(*geometry_function)(gdouble *, gint *, gint) ;
nbi_surface_t *geometry_ellipsoid_ico(gdouble argd[], gint argi[], gint nq) ;
nbi_surface_t *geometry_stellarator(gdouble argd[], gint argi[], gint nq) ;
nbi_surface_t *geometry_grid_xy(gdouble argd[], gint argi[], gint nq) ;
nbi_surface_t *geometry_grid_yz(gdouble argd[], gint argi[], gint nq) ;
nbi_surface_t *geometry_grid_zx(gdouble argd[], gint argi[], gint nq) ;

gpointer geometries[] = {"ellipsoid-ico",
			 geometry_ellipsoid_ico,
			 "subdivision depth",
			 "semi-axis x, semi-axis y, semi-axis z",

			 "stellarator",
			 geometry_stellarator,
			 "",
			 "",

			 "grid-xy",
			 geometry_grid_xy,
			 "nu, nv",
			 "xmin, ymin, xmax, ymax, z",

			 "grid-yz",
			 geometry_grid_yz,
			 "nu, nv",
			 "ymin, zmin, ymax, zmax, x",

			 "grid-zx",
			 geometry_grid_zx,
			 "nu, nv",
			 "zmin, xmin, zmax, xmax, y",
			 NULL, NULL, NULL, NULL} ;

nbi_surface_t *geometry_ellipsoid_ico(gdouble argd[], gint argi[], gint nq)

{
  nbi_surface_t *s = NULL ;
  gint ico, np ;
  gdouble rx, ry, rz ;
  
  ico = argi[0] ;
  rx = argd[0] ; ry = argd[1] ; rz = argd[2] ; 

  np = 20*(1 << 2*ico) ;
  
  fprintf(stderr,
	  "%s: icosohedral ellipsoid\n"
	  "  subdivision depth: %d\n"
	  "  points per patch:  %d\n"
	  "  semi-axis x: %lg\n"
	  "  semi-axis y: %lg\n"
	  "  semi-axis z: %lg\n",
	  progname, ico, nq, rx, ry, rz) ;

  s = nbi_surface_alloc(np*nq, np) ;    
  nbi_geometry_ellipsoid_ico(s, rx, ry, rz, ico, nq) ;
  
  return s ;
}

nbi_surface_t *geometry_stellarator(gdouble argd[], gint argi[], gint nq)

{
  nbi_surface_t *s = NULL ;

  return s ;
}

nbi_surface_t *geometry_grid_xy(gdouble argd[], gint argi[], gint nq)

{
  nbi_surface_t *s = NULL ;
  gint nu, nv, i, xstr ;
  gdouble xmin, xmax, ymin, ymax, z, *x, u, v ;

  nu = argi[0] ; nv = argi[1] ;
  xmin = argd[0] ; ymin = argd[1] ; 
  xmax = argd[2] ; ymax = argd[3] ; 
  z = argd[4] ;
  
  s = nbi_surface_alloc(2*nu*nv*nq, 2*nu*nv) ;

  nbi_geometry_grid(s, 0, 1, nu, 0, 1, nv, nq) ;

  x = (NBI_REAL *)nbi_surface_node(s, 0) ;
  xstr = NBI_SURFACE_NODE_LENGTH ;

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    u = x[i*xstr+0] ; v = x[i*xstr+1] ;
    x[i*xstr+0] = xmin + (xmax - xmin)*u ;
    x[i*xstr+1] = ymin + (ymax - ymin)*v ;
    x[i*xstr+2] = z ;
    x[i*xstr+3] = 0 ; x[i*xstr+4] = 0 ; x[i*xstr+5] = 1 ; 
  }    
  
  return s ;
}

nbi_surface_t *geometry_grid_yz(gdouble argd[], gint argi[], gint nq)

{
  nbi_surface_t *s = NULL ;
  gint nu, nv, i, xstr ;
  gdouble ymin, ymax, zmin, zmax, xp, *x, u, v ;

  nu = argi[0] ; nv = argi[1] ;
  ymin = argd[0] ; zmin = argd[1] ; 
  ymax = argd[2] ; zmax = argd[3] ; 
  xp = argd[4] ;
  
  s = nbi_surface_alloc(2*nu*nv*nq, 2*nu*nv) ;

  nbi_geometry_grid(s, 0, 1, nu, 0, 1, nv, nq) ;

  x = (NBI_REAL *)nbi_surface_node(s, 0) ;
  xstr = NBI_SURFACE_NODE_LENGTH ;

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    u = x[i*xstr+0] ; v = x[i*xstr+1] ;
    x[i*xstr+0] = xp ;
    x[i*xstr+1] = ymin + (ymax - ymin)*u ;
    x[i*xstr+2] = zmin + (zmax - zmin)*v ;
    x[i*xstr+3] = 0 ; x[i*xstr+4] = 0 ; x[i*xstr+5] = 1 ; 
  }    
  
  return s ;
}

nbi_surface_t *geometry_grid_zx(gdouble argd[], gint argi[], gint nq)

{
  nbi_surface_t *s = NULL ;
  gint nu, nv, i, xstr ;
  gdouble zmin, zmax, xmin, xmax, y, *x, u, v ;

  nu = argi[0] ; nv = argi[1] ;
  zmin = argd[0] ; xmin = argd[1] ; 
  zmax = argd[2] ; xmax = argd[3] ; 
  y = argd[4] ;
  
  s = nbi_surface_alloc(2*nu*nv*nq, 2*nu*nv) ;

  nbi_geometry_grid(s, 0, 1, nu, 0, 1, nv, nq) ;

  x = (NBI_REAL *)nbi_surface_node(s, 0) ;
  xstr = NBI_SURFACE_NODE_LENGTH ;

  for ( i = 0 ; i < nbi_surface_node_number(s) ; i ++ ) {
    u = x[i*xstr+0] ; v = x[i*xstr+1] ;
    x[i*xstr+0] = xmin + (xmax - xmin)*u ;
    x[i*xstr+1] = y ;
    x[i*xstr+2] = zmin + (zmax - zmin)*v ;
    x[i*xstr+3] = 0 ; x[i*xstr+4] = 0 ; x[i*xstr+5] = 1 ; 
  }    
  
  return s ;
}

static geometry_function parse_geometry(gchar *str)

{
  gint i ;

  for ( i = 0 ; geometries[i] != NULL ; i += 4 ) {
    if ( strcmp(str, (gchar *)geometries[i]) == 0 )
      return geometries[i+1] ;
  }
  
  return NULL ;
}

static gint list_geometries(FILE *f)

{
  gint i ;
  
  for ( i = 0 ; geometries[i] != NULL ; i += 4 ) {
    fprintf(f, "%s\n", (gchar *)geometries[i]) ;
    fprintf(f, "  integer arguments: %s\n", (gchar *)geometries[i+2]) ;
    fprintf(f, "     real arguments: %s\n", (gchar *)geometries[i+3]) ;
  }

  return 0 ;
}

static void print_help_text(FILE *output, gint nqp)

{
  fprintf(output,
	  "%s: generate and discretize NBI surfaces\n\n"
	  "Usage:\n\n"
	  "  %s <options> > <output>\n\n",
	  progname, progname) ;

  fprintf(output,
	  "Options:\n\n"
	  "  -h print this message and exit\n"
#ifdef HAVE_AGG
	  "  -a # AGG geometry file\n"
#endif /*HAVE_AGG*/
	  "  -d # append a real argument for geometry specification\n"
	  "  -f write list of geometry formats handled, and exit\n"
	  "  -g # select a built-in geometry\n"
	  "  -G list available built-in geometries\n"
	  "  -i # append an integer argument for geometry specification\n"
#ifdef HAVE_LIBGMSH
	  "  -m # GMSH geometry file\n"
#endif /*HAVE_LIBGMSH*/
	  "  -o # output file\n"
	  "  -q # number of quadrature points per surface patch (%d)\n",
	  nqp
	  ) ;
  return ;
}

static void list_input_formats(FILE *f)

{
  fprintf(f, "geometry formats handled by %s\n\n", progname) ;

  fprintf(f, "  internal geometry generation (-G, -g)\n") ;
#ifdef HAVE_LIBGMSH
  fprintf(f, "  mesh based on GMSH geometry file (-m)\n") ;
#endif /*HAVE_LIBGMSH*/
  
  return ; 
}

gint main(gint argc, gchar **argv)

{
  geometry_function gfunc ;
  gchar ch, *opfile ;
#ifdef HAVE_LIBGMSH
  gchar *gmshfile ;
#endif /*HAVE_LIBGMSH*/
#ifdef HAVE_AGG
  gchar *aggfile ;
#endif /*HAVE_AGG*/
  gdouble argd[16] ;
  gint argi[16], nq, argci, argcd ;
  nbi_surface_t *s ;
  FILE *output ;
  
  opfile = NULL ;
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  gfunc = NULL ;
#ifdef HAVE_LIBGMSH
  gmshfile = NULL ;
#endif /*HAVE_LIBGMSH*/
#ifdef HAVE_AGG
  aggfile = NULL ;
#endif /*HAVE_AGG*/

  s = NULL ;
  nq = 7 ;
  
  argd[0] = 1.0 ; argd[1] = 1.0 ; argd[2] = 1.0 ;
  argi[0] = 1 ;

  argci = argcd = 0 ;
  
  while ( (ch = getopt(argc, argv, "ha:d:fGg:i:m:o:q:")) != EOF ) {
    switch (ch ) {
    default: g_assert_not_reached() ; break ;
    case 'h': print_help_text(stderr, nq) ; return 0 ; break ;
#ifndef HAVE_AGG
    case 'a':
      fprintf(stderr, "%s: AGG not supported\n", progname) ;
      return 1 ;
      break ;
#endif /*HAVE_AGG*/
#ifdef HAVE_AGG
    case 'a':
      if ( gfunc != NULL
#ifdef HAVE_LIBGMSH
	   || gmshfile != NULL
#endif /*HAVE_LIBGMSH*/	   
	   ) {
	fprintf(stderr, "%s: use only one geometry specification\n",
		progname) ;
	return 1 ;
      }
      aggfile = g_strdup(optarg) ; break ;
      break ;
#endif /*HAVE_AGG*/
    case 'd': argd[argcd] = atof(optarg) ; argcd ++ ; break ;
    case 'f': list_input_formats(stderr) ; return 0 ; break ;
    case 'g':
      if ( gmshfile != NULL ) {
	fprintf(stderr, "%s: use -g or -m but not both\n", progname) ;
	return 1 ;
      }
      gfunc = parse_geometry(optarg) ; break ;
    case 'G':
      list_geometries(stderr) ;
      return 0 ;
      break ;
    case 'i': argi[argci] = atoi(optarg) ; argci ++ ; break ;
#ifdef HAVE_LIBGMSH
    case 'm':
      if ( gfunc != NULL
#ifdef HAVE_AGG
	   || aggfile != NULL
#endif /*HAVE_AGG*/	   
	   ) {
	fprintf(stderr, "%s: use only one geometry specification\n",
		progname) ;
	return 1 ;
      }
      gmshfile = g_strdup(optarg) ; break ;
#endif /*HAVE_LIBGMSH*/
    case 'o': opfile = g_strdup(optarg) ; break ;
    case 'q': nq = atoi(optarg) ; break ;
    }
  }

  if ( gfunc == NULL
#ifdef HAVE_LIBGMSH
       && gmshfile == NULL
#endif /*HAVE_LIBGMSH*/
#ifdef HAVE_AGG
       && aggfile == NULL
#endif /*HAVE_AGG*/
       ) {
    fprintf(stderr, "%s: no geometry specified\n", progname) ;

    return 1 ;
  }

#ifdef HAVE_LIBGMSH
  gint gmsh_err ;
  if ( gmshfile != NULL ) {
    /* gmshInitialize(argc, argv, 1, 0, &gmsh_err) ; */
    /*don't pass the nbi options to GMSH (who knows what might happen?)*/
#if GMSH_API_VERSION_MINOR < 10
    gmshInitialize(0, NULL, 1,    &gmsh_err) ;
#else /*GMSH_API_VERSION_MINOR < 10*/
    gmshInitialize(0, NULL, 1, 0, &gmsh_err) ;
#endif /*GMSH_API_VERSION_MINOR < 10*/
    gmshOptionSetNumber("General.Verbosity", 0, &gmsh_err) ;

    s = nbi_gmsh_mesh(gmshfile, nq) ;
    
    gmshFinalize(&gmsh_err) ;
  }
#endif /*HAVE_LIBGMSH*/

#ifdef HAVE_AGG
  if ( aggfile != NULL ) {
    s = nbi_agg_mesh(aggfile, nq) ;
  }
#endif /*HAVE_AGG*/

  if ( gfunc != NULL ) {
    s = gfunc(argd, argi, nq) ;
  }

  if ( s != NULL ) {
    fprintf(stderr,
	    "%s: patches: %d\n"
	    "%s: nodes:   %d\n",
	    progname, nbi_surface_patch_number(s),
	    progname, nbi_surface_node_number(s)) ;
  }
  
  if ( opfile == NULL ) {
    output = stdout ;
  } else {
    output = fopen(opfile, "w") ;
    if ( output == NULL ) {
      fprintf(stderr, "%s: cannot open file %s\n", progname, opfile) ;
      exit(1) ;
    }
  }

  if ( s != NULL ) nbi_surface_write(s, output) ;

  if ( output != stdout ) fclose(output) ;
  
  return 0 ;
}
