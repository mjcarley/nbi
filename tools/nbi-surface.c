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

#include "nbi-private.h"

GTimer *timer ;
gchar *progname ;

typedef nbi_surface_t *(*geometry_function)(gdouble *, gint *, gint) ;
nbi_surface_t *geometry_ellipsoid_ico(gdouble argd[], gint argi[], gint nq) ;
nbi_surface_t *geometry_stellarator(gdouble argd[], gint argi[], gint nq) ;

gpointer geometries[] = {"ellipsoid-ico",
			 geometry_ellipsoid_ico,
			 "subdivision depth",
			 "semi-axis x, semi-axis y, semi-axis z",
			 "stellarator",
			 geometry_stellarator,
			 "",
			 "",
			 NULL} ;

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

static void print_help_text(FILE *output)

{
  fprintf(output,
	  "Usage:\n\n"
	  "  %s <options>\n\n",
	  progname) ;

  fprintf(output,
	  "Options:\n\n"
	  "  -h print this message and exit\n"
#ifdef HAVE_AGG
	  "  -a # AGG input file\n"
#endif /*HAVE_AGG*/
	  "  -d # append a real argument for geometry specification\n"
	  "  -g # select a geometry\n"
	  "  -G list available built-in geometries\n"
	  "  -i # append an integer argument for geometry specification\n"
	  "  -o # output file\n"
	  "  -q # number of quadrature points per surface patch\n"
	  ) ;
  return ;
}

gint main(gint argc, gchar **argv)

{
  geometry_function gfunc ;
  gchar ch, *opfile, *aggfile ;
  gdouble argd[16] ;
  gint argi[16], nq, argci, argcd ;
  nbi_surface_t *s ;
  FILE *output, *input ;
  
  opfile = NULL ;
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  gfunc = NULL ; aggfile = NULL ;
  s = NULL ;
  nq = 7 ;
  
  argd[0] = 1.0 ; argd[1] = 1.0 ; argd[2] = 1.0 ;
  argi[0] = 1 ;

  argci = argcd = 0 ;
  
  while ( (ch = getopt(argc, argv, "ha:d:Gg:i:o:q:")) != EOF ) {
    switch (ch ) {
    default: g_assert_not_reached() ; break ;
    case 'a':
#ifdef HAVE_AGG
      if ( gfunc != NULL ) {
	fprintf(stderr, "%s: use -a or -g but not both\n", progname) ;
	return 1 ;
      }
      aggfile = g_strdup(optarg) ;
#else /*HAVE_AGG*/
      fprintf(stderr,
	      "AGG geometry support not implemented: install AGG and "
	      "recompile\n") ;
      exit(1) ;
#endif /*HAVE_AGG*/
      break ;
    case 'h': print_help_text(stderr) ; return 0 ; break ;
    case 'd': argd[argcd] = atof(optarg) ; argcd ++ ; break ;
    case 'g':
      if ( aggfile != NULL ) {
	fprintf(stderr, "%s: use -a or -g but not both\n", progname) ;
	return 1 ;
      }
      gfunc = parse_geometry(optarg) ; break ;
    case 'G':
      list_geometries(stderr) ;
      return 0 ;
      break ;
    case 'i': argi[argci] = atoi(optarg) ; argci ++ ; break ;
    case 'o': opfile = g_strdup(optarg) ; break ;
    case 'q': nq = atoi(optarg) ; break ;
    }
  }

  if ( aggfile == NULL && gfunc == NULL ) {
    fprintf(stderr, "%s: no geometry specified\n",
	    progname) ;

    return 1 ;
  }

#ifdef HAVE_AGG
  if ( aggfile != NULL ) {
    gint fid ;
    input = fopen(aggfile, "r") ;
    if ( input == NULL ) {
      fprintf(stderr, "%s: cannot open file %s\n", progname, aggfile) ;
      return 1 ;
    }
    fid = fileno(input) ;
    s = nbi_agg_mesh(fid, nq) ;
    fclose(input) ;
    close(fid) ;
  }
#endif /*HAVE_AGG*/

  
  if ( gfunc != NULL ) {
    s = gfunc(argd, argi, nq) ;
  }

  if ( s != NULL ) {
    fprintf(stderr,
	    "  patches: %d\n"
	    "  nodes:   %d\n",
	    nbi_surface_patch_number(s), nbi_surface_node_number(s)) ;
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
