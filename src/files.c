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

#include "nbi-private.h"

/**
 * @ingroup export
 *
 * @{
 * 
 */

/** 
 * Insert a string into the header of an output data file
 * 
 * @param header header string;
 * @param i index of insertion point in \a header;
 * @param len length of \a header;
 * @param str string to insert.
 * 
 * @return 0 on success.
 */

gint nbi_header_insert_string(gchar *header, gint i, gint len, gchar *str)

{
  gint slen, j ;
  
  if ( (slen = strlen(str)) > len )
    g_error("%s: string %s is too long (%d) for header (%d)", __FUNCTION__,
	    str, slen, len) ;

  if ( slen + i > NBI_HEADER_LENGTH )
    g_error("%s: not enough space in header for string ""%s"" of length %d "
	    "at position %d",  __FUNCTION__, str, slen, i) ;
  
  for ( j = 0 ; j < slen ; j ++ ) {
    header[i+j] = str[j] ;
  }
  
  return 0 ;
}

/** 
 * Initialize a header for an output file
 * 
 * @param header header string;
 * @param id identification of output file type;
 * @param version version number of output format;
 * @param type type of file;
 * @param format ASCII or binary.
 * 
 * @return 0 on success.
 */

gint nbi_header_init(gchar *header,
		     gchar *id,
		     gchar *version,
		     gchar *type,
		     gchar *format)

{
  gint i ;

  for ( i = 0 ; i < NBI_HEADER_LENGTH ; i ++ ) header[i] = ' ' ;
  header[NBI_HEADER_LENGTH] = '\0' ;

  nbi_header_insert_string(header, NBI_HEADER_ID,      3, id) ;
  nbi_header_insert_string(header, NBI_HEADER_VERSION, 3, version) ;
  nbi_header_insert_string(header, NBI_HEADER_TYPE,    3, type) ;
  nbi_header_insert_string(header, NBI_HEADER_FORMAT,  1, format) ;

  if ( format[0] != 'A' && format[0] != 'B' )
    g_error("%s: data format must be `A' or `B' (ASCII or binary)",
	    __FUNCTION__) ;

  return 0 ;
}

/** 
 * Read a header from input file
 * 
 * @param f input file stream;
 * @param header header string.
 * 
 * @return 0 on success.
 */

gint nbi_header_read(FILE *f, gchar header[])

{
  gchar c ;
  
  fread(header, sizeof(gchar), NBI_HEADER_LENGTH, f) ;
  fread(&c    , sizeof(gchar), 1,                 f) ;

  return 0 ;
}

/** 
 * Write a header to output
 * 
 * @param f output file stream;
 * @param header header string.
 * 
 * @return 0 on success.
 */

gint nbi_header_write(FILE *f, gchar header[])

{
  gchar c = '\n' ;

  fwrite(header, sizeof(gchar), NBI_HEADER_LENGTH, f) ;
  fwrite(&c    , sizeof(gchar), 1,                 f) ;

  return 0 ;
}

/** 
 * Problem type string to be used in file headers
 * 
 * @param p problem type;
 * 
 * @return string identifying \a p.
 */

gchar *nbi_problem_type_string(nbi_problem_t p)

{
  switch ( p ) {
  default: return NULL ; break ;
  case NBI_PROBLEM_UNDEFINED: return "UNDEFINED" ; break ;
  case NBI_PROBLEM_LAPLACE  : return "LAPLACE" ; break ;
  case NBI_PROBLEM_HELMHOLTZ: return "HELMHOLTZ" ; break ;
  }
  
  return NULL ;
}

/** 
 * Parse a problem identifier
 * 
 * @param p problem string, from a file header.
 * 
 * @return ::nbi_problem_t corresponding to \a p, or ::NBI_PROBLEM_UNDEFINED.
 */

nbi_problem_t nbi_problem_from_string(gchar *p)

{
  if ( strcmp(p, "UNDEFINED") == 0 ) return NBI_PROBLEM_UNDEFINED ;
  if ( strcmp(p, "LAPLACE") == 0   ) return NBI_PROBLEM_LAPLACE ;
  if ( strcmp(p, "HELMHOLTZ") == 0 ) return NBI_PROBLEM_HELMHOLTZ ;
  
  return NBI_PROBLEM_UNDEFINED ;
}


/**
 *
 * @}
 *
 */
