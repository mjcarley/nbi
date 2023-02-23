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

gdouble nbi_function_gfunc_laplace_G(gdouble x, gdouble y, gdouble z)

{
  gdouble G, R ;

  R = x*x + y*y + z*z ;

  if ( R == 0.0 ) return G_MAXDOUBLE ;

  R = sqrt(R) ;

  G = 0.25*M_1_PI/R ;
  
  return G ;
}

gdouble nbi_function_gfunc_laplace_dG(gdouble x, gdouble y, gdouble z,
				     gdouble nx, gdouble ny, gdouble nz)

{
  gdouble dG, R ;

  R = x*x + y*y + z*z ;

  if ( R == 0.0 ) return G_MAXDOUBLE ;

  R = sqrt(R) ;

  dG = -0.25*M_1_PI*(x*nx + y*ny + z*nz)/R/R/R ;
  
  return dG ;
}

gdouble nbi_function_gfunc_helmholtz_G_real(gdouble k,
					    gdouble x, gdouble y, gdouble z)

{
  gdouble G, R ;

  R = x*x + y*y + z*z ;

  if ( R == 0.0 ) return G_MAXDOUBLE ;

  R = sqrt(R) ;

  G = 0.25*M_1_PI*cos(k*R)/R ;
  
  return G ;
}

gdouble nbi_function_gfunc_helmholtz_G_imag(gdouble k,
					    gdouble x, gdouble y, gdouble z)

{
  gdouble G, R ;

  R = x*x + y*y + z*z ;

  if ( R == 0.0 ) return G_MAXDOUBLE ;

  R = sqrt(R) ;

  G = 0.25*M_1_PI*sin(k*R)/R ;
  
  return G ;
}

gdouble nbi_function_gfunc_helmholtz_dG_real(gdouble k,
					     gdouble x, gdouble y, gdouble z,
					     gdouble nx, gdouble ny, gdouble nz)

{
  gdouble dG, R ;

  R = x*x + y*y + z*z ;

  if ( R == 0.0 ) return G_MAXDOUBLE ;

  R = sqrt(R) ;

  dG = 0.25*M_1_PI*(x*nx + y*ny + z*nz)/R/R/R*(-cos(k*R) - k*R*sin(k*R)) ;
  
  return dG ;
}

gdouble nbi_function_gfunc_helmholtz_dG_imag(gdouble k,
					     gdouble x, gdouble y, gdouble z,
					     gdouble nx, gdouble ny, gdouble nz)

{
  gdouble dG, R ;

  R = x*x + y*y + z*z ;

  if ( R == 0.0 ) return G_MAXDOUBLE ;

  R = sqrt(R) ;

  dG = 0.25*M_1_PI*(x*nx + y*ny + z*nz)/R/R/R*(-sin(k*R) + k*R*cos(k*R)) ;
  
  return dG ;
}
