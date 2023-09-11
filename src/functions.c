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

#define COMPLEX_MUL_REAL(_a,_b) (((_a)[0])*((_b)[0]) - ((_a)[1])*((_b)[1]))
#define COMPLEX_MUL_IMAG(_a,_b) (((_a)[1])*((_b)[0]) + ((_a)[0])*((_b)[1]))

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

static void gfunc_helmholtz(gdouble k, gdouble x, gdouble y, gdouble z,
			    gdouble *R, gdouble *G)

/*
 * Helmholtz Green's function for wavenumber k
 */
  
{
  *R = sqrt(x*x + y*y + z*z) ;
  G[0] = cos(k*(*R))/4.0/M_PI/(*R) ;
  G[1] = sin(k*(*R))/4.0/M_PI/(*R) ;
  
  return ;
}

static void dgfunc_helmholtz(gdouble k, gdouble x, gdouble y, gdouble z,
			     gdouble *dR, gdouble *dG)

/*
 * Helmholtz Green's function and its gradient for wavenumber k
 */

{
  gdouble jkR[2] ;
  
  dR[0] = sqrt(x*x + y*y + z*z) ;
  dG[0] = cos(k*dR[0])/4.0/M_PI/dR[0] ;
  dG[1] = sin(k*dR[0])/4.0/M_PI/dR[0] ;

  jkR[0] = -1 ; jkR[1] = k*dR[0] ;

  dG[2] = dG[4] = dG[6] = (dG[0]*jkR[0] - dG[1]*jkR[1])/dR[0] ;
  dG[3] = dG[5] = dG[7] = (dG[0]*jkR[1] + dG[1]*jkR[0])/dR[0] ;

  dR[1] = x/dR[0] ; dR[2] = y/dR[0] ; dR[3] = z/dR[0] ; 

  dG[2] *= dR[1] ; dG[3] *= dR[1] ;
  dG[4] *= dR[2] ; dG[5] *= dR[2] ;
  dG[6] *= dR[3] ; dG[7] *= dR[3] ;
  
  return ;
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
  gdouble dG[8], dR[4] ;

  dgfunc_helmholtz(k, x, y, z, dR, dG) ;
  
  return dG[2]*nx + dG[4]*ny + dG[6]*nz ;
}

gdouble nbi_function_gfunc_helmholtz_dG_imag(gdouble k,
					     gdouble x, gdouble y, gdouble z,
					     gdouble nx, gdouble ny, gdouble nz)

{
  gdouble dG[8], dR[4] ;

  dgfunc_helmholtz(k, x, y, z, dR, dG) ;
  
  return dG[3]*nx + dG[5]*ny + dG[7]*nz ;
}

static void gfunc_helmholtz_ring(gdouble a, gdouble n, gdouble k,
				 gdouble x, gdouble y, gdouble z,
				 gint ngp, gdouble *G)

{
  gdouble E[2], R, g[2], th1, y1, z1 ;
  gint i ;

  G[0] = G[1] = 0 ;
  for ( i = 0 ; i < ngp ; i ++ ) {
    th1 = 2.0*M_PI*i/ngp ;
    E[0] = cos(n*th1) ; E[1] = sin(n*th1) ;
    y1 = a*cos(th1) ; z1 = a*sin(th1) ;
    gfunc_helmholtz(k, x, y-y1, z-z1, &R, g) ;
    G[0] += COMPLEX_MUL_REAL(E, g) ;
    G[1] += COMPLEX_MUL_IMAG(E, g) ;
  }
  
  G[0] *= 2.0*M_PI/ngp ;
  G[1] *= 2.0*M_PI/ngp ;

  return ;
}

static void gfunc_helmholtz_ring_gradient(gdouble a, gdouble n, gdouble k,
					  gdouble x, gdouble y, gdouble z,
					  gint ngp, gdouble *dG)

{
  gdouble E[2], dR[4], dg[8], th1, y1, z1, sc ;
  gint i ;

  memset(dG, 0, 8*sizeof(gdouble)) ;
  sc = 2.0*M_PI/ngp ;  
  for ( i = 0 ; i < ngp ; i ++ ) {
    th1 = 2.0*M_PI*i/ngp ;
    E[0] = cos(n*th1) ; E[1] = sin(n*th1) ;
    y1 = a*cos(th1) ; z1 = a*sin(th1) ;
    dgfunc_helmholtz(k, x, y-y1, z-z1, dR, dg) ;
    dG[0] += COMPLEX_MUL_REAL(E, &(dg[0])) ;
    dG[1] += COMPLEX_MUL_IMAG(E, &(dg[0])) ;
    dG[2] += COMPLEX_MUL_REAL(E, &(dg[2])) ;
    dG[3] += COMPLEX_MUL_IMAG(E, &(dg[2])) ;
    dG[4] += COMPLEX_MUL_REAL(E, &(dg[4])) ;
    dG[5] += COMPLEX_MUL_IMAG(E, &(dg[4])) ;
    dG[6] += COMPLEX_MUL_REAL(E, &(dg[6])) ;
    dG[7] += COMPLEX_MUL_IMAG(E, &(dg[6])) ;
  }

  dG[0] *= sc ; dG[1] *= sc ; dG[2] *= sc ; dG[3] *= sc ;
  dG[4] *= sc ; dG[5] *= sc ; dG[6] *= sc ; dG[7] *= sc ;
  
  return ;
}

gdouble nbi_function_gfunc_helmholtz_ring_real(gdouble a, gdouble n, gdouble k,
					       gdouble x, gdouble y, gdouble z)

{
  gdouble G[2] ;
  gint ngp ;

  ngp = 1024 ;

  gfunc_helmholtz_ring(a, n, k, x, y, z, ngp, G) ;
  
  return G[0] ;
}
					       
gdouble nbi_function_gfunc_helmholtz_ring_imag(gdouble a, gdouble n, gdouble k,
					       gdouble x, gdouble y, gdouble z)

{
  gdouble G[2] ;
  gint ngp ;

  ngp = 1024 ;
  gfunc_helmholtz_ring(a, n, k, x, y, z, ngp, G) ;
  
  return G[1] ;
}
					       
gdouble nbi_function_gfunc_helmholtz_ring_real_x(gdouble a, gdouble n,
						 gdouble k,
						 gdouble x, gdouble y,
						 gdouble z)

{
  gdouble G[8] ;
  gint ngp ;

  ngp = 1024 ;
  gfunc_helmholtz_ring_gradient(a, n, k, x, y, z, ngp, G) ;
  
  return G[2] ;
}

gdouble nbi_function_gfunc_helmholtz_ring_real_y(gdouble a, gdouble n,
						 gdouble k,
						 gdouble x, gdouble y,
						 gdouble z)

{
  gdouble G[8] ;
  gint ngp ;

  ngp = 1024 ;
  gfunc_helmholtz_ring_gradient(a, n, k, x, y, z, ngp, G) ;
  
  return G[4] ;
}

gdouble nbi_function_gfunc_helmholtz_ring_real_z(gdouble a, gdouble n,
						 gdouble k,
						 gdouble x, gdouble y,
						 gdouble z)

{
  gdouble G[8] ;
  gint ngp ;

  ngp = 1024 ;
  gfunc_helmholtz_ring_gradient(a, n, k, x, y, z, ngp, G) ;
  
  return G[6] ;
}
					       
gdouble nbi_function_gfunc_helmholtz_ring_imag_x(gdouble a, gdouble n,
						 gdouble k,
						 gdouble x, gdouble y,
						 gdouble z)

{
  gdouble G[8] ;
  gint ngp ;

  ngp = 1024 ;
  gfunc_helmholtz_ring_gradient(a, n, k, x, y, z, ngp, G) ;
  
  return G[3] ;
}

gdouble nbi_function_gfunc_helmholtz_ring_imag_y(gdouble a, gdouble n,
						 gdouble k,
						 gdouble x, gdouble y,
						 gdouble z)

{
  gdouble G[8] ;
  gint ngp ;

  ngp = 1024 ;
  gfunc_helmholtz_ring_gradient(a, n, k, x, y, z, ngp, G) ;
  
  return G[5] ;
}

gdouble nbi_function_gfunc_helmholtz_ring_imag_z(gdouble a, gdouble n,
						 gdouble k,
						 gdouble x, gdouble y,
						 gdouble z)

{
  gdouble G[8] ;
  gint ngp ;

  ngp = 1024 ;
  gfunc_helmholtz_ring_gradient(a, n, k, x, y, z, ngp, G) ;
  
  return G[7] ;
}

static void spherical_scattered(gdouble k, gdouble a, gdouble r, gdouble th,
				gint M, gdouble *p)

/*
 * scattered pressure from sphere of radius a due to incident wave
 * \exp[j k z]
 */

{
  gdouble ha[258], hr[258], P[129], ka, kr, C, jm[2], Am[2], An[2], Ad[2] ;
  gint m ;
  
  g_assert(M < 128) ;

  /*generate Bessel functions and Legendre polynomials*/
  ka = k*a ; kr = k*r ;
  ha[2*0+0] = sin(ka)/ka ; ha[2*0+1] = -cos(ka)/ka ;
  ha[2*1+0] =  sin(ka)/ka/ka - cos(ka)/ka ;
  ha[2*1+1] = -cos(ka)/ka/ka - sin(ka)/ka ;
  hr[2*0+0] = sin(kr)/kr ; hr[2*0+1] = -cos(kr)/kr ;
  hr[2*1+0] =  sin(kr)/kr/kr - cos(kr)/kr ;
  hr[2*1+1] = -cos(kr)/kr/kr - sin(kr)/kr ;
  C = cos(th) ;
  P[0] = 1.0 ;
  P[1] = C ;

  for ( m = 0 ; m <= M-1 ; m ++ ) {
    P[m+2] = ((2*m+3)*C*P[m+1] - (m+1)*P[m])/(m+2) ;
    ha[2*(m+2)+0] = (2*m+3)/ka*ha[2*(m+1)+0] - ha[2*m+0] ;
    ha[2*(m+2)+1] = (2*m+3)/ka*ha[2*(m+1)+1] - ha[2*m+1] ;
    hr[2*(m+2)+0] = (2*m+3)/kr*hr[2*(m+1)+0] - hr[2*m+0] ;
    hr[2*(m+2)+1] = (2*m+3)/kr*hr[2*(m+1)+1] - hr[2*m+1] ;
  }

  p[0] = p[1] = 0 ;

  jm[0] = 1 ; jm[1] = 0 ;
  for ( m = 0 ; m <= M-1 ; m ++ ) {
    Ad[0] = m*ha[2*(m-1)+0] - (m+1)*ha[2*(m+1)+0] ;
    Ad[1] = m*ha[2*(m-1)+1] - (m+1)*ha[2*(m+1)+1] ;

    An[0] =  Ad[0]*Ad[0] ; An[1] = -Ad[0]*Ad[1] ;
    An[0] /= Ad[0]*Ad[0] + Ad[1]*Ad[1] ;
    An[1] /= Ad[0]*Ad[0] + Ad[1]*Ad[1] ;
    
    Am[0] = COMPLEX_MUL_REAL(jm, An)*(2*m+1) ;
    Am[1] = COMPLEX_MUL_IMAG(jm, An)*(2*m+1) ;
    p[0] -= COMPLEX_MUL_REAL(Am, &(hr[2*m]))*P[m] ;
    p[1] -= COMPLEX_MUL_IMAG(Am, &(hr[2*m]))*P[m] ;

    ka = jm[0] ; jm[0] = jm[1] ; jm[1] = ka ; jm[0] *= -1 ;
  }

  return ;
}

gdouble nbi_function_sphere_scattered_r(gdouble a, gdouble k,
					gdouble x, gdouble y,
					gdouble z)

{
  gdouble p[2], th, r ;
  gint M ;

  M = 8 ;
  r = sqrt(x*x + y*y + z*z) ;
  if ( z > r ) {
    th = 0 ;
  } else {
    if ( z < -r ) {
      th = M_PI ;
    } else {
      th = acos(z/r) ;
    }
  }

  spherical_scattered(k, a, r, th, M, p) ;
  
  return p[0] ;
}

gdouble nbi_function_sphere_scattered_i(gdouble a, gdouble k,
					gdouble x, gdouble y,
					gdouble z)

{
  gdouble p[2], th, r ;
  gint M ;

  M = 8 ;
  r = sqrt(x*x + y*y + z*z) ;
  if ( z > r ) {
    th = 0 ;
  } else {
    if ( z < -r ) {
      th = M_PI ;
    } else {
      th = acos(z/r) ;
    }
  }

  spherical_scattered(k, a, r, th, M, p) ;
  
  return p[1] ;
}
