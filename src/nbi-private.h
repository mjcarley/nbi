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

#ifndef NBI_PRIVATE_H_INCLUDED
#define NBI_PRIVATE_H_INCLUDED

#ifdef NBI_SINGLE_PRECISION

#define NBI_REAL gfloat

#define NBI_FUNCTION_NAME(NBI_func) NBI_func##_f

#define SQRT(NBI_x) sqrtf((NBI_x))
#define CBRT(NBI_x) cbrtf((NBI_x))
#define SIN(NBI_x) sinf((NBI_x))
#define COS(NBI_x) cosf((NBI_x))
#define ACOS(NBI_x) acosf((NBI_x))
#define ATAN(NBI_x) atanf((NBI_x))
#define ATAN2(NBI_y,NBI_x) atan2f((NBI_y),(NBI_x))
#define LOG(NBI_x) logf((NBI_x))
#define EXP(NBI_x) expf(NBI_x)

#else

#define NBI_REAL gdouble

#define NBI_FUNCTION_NAME(NBI_func) NBI_func

#define SQRT(NBI_x) sqrt((NBI_x))
#define CBRT(NBI_x) cbrt((NBI_x))
#define SIN(NBI_x) sin((NBI_x))
#define COS(NBI_x) cos((NBI_x))
#define ACOS(NBI_x) acos((NBI_x))
#define ATAN(NBI_x) atan((NBI_x))
#define ATAN2(NBI_y,NBI_x) atan2((NBI_y),(NBI_x))
#define LOG(NBI_x) log((NBI_x))
#define EXP(NBI_x) exp(NBI_x)

#endif /*NBI_SINGLE_PRECISION*/

#define nbi_vector_scalar(NBI_A,NBI_B)				\
  (((NBI_A)[0])*((NBI_B)[0])+					\
   ((NBI_A)[1])*((NBI_B)[1])+					\
   ((NBI_A)[2])*((NBI_B)[2]))
#define nbi_vector_length(NBI_A)				\
  (SQRT(((NBI_A)[0])*((NBI_A)[0])+				\
	((NBI_A)[1])*((NBI_A)[1]) +				\
	((NBI_A)[2])*((NBI_A)[2])))
#define nbi_vector_cross(NBI_C,NBI_A,NBI_B)				\
  ((NBI_C)[0] = (NBI_A)[1]*(NBI_B)[2] - (NBI_A)[2]*(NBI_B)[1],		\
   (NBI_C)[1] = (NBI_A)[2]*(NBI_B)[0] - (NBI_A)[0]*(NBI_B)[2],		\
   (NBI_C)[2] = (NBI_A)[0]*(NBI_B)[1] - (NBI_A)[1]*(NBI_B)[0])
#define nbi_vector_distance2(NBI_A,NBI_B)		\
  ( ((NBI_A)[0]-(NBI_B)[0])*((NBI_A)[0]-(NBI_B)[0]) +	\
    ((NBI_A)[1]-(NBI_B)[1])*((NBI_A)[1]-(NBI_B)[1]) +	\
    ((NBI_A)[2]-(NBI_B)[2])*((NBI_A)[2]-(NBI_B)[2]) )

#define nbi_vector_distance(NBI_A,NBI_B)	\
  (SQRT((nbi_vector_distance2(NBI_A,NBI_B))))
#define nbi_vector_diff_scalar(NBI_A,NBI_B,NBI_C)			\
  (((NBI_A)[0]-(NBI_B)[0])*((NBI_C)[0]) +				\
   ((NBI_A)[1]-(NBI_B)[1])*((NBI_C)[1]) +				\
   ((NBI_A)[2]-(NBI_B)[2])*((NBI_C)[2]))

#define nbi_vector_diff(NBI_A,NBI_B,NBI_C)	\
  do {						\
    (NBI_A)[0] = (NBI_B)[0] - (NBI_C)[0] ;		\
    (NBI_A)[1] = (NBI_B)[1] - (NBI_C)[1] ;		\
    (NBI_A)[2] = (NBI_B)[2] - (NBI_C)[2] ;		\
  } while (0)

#endif /*NBI_PRIVATE_H_INCLUDED*/
