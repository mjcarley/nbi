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

#define NBI_LOCAL_CUTOFF_RADIUS 1e-6

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

#define NBI_THREAD_DATA_SIZE         8
#define NBI_THREAD_DATA_MATRIX       0
#define NBI_THREAD_DATA_INT          1
#define NBI_THREAD_DATA_INT_POINTER  2
#define NBI_THREAD_DATA_REAL         3
#define NBI_THREAD_DATA_REAL_POINTER 4

#define NBI_THREAD_DATA_INT_SIZE     16
#define NBI_THREAD_DATA_INT_PSTR     0
#define NBI_THREAD_DATA_INT_NSTR     1
#define NBI_THREAD_DATA_INT_FSTR     2
#define NBI_THREAD_DATA_INT_OFFSET   3

#define NBI_THREAD_DATA_REAL_SIZE    8
#define NBI_THREAD_DATA_REAL_WT1     0
#define NBI_THREAD_DATA_REAL_WT2     2
#define NBI_THREAD_DATA_REAL_SIGN    4
/* #define NBI_THREAD_DATA_REAL_ */

#define NBI_THREAD_DATA_REAL_PTR_SIZE 8
#define NBI_THREAD_DATA_REAL_PTR_P    0
#define NBI_THREAD_DATA_REAL_PTR_PN   1
#define NBI_THREAD_DATA_REAL_PTR_F    2

#define NBI_THREAD_MAIN_DATA_SIZE    5
#define NBI_THREAD_MAIN_DATA_THREAD  0
#define NBI_THREAD_MAIN_DATA_NTHREAD 1
#define NBI_THREAD_MAIN_DATA_DATA    2
#define NBI_THREAD_MAIN_DATA_WORK    3
#define NBI_THREAD_MAIN_DATA_CONV    4

#define NBI_THREAD_NUMBER_MAX     8

#define nbi_vector_scalar(NBI_A,NBI_B)				\
  (((NBI_A)[0])*((NBI_B)[0])+					\
   ((NBI_A)[1])*((NBI_B)[1])+					\
   ((NBI_A)[2])*((NBI_B)[2]))
#define nbi_vector_length2(NBI_A)				\
  (((NBI_A)[0])*((NBI_A)[0]) +					\
   ((NBI_A)[1])*((NBI_A)[1]) +					\
   ((NBI_A)[2])*((NBI_A)[2]))
#define nbi_vector_length(NBI_A) (SQRT(nbi_vector_length2(NBI_A)))
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


gdouble nbi_function_gfunc_laplace_G(gdouble x, gdouble y, gdouble z) ;
gdouble nbi_function_gfunc_laplace_dG(gdouble x, gdouble y, gdouble z,
				      gdouble nx, gdouble ny, gdouble nz) ;

#endif /*NBI_PRIVATE_H_INCLUDED*/
