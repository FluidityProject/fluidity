! Copyright (C) 2006 Imperial College London and others.
! 
! Please see the AUTHORS file in the main source directory for a full list
! of copyright holders.
! 
! Adrian Umpleby
! Applied Modelling and Computation Group
! Department of Earth Science and Engineering
! Imperial College London
! 
! adrian@Imperial.ac.uk
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA

#include "confdefs.h"
#include "ewrite.h"

! #define FORTRAN_DISALLOWS_LONG_LINES

#ifndef __FILE__
#error __FILE__ does not work
#endif

#ifndef __LINE__
#error __LINE__ does not work
#endif

#define adabort(X) call adabort_pinpoint(X, __FILE__, __LINE__)
#define adabort(X) call adabort_pinpoint(X, __FILE__, __LINE__)

#ifdef NDEBUG
#define ASSERT(X)
#else
#ifdef FORTRAN_DISALLOWS_LONG_LINES
#define ASSERT(X) IF(.NOT.(X)) adabort('Failed assertion ')
#else
#define ASSERT(X) IF(.NOT.(X)) adabort('Failed assertion '//'X')
#endif
#endif
#define assert(X) ASSERT(X)
