/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    C.Pain@Imperial.ac.uk
    
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation,
    version 2.1 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
    USA
*/
#include "confdefs.h"

#include <iostream>
#include <string.h>

#include "InterpolationTri3Tri3.h"

#ifdef DOUBLEP
typedef double gigreal_t;
#else
typedef float gigreal_t;
#endif

extern "C" {
#define fltri3totri3_fc F77_FUNC(fltri3totri3, FLTRI3TOTRI3)

  void fltri3totri3_fc(int *NNode1,  int *NElems1, gigreal_t *X1, gigreal_t *Y1,
           int *ENLIST1, gigreal_t *RMEM1,  int *FIELDS1, int *NFIELDS, 
           int *NNode2,  int *NElems2, gigreal_t *X2, gigreal_t *Y2,
           int *ENLIST2, gigreal_t *RMEM2,  int *FIELDS2, int *IERROR){
    // reset error
    *IERROR = 0;  
    InterpolationTri3Tri3<gigreal_t> interpolator(1);
    interpolator.SetSource(ENLIST1, X1, Y1, *NNode1, *NElems1);
    interpolator.SetDest(ENLIST2, X2, Y2, *NNode2, *NElems2); 
    
    for(int i=0; i<*NFIELDS; i++) 
      interpolator.Interpolate(&(RMEM1[FIELDS1[i]-1]), &(RMEM2[FIELDS2[i]-1]));
    
    return;  
  }

}
