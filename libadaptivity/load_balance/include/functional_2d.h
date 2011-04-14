/*
 *      Copyright (c) 2004-2006 by Gerard Gorman
 *      Copyright (c) 2006- Imperial College London
 *      See COPYING file for copying and redistribution conditions.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; version 2 of the License.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      Contact info: gerard.j.gorman@gmail.com/g.gorman@imperial.ac.uk
 */
#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H
#include <cfloat>
#include <cmath>

#include "samtypes.h"

class functional_2d{
 public:
  functional_2d();
  ~functional_2d();
    
  samfloat_t standard(samfloat_t x0, samfloat_t y0, const samfloat_t *m0,
		  samfloat_t x1, samfloat_t y1, const samfloat_t *m1,
		  samfloat_t x2, samfloat_t y2, const samfloat_t *m2);
  
  samfloat_t oddy(samfloat_t x0, samfloat_t y0, const samfloat_t *m0,
	      samfloat_t x1, samfloat_t y1, const samfloat_t *m1,
	      samfloat_t x2, samfloat_t y2, const samfloat_t *m2);
  
 private:
  samfloat_t dot(samfloat_t x0, samfloat_t y0,
	     samfloat_t x1, samfloat_t y1);

  samfloat_t len2(samfloat_t x0, samfloat_t y0,
	      samfloat_t x1, samfloat_t y1);
  
  samfloat_t GetInRadius(samfloat_t x0, samfloat_t y0,
		     samfloat_t x1, samfloat_t y1,
		     samfloat_t x2, samfloat_t y2);
  
  samfloat_t GetShape(samfloat_t x0, samfloat_t y0,
		  samfloat_t x1, samfloat_t y1,
		  samfloat_t x2, samfloat_t y2);
  
  samfloat_t GetSize(samfloat_t x0, samfloat_t y0,
		 samfloat_t x1, samfloat_t y1,
		 samfloat_t x2, samfloat_t y2);

  samfloat_t GetArea(samfloat_t r0_x, samfloat_t r0_y,
					  samfloat_t r1_x, samfloat_t r1_y,
					  samfloat_t r2_x, samfloat_t r2_y);
  
  const samfloat_t alpha;
  samfloat_t m[3];
};
#endif
