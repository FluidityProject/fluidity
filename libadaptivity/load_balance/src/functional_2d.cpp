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

#include "functional_2d.h"
#include <iostream>

using namespace std;

functional_2d::functional_2d(): alpha(0.5/sqrt(3.0)){}
functional_2d::~functional_2d(){}

samfloat_t functional_2d::dot(samfloat_t x0, samfloat_t y0,
		       samfloat_t x1, samfloat_t y1){
  return x0*(x1*m[0]+y1*m[1]) + y0*(x1*m[1]+y1*m[2]);
}

samfloat_t functional_2d::len2(samfloat_t x0, samfloat_t y0,
			samfloat_t x1, samfloat_t y1){
  samfloat_t vx = x0-x1;
  samfloat_t vy = y0-y1;
  return vx*vx*m[0] + 2*vx*vy*m[1] + vy*vy*m[2];
}

samfloat_t functional_2d::GetInRadius(samfloat_t x0, samfloat_t y0,
			       samfloat_t x1, samfloat_t y1,
			       samfloat_t x2, samfloat_t y2){
  samfloat_t r0 = sqrt(len2(x1, y1, x2, y2));
  samfloat_t r1 = sqrt(len2(x2, y2, x0, y0));
  samfloat_t r2 = sqrt(len2(x0, y0, x1, y1));
  
  samfloat_t s=0.5*(r0+r1+r2);
  return sqrt((s-r0)*(s-r1)*(s-r2)/s);  
}

samfloat_t functional_2d::GetArea(samfloat_t r0_x, samfloat_t r0_y,
					  samfloat_t r1_x, samfloat_t r1_y,
					  samfloat_t r2_x, samfloat_t r2_y){
  samfloat_t x1 = r1_x - r0_x;
  samfloat_t y1 = r1_y - r0_y;
  
  samfloat_t x2 = r2_x - r0_x;
  samfloat_t y2 = r2_y - r0_y;
  
  return -0.5*(x1*y2 - x2*y1);
}


samfloat_t functional_2d::GetShape(samfloat_t x0, samfloat_t y0,
			    samfloat_t x1, samfloat_t y1,
			    samfloat_t x2, samfloat_t y2){
  
  samfloat_t inradius = GetInRadius(x0, y0,
				x1, y1,
				x2, y2);
  
  samfloat_t shape=0.0;
  if(inradius>10e-6)
    shape = (alpha/inradius - 1)*(alpha/inradius - 1);
  else
    shape = 1.0e6;
  
  return shape;
}

samfloat_t functional_2d::GetSize(samfloat_t x0, samfloat_t y0,
			   samfloat_t x1, samfloat_t y1,
			   samfloat_t x2, samfloat_t y2){
  samfloat_t size=0.0;
  
  // edge a
  samfloat_t r0 = sqrt(len2(x0, y0, x1, y1));
  size += ((r0-1)*(r0-1));
  
  // edge b
  samfloat_t r1 = sqrt(len2(x0, y0, x2, y2));
  size += ((r1-1)*(r1-1));
  
  // edge c
  samfloat_t r2 = sqrt(len2(x1, y1, x2, y2));
  size += ((r2-1)*(r2-1));
  
  return size;
}

samfloat_t functional_2d::standard(samfloat_t x0, samfloat_t y0, const samfloat_t *m0,
                                   samfloat_t x1, samfloat_t y1, const samfloat_t *m1,
                                   samfloat_t x2, samfloat_t y2, const samfloat_t *m2){
  
  m[0] = (m0[0]+m1[0]+m2[0])/3.0;
  m[1] = (m0[1]+m1[1]+m2[1])/3.0;
  m[2] = (m0[2]+m1[2]+m2[2])/3.0;
  
  samfloat_t shape = GetShape(x0, y0,
                              x1, y1,
                              x2, y2);
  
  samfloat_t size = GetSize(x0, y0,
                            x1, y1,
                            x2, y2);
  
  samfloat_t functional = 10.0*size + shape;
  if(!finite(functional)){
    return 1.0e6;
  }

  return functional;
}

samfloat_t functional_2d::oddy(samfloat_t x0, samfloat_t y0, const samfloat_t *m0,
			samfloat_t x1, samfloat_t y1, const samfloat_t *m1,
			samfloat_t x2, samfloat_t y2, const samfloat_t *m2){

  m[0] = (m0[0]+m1[0]+m2[0])/3.0;
  m[1] = (m0[1]+m1[1]+m2[1])/3.0;
  m[2] = (m0[2]+m1[2]+m2[2])/3.0;
  
  if(GetArea(x0, y0, x1, y1, x2, y2)<0.0)
    return 1.0e6;
  
  samfloat_t l2[3];
  l2[0] = len2(x1, y1, x2, y2);
  l2[1] = len2(x2, y2, x0, y0);
  l2[2] = len2(x0, y0, x1, y1);
  
  samfloat_t b2[3];
  b2[0] = dot(x2-x1, y2-y1, x0-x2, y0-y2);
  b2[1] = dot(x0-x2, y0-y2, x1-x0, y1-y0);
  b2[2] = dot(x1-x0, y1-y0, x2-x1, y2-y1);
  
  samfloat_t g[3];
  for(size_t i=0;i<3;i++){
    b2[i]*=b2[i];
    g[i] = l2[i]*l2[(i+1)%3] - b2[i];
  }
  
  samfloat_t functional = 0;
  for(size_t i=0;i<3;i++){
    samfloat_t a = (l2[i]-l2[(i+1)%3]);
    functional += (a*a + 4*b2[i])/g[i];
    // cout<<"functional "<<functional<<" "<<a<<" "<<b2[i]<<" "<<g[i]<<endl;
  }
  
  return 0.5*functional;
}
