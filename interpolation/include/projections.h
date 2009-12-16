/*
  
  This code is taken from Terreno with the original copyright notice:

  *
  *      Copyright (c) 2004-2006 by Gerard Gorman
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
  *      Contact info: gerard.j.gorman@gmail.com
  *
  
  It was re-licensed by Gerard Gorman on 26th December 2007 to
  Fluidity with the new license:
  
  Copyright (C) 2006 Imperial College London and others.
  
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

#ifndef PROJECTIONS_H
#define PROJECTIONS_H

#include <cmath>

#ifndef PI
#define PI 3.1415926535897931
#endif

#ifndef RAD_TO_DEG
#define RAD_TO_DEG 57.295779523
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD (1.0/RAD_TO_DEG)
#endif

#ifndef RADIUS_OF_EARTH
#define RADIUS_OF_EARTH 6378000.0
#endif

namespace projection{
  template<class R>
    R degrees(R radians){
    return radians*RAD_TO_DEG;
  }
  
  template<class R>
    R radians(R degrees){
    return degrees*DEG_TO_RAD;
  }

  template<class R>
    int cartesian2spherical(R x, R y, R z, R &latitude, R &longitude, R &depth){
    double r = sqrt(x*x+y*y+z*z);
    
    depth = r - RADIUS_OF_EARTH;
    // latitude  = 90.0 - theta;
    latitude  = 90.0 - projection::degrees(acos(z/r));
    longitude = projection::degrees(atan2(y, x));
    
    return 0;
  }
  
  template<class R>
    int spherical2cartesian(R latitude, R longitude, R depth, R &x, R &y, R &z){
    // latitude  = 90.0 - theta;
    latitude  = projection::radians(latitude);
    longitude = projection::radians(longitude);
    
    x = (RADIUS_OF_EARTH+depth)*cos(longitude)*sin(PI*0.5-latitude);
    y = (RADIUS_OF_EARTH+depth)*sin(longitude)*sin(PI*0.5-latitude);
    z = (RADIUS_OF_EARTH+depth)*cos(PI*0.5-latitude);
    
    return 0;
  }

  template<class R>
    int spherical2stereographic(R longitude_0, R latitude_0,
                                R longitude, R latitude,
                                R &x, R &y){
    // http://mathworld.wolfram.com/StereographicProjection.html
    longitude = radians(longitude);
    latitude  = radians(latitude);
        
    if(latitude_0==latitude){
      x = 0.0;
      y = 0.0;
      return -1;
    }
    
    double k = 1.0/
    (1.0 + 
     sin(latitude_0)*sin(latitude) + 
     cos(latitude_0)*cos(latitude)*cos(longitude-longitude_0));
    
    x = k*cos(latitude)*sin(longitude-longitude_0);
    y = k*(cos(latitude_0)*sin(latitude) - 
	   sin(latitude_0)*cos(latitude)*cos(longitude-longitude_0));
    
    return 0;
  }
  
  template<class R>
    int cartesian2stereographic(R longitude_0, R latitude_0,
                                R x, R y, R z,
                                R &xx, R &yy){
    R longitude, latitude, depth;
    
    cartesian2spherical(x, y, z, latitude, longitude, depth);
    spherical2stereographic(longitude_0, latitude_0,
                            longitude, latitude,
                            xx, yy);

    return 0;
  }

  template<class R>
    int cartesian2stereographic(R longitude_0, R latitude_0,
                                R x, R y, R z,
                                R &xx, R &yy, R &depth){
    R longitude, latitude;
    
    cartesian2spherical(x, y, z, latitude, longitude, depth);
    spherical2stereographic(longitude_0, latitude_0,
                            longitude, latitude,
                            xx, yy);
    
    return 0;
  }

  template<class R>
    int stereographic2cartesian(R longitude_0, R latitude_0,
                                R xx, R yy,
                                R &x, R &y, R &z){
    R longitude, latitude;
    stereographic2spherical(longitude_0, latitude_0,
                            xx, yy,
                            longitude, latitude);
    spherical2cartesian(longitude, latitude, x, y, z);
    return 0;
  }
  
  template<class R>
    int stereographic2spherical(R longitude_0, R latitude_0,
                                R x, R y,
                                R &longitude, R &latitude){
    // http://mathworld.wolfram.com/StereographicProjection.html
    if((x==0.0)&&(y==0.0)){
      latitude  = degrees(latitude_0);
      longitude = degrees(longitude_0);
      return 0;
    }
    
    double rho = sqrt(x*x + y*y);
    double c = 2.0*atan2(rho, 1.0);
    
    latitude = asin(cos(c)*sin(latitude_0) + y*sin(c)*cos(latitude_0)/rho);
    longitude = longitude_0 + atan2(x*sin(c), 
				    rho*cos(latitude_0)*cos(c) -
				    y*sin(latitude_0)*sin(c));
    
    latitude  = degrees(latitude);
    longitude = degrees(longitude);
  
    return 0;
  }
}

#endif
