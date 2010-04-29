/*  Copyright (C) 2009 Imperial College London and others.
    
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


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include "SampleNetCDF2.h"
#include "confdefs.h"

using namespace std;

void get_lldepth(double x, double y, double z, double &longitude, double &latitude, double &depth){

  const double pi=4.0*atanl(1.0);
  const double earth_radius=6378100.0;
  const double rad_to_deg=180.0/pi;
  
  double r = sqrt(x*x+y*y+z*z);
  longitude = atan2(y, x);
  latitude = acos(z/r);
  depth = earth_radius - r;
  longitude*=rad_to_deg;
  latitude = 90 - latitude*rad_to_deg;
}

extern "C" {
#define set_from_map_fc F77_FUNC(set_from_map, SET_FROM_MAP)
        void set_from_map_fc(const char* filename, double *x, double *y, double *z, double *depth);
}

void set_from_map_fc(const char* filename, double *x, double *y, double *z, double *depth){

        string file=string(filename);
        SampleNetCDF2 map(filename);
        double longitude, latitude, depth_tmp;
        get_lldepth(*x, *y, *z, longitude, latitude, depth_tmp);
        if(map.HasPoint(longitude, latitude)){
                *depth=map.GetValue(longitude, latitude);
        }else{
                cerr<<"Point not found in netcdf file\n";
                *depth=0.0;
        }
}

#ifdef DEPTH_UNITTEST
// This is a simple test program to check the depth returned for three given lon,lat locations
int main(int argc, char **argv){

  double pi=4.0*atanl(1.0);
  double deg_to_rad=pi/180.0;
  double depth=0.0;

  bool open_file=true;

  string data_file="/data/maps/gridone.grd";

  double r=6378100.0;
  double longitude=5.974731;
  double latitude=38.387311;

  double xloc=r*cos(latitude*deg_to_rad)*cos(longitude*deg_to_rad);
  double yloc=r*cos(latitude*deg_to_rad)*sin(longitude*deg_to_rad);
  double zloc=r*sin(latitude*deg_to_rad);

  set_from_map_fc(data_file.c_str(), &xloc, &yloc, &zloc, &depth);

  cout << "The depth at ( " << longitude << " , " << latitude << " ) is " << depth << endl;

  open_file=false;

  longitude=6.974731;
  latitude=39.387311;

  xloc=r*cos(latitude*deg_to_rad)*cos(longitude*deg_to_rad);
  yloc=r*cos(latitude*deg_to_rad)*sin(longitude*deg_to_rad);
  zloc=r*sin(latitude*deg_to_rad);

  set_from_map_fc(data_file.c_str(), &xloc, &yloc, &zloc, &depth);

  cout << "The depth at ( " << longitude << " , " << latitude << " ) is " << depth << endl;

  longitude=5.974731;
  latitude=40.387311;

  xloc=r*cos(latitude*deg_to_rad)*cos(longitude*deg_to_rad);
  yloc=r*cos(latitude*deg_to_rad)*sin(longitude*deg_to_rad);
  zloc=r*sin(latitude*deg_to_rad);

  set_from_map_fc(data_file.c_str(), &xloc, &yloc, &zloc, &depth);

  cout << "The depth at ( " << longitude << " , " << latitude << " ) is " << depth << endl;
  cout << "End of unit test.\n";  
  
}
#endif

