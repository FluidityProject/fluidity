/*  Copyright (C) 2009 Imperial College London and others.
    
Please see the AUTHORS file in the main source directory for a full list
of copyright holders.
        
Prof. C Pain
Applied Modelling and Computation Group
Department of Earth Science and Engineering
Imperial College London
                                    
amcgsoftware@imperial.ac.uk
                                            
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
  latitude = 90.0 - latitude*rad_to_deg;
}

extern "C" {
#define set_from_map_fc F77_FUNC(set_from_map, SET_FROM_MAP)
        void set_from_map_fc(const char* filename, const double *X, const double *Y, const double *Z, double *depth, int *ncolumns);

#define set_from_map_beta_fc F77_FUNC(set_from_map_beta, SET_FROM_MAP_BETA)
        void set_from_map_beta_fc(const char* filename, const double *X, const double *Y, double *depth, int *ncolumns, double *surf_h);

}

void set_from_map_fc(const char* filename, const double *X, const double *Y, const double *Z, double *depth, int *n){

        string file=string(filename);
        SampleNetCDF2 map(file);

        const int ncolumns = *n;
        double *x = new double[ncolumns];
        double *y = new double[ncolumns];
        double *z = new double[ncolumns];
        double height[ncolumns];
        for (int i = 0; i < ncolumns; i++) {
          x[i] = X[i];
          y[i] = Y[i];
          z[i] = Z[i];
          double longitude, latitude, depth_tmp;
          get_lldepth(x[i], y[i], z[i], longitude, latitude, depth_tmp);
          if(map.HasPoint(longitude, latitude)){
                height[i]=map.GetValue(longitude, latitude);
          }else{
                cerr<<"Point not found in netcdf file\n";
                height[i]=0.0;
          }
        }

        delete [] x;
        delete [] y;
        delete [] z;

        for (int i=0; i<ncolumns; i++) {

          depth[i]=-height[i];

        }
}

void set_from_map_beta_fc(const char* filename, const double *X, const double *Y, double *depth, int *n, double *surf_h){

        string file=string(filename);
        SampleNetCDF2 map(file);

        const int ncolumns = *n;
        const double sh = *surf_h;
        double *x = new double[ncolumns];
        double *y = new double[ncolumns];
        double height[ncolumns];
        for (int i = 0; i < ncolumns; i++) {
          x[i] = X[i];
          y[i] = Y[i];
          if(map.HasPoint(x[i], y[i])){
                height[i]=map.GetValue(x[i], y[i]);
          }else{
                cerr<<"Point not found in netcdf file\n";
                height[i]=0.0;
          }
        }

        delete [] x;
        delete [] y;

        for (int i=0; i<ncolumns; i++) {

          depth[i]=sh-height[i];

        }
}

#ifdef DEPTH_UNITTEST
// This is a simple test program to check the depth returned for three given lon,lat locations
int main(int argc, char **argv){

  double pi=4.0*atanl(1.0);
  double deg_to_rad=pi/180.0;
  double depth[3];
  int columns=3;

  string data_file="/data/maps/gridone.grd";

  double xloc[3], yloc[3], zloc[3], longitude[3], latitude[3];

  double r=6378100.0; // Approx. radius of Earth

  // 1st location
  longitude[0]=5.974731;
  latitude[0]=38.387311;

  xloc[0]=r*cos(latitude[0]*deg_to_rad)*cos(longitude[0]*deg_to_rad);
  yloc[0]=r*cos(latitude[0]*deg_to_rad)*sin(longitude[0]*deg_to_rad);
  zloc[0]=r*sin(latitude[0]*deg_to_rad);

  // 2nd location
  longitude[1]=6.974731;
  latitude[1]=39.387311;

  xloc[1]=r*cos(latitude[1]*deg_to_rad)*cos(longitude[1]*deg_to_rad);
  yloc[1]=r*cos(latitude[1]*deg_to_rad)*sin(longitude[1]*deg_to_rad);
  zloc[1]=r*sin(latitude[1]*deg_to_rad);

  // 3rd location
  longitude[2]=5.974731;
  latitude[2]=40.387311;

  xloc[2]=r*cos(latitude[2]*deg_to_rad)*cos(longitude[2]*deg_to_rad);
  yloc[2]=r*cos(latitude[2]*deg_to_rad)*sin(longitude[2]*deg_to_rad);
  zloc[2]=r*sin(latitude[2]*deg_to_rad);

  cout << "Retrieving bathymetry data.\n";
  set_from_map_fc(data_file.c_str(), xloc, yloc, zloc, depth, &columns);

  for (int i = 0; i < columns; i++) {
    cout << "The depth at ( " << longitude[i] << " , " << latitude[i] << " ) is " << depth[i] << endl;
  }
  cout << "End of unit test.\n";  
  
}
#endif

#ifdef DEPTH_UNITTEST2
// This is a simple test program to check the depth returned for three given x,y beta plane locations
int main(int argc, char **argv){

  double depth[3];
  int columns=3;

  string data_file="/data/myfluidity/tests/MonaiValley/outfileMonaiValley.grd";

  double xloc[3], yloc[3];
  
  // The Monaia Valley data bounds are
  // x:actual_range = 0., 5.4916666667
  // y:actual_range = 0., 3.4083333333

  // 1st location
  xloc[0]=1.0;
  yloc[0]=1.0;

  // 2nd location
  xloc[1]=2.0;
  yloc[1]=2.0;

  // 3rd location
  xloc[2]=3.0;
  yloc[2]=3.0;

  cout << "Retrieving bathymetry data.\n";
  set_from_map_beta_fc(data_file.c_str(), xloc, yloc, depth, &columns);

  for (int i = 0; i < columns; i++) {
    cout << "The depth at ( " << xloc[i] << " , " << yloc[i] << " ) is " << depth[i] << endl;
  }
  cout << "End of unit test.\n";  
  
}
#endif

