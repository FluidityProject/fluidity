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

void get_ll(double x, double y, double z, double &longitude, double &latitude){

  const double pi=4.0*atanl(1.0);
  const double rad_to_deg=180.0/pi;
  
  double r = sqrt(x*x+y*y+z*z);
  longitude = atan2(y, x);
  latitude = acos(z/r);
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
          double longitude, latitude;
          get_ll(x[i], y[i], z[i], longitude, latitude);
          if(map.HasPoint(longitude, latitude)){
                height[i]=map.GetValue(longitude, latitude);
          }else{
                cerr << "Point [" << x[i] << ", " << y[i] << ", " << z[i] << "] with longitude and latitude [" << longitude << ", " << latitude <<  "] not found in netCDF file, " << file << ".\n";
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
                cerr << "Point [" << x[i] << ", " << y[i] << "] not found in netCDF file, " << file << ".\n";
                height[i]=0.0;
          }
        }

        delete [] x;
        delete [] y;

        for (int i=0; i<ncolumns; i++) {

          depth[i]=sh-height[i];

        }
}

