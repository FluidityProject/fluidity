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

extern "C" {
#define get_field_values_fc F77_FUNC(get_field_values, GET_FIELD_VALUES)
        void get_field_values_fc(const char* filename, const double *X, const double *Y, double *Z, int *nodes);

}

void get_field_values_fc(const char* filename, const double *X, const double *Y, double *Z, int *n){

        string file=string(filename);
        SampleNetCDF2 map(file);

        const int nodes = *n;
        double *x = new double[nodes];
        double *y = new double[nodes];
        for (int i = 0; i < nodes; i++) {
          x[i] = X[i];
          y[i] = Y[i];
          if(map.HasPoint(x[i], y[i])){
                Z[i]=map.GetValue(x[i], y[i]);
          }else{
                cerr << "Point [" << x[i] << ", " << y[i] << "] not found in netCDF file, " << file << ".\n";
                Z[i]=0.0;
          }
        }

        delete [] x;
        delete [] y;

}
