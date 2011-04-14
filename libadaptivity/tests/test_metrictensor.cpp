/* Copyright (C) 2006 Imperial College London and others.

 Please see the AUTHORS file in the main source directory for a full list
 of copyright holders.

 Dr Gerard J Gorman
 Applied Modelling and Computation Group
 Department of Earth Science and Engineering
 Imperial College London

 g.gorman@imperial.ac.uk

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
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

#include "MetricTensor.h"

using namespace std;

int main(){
  {
    MetricTensor Ma(82800,
                    -30857.1,    28800);
    Ma.write_vtk("M");
    cout<<"Ma =\n"<<Ma<<endl;
    
    double D[2], V[4];
    Ma.eigen_decomp(D, V);
    Ma.eigen_undecomp(D, V);
    cout<<"Ma (after eigen decomposition and reforming matrix) =\n"<<Ma<<endl;
    
    MetricTensor Mb(Ma);
    
    MetricTensor limit(10000,
                        0, 30000);
    cout<<"limit =\n"<<limit<<endl;
    limit.write_vtk("limit");
    
    Ma.circumscribe(limit);
    cout<<"circumscribed =\n"<<Ma<<endl;
    Ma.write_vtk("circumscribe");
    
    Mb.inscribe(limit);
    cout<<"inscribe =\n"<<Mb<<endl;
    Mb.write_vtk("inscribe");
  }
  
  {
    MetricTensor Ma(1.00467,
                    0.00162196, 6.2432,
                    -0.24383, -1.98199, 14.4226);
    Ma.write_vtk("M");
    cout<<"Ma =\n"<<Ma<<endl;
    
    double D[3], V[9];
    Ma.eigen_decomp(D, V);
    Ma.eigen_undecomp(D, V);
    cout<<"Ma (after eigen decomposition and reforming matrix) =\n"<<Ma<<endl;
    
    MetricTensor Mb(Ma);
    
    MetricTensor limit(4.05543,
                        0.0983953, 1.34818,
                        0.643409, -0.989943, 4.09606);
    
    cout<<"limit =\n"<<limit<<endl;
    limit.write_vtk("limit");
    
    Ma.circumscribe(limit);
    cout<<"circumscribed =\n"<<Ma<<endl;
    Ma.write_vtk("circumscribe");
    
    Mb.inscribe(limit);
    cout<<"inscribe =\n"<<Mb<<endl;
    Mb.write_vtk("inscribe");
  }

  {
    MetricTensor Ma(1.0,
                     0.0, 1.0,
                     0.0, 0.0, 1.0);
    Ma.write_vtk("M");
    cout<<"Ma =\n"<<Ma<<endl;
    
    double D[3], V[9];
    Ma.eigen_decomp(D, V);
    Ma.eigen_undecomp(D, V);
    cout<<"Ma (after eigen decomposition and reforming matrix) =\n"<<Ma<<endl;
    
    MetricTensor Mb(Ma);
    
    MetricTensor limit(0.5,
                        0.0, 0.5,
                        0.0, 0.0, 0.5);
    
    cout<<"limit =\n"<<limit<<endl;
    limit.write_vtk("limit");
    
    Ma.circumscribe(limit);
    cout<<"circumscribed =\n"<<Ma<<endl;
    Ma.write_vtk("circumscribe");
    
    Mb.inscribe(limit);
    cout<<"inscribe =\n"<<Mb<<endl;
    Mb.write_vtk("inscribe");
  }

  return(0);
}
