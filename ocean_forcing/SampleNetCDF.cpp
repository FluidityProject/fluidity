/*
 *    Copyright (c) 2004-2006 by Gerard Gorman
 *    Copyright (c) 2006- Imperial College London
 *    See COPYING file for copying and redistribution conditions.
 *   
 *    This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation,
 *    version 2.1 of the License.
 *
 *    This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free Software
 *    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
 *    USA
 *
 *    Contact info: gerard.j.gorman@gmail.com/g.gorman@imperial.ac.uk
 */

#include <cassert>
#include <iostream>
#include <map>
#include <vector>

#include <vtk.h>

#include "SampleNetCDF.h"

using namespace std;

SampleNetCDF::SampleNetCDF(){
  VerboseOff();
}

SampleNetCDF::~SampleNetCDF(){}


SampleNetCDF::SampleNetCDF(string filename){
  VerboseOff();
  SetFile(filename);
  return;
}

void SampleNetCDF::Close(){
  if(verbose)
    cout<<"void SampleNetCDF::Close()\n";
  
  reader.Close();
}

const SampleNetCDF& SampleNetCDF::operator=(const SampleNetCDF &in){
  if(verbose)
    cout<<"const SampleNetCDF& SampleNetCDF::operator=(const SampleNetCDF &in)\n";

  for(int i=0;i<2;i++){
    x_range[i] = in.x_range[i];
    y_range[i] = in.y_range[i];
  
    spacing[i] = in.spacing[i];
    dimension[i] = in.dimension[i];
  }

  data = in.data;
  verbose = in.verbose;
  
  reader = in.reader;
  
  return *this;
}

bool SampleNetCDF::HasPoint(double longitude, double latitude) const{
  if(verbose)
    cout<<"double SampleNetCDF::HasValue("<<longitude<<", "<<latitude<<") const\n";
  
  if(longitude<x_range[0])
    return false;
  if(longitude>x_range[1])
    return false;
  
  if(latitude<y_range[0])
    return false;
  if(latitude>y_range[1])
    return false;
  
  return true;
}

double SampleNetCDF::GetValue(int ix, int iy) const{
  if(verbose)
    cout<<"double SampleNetCDF::GetValue("<<ix<<", "<<iy<<") const";

  //assert(ix>=0);
  ix = std::max(ix, 0);

  //assert(ix<dimension[0]);
  ix = std::min(ix, dimension[0]-1);

  //assert(iy>=0);
  iy = std::max(iy, 0);

  //assert(iy<dimension[1]);
  iy = std::min(iy, dimension[1]-1);

  double val = data[iy*dimension[0]+ix];
  if(verbose)
    cout<<" = "<<val<<endl;
  
  return val;
}

double SampleNetCDF::GetValue(double longitude, double latitude) const{
  if(verbose)
    cout<<"double SampleNetCDF::GetValue("<<longitude<<", "<<latitude<<") const\n";

  longitude = std::min(std::max(x_range[0], longitude), x_range[1]);
  latitude = std::min(std::max(y_range[0], latitude), y_range[1]);
  if(verbose)
    cout<<"Shifted coordinate = "<<longitude<<" "<<latitude<<endl;

  int i0 = (int)floor((longitude - x_range[0])/spacing[0]);
  int i1 =  (int)ceil((longitude - x_range[0])/spacing[0]);
  
  int j0 = (int)floor((latitude - y_range[0])/spacing[1]);
  int j1 =  (int)ceil((latitude - y_range[0])/spacing[1]);
  
  double x0 = x_range[0] + i0*spacing[0];
  double x1 = x_range[0] + i1*spacing[0];

  double y0 = y_range[0] + j0*spacing[1];
  double y1 = y_range[0] + j1*spacing[1];

  if(verbose)
    cout<<"longitude "<<x0<<" < "<<longitude<<" < "<<x1<<endl
        <<"latitude  "<<y0<<" < "<<latitude <<" < "<<y1<<endl
        <<"spacing   "<<spacing[0]<<", "<<spacing[1]<<endl
        <<"x range   "<<x_range[0]<<", "<<x_range[1]<<endl
        <<"y range   "<<y_range[0]<<", "<<y_range[1]<<endl;
  
  double z00 = GetValue(i0, j0);
  double z01 = GetValue(i0, j1);
  double z10 = GetValue(i1, j0);
  double z11 = GetValue(i1, j1);
    
  if(i0==i1){ // No interpolation along longitude
    if(verbose)
      cout<<"Case 1\n";
    if(j0==j1){
      return z00;
    }else{
      return z01 - (y1-latitude)*(z01-z00)/(y1-y0);
    }
  }else if(j0==j1){ // No interpolation along latitude
    if(verbose)
      cout<<"Case 2\n";
    return z11 - (x1-longitude)*(z11-z10)/(x1-x0);
  }
  //else{ // Bi-linear interpolation
  if(verbose)
    cout<<"Case 3"<<endl;
  double dx = x1-x0;
  double dy = y1-y0;
  
  double val = (z00*(x1-longitude)*(y1-latitude) +
                z10*(longitude-x0)*(y1-latitude) +
                z01*(x1-longitude)*(latitude-y0) +
                z11*(longitude-x0)*(latitude-y0))/(dx*dy);

  return val; 
}

void SampleNetCDF::SetFile(string filename){
  if(verbose)
    cout<<"void SampleNetCDF::SetFile("<<filename<<")\n";

  reader.SetFile(filename.c_str());
  
  reader.GetXRange(x_range[0], x_range[1]);
  reader.GetYRange(y_range[0], y_range[1]);
  
  reader.GetSpacing(spacing[0], spacing[1]);
  reader.GetDimension(dimension[0], dimension[1]);
}

void SampleNetCDF::SetVariable(string varname){
  if(verbose)
    cout<<"void SampleNetCDF::SetVariable("<<varname<<")\n";

  reader.Read(varname, data);
}

void SampleNetCDF::VerboseOff(){
  verbose=false;
  reader.VerboseOff();
}

void SampleNetCDF::VerboseOn(){
  verbose=true;
  reader.VerboseOn();
}


// Fortran interface
map<int, SampleNetCDF> netcdf_sampler;

extern "C" {
#define samplenetcdf_open_fc F77_FUNC_(samplenetcdf_open_c, SAMPLENETCDF_OPEN_C)
  void samplenetcdf_open_fc(const char *name, const int *len, int *id){
    *id=0;
    for(;;(*id)++)
      if(netcdf_sampler.find(*id)==netcdf_sampler.end())
        break;
    
    netcdf_sampler[*id] = SampleNetCDF(string(name, *len));
  }
  
#define samplenetcdf_setvariable_fc F77_FUNC_(samplenetcdf_setvariable_c, SAMPLENETCDF_SETVARIABLE_C)
  void samplenetcdf_setvariable_fc(const int *id, const char *varname, const int *len){
    if(netcdf_sampler.find(*id)==netcdf_sampler.end()){
      cerr<<"ERROR: netcdf file has not been opened\n";
      exit(-1);
    }
    
    netcdf_sampler[*id].SetVariable(string(varname, *len));
  }
  
#define samplenetcdf_getvalue_fc F77_FUNC_(samplenetcdf_getvalue_c, SAMPLENETCDF_GETVALUE_C)
  void samplenetcdf_getvalue_fc(const int *id, const double *longitude, const double *latitude, double *val){
    if(netcdf_sampler.find(*id)==netcdf_sampler.end()){
      cerr<<"ERROR: netcdf file has not been opened\n";
      exit(-1);
    }
    
    *val = netcdf_sampler[*id].GetValue(*longitude, *latitude);
  }
  
#define samplenetcdf_close_fc F77_FUNC_(samplenetcdf_close_c, SAMPLENETCDF_CLOSE_C)
  void samplenetcdf_close_fc(const int *id){
    map<int, SampleNetCDF>::iterator it = netcdf_sampler.find(*id);
    if(it!=netcdf_sampler.end()){
      netcdf_sampler[*id].Close();
      netcdf_sampler.erase(it);
    }
  }
}
