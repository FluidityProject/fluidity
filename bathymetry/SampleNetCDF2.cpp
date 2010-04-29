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

#include <cassert>
#include <iostream>
#include <map>
#include <vector>

#ifdef USING_VTK
#include <vtkImageData.h>
#include <vtkShortArray.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPointData.h>
#include <vtkZLibDataCompressor.h>
#include <vtkBMPWriter.h>
#endif

#include "SampleNetCDF2.h"

using namespace std;

SampleNetCDF2::SampleNetCDF2(){
  VerboseOff();
  is_constant = true;
}

SampleNetCDF2::~SampleNetCDF2(){}

SampleNetCDF2::SampleNetCDF2(string filename){
  VerboseOff();
  SetFile(filename);
  return;
}

bool SampleNetCDF2::HasPoint(double longitude, double latitude) const{
  if(verbose)
    cout<<"double SampleNetCDF2::HasValue("<<longitude<<", "<<latitude<<") const\n";
  
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

double SampleNetCDF2::GetValue(int ix, int iy) const{
  if(verbose)
    cout<<"double SampleNetCDF2::GetValue("<<ix<<", "<<iy<<") const";

  assert(!is_constant);

  //assert(ix>=0);
  ix = std::max(ix, 0);

  //assert(ix<dimension[0]);
  ix = std::min(ix, dimension[0]-1);

  //assert(iy>=0);
  iy = std::max(iy, 0);

  //assert(iy<dimension[1]);
  iy = std::min(iy, dimension[1]-1);

  if(verbose)
    cout<<" = "<<data[iy*dimension[0]+ix]<<endl;
  
  return data[iy*dimension[0]+ix];
}

double SampleNetCDF2::GetValue(double longitude, double latitude) const{
  if(verbose)
    cout<<"double SampleNetCDF2::GetValue("<<longitude<<", "<<latitude<<") const\n";

  if(is_constant)
    return 1.0;

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

  double val;
  if(i0==i1){ // No interpolation along longitude
    if(verbose)
      cout<<"Case 1\n";
    if(j0==j1){
      val = z00;
    }else{
      val = z00 + (latitude-y0)*(z01-z00)/(y1-y0);
    }
    if(verbose)
      cout<<"z00, z, z01 = "<<z00<<", "<<val<<", "<<z01<<endl;
  }else if(j0==j1){ // No interpolation along latitude
    if(verbose)
      cout<<"Case 2\n";
    val =  z00 + (longitude-x0)*(z10-z00)/(x1-x0);
    
    if(verbose)
      cout<<"z00, z, z10 = "<<z00<<", "<<val<<", "<<z10<<endl;
  }else{ // Bi-linear interpolation
    if(verbose)
      cout<<"Case 3"<<endl;
    double dx = x1-x0;
    double dy = y1-y0;
    
    val = (z00*(x1-longitude)*(y1-latitude) +
           z10*(longitude-x0)*(y1-latitude) +
           z01*(x1-longitude)*(latitude-y0) +
           z11*(longitude-x0)*(latitude-y0))/(dx*dy);
    
    if(verbose)
      cout<<"z00, z10, z01, z11, val = "
          <<z00<<", "<<z10<<", "<<z01<<", "<<z11<<", "<<val<<endl;
  }

  return val; 
}

void SampleNetCDF2::SetFile(string filename){
  if(filename.empty()){
    is_constant=true;
    return;
  }
  
  NetCDF_reader reader(filename.c_str(), verbose);
  
  reader.GetXRange(x_range[0], x_range[1]);
  reader.GetYRange(y_range[0], y_range[1]);
  
  reader.GetSpacing(spacing[0], spacing[1]);
  reader.GetDimension(dimension[0], dimension[1]);
  
  reader.Read(data);
  
  is_constant = false;
}

void SampleNetCDF2::VerboseOff(){
  verbose=false;
}

void SampleNetCDF2::VerboseOn(){
  verbose=true;
}
