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
#ifndef SAMPLENETCDF2_H
#define SAMPLENETCDF2_H
#include <cmath>
#include <string>
#include <vector>

#include <string.h>
extern "C" {
#include <netcdf.h>
}
// #include "tools.h" // Experimenting to see whether this is needed
#include "NetCDF_reader.h"

class SampleNetCDF2{
 public:
  SampleNetCDF2();
  SampleNetCDF2(std::string);
  ~SampleNetCDF2();
  
  double GetValue(double, double) const;

  bool HasPoint(double, double) const;

  void SetFile(std::string);
  
  void VerboseOff();
  void VerboseOn();

 private:
  void Debug(char *) const;
  double GetValue(int, int) const;

  bool is_constant;
  
  double x_range[2];          // Longitude
  double y_range[2];          // Latitude
  
  double spacing[2];          // Grid spacing
  int dimension[2];           // Grid dimensions
  
  std::vector<double> data;
  bool verbose;
};
#endif
