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
#ifndef NETCDF_READER_H
#define NETCDF_READER_H

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifdef USING_VTK
#include <vtkImageData.h>
#include <vtkShortArray.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPointData.h>
#include <vtkZLibDataCompressor.h>
#include <vtkBMPWriter.h>
#endif

extern "C" {
#include <netcdf.h>
}

class NetCDF_reader{
 public:
  NetCDF_reader(const char *, bool);
  ~NetCDF_reader();
  
  int GetDimension(int &ilen, int &jlen) const;
  int GetSpacing(double &dx, double &dy) const;
  int GetXRange(double &xmin, double &xmax) const;
  int GetYRange(double &ymin, double &ymax) const;
  bool IsColatitude() const;
  int Read(std::vector<double> &z) const;

  void NetCDFErrorCheckingOff();
  void NetCDFErrorCheckingOn();

  void VerboseOff();
  void VerboseOn();
  
 private:
  int ncid;                   // netCDF file handle
  enum {gebco , coards, coards_cf10, unknown} data_source;
                              // 1 - "GEBCO One Minute Grid"
  double x_range[2];          // Longitude
  double y_range[2];          // Latitude

  double spacing[2];          // Grid spacing
  int dimension[2];           // Grid dimensions

  long xysize;                // xysize

  bool verbose;
};
#endif
