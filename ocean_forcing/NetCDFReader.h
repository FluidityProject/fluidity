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

#ifndef NETCDF_READER_H
#define NETCDF_READER_H

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "confdefs.h"

#ifdef USING_VTK
#include <vtk.h>
#endif

#ifdef HAVE_NETCDF
extern "C" {
#include <netcdf.h>
}
#endif

class NetCDFReader{
 public:
  NetCDFReader();
  NetCDFReader(const char *);
  ~NetCDFReader();
  
  // Overloaded operators.
  const NetCDFReader& operator=(const NetCDFReader &in);

  void Close();

  int GetDimension(int &ilen, int &jlen) const;  
  int GetSpacing(double &dx, double &dy) const;  
  int GetXRange(double &xmin, double &xmax) const;  
  int GetYRange(double &ymin, double &ymax) const; 
  int Read(std::string varname, std::vector<double> &var);

  void SetFile(const char *filename);

  void VerboseOff();
  void VerboseOn();

 private:
  void GetGrid();

  int ncid;                  // netCDF file handle
  bool fileOpen;             // true if ncid refers to an open netCDF file

  double x_range[2];          // Longitude
  double y_range[2];          // Latitude

  double spacing[2];         // Grid spacing
  long dimension[2];         // Grid dimensions

  bool verbose;
};
#endif
