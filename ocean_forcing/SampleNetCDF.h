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

#ifndef SAMPLENETCDF_H
#define SAMPLENETCDF_H
#include <cmath>
#include <string>
#include <vector>

#include <string.h>
#ifdef HAVE_NETCDF
extern "C" {
#include <netcdf.h>
}
#endif
#include "NetCDFReader.h"
#include "confdefs.h"

class SampleNetCDF{
 public:
  SampleNetCDF();
  SampleNetCDF(std::string);
  ~SampleNetCDF();

  // Overloaded operators.
  const SampleNetCDF& operator=(const SampleNetCDF &in);
  
  void Close();
  double GetValue(double, double) const;

  bool HasPoint(double, double) const;

  void SetFile(std::string);
  void SetVariable(std::string);
  
  void VerboseOff();
  void VerboseOn();
  
 private:
  void Debug(char *) const;
  double GetValue(int, int) const;

  double x_range[2];          // Longitude
  double y_range[2];          // Latitude
  
  double spacing[2];          // Grid spacing
  int dimension[2];           // Grid dimensions
  
  std::vector<double> data;
  bool verbose;

  NetCDFReader reader;
};
#endif
