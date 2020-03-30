/*  Copyright (C) 2006 Imperial College London and others.
    
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

#ifndef FLUXESREADER_H
#define FLUXESREADER_H

#include "fmangle.h"
#include "c++debug.h"

#include <cassert>
#include <cmath>
#ifdef __SUN
double round(double);
#endif
#include <deque>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifdef HAVE_LIBNETCDF
extern "C" {
#include <netcdf.h>
}
#endif

#include "flmpi.h"

#include <stdio.h>
#include <stdlib.h>

#include "Calendar.h"

class FluxesReader{
 public:
  FluxesReader();
  ~FluxesReader();
  
  void AddFieldOfInterest(std::string);
  void ClearFields();
  bool Enabled() const;
  int GetScalars(double, double, double*);
  int RegisterDataFile(std::string);
  int SetSimulationTimeUnits(std::string);
  int SetTimeSeconds(double);
  void VerboseOff();
  void VerboseOn();
  double GetScalar(std::string, double, double);

 private:
  bool verbose, modified;
  double time_set;

  std::pair<size_t, size_t> GetInterval(double value, const std::map<double, int>&lut) const;

  double GetScaleFactor(std::string) const;
  int Read(int);
  int Update();
  inline double SolveLine(double y, double x0, double x1, double y0, double y1) const;

  int NProcs, MyRank, NParts;
  std::deque<std::string> ERA_data_files;
  
  //
  // data file specification
  //
  int ncid; // fluxes file handle
  std::map<std::string, long> spec;
  std::vector<double> longitude, latitude;
  std::map<double, int> lut_longitude, lut_latitude;
  double dlong, dlat;
  std::vector<double> time;
  std::map<double, int> lut_time;
  
  std::string simulation_time_units, data_time_units;
  Calendar *calendar;

  //
  // Data
  //
  int time_fore, time_aft;
  int hour_index;
  
  // Fields of interest
  std::deque<std::string> fields_of_interest;

  //    time index    |  name of field  | field values
  std::map<int, std::map<std::string, std::vector<double> > > fields;
};

extern FluxesReader FluxesReader_global;

#endif
