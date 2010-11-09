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
#ifndef CLIMATEREADER_H
#define CLIMATEREADER_H
#include "confdefs.h"

#include "Calendar.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_LIBNETCDF
#include <netcdf.h>
#endif

class ClimateReader{
 public:
  ClimateReader();
  ~ClimateReader();

  double GetValue(std::string name, double x, double y, double z);
  double GetValue(std::string name, double x, double y, double z, int ilevel);
  int SetClimatology(std::string filename);
  int SetSimulationTimeUnits(std::string str);
  int SetTimeSeconds(double sim_seconds);
  void VerboseOff();
  void VerboseOn();
  
 private:
  int Cartesian2Grid(double x, double y, double z, int &ilong, int &ilat, int &ilevel);
  int Cartesian2Spherical(double x, double y, double z, double &longitude, double &latitude, double &depth);
  double GetValue(std::string name, int ilong, int ilat, int ilevel);
  double SolveLine(double x0, double x1, double t0, double t1, double t) const;
  int Spherical2Grid(double longitude, double latitude, double depth, int &ilong, int &ilat, int &ilevel);
  double Uncompress(short svar);

  double earth_radius, pi, rad_to_deg, deg_to_rad;
  int idim0, jdim0;
  double space;
  std::vector<float> levels, months, seasons;
  int ncid;
  
  double scale_factor, add_offset;
  short fill_value, missing_value;
  bool verbose, have_ref_date;

  std::string simulation_time_units, data_time_units;
  Calendar *calendar;
  double time_set; // seriously needs revisiting
  std::pair<int, double> time_month, time_season;
};

extern ClimateReader ClimateReader_global;
#endif
