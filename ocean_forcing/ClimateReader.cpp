/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    C.Pain@Imperial.ac.uk
    
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

#include "ClimateReader.h"
using namespace std;

ClimateReader::ClimateReader(){
  earth_radius=6378000.0;
  pi = 4.0*atan(1.0);
  rad_to_deg = 180.0/pi;
  deg_to_rad = pi/180;
  
  idim0=1440;
  jdim0=720;
  space = 0.25;
  
  ncid = -1;
  
  VerboseOff();
  have_ref_date = false;

  time_month = pair<int, double>(0, 0.0);
  time_season = pair<int, double>(0, 0.0);
  
  // Standard levels
  levels.push_back(0.);
  levels.push_back(-10.);
  levels.push_back(-20.);
  levels.push_back(-30.);
  levels.push_back(-50.); 
  levels.push_back(-75.);
  levels.push_back(-100.); 
  levels.push_back(-125.); 
  levels.push_back(-150.); 
  levels.push_back(-200.); 
  levels.push_back(-250.); 
  levels.push_back(-300.); 
  levels.push_back(-400.); 
  levels.push_back(-500.);
  levels.push_back(-600.); 
  levels.push_back(-700.); 
  levels.push_back(-800.); 
  levels.push_back(-900.); 
  levels.push_back(-1000.); 
  levels.push_back(-1100.); 
  levels.push_back(-1200.); 
  levels.push_back(-1300.); 
  levels.push_back(-1400.);
  levels.push_back(-1500.);
  levels.push_back(-1750.); 
  levels.push_back(-2000.);
  levels.push_back(-2500.);
  levels.push_back(-3000.);
  levels.push_back(-3500.);
  levels.push_back(-4000.);
  levels.push_back(-4500.);
  levels.push_back(-5000.);
  levels.push_back(-5500.);
  
  // Months in seconds
  months.push_back(24*60*60*31.0);
  months.push_back(24*60*60*28.25);
  months.push_back(24*60*60*31.0);
  months.push_back(24*60*60*30.0);
  months.push_back(24*60*60*31.0);
  months.push_back(24*60*60*30.0);
  months.push_back(24*60*60*31.0);
  months.push_back(24*60*60*31.0);
  months.push_back(24*60*60*30.0);
  months.push_back(24*60*60*31.0);
  months.push_back(24*60*60*30.0);
  months.push_back(24*60*60*31.0);
  
  // Seasons in seconds
  seasons.push_back(months[0]+months[1]+months[2]);
  seasons.push_back(months[3]+months[4]+months[5]);
  seasons.push_back(months[6]+months[7]+months[8]);
  seasons.push_back(months[9]+months[10]+months[11]);

  calendar = NULL;
}

ClimateReader::~ClimateReader(){
  if(calendar!=NULL)
    delete calendar;
}

int ClimateReader::Cartesian2Grid(double x, double y, double z, int &ilong, int &ilat, int &ilevel){
  if(verbose)
    cout<<"int Cartesian2Grid("<<x<<", "<<y<<", "<<z<<", "<<"int &ilong, int &ilat, int &ilevel)\n";
  
  double longitude, latitude, depth;
  Cartesian2Spherical(x, y, z, longitude, latitude, depth);
  Spherical2Grid(longitude, latitude, depth, ilong, ilat, ilevel);
  
  return 0;
}
  
int ClimateReader::Cartesian2Spherical(double x, double y, double z, double &longitude, double &latitude, double &depth){
  if(verbose)
    cout<<"int Cartesian2Spherical("<<x<<", "<<y<<", "<<z<<", double &longitude, double &latitude, double &depth)\n";
  
  double r = sqrt(x*x+y*y+z*z);
  
  depth = r - earth_radius;
  
  longitude = rad_to_deg*atan2(y, x);
  
  latitude = 90.0 - rad_to_deg*acos(z/r);
  
  return 0;
}

double ClimateReader::Uncompress(short svar){
  if((svar==fill_value)||(svar==missing_value)){
    if(verbose)
      cout<<"WARNING: missing value\n";
    return 0.0;
  }
  return svar*scale_factor + add_offset;
}

double ClimateReader::GetValue(string name, double x, double y, double z){
  if(verbose)
    cout<<"double GetValue("<<name<<", "<<x<<", "<<y<<", "<<z<<")\n";
  
  assert(ncid>0);
  
  int ilong, ilat, ilevel;
  Cartesian2Grid(x, y, z, ilong, ilat, ilevel);
  
  return GetValue(name, ilong, ilat, ilevel);
}

double ClimateReader::GetValue(string name, double x, double y, double z, int _ilevel){
  if(verbose)
    cout<<"double GetValue("<<name<<", "<<x<<", "<<y<<", "<<z<<")\n";
  
  assert(ncid>0);
  
  int ilong, ilat, ilevel;
  Cartesian2Grid(x, y, z, ilong, ilat, ilevel);
  
  // A specific level is required so clobber calculated level
  ilevel = _ilevel;
  
  return GetValue(name, ilong, ilat, ilevel);
}
  
double ClimateReader::GetValue(string name, int ilong, int ilat, int ilevel){
#ifdef HAVE_NETCDF
  if(verbose)
    cout<<"double GetValue("<<name<<", "<<ilong<<", "<<ilat<<", "<<ilevel<<")\n";
  
  assert(ncid>0);
  
  string varname(name);
  int t0, t1;
  double s, s0, s1;
  if(ilevel<24){
    varname+=string("_monthly");
    t0 = time_month.first;
    s = time_month.second;
    s0 = months[t0]*0.5; // mid-month
    if(s<s0){ // the middle of the previous month
      t1 = (t0+11)%12;
      s1 = -months[t1]*0.5;
    }else{ // the middle of the next month
      t1 = (t0+1)%12;
      s1 = months[t0]+months[t1]*0.5; 
    }
  }else{
    ilevel-=24;
    varname+=string("_seasonal");
    t0 = time_season.first;
    s = time_season.second;
    s0 = seasons[t0]*0.5; // mid-season
    if(s<s0){ // the middle of the previous season
      t1 = (t0+3)%4;
      s1 = -seasons[t1]*0.5; 
    }else{ // the middle of the next season
      t1 = (t0+1)%4;
      s1 = seasons[t0]+seasons[t1]*0.5; 
    }
  }
  
  int varid = ncvarid(ncid, varname.c_str());
  
  // Attributes to read
  nc_get_att_double(ncid, varid, "scale_factor",  &scale_factor);
  nc_get_att_double(ncid, varid, "add_offset",    &add_offset);
  nc_get_att_short(ncid,  varid, "_FillValue",    &fill_value);
  nc_get_att_short(ncid,  varid, "missing_value", &missing_value);
  
  long mindex[4];
  mindex[0] = t0;
  mindex[1] = ilevel;
  mindex[2] = ilat;
  mindex[3] = ilong;
  if(verbose)
    cout<<"Reading "<<varname<<" 0: "<<mindex[0]<<", "<<mindex[1]<<", "<<mindex[2]<<", "<<mindex[3]<<": "<<scale_factor<<", "<<add_offset<<endl;
  
  short value;
  int err = ncvarget1(ncid, varid, mindex, &value);
  assert(err<=0);
  double rval0 = Uncompress(value);
  
  mindex[0] = t1;
  if(verbose)
    cout<<"Reading "<<varname<<" 1: "<<mindex[0]<<", "<<mindex[1]<<", "<<mindex[2]<<", "<<mindex[3]<<": "<<scale_factor<<", "<<add_offset<<endl;
  
  err = ncvarget1(ncid, varid, mindex, &value);
  assert(err<=0);
  double rval1 = Uncompress(value);
  
  return SolveLine(rval0, rval1, s0, s1, s);
#else
  cerr<<"ERROR: no NetCDF support compiled\n";
  exit(-1);

  // This is just to keep the compiler quiet.
  return 0.0;
#endif
}

int ClimateReader::SetClimatology(string filename){
#ifdef HAVE_NETCDF
  if(verbose)
    cout<<"int set_climatology("<<filename<<")\n";

  ncid = ncopen(filename.c_str(), NC_NOWRITE);    
  return ncid;
#else
  cerr<<"ERROR: no NetCDF support compiled\n";
  exit(-1);

  // This is just to keep the compiler quiet.
  return 0;
#endif
}

int ClimateReader::SetSimulationTimeUnits(std::string str){
  simulation_time_units = str;
  have_ref_date = true;
  return 0;
}

int ClimateReader::SetTimeSeconds(double seconds){
  if(verbose)
    cout<<"int SetTimeSeconds("<<seconds<<")\n";

  // Convert from simulation time units to data time units
  calendar->SetTransformation(simulation_time_units, data_time_units, "gregorian");
  
  int err = calendar->Convert(seconds, time_set);
  if(err){
    cerr<<"ERROR ("<<__FILE__<<"): Conversion between time units has failed.\n";
    exit(-1);
  }
    
  return (0);
}

double ClimateReader::SolveLine(double x0, double x1, double t0, double t1, double t) const{
  if(verbose)
    cout<<"double ClimateReader::solve_line("<<x0<<", "<<x1<<", "<<t0<<", "<<t1<<", "<<t<<") const\n";
  
  if(fabs(t1-t0)<1.0) // This is an instant as far as climatology is concerned
    return 0.5*(x0+x1);
  return x0 + (x1-x0)*(t-t0)/(t1-t0);
}
  
int ClimateReader::Spherical2Grid(double longitude, double latitude, double depth, int &ilong, int &ilat, int &ilevel){
  if(verbose)
    cout<<"int spherical2grid("<<longitude<<", "<<latitude<<", "<<depth<<", int &ilong, int &ilat, int &ilevel)\n";
  
  ilat= (int)(((latitude+90.)/(space+1.E-7))+0.5);
  if(longitude<0)
    ilong=(int)(((longitude+360.)/(space+1.E-7))+0.5);
  else
    ilong=(int)((longitude/(space+1.E-7))+0.5);
  
  ilong=ilong%idim0;

  ilevel=levels.size()-1;
  for(int i=1;i<ilevel;i++){
    if(depth>levels[i]){
      ilevel = i-1;
      break;
    }
  }
  
  return 0;
}

void ClimateReader::VerboseOff(){
  verbose = false;
}

void ClimateReader::VerboseOn(){
  verbose = true;
}

// Global
ClimateReader ClimateReader_global;

// Fortran interface
extern "C" {
#define climatology_settimeseconds_fc F77_FUNC_(climatology_settimeseconds_c, CLIMATOLOGY_SET_TIME_SECONDS_C)
  void climatology_settimeseconds_fc(const double *_time){
    ClimateReader_global.SetTimeSeconds(*_time);  
    return;
  }
  
#define climatology_get_value_fc F77_FUNC_(climatology_get_value_c, CLIMATOLOGY_GET_VALUE_C)
  void climatology_get_value_fc(const char *name, const int *len, const double *x, const double *y, const double *z, double *value){
    *value = ClimateReader_global.GetValue(string(name, *len), *x, *y, *z);
    return;
  }

#define climatology_get_surface_value_fc F77_FUNC_(climatology_get_surface_value_c, CLIMATOLOGY_GET_SURFACE_VALUE_C)
  void climatology_get_surface_value_fc(const char *name, const int *len, const double *x, const double *y, const double *z, double *value){
    *value = ClimateReader_global.GetValue(string(name, *len), *x, *y, *z, 0);
    return;
  }
}
