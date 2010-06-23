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
#include "NetCDF_reader.h"

#include <fstream>

using namespace std;

NetCDF_reader::NetCDF_reader(const char *filename, bool _verbose){
#ifdef HAVE_NETCDF
  verbose = _verbose;
 
  if(verbose)
    cout<<"NetCDF_reader::NetCDF_reader()\n";
  
  NetCDFErrorCheckingOn();

  // Check that the file exists.
  fstream ncfile;
  ncfile.open(filename, ios::in);
  if(!ncfile.is_open()){
    cerr<<"ERROR: NetCDF file, "<<filename<<", cannot be opened. Does it exist? Have you read permission?\n";
    exit(-1);
  }
  ncfile.close();
  
  // Open the netCDF file.
  nc_open(filename, NC_NOWRITE, &ncid);
  if(ncerr!=NC_NOERR){
    cout.flush();
    cerr<<__FILE__<<", "<<__LINE__<<": Failed to open file "<<filename<<endl;
    exit(-1);
  }

  // Is this file:
  // - "GEBCO One Minute Grid"
  // - "COARDS" convention
  // - unknown
  data_source = unknown;
  {
    // Check if this is GEBCO 1 minuite data
    size_t lenp;
    if(nc_inq_attlen(ncid, NC_GLOBAL, "title", &lenp)==NC_NOERR){
      char *title = new char[lenp+1];
      nc_get_att_text(ncid, NC_GLOBAL, "title", title);
      title[lenp] = '\0';
      if(string("GEBCO One Minute Grid").compare(title)==0){
        data_source = gebco;
      }
      delete [] title;
    }
    
    if(data_source==unknown){
      if(nc_inq_attlen(ncid, NC_GLOBAL, "Conventions", &lenp)==NC_NOERR){
        char *Conventions = new char[lenp+1];
        nc_get_att_text(ncid, NC_GLOBAL, "Conventions", Conventions);
        Conventions[lenp] = '\0';
        if(string("COARDS").compare(Conventions)==0){
          data_source = coards;
        }else if(string("COARDS/CF-1.0").compare(Conventions)==0){
          data_source = coards_cf10;
        }
        delete [] Conventions;
      }
    }
  }
  
  if(verbose){
    cout<<"NetCDF dataset - ";
    if(data_source==gebco) cout<<"GEBCO"<<endl;
    if(data_source==coards) cout<<"COARDS"<<endl;
    if(data_source==coards_cf10) cout<<"COARDS/CF-1.0"<<endl;
    if(data_source==unknown) cout<<"unknown"<<endl;
  }
  
  // 
  long start=0;
  nc_type xtypep;                 /* variable type */
  int ndims;                      /* number of dims */

  // Get dimensions
  if(data_source==coards){
    int id = ncdimid(ncid, "x");
    if(ncerr!=NC_NOERR){
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<": ERROR - dimension variable x does not exist\n";
      exit(-1);
    }
    long x;
    ncdiminq(ncid, id, (char *)0, &x);
    
    id = ncdimid(ncid, "y");
    if(ncerr!=NC_NOERR){
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<": ERROR - dimension variable y does not exist\n";
      exit(-1);
    }
    long y;
    ncdiminq(ncid, id, (char *)0, &y);
    
    dimension[0] = x;
    dimension[1] = y;
    xysize = x*y;
  }else if(data_source==coards_cf10){
    NetCDFErrorCheckingOff();
    
    int id = ncdimid(ncid, "lon");
    if(ncerr!=NC_NOERR){
      ncerr=NC_NOERR;
      id = ncdimid(ncid, "x");
      if(ncerr!=NC_NOERR){
        cout.flush();
        cerr<<__FILE__<<", "<<__LINE__<<": ERROR - neither dimension variable lon nor x exist\n";
        exit(-1);
      }
    }
    long x;
    ncdiminq(ncid, id, (char *)0, &x);
    
    id = ncdimid(ncid, "lat");
    if(ncerr!=NC_NOERR){
      ncerr=NC_NOERR;
      id = ncdimid(ncid, "y");
      if(ncerr!=NC_NOERR){
        cout.flush();
        cerr<<__FILE__<<", "<<__LINE__<<": ERROR - dimension variable lat does not exist\n";
        exit(-1);
      }
    }
    long y;
    ncdiminq(ncid, id, (char *)0, &y);
    
    dimension[0] = x;
    dimension[1] = y;
    xysize = x*y;
    NetCDFErrorCheckingOn();
  }else{
    int varid = ncdimid(ncid, "xysize");
    if(ncerr!=NC_NOERR){
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<": ERROR - xysize doesn't exist\n";
      exit(-1);
    }
    ncdiminq(ncid, varid, (char *)0, &xysize);
    
    long side;
    varid = ncdimid(ncid, "side");
    if(ncerr!=NC_NOERR){
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<": ERROR - side doesn't exist\n";
      exit(-1);
    }
    ncdiminq(ncid, varid, (char *)0, &side);
    
    varid = ncvarid(ncid, "dimension");
    ncvarinq(ncid, varid, 0, &xtypep, &ndims, NULL, NULL);
    assert(ndims==1);
    
    if(xtypep==NC_SHORT){
      short var[2];
      ncvarget(ncid, varid, &start, &side, var);
      dimension[0] = var[0];
      dimension[1] = var[1];
    }else if(xtypep==NC_INT){
      ncvarget(ncid, varid, &start, &side, dimension);
    }else{
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<" - ERROR: unexpected data type for dimension\n";
      exit(-1);
    }
    if(verbose) 
      cout<<"dimension = "<<dimension[0]<<", "<<dimension[1]<<endl;    
  }
  
  if(verbose)
    cout<<"xysize = "<<xysize<<endl;
  
  if(data_source==coards){
    // x variable
    int varid = ncvarid(ncid, "x");
    if(ncerr!=NC_NOERR){
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<": ERROR - variable x does not exist.\n";
      exit(-1);
    }
    
    int len;
    ncattinq(ncid, varid, "actual_range", &xtypep, &len);
    assert(len==2);
    if(xtypep==NC_FLOAT){ // single precision floating point number
      float var[2];
      ncattget(ncid, varid, "actual_range", var);
      x_range[0] = var[0];
      x_range[1] = var[1];
    }else if(xtypep==NC_DOUBLE){         // double precision floating point number
      ncattget(ncid, varid, "actual_range", x_range);
    }else{
      cerr<<__FILE__<<", "<<__LINE__<<" - ERROR: unexpected data type for x_range\n";
      exit(-1);
    }

    // y variable
    varid = ncvarid(ncid, "y");
    if(ncerr!=NC_NOERR){
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<": ERROR - variable y does not exist.\n";
      exit(-1);
    }
    
    ncattinq(ncid, varid, "actual_range", &xtypep, &len);
    assert(len==2);
    if(xtypep==NC_FLOAT){ // single precision floating point number
      float var[2];
      ncattget(ncid, varid, "actual_range", var);
      y_range[0] = var[0];
      y_range[1] = var[1];
    }else if(xtypep==NC_DOUBLE){         // double precision floating point number
      ncattget(ncid, varid, "actual_range", y_range);
    }else{
      cerr<<__FILE__<<", "<<__LINE__<<" - ERROR: unexpected data type for y_range\n";
      exit(-1);
    }
  }else if(data_source==coards_cf10){
    NetCDFErrorCheckingOff();
    // lon variable
    int varid = ncvarid(ncid, "lon");
    if(ncerr!=NC_NOERR){
      ncerr=NC_NOERR;
      varid = ncvarid(ncid, "x");
      if(ncerr!=NC_NOERR){
        cout.flush();
        cerr<<__FILE__<<", "<<__LINE__<<": ERROR - variable lon does not exist.\n";
        exit(-1);
      }
    }
    int len;
    ncattinq(ncid, varid, "actual_range", &xtypep, &len);
    assert(len==2);
    if(xtypep==NC_FLOAT){ // single precision floating point number
      float var[2];
      ncattget(ncid, varid, "actual_range", var);
      x_range[0] = var[0];
      x_range[1] = var[1];
    }else if(xtypep==NC_DOUBLE){         // double precision floating point number
      ncattget(ncid, varid, "actual_range", x_range);
    }else{
      cerr<<__FILE__<<", "<<__LINE__<<" - ERROR: unexpected data type for x_range\n";
      exit(-1);
    }
    // y variable
    varid = ncvarid(ncid, "lat");
    if(ncerr!=NC_NOERR){
      ncerr=NC_NOERR;
      varid = ncvarid(ncid, "y");
      if(ncerr!=NC_NOERR){
        cout.flush();
        cerr<<__FILE__<<", "<<__LINE__<<": ERROR - variable lat does not exist.\n";
        exit(-1);
      }
    }
    ncattinq(ncid, varid, "actual_range", &xtypep, &len);
    assert(len==2);
    if(xtypep==NC_FLOAT){ // single precision floating point number
      float var[2];
      ncattget(ncid, varid, "actual_range", var);
      y_range[0] = var[0];
      y_range[1] = var[1];
    }else if(xtypep==NC_DOUBLE){         // double precision floating point number
      ncattget(ncid, varid, "actual_range", y_range);
    }else{
      cerr<<__FILE__<<", "<<__LINE__<<" - ERROR: unexpected data type for lat_range\n";
      exit(-1);
    }
    NetCDFErrorCheckingOn();
  }else{    
    // -
        long side;
        int varid = ncdimid(ncid, "side");
    if(ncerr!=NC_NOERR){
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<": ERROR - side doesn't exist\n";
      exit(-1);
    }
    ncdiminq(ncid, varid, (char *)0, &side);
    
    varid = ncvarid(ncid, "x_range");
    ncvarinq(ncid, varid, 0, &xtypep, &ndims, NULL, NULL);
    assert(ndims==1);

    if(xtypep==NC_FLOAT){ // single precision floating point number
      float var[2];
      ncvarget(ncid, varid, &start, &side, var);
      x_range[0] = var[0];
      x_range[1] = var[1];
    }else if(xtypep==NC_DOUBLE){         // double precision floating point number
      ncvarget(ncid, varid, &start, &side, x_range);
    }else{
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<" - ERROR: unexpected data type for x_range\n";
      exit(-1);
    }
    if(verbose)
      cout<<"x_range = "<<x_range[0]<<", "<<x_range[1]<<endl;
    
    // -
    varid = ncvarid(ncid, "y_range");
    ncvarinq(ncid, varid, 0, &xtypep, &ndims, NULL, NULL);
    assert(ndims==1);
    
    if(xtypep==NC_FLOAT){ // single precision floating point number
      float var[2];
      ncvarget(ncid, varid, &start, &side, var);
      y_range[0] = var[0];
      y_range[1] = var[1];
    }else if(xtypep==NC_DOUBLE){ // double precision floating point number
      ncvarget(ncid, varid, &start, &side, y_range);
    }else{
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<" - ERROR: unexpected data type for y_range\n";
      exit(-1);
    }
    if(verbose)
      cout<<"y_range = "<<y_range[0]<<", "<<y_range[1]<<endl;
  }

  spacing[0] = (x_range[1] - x_range[0])/(dimension[0]-1);
  spacing[1] = (y_range[1] - y_range[0])/(dimension[1]-1);
  
  return;
#else
  cerr << "No NetCDF support" << endl;
  exit(-1);
#endif
}

NetCDF_reader::~NetCDF_reader(){
#ifdef HAVE_NETCDF
  // Close the netCDF file.
  ncclose(ncid);
#else
  cerr << "No NetCDF support" << endl;
#endif
}

int NetCDF_reader::GetDimension(int &ilen, int &jlen) const{
  ilen = dimension[0];
  jlen = dimension[1];
  return 0;
}

int NetCDF_reader::GetSpacing(double &dx, double &dy) const{
  dx = spacing[0];
  dy = spacing[1];
  return 0;
}

int NetCDF_reader::GetXRange(double &xmin, double &xmax) const{
  xmin = x_range[0];
  xmax = x_range[1];
  return 0;
}

int NetCDF_reader::GetYRange(double &ymin, double &ymax) const{
  ymin = y_range[0];
  ymax = y_range[1];
  return 0;
}

bool NetCDF_reader::IsColatitude() const{
  if(data_source==coards_cf10)
    return false;
  // else
  return data_source!=coards;
}

int NetCDF_reader::Read(vector<double> &z) const{
#ifdef HAVE_NETCDF
  if(verbose)
    cout<<"int NetCDF_reader::Read(vector<double> &z) const\n";

  long start[]={0,0}, count[2];
  if((data_source==coards)||(data_source==coards_cf10)){
    count[0] = dimension[1];
    count[1] = dimension[0];
  }else{
    count[0] = xysize;
    count[1] = 0;
  }

  nc_type xtypep;                 /* variable type */
  int ndims;                      /* number of dims */
  int dims[MAX_VAR_DIMS];         /* variable shape */
  
  int varid = ncvarid(ncid, "z");
  ncvarinq(ncid, varid, 0, &xtypep, &ndims, dims, NULL);

  z.resize(xysize);
  if(xtypep==NC_BYTE){ // signed 1 byte integer
    vector<char> var(xysize);
    ncvarget(ncid, varid, start, count, &(var[0]));
    for(long i=0;i<xysize;i++)
      z[i] = var[i];
  }else if(xtypep==NC_CHAR){ // ISO/ASCII character
    vector<char> var(xysize);
    ncvarget(ncid, varid, start, count, &(var[0]));
    for(long i=0;i<xysize;i++)
      z[i] = var[i];
  }else if(xtypep==NC_SHORT){ // signed 2 byte integer
    vector<short> var(xysize);
    ncvarget(ncid, varid, start, count, &(var[0]));
    for(long i=0;i<xysize;i++)
      z[i] = var[i];
  }else if(xtypep==NC_INT){ // signed 4 byte integer
    vector<int> var(xysize);
    ncvarget(ncid, varid, start, count, &(var[0]));
    for(long i=0;i<xysize;i++)
      z[i] = var[i];
  }else if(xtypep==NC_FLOAT){ // single precision floating point number
    vector<float> var(xysize);
    ncvarget(ncid, varid, start, count, &(var[0]));
    for(long i=0;i<xysize;i++){
      z[i] = var[i];
    }
  }else if(xtypep==NC_DOUBLE){ // double precision floating point number
    ncvarget(ncid, varid, start, count, &(z[0]));
  }else{
    cerr<<"ERROR: unknown data type for z\n";
    exit (-1);
  }

  return 0;
#else
  cerr << "No NetCDF support" << endl;
  exit(-1);
#endif
}

void NetCDF_reader::NetCDFErrorCheckingOff(){
#ifdef HAVE_NETCDF
  ncopts = 0;
#else
  cerr << "No NetCDF support" << endl;
  exit(-1);
#endif
}

void NetCDF_reader::NetCDFErrorCheckingOn(){
#ifdef HAVE_NETCDF
  // Make NetCDF errors verbose and fatal
  ncopts = NC_VERBOSE | NC_FATAL;
#else
  cerr << "No NetCDF support" << endl;
  exit(-1);
#endif
}

void NetCDF_reader::VerboseOff(){
  verbose = false;
}

void NetCDF_reader::VerboseOn(){
  verbose = true;
}
