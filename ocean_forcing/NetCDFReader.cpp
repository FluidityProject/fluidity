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

#include "NetCDFReader.h"

#include <fstream>

using namespace std;

NetCDFReader::NetCDFReader() : fileOpen(false) {
#ifdef HAVE_NETCDF
  VerboseOff();

  // Make NetCDF errors verbose and fatal
  ncopts = NC_VERBOSE;
#endif
}

NetCDFReader::NetCDFReader(const char *filename){
#ifdef HAVE_NETCDF
  VerboseOff();

  // Make NetCDF errors verbose and fatal
  ncopts = NC_VERBOSE;

  SetFile(filename);
#endif
}

const NetCDFReader& NetCDFReader::operator=(const NetCDFReader &in){
  if(verbose)
    cout<<"const NetCDFReader& NetCDFReader::operator=(const NetCDFReader &in)\n";

  ncid = in.ncid;
  fileOpen = in.fileOpen;

  for(int i=0;i<2;i++){
    x_range[i] = in.x_range[i];
    y_range[i] = in.y_range[i];

    spacing[i] = in.spacing[i];
    dimension[i] = in.dimension[i];
  }

  verbose = in.verbose;

  return *this;
}

void NetCDFReader::Close(){
  if(verbose)
    cout<<"void NetCDFReader::Close()\n";
#ifdef HAVE_NETCDF
  if(fileOpen){
    // Close the netCDF file.
    ncclose(ncid);  
  }
  fileOpen = false;
#endif
}

void NetCDFReader::SetFile(const char *filename){
  if(verbose)
    cout<<"void NetCDFReader::SetFile(const char *filename)\n";
#ifdef HAVE_NETCDF
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

  fileOpen = true;

  GetGrid();
#else
  cerr<<"WARNING ("<<__FILE__<<", "<<__LINE__<<"): "
      <<"No NetCDF support available. Recompilation required."<<endl;
#endif
  return;
}

void NetCDFReader::GetGrid(){
  if(verbose)
    cout<<"NetCDFReader::GetGrid()\n";
#ifdef HAVE_NETCDF
  // Get dimensions -- longitude
  int id = ncdimid(ncid, "longitude");
  if(ncerr!=NC_NOERR){
    cout.flush();
    cerr<<__FILE__<<", "<<__LINE__<<": ERROR - dimension variable \"longitude\" does not exist\n";
    exit(-1);
  }
  ncdiminq(ncid, id, (char *)0, &(dimension[0]));
    
  // Get dimensions -- latitude
  id = ncdimid(ncid, "latitude");
  if(ncerr!=NC_NOERR){
    cout.flush();
    cerr<<__FILE__<<", "<<__LINE__<<": ERROR - dimension variable \"latitude\" does not exist\n";
    exit(-1);
  }
  ncdiminq(ncid, id, (char *)0, &(dimension[1]));

  // Longitude range
  int varid = ncvarid(ncid, "longitude");
  if(ncerr!=NC_NOERR){
    cout.flush();
    cerr<<__FILE__<<", "<<__LINE__<<": ERROR - variable \"longitude\" does not exist.\n";
    exit(-1);
  }
  
  nc_type xtypep;
  int len;
  ncattinq(ncid, varid, "actual_range", &xtypep, &len);
  assert(len==2);
  if(xtypep==NC_DOUBLE){
    ncattget(ncid, varid, "actual_range", x_range);
  }else if(xtypep==NC_FLOAT){
    float val[2];
    ncattget(ncid, varid, "actual_range", val);
    x_range[0] = val[0];
    x_range[1] = val[1];
  }else{
    cerr<<__FILE__<<", "<<__LINE__<<" - ERROR: unexpected data type for x_range\n";
    exit(-1);
  }

  // Latitude range
  varid = ncvarid(ncid, "latitude");
  if(ncerr!=NC_NOERR){
    cout.flush();
    cerr<<__FILE__<<", "<<__LINE__<<": ERROR - variable latitude does not exist.\n";
    exit(-1);
  }

  ncattinq(ncid, varid, "actual_range", &xtypep, &len);
  assert(len==2);
  if(xtypep==NC_DOUBLE){
    ncattget(ncid, varid, "actual_range", y_range);
  }else if(xtypep==NC_FLOAT){
    float val[2];
    ncattget(ncid, varid, "actual_range", val);
    y_range[0] = val[0];
    y_range[1] = val[1];
  }else{
    cerr<<__FILE__<<", "<<__LINE__<<" - ERROR: unexpected data type for y_range\n";
    exit(-1);
  }
  
  spacing[0] = (x_range[1] - x_range[0])/(dimension[0]-1);
  spacing[1] = (y_range[1] - y_range[0])/(dimension[1]-1);
#endif
  return;
}

NetCDFReader::~NetCDFReader(){
  // Close();
}

int NetCDFReader::GetDimension(int &ilen, int &jlen) const{
  if(verbose)
    cout<<"int NetCDFReader::GetDimension(int &ilen, int &jlen) const\n";

  ilen = dimension[0];
  jlen = dimension[1];

  return 0;
}

int NetCDFReader::GetSpacing(double &dx, double &dy) const{
  if(verbose)
    cout<<"int NetCDFReader::GetSpacing(double &dx, double &dy) const\n";

  dx = spacing[0];
  dy = spacing[1];

  return 0;
}

int NetCDFReader::GetXRange(double &xmin, double &xmax) const{
  if(verbose)
    cout<<"int NetCDFReader::GetXRange(double &xmin, double &xmax) const\n";

  xmin = x_range[0];
  xmax = x_range[1];

  return 0;
}

int NetCDFReader::GetYRange(double &ymin, double &ymax) const{
  if(verbose)
    cout<<"int NetCDFReader::GetYRange(double &ymin, double &ymax) const\n";

  ymin = y_range[0];
  ymax = y_range[1];

  return 0;
}

int NetCDFReader::Read(string varname, vector<double> &var){
  if(verbose)
    cout<<"int NetCDFReader::Read("<<varname<<", vector<double> &) const\n";
#ifdef HAVE_NETCDF
  nc_type xtypep;                 /* variable type */
  int ndims;                      /* number of dims */
  int dims[MAX_VAR_DIMS];         /* variable shape */
  
  int varid = ncvarid(ncid, varname.c_str());
  ncvarinq(ncid, varid, 0, &xtypep, &ndims, dims, NULL);

  long start[]={0,0}, count[2];
  count[0] = dimension[1];
  count[1] = dimension[0];

  long len = dimension[0]*dimension[1];
  var.resize(len);
  if(xtypep==NC_FLOAT){ // single precision floating point number
    vector<float> fvar(len);
    ncvarget(ncid, varid, start, count, &(fvar[0]));
    for(long i=0;i<len;i++)
      var[i] = fvar[i];
    fvar.clear();
  }else if(xtypep==NC_DOUBLE){ // double precision floating point number
    ncvarget(ncid, varid, start, count, &(var[0]));
  }else{
    cerr<<"ERROR: unknown data type for z\n";
    exit (-1);
  }
#endif
  return 0;
}

void NetCDFReader::VerboseOff(){
  verbose = false;
}

void NetCDFReader::VerboseOn(){
  verbose = true;
}
