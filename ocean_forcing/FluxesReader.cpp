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
#include "FluxesReader.h"
#include <limits>
#include <string.h>

#ifdef __SUN
double round(double x)
{
  return floor(x+0.5);
}
#endif

using namespace std;

FluxesReader::FluxesReader(){
  verbose = false;
  MyRank = 0;
  NProcs = 1;
#ifdef HAVE_MPI
  int init_flag;
  MPI_Initialized(&init_flag);
  if(init_flag){
    int MyRank, NProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
    MPI_Comm_size(MPI_COMM_WORLD, &NProcs);
  }
#endif
  if(NProcs>1){
    NParts=NProcs*10;
  }

  modified = false;
  time_set = -1.0;

  calendar = NULL;
  return;
}

FluxesReader::~FluxesReader(){
  if(calendar!=NULL)
    delete calendar;
}

// Add field to the list of fields to be interpolated - the order that
// they are added is the order in which they will be returned on
// interpolation.
void FluxesReader::AddFieldOfInterest(std::string field_name){
  if(verbose)
    cout<<"void FluxesReader::AddFieldOfInterest("<<field_name<<")\n";
  
  for(deque<string>::const_iterator ifield=fields_of_interest.begin(); ifield!=fields_of_interest.end(); ifield++)
    if(field_name==*ifield)
      cerr<<"ERROR: void FluxesReader::AddFieldOfInterest("<<field_name<<")\n Field added multiple times\n";
  
  fields_of_interest.push_back(field_name);
  
  modified = true;
  
  return;
}

void FluxesReader::ClearFields(){
  if(verbose)
    cout<<"void FluxesReader::ClearFields()\n";
  
  fields_of_interest.clear();
  return;
}

bool FluxesReader::Enabled() const{
  if(verbose)
    cout<<"void FluxesReader::Enabled() const\n";
  
  return ERA_data_files.size()!=0;
}

pair<size_t, size_t> FluxesReader::GetInterval(double value, const map<double, int>&lut) const{
  map<double, int>::const_iterator lower_bound = lut.lower_bound(value);
  assert(lower_bound!=lut.end());

  map<double, int>::const_iterator upper_bound = lower_bound;
  if(lower_bound!=lut.begin())
    lower_bound--;

  return pair<size_t, size_t>(lower_bound->second, upper_bound->second);
}

double FluxesReader::GetScalar(string scalar, double xlong, double ylat){
  if(verbose)
    cout<<"void FluxesReader::GetScalar("<<scalar<<", "<<xlong<<", "<<ylat<<")\n";
  
  // Ensure that time has been set and is within range.
  if((time_set<*time.begin())||(time_set>*time.rbegin())){
    cerr<<"ERROR: int FluxesReader::GetScalar( ... )\n"
        <<"time range is "<<*time.begin()<<" --> "<<*time.rbegin()<<endl;
    exit(-1);
  }
  
  if(modified)
    Update();
  
  // Ensure that the input is sensiable
  while(xlong<0)
    xlong+=360.0;

  while(xlong>=360.0)
    xlong-=360.0;
  
  assert(ylat>=-90);
  assert(ylat<=90);
  
  map<double, int>::const_iterator long1 = lut_longitude.lower_bound(xlong);
  assert(long1!=lut_longitude.end());

  map<double, int>::const_iterator long0 = long1;
  if((long0!=lut_longitude.begin())&&(fabs(long0->first-xlong)>0.001))
    long0--;
  size_t i0=long0->second;
  size_t i1=long1->second;
  if(verbose)
    cout<<"longitude "<<xlong<<" Use data from "<<long0->first<<" --> "<<long1->first<<endl;
  
  map<double, int>::const_iterator lat1 = lut_latitude.lower_bound(ylat);
  assert(lat1!=lut_latitude.end());
  if(verbose)
    cout<<"latitude lower bound "<<ylat<<" --> "<<lat1->first<<endl;

  map<double, int>::const_iterator lat0 = lat1;
  if((lat0!=lut_latitude.begin())&&(fabs(lat0->first-ylat)>0.001))
    lat0--;
  
  size_t j0=lat0->second;
  size_t j1=lat1->second;
  
  if(verbose)
    cout<<"Selecting grid ("<<i0<<", "<<j0<<"), ("<<i1<<", "<<j1<<")\n";

  assert(fields.begin()!=fields.end());
  size_t t0 = fields.begin()->first;
  size_t t1 = fields.rbegin()->first;
  
  double values[2];
  size_t time_level=0;
  size_t nvalues=0, stride=longitude.size();

  for(map<int, map<string, vector<double> > >::const_iterator itime=fields.begin(); itime!=fields.end(); itime++){
    const vector<double> &fld=itime->second.find(scalar)->second;
    
    if(i0==i1){ // No interpolation along longitude
      if(verbose)
        cout<<"Case 1\n";
      
      if(j0==j1){
        double x = fld[j0*stride+i0];
        values[time_level] = x;
      
      }else{
      
        double x0 = fld[j0*stride+i0];
        double x1 = fld[j1*stride+i0];

        double x = SolveLine(x0, x1, lat0->first, lat1->first, ylat);
        
        values[time_level] = x;
      }
    
    } else if(j0==j1){ // No interpolation along latitude
      if(verbose)
        cout<<"Case 2\n";
      
      double x0 = fld[j0*stride+i0];
      double x1 = fld[j0*stride+i1];

      double long_temp0=long0->first;
      double long_temp1=long1->first;
      if (long_temp1 < long_temp0) {
        long_temp1 += 360.0;
      }
       
      double x = SolveLine(x0, x1, long_temp0, long_temp1, xlong);

      values[time_level] = x;
    }else{ // Bi-linear interpolation
      if(verbose)
        cout<<"Case 3 -- "<<fields.size()<<endl;
      
      double x00 = fld[j0*longitude.size()+i0];
      double x10 = fld[j0*longitude.size()+i1];
      double x01 = fld[j1*longitude.size()+i0];
      double x11 = fld[j1*longitude.size()+i1];

      // Calculating the latitude differences between the point of interpolation and the data points
      double dl0=xlong-longitude[i0];
      double dl1=xlong-longitude[i1];

      // xlong should be greater than longitude[i0]; fix this if not
      if (dl0 < 0.0) {
        dl0+=360.0;
      }

      // xlong should be less than longitude[i1]; fix this if not
      if (dl1 > 0.0) {
        dl1-=360.0;
      }

      double x =  (x00*(dl1)*(ylat-latitude[j1]) -
                   x10*(dl0)*(ylat-latitude[j1]) -
                   x01*(dl1)*(ylat-latitude[j0]) + 
                   x11*(dl0)*(ylat-latitude[j0]))
                   /(dlong*dlat);
      values[time_level] = x;
    }
    time_level++;
  }
  
  double val;
  if(fields.size()==1){
    val = values[0];
  }else{
    assert(fields.size()==2);
    val = SolveLine(values[0], values[1], (double)(time[t0]), (double)(time[t1]), time_set);
  }

  if(verbose)
    cout<<"Value = "<<val/GetScaleFactor(scalar)<<endl;

  return val/GetScaleFactor(scalar);
}

int FluxesReader::GetScalars(double xlong, double ylat, double *scalars){
  if(verbose)
    cout<<"int FluxesReader::GetScalars("<<xlong<<", "<<ylat<<", scalars)\n";
  
  if(modified)
    Update();
  
  // Ensure that the input is sensible
  while(xlong<0)
    xlong+=360.0;
  while(xlong>=360.0)
    xlong-=360.0;
  
  assert(ylat>=-90);
  assert(ylat<=90);
  
  map<double, int>::const_iterator long1 = lut_longitude.lower_bound(xlong);
  assert(long1!=lut_longitude.end());
  map<double, int>::const_iterator long0 = long1;
  if((long0!=lut_longitude.begin())&&(fabs(long0->first-xlong)>0.001))
    long0--;
  size_t i0=long0->second;
  size_t i1=long1->second;
  if(verbose)
    cout<<"longitude "<<xlong<<" will be calculated from "<<long0->first<<" --> "<<long1->first<<endl;

  
  map<double, int>::const_iterator lat1 = lut_latitude.lower_bound(ylat);
  assert(lat1!=lut_latitude.end());
  map<double, int>::const_iterator lat0 = lat1;
  if((lat0!=lut_latitude.begin())&&(fabs(lat0->first-ylat)>0.001))
    lat0--;
  if(verbose)
    cout<<"latitude "<<ylat<<" will be calculated from "<<lat0->first<<" --> "<<lat1->first<<endl;

  
  size_t j0=lat0->second;
  size_t j1=lat1->second;
  
  if(verbose)
    cout<<"Selecting grid ("<<i0<<", "<<j0<<"), ("<<i1<<", "<<j1<<")\n";

  assert(fields.begin()!=fields.end());
  size_t t0 = fields.begin()->first;
  size_t t1 = fields.rbegin()->first;

  double values[1024]; // This is way more than the number of
  // variables calculated in ERA-40
  size_t nvalues=0, stride=longitude.size();
  
  for(map<int, map<string, vector<double> > >::const_iterator itime=fields.begin(); itime!=fields.end(); itime++){
    for(deque<string>::const_iterator ifield=fields_of_interest.begin(); ifield!=fields_of_interest.end(); ifield++){
      const vector<double> &fld=itime->second.find(*ifield)->second;
      
      if(i0==i1){ // No interpolation along longitude
        if(verbose)
          cout<<"Case 1\n";
        if(j0==j1){
          double x = fld[j0*stride+i0];
          values[nvalues] = x;
        }else{
          double x0 = fld[j0*stride+i0];
          double x1 = fld[j1*stride+i0];
          
          double x = SolveLine(x0, x1, lat0->first, lat1->first, ylat);
          values[nvalues] = x;
        }
        // }
      }else if(j0==j1){ // No interpolation along latitude
        if(verbose)
          cout<<"Case 2\n";
        double x0 = fld[j0*stride+i0];
        double x1 = fld[j0*stride+i1];

        double long_temp0=long0->first;
        double long_temp1=long1->first;

        if (long_temp1 < long_temp0) {
          long_temp1 += 360.0;
        }

        
        double x = SolveLine(x0, x1, long_temp0, long_temp1, xlong);
        values[nvalues] = x;
        // }
      }else{ // Bi-linear interpolation
        if(verbose) {
          cout<<"Case 3 -- "<<fields.size()<<endl;
        }
        double x00 = fld[j0*stride+i0];
        double x10 = fld[j0*stride+i1];
        double x01 = fld[j1*stride+i0];
        double x11 = fld[j1*stride+i1];

        // Calculating the latitude differences between the point of interpolation and the data points
        double dl0=xlong-longitude[i0];
        double dl1=xlong-longitude[i1];

        // xlong should be greater than longitude[i0]; fix this if not
        if (dl0 < 0.0) {
          dl0+=360.0;
        }

        // xlong should be less than longitude[i1]; fix this if not
        if (dl1 > 0.0) {
          dl1-=360.0;
        }

        double x =  (x00*(dl1)*(ylat-latitude[j1]) -
                     x10*(dl0)*(ylat-latitude[j1]) -
                     x01*(dl1)*(ylat-latitude[j0]) + 
                     x11*(dl0)*(ylat-latitude[j0]))
                     /(dlong*dlat);
        values[nvalues] = x;
      }
      nvalues++;
    }
  }
  
  if(fields.size()==1){
    assert(nvalues==fields_of_interest.size());
    for(size_t i=0;i<nvalues; i++)
      scalars[i] = values[i];
  }else{
    assert(nvalues==(2*fields_of_interest.size()));
    double *val2=&(values[fields_of_interest.size()]);
    nvalues/=2;
    for(size_t i=0;i<nvalues; i++)
      scalars[i] = SolveLine(values[i], val2[i], (double)(time[t0]), (double)(time[t1]), time_set);
  }
  
  for(size_t i=0;i<nvalues; i++){
    scalars[i]/=GetScaleFactor(fields_of_interest[i]);
    if(verbose)
      cout<<"Value "<<i<<": "<<scalars[i]<<endl;
  }
  return 0;
}

double FluxesReader::GetScaleFactor(std::string variable_name) const{
  return 1.0;  
}

int FluxesReader::Read(int time_index){
  if(verbose)
    cout<<"int Read("<<time_index<<")\n";
  
  // Make sure that a file has been registered
  if(ERA_data_files.size()==0){
    cerr<<"ERROR: int Read("<<time_index<<") -- no file has been registered\n";
    exit(-1);
  }
  
#ifdef HAVE_LIBNETCDF
  for(deque<string>::iterator ifield=fields_of_interest.begin(); ifield!=fields_of_interest.end(); ifield++){
    // Check if we already have a copy of this field
    if(fields.find(time_index)!=fields.end())
      if(fields[time_index].find(*ifield)!=fields[time_index].end())
        if(fields[time_index][*ifield].empty())
          continue;
    
    long id = ncvarid(ncid, ifield->c_str());
    nc_type xtypep;                 /* variable type */
    int ndims;                      /* number of dims */
    int dims[MAX_VAR_DIMS];         /* variable shape */
    int natts;                      /* number of attributes */
    
    ncvarinq(ncid, id, 0, &xtypep, &ndims, dims, &natts);
    assert(xtypep==NC_SHORT);
    
    long start[]={time_index, 0, 0}, count[]={1, spec["latitude"], spec["longitude"]};
    size_t len = spec["latitude"]*spec["longitude"];
    vector<short> field(len);
    ncvarget(ncid, id, start, count, &(field[0]));
    
    // Attributes to read
    int err;
    double scale_factor;
    err = nc_get_att_double(ncid, id, "scale_factor", &scale_factor);
    
    double add_offset;
    err = nc_get_att_double(ncid, id, "add_offset", &add_offset);
    
    short fill_value;
    err = nc_get_att_short(ncid, id, "_FillValue", &fill_value);
    if (err != NC_NOERR) fill_value = std::numeric_limits<short>::quiet_NaN();
    
    short missing_value;
    err = nc_get_att_short(ncid, id, "missing_value", &missing_value);
    if (err != NC_NOERR) missing_value = std::numeric_limits<short>::quiet_NaN();
    
    if(verbose){
      cout<<*ifield<<" scale_factor = "<<scale_factor<<endl
          <<*ifield<<" add_offset   = "<<add_offset<<endl
          <<*ifield<<"  _FillValue   = "<<fill_value<<endl
          <<*ifield<<"  _missing_value   = "<<missing_value<<endl;
    }
    
    fields[time_index][*ifield].resize(len);
    for(size_t i=0; i<len; i++){
      if((field[i]==fill_value)||
         (field[i]==missing_value)){
            cerr<<"FluxesReader: Missing value in input file:\n";
            cerr<<"Time: "<<time_index<<". NCID: "<<ncid<<endl;
            exit(-1);
      }else{
        fields[time_index][*ifield][i] = scale_factor*field[i] + add_offset;
      }
    }
  }
#else
  cerr<<"ERROR: No fluxes support compiled\n";
  exit(-1);
#endif
  return 0;
}

// Returns 0 on success, negitive if error.
int FluxesReader::RegisterDataFile(string file){
#ifdef HAVE_LIBNETCDF
  if(verbose)
    cout<<"int FluxesReader::RegisterDataFile("<<file<<")\n";
  
  // Check what has gone before
  if(ERA_data_files.size()==1){
    assert(ERA_data_files[0]==file);
    return 0;
  }
  // else
  
  ERA_data_files.push_back(file);
  
  //
  // Get the basic information from the fluxes file
  //
  
  // Make fluxes errors verbose and fatal
  ncopts = NC_VERBOSE | NC_FATAL;
  
  // Open the netCDF file.
  if(file.c_str()[0]!='/')
    file = string(getenv("PWD"))+string("/")+file;
  ncid = ncopen(file.c_str(), NC_NOWRITE);
  if(ncid==-1){
    cerr<<"ERROR: could not open netcdf file: "<<file<<endl;
    FLAbort("int FluxesReader::RegisterDataFile(string file)", __FILE__, __LINE__);
  }
  
  // If this is the first time we're opened one of these files then we
  // need to store it's basic specification (dimensions
  // etc). Otherwise we read in the new specification and reassure
  // ourselves that it's the same as the original.
  
  //
  // Get dimensions
  //
  long nlatitude;
  long nlongitude;
  { // longitude
    int id = ncdimid(ncid, "longitude");
    ncdiminq(ncid, id, (char *)0, &nlongitude);
    
    assert(spec.find("longitude")==spec.end());
    spec["longitude"] = nlongitude;

    if(verbose)
      cout<<"("<<__FILE__<<", "<<__LINE__<<"): longitude = "<<nlongitude<<endl;
  }
  
  { // latitude
    int id = ncdimid(ncid, "latitude");
    ncdiminq(ncid, id, (char *)0, &nlatitude);

    assert(spec.find("latitude")==spec.end());
    spec["latitude"] = nlatitude;
    dlat = 180.0/(nlatitude-1);

    if(verbose)
      cout<<"("<<__FILE__<<", "<<__LINE__<<"): latitude = "<<nlatitude<<endl;
  }

  { // time
    int id = ncdimid(ncid, "time");
    long ntime;
    ncdiminq(ncid, id, (char *)0, &ntime);

    assert(spec.find("time")==spec.end());
    spec["time"] = ntime;
    
    if(verbose)
      cout<<"("<<__FILE__<<", "<<__LINE__<<"): time = "<<ntime<<endl;
  }

  // Get variables
  nc_type xtypep;                 /* variable type */
  int ndims;                      /* number of dims */
  int dims[MAX_VAR_DIMS];         /* variable shape */
  int natts;                      /* number of attributes */
  
  { // longitude
    long id = ncvarid(ncid, "longitude");
    ncvarinq(ncid, id, 0, &xtypep, &ndims, dims, &natts);
    
    long start=0, count=spec["longitude"];
    longitude.resize(count);
    
    if(xtypep==NC_INT){
      vector<int> longitude_tmp(count);
      ncvarget(ncid, id, &start, &count, &(longitude_tmp[0]));
      for(size_t i=0;i<(size_t)count;i++)
        longitude[i] = longitude_tmp[i];
    }else if(xtypep==NC_FLOAT){
      vector<float> longitude_tmp(count);
      ncvarget(ncid, id, &start, &count, &(longitude_tmp[0]));
      for(size_t i=0;i<(size_t)count;i++)
        longitude[i] = longitude_tmp[i];
    }else if(xtypep==NC_DOUBLE){
      ncvarget(ncid, id, &start, &count, &(longitude[0]));
    } else if(xtypep==NC_SHORT) {
       double scale_factor;
       nc_get_att_double(ncid, id, "scale_factor", &scale_factor); 
       double add_offset;
       nc_get_att_double(ncid, id, "add_offset", &add_offset);
       vector<short> longitude_tmp(count);
       ncvarget(ncid, id, &start, &count, &(longitude_tmp[0]));
       for (size_t i=0;i<(size_t)count; i++)
           longitude[i] = longitude_tmp[i]*scale_factor + add_offset;
    }else{
      cerr<<"ERROR: ("<<__FILE__<<"): longitude has unexpected type.\n";
      exit(-1);
    }
    // working out Dlong needs to be done after we know how large the grid is...
    double latLength = abs(longitude[0]-longitude[count-1]);
    if (latLength<360.0) {
        dlong = latLength/(nlongitude-1);
    } else {
        dlong = latLength/nlongitude;
    }
    // Note that longitude 0 is repeated at 360.0 - therefore no -1!


    for(size_t i=0;i<(size_t)count;i++)
      lut_longitude[longitude[i]] = i;

    if(lut_longitude.begin()->first<0.001)
      lut_longitude[360.0] = lut_longitude.begin()->second;
    
    if(verbose){
      cout<<"longitude =";
      for(vector<double>::iterator ilong=longitude.begin(); ilong!=longitude.end(); ilong++)
        cout<<*ilong<<" ";
      cout<<endl;
    }
  }

  { // latitude
    long id = ncvarid(ncid, "latitude");
    ncvarinq(ncid, id, 0, &xtypep, &ndims, dims, &natts);
    
    long start=0, count=spec["latitude"];
    latitude.resize(count);

    if(xtypep==NC_INT){
      vector<int> latitude_tmp(count);
      ncvarget(ncid, id, &start, &count, &(latitude_tmp[0]));
      for(size_t i=0;i<(size_t)count;i++)
        latitude[i] = latitude_tmp[i];
    }else if(xtypep==NC_FLOAT){
      vector<float> latitude_tmp(count);
      ncvarget(ncid, id, &start, &count, &(latitude_tmp[0]));
      for(size_t i=0;i<(size_t)count;i++)
        latitude[i] = latitude_tmp[i];
    }else if(xtypep==NC_DOUBLE){
      ncvarget(ncid, id, &start, &count, &(latitude[0]));
    } else if (xtypep==NC_SHORT){
       double scale_factor;
       nc_get_att_double(ncid, id, "scale_factor", &scale_factor); 
       double add_offset;
       nc_get_att_double(ncid, id, "add_offset", &add_offset);
       vector<short> latitude_tmp(count);
       ncvarget(ncid, id, &start, &count, &(latitude_tmp[0]));
       for (size_t i=0;i<(size_t)count; i++){
           latitude[i] = latitude_tmp[i]*scale_factor + add_offset;
       }
    }else{
      cerr<<"ERROR: ("<<__FILE__<<"): latitude has unexpected type.\n";
      exit(-1);
    }
    // working out Dlat needs to be done after we know how large the grid is...
    dlat = (abs(latitude[0]-latitude[count-1]))/(nlatitude-1);

    
    for(size_t i=0;i<(size_t)count;i++)
      lut_latitude[latitude[i]] = i;
    
    if(verbose){
      cout<<"latitude =";
      for(vector<double>::iterator ilat=latitude.begin(); ilat!=latitude.end(); ilat++)
        cout<<*ilat<<" ";
      cout<<endl;
    }
  }
  
  { // time
    long id = ncvarid(ncid, "time");
    ncvarinq(ncid, id, 0, &xtypep, &ndims, dims, &natts);
    
    long start=0, count=spec["time"];
    time.resize(count);
    if(xtypep==NC_INT){
      vector<int> time_tmp(count);
      ncvarget(ncid, id, &start, &count, &(time_tmp[0]));
      for(size_t i=0;i<(size_t)count;i++)
        time[i] = time_tmp[i];
    }else if(xtypep==NC_FLOAT){
      vector<float> time_tmp(count);
      ncvarget(ncid, id, &start, &count, &(time_tmp[0]));
      for(size_t i=0;i<(size_t)count;i++)
        time[i] = time_tmp[i];
    }else if(xtypep==NC_DOUBLE){
      ncvarget(ncid, id, &start, &count, &(time[0]));
    } else if(xtypep==NC_SHORT) {
       double scale_factor;
       nc_get_att_double(ncid, id, "scale_factor", &scale_factor); 
       double add_offset;
       nc_get_att_double(ncid, id, "add_offset", &add_offset);
       vector<short> time_tmp(count);
       ncvarget(ncid, id, &start, &count, &(time_tmp[0]));
       for (size_t i=0;i<(size_t)count; i++)
           time[i] = time_tmp[i]*scale_factor + add_offset;
    }else{
      cerr<<"ERROR: ("<<__FILE__<<"): time has unexpected type.\n";
      exit(-1);
    }
    
    lut_time[time[0]] = 0;
    for(size_t i=1;i<(size_t)count;i++)
      if(time[i]>time[i-1])
        lut_time[time[i]] = i;
    
    if(verbose){
      for(vector<double>::iterator it=time.begin(); it!=time.end(); it++)
        cout<<*it<<" ";
      cout<<endl;
    }
    
    // Attributes to read
    char units[1024];
    memset(units, '\0', 1024);
    nc_get_att_text(ncid, id, "units", units);
    
    if(calendar==NULL)
      calendar = new Calendar;
    if (verbose) {
        cout << "Time in file has units "<<units<<endl;
    }
    data_time_units = string(units);
  }

#else
  cerr<<"ERROR: No fluxes support compiled\n";
  exit(-1);
#endif  
  return 0;
}

int FluxesReader::SetSimulationTimeUnits(std::string units){
  simulation_time_units = units;

  return 0;
}

int FluxesReader::SetTimeSeconds(double seconds){
  if(verbose) {
    cout<<"int FluxesReader::SetTimeSeconds("<<seconds<<")\n";
    cout<<"will convert from: "<<simulation_time_units<<" to "<<data_time_units<<endl;
  }

  // Convert from simulation time units to data time units
  int err = calendar->SetTransformation(simulation_time_units, data_time_units, "gregorian");
  if(err){
    cerr<<"ERROR ("<<__FILE__<<"): Setting the time transformation has failed.\n";
    exit(-1);
  }
  
  err = calendar->Convert(seconds, time_set);
  if(err){
    cerr<<"ERROR ("<<__FILE__<<"): Conversion between time units has failed.\n";
    exit(-1);
  }

  if (verbose)
    cout<<"New time for data="<<time_set<<"\n";

  modified = true;
  
  Update();
  
  return (0);
}

double FluxesReader::SolveLine(double x0, double x1, double y0, double y1, double y) const{
  //cout<<x0<<" , "<<x1<<" , "<<y0<<" , "<<y1<<" , "<<y<<endl;
  return x0 + (x1-x0)*(y-y0)/(y1-y0);
}

int FluxesReader::Update(){
  if(!modified)
    return 0;
  
  // Make sure that a file has been registered
  if(ERA_data_files.size()==0)
    return 0;
  
  // Get an iterator to the first element which has a value greater than or equal to key
  pair<size_t, size_t> time_interval = GetInterval(time_set, lut_time);
  size_t t0=time_interval.first;
  size_t t1=time_interval.second;
  
  bool need_second_frame = (t0!=t1);
  
  // Deleted data which is no longer needed
  if(!fields.empty())
    if((fields.begin()->first != (int)t0)&&
       (fields.begin()->first != (int)t1))
      fields.erase(fields.begin());
  
  if(!fields.empty())
    if((fields.rbegin()->first != (int)t0)&&
       (fields.rbegin()->first != (int)t1))
      fields.erase(fields.rbegin()->first);
  
  // Read in the new data
  Read(t0);
  Read(t1);
  
  modified = false;
  
  return 0;
}

void FluxesReader::VerboseOff(){
  verbose = false;
  return;
}

void FluxesReader::VerboseOn(){
  verbose = true;
  return;
}

//
// FORTRAN interface
//

FluxesReader FluxesReader_global;

extern "C" {
  void fluxes_addfieldofinterest_c(const char *scalar){
    FluxesReader_global.AddFieldOfInterest(string(scalar));
    return;
  }
  
  void fluxes_clearfields_c(){
    FluxesReader_global.ClearFields();
    return;
  }
  
  void fluxes_getscalars_c(const double *longitude, const double *latitude, double *scalars){
    FluxesReader_global.GetScalars(*longitude, *latitude, scalars);
    return;
  }

  void fluxes_getscalar_c(const char *name, const double *longitude, const double *latitude, double *scalar){
    *scalar = FluxesReader_global.GetScalar(string(name), *longitude, *latitude);
    return;
  }
  
  void fluxes_registerdatafile_c(const char *filename){
    FluxesReader_global.RegisterDataFile(string(filename));
    return;
  }

  void fluxes_setsimulationtimeunits_c(const char *units){
    FluxesReader_global.SetSimulationTimeUnits(string(units));  
    return;
  }
  
  void fluxes_settimeseconds_c(const double *time){
    FluxesReader_global.SetTimeSeconds(*time);  
    return;
  }
}
