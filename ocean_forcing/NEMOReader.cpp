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
#include "NEMOReader.h"
#include <limits>
#include <string.h>

#ifdef __SUN
double round(double x)
{
  return floor(x+0.5); // floor: Returns the largest integral value that is not greater than x (i.e. floor(2.3)=2).
}
#endif

using namespace std;

NEMOReader::NEMOReader(){
  verbose = false;
  MyRank = 0;
  NProcs = 1;
#ifdef HAVE_MPI
  if(MPI::Is_initialized()){
    MyRank = MPI::COMM_WORLD.Get_rank();
    NProcs = MPI::COMM_WORLD.Get_size();
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

NEMOReader::~NEMOReader(){
  if(calendar!=NULL)
    delete calendar;
}

// Add field to the list of fields to be interpolated - the order that
// they are added is the order in which they will be returned on
// interpolation.
void NEMOReader::AddFieldOfInterest(std::string field_name){
  if(verbose)
    cout<<"void NEMOReader::AddFieldOfInterest("<<field_name<<")\n";
  
  for(deque<string>::const_iterator ifield=fields_of_interest.begin(); ifield!=fields_of_interest.end(); ifield++)
    if(field_name==*ifield)
      cerr<<"ERROR: void NEMOReader::AddFieldOfInterest("<<field_name<<")\n Field added multiple times\n";
  
  fields_of_interest.push_back(field_name);
  
  modified = true;
  
  return;
}

void NEMOReader::ClearFields(){
  if(verbose)
    cout<<"void NEMOReader::ClearFields()\n";
  
  fields_of_interest.clear();
  return;
}

bool NEMOReader::Enabled() const{
  if(verbose)
    cout<<"void NEMOReader::Enabled() const\n";
  
  return NEMO_data_files.size()!=0;
}

pair<size_t, size_t> NEMOReader::GetInterval(double value, const map<double, int>&lut) const{
  map<double, int>::const_iterator lower_bound = lut.lower_bound(value);
//   assert(lower_bound!=lut.end());

  map<double, int>::const_iterator upper_bound = lower_bound;
  if(lower_bound!=lut.begin())
    lower_bound--;

  return pair<size_t, size_t>(lower_bound->second, upper_bound->second);
}

int NEMOReader::GetScalars(double xlong, double ylat, double p_depth, double *scalars){
  if(verbose)
    cout<<"int NEMOReader::GetScalars("<<xlong<<", "<<ylat<<", "<<p_depth<<", scalars)\n";
  
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
  if(verbose)
    cout<<"longitude lower bound "<<xlong<<" --> "<<long1->first<<endl;

  map<double, int>::const_iterator long0 = long1;
  if((long0!=lut_longitude.begin())&&(fabs(long0->first-xlong)>0.001))
    long0--;
  size_t i0=long0->second;
  size_t i1=long1->second;
  
  map<double, int>::const_iterator lat1 = lut_latitude.lower_bound(ylat);
  assert(lat1!=lut_latitude.end());
  if(verbose)
    cout<<"latitude lower bound "<<ylat<<" --> "<<lat1->first<<endl;

  map<double, int>::const_iterator lat0 = lat1;
  if((lat0!=lut_latitude.begin())&&(fabs(lat0->first-ylat)>0.001))
    lat0--;
  
  size_t j0=lat0->second;
  size_t j1=lat1->second;

  map<double, int>::const_iterator depth1 = lut_depth.lower_bound(p_depth);
  assert(depth1!=lut_depth.end());
  if(verbose)
    cout<<"depth lower bound "<<p_depth<<" --> "<<depth1->first<<endl;

  map<double, int>::const_iterator depth0 = depth1;
  if((depth0!=lut_depth.begin())&&(fabs(depth0->first-p_depth)>0.001))
    depth0--;
  size_t k0=depth0->second;
  size_t k1=depth1->second;
  
  if(verbose)
    cout<<"Selecting grid ("<<i0<<", "<<j0<<", "<<k0<<"), ("<<i1<<", "<<j1<<", "<<k1<<")\n";


  assert(fields.begin()!=fields.end());
  size_t t0 = fields.begin()->first;
  size_t t1 = fields.rbegin()->first;

  double values[1024]; // This is way more than the number of
                       // variables calculated in NEMO
  size_t nvalues=0, stride=longitude.size(), vstride=longitude.size()*latitude.size();
  
  for(map<int, map<string, vector<double> > >::const_iterator itime=fields.begin(); itime!=fields.end(); itime++){
    for(deque<string>::const_iterator ifield=fields_of_interest.begin(); ifield!=fields_of_interest.end(); ifield++){
      const vector<double> &fld=itime->second.find(*ifield)->second;
      
      string dtest="ssh";
      if(ifield->c_str()==dtest){ // Use ssh to set prssure
        if(i0==i1){ // No interpolation along longitude
          if(verbose)
            cout<<"Case 1.1\n";
          if(j0==j1){
            double x = fld[j0*stride+i0];
            values[nvalues] = x;
          }else{
            double x0 = fld[j0*stride+i0];
            double x1 = fld[j1*stride+i0];
          
            double x = SolveLine(x0, x1, ylat-lat0->first, lat1->first, ylat);
            values[nvalues] = x;
          }
        }else if(j0==j1){ // No interpolation along latitude
          if(verbose)
            cout<<"Case 1.2\n";
          double x0 = fld[j0*stride+i0];
          double x1 = fld[j0*stride+i1];

          double long_temp0=long0->first;
          double long_temp1=long1->first;

          if (long_temp1 < long_temp0) {
            long_temp1 += 360.0;
          }
        
          double x = SolveLine(x0, x1, long_temp0, long_temp1, xlong);
          values[nvalues] = x;
        }else{ // Bi-linear interpolation
          if(verbose) {
            cout<<"Case 1.3 -- "<<fields.size()<<endl;
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
      }else if(k0==k1){ // No interpolation in the radial direction
        if(i0==i1){ // No interpolation along longitude
          if(verbose)
            cout<<"Case 2.1\n";
          if(j0==j1){
            double x = fld[k0*vstride+j0*stride+i0];
            values[nvalues] = x;
          }else{
            double x0 = fld[k0*vstride+j0*stride+i0];
            double x1 = fld[k0*vstride+j1*stride+i0];
          
            double x = SolveLine(x0, x1, ylat-lat0->first, lat1->first, ylat);
            values[nvalues] = x;
          }
        }else if(j0==j1){ // No interpolation along latitude
          if(verbose)
            cout<<"Case 2.2\n";
          double x0 = fld[k0*vstride+j0*stride+i0];
          double x1 = fld[k0*vstride+j0*stride+i1];

          double long_temp0=long0->first;
          double long_temp1=long1->first;

          if (long_temp1 < long_temp0) {
            long_temp1 += 360.0;
          }
        
          double x = SolveLine(x0, x1, long_temp0, long_temp1, xlong);
          values[nvalues] = x;
        }else{ // Bi-linear interpolation
          if(verbose) {
            cout<<"Case 2.3 -- "<<fields.size()<<endl;
          }
          double x00 = fld[k0*vstride+j0*stride+i0];
          double x10 = fld[k0*vstride+j0*stride+i1];
          double x01 = fld[k0*vstride+j1*stride+i0];
          double x11 = fld[k0*vstride+j1*stride+i1];

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
      }else{
        if(i0==i1){ // No interpolation along longitude
          if(verbose)
            cout<<"Case 3.1\n";
          if(j0==j1){
            double x1 = fld[k0*vstride+j0*stride+i0];
            double x2 = fld[k1*vstride+j0*stride+i0];
            double ddepth=depth[k1]-depth[k0];
            double pdepth=p_depth-depth[k0];
            double x=x2*(pdepth/ddepth)+x1*(ddepth-pdepth)/ddepth;
            values[nvalues] = x;
          }else{
            double x0 = fld[k0*vstride+j0*stride+i0];
            double x1 = fld[k0*vstride+j1*stride+i0];
            double x2 = fld[k1*vstride+j0*stride+i0];
            double x3 = fld[k1*vstride+j1*stride+i0];
          
            double xk0 = SolveLine(x0, x1, ylat-lat0->first, lat1->first, ylat);
            double xk1 = SolveLine(x2, x3, ylat-lat0->first, lat1->first, ylat);
            double ddepth=depth[k1]-depth[k0];
            double pdepth=p_depth-depth[k0];
            double x=xk1*(pdepth/ddepth)+xk0*(ddepth-pdepth)/ddepth;
            values[nvalues] = x;
          }
        }else if(j0==j1){ // No interpolation along latitude
          if(verbose)
            cout<<"Case 3.2\n";
          double x0 = fld[k0*vstride+j0*stride+i0];
          double x1 = fld[k0*vstride+j0*stride+i1];
          double x2 = fld[k1*vstride+j0*stride+i0];
          double x3 = fld[k1*vstride+j0*stride+i1];

          double long_temp0=long0->first;
          double long_temp1=long1->first;

          if (long_temp1 < long_temp0) {
            long_temp1 += 360.0;
          }
        
          double xk0 = SolveLine(x0, x1, long_temp0, long_temp1, xlong);
          double xk1 = SolveLine(x2, x3, long_temp0, long_temp1, xlong);
          double ddepth=depth[k1]-depth[k0];
          double pdepth=p_depth-depth[k0];
          double x=xk1*(pdepth/ddepth)+xk0*(ddepth-pdepth)/ddepth;
          values[nvalues] = x;
        }else{ // Tri-linear interpolation
          if(verbose) {
            cout<<"Case 3.3 -- "<<fields.size()<<endl;
          }
          double x000 = fld[k0*vstride+j0*stride+i0];
          double x100 = fld[k0*vstride+j0*stride+i1];
          double x010 = fld[k0*vstride+j1*stride+i0];
          double x110 = fld[k0*vstride+j1*stride+i1];

          double x001 = fld[k1*vstride+j0*stride+i0];
          double x101 = fld[k1*vstride+j0*stride+i1];
          double x011 = fld[k1*vstride+j1*stride+i0];
          double x111 = fld[k1*vstride+j1*stride+i1];

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

          double x1 =  (x000*(dl1)*(ylat-latitude[j1]) -
                        x100*(dl0)*(ylat-latitude[j1]) -
                        x010*(dl1)*(ylat-latitude[j0]) + 
                        x110*(dl0)*(ylat-latitude[j0]))
                        /(dlong*dlat);
          double x2 =  (x001*(dl1)*(ylat-latitude[j1]) -
                        x101*(dl0)*(ylat-latitude[j1]) -
                        x011*(dl1)*(ylat-latitude[j0]) + 
                        x111*(dl0)*(ylat-latitude[j0]))
                        /(dlong*dlat);
          double ddepth=depth[k1]-depth[k0];
          double pdepth=p_depth-depth[k0];
          double x=x2*(pdepth/ddepth)+x1*(ddepth-pdepth)/ddepth;
          values[nvalues] = x;
        }
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

double NEMOReader::GetScaleFactor(std::string variable_name) const{
  return 1.0;  
}

int NEMOReader::Read(int time_index){
  if(verbose)
    cout<<"int Read("<<time_index<<")\n";
  
  // Make sure that a file has been registered
  if(NEMO_data_files.size()==0){
    cerr<<"ERROR: int Read("<<time_index<<") -- no file has been registered\n";
    exit(-1);
  }
  
#ifdef HAVE_NETCDF
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
//    assert(xtypep==NC_SHORT);

    size_t len;
    vector<float> field;
   
    string dtest="ssh";
    if (ifield->c_str()==dtest){
      len = vdimension[1]*vdimension[0];
      field.resize(len);
      long start[]={time_index, 0, 0}, count[]={1,vdimension[1],vdimension[0]};
      ncvarget(ncid, id, start, count, &(field[0]));
    }else{
      len = ndepth*vdimension[1]*vdimension[0];
      field.resize(len);
      long start[]={time_index, 0, 0, 0}, count[]={1,ndepth,vdimension[1],vdimension[0]};
      ncvarget(ncid, id, start, count, &(field[0]));
    }
    
//     string extension=".dat";
//     string fname=ifield->c_str();
//     string dfile=fname+extension;
//     ofstream myfile(dfile.c_str());
// 
//     for(size_t i=0;i<(size_t)len;i++){
//       myfile << i << '\t' << field[i] << endl;
//     }

    fields[time_index][*ifield].resize(len);
    for(size_t i=0; i<len; i++){
      fields[time_index][*ifield][i] = field[i];
    }

  }
#else
  cerr<<"ERROR: No NetCDF support compiled\n";
  exit(-1);
#endif
  return 0;
}

// Returns 0 on success, negitive if error.
int NEMOReader::RegisterDataFile(string file){
#ifdef HAVE_NETCDF
  if(verbose)
    cout<<"int NEMOReader::RegisterDataFile("<<file<<")\n";
  
  // Check what has gone before
  if(NEMO_data_files.size()==1){
    assert(NEMO_data_files[0]==file);
    return 0;
  }
  // else
  
  NEMO_data_files.push_back(file);
  
  //
  // Get the basic information from the NetCDF file
  //
  
  // Make NetCDF errors verbose and fatal
  ncopts = NC_VERBOSE | NC_FATAL;
  
  // Open the netCDF file.
  if(file.c_str()[0]!='/')
    file = string(getenv("PWD"))+string("/")+file;
  ncid = ncopen(file.c_str(), NC_NOWRITE);

  if(ncid==-1){
    cerr<<"ERROR: could not open netcdf file: "<<file<<endl;
    FLAbort("int NEMOReader::RegisterDataFile(string file)", __FILE__, __LINE__);
  }

  // If this is the first time we're opened one of these files then we
  // need to store it's basic specification (dimensions
  // etc). Otherwise we read in the new specification and reassure
  // ourselves that it's the same as the original.

  { // Get dimensions -- longitude
    int id = ncdimid(ncid, "longitude");
    if(ncerr!=NC_NOERR){
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<": ERROR - dimension variable \"longitude\" does not exist\n";
      exit(-1);
    }
  ncdiminq(ncid, id, (char *)0, &(vdimension[0]));
  
  if(verbose)
      cout<<"("<<__FILE__<<", "<<__LINE__<<"): longitude = "<<vdimension[0]<<endl;
  }

  { // Get dimensions -- latitude
    int id = ncdimid(ncid, "latitude");
    if(ncerr!=NC_NOERR){
      cout.flush();
      cerr<<__FILE__<<", "<<__LINE__<<": ERROR - dimension variable \"latitude\" does not exist\n";
      exit(-1);
    }
    ncdiminq(ncid, id, (char *)0, &(vdimension[1]));

    if(verbose)
      cout<<"("<<__FILE__<<", "<<__LINE__<<"): latitude = "<<vdimension[1]<<endl;
  }

  { // Get the number of layers
    int id = ncdimid(ncid, "depth_level");
    ncdiminq(ncid, id, (char *)0, &ndepth);

    assert(spec.find("depth_level")==spec.end());
    spec["depth_level"] = ndepth;
    
    if(verbose)
      cout<<"("<<__FILE__<<", "<<__LINE__<<"): depth_level = "<<ndepth<<endl;
  }

  { // time ('time_counter' in NEMO datafiles)
    int id = ncdimid(ncid, "time_counter");
    ncdiminq(ncid, id, (char *)0, &ntime);

    assert(spec.find("time_counter")==spec.end());
    spec["time_counter"] = ntime;
    
    if(verbose)
      cout<<"("<<__FILE__<<", "<<__LINE__<<"): time_counter = "<<ntime<<endl;
  }

  {// Load longitude values
    nc_type xtypep;                 /* variable type */
    int ndims;                      /* number of dims */
    int dims[MAX_VAR_DIMS];         /* variable shape */
  
    int id = ncvarid(ncid, "longitude");
    ncvarinq(ncid, id, 0, &xtypep, &ndims, dims, NULL);

    long start[]={0}, count[1];
    count[0] = vdimension[0];

    long len = vdimension[0];
    longitude.resize(len);
    if(xtypep==NC_FLOAT){ // single precision floating point number
      vector<float> fvar(len);
      ncvarget(ncid, id, start, count, &(fvar[0]));
      for(long i=0;i<len;i++)
        longitude[i] = fvar[i];
      fvar.clear();
    }else if(xtypep==NC_DOUBLE){ // double precision floating point number
      ncvarget(ncid, id, start, count, &(longitude[0]));
    }else{
      cerr<<"ERROR: unknown data type\n";
      exit (-1);
    }

    // Currently, can calculate dlong as this:
    dlong=longitude[1]-longitude[0];

    // Set Longitude to be between 0 and 360
    for(long i=0;i<len;i++){
      if(longitude[i]<0.0){
        longitude[i]+=360.0;
      }
    }

    for(long i=0;i<len;i++)
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
  

  {// Load latitudes values
    nc_type xtypep;                 /* variable type */
    int ndims;                      /* number of dims */
    int dims[MAX_VAR_DIMS];         /* variable shape */
  
    int id = ncvarid(ncid, "latitude");
    ncvarinq(ncid, id, 0, &xtypep, &ndims, dims, NULL);

    long start[]={0}, count[1];
    count[0] = vdimension[1];

    long len = vdimension[1];
    latitude.resize(len);
    if(xtypep==NC_FLOAT){ // single precision floating point number
      vector<float> fvar(len);
      ncvarget(ncid, id, start, count, &(fvar[0]));
      for(long i=0;i<len;i++)
        latitude[i] = fvar[i];
      fvar.clear();
    }else if(xtypep==NC_DOUBLE){ // double precision floating point number
      ncvarget(ncid, id, start, count, &(latitude[0]));
    }else{
      cerr<<"ERROR: unknown data type\n";
      exit (-1);
    }

    // Currently, can calculate dat as this:
    dlat=latitude[1]-latitude[0];

    for(long i=0;i<len;i++)
      lut_latitude[latitude[i]] = i;
    
    if(verbose){
      cout<<"latitude =";
      for(vector<double>::iterator ilat=latitude.begin(); ilat!=latitude.end(); ilat++)
        cout<<*ilat<<" ";
      cout<<endl;
    }
  
  }

//   {
//      ofstream myfile("longlat.dat");
// 
//      long leni = vdimension[0];
//      long lenj = vdimension[1];
//      for(long i=0;i<leni;i++){
//        for(long j=0;j<lenj;j++){
//          myfile << i << '\t' << longitude[i] << '\t' << j << '\t' << latitude[j] << endl;
//        }
//      }
//   }

  { // load depth levels
    nc_type xtypep;                 /* variable type */
    int ndims;                      /* number of dims */
    int dims[MAX_VAR_DIMS];         /* variable shape */
  
    int id = ncvarid(ncid, "depth_level");
    ncvarinq(ncid, id, 0, &xtypep, &ndims, dims, NULL);

    long start[]={0}, count[1];
    count[0] = ndepth;

    long len = ndepth;
    depth.resize(len);
    if(xtypep==NC_FLOAT){ // single precision floating point number
      vector<float> fvar(len);
      ncvarget(ncid, id, start, count, &(fvar[0]));
      for(long i=0;i<len;i++)
        depth[i] = fvar[i];
      fvar.clear();
    }else if(xtypep==NC_DOUBLE){ // double precision floating point number
      ncvarget(ncid, id, start, count, &(depth[0]));
    }else{
      cerr<<"ERROR: unknown data type\n";
      exit (-1);
    }

    for(long i=0;i<len;i++)
      lut_depth[depth[i]] = i;
    
    if(verbose){
      cout<<"depth =";
      for(vector<double>::iterator idepth=depth.begin(); idepth!=depth.end(); idepth++)
        cout<<*idepth<<" ";
      cout<<endl;
    }
  
  }
  
  { // time
    nc_type xtypep;                 /* variable type */
    int ndims, natts;               /* number of dims  and attributes */
    int dims[MAX_VAR_DIMS];         /* variable shape */

    long id = ncvarid(ncid, "time_counter");
    ncvarinq(ncid, id, 0, &xtypep, &ndims, dims, &natts);
    
    long start=0, count=spec["time_counter"];
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
  cerr<<"ERROR: No NetCDF support compiled\n";
  exit(-1);
#endif  
  return 0;
}

int NEMOReader::SetSimulationTimeUnits(std::string units){
  simulation_time_units = units;

  return 0;
}

int NEMOReader::SetTimeSeconds(double seconds){
  if(verbose) {
    cout<<"int NEMOReader::SetTimeSeconds("<<seconds<<")\n";
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

double NEMOReader::SolveLine(double x0, double x1, double y0, double y1, double y) const{
  return x0 + (x1-x0)*(y-y0)/(y1-y0);
}

int NEMOReader::Update(){
  if(!modified)
    return 0;
  
  // Make sure that a file has been registered
  if(NEMO_data_files.size()==0)
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

void NEMOReader::VerboseOff(){
  verbose = false;
  return;
}

void NEMOReader::VerboseOn(){
  verbose = true;
  return;
}

//
// FORTRAN interface
//

NEMOReader NEMOReader_global;

extern "C" {
#define nemo_addfieldofinterest_fc F77_FUNC_(nemo_addfieldofinterest, NEMO_ADDFIELDOFINTEREST)
  void nemo_addfieldofinterest_fc(char *_scalar, int *len){
    char scalar[1024];
    assert(*len<1023);
    strncpy(scalar, _scalar, *len);
    scalar[*len] = '\0';
    NEMOReader_global.AddFieldOfInterest(string(scalar));
    return;
  }
  
#define nemo_clearfields_fc F77_FUNC_(nemo_clearfields, NEMO_CLEARFIELDS)
  void nemo_clearfields_fc(){
    NEMOReader_global.ClearFields();
    return;
  }
  
#define nemo_getscalars_fc F77_FUNC_(nemo_getscalars, NEMO_GETSCALARS)
  void nemo_getscalars_fc(double *longitude, double *latitude, double *p_depth, double *scalars){
    NEMOReader_global.GetScalars(*longitude, *latitude, *p_depth, scalars);
    return;
  }
  
#define nemo_registerdatafile_fc F77_FUNC_(nemo_registerdatafile, NEMO_REGISTERDATAFILE)
  void nemo_registerdatafile_fc(char *_filename, int *len){
    char filename[4096];
    assert(*len<4095);
    strncpy(filename, _filename, *len);
    filename[*len] = '\0';
    NEMOReader_global.RegisterDataFile(string(filename));
    return;
  }
  
#define nemo_settimeseconds_fc F77_FUNC_(nemo_settimeseconds, NEMO_SETTIMESECONDS)
  void nemo_settimeseconds_fc(double *_time){
    NEMOReader_global.SetTimeSeconds(*_time);  
    return;
  }
}