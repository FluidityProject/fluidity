/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    C.Pain@imperial.ac.uk
    
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

    This creates a climatology atlas, for use with ICOM, using "High
    resolution (1/4 degree) Temperature and Salinity Analyses of the
    World's Oceans. Version 2"

    The 1/4 degree grid climatological mean fields of in situ
    temperature (degrees Celsius) and salinity (practical salinity
    scale) for the annual, seasonal, and monthly time periods were
    calculated by Boyer2005 using objective analysis techniques. The
    data and associated metadata was obtained from the NODC (
    http://www.nodc.noaa.gov/OC5/WOA01/qd_ts01.html). All the data
    files are gzipped ASCII 1/4 gridded data files and contain
    the DOS end-of-line character (^M). There are 12 monthly averages
    of temperature (t001_v2.gz, t002_v2.gz, t003_v2.gz, t004_v2.gz,
    t005_v2.gz, t006_v2.gz, t007_v2.gz, t008_v2.gz, t009_v2.gz,
    t010_v2.gz, t011_v2.gz, t012_v2.gz) and 12 monthly averages of
    salinity (s001_v2.gz, s002_v2.gz, s003_v2.gz, s004_v2.gz,
    s005_v2.gz, s006_v2.gz, s007_v2.gz, s008_v2.gz, s009_v2.gz,
    s010_v2.gz, s011_v2.gz, s012_v2.gz) corresponding to January,
    February, March, April, May, June, July, August, September,
    October, November and December respectively. In addition, there 4
    seasonal averages of temperature (t013AV_v2.gz, t014AV_v2.gz,
    t053AV_v2.gz, t016AV_v2.gz) and 4 seasonal averages of salinity
    (s013AV_v2.gz, s014AV_v2.gz, s053AV_v2.gz, s016AV_v2.gz)
    corresponding to winter (defined as January, February and March),
    spring summer and autumn.

    The first value is the value for the 1/4 degree square centered at
    (0.125E, -89.875S) representing the area from 0 to 0.25E longitude
    and -90. to -89.75 latitude. The indices increment in the
    longitude direction first, so the second value represents the 1/4
    degree square centered at (0.375E, -89.875S). There are 1440
    longitude boxes for each latitude, so the 1441st value represent
    the 1/4 degree grid box centered at (0.125E,-89.625S). There are
    24 depth levels in the monthly fields. The levels correspond to
    the standard depths (0., 10., 20., 30., 50., 75., 100., 125.,
    150., 200., 250., 300., 400., 500., 600., 700., 800., 900., 1000.,
    1100., 1200., 1300., 1400., 1500. meters). For the seasonal
    averages there are an additional 9 standard depth levels (1750.,
    2000., 2500., 3000., 3500., 4000., 4500., 5000.,
    5500. meters). The land is masked using a value of -99.9999.

    For use with ICOM (relaxing to climatology at the boundaries), a
    NetCDF data was created containing the monthly means, and the
    additional 9 standard levels from the seasonal means in order to
    provide information below the 1500m level. In addition, Killworth
    correction (Killworth 1995) is applied to the data in order to
    facilitate accutate time interpolation.

    To use this just execute it in the directory which contains all
    the files listed above gunzipped.
*/
#include "confdefs.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <iostream>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_LIBNETCDF
#include <netcdf.h>
#endif

using namespace std;

const short fillvalue=32767;
const int idim0=1440, jdim0=720;
const float land_val=-99.9999;

extern "C" {
  void dgesv_(int *N, int *NRHS, double A[], int *LDA, int IPIVOT[], double B[], int *LDB, int *INFO);
}

#if defined(__SUN) || defined (__PGI)
long lround(double x){
  if(x<0.0)
    return (long)(-floor(-x+0.5));
  else
    return (long)(floor(x+0.5));
}
#endif

int killworth_matrix_seasonal(double *A){
  for(size_t i=0;i<16;i++)
    A[i] = 0.0;
  
  double l[] = {90.25, 91.0, 91.0, 92.0};
  for(size_t n=0;n<4;n++){
    double e = 0.25*l[n]/(l[(n+4-1)%4]+l[n]);
    double g = 0.25*l[n]/(l[n]+l[(n+1)%4]);
    double f = 1.0 - e - g;
    
    A[n*4+(n+4-1)%4] = e;
    A[n*4+n]           = f;
    A[n*4+(n+1)%4]    = g;
  }
  
  return 0;
}

int killworth_matrix_monthly(double *A){
  for(size_t i=0;i<144;i++)
    A[i] = 0.0;
  
  double l[] = {31.0, 28.25, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0};
  for(size_t n=0;n<12;n++){
    double e = 0.25*l[n]/(l[(n+12-1)%12]+l[n]);
    double g = 0.25*l[n]/(l[n]+l[(n+1)%12]);
    double f = 1.0 - e - g;
    
    A[n*12+(n+12-1)%12] = e;
    A[n*12+n]           = f;
    A[n*12+(n+1)%12]    = g;
  }
  
  return 0;
}

int kilworth_correction(vector<double> &A, vector<double> &d){
  int n=d.size();
  int nrhs=1,ldb=n;
  vector<int> ipivot(n);
  int info;
  vector<double> AA(A);
  dgesv_(&n, &nrhs, &(AA[0]), &n, &(ipivot[0]), &(d[0]), &ldb, &info);
  return 0;
}

int killworth_correction(deque< deque< vector<float> > >& variable, float &max_v, float &min_v){
  const int n=variable.size();
  const size_t nlevels=variable[0].size();
  vector<double> A(n*n, 0.0);
  if(n==12){
    killworth_matrix_monthly(&(A[0]));
  }else if(n==4){
    killworth_matrix_seasonal(&(A[0]));
  }else{
    cerr<<"ERROR: time not reckonised\n";
    exit(-1);
  }
  
  //cout<<"again A="<<endl;
  //for(size_t ii=0;ii<n;ii++){
  //  for(size_t jj=0;jj<n;jj++)
  //    cout<<A[ii*n+jj]<<" ";
  //  cout<<endl;
  //}
 
  max_v=-1000.0;
  min_v=1000.0;
  
  vector<double> d(n);
  for(size_t k=0;k<nlevels;k++){
    for(size_t j=0; j<(size_t)jdim0; j++){
      for(size_t i=0; i<(size_t)idim0; i++){
    if(fabs(variable[0][k][j*idim0+i]-land_val)>1.0){
      for(size_t m=0;m<variable.size();m++)
        d[m] = variable[m][k][j*idim0+i];
      
      kilworth_correction(A, d);
      
      for(size_t m=0;m<variable.size();m++){
      if(d[m]>1000.0){
        cout<<"A="<<endl;
        for(size_t ii=0;ii<(size_t)n;ii++){
        for(size_t jj=0;jj<(size_t)n;jj++)
          cout<<A[ii*n+jj]<<" ";
        cout<<endl;
        }
        cout<<"d = ";
        for(size_t m=0;m<variable.size();m++)
        cout<<variable[m][k][j*idim0+i]<<" ";
        cout<<endl;
        cout<<"d' = ";
        for(size_t m=0;m<variable.size();m++)
        cout<<d[m]<<" ";
        cout<<endl;
      }
      max_v = max(max_v, (float)d[m]);
      min_v = min(min_v, (float)d[m]);
      variable[m][k][j*idim0+i] = d[m];
      }
    }
      }
    }
  }

 return 0;
}

int diffuse_boundaries(deque< deque< vector<float> > > &dat, float minv, int iterations){
  // Iterate through time
  for(deque< deque< vector<float> > >::iterator tdat=dat.begin(); tdat!=dat.end(); tdat++){
  // Iterate through levels
    for(size_t k=0;k<tdat->size();k++){
    
    // Find which points are null
    vector<bool> null_points(idim0*jdim0, false);
    for(size_t j=0; j<(size_t)jdim0; j++){
    for(size_t i=0; i<(size_t)idim0; i++){
          if(fabs((*tdat)[k][j*idim0+i]-land_val)<1.0){
      null_points[j*idim0+i] = true;
      (*tdat)[k][j*idim0+i] = minv;
          }
    }
      }
    
    // Start diffusing
    for(size_t its=0;its<iterations;its++){
    for(size_t j=0; j<(size_t)jdim0; j++){ 
       for(size_t i=0; i<(size_t)idim0; i++){
      if(null_points[j*idim0+i]){
        float sum=(*tdat)[k][j*idim0+i];
        int npts=1;
        
        if(j>0){
        sum+=(*tdat)[k][(j-1)*idim0+i];
        npts++;
        }

         if(j+1<jdim0){
        sum+=(*tdat)[k][(j+1)*idim0+i];
        npts++;
        }

        sum+=(*tdat)[k][j*idim0+(i+idim0-1)%idim0];
        npts++;

        sum+=(*tdat)[k][j*idim0+(i+1)%idim0];
        npts++;
        
        if(k>0){
        sum+=(*tdat)[k-1][j*idim0+i];
        npts++;
        }

        (*tdat)[k][j*idim0+i] = sum/npts;
      }
      }
    }
    }
  }
  }

  return 0;
}

int write_netcdf_variable(int ncid, const int* dimids, const char* name, const char* long_name, const char *units,
        deque< deque< vector<float> > >& variable){
#ifdef HAVE_LIBNETCDF

  float max_v, min_v;
  killworth_correction(variable, max_v, min_v);

  cout<<"Diffusing values into interior\n";
  diffuse_boundaries(variable, min_v, 1000);
  
  cout<<"Writing out "<<name<<" with vals in the range "<<min_v<<" -> "<<max_v<<endl;

  // Put file into define mode
  ncredef(ncid);
  
  int varid = ncvardef(ncid, name, NC_SHORT, 4, dimids);
  
  ncattput(ncid, varid, "long_name",     NC_CHAR,  strlen(long_name), long_name);
  ncattput(ncid, varid, "units",         NC_CHAR,  strlen(units),     units);
  ncattput(ncid, varid, "_FillValue",    NC_SHORT, 1,                 &fillvalue);
  ncattput(ncid, varid, "missing_value", NC_SHORT, 1,                 &fillvalue);

  // Packing algorithm
  float add_offset = (max_v + min_v)*0.5;
  ncattput(ncid, varid, "add_offset", NC_FLOAT, 1, &add_offset);
  
  // 2^16
  float scale_factor = (max_v - min_v)/(65536 - 5);
  ncattput(ncid, varid, "scale_factor", NC_FLOAT, 1, &scale_factor);
  
  // Finish defining metadata for salinity
  ncendef(ncid);
  
  // Write out variable
  varid = ncvarid(ncid, name);
  
  long start[]={0,0,0,0}, count[]={1,1,jdim0,idim0};
  vector<short> frame(idim0*jdim0);
  for(size_t m=0;m<variable.size();m++){
    start[0] = m;
    for(size_t k=0;k<variable[m].size();k++){
      start[1] = k;
      for(size_t j=0; j<(size_t)jdim0; j++)
    for(size_t i=0; i<(size_t)idim0; i++){
      // After applying a diffusion filter we need to ensure the
      // values are still within bounds.
      float v = min(max_v, variable[m][k][j*idim0+i]);
      v = max(v, min_v);
      
      frame[j*idim0+i] = lround((v - add_offset)/scale_factor);
    }
      ncvarput(ncid, varid, start, count, &(frame[0]));
    }
  }
#else
  cerr<<"ERROR: No NetCDF support compiled.\n";
  exit(-1);
#endif
  return 0;
}

int main(){
#ifdef HAVE_LIBNETCDF
  // Start to write the climatology data
  int ncid = nccreate("climatology.nc", NC_CLOBBER);
  
  //
  // Define the dimensions
  //
  int dimids[4];
  int time0 = ncdimdef(ncid, "time0",     12);
  int time1 = ncdimdef(ncid, "time1",     4);

  int levels0 = ncdimdef(ncid, "levels0",   24);
  int varid = ncvardef(ncid, "levels0", NC_FLOAT, 1, &levels0);
  ncattput(ncid, varid, "long_name", NC_CHAR, 19, "Top standard levels");
  ncattput(ncid, varid, "units", NC_CHAR, 6, "Meters");

  int levels1 = ncdimdef(ncid, "levels1",   9);
  varid = ncvardef(ncid, "levels1", NC_FLOAT, 1, &levels1);
  ncattput(ncid, varid, "long_name", NC_CHAR, 22, "Bottom standard levels");
  ncattput(ncid, varid, "units", NC_CHAR, 6, "Meters");
  
  int latitude = ncdimdef(ncid, "latitude",  jdim0);
  varid = ncvardef(ncid, "latitude", NC_FLOAT, 1, &latitude);
  ncattput(ncid, varid, "long_name", NC_CHAR, 16, "Degrees latitude");
  ncattput(ncid, varid, "units", NC_CHAR, 7, "Degrees");

  int longitude = ncdimdef(ncid, "longitude", idim0);
  varid = ncvardef(ncid, "longitude", NC_FLOAT, 1, &longitude);
  ncattput(ncid, varid, "long_name", NC_CHAR, 17, "Degrees longitude");
  ncattput(ncid, varid, "units", NC_CHAR, 7, "Degrees");

  // Finish writing mata data for now
  ncendef(ncid);

  // Write out dimension markers

  vector<float> markers;
  long start=0, count;
  varid = ncvarid(ncid, "levels0");
  markers.push_back(0.);
  markers.push_back(10.);
  markers.push_back(20.);
  markers.push_back(30.);
  markers.push_back(50.); 
  markers.push_back(75.);
  markers.push_back(100.); 
  markers.push_back(125.); 
  markers.push_back(150.); 
  markers.push_back(200.); 
  markers.push_back(250.); 
  markers.push_back(300.); 
  markers.push_back(400.); 
  markers.push_back(500.);
  markers.push_back(600.); 
  markers.push_back(700.); 
  markers.push_back(800.); 
  markers.push_back(900.); 
  markers.push_back(1000.); 
  markers.push_back(1100.); 
  markers.push_back(1200.); 
  markers.push_back(1300.); 
  markers.push_back(1400.);
  markers.push_back(1500.);
  count = markers.size();
  ncvarput(ncid, varid, &start, &count, &(markers[0]));
  
  varid = ncvarid(ncid, "levels1");
  markers.clear();
  markers.push_back(1750.); 
  markers.push_back(2000.);
  markers.push_back(2500.);
  markers.push_back(3000.);
  markers.push_back(3500.);
  markers.push_back(4000.);
  markers.push_back(4500.);
  markers.push_back(5000.);
  markers.push_back(5500.);
  count = markers.size();
  ncvarput(ncid, varid, &start, &count, &(markers[0]));
  
  varid = ncvarid(ncid, "latitude");
  markers.clear();
  for(int i=0;i<jdim0;i++)
    markers.push_back(-90.0 + i*0.25);
  count = markers.size();
  ncvarput(ncid, varid, &start, &count, &(markers[0]));
  
  varid = ncvarid(ncid, "longitude");
  markers.clear();
  for(int i=0;i<idim0;i++)
    markers.push_back(i*0.25);
  count = markers.size();
  ncvarput(ncid, varid, &start, &count, &(markers[0]));
      
  // Set up 12 monthly means
  deque< deque< vector<float> > > dat(12);
  
  // Allocate 24 levels for each month
  for(size_t i=0;i<12;i++){
    dat[i].resize(24);
    // Allocate 1440x720 frame for each level
    for(size_t j=0;j<24;j++)
      dat[i][j].resize(idim0*jdim0);
  }
  
  // Monthly salinity
  const char *salinity_monthly[] = {"s001_v2","s002_v2","s003_v2","s004_v2","s005_v2","s006_v2",
                                    "s007_v2","s008_v2","s009_v2","s010_v2","s011_v2","s012_v2"};
  
  for(int m=0;m<12;m++){
    FILE *fp;
    
    if ((fp = fopen(salinity_monthly[m],"rb+\0")) == NULL)
      cerr<<"ERROR: UNABLE TO OPEN FILE\n";
    
    cout<<"Reading file "<<salinity_monthly[m]<<endl;
    
    // Read data
    for(size_t k=0; k<24; k++){
      // cout<<"Reading level "<<k<<endl;
      for(size_t j=0; j<(size_t)jdim0; j++)
  for(size_t i=0; i<(size_t)idim0; i++){
    int ierr = fscanf(fp, "%f", &(dat[m][k][j*idim0+i]));
    if (ierr<1) cerr<<"File read error";
  }
    }
    fclose(fp);
  }
  
  dimids[0] = time0;
  dimids[1] = levels0;
  dimids[2] = latitude;
  dimids[3] = longitude;
  write_netcdf_variable(ncid, dimids, "salinity_monthly", "Monthly mean salinity", "Practical salinity units", dat);

  // Monthly temperature
  const char *temperature_monthly[] = {"t001_v2","t002_v2","t003_v2","t004_v2","t005_v2","t006_v2",
                                       "t007_v2","t008_v2","t009_v2","t010_v2","t011_v2","t012_v2"};
    
  for(int m=0;m<12;m++){
    FILE *fp;
    
    if ((fp = fopen(temperature_monthly[m],"rb+\0")) == NULL)
      cerr<<"ERROR: UNABLE TO OPEN FILE\n";
    
    cout<<"Reading file "<<temperature_monthly[m]<<endl;
    
    // Read data
    for(size_t k=0; k<24; k++){
      // cout<<"Reading level "<<k<<endl;
      for(size_t j=0; j<(size_t)jdim0; j++)
  for(size_t i=0; i<(size_t)idim0; i++){
    int ierr = fscanf(fp, "%f", &(dat[m][k][j*idim0+i]));
    if (ierr<1) cerr<<"File read error";
  }
    }
    fclose(fp);
  }
  
  dimids[0] = time0;
  dimids[1] = levels0;
  dimids[2] = latitude;
  dimids[3] = longitude;
  write_netcdf_variable(ncid, dimids, "temperature_monthly", "Monthly mean in situ temperature", "Degrees Celsius", dat);
  
  // Seasonal salinity
  dat.resize(4);
  
  // Allocate 9 levels for each season
  for(size_t i=0;i<4;i++){
    dat[i].resize(9);
    // Allocate 1440x720 frame for each level
    for(size_t j=0;j<9;j++)
      dat[i][j].resize(idim0*jdim0);
  }
  const char *salinity_seasonal[] = {"s013AV_v2","s014AV_v2","s015AV_v2","s016AV_v2"};
    
  for(int m=0;m<4;m++){
    FILE *fp;
    
    if ((fp = fopen(salinity_seasonal[m],"rb+\0")) == NULL)
      cerr<<"ERROR: UNABLE TO OPEN FILE\n";
    
    cout<<"Reading file "<<salinity_seasonal[m]<<endl;
    
    // Read data
    for(size_t k=0; k<33; k++){
      // cout<<"Reading level "<<k<<endl;
      for(size_t j=0; j<(size_t)jdim0; j++)
  for(size_t i=0; i<(size_t)idim0; i++){
    float v;
    int ierr = fscanf(fp, "%f", &v);
    if (ierr<1) cerr<<"File read error";
    if(k>=24){
      dat[m][k-24][j*idim0+i] = v;
    }
  }
    }
    fclose(fp);
  }
  
  dimids[0] = time1;
  dimids[1] = levels1;
  dimids[2] = latitude;
  dimids[3] = longitude;
  write_netcdf_variable(ncid, dimids, "salinity_seasonal", "Seasonal mean salinity", "Practical salinity units", dat);

  // Seasonal temperature
  const char *temperature_seasonal[] = {"t013AV_v2","t014AV_v2","t015AV_v2","t016AV_v2"};
    
  for(int m=0;m<4;m++){
    FILE *fp;
    
    if ((fp = fopen(temperature_seasonal[m],"rb+\0")) == NULL)
      cerr<<"ERROR: UNABLE TO OPEN FILE\n";
    
    cout<<"Reading file "<<temperature_seasonal[m]<<endl;
    
    // Read data
    for(size_t k=0; k<33; k++){
      // cout<<"Reading level "<<k<<endl;
      for(size_t j=0; j<(size_t)jdim0; j++)
  for(size_t i=0; i<(size_t)idim0; i++){
    float v;
    int ierr = fscanf(fp, "%f", &v);
    if (ierr<1) cerr<<"File read error";
    if(k>=24){
      dat[m][k-24][j*idim0+i] = v;
    }
  }
    }
    fclose(fp);
  }
  
  dimids[0] = time1;
  dimids[1] = levels1;
  dimids[2] = latitude;
  dimids[3] = longitude;
  write_netcdf_variable(ncid, dimids, "temperature_seasonal", "Seasonal mean temperature", "Degrees Celsius", dat);
  
  ncclose(ncid);
#else
  cerr<<"ERROR: No NetCDF support compiled.\n";
  exit(-1);
#endif

return 0;
}
