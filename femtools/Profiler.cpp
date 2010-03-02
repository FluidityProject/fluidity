/* Copyright (C) 2010- Imperial College London and others.
   
   Please see the AUTHORS file in the main source directory for a full
   list of copyright holders.
   
   Dr Gerard J Gorman
   Applied Modelling and Computation Group
   Department of Earth Science and Engineering
   Imperial College London
   
   g.gorman@imperial.ac.uk
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
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

#include "Profiler.h"

using namespace std;

Profiler::Profiler(){}

Profiler::~Profiler(){}

double Profiler::get(const std::string &key) const{
  double time = timings.find(key)->second.second;
  if(MPI::Is_initialized()){
    double gtime;
    MPI::COMM_WORLD.Reduce(&time, &gtime, 1, MPI::DOUBLE, MPI::SUM, 0);
    return gtime;
  }
  
  return time; 
}

void Profiler::print() const{
  for(map< string, pair<double, double> >::const_iterator it=timings.begin();it!=timings.end();++it){
    cout<<it->first<<" :: "<<it->second.second<<endl;
  }
}

void Profiler::tic(const std::string &key){
  timings[key].first = wall_time();
}

void Profiler::toc(const std::string &key){
  timings[key].second += wall_time() - timings[key].first;
}

double Profiler::wall_time() const{
#ifdef HAVE_MPI
  return MPI::Wtime();
#else
  return 0.0;
#endif
}

void Profiler::zero(){
  for(map< string, pair<double, double> >::iterator it=timings.begin();it!=timings.end();++it){
    it->second.second = 0.0;
  }
}

void Profiler::zero(const std::string &key){
  timings[key].second = 0.0;
}

// Opaque instance of profiler.
Profiler flprofiler;

// Fortran interface
extern "C" {
#define cprofiler_get_fc F77_FUNC(cprofiler_get, CPROFILER_GET)
  void cprofiler_get_fc(const char *key, const int *key_len, double *time){
    *time = flprofiler.get(string(key, *key_len));
  }

#define cprofiler_tic_fc F77_FUNC(cprofiler_tic, CPROFILER_TIC)
  void cprofiler_tic_fc(const char *key, const int *key_len){
    flprofiler.tic(string(key, *key_len));
  }

#define cprofiler_toc_fc F77_FUNC(cprofiler_toc, CPROFILER_TOC)
  void cprofiler_toc_fc(const char *key, const int *key_len){
    flprofiler.toc(string(key, *key_len));
  }

#define cprofiler_zero_fc F77_FUNC(cprofiler_zero, CPROFILER_ZERO)
  void cprofiler_zero_fc(){
    flprofiler.zero();
  }
}
