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
#include "ParticleList.h"
#include <fstream>
#include "confdefs.h"
#include <cstdlib>
#include <map>
#include <iomanip>

using namespace std;

#ifndef ParticleWriterStat_H
#define ParticleWriterStat_H
class ParticleWriterStat
{
 public:
  ParticleWriterStat(ParticleList *plist, char *filename, bool binary);
  ~ParticleWriterStat();

  void register_field(char *name, void *ptr, char *mphase, int components);
  void write_header();
  void write_particle_state(double time, double dt);

 private:
  string filename;
  ofstream file;
  bool is_binary;

  ParticleList *plist;
  map<string, void*> field_ptrs;
  map<string, int> field_comps;
  map<string, string> field_mphase;
  int column_cnt;

  void write_header_constants();
  void write_header_fields();
  string constant_tag(string name, string type, string value);
  string constant_tag(string name, string type, int value);
  string field_tag(string name, int column, string statistic,
                   string material_phase, int components);

  void write_scalar(double value);
  void write_vector(vector<double> values);
  void write_vector(double *values, size_t size);
};
#endif // ParticleWriterStat_H
