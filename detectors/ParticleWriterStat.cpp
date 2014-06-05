/*  Copyright (C) 2009 Imperial College London and others.
    
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
#include "ParticleWriterStat.h"

ParticleWriterStat::ParticleWriterStat(ParticleList *plist, char *filename, bool binary) {
  this->filename = string(filename);
  this->is_binary = binary;
  this->plist = plist;
  column_cnt = 1;

  file.open(filename, ios::out);
}

ParticleWriterStat::~ParticleWriterStat()
{
  file.close();
}

void ParticleWriterStat::register_field(char *name, void *ptr,
                                        char *mphase, int components)
{
  string fieldname = string(name);
  field_ptrs.insert( pair<string, void*>(fieldname, ptr) );
  field_comps.insert( pair<string, int>(fieldname, components) );
  field_mphase.insert( pair<string, string>(fieldname, string(mphase)) );
}

void ParticleWriterStat::write_header()
{
  file << "<header>" << endl;
  write_header_constants();
  write_header_fields();
  file << "</header>" << endl;

  /* Re-open file for binary output */
  if (is_binary) {
    file.close();
    string fname = filename.append(".dat");
    file.open(fname.c_str(), ios::out | ios::binary);
  }
}

void ParticleWriterStat::write_header_constants()
{
#ifdef __FLUIDITY_VERSION__
  string version = __FLUIDITY_VERSION__;
#else
  string version = "";
#endif
  file << constant_tag("FluidityVersion", "string", version) << endl;

  string datetime = string(__DATE__).append(" ").append(__TIME__);
  file << constant_tag("CompileTime", "string", datetime) << endl;

  time_t t = time(0);
  string time = ctime(&t);
  time.erase(time.find("\n", 0), 1);  // Remove \n from ctime() string
  file << constant_tag("StartTime", "string", time) << endl;

  string hostname = "Unknown";
  char *host = getenv("HOSTNAME");

  if (host != NULL) hostname = string(host);
  file << constant_tag("HostName", "string", hostname) << endl;

  if (is_binary) {
    // TODO: Binary data types...
  } else {
    file << constant_tag("format", "string", "plain_text") << endl;
  }
}

void ParticleWriterStat::write_header_fields()
{
  /* Header for timekeeping */
  file << field_tag("ElapsedTime", column_cnt++, "value", string(), 1) << endl;
  file << field_tag("dt", column_cnt++, "value", string(), 1) << endl;

  /* Header for detector position */
  for (list<Particle*>::iterator p=plist->begin(); p != plist->end(); ++p) {
    string name = (*p)->get_name();
    int components = (*p)->get_position().size();
    file << field_tag(name, column_cnt, "position", string(), components) << endl;
    column_cnt += components;
  }

  /* Header for detector fields */
  for (map<string, string>::iterator field=field_mphase.begin();
       field!=field_mphase.end(); ++field) {

    for (list<Particle*>::iterator p=plist->begin(); p != plist->end(); ++p) {
      string name = field->first;
      string mphase = field->second;
      string statistic = (*p)->get_name();
      map<string, int>::iterator comp = field_comps.find(name);
      int components = comp->second;
      file << field_tag(name, column_cnt, statistic, mphase, components) << endl;
      column_cnt += components;
    }
  }
}

string ParticleWriterStat::constant_tag(string name, string type, string value)
{
  stringstream tag;
  tag << "<constant name=\"" << string(name) << "\" type=\"" << type;
  tag << "\" value=\"" << value << "\" />";
  return tag.str();
}

string ParticleWriterStat::field_tag(string name, int column, string statistic,
                                     string material_phase, int components)
{
  stringstream tag;
  tag << "<field column=\"" << column << "\" name=\"" << name;
  tag << "\" statistic=\"" << statistic << "\"";
  if (!material_phase.empty()) tag << " material_phase=\"" << material_phase << "\"";
  if (components > 1) tag << " components=\"" << components << "\"";
  tag << "/>";
  return tag.str();
}

void ParticleWriterStat::write_particle_state(double time, double dt)
{
  /* Timekeeping */
  write_scalar(time);
  write_scalar(dt);

  /* Particle position */
  for (list<Particle*>::iterator p=plist->begin(); p != plist->end(); ++p) {
    write_vector((*p)->get_position());
  }

  /* Sampled field data */
  for (map<string, void*>::iterator field=field_ptrs.begin();
       field!=field_ptrs.end(); ++field) {
    void *ptr = field->second;
    int comps = field_comps.find(field->first)->second;
    if (comps > 1) {
      /* Vector fields */
      double *vector_value = new double[comps];
      for (list<Particle*>::iterator p=plist->begin(); p != plist->end(); ++p) {
        (*p)->vector_field_value(ptr, vector_value);
        write_vector(vector_value, comps);
      }
    } else {
      /* Scalar fields */
      double scalar_value = 0.;
      for (list<Particle*>::iterator p=plist->begin(); p != plist->end(); ++p) {
        (*p)->scalar_field_value(ptr, &scalar_value);
        write_scalar(scalar_value);
      }
    }
  }
  file << endl;
}

void ParticleWriterStat::write_scalar(double value)
{
  file << "  " << std::fixed << std::setprecision(6) << value;
}

void ParticleWriterStat::write_vector(vector<double> values)
{
  for (size_t i=0; i<values.size(); ++i) write_scalar(values[i]);
}

void ParticleWriterStat::write_vector(double *values, size_t size)
{
  for (size_t i=0; i<size; ++i) write_scalar(values[i]);
}

extern "C" {
  ParticleWriterStat *create_particle_writer_stat(ParticleList *plist,
                                                  char *filename, int binary)
  {
    if (binary==0) return new ParticleWriterStat(plist, filename, true);
    else           return new ParticleWriterStat(plist, filename, false);
  }

  void particle_writer_register_field(ParticleWriterStat *writer, char *field, 
                                      void *field_ptr, char *mphase, int components)
  {
    writer->register_field(field, field_ptr, mphase, components);
  }

  void particle_writer_write_header(ParticleWriterStat *writer)
  {
    writer->write_header();
  }

  void particle_writer_write_state(ParticleWriterStat *writer, double time, double dt)
  {
    writer->write_particle_state(time, dt);
  }
}
