/* Copyright (C) 2006 Imperial College London and others.

 Please see the AUTHORS file in the main source directory for a full list
 of copyright holders.

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
#include <mpi.h>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <set>
#include <cassert>
#include <cstdlib>

#include "c++debug.h"
#include "Element.h"

using namespace std;

Element::Element(){
}
Element::Element(const Element& element){
  *this = element;
}
Element::~Element(){}

// Assignment operator for Element.
Element &Element::operator=( const Element &e ){
  eid      = e.eid;
  flags    = e.flags;
  nodes    = e.nodes;
  MFnodes  = e.MFnodes;
  ifields  = e.ifields;
  fields   = e.fields;
  
  return( *this );
}

// Equivalence operator for elements
// Two elements are considered to be equivalent if they contain the same nodes.
bool Element::operator==(const Element& elem) const{
  unsigned len = nodes.size();
  
  // If these are not the same length then obviously...
  if(len != elem.nodes.size())
    return(false);

  // Make no assumption about the order
  for(unsigned i=0; i<len; i++){
    bool found = false;
    for(unsigned j=0; j<len; j++){
      found = (nodes[i] == elem.nodes[j]);
      if(found)
	break;
    }
    if(!found)
      return(false);
  }

  len = MFnodes.size();
  { // Mixed-Formulation
  
    // If these are not the same length then obviously...
    if(len != elem.MFnodes.size())
      return(false);
    
    // Make no assumption about the order
    for(unsigned j=0; j<len; j++){
      bool found = false;
      for(unsigned k=0; k<len; k++){
	found = (MFnodes[j] == elem.MFnodes[k]);
	if(found)
	  break;
      }
      if(!found)
	return(false);
    }
  }
  return(true);
}

bool Element::operator!=(const Element& elem) const{
  return( !(*this == elem) );
}

bool Element::operator<(const Element& in) const{
  
  if(*this == in) 
    return false;
  
  unsigned len = nodes.size();
  if(len != in.nodes.size())
    return(len<in.nodes.size());
  
  {
    set<unsigned> set1;
    set<unsigned> set2;
    
    for(unsigned i=0; i<len; i++){
      set1.insert(nodes[i]);
      set2.insert(in.nodes[i]);
    }
    
    if(set1!=set2)
      return(set1<set2);
  }
  
  len = MFnodes.size();
  if(len != in.MFnodes.size())
    return(len<in.MFnodes.size());
  
  {
    set<unsigned> set1;
    set<unsigned> set2;
    
    for(unsigned i=0; i<len; i++){
      set1.insert( MFnodes[i] );
      set2.insert( in.MFnodes[i] );
    }
    
    if(set1!=set2)
      return(set1<set2);
  }

  cerr<<"ERROR: The programmer is a brainless git"<<endl
      <<__FILE__<<" "<<__LINE__<<endl
      <<*this<<endl
      <<in<<endl;
  if(*this==in){
      cerr<<"they are the same"<<endl;
  }else{
      cerr<<"they are different"<<endl;
  }
  
  exit(-1);
  
  return false;
}

void Element::append_field(const std::vector<samfloat_t>& _field){
  fields.insert(fields.end(), _field.begin(), _field.end());
  
}

void Element::append_field(const samfloat_t* _field, const size_t block_size){
  for(size_t i=0;i<block_size;i++)
    fields.push_back(_field[i]);
}

void Element::set_eid(const eid_t _eid){ 
  eid = _eid; 
}
eid_t Element::get_eid() const{ 
  return eid; 
}

void Element::set_flags(const unsigned char _flags){ 
  flags = _flags; 
}
unsigned char Element::get_flags() const{ 
  return flags; 
}

void Element::set_enlist(const vector<unn_t>& in){ 
  nodes = in; 
}
void Element::set_MFenlist(const vector<unn_t>& in){
  MFnodes = in; 
}
//void          Element::set_enlist(const unn_t* _unn, const unsigned len){
//  nodes.resize(len);
//  for(unsigned i=0; i<len; i++)
//    nodes[i] = _unn[i];
//}
const vector<unn_t>& Element::get_enlist() const{
  return nodes; 
}
//const unn_t*              Element::get_cptr_enlist() const{ 
//  return nodes.begin(); 
//}
const vector<unn_t>& Element::get_MFenlist() const{
  return( MFnodes );
}

unsigned Element::get_size_enlist() const{ 
  return nodes.size(); 
}
unsigned Element::get_size_MFenlist() const{
  return MFnodes.size(); 
}

void Element::set_ifields(const vector<int>& in){ 
  ifields = in; 
}
const vector<int>& Element::get_ifields() const{ 
  return ifields; 
}
size_t Element::get_size_ifields() const{ 
  return ifields.size(); 
}

void Element::set_fields(const vector<samfloat_t>& in){ 
  fields = in; 
}
const vector<samfloat_t>& Element::get_fields() const{ 
  return fields;
}
size_t Element::get_size_fields() const{ 
  return fields.size(); 
}

void Element::pack(char *buffer, int& bsize, int& offset) const{  
  ECHO("packing..." << *this);
  
  //  MPI::EID_T.Pack(&eid,           1, buffer, bsize, offset, MPI::COMM_WORLD);
  MPI::UNSIGNED_CHAR.Pack(&flags, 1, buffer, bsize, offset, MPI::COMM_WORLD);
  
  { // Pack element-node list for this element
    unsigned ncnt = nodes.size();
    MPI::UNSIGNED.Pack(&ncnt,         1, buffer, bsize, offset, MPI::COMM_WORLD);
    MPI::UNN_T.Pack(&(nodes[0]), ncnt, buffer, bsize, offset, MPI::COMM_WORLD);
  }
  
  { // Pack Mixed-Formulation based element-node lists for this element
    unsigned ncnt = MFnodes.size();
    MPI::UNSIGNED.Pack(&ncnt,          1, buffer, bsize, offset, MPI::COMM_WORLD);
    if( ncnt > 0 )
      MPI::UNN_T.Pack(&(MFnodes[0]), ncnt, buffer, bsize, offset, MPI::COMM_WORLD);
  }
  
  { // Pack integer fields
    unsigned icnt = ifields.size();
    MPI::UNSIGNED.Pack(&icnt,         1, buffer, bsize, offset, MPI::COMM_WORLD);
    MPI::INT.Pack(&(ifields[0]), icnt, buffer, bsize, offset, MPI::COMM_WORLD);
  }

  { // Pack real fields field
    unsigned fcnt = fields.size();
    MPI::UNSIGNED.Pack(&fcnt,            1, buffer, bsize, offset, MPI::COMM_WORLD);
    MPI::SAMFLOAT.Pack(&(fields[0]), fcnt, buffer, bsize, offset, MPI::COMM_WORLD);
  }

}

void Element::unpack(char *buffer, int& bsize, int& offset){
  
  // MPI::EID_T.Unpack(buffer,         bsize, &eid,  1, offset, MPI::COMM_WORLD);
  MPI::UNSIGNED_CHAR.Unpack(buffer, bsize, &flags,1, offset, MPI::COMM_WORLD);
  
  { // unpack element-node list
    unsigned ncnt;
    MPI::UNSIGNED.Unpack(buffer, bsize, &ncnt, 1, offset, MPI::COMM_WORLD);
    nodes.resize(ncnt);
    MPI::UNN_T.Unpack(buffer, bsize, &(nodes[0]), ncnt, offset, MPI::COMM_WORLD);
  }
  
  { // unpack pressure element-node list
    unsigned ncnt;
    MPI::UNSIGNED.Unpack(buffer, bsize, &ncnt, 1, offset, MPI::COMM_WORLD);
    
    if( ncnt > 0 ){
      MFnodes.resize( ncnt );
      MPI::UNN_T.Unpack(buffer, bsize, &(MFnodes[0]), ncnt, offset, MPI::COMM_WORLD);
    }
  }

  { // unpack the integer fields
    unsigned icnt;
    MPI::UNSIGNED.Unpack(buffer, bsize, &icnt, 1, offset, MPI::COMM_WORLD);
    ifields.resize(icnt);
    MPI::INT.Unpack(buffer, bsize, &(ifields[0]), icnt, offset, MPI::COMM_WORLD);
  }

  { // unpack real field values
    unsigned fcnt;
    MPI::UNSIGNED.Unpack(buffer, bsize, &fcnt, 1, offset, MPI::COMM_WORLD);
    fields.resize(fcnt);
    MPI::SAMFLOAT.Unpack(buffer, bsize, &(fields[0]), fcnt, offset, MPI::COMM_WORLD);
  }

}

//
// Get an estimate of the number of bytes required to pack an element.
//
unsigned Element::pack_size() const{  
  int total;
  
  total  = MPI::EID_T.Pack_size(1, MPI::COMM_WORLD);
  total += MPI::UNSIGNED_CHAR.Pack_size(1, MPI::COMM_WORLD);
  
  { // Pack_size element-node list for this element
    unsigned ncnt = nodes.size();
    total += MPI::UNSIGNED.Pack_size(1, MPI::COMM_WORLD);
    total += MPI::UNN_T.Pack_size(ncnt, MPI::COMM_WORLD);
  }
  
  { // Pack_size Mixed-Formulation based element-node lists for this element
    unsigned ncnt = MFnodes.size();
    total += MPI::UNSIGNED.Pack_size(1, MPI::COMM_WORLD);
    if( ncnt > 0 )
      total += MPI::UNN_T.Pack_size(ncnt, MPI::COMM_WORLD);
  }
  
  { // Pack_size the ifields
    unsigned icnt = ifields.size();
    total += MPI::UNSIGNED.Pack_size(1, MPI::COMM_WORLD);
    if( icnt>0 )
      total += MPI::INT.Pack_size(icnt, MPI::COMM_WORLD);
  }
  
  { // Pack_size element-wise field values
    unsigned fcnt = fields.size();
    total += MPI::UNSIGNED.Pack_size(1, MPI::COMM_WORLD);
    if( fcnt>0 )
      total += MPI::SAMFLOAT.Pack_size(fcnt, MPI::COMM_WORLD);
  }
  
  return total;
}

// printing operator for Element
ostream& operator<<(ostream& s, const Element& e){
  s << "Element Profile" << endl;
  
  s << "ID number\t\t" << e.eid << endl;
  s << "Flags\t\t\t" << (unsigned)(e.flags) << endl;
  
  { // Node info.
    s << "Number of regular nodes in element\t" << e.nodes.size() << endl;
    s << "{";
    for(vector<unn_t>::const_iterator it = e.nodes.begin(); it != e.nodes.end(); ++it)
      s << *it << " ";
    s << "}" << endl;
  }

  { // Mixed-formulation node.
    unsigned len = e.MFnodes.size();
    if(len>0){
      s << "Mixed-formulation element:";
      for(unsigned i = 0; i<len; i++)
	s << e.MFnodes[i] << " ";
      s << "}" << endl;
    }
  }

  { // integer fields
    if(e.ifields.size() ){
      s << "Element integer fields: ";
      for(vector<int>::const_iterator it = e.ifields.begin(); it != e.ifields.end(); ++it){
	s << *it <<"\t";
      }
      s<<endl;
    }
  }
  
  { // real fields
    if(e.fields.size() ){
      s << "Element real fields: ";
      for(vector<samfloat_t>::const_iterator it = e.fields.begin(); it != e.fields.end(); ++it){
	s << *it <<"\t";
      }
      s<<endl;
    }
  }
  
  return(s);
}
