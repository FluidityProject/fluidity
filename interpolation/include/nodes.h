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

#ifndef H_NODES
#define H_NODES
// Node class

#include <cassert>
#include <vector>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <math.h>
#include <iostream>

template<class R>
class Node{

 public:
  Node(void);                // "No argument" constructor
  Node(const Node& in);      // Copy constructor
  ~Node(void);               // destructer
  
  // Member functions.
  const Node& operator=(const Node &in); 
  bool operator==(const Node& in) const;
  bool operator!=(const Node& in) const;
  bool operator<(const Node& in) const;
  Node<R> operator-(const Node<R>& rhs) const;

  void set_gnn(const unsigned eid);
  unsigned get_gnn() const;

  void set_flags(const unsigned char);
  unsigned char get_flags() const;

  void set_start_element(const unsigned eid);
  unsigned get_start_element() const;

  void set_fields(const std::vector<R>& in);
  void set_fields(const R* in, const unsigned inlen);
  void add_field(R in);
  
  const std::vector<R>& get_fields() const;
  unsigned get_size_fields() const;

  void set_coord(const std::vector<R>&);
  void set_coord(const R, const R);
  void set_coord(const R, const R, const R);

  const std::vector<R>& get_coord() const;
  R get_x() const;
  R get_y() const;
  R get_z() const;
  
 private:
  unsigned gnn;
  std::vector<R> x;         // position
  unsigned char  flags;      
  std::vector<R> fields;
  unsigned eguess;
};

template<class R>
Node<R>::Node(void){
  flags=0;
  eguess=0;
  return;
}

template<class R>
Node<R>::Node(const Node& in){
  *this = in;
  return;
}

template<class R>
Node<R>::~Node(void){
  return;
}

template<class R>
const Node<R> &Node<R>::operator=(const Node<R>& in){
  gnn    = in.gnn;
  x      = in.x;     
  flags  = in.flags;      
  fields = in.fields;
  eguess = in.eguess;
  return *this;
}

template<class R>
bool Node<R>::operator==(const Node<R>& in) const{
  unsigned len = x.size();
  if(len != in.x.size()) return false;
  for(unsigned i=0; i<len; i++)
    if(x[i] != in.x[i]) return false;
  
  if(flags != in.flags) return false;

  len = fields.size();
  if(len != in.fields.size()) return false;
  for(unsigned i=0; i<len; i++)
    if(fields[i] != in.fields[i]) return false;

  return true;
}

template<class R>
bool Node<R>::operator!=(const Node<R>& in) const{
  return !(*this == in);
}

template<class R>
bool Node<R>::operator<(const Node<R>& in) const{

  unsigned len = x.size();
  if(len != in.x.size())
    return len < in.x.size();
  
  for(unsigned i=0; i<len; i++)
    if(x[i] != in.x[i]) 
      return x[i] < in.x[i];
  
  if(flags != in.flags) 
    return flags < in.flags;
  
  len = fields.size();
  if(len != in.fields.size()) 
    return len < in.fields.size();
  
  for(unsigned i=0; i<len; i++)
    if(fields[i] != in.fields[i]) 
      return fields[i] != in.fields[i];
  
  // if we get here it means that they're equal!
  return false;
}

template<class R>
void Node<R>::set_flags(const unsigned char in){
  flags = in;
  return;
}

template<class R>
unsigned char Node<R>::get_flags() const{
  return flags;
}

template<class R>
void Node<R>::set_fields(const std::vector<R>& in){
  fields = in;
  return;
}

template<class R>
void Node<R>::set_fields(const R* in, const unsigned inlen){
  fields.resize(inlen);
  for(unsigned i=0; i<inlen; i++)
    fields[i] = in[i];

  return;
}

template<class R>
void Node<R>::add_field(R in){
  fields.push_back(in);
  return;
}

template<class R>
const std::vector<R>& Node<R>::get_fields() const{
  return fields;
}

//template<class R>
//std::ostream & Node<R>::operator<<(std::ostream& out, const Node& in){
//  out << "x = ";
//  for(std::vector<R> it=nodes.begin(); it!=nodes.end(); ++it)
//    out<<*it<<" ";
//  out <<endl;
//
//  out << "fields = ";
//  for(std::vector<R> it=fields.begin(); it!=fields.end(); ++it)
//    out<<*it<<" ";
//  out << endl;
//
//  return out;
//}

template<class R>
unsigned Node<R>::get_size_fields() const{
  return fields.size();
}

template<class R>
void Node<R>::set_coord(const std::vector<R>& in){
  x = in;
}

template<class R>
void Node<R>::set_coord(const R x0, const R x1){
  x.resize(2);
  x[0] = x0;
  x[1] = x1;

  return;
}

template<class R>
void Node<R>::set_coord(const R x0, const R x1, const R x2){
  x.resize(3);
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
  
  return;
}

template<class R>
const std::vector<R>& Node<R>::get_coord() const{
  return x;
}

template<class R>
R Node<R>::get_x() const{
  assert(x.size() > 0);
  return x[0];
}

template<class R>
R Node<R>::get_y() const{
  assert(x.size() > 1);
  return x[1];
}

template<class R>
R Node<R>::get_z() const{
  assert(x.size() > 2);
  return x[2];
}

template<class R>
void Node<R>::set_start_element(const unsigned eid){
  eguess = eid;
  return;
}

template<class R>
unsigned Node<R>::get_start_element() const {
  return(eguess);
}

template<class R>
void Node<R>::set_gnn(const unsigned in){
  gnn = in;
  return;
}

template<class R>
unsigned Node<R>::get_gnn() const {
  return(gnn);
}

template<class R>
Node<R> Node<R>::operator-(const Node<R> &rhs) const{
  unsigned len = fields.size();
  assert(len == rhs.fields.size());
  
  Node<R> node( *this );
  for(unsigned i=0; i<len; i++)
    node.fields[i] -= rhs.fields[i];
  
  return(node);
}

#endif











