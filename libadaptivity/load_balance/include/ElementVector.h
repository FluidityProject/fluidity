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
#ifndef ELEMENTVECTOR_H
#define ELEMENTVECTOR_H

#include "confdefs.h"

#include <vector>
#include <deque>
#include <iostream>
#include "samtypes.h"

//
// This class describes a template container for storing elements.
//
template<class ElementType>
class ElementVector: public std::deque<ElementType>{
 public:
  ElementVector();
  ElementVector(const ElementVector&);
  ~ElementVector();
  
  // Overloaded operators.
  ElementVector &operator=(const ElementVector<ElementType> &in);
#ifndef NDEBUG
  const ElementType &operator[](const unsigned in) const;
  ElementType &operator[](const unsigned in);
#endif

  // Other member functions.
  void shrink(const unsigned);
  
  void push_back(const ElementType& in);
  unsigned size() const;
  void clear();
  
 private:
  unsigned len;
};

template<class ElementType>
ElementVector<ElementType>::ElementVector(){
  len=0;
}

template<class ElementType>
ElementVector<ElementType>::ElementVector(const ElementVector& newkid){
  *this = newkid;
}

template<class ElementType>
ElementVector<ElementType>::~ElementVector(){
  clear();
}

// Overloaded operators.
template<class ElementType>
ElementVector<ElementType>& ElementVector<ElementType>::operator=(const ElementVector<ElementType>& in){
  len = in.len;
  insert( (*this).end(), in.begin(), in.end() );
  
  return *this;
}

#ifndef NDEBUG
template<class ElementType>
ElementType& ElementVector<ElementType>::operator[](const unsigned i){
  // Don't forgive out of bounds referencing.
  assert(i<len); 
  return *(this->begin() + i);
}

template<class ElementType>
const ElementType& ElementVector<ElementType>::operator[](const unsigned i) const{
  
  // Don't forgive out of bounds referencing.
  assert( (i>=0)&&(i<len)); 
  
  return *(this->begin() + i);
}
#endif

// Other member functions
template<class ElementType>
void ElementVector<ElementType>::push_back(const ElementType& in){
  insert(this->end(), in);
  len++;
}

template<class ElementType>
void ElementVector<ElementType>::shrink(const unsigned newlen){
  if(newlen == len)
    return;
  
  assert(newlen < len);
  this->resize( newlen );
  len = newlen;
}

template<class ElementType>
unsigned ElementVector<ElementType>::size() const{
  return len;
}

template<class ElementType>
void ElementVector<ElementType>::clear(){
  len = 0;
  erase(this->begin(), this->end());
}

#endif


