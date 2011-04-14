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
#include <string>
using std::string;

#include <assert.h>
#include "unn_t.h"

// default constructor
unn_t::unn_t(){
  s.resize(UNN_LEN+1);
}

// copy constructor
unn_t::unn_t(const unn_t& unn){
  *this = unn;
}

// destructor
unn_t::~unn_t(){}

// return the unn as a c type string
const char *unn_t::number(){
  return(s.c_str());
}

// Overloaded operators.
unn_t & unn_t::operator=(const unn_t &a){
  s = a.s;
  return(*this);
}
unn_t & unn_t::operator=(const string &a){
  char num[UNN_LEN+1];
  
  if( a.size() < UNN_LEN ){
    sprintf(num, "%*s", UNN_LEN, a.c_str());
    for(int i=0;i<UNN_LEN;i++){
      if(num[i] == ' ')
	num[i] = '0';
      else
	break;
    }
  }else{
    strncpy(num, a.c_str(), UNN_LEN);
  }
  num[UNN_LEN] = '\0';
  
  s = num;
  return(*this);
}
unn_t & unn_t::operator=(const int &a){
  char num[UNN_LEN+1];
  sprintf(num, "%0*d", UNN_LEN, a);
  s = num;
  return(*this);
}
