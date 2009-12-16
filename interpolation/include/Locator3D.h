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

#ifndef H_LOCATOR
#define H_LOCATOR
#include "confdefs.h"

#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef DOUBLEP
typedef double gigreal_t;
#else
typedef float gigreal_t;
#endif

#include "nodes.h"
#include "tetra4.h"

class Locator3D{
 public:
  Locator3D();

  int get_node2element(int) const;
  
  int interpolate_linear(gigreal_t *, int *, int, gigreal_t *,  int *);
  
  int t4_to_t4_Search(int, int, const gigreal_t [], const gigreal_t [], const gigreal_t [], const int[],
		      int, int, const gigreal_t [], const gigreal_t [], const gigreal_t [], const int[]);
  
  int test_gid(std::string, std::string);
  
  void verbose_off();
  
  void verbose_on();
  
 private:  
  size_t bruteforceSearch(const std::deque< Tetra4<gigreal_t> >& elements, 
			  const std::vector<gigreal_t> x);
  
  int GiD_read_msh(const std::string& fileName, 
		   std::deque< Tetra4<gigreal_t> >& ElementList, 
		   std::deque< Node<gigreal_t> >& nodeList);
  
  int interpolate_linear();
  
  int mkEElist(const std::deque< Tetra4<gigreal_t> >& elements,
	       std::deque< std::vector<int> >& EElist);
  
  int mkNNlist(const std::deque< Tetra4<gigreal_t> >& elements,
	       std::deque< std::vector<int> >& NNlist);
  
  void Tokenize(const std::string& str,
		std::vector<std::string>& tokens);
  
  int vicinitySearch(const std::deque< Tetra4<gigreal_t> >& elements,
		     const std::deque< std::vector<int> >& EElist,
		     const std::vector<gigreal_t>& x,
		     const size_t start,
		     const int NTries,
		     size_t& last);
  
  int t4_to_t4_Search();
  
  std::vector<int> node2element;
  std::deque< Node<gigreal_t> > nodes_1;
  std::deque< Tetra4<gigreal_t> > elements_1;
  
  std::deque< Node<gigreal_t> > nodes_2;
  std::deque< Tetra4<gigreal_t> > elements_2;
  
  bool verbose;
};

#endif
