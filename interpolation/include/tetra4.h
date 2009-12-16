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

#ifndef H_TETRA4
#define H_TETRA4
// 4 node tetrahedron class.

#include <algorithm>
#include <vector>
#include <deque>
#include <string>
#include <math.h>
#include <iostream>

#include <assert.h>
#include <time.h>

template<class R>
class Tetra4{
 public:
  
  // Constructers/destructer
  Tetra4();
  Tetra4(const Tetra4& in); 
  ~Tetra4();
  
  // Overloaded operators.
  const Tetra4 &operator=(const Tetra4& in);
  bool operator==(const Tetra4& in) const;
  bool operator!=(const Tetra4& in) const;
  bool operator<(const Tetra4& in) const;
  
  void set_flags(const unsigned char flags);
  unsigned char get_flags() const;
  
  void set_nodeList(const std::vector<size_t>& in);
  const std::vector<size_t>& get_nodeList() const;
  
  unsigned get_size_enlist() const;
  int calcShapeFxn(const std::vector<R>& x0, const std::vector<R>& x1, 
		   const std::vector<R>& x2, const std::vector<R>& x3);

  R IO(const std::vector<R>& x) const;
  R IO(const std::vector<R>& x, size_t& Nn) const;

  int interpolate(const std::vector<R>& f0, const std::vector<R>& f1, 
		  const std::vector<R>& f2, const std::vector<R>& f3, 
		  const std::vector<R>& x,  std::vector<R>& f);

 private:
  bool GotShapeFxn;
  
  unsigned char flags;  
  std::vector<size_t> nodes;

  // Shape function stuff.
  double Volume, DET;
  std::vector<R> x0, x1, x2, x3;

  // know your orientation
  float flip;

  double get_L1(const std::vector<R>& x) const;
  double get_L2(const std::vector<R>& x) const;
  double get_L3(const std::vector<R>& x) const;
  
  double calc_vol(const std::vector<R>& r0, const std::vector<R>& r1, 
		  const std::vector<R>& r2, const std::vector<R>& r3) const;
};

template<class R>
Tetra4<R>::Tetra4(){
  GotShapeFxn = false;
  flip = 1.0;
  return;
}

template<class R>
Tetra4<R>::Tetra4(const Tetra4& in){
  *this = in;
  flip = 1.0;
  return;
}

template<class R>
Tetra4<R>::~Tetra4(){
  return;
}

template<class R>
const Tetra4<R>& Tetra4<R>::operator=(const Tetra4<R>& in){

  GotShapeFxn = in.GotShapeFxn;
  flags = in.flags;
  nodes = in.nodes;

  flip = in.flip;

  // Shape function stuff.
  Volume = in.Volume;
  DET = in.DET;
  x0 = in.x0;
  x1 = in.x1;
  x2 = in.x2;
  x3 = in.x3;

  return *this;
}

template<class R>
bool Tetra4<R>::operator==(const Tetra4<R>& in) const{
  
  unsigned len = nodes.size();
  if(len != in.nodes.size()) return false;
  
  for(unsigned i=0; i<len; i++)
    if(nodes[i] != in.nodes[i]) return false;
  
  return true;
}

template<class R>
bool Tetra4<R>::operator!=(const Tetra4<R>& in) const{
  return !(*this==in);
}

template<class R>
bool Tetra4<R>::operator<(const Tetra4<R>& in) const{
  unsigned len = nodes.size();
  if(len < in.nodes.size()) return true;
  if(len > in.nodes.size()) return false;
  
  for(unsigned i=0; i<len; i++)
    if(nodes[i] < in.nodes[i]) return true;
  
  return false; 
}

template<class R>
int Tetra4<R>::calcShapeFxn(const std::vector<R>& r0, const std::vector<R>& r1, 
			    const std::vector<R>& r2, const std::vector<R>& r3){
  // http://huizen.dto.tudelft.nl/deBruijn/programs/suna57.htm

  // The volume of an tetrahedral element is defined as:
  // V = |r_a r_b r_c| / 6
  // where r_a, r_b, r_c are vectors pointing from one the the
  // verticies to the other three verticies.
  
  DET = calc_vol(r0, r1, r2, r3);
  
  if(DET == 0.0){
    std::cerr<< "ERROR: Tetrahedral element has zero volume. "
	     << "Can't calculate shape fxn"<<std::endl;
    std::cerr<<r0[0]<<" "<<r0[1]<<" "<<r0[2]<<std::endl;
    std::cerr<<r1[0]<<" "<<r1[1]<<" "<<r1[2]<<std::endl;
    std::cerr<<r2[0]<<" "<<r2[1]<<" "<<r2[2]<<std::endl;
    std::cerr<<r3[0]<<" "<<r3[1]<<" "<<r3[2]<<std::endl;
    return(-1);
  }
  
  // fix orientation
  flip = DET/fabs(DET);
  DET *= flip;
  Volume = flip*DET/6.0;
  
  x0 = r0;
  x1 = r1;
  x2 = r2;
  x3 = r3;
  
  //   N(0) = 1-L1-L2-L3;  N(1) = L1;  N(2) = L2;  N(3) = L3.
  GotShapeFxn = true;
  return(0);
}

template<class R>
double Tetra4<R>::get_L1(const std::vector<R>& x) const{
  assert(GotShapeFxn);
  return calc_vol(x0, x, x2, x3)/DET;
}

template<class R>
double Tetra4<R>::get_L2(const std::vector<R>& x) const{
  assert(GotShapeFxn);
  return calc_vol(x0, x1, x, x3)/DET;
}

template<class R>
double Tetra4<R>::get_L3(const std::vector<R>& x) const{
  assert(GotShapeFxn);
  return calc_vol(x0, x1, x2, x)/DET;
}

template<class R>
R Tetra4<R>::IO(const std::vector<R>& x) const{
  assert(GotShapeFxn);

  // For IO(h,k,l) > 0 a point (h,k,l) is inside the element;
  // For IO(h,k,l) < 0 a point (h,k,l) is outside the element;
  // For IO(h,k,l) = 0 a point (h,k,l) is at the element boundary.
  // IO = min( h , k , l , 1-h-k-l ) .
  double L1 = get_L1(x);
  double L2 = get_L2(x);
  double L3 = get_L3(x);

  double L0 = 1.0-L1-L2-L3;
  
  double IO = L0;
  IO = (IO<L1)?IO:L1;
  IO = (IO<L2)?IO:L2;
  IO = (IO<L3)?IO:L3;
  
  // IO  > 0.0 => point inside element
  // IO == 0.0 => point on element boundary
  // IO  < 0.0 => point outside element
  return IO;
}

template<class R>
R Tetra4<R>::IO(const std::vector<R>& x, size_t& Nn) const{
  assert(GotShapeFxn);

  // For IO(h,k,l) > 0 a point (h,k,l) is inside the element;
  // For IO(h,k,l) < 0 a point (h,k,l) is outside the element;
  // For IO(h,k,l) = 0 a point (h,k,l) is at the element boundary.
  // IO = min( h , k , l , 1-h-k-l ) .
  double L1 = get_L1(x);
  double L2 = get_L2(x);
  double L3 = get_L3(x);

  double L0 = 1.0-L1-L2-L3;
  
  double IO = L0;
  Nn = 0;

  if(L1<IO){
    IO = L1;
    Nn = 1;
  }else if(L1==IO){
    // introduce sudo-randomness
    if(time(NULL)%2){
      IO = L1;
      Nn = 1; 
    }
  }

  if(L2<IO){
    IO = L2;
    Nn = 2;
  }else if(L2==IO){
    if(time(NULL)%2){
      IO = L2;
      Nn = 2;
    }
  }

  if(L3<IO){
    IO = L3;
    Nn = 3;
  }else if(L3==IO){
    if(time(NULL)%2){
      IO = L3;
      Nn = 3;
    }
  }
  
  // IO  > 0.0 => point inside element
  // IO == 0.0 => point on element boundary
  // IO  < 0.0 => point outside element
  return IO;
}

template<class R>
int Tetra4<R>::interpolate(const std::vector<R>& f0, const std::vector<R>& f1, 
			   const std::vector<R>& f2, const std::vector<R>& f3, 
			   const std::vector<R>& x,  std::vector<R>& f){
  unsigned len = f0.size();
  assert(len == f1.size());
  assert(len == f2.size());
  assert(len == f3.size());
  assert(x.size() == 3);
  
  // N(1) = xi  ;  N(2) = eta  ;  N(3) = zeta  ;  N(4) = 1-xi-eta-zeta.
  double L1 = get_L1(x);
  double L2 = get_L2(x);
  double L3 = get_L3(x);
  double N0 = 1.0-L1-L2-L3;
  double N1 = L1;
  double N2 = L2;
  double N3 = L3;

  f.resize(len);
  for(unsigned i=0; i<len; i++){
    f[i] = N0*f0[i] + N1*f1[i] + N2*f2[i] + N3*f3[i];
  }
  
  return(0);
}

template<class R>
const std::vector<size_t>& Tetra4<R>::get_nodeList() const{
  return nodes;
}

template<class R>
void Tetra4<R>::set_nodeList(const std::vector<size_t>& in){
  nodes = in;
  return;
}

template<class R>
double Tetra4<R>::calc_vol(const std::vector<R>& r0, const std::vector<R>& r1, 
						   const std::vector<R>& r2, const std::vector<R>& r3) const{
  double a1 = r1[0]-r0[0];
  double a2 = r1[1]-r0[1];
  double a3 = r1[2]-r0[2];
  
  double b1 = r2[0]-r0[0];
  double b2 = r2[1]-r0[1];
  double b3 = r2[2]-r0[2];
  
  double c1 = r3[0]-r0[0];
  double c2 = r3[1]-r0[1];
  double c3 = r3[2]-r0[2];
  
  // volume = | r_a r_b r_c | / 6
  double det = a1*(b2*c3 - b3*c2) - b1*(a2*c3 - a3*c2) + c1*(a2*b3 - a3*b2);  
  return det*flip;
}

#endif
