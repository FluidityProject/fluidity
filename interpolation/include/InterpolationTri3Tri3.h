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

#ifndef INTERPOLATIONTRI3TRI3_H
#define INTERPOLATIONTRI3TRI3_H

template<class R>
class InterpolationTri3Tri3;

#include <cassert>
#include <vector>
#include <set>
#include <utility>

template<class R>
class InterpolationTri3Tri3{
 public:
  InterpolationTri3Tri3(int _offset){offset = _offset;NBruteForce=0;}
  void Interpolate(const R *source, R* dest);
  void SetDest(const int *enlist, const R *X, const R *Y,
		 int NNodes, int NElements);
  void SetSource(const int *enlist, const R *X, const R *Y,
		 int NNodes, int NElements);
 private:
  InterpolationTri3Tri3(const InterpolationTri3Tri3&);
  const InterpolationTri3Tri3& operator=(const InterpolationTri3Tri3&);

  void advancing_front_search();
  double area(int eid);
  double area(int n1, int n2, int n3);
  double area(int n1, int n2, R x, R y);
  double area(R x1, R y1,
	      R x2, R y2,
	      R x3, R y3);
  int bruteForce(R x, R y);
  void find_elements();
  std::pair<int, double> IO(int eid, R x, R y);
  void shapeFunction(int eid, R x, R y, double *L);
  std::pair<int, double> vicinity(int start, R x, R y, int tries);
  int offset,NBruteForce;

  int sourceNNodes, sourceNElements;
  std::vector<int> sourceENList;
  std::vector<R> sourceX, sourceY;
  std::vector<int> sourceEEList;

  std::vector<double> areas;

  int destNNodes, destNElements;
  std::vector<int> destENList;
  std::vector<R> destX, destY;
  std::vector< std::set<int> > destNNList;

  std::vector<int> element_containing_node;
  std::vector<double> L012;
};

template<class R>
void InterpolationTri3Tri3<R>::advancing_front_search(){
  
  { // Create the Element-Element adjacency lists
    std::vector< std::set<int> > NEList(sourceNNodes);
    for(int i=0; i<sourceNElements; i++){
      for(int j=0; j<3; j++){
	int n = sourceENList[i*3+j]-offset;
	assert(n>=0); assert(n<sourceNNodes);
	NEList[n].insert(i);
      }
    }
    
    // We want the n'th element in the element-element adjacency list to
    // be opposit the n'th node in the element
    sourceEEList.resize(sourceNElements*3);
    for(int i=0; i<sourceNElements; i++){
      for(int j=0; j<3; j++){
	int elm=-1;
	int n1 = sourceENList[i*3 + ((j+1)%3)] - offset;
	int n2 = sourceENList[i*3 + ((j+2)%3)] - offset;
	for(std::set<int>::iterator it=NEList[n1].begin(); it!=NEList[n1].end(); ++it){
	  if((*it)==i)
	    continue;
	  
	  if(NEList[n2].find(*it)!=NEList[n2].end()){
	    elm = *it;
	    break;
	  }
	}
	sourceEEList[i*3+j] = elm;
      }
    }
  }

  { // create node-node adjancy list
    destNNList.resize(destNNodes);
    for(int i=0; i<destNElements; i++){
      for(int j=0; j<3; j++){
	int n1 = destENList[i*3 + j]-offset;
	assert(n1>=0); assert(n1<destNNodes);
	for(int k=0; k<3; k++){
	  if(j!=k){
	    int n2 = destENList[i*3 + k]-offset;
	    assert(n2>=0); assert(n2<destNNodes);
	    destNNList[n1].insert(n2);
	  }
	}
      }
    }
  }

  element_containing_node.resize(destNNodes);
  for(int i=0;i<destNNodes;i++) 
    element_containing_node[i] = -1;
  
  std::vector<int> start(destNNodes, -1);

  for(int i=0; i<destNNodes; i++){
    if(element_containing_node[i]<0){
      // cerr<<"searching for node "<<i<<endl;

      // Get the ball rolling
      // cerr<<"Get the ball rolling\n";
      element_containing_node[i] = bruteForce(destX[i], destY[i]);
      
      // Form initial front
      // cerr<<"Form initial front\n";
      std::set<int> front;
      for(std::set<int>::iterator it=destNNList[i].begin(); it!=destNNList[i].end(); ++it){
	if(start[*it]<0){
	  front.insert(*it);
	  start[*it] = element_containing_node[i];
	}
      }

      while(!front.empty()){
	// cerr<<"size of new front = "<<front.size()<<endl;
	std::set<int> newFront;
	for(std::set<int>::const_iterator it=front.begin(); it!=front.end(); ++it){
	  // First try vacinity search
	  assert(start[*it]>=0);
	  // cerr<<"starting point for "<<*it<<" is "<<start[*it]<<endl;

	  std::pair<int, double> guess( vicinity(start[*it], destX[*it], destY[*it], 100) );
	  if(guess.second>=0.0){
	    element_containing_node[*it] = guess.first;
	  }else{
	    element_containing_node[*it] = bruteForce(destX[*it], destY[*it]);
	  }
	  // cerr<<"found "<<*it<<" in element "<<element_containing_node[*it]<<endl;
	  for(std::set<int>::iterator jt=destNNList[*it].begin(); jt!=destNNList[*it].end(); ++jt){
	    if(start[*jt]<0){
	      newFront.insert(*jt);
	      start[*jt] = element_containing_node[*it];
	    }
	  } 
	}
	front.swap(newFront);
      }
    }
  }

  // finally calculate all the correct shape functions so we're ready
  // to interpolate
  L012.resize(destNNodes*3);
  for(int i=0; i<destNNodes; i++)
    shapeFunction(element_containing_node[i], destX[i], destY[i], &(L012[i*3]));

  // Don't need the adjancy lists any more so lets delete them for memory
  sourceEEList.clear();
  destNNList.clear();
}

template<class R>
double InterpolationTri3Tri3<R>::area(int eid){
  return area(sourceENList[eid*3  ],
	      sourceENList[eid*3+1],
	      sourceENList[eid*3+2]);
}

template<class R>
void InterpolationTri3Tri3<R>::Interpolate(const R *source, R* dest){
  if(element_containing_node.empty())
    advancing_front_search();

  for(int i=0; i<destNNodes; i++){
    int n0 = sourceENList[element_containing_node[i]*3  ] - offset;
    int n1 = sourceENList[element_containing_node[i]*3+1] - offset;
    int n2 = sourceENList[element_containing_node[i]*3+2] - offset;

    dest[i] = 
      source[n0]*L012[i*3  ] + 
      source[n1]*L012[i*3+1] + 
      source[n2]*L012[i*3+2];
  }
}

template<class R>
double InterpolationTri3Tri3<R>::area(int n1, int n2, int n3){
  n1 = n1 - offset;
  n2 = n2 - offset;
  n3 = n3 - offset;
  
  return area(sourceX[n1], sourceY[n1],
	      sourceX[n2], sourceY[n2],
	      sourceX[n3], sourceY[n3]);
}

template<class R>
double InterpolationTri3Tri3<R>::area(int n1, int n2, R x, R y){
  n1 = n1 - offset;
  n2 = n2 - offset;
  
  return area(sourceX[n1], sourceY[n1],
	      sourceX[n2], sourceY[n2],
	      x, y);
}

template<class R>
double InterpolationTri3Tri3<R>::area(R x1, R y1,
				      R x2, R y2,
				      R x3, R y3){
  double vx1 = x2 - x1;
  double vy1 = y2 - y1;
  
  double vx2 = x3 - x1;
  double vy2 = y3 - y1;
  
  return 0.5*(vx1*vy2 - vx2*vy1);
}

template<class R>
int InterpolationTri3Tri3<R>::bruteForce(R x, R y){
  NBruteForce++;
  std::cerr<<"brute force "<<NBruteForce<<std::endl;

  std::pair<int, double> best( IO(0, x, y) );
  best.first = 0;
  
  if(best.second>=0.0) 
    return 0;
  
  for(int i=1; i<sourceNElements; i++){
    std::pair<int, double> guess( IO(i, x, y) );
    if(best.second<guess.second){
      best.first  = i;
      best.second = guess.second;

      if(best.second>=0.0) 
	return best.first;
    }
  }
  return best.first;
}

template<class R>
void InterpolationTri3Tri3<R>::SetDest(const int *enlist, const R *X, const R *Y,
				       int NNodes, int NElements){
  element_containing_node.clear();

  destNNodes    = NNodes;
  destNElements = NElements;

  destENList.resize(NElements*3);
  memcpy(&(destENList[0]), enlist, NElements*3*sizeof(int));
  
  destX.resize(NNodes);
  memcpy(&(destX[0]), X, NNodes*sizeof(R));

  destY.resize(NNodes);
  memcpy(&(destY[0]), Y, NNodes*sizeof(R));
}

template<class R>
void InterpolationTri3Tri3<R>::SetSource(const int *enlist, const R *X, const R *Y,
					 int NNodes, int NElements){
  element_containing_node.clear();

  sourceNNodes    = NNodes;
  sourceNElements = NElements;
  sourceENList.resize(NElements*3);
  memcpy(&(sourceENList[0]), enlist, NElements*3*sizeof(int));
  sourceX.resize(NNodes);
  memcpy(&(sourceX[0]), X, NNodes*sizeof(R));
  sourceY.resize(NNodes);
  memcpy(&(sourceY[0]), Y, NNodes*sizeof(R));

  areas.resize(NElements);
  for(int i=0; i<NElements; i++)
    areas[i] = area(i);
}

template<class R>
std::pair<int, double> InterpolationTri3Tri3<R>::IO(int eid, R x, R y){
  double L[3];
  shapeFunction(eid, x, y, L);

  std::pair<int, double> minL(0, L[0]);
  if(L[1]<minL.second){
    minL.first  = 1;
    minL.second = L[1];
  }
  if(L[2]<minL.second){
    minL.first  = 2;
    minL.second = L[2];
  }
  
  // send back the weight and element id
  if(minL.second>=0.0)
    minL.first = eid;
  else
    minL.first = sourceEEList[eid*3 + minL.first];

  // cerr<<minL.second<<endl;
  return minL;
}

template<class R>
void InterpolationTri3Tri3<R>::shapeFunction(int eid, R x, R y, double *L){
  L[0] = area(sourceENList[eid*3 + 1], sourceENList[eid*3 + 2], x, y)/areas[eid];
  L[1] = area(sourceENList[eid*3 + 2], sourceENList[eid*3 + 0], x, y)/areas[eid];
  L[2] = area(sourceENList[eid*3 + 0], sourceENList[eid*3 + 1], x, y)/areas[eid];
}

template<class R>
std::pair<int, double> InterpolationTri3Tri3<R>::vicinity(int start, R x, R y, int tries){
  std::pair<int, double> guess( IO(start, x, y) );
  if((guess.first<0)||(guess.second>=0.0)) // if dead-end or success then return
    return guess;
  
  for(int i=1;i<tries; i++){
    guess = IO(guess.first, x, y);
    if((guess.first<0)||(guess.second>=0.0)) // if dead-end or success then return
      return guess;
  }

  return guess;
}

#endif
