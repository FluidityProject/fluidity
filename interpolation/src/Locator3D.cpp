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
#include "Locator3D.h"

using namespace std;

Locator3D::Locator3D(){
  verbose_off();
}

size_t Locator3D::bruteforceSearch(const deque< Tetra4<gigreal_t> >& elements, 
           const vector<gigreal_t> x){
  if(verbose)
    cout<<"size_t Locator3D::bruteforceSearch(const deque< Tetra4<gigreal_t> >& elements, const vector<gigreal_t> x)\n";
  
  // initalise best_XX
  assert(!elements.empty());
  gigreal_t best_IO = elements[0].IO(x);
  size_t best_eid = 0;
  
  // Search through elements
  size_t eid = 0;
  for(deque< Tetra4<gigreal_t> >::const_iterator j=elements.begin(); j!=elements.end(); ++j){
    gigreal_t IO = (*j).IO(x);
    
    if(IO>=0.0){
      // eid is ok
      return eid;
    }else{
      if(IO>best_IO){
        best_IO   = IO;
        best_eid  = eid;
      }
    }
    eid++;      
  }
  
  return best_eid;
}

int Locator3D::get_node2element(int nid) const{
  return node2element[nid];
}

int Locator3D::GiD_read_msh(const string& fileName, 
          deque< Tetra4<gigreal_t> >& elementList, 
          deque< Node<gigreal_t> >& nodeList){
  if(verbose)
    cout<<"int Locator3D::GiD_read_msh(const string& fileName, deque< Tetra4<gigreal_t> >& elementList, deque< Node<gigreal_t> >& nodeList)\n";
  
  ifstream GiDFile( fileName.c_str() );
  string line;
  
  if (! GiDFile.is_open() ){
    cerr << "Error opening file";
    return(-1); 
  }

  // Get the header
  getline(GiDFile, line);
  vector<string> tokens;
  Tokenize(line, tokens);
  assert( tokens.size() );
  assert(tokens[0] == "MESH");
  assert(tokens[1] == "dimension");
  int ndim;
  stringstream(tokens[2]) >> ndim;
  assert(ndim == 3);
  assert(tokens[3] == "ElemType");
  string ElemType = tokens[4];
  assert(tokens[5] == "Nnode");
  int Nnode;
  stringstream(tokens[6]) >> Nnode;
  assert(Nnode == 4);
  
  // Read nodes
  getline(GiDFile, line);
  assert(line == "Coordinates");
  size_t gnn=0;
  for(;;){
    getline(GiDFile, line);
    Tokenize(line, tokens);
    
    assert( !tokens.empty() );
    if(tokens[0] == "end"){
      break;
    }else{
      Node<gigreal_t> nod;
    gigreal_t x,y,z;
    stringstream(tokens[1]) >> x;
    stringstream(tokens[2]) >> y;
    stringstream(tokens[3]) >> z;
    
      nod.set_coord(x, y, z);
      nod.set_gnn(gnn);
      nodeList.push_back( nod );
      gnn++;
    }
  }
  
  getline(GiDFile, line);  

  // Read elements
  getline(GiDFile, line);  
  assert(line == "Elements");
  for(;;){
    getline(GiDFile, line);
    Tokenize(line, tokens);
    
    if(tokens[0] == "end") break;
    
    assert(tokens.size() == 5);
    Tetra4<gigreal_t> tet;
    vector<size_t> nods(4);
  for(size_t j=0;j<4;j++){
    stringstream(tokens[j+1]) >> nods[j];
    nods[j]--;
  }
    tet.set_nodeList(nods);

    elementList.push_back(tet);
  }

  GiDFile.close();

  return 0;
}

int Locator3D::interpolate_linear(gigreal_t *RMEM1,  int *FIELDS1, int NFIELDS,
                  gigreal_t *RMEM2,  int *FIELDS2){
  if(verbose)
  cout<<"int Locator3D::interpolate_linear(gigreal_t *RMEM1,  int *FIELDS1, int NFIELDS, gigreal_t *RMEM2,  int *FIELDS2)\n";
  
  for(int i=0; i<nodes_1.size(); i++){
  for(int j=0; j<NFIELDS; j++){
    nodes_1[i].add_field(RMEM1[FIELDS1[j]-1+i] );
  }
  }
  
  interpolate_linear();
  
  for(int i=0; i<nodes_2.size(); i++){
  const vector<gigreal_t>& flds = nodes_2[i].get_fields();
  assert(flds.size()==(size_t)NFIELDS);
  for(int j=0; j<NFIELDS; j++){
    RMEM2[FIELDS2[j]-1+i] = flds[j];
  }
  }
  
  return 0;
}

int Locator3D::interpolate_linear(){
  if(verbose)
  cout<<"int Locator3D::interpolate_linear()\n";

  for(int i=0; i<nodes_2.size(); i++){
  int eid = node2element[i];
  Node<gigreal_t>& node = nodes_2[i];
  const vector<gigreal_t>& x = node.get_coord();

  // interpolation
  if(verbose)
    cout<<eid<<" -- "<<elements_1.size()<<endl;
  
  assert(eid < elements_1.size());
  const vector<size_t>& nodes = elements_1[eid].get_nodeList();
  const vector<gigreal_t>& f0 = nodes_1[ nodes[0] ].get_fields();
  const vector<gigreal_t>& f1 = nodes_1[ nodes[1] ].get_fields();
  const vector<gigreal_t>& f2 = nodes_1[ nodes[2] ].get_fields();
  const vector<gigreal_t>& f3 = nodes_1[ nodes[3] ].get_fields();

  if(verbose)
    cout<<" -- "<<nodes[0]<<" "<<nodes[1]<<" "<<nodes[2]<<" "<<nodes[3]<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
  
  vector<gigreal_t> f;
  elements_1[eid].interpolate(f0, f1, f2, f3, x, f);
  node.set_fields(f);
  }
  
  return 0;
}

int Locator3D::mkEElist(const deque< Tetra4<gigreal_t> >& elements,
            deque< vector<int> >& EElist){
  if(verbose)
  cout<<"int Locator3D::mkEElist(const deque< Tetra4<gigreal_t> >& elements, deque< vector<int> >& EElist)\n";
  
  // First step: create a node->element adjancy list
  map<size_t, set<size_t> > nodeE;
  size_t eid=0;
  for(deque< Tetra4<gigreal_t> >::const_iterator e=elements.begin(); e!=elements.end(); ++e){
    const vector<size_t>& nodes=(*e).get_nodeList();
    for(vector<size_t>::const_iterator n=nodes.begin(); n!=nodes.end(); ++n){
      nodeE[ *n ].insert( eid );
    }
    eid++;
  }
  
  // The EE list must have the following property. The i'th element in
  // the adjancy list must share the face opposit the i'th node of the
  // element.
  size_t ecnt = elements.size();
  EElist.clear();
  EElist.resize( ecnt );
  for(size_t i=0; i<ecnt; i++){
    const vector<size_t>& nodes=elements[i].get_nodeList();
    assert(nodes.size() == 4);
    
    EElist[i].resize(4);
    for(size_t j=0; j<4; j++){
      
      // Make the correct facet.
      size_t facet[3];
      if     (j==0){facet[0]=nodes[1];  facet[1]=nodes[2];  facet[2]=nodes[3];}
      else if(j==1){facet[0]=nodes[0];  facet[1]=nodes[2];  facet[2]=nodes[3];}
      else if(j==2){facet[0]=nodes[0];  facet[1]=nodes[1];  facet[2]=nodes[3];}
      else         {facet[0]=nodes[0];  facet[1]=nodes[1];  facet[2]=nodes[2];}
      
      int partner=-1;
      for(set<size_t>::const_iterator it=nodeE[ facet[0] ].begin(); it!=nodeE[ facet[0] ].end(); ++it){
  if(*it == i)
    continue;
  
  if( (nodeE[ facet[1] ].find( *it )!=nodeE[ facet[1] ].end())&&
      (nodeE[ facet[2] ].find( *it )!=nodeE[ facet[2] ].end()) ){
    partner = *it;
    break;
  }
      }
      // partner==-1 => boundary element
      // cout<<"partner - "<<partner<<endl;
      EElist[i][j] = partner;
    }
    
  }

  return 0;
}

int Locator3D::t4_to_t4_Search(int NNode1,  int NElems1, const gigreal_t X1[], const gigreal_t Y1[], const gigreal_t Z1[], const int ENList1[],
             int NNode2,  int NElems2, const gigreal_t X2[], const gigreal_t Y2[], const gigreal_t Z2[], const int ENList2[]){
  if(verbose)
    cout<<"int Locator3D::t4_to_t4_Search(int NNode1,  int NElems1, const gigreal_t X1[], const gigreal_t Y1[], const gigreal_t Z1[], const int ENList1[], "
  <<"int NNode2,  int NElems2, const gigreal_t X2[], const gigreal_t Y2[], const gigreal_t Z2[], const int ENList2[])\n";
  
  // Put stuff into internal data-structures. This really should go.
  nodes_1.resize(NNode1);
  for(int i=0; i<NNode1; i++){
    nodes_1[i].set_gnn(i);
    nodes_1[i].set_coord(X1[i], Y1[i], Z1[i]);
  }
  
  elements_1.resize(NElems1);
  for(int i=0; i<NElems1; i++){
    vector<size_t> nids(4);
    for(size_t j=0; j<4; j++){
      nids[j] = ENList1[4*i+j]-1;
    }
    elements_1[i].set_nodeList(nids);
  }
  
  nodes_2.resize(NNode2);
  for(int i=0; i<NNode2; i++){
    nodes_2[i].set_gnn(i);
    nodes_2[i].set_coord(X2[i], Y2[i], Z2[i]);
  }
  
  elements_2.resize(NElems2);
  for(int i=0; i<NElems2; i++){
    vector<size_t> nids(4);
    for(int j=0; j<4; j++){
      nids[j] = ENList2[4*i+j]-1;
    }
    elements_2[i].set_nodeList(nids);
  }
  
  t4_to_t4_Search();
  
  return(0);
}

int Locator3D::t4_to_t4_Search(){
  if(verbose)
  cout<<"int Locator3D::t4_to_t4_Search()\n";

  //
  // Fast grid-2-grid interpolation
  //
  node2element.resize(nodes_2.size());

  // Calculate shape functions.
  for(deque< Tetra4<gigreal_t> >::iterator j=elements_1.begin(); j!=elements_1.end(); ++j){
    const vector<size_t>& nods = (*j).get_nodeList();
    assert(nods.size() == 4);
    
  if(verbose)
    cout<<nods[0]<<" "<<nods[1]<<" "<<nods[2]<<" "<<nods[3]<<endl;

    int err = j->calcShapeFxn(nodes_1[ nods[0] ].get_coord(), nodes_1[ nods[1] ].get_coord(),
                nodes_1[ nods[2] ].get_coord(), nodes_1[ nods[3] ].get_coord());

    if(err){
      exit(-1);
    }
  }
  
  // Make element-element adjancy list.
  deque< vector<int> > EElist;
  mkEElist(elements_1, EElist);
  
  // Make node-node adjancy list.
  deque< vector<int> > NNlist;
  mkNNlist(elements_2, NNlist);
  
  size_t bcnt=0;
  for(;;){
    // The advancing front
    deque< size_t > frontPoints;
    
    {// Initalise advancing front list
      size_t nid=0;
      for(deque< Node<gigreal_t> >::iterator i=nodes_2.begin(); i!=nodes_2.end(); ++i){
    const unsigned char flag = (*i).get_flags();
    
    // Get the first point in the array that hasn't already been
    // done
    if( !(flag&0x2) ){
      frontPoints.push_back( nid );
      break;
    }
    nid++;
      }
    }
  
    if(frontPoints.empty()) // Finished :-))
      break;

  if(verbose)
    cout<<"Starting new front!\n";
    
    while(!frontPoints.empty()){
      Node<gigreal_t>& node = nodes_2[ frontPoints[0] ];
      const vector<gigreal_t>& x = node.get_coord();
      
      // Find the element that contains this node or the nearest element
      // to this node.
      size_t eid=0;
      size_t startElement = node.get_start_element();
      int findResult=vicinitySearch(elements_1, EElist, x, startElement, 10, eid);
      if(findResult<0){
    if(verbose)
      cout<<"Trying brute force.\n";
    // Vicinity search failed :-< Try brute force.
    eid = bruteforceSearch(elements_1, x);
    bcnt++;
      }
    
    node2element[node.get_gnn()] = eid;
      
      // Flag as interpolated
      node.set_flags( node.get_flags()|0x2 );
    
      frontPoints.pop_front();
    
      size_t gnn = node.get_gnn();
      assert(gnn<NNlist.size());
      for(vector<int>::const_iterator in=NNlist[gnn].begin(); in!=NNlist[gnn].end(); ++in){
    const unsigned char flag = nodes_2[*in].get_flags();
    if( !(flag&0x2) ){ // non-interpolated point
      if( !(flag&0x4) ){ // not already in front list
      frontPoints.push_back( *in );
      nodes_2[*in].set_start_element(eid); 
      nodes_2[*in].set_flags( flag|0x4 );
      }
    }
      }
    }
  } 
  
  return(0);
}

int Locator3D::test_gid(string file1, string file2){
  if(verbose)
  cout<<"int Locator3D::test_gid(string file1, string file2)\n";
  
  int err = GiD_read_msh(file1, elements_1, nodes_1);
  if(err){
    cerr<<"GiD_read_msh() failed!!"<<endl;
    exit(err);
  }
  
  err = GiD_read_msh(file2, elements_2, nodes_2);
  if(err){
    cerr<<"GiD_read_msh() failed!!"<<endl;
    exit(err);
  }
  
  // Choose an analitical function and use it to calculate field
  // values on mesh 1
  for(deque< Node<gigreal_t> >::iterator i=nodes_1.begin(); i!=nodes_1.end(); ++i){
    gigreal_t x = (*i).get_x();
    gigreal_t y = (*i).get_y();
    gigreal_t z = (*i).get_z();
    // gigreal_t f = x*x + y*y + z*z;
  gigreal_t f = x + 2.0*y + 3.0*z;
    (*i).set_fields(&f, 1);
  }
  
  t4_to_t4_Search();
  
  interpolate_linear();
  
  // Calculate the interpolation error
  gigreal_t e = 0.0;
  int cnt = 0;
  for(deque< Node<gigreal_t> >::iterator i=nodes_2.begin(); i!=nodes_2.end(); ++i){
    gigreal_t x = (*i).get_x();
    gigreal_t y = (*i).get_y();
    gigreal_t z = (*i).get_z();
    // gigreal_t f = x*x + y*y + z*z;
  gigreal_t f = x + 2.0*y + 3.0*z;

    const vector<gigreal_t>& fld = (*i).get_fields();

    cnt++;
    e += pow(f - fld[0], 2)/f;
  }
  gigreal_t Chi2 = e/cnt;
  cerr << "Chi2 = " << Chi2 << endl;

  return 0;
}

int Locator3D::mkNNlist(const deque< Tetra4<gigreal_t> >& elements,
            deque< vector<int> >& NNlist){
  if(verbose)
  cout<<"int Locator3D::mkNNlist(const deque< Tetra4<gigreal_t> >& elements, deque< vector<int> >& NNlist)\n";
  
  // First step: create a node->node using the map and set template
  map<size_t, set<size_t> > nodeN;
  for(deque< Tetra4<gigreal_t> >::const_iterator e=elements.begin(); e!=elements.end(); ++e){
    const vector<size_t>& nodes=(*e).get_nodeList();
    for(vector<size_t>::const_iterator n1=nodes.begin(); n1!=nodes.end(); ++n1){
      for(vector<size_t>::const_iterator n2=nodes.begin(); n2!=nodes.end(); ++n2){
  nodeN[ (*n1) ].insert( *n2 );
      }
    }
  }
  
  // Copy this into the real thing
  NNlist.resize( nodeN.size() );
  size_t ncnt=0;
  for(map<size_t, set<size_t> >::const_iterator it=nodeN.begin(); it!=nodeN.end(); ++it){
    assert((*it).first == ncnt);
    for(set<size_t>::const_iterator in=(*it).second.begin(); in!=(*it).second.end(); ++in){
      NNlist[ncnt].push_back( *in );
    }
    ncnt++;
  }
  
  return 0;
}

void Locator3D::Tokenize(const string& str, vector<string>& tokens){
  string delimiters = " ";
  
  tokens.clear();
  
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos){
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);

    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }

  return;
}

void Locator3D::verbose_off(){
  verbose = false;
}

void Locator3D::verbose_on(){
  verbose = true;
}

// vicinitySearch does a neighbour-to-neighbour vicinity search. 
// Return values: 
// 0 => found element 
//-1 => reached edge boundary/dead-end
//-2 => didn't find element with NTries


// The EE list must have the following property. The i'th element in
// the adjancy list must share the face opposit the i'th node of the
// element.
int Locator3D::vicinitySearch(const deque< Tetra4<gigreal_t> >& elements,
                const deque< vector<int> >& EElist,
                const vector<gigreal_t>& x,
                const size_t start, 
                const int NTries,
                size_t& last){
  if(verbose)
  cout<<"int Locator3D::vicinitySearch(const deque< Tetra4<gigreal_t> >& elements, const deque< vector<int> >& EElist, const vector<gigreal_t>& x, "
    <<"const size_t start, const int NTries, size_t& last)\n";
  
  last = start;
  for(size_t i=0; i<NTries; i++){
    size_t node;
    
    gigreal_t IO = elements[ last ].IO(x, node);
    
    // Successful
    if(IO>=0.0)
      return(0);
    
    // Check for dead-end
    if(EElist[last][node] < 0)
      return(-1);
    
    last = EElist[last][node];
  }

  return(-2);
}

// Fortran interface
extern "C" {  
  // 1 is the source mesh, 2 is the target mesh
#define fltetra4totetra4_fc F77_FUNC(fltetra4totetra4, FLTETRA4TOTETRA4)
  void fltetra4totetra4_fc(int *NNode1,  int *NElems1, gigreal_t *X1, gigreal_t *Y1, gigreal_t *Z1, 
         int *ENList1, gigreal_t *RMEM1,  int *FIELDS1, int *NFIELDS, 
         int *NNode2,  int *NElems2, gigreal_t *X2, gigreal_t *Y2, gigreal_t *Z2, 
         int *ENList2, gigreal_t *RMEM2,  int *FIELDS2, int *IERROR){
    // reset error
    *IERROR = 0;
    
    Locator3D mapper;
    mapper.t4_to_t4_Search(*NNode1,  *NElems1, X1, Y1, Z1, ENList1,
         *NNode2,  *NElems2, X2, Y2, Z2, ENList2);
    
    mapper.interpolate_linear(RMEM1,  FIELDS1, *NFIELDS, RMEM2, FIELDS2);
    
    return;
  }
  
#define fl_node2tetra_locator_fc F77_FUNC(fl_node2tetra_locator, FL_NODE2TETRA_LOCATOR)
  void fl_node2tetra_locator_fc(const int *NNode1, const int *NElems1, const gigreal_t *X1, const gigreal_t *Y1, const gigreal_t *Z1, 
        const int *ENList1,
        const int *NNode2, const int *NElems2, const gigreal_t *X2, const gigreal_t *Y2, const gigreal_t *Z2, 
        const int *ENList2, int *node2element, int *IERROR){
    // reset error
    *IERROR = 0;
    
    Locator3D mapper;
    mapper.t4_to_t4_Search(*NNode1,  *NElems1, X1, Y1, Z1, ENList1,
         *NNode2,  *NElems2, X2, Y2, Z2, ENList2);
    
    for(int i=0;i<(*NNode2);i++)
      node2element[i] = mapper.get_node2element(i);
    
    return;
  }
}
