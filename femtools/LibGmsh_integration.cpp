/*  Copyright (C) 2006 Imperial College London and others.
    
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
#include "fmangle.h"

#ifdef HAVE_LIBGMSH
#include "gmsh/Gmsh.h"
#include "gmsh/GModel.h"
#include "gmsh/GEntity.h"
#include "gmsh/MVertex.h"
#include "gmsh/Context.h"
#include "gmsh/PView.h"
#include "gmsh/discreteVertex.h"
#include "gmsh/discreteEdge.h"
#include "gmsh/discreteFace.h"
#include "gmsh/discreteRegion.h"

#include <vector>
#include <map>
#include <iterator>

typedef std::vector<GEntity*> EntityVector;
typedef std::vector<MVertex*> VertexVector;

int dim;

int FluidityNodeOrdering(int type, int node) {
  int quads[]= {0,1,3,2};
  int hexes[]= {0,1,3,2,4,5,7,6};
  switch(type) {
  case TYPE_QUA:
    return quads[node];
  case TYPE_HEX:
    return hexes[node];
  }
  return node;
}

extern "C" {

  void cgmsh_initialise() {
    GmshInitialize();
  }

  void cgmsh_finalise(GModel *&gm) {
    delete gm;
    GmshFinalize();
  }

  void cread_gmsh_file(GModel *&gm, const char* filename) {
      GmshInitialize();
      gm = new GModel;
      CTX::instance()->mesh.switchElementTags=1;
      gm->readMSH(filename);
      PView::readMSH(filename);
  }

  void cread_gmsh_sizes(GModel *&gm, int &numNodes, int &numFaces,
			int &numElements, bool &haveRegionIDs,
			bool &haveBounds, bool &haveElementOwners,
			int &gdim, int &loc, int &sloc) {

    numNodes = gm->getNumMeshVertices();
    
    dim = 0;
    if (gm->getNumRegions()>0) {
      dim = 3;
    } else if (gm->getNumFaces()>0) {
      dim = 2;
    } else if (gm->getNumEdges()>0) {
      dim = 1;
    }
    gdim = dim;

    numElements = 0;
    numFaces = 0;
    haveRegionIDs = false;
    haveBounds = false;
    haveElementOwners = false;
    std::vector<GEntity*> eles;
    std::vector<GEntity*> faces;
    gm->getEntities(eles, dim);
    gm->getEntities(faces, dim-1);
    if (eles.size()>0) {
      loc = eles[0]->getMeshElement(0)->getNumVertices();
    } else {
      loc = 0;
    }
    if (faces.size()>0) {
      sloc = faces[0]->getMeshElement(0)->getNumVertices();
    } else {
      sloc = 0;
    }
    for (EntityVector::iterator it = eles.begin(); it != eles.end(); ++it) {
      numElements = numElements + (*it)->getNumMeshElements();
      haveRegionIDs = haveRegionIDs || (*it)->tag()>0;
    }
    for (EntityVector::iterator it = faces.begin(); it != faces.end(); ++it) {
      numFaces = numFaces + (*it)->getNumMeshElements();
      haveBounds = haveBounds || (*it)->tag()>0;
      haveElementOwners = haveElementOwners || (*it)->getPhysicalEntities().size()==3;
    }
  }

  void cread_gmsh_element_connectivity(void *&gmv, int &numElements,
				       int &loc, int *ndglno, int *regionIDs) {

    GModel *gm = (GModel*) gmv; 
    int k=0;
    std::vector<GEntity*> ents;
    gm->getEntities(ents, dim);
    for (EntityVector::iterator ent = ents.begin(); ent != ents.end(); ++ent) {
      for (unsigned int i=0; i< (*ent)->getNumMeshElements(); ++i) {
	MElement* e = (*ent)->getMeshElement(i);
	for (int j=0; j<e->getNumVertices(); ++j) {
	  MVertex* v = e->getVertex(FluidityNodeOrdering(e->getType(),j));
	  ndglno[loc*k+j] = v->getIndex();
	}
	if ((*ent)->tag()>0) regionIDs[k] = (*ent)->tag();
	++k;
      }
    }
  }

  void cread_gmsh_points(GModel *&gm, int &dim, int &numNodes, double *val) {
    for (int i=0; i<numNodes; ++i) {
      MVertex *v = gm->getMeshVertexByTag(i+1);
      switch(dim) {
      case 3:
	val[dim*i+2] = v->z(); 
      case 2:
        val[dim*i+1] = v->y(); 
      case 1:
	val[dim*i+0] = v->x(); 
      }
    }
  }

  void cread_gmsh_face_connectivity(GModel *&gm, int *numFaces,
				    int &sloc, int *sndglno,
				    int *boundaryIDs, int *faceOwner) {
    int k=0;
    std::vector<GEntity*> ents;
    gm->getEntities(ents, dim-1);
    for (EntityVector::iterator ent = ents.begin(); ent != ents.end(); ++ent) {
      for (unsigned int i=0; i< (*ent)->getNumMeshElements(); ++i) {
	MElement* e = (*ent)->getMeshElement(i);
	for (int j=0; j<e->getNumVertices(); ++j) {
	  MVertex* v = e->getVertex(FluidityNodeOrdering(e->getType(),j));
	  sndglno[sloc*k+j] = v->getIndex();
	}
	if ((*ent)->tag()>0) boundaryIDs[k] = (*ent)->tag();
	if ((*ent)->getPhysicalEntities().size()==3)
	  faceOwner[k] = (*ent)->getPhysicalEntities()[2];
	++k;
      }
    }
  }

  void cmesh_to_gmodel(GModel *&gm, int &binary,
		       int &numNodes, int &numElements, int &numFaces,
		       int &loc, int &sloc,
		       int &gdim, double* val, 
		       int &etype, int* eles,
		       int &ftype, int* faces,
		       int *ele_ids, int *face_ids) {
    
    GmshInitialize();
    gm = new GModel;

    GEntity *e;
    std::map<int,GVertex*> point;
    std::map<int,GEdge*> edge;
    std::map<int,GFace*> face;
    std::map<int,GRegion*> region;
    
    std::set<int> uele_ids;
    std::set<int> uface_ids;
    
    uele_ids.insert(0);
    if (ele_ids) {
      for (int i=0;i<numElements;++i) {
	uele_ids.insert(ele_ids[i]);
      }
    }
    uface_ids.insert(0); 
    if(face_ids) {
      for (int i=0;i<numFaces;++i) {
	uface_ids.insert(face_ids[i]);
      }
    }

    switch(gdim) {
    case 1:
      for (std::set<int>::iterator it=uele_ids.begin();
	   it != uele_ids.end(); ++it) {
	edge[*it] = new discreteEdge(gm,0,NULL,NULL);
	edge[*it]->addPhysicalEntity(*it);
	edge[*it]->setTag(*it);
	gm->add(edge[*it]);
      }
      for (std::set<int>::iterator it=uface_ids.begin();
	   it != uface_ids.end(); ++it) {
	point[*it] = new discreteVertex(gm,0);
	point[*it]->addPhysicalEntity(*it);
	point[*it]->setTag(*it);
	gm->add(point[*it]);
      }
      e = edge[0];
      break;
    case 2:
      for (std::set<int>::iterator it=uele_ids.begin();
	   it != uele_ids.end(); ++it) {
	face[*it] =  new discreteFace(gm,0);
	face[*it]->addPhysicalEntity(*it);
	face[*it]->setTag(*it);
	gm->add(face[*it]);
      }
      for (std::set<int>::iterator it=uface_ids.begin();
	   it != uface_ids.end(); ++it) {
	edge[*it] = new discreteEdge(gm,0,NULL,NULL);
	edge[*it]->addPhysicalEntity(*it);
	edge[*it]->setTag(*it);
	gm->add(edge[*it]);
      }
      e=face[0];
      break;
    case 3:
      for (std::set<int>::iterator it=uele_ids.begin();
	   it != uele_ids.end(); ++it) {
	region[*it] = new discreteRegion(gm,0);
	region[*it]->addPhysicalEntity(*it);
	region[*it]->setTag(*it);
	gm->add(region[*it]);
      }
      for (std::set<int>::iterator it=uface_ids.begin();
	   it != uface_ids.end(); ++it) {
	face[*it] = new discreteFace(gm,0);
	face[*it]->addPhysicalEntity(*it);
	face[*it]->setTag(*it);
	gm->add(face[*it]);
      }
      e=region[0];
      break;
    }

    MElementFactory f;

    double x[3] = {0,0,0};
    for (int n=0;n<numNodes;++n) {
      for (int i=0;i<gdim;++i) {
	x[i] = val[gdim*n+i];
      }
      MVertex* v =  new MVertex(x[0],x[1],x[2],NULL,n+1);
      v->setIndex(n+1);
      e->addMeshVertex(v);
    }

    for (int n=0;n<numElements;++n) {
      int id = 0;
      if (ele_ids) id = ele_ids[n];

      VertexVector vertices;
      for (int i=0;i<loc;++i) {
	vertices.push_back(e->getMeshVertex(eles[loc*n+i]-1));
      }
      MElement *ele = f.create(etype, vertices);
      switch(gdim) {
      case 1:
	edge[id]->addLine((MLine*) ele);
	break;
      case 2:
	face[id]->addTriangle((MTriangle*) ele);
	break;
      case 3:
	region[id]->addTetrahedron((MTetrahedron*) ele);
	break;
      }
    }

    for (int n=0;n<numFaces;++n) {
      int id = 0;
      if (face_ids) id = face_ids[n];

      VertexVector vertices;
      for (int i=0;i<sloc;++i) {
	vertices.push_back(e->getMeshVertex(faces[sloc*n+i]-1));
      }
      MElement *ele = f.create(ftype, vertices);
      switch(gdim) {
      case 1:
	point[id]->addPoint((MPoint*) ele);
	break;
      case 2:
	edge[id]->addLine((MLine*) ele);
	break;
      case 3:
	face[id]->addTriangle((MTriangle*) ele);
	break;
      }
    }
    
  }

  void cwrite_gmsh_file(GModel *&gm, char* filename) {
    gm->save(filename);
  }

}

#endif
