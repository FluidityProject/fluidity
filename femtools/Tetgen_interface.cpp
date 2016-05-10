#ifdef HAVE_LIBTET
#include "tetgen.h"
#include <iostream>


extern "C" {
  void tetgen_binding(int nnondes, double* nodes, int nfacets, int* facets, int* neles, void** context);
  void set_from_tetgenio(int neles, int* element_list, void* context);
  void tetgen_cleanup(void* context);
    }



void tetgen_binding(int nnodes, double* nodes, int nfacets, int* facets, int* neles, void** context) {
  tetgenio in;
  tetgenio* out;

  
  std::cout<< "hello!"<<std::endl;
  out = new tetgenio;
  *context = (void*) out;
  tetgenio::facet *f;
  tetgenio::polygon *p;

  // Fortran indices start from 1
  in.firstnumber = 1;

  in.numberofpoints = nnodes;
  in.pointlist = nodes;

  in.numberoffacets=nfacets;

  std::cout<< "facetlist"<<std::endl;
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];

 std::cout<< "facets"<<std::endl;
  for (int i=0;i<nfacets;i++) {
    f = &in.facetlist[i];
    f->numberofpolygons=1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 3;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = facets[3*i+0];
    p->vertexlist[1] = facets[3*i+1];
    p->vertexlist[2] = facets[3*i+2];
    in.facetmarkerlist[i]=1;
  }

  in.numberofholes = 1;
  in.holelist = new double[3];
  in.holelist[0]=0.0;
  in.holelist[1]=0.0;
  in.holelist[2]=3.0;


  tetgenbehavior b;
  b.parse_commandline("pYS0CO2/1");

  std::cout<< "tetrahedralize"<<std::endl;
    
  tetrahedralize(&b, &in, out);

  *neles = out->numberoftetrahedra;

  in.numberofpoints = 0;
  in.pointlist = NULL;


  std::cout<< "facets"<<std::endl;
  for (int i=0;i<nfacets;i++) {
    f = &in.facetlist[i];
    f->numberofpolygons=0;
    p = &f->polygonlist[0];
    p->numberofvertices = 0;
    delete p->vertexlist;
    p->vertexlist=NULL;
    delete f->polygonlist;
    f->polygonlist=NULL;
  }

  in.numberofholes=0;
  delete in.holelist;
  in.holelist=NULL;

  delete in.facetlist;
  in.facetlist=NULL;
  delete in.facetmarkerlist;
  in.facetmarkerlist=NULL;

  std::cout<< "out"<<std::endl;
}

void set_from_tetgenio(int neles, int* element_list, void* context){
  tetgenio* out;
  out= (tetgenio*) context;
  
  for(int i=0; i<4*neles; i++){
    element_list[i] = out->tetrahedronlist[i];
  }
}

void tetgen_cleanup(void* context)
{
  tetgenio* out;
  out= (tetgenio*) context;

  std::cout<< "cleanup"<<std::endl;

  delete out;
}

#endif
