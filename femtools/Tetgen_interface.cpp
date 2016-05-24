#include "confdefs.h"
#include <iostream>
#include <stdlib.h>


extern "C" {
  struct mesh_data {
    int nnodes, nelements, nfacets, nholes;
    int lnode_ids, lregion_ids, lface_ids;
    double* nodes;
    int* node_ids;
    int* ndglno;
    int* region_ids;
    int* facets;
    int* face_ids;
    double* holes;
  };
void mesh_data_cleanup(struct mesh_data* mesh, int* dim);
}

void mesh_data_cleanup(struct mesh_data* mesh, int* dim){
  switch(*dim){
  case(2):{
    // triangle uses malloc
    free(mesh->nodes);
    mesh->nodes = NULL;
    free(mesh->ndglno);
    mesh->ndglno = NULL;
    free(mesh->facets);
    mesh->facets = NULL;
    break;
  }
  case(3):{
    // tetgen uses new
  if (mesh->nodes) {
    delete mesh->nodes;
    mesh->nodes = NULL;
  }
  if (mesh->ndglno) {
    delete mesh->ndglno;
    mesh->ndglno = NULL;
  }
  if (mesh->facets) {
    delete mesh->facets;
    mesh->facets = NULL;
  }
  break;
  }
  }
}

#ifdef HAVE_LIBTET
#include "tetgen.h"

extern "C" {
  struct mesh_data tetgen(struct mesh_data* input_mesh,char* command);
  void* mesh_to_tetgenio(struct mesh_data* mesh);
  struct mesh_data tetgenio_to_mesh(void* out);
    }

void* mesh_to_tetgenio(struct mesh_data* mesh) {
  tetgenio* out;
  out = new tetgenio;
  tetgenio::facet *f;
  tetgenio::polygon *p;

   // Fortran indices start from 1
  out->firstnumber = 1;

  out->numberofpoints = mesh->nnodes;
  out->pointlist = mesh->nodes;

  out->numberoffacets=mesh->nfacets;
  out->facetlist = new tetgenio::facet[out->numberoffacets];
  out->facetmarkerlist = new int[out->numberoffacets];

  for (int i=0;i<mesh->nfacets;i++) {
    f = &out->facetlist[i];
    f->numberofpolygons=1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 3;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = ((int*) mesh->facets)[3*i+0];
    p->vertexlist[1] = ((int*) mesh->facets)[3*i+1];
    p->vertexlist[2] = ((int*) mesh->facets)[3*i+2];
    if (!mesh->lface_ids) {
      out->facetmarkerlist[i]=i;
    }
  }

  out->numberoftetrahedra = mesh->nelements;
  out->tetrahedronlist = mesh->ndglno;

  
  out->numberofholes = mesh->nholes;
  out->holelist=mesh->holes;

  return (void*) out;
}

struct mesh_data tetgenio_to_mesh(void* context){
  tetgenio* in = (tetgenio*) context;
  struct mesh_data out;

  out.nnodes = in->numberofpoints;
  out.nodes = in->pointlist;
  out.nfacets = in->numberoffacets;
  out.facets = new int[3*in->numberoffacets];
  for (int i=0;i<out.nfacets;i++) {
    for(int j=0;j<3;j++) {
      out.facets[3*i+j] = in->facetlist[i].polygonlist[0].vertexlist[j];
    }
  }
  out.nelements = in->numberoftetrahedra;
  out.ndglno = in->tetrahedronlist;

  //  We ignore the holes
  //  out.nholes=in->numberofholes;
  //  out.holes=in->holelist;

  return out;
}
 
struct mesh_data tetgen(struct mesh_data* input_mesh,char* command){
  tetgenio input = *(tetgenio*) mesh_to_tetgenio(input_mesh);
  tetgenio output;
  struct mesh_data output_mesh;

  tetgenbehavior b;
  b.parse_commandline(command);

  tetrahedralize(&b, &input, &output);
  
  output_mesh = tetgenio_to_mesh(&output);

  // prevent the tetgenio destructor from wiping out the data
  // now held in our mesh struct

  input.numberofpoints = 0;
  input.pointlist = NULL;
  input.numberoftetrahedra = 0;
  input.tetrahedronlist = NULL;
  input.numberofholes = 0;
  input.holelist = NULL;

  output.numberofpoints = 0;
  output.pointlist = NULL;
  output.numberoftetrahedra = 0;
  output.tetrahedronlist = NULL;

  return output_mesh;
}

#endif

#ifdef HAVE_LIBTRIANGLE
#define VOID void
#define ANSI_DECLARATORS 1
extern "C" {
  #include "triangle.h"
}

extern "C" {
  struct mesh_data triangle(struct mesh_data* input_mesh,char* command);
   void* mesh_to_triangleio(struct mesh_data* mesh);
  struct mesh_data triangleio_to_mesh(void* out);
    }

void* mesh_to_triangleio(struct mesh_data* mesh) {
  struct triangulateio* out;
  out = new struct triangulateio;

   // Fortran indices start from 1

  out->numberofpoints = mesh->nnodes;
  out->pointlist = mesh->nodes;

  out->numberofsegments=mesh->nfacets;
  out->segmentlist = mesh->facets;
  out->segmentmarkerlist = NULL;

  out->numberoftriangles = mesh->nelements;
  out->trianglelist = mesh->ndglno;
  out->numberofcorners=3;

  out->numberofholes = mesh->nholes;
  out->holelist=mesh->holes;

  return (void*) out;
}

struct mesh_data triangleio_to_mesh(void* context){
  struct triangulateio* in = (struct triangulateio*) context;
  struct mesh_data out;

  out.nnodes = in->numberofpoints;
  out.nodes = in->pointlist;
  out.nfacets = in->numberofsegments;
  out.facets = in->segmentlist;
  out.nelements = in->numberoftriangles;
  out.ndglno = in->trianglelist;

  //  We ignore the holes
  //  out.nholes=in->numberofholes;
  //  out.holes=in->holelist;

  return out;
}

struct mesh_data triangle(struct mesh_data* input_mesh,char* command){
  struct triangulateio input = *(struct triangulateio*) mesh_to_tetgenio(input_mesh);
  struct triangulateio output, varout;
  struct mesh_data output_mesh;

  triangulate(command, &input, &output, &varout);
  
  output_mesh = triangleio_to_mesh(&output);

 // triangle makes us clean up our own mess, so no clearup heer (I think)

  return output_mesh;
}

#endif
