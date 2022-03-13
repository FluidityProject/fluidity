#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "H5hut.h"

#ifndef H5_HAVE_PARALLEL
#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
#endif
#endif

const char* FNAME = "simple_tet1.h5";

typedef struct vertex {
	h5_float64_t P[3];
} vertex_t;

typedef struct elem {
	h5_loc_idx_t vids[4];
} elem_t;

vertex_t Vertices[] = {
	{{-1.0,  0.0,  0.0}},
	{{ 0.0,  0.0,  1.0}},
	{{ 0.0,  1.0,  0.0}},
	{{ 1.0,  0.0,  0.0}},
};

elem_t Elems[] = {
	{{ 0, 1, 2, 3 }},
};

const int num_vertices = sizeof (Vertices) / sizeof (Vertices[0]);
const int num_elems = sizeof (Elems) / sizeof (Elems[0]);

int
main (
	int argc,
	char* argv[]
	) {
	/* abort program on errors in library */
	H5SetErrorHandler (H5AbortErrorhandler);
	H5SetVerbosityLevel (H5_DEBUG_ALL);

	/* open file and add mesh */
	h5_file_t const f = H5OpenFile (FNAME, H5_O_WRONLY, 0);
	h5t_mesh_t* mesh;
	H5FedAddTetrahedralMesh (f, "0", &mesh);

	/* store vertices */
	H5FedBeginStoreVertices (mesh, num_vertices);
	int i;
	for (i = 0; i < num_vertices; i++) {
		H5FedStoreVertex (mesh, -1, Vertices[i].P);
	}
	H5FedEndStoreVertices (mesh);

	/* store elements */
	H5FedBeginStoreElements (mesh, num_elems);
	for (i = 0; i < num_elems; i++) {
		H5FedStoreElement (mesh, Elems[i].vids);
	}
	H5FedEndStoreElements (mesh);

	/* add 1. Level */
	H5FedBeginRefineElements (mesh);
	H5FedRefineElement (mesh, 0);
	H5FedEndRefineElements (mesh);

	H5FedCloseMesh (mesh);
	H5CloseFile (f);
	return 0;
}
