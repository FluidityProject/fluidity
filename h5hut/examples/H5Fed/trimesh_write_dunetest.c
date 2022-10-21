#include <stdio.h>
#include <stdlib.h>

#include "H5hut.h"

const char* FNAME = "simple_dunetest.h5";

typedef struct vertex {
	h5_float64_t P[3];
} vertex_t; 

typedef struct elem {
	h5_loc_idx_t vids[3];
} elem_t;
	       
vertex_t Vertices[] = {
	{ { 0.0,  0.0,  0.0} },
	{ {-0.5,  1.0,  0.0} },
	{ { 0.5,  1.0,  0.0} },
	{ { 1.0,  0.0,  0.0} },
	{ { 0.5, -1.0,  0.0} },
	{ {-0.5, -1.0,  0.0} },
	{ {-1.0,  0.0,  0.0} }
};

elem_t Elems[] = {
	{ { 0, 1, 2 } },
	{ { 0, 2, 3 } },
	{ { 0, 3, 4 } },
	{ { 0, 4, 5 } },
	{ { 0, 5, 6 } },
	{ { 0, 1, 6 } }
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
	H5SetVerbosityLevel (2);

	/* open file and add mesh */
	const h5_file_t f = H5OpenFile (FNAME, H5_O_WRONLY, 0);
	h5t_mesh_t* m;
	H5FedAddTriangleMesh (f, "0", &m);

	/* store vertices */
	H5FedBeginStoreVertices (m, num_vertices);
	int i;
	for (i = 0; i < num_vertices; i++) {
		H5FedStoreVertex (m, -1, Vertices[i].P);
	}
	H5FedEndStoreVertices (m);

	/* store elements */
	H5FedBeginStoreElements (m, num_elems);
	for (i = 0; i < num_elems; i++) {
		H5FedStoreElement (m, Elems[i].vids);
	}
	H5FedEndStoreElements (m);

	H5FedCloseMesh (m);
	H5CloseFile (f);
	return 0;
}
