#include <stdio.h>
#include <stdlib.h>

#include "H5hut.h"

#ifndef H5_HAVE_PARALLEL
#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
#endif
#endif

struct vertex {
	h5_glb_id_t global_id;
	h5_float64_t P[3];
};

typedef struct vertex vertex_t; 

struct tet {
	h5_glb_id_t global_id;
	h5_glb_id_t parent_id;
	h5_glb_id_t vids[4];
};
typedef struct tet tet_t;

struct boundary {
	h5_glb_id_t vids[3];
};
typedef struct boundary boundary_t;

vertex_t V0[5] = {
	{ 0, {-1.0,  0.0,  0.0} },
	{ 1, { 1.0,  0.0,  0.0} },
	{ 2, { 0.0,  1.0,  0.0} },
	{ 3, { 0.0,  0.0,  1.0} },
	{ 4, { 0.0, -1.0,  0.0} }
};

vertex_t V1[1] = {
	{ 5, {0.0,  0.0,  0.0 } }
};

// sorted vertices: 0, 4, 5, 3, 2, 1

tet_t T0[2] = {
	{ 1, -1, { 0, 1, 2, 3 } },	// 0, 3, 2, 1
	{ 0, -1, { 0, 1, 3, 4 } }	// 0, 4, 3, 1
};

tet_t T1[2] = {
	{ 2, 0, { 0, 3, 4, 5 } },	// 0, 4, 5, 3
	{ 3, 0, { 1, 3, 4, 5 } }	// 4, 5, 3, 1
};

boundary_t  boundary[2] = {
	{{ 0, 1, 2 }},
	{{ 2, 3, 4 }}
};

int
main (
	int argc,
	char *argv[]
	) {
	H5PartSetVerbosityLevel ( 4 );

	h5_file_t *f = H5OpenFile ( "simple_tet.h5", 0 );
	if ( f == NULL ) {
		fprintf ( stderr, "!!! Can't open file.\n" );
		return -1;
	}

	h5_err_t h5err = H5FedOpenMesh ( f, 0, TETRAHEDRAL_MESH );
	if ( h5err < 0 ) {
		fprintf ( stderr, "!!! Can't open mesh %d\n", 0 );
		return -1;
	}
	h5err = H5FedAddBoundary ( f );
	if ( h5err < 0 ) {
		fprintf ( stderr, "!!! Can't add boundary.\n" );
		return -1;
	}
	h5err = H5FedAddNumBoundaryfaces ( f, 1 );
	if ( h5err < 0 ) {
		fprintf ( stderr, "!!! Can't add boundary.\n" );
		return -1;
	}
	h5err = H5FedStoreBoundaryface ( f, boundary[0].vids );
	if ( h5err < 0 ) {
		fprintf ( stderr, "!!! Can't write boundary.\n" );
		return -1;
	}

	h5err = H5CloseFile ( f );
	if ( h5err < 0 ) {
		fprintf ( stderr, "!!! Can't close file.\n" );
		return -1;
	}
	return 0;
}
