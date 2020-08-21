#include <stdio.h>
#include <stdlib.h>

#include "H5hut.h"
#if defined (H5_HAVE_PARALLEL)
#include <mpi.h>
#endif

const char* FNAME = "simple_triangle.h5";

static h5_err_t
traverse_vertices (
	h5t_mesh_t* const m
	) {
	/* get number of vertices we have to expect */
	h5_size_t num_vertices_expect = H5FedGetNumVerticesTotal (m);

	/* get iterator for co-dim 2 entities, i.e. vertices */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 2);

	/* iterate */
	h5_loc_id_t local_id;
	h5_size_t num_vertices = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		h5_float64_t P[3];
		H5FedGetVertexCoordsByID (m, local_id, P);
		char v[256];
		snprintf (v, sizeof(v), "=%llx=", (long long)local_id);
		printf ("| %-18s | (%f, %f, %f) |\n",
			 v, P[0], P[1], P[2]);
		num_vertices++;
	}

	/* done */
	H5FedEndTraverseEntities (iter);

	/* report error if we got a different number then expected */
	if (num_vertices != num_vertices_expect) {
		fprintf (stderr, "!!! Got %lld vertices, but expected %lld.\n",
			 (long long)num_vertices, (long long)num_vertices_expect);
		exit (1);
	}

	printf ("    Number of vertices on level: %lld\n", (long long)num_vertices);
	return H5_SUCCESS;
}

static h5_err_t
traverse_edges (
	h5t_mesh_t* const m
	) {
	printf ( "Travering edges on level %lld:\n", (long long)H5FedGetLevel(m) );

	/* get iterator for co-dim 1 entities, i.e. edges */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 1);

	/* iterate */
	h5_loc_id_t local_id;
	h5_size_t num_edges = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		char v[256];
		char k[256];
		h5_loc_id_t local_vids[4];
		snprintf ( k, sizeof(k), "=%llx=", (long long)local_id );
		H5FedGetVertexIndicesOfEntity ( m, local_id, local_vids );
		snprintf ( v, sizeof(v), "=[%lld,%lld]=",
			   (long long)local_vids[0], (long long)local_vids[1] );
		printf ( "| %-18s | %-18s |\n", k, v );
		num_edges++;
	}

	/* done */
	H5FedEndTraverseEntities (iter);

	printf ("    Number of edges: %lld\n", (long long)num_edges);
	return H5_SUCCESS;
}

static h5_err_t
traverse_boundary_edges (
	h5t_mesh_t* const m
	) {
	printf ("Travering boundary edges on level %lld:\n", (long long)H5FedGetLevel(m));

	/* get iterator for co-dim 1 entities, i.e. edges */
	h5t_iterator_t* iter = H5FedBeginTraverseBoundaryFaces (m, 1);

	/* iterate */
	h5_loc_id_t local_id;
	h5_size_t num_edges = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		char v[256];
		char k[256];
		h5_loc_id_t local_vids[4];
		snprintf ( k, sizeof(k), "=%llx=", (long long)local_id);
		H5FedGetVertexIndicesOfEntity (m, local_id, local_vids);
		snprintf ( v, sizeof(v), "=[%lld,%lld]=",
			   (long long)local_vids[0], (long long)local_vids[1] );
		printf ( "| %-18s | %-18s |\n", k, v );
		num_edges++;
	}

	/* done */
	H5FedEndTraverseEntities (iter);

	printf ("    Number of edges: %lld\n", (long long)num_edges);
	return H5_SUCCESS;
}

static h5_err_t
traverse_elems (
	h5t_mesh_t* const m
	) {
	/* get number of elements we have to expect */
	h5_size_t num_elems_expect = H5FedGetNumElementsTotal (m);

	/* get iterator for co-dim 0 */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 0);

	/* iterate over all co-dim 0 entities, i.e. elements */
	h5_loc_id_t local_id;
	h5_size_t num_elems = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		char v[256];
		char t[256];
		h5_loc_id_t local_vids[4];
		snprintf (t, sizeof(t), "=%llx=", (long long)local_id);
		H5FedGetVertexIndicesOfEntity (m, local_id, local_vids);
		snprintf (v, sizeof(v), "=[%lld,%lld,%lld]=",
			  (long long)local_vids[0],
			  (long long)local_vids[1],
			  (long long)local_vids[2]);
		printf ("| %-18s | %-18s |\n", t, v);
		num_elems++;
	}

	/* done */
	H5FedEndTraverseEntities (iter);

	/* report error if we got a different number then expected */
	if (num_elems != num_elems_expect) {
		fprintf (stderr, "!!! Got %lld elements, but expected %lld.\n",
			 (long long)num_elems, (long long)num_elems_expect);
		exit(1);


	}

	printf ("    Number of elements on level: %lld\n", (long long)num_elems);
	return H5_SUCCESS;
}

/*
  Traverse elements and output coordinates
 */
static h5_err_t
traverse_elems2 (
	h5t_mesh_t* const m
	) {
	/* get number of elements we have to expect */
	h5_size_t num_elems_expect = H5FedGetNumElementsTotal (m);

	/* get iterator for co-dim 0 */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 0);

	/* iterate over all co-dim 0 entities, i.e. elements */
	h5_loc_id_t local_id;
	h5_size_t num_elems = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		printf ("%05llu", (unsigned long long)num_elems);
		h5_loc_id_t local_vids[4];
		H5FedGetVertexIndicesOfEntity (m, local_id, local_vids);
		int i;
		for (i = 0; i < 3; i++) {
			h5_float64_t P[3];
			H5FedGetVertexCoordsByIndex (m, local_vids[i], P);
			printf (" %8.6f %8.6f %8.6f", P[0], P[1], P[2]);
		}
		printf ("\n");
		num_elems++;
	}

	/* done */
	H5FedEndTraverseEntities (iter);

	/* report error if we got a different number then expected */
	if (num_elems != num_elems_expect) {
		fprintf (stderr, "!!! Got %lld elements, but expected %lld.\n",
			 (long long)num_elems, (long long)num_elems_expect);
		exit(1);
	}
	return H5_SUCCESS;
}

static h5_err_t
traverse_level (
	h5t_mesh_t* const m,
	const h5_loc_id_t level_id
	) {
	printf ("    Setting level to %d\n", level_id);
	H5FedSetLevel (m, level_id);
	traverse_vertices (m);
	traverse_edges (m);
	traverse_boundary_edges (m);
	traverse_elems (m);
	traverse_elems2 (m);
	return H5_SUCCESS;
}

static h5_err_t
traverse_mesh (
	const h5_file_t f,
	const h5_id_t mesh_id
	) {
	h5t_mesh_t* m;
	/* open mesh and get number of levels */
	printf ("    Opening mesh with id %lld\n", (long long)mesh_id);
	H5FedOpenTriangleMeshByIndex (f, mesh_id, &m);
	h5_size_t num_levels = H5FedGetNumLevels (m);
	printf ("    Number of levels in mesh: %lld\n", (long long)num_levels);

	/* loop over all levels */
	h5_lvl_idx_t level_id;
	for (level_id = 0; level_id < num_levels; level_id++) {
		traverse_level (m, level_id);
	}
	/* done */
	H5FedCloseMesh (m);
	return H5_SUCCESS;
}


int
main (
	int argc,
	char* argv[]
	) {

	/* abort program on error, so we don't have to handle them */
	H5SetErrorHandler (H5AbortErrorhandler);
	H5SetVerbosityLevel (H5_DEBUG_ALL);

	/* open file and get number of meshes */
	h5_file_t f = H5OpenFile (FNAME, H5_O_RDONLY, H5_PROP_DEFAULT);
	h5_size_t num_meshes = H5FedGetNumTriangleMeshes (f);
	printf ("    Number of meshes: %lld\n", (long long)num_meshes);

	/* loop over all meshes */
	h5_id_t mesh_id;
	for (mesh_id = 0; mesh_id < num_meshes; mesh_id++) {
		traverse_mesh (f, mesh_id);
	}

	/* done */
	H5CloseFile (f);
#if defined (H5_HAVE_PARALLEL)
	MPI_Finalize ();
#endif

	return 0;
}
