#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <strings.h>

#include "H5hut.h"

const char* FNAME = "simple_tet.h5";


typedef struct timer {
	clock_t _start;
	clock_t _end;
	struct timer* (*new)(void);
	void (*delete)(struct timer*);
	void (*start)(struct timer*);
	void (*stop)(struct timer*);
	double (*elapsed)(struct timer*);
} Timer;

extern Timer Timer_;
static Timer*
new (
	void
	) {
	Timer* this = (Timer*)malloc (sizeof (Timer));
	*this = Timer_;
	return this;
}

static void
delete (Timer* this) {
	free (this);
}

static void
start (Timer* this) {
	this->_start = clock();
}

static void
stop (Timer* this) {
	this->_end = clock();
}

static double
elapsed (Timer* this) {
	return (double)(this->_end - this->_start)/CLOCKS_PER_SEC;
}

Timer Timer_ = {
	0,
	0,
	new,
	delete,
	start,
	stop,
	elapsed
};

static h5_err_t
traverse_vertices (
	h5t_mesh_t* const m,
	int dumpit,
	Timer* timer
	) {
	printf ( "  %-32s", "Traversing vertices: ");

	/* get number of vertices we have to expect */
	h5_size_t num_vertices_expect = H5FedGetNumVerticesTotal (m);

	/* get iterator for co-dim 3 entities, i.e vertices */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 3);

	/* iterate  */
	h5_loc_id_t local_id;
	h5_size_t num_vertices = 0;
	timer->start(timer);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		if (dumpit) {
			printf ("\n");
			h5_float64_t P[3];
			H5FedGetVertexCoordsByID (m, local_id, P);
			char v[256];
			snprintf (v, sizeof(v), "=%llx=", (long long)local_id);
			printf ("| %-18s | (%f, %f, %f) |\n",
				v, P[0], P[1], P[2]);
		}
		num_vertices++;
	}

	/* done */
	H5FedEndTraverseEntities (iter);
	timer->stop(timer);

	/* report error if we got a different number then expected */
	if (num_vertices != num_vertices_expect) {
		fprintf (stderr, "!!! Got %lud vertices, but expected %lud.\n",
			 (unsigned long)num_vertices,
			 (unsigned long)num_vertices_expect);
	}

	printf ("  %fsec, number of vertices: %lu\n",
		timer->elapsed(timer),
		(unsigned long)num_vertices);
	return H5_SUCCESS;
}

static h5_err_t
traverse_edges (
	h5t_mesh_t* const m,
	int dumpit,
	Timer* timer
	) {
	printf ( "  %-32s", "Traversing edges: ");
	timer->start(timer);
	/* get iterator for co-dim 2 entities, i.e. edges */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 2);

	/* iterate */
	h5_loc_id_t local_id;
	h5_size_t num_edges = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		if (dumpit) {
			printf ("\n");
			char v[256];
			char k[256];
			h5_loc_id_t local_vids[4];
			snprintf ( k, sizeof(k), "=%llx=", (long long)local_id );
			H5FedGetVertexIndicesOfEntity ( m, local_id, local_vids );
			snprintf ( v, sizeof(v), "=[%lld,%lld]=",
				   (long long)local_vids[0], (long long)local_vids[1] );
			printf ( "| %-18s | %-18s |\n", k, v );
		}
		num_edges++;
	}

	/* done */
	H5FedEndTraverseEntities (iter);
	timer->stop(timer);
	printf ("  %fsec, number of edges: %lu\n",
		timer->elapsed(timer),
		(unsigned long)num_edges);
	return H5_SUCCESS;
}

static h5_err_t
traverse_triangles (
	h5t_mesh_t* const m,
	int dumpit,
	Timer* timer
	) {
	printf ( "  %-32s", "Traversing triangles: ");

	timer->start (timer);
	/* get iterator for co-dim 1 entities, i.e. triangles */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 1);

	/* iterate */
	h5_loc_id_t local_id;
	h5_size_t num_triangles = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		if (dumpit) {
			printf ("\n");
			char v[256];
			char d[256];
			h5_loc_id_t local_vids[4];
			snprintf ( d, sizeof(d), "=%llx=", (long long)local_id );
			H5FedGetVertexIndicesOfEntity ( m, local_id, local_vids );
			snprintf ( v, sizeof(v), "=[%lld,%lld,%lld]=",
				   (long long)local_vids[0],
				   (long long)local_vids[1],
				   (long long)local_vids[2] );
			printf ( "| %-18s | %-18s |\n", d, v );
		}
		num_triangles++;
	}

	/* done */
	H5FedEndTraverseEntities (iter);
	timer->stop(timer);
	printf ("  %fsec, number of triangles: %lu\n",
		timer->elapsed(timer),
		(unsigned long)num_triangles);
	return H5_SUCCESS;
}

static h5_err_t
traverse_boundary_triangles (
	h5t_mesh_t* const m,
	int dumpit,
	Timer* timer
	) {
	printf ( "  %-32s", "Traversing boundary triangles: ");

	timer->start(timer);
	/* get iterator for co-dim 1 entities, i.e. triangles */
	h5t_iterator_t* iter = H5FedBeginTraverseBoundaryFaces (m, 1);

	/* iterate */
	h5_loc_id_t local_id;
	h5_size_t num_triangles = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		if (dumpit) {
			printf ("\n");
			char v[256];
			char d[256];
			h5_loc_id_t local_vids[4];
			snprintf ( d, sizeof(d), "=%llx=", (long long)local_id );
			H5FedGetVertexIndicesOfEntity ( m, local_id, local_vids );
			snprintf ( v, sizeof(v), "=[%lld,%lld,%lld]=",
				   (long long)local_vids[0],
				   (long long)local_vids[1],
				   (long long)local_vids[2] );
			printf ( "| %-18s | %-18s |\n", d, v );
		}
		num_triangles++;
	}

	/* done */
	H5FedEndTraverseEntities (iter);
	timer->stop(timer);
	printf ("  %fsec, number of boundary triangles: %lu\n",
		timer->elapsed(timer),
		(unsigned long)num_triangles);
	return H5_SUCCESS;
}

static h5_err_t
traverse_elems (
	h5t_mesh_t* const m,
	int dumpit,
	Timer* timer
	) {
	printf ( "  %-32s", "Traversing elements: ");
	timer->start(timer);
	/* get number of elements we have to expect */
	h5_size_t num_elems_expect = H5FedGetNumElementsTotal (m);

	/* get iterator for co-dim 0 */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 0);

	/* iterate over all co-dim 0 entities, i.e. elements */
	h5_loc_id_t local_id;
	h5_size_t num_elems = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		if (dumpit) {
			printf ("\n");
			char v[256];
			char t[256];
			h5_loc_id_t local_vids[4];
			snprintf ( t, sizeof(t), "=%llx=", (long long)local_id );
			H5FedGetVertexIndicesOfEntity (m, local_id, local_vids);
			snprintf ( v, sizeof(v), "=[%lld,%lld,%lld,%lld]=",
				   (long long)local_vids[0], (long long)local_vids[1],
				   (long long)local_vids[2], (long long)local_vids[3] );
			printf ( "| %-18s | %-18s |\n", t, v );
		}
		num_elems++;
	}

	/* done */
	H5FedEndTraverseEntities (iter);
	timer->stop(timer);
	/* report error if we got a different number then expected */
	if (num_elems != num_elems_expect) {
		fprintf (stderr, "!!! Got %lld elements, but expected %lld.\n",
			 (long long)num_elems, (long long)num_elems_expect);
		exit(1);
	}

	printf ("  %fsec, number of elements: %lld\n",
		timer->elapsed(timer),
		(long long)num_elems);
	return H5_SUCCESS;
}


static h5_err_t
traverse_level (
	h5t_mesh_t* const m,
	const h5_lvl_idx_t level_id,
	int dumpit,
	Timer* timer
	) {
	printf ("Setting level to %d\n", level_id);
	H5FedSetLevel (m, level_id);
	traverse_vertices (m, dumpit, timer);
	traverse_edges (m, dumpit, timer);
	traverse_triangles (m, dumpit, timer);
	traverse_boundary_triangles (m, dumpit, timer);
	traverse_elems (m, dumpit, timer);
	return H5_SUCCESS;
}

static h5_err_t
traverse_mesh (
	h5_file_t const f,
	const h5_id_t mesh_id,
	int dumpit
	) {
	h5t_mesh_t* mesh;
	Timer* timer = Timer_.new();

	timer->start(timer);
	/* open mesh and get number of levels */
	printf ("    Opening mesh with id %lld\n", (long long)mesh_id);
	H5FedOpenTetrahedralMeshByIndex (f, mesh_id, &mesh);
	timer->stop(timer);
	printf ("  %fsec\n", timer->elapsed(timer));
	h5_size_t num_levels = H5FedGetNumLevels (mesh);
	printf ("    Number of levels in mesh: %lld\n", (long long)num_levels);

	/* loop over all levels */
	h5_lvl_idx_t level_id;
	for (level_id = 0; level_id < num_levels; level_id++) {
		traverse_level (mesh, level_id, dumpit, timer);
	}
	/* done */
	H5FedCloseMesh (mesh);
	timer->delete(timer);
	return H5_SUCCESS;
}

int
main (
	int argc,
	char* argv[]
	) {

	/* abort program on error, so we don't have to handle them */
	H5SetErrorHandler (H5AbortErrorhandler);
	H5SetVerbosityLevel (0);

	/* open file and get number of meshes */
	h5_file_t f = H5OpenFile (FNAME, H5_O_RDONLY, H5_PROP_DEFAULT);
	h5_size_t num_meshes = H5FedGetNumTetrahedralMeshes (f);
	printf ("    Number of meshes: %lld\n", (long long)num_meshes);

	/* loop over all meshes */
	h5_id_t mesh_id;
	int dumpit = 0;
	for (mesh_id = 0; mesh_id < num_meshes; mesh_id++) {
		traverse_mesh (f, mesh_id, dumpit);
	}

	/* done */
	H5CloseFile (f);
	return 0;
}
