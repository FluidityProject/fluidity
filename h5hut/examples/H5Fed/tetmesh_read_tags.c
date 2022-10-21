#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "H5hut.h"

// name of input file
const char* FNAME = "simple_tet.h5";

// H5hut verbosity/debug level
const h5_int64_t h5_verbosity = H5_VERBOSE_DEFAULT;
const h5_int64_t h5_debugmsk = H5_DEBUG_ALL;

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
	Timer* self = (Timer*)malloc (sizeof (Timer));
	*self = Timer_;
	return self;
}

static void
delete (Timer* self) {
	free (self);
}

static void
start (Timer* self) {
	self->_start = clock();
}

static void
stop (Timer* self) {
	self->_end = clock();
}

static double
elapsed (Timer* self) {
	return (double)(self->_end - self->_start)/CLOCKS_PER_SEC;
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

struct vertex {
	h5_float64_t P[3];
};

typedef struct vertex vertex_t; 

struct tet {
	h5_glb_id_t global_id;
	h5_glb_id_t parent_id;
	h5_glb_id_t vids[4];
};
typedef struct tet tet_t;

static h5_err_t
traverse_vertices (
	h5t_mesh_t* const m,
	h5t_tagset_t* tagset,
	int dumpit,
	Timer* timer
	) {
	printf ( "  %-32s", "Traversing vertices: ");

	/* get iterator for co-dim 3 entities, i.e vertices */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 3);

	/* iterate  */
	h5_loc_id_t local_id;
	h5_size_t num_entities = 0;
	h5_size_t num_entities_tagged = 0;
	timer->start(timer);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		num_entities++;
		if (dumpit) {
			printf ("\n");
			h5_float64_t P[3];
			H5FedGetVertexCoordsByID (m, local_id, P);
			char v[256];
			snprintf (v, sizeof(v), "=%llx=", (long long)local_id);
			printf ("| %-18s | (%f, %f, %f) |",
				v, P[0], P[1], P[2]);
		}

		h5_size_t size = 3;
		h5_int64_t tval[3];
		if ((local_id = H5FedGetTag (tagset, local_id, &size, tval)) < 0) {
			continue;	// not tagged
		}
		if (dumpit) {
			printf (" (%llx, %llx, %llx) |",
				(long long)tval[0], (long long)tval[1], (long long)tval[2]);
		}
		num_entities_tagged++;
		if (tval[0] != local_id ||
		    tval[1] != local_id+1 ||
		    tval[2] != local_id+2) {
			fprintf (stderr,
				 "!!! Wrong tag values for entity %lld\n",
				 (long long)local_id);
			exit (1);
		}
	}
	/* done */
	H5FedEndTraverseEntities (iter);
	timer->stop(timer);

	if (dumpit) {
		printf ("\n\nTime elapsed:");
	}
	printf ("  %fsec, total number of entities: %lu; tagged: %lu\n",
		timer->elapsed(timer),
		(unsigned long)num_entities, (unsigned long)num_entities_tagged);
	return H5_SUCCESS;
}

static h5_err_t
traverse_edges (
	h5t_mesh_t* const m,
	h5t_tagset_t* tagset,
	int dumpit,
	Timer* timer
	) {
	printf ( "  %-32s", "Traversing edges: ");
	timer->start(timer);

	/* get iterator for co-dim 2 entities, i.e. edges */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 2);

	/* iterate */
	h5_loc_id_t local_id;
	h5_size_t num_entities = 0;
	h5_size_t num_entities_tagged = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		num_entities++;
		if (dumpit) {
			printf ("\n");
			char v[256];
			char k[256];
			h5_loc_id_t local_vids[4];
			snprintf ( k, sizeof(k), "=%llx=", (long long)local_id );
			H5FedGetVertexIndicesOfEntity ( m, local_id, local_vids );
			snprintf ( v, sizeof(v), "=[%lld,%lld]=",
				   (long long)local_vids[0], (long long)local_vids[1] );
			printf ( "| %-18s | %-18s |", k, v );
		}
		h5_size_t size = 3;
		h5_int64_t tval[3];
		if ((local_id = H5FedGetTag (tagset, local_id, &size, tval)) < 0) 
			continue;	// not tagged
		if (dumpit) {
			printf (" (%llx, %llx, %llx) |",
				(long long)tval[0], (long long)tval[1], (long long)tval[2]);
		}
		num_entities_tagged++;

		if (tval[0] != local_id ||
		    tval[1] != local_id+1 ||
		    tval[2] != local_id+2) {
			fprintf (stderr,
				 "!!! Wrong tag values for entity %lld\n",
				 (long long)local_id);
			exit (1);
		}
	}

	/* done */
	H5FedEndTraverseEntities (iter);
	timer->stop(timer);

	if (dumpit) {
		printf ("\n\nTime elapsed:");
	}
	printf ("  %fsec, total number of entities: %lu; tagged: %lu\n",
		timer->elapsed(timer),
		(unsigned long)num_entities, (unsigned long)num_entities_tagged);
	return H5_SUCCESS;
}

static h5_err_t
traverse_triangles (
	h5t_mesh_t* const m,
	h5t_tagset_t* tagset,
	int dumpit,
	Timer* timer
	) {
	printf ( "  %-32s", "Traversing triangles: ");
	timer->start(timer);
	/* get iterator for co-dim 1 entities, i.e. triangles */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 1);

	/* iterate */
	h5_loc_id_t local_id;
	h5_size_t num_entities = 0;
	h5_size_t num_entities_tagged = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		num_entities++;
		if (dumpit) {
			printf ("\n");
			char v[256];
			char k[256];
			h5_loc_id_t local_vids[4];
			snprintf ( k, sizeof(k), "=%llx=", (long long)local_id );
			H5FedGetVertexIndicesOfEntity ( m, local_id, local_vids );
			snprintf ( v, sizeof(v), "=[%lld,%lld,%lld]=",
				   (long long)local_vids[0],
				   (long long)local_vids[1],
				   (long long)local_vids[2] );
			printf ( "| %-18s | %-18s |", k, v );
		}
		h5_size_t size = 3;
		h5_int64_t tval[3];
		if ((local_id = H5FedGetTag (tagset, local_id, &size, tval)) < 0) 
			continue;	// not tagged
		if (dumpit) {
			printf (" (%llx, %llx, %llx) |",
				(long long)tval[0], (long long)tval[1], (long long)tval[2]);
		}
		num_entities_tagged++;

		if (tval[0] != local_id || tval[1] != local_id+1 || tval[2] != local_id+2) {
			fprintf (stderr,
				 "!!! Wrong tag values for entity %lld\n",
				 (long long)local_id);
			exit (1);
		}
	}

	/* done */
	H5FedEndTraverseEntities (iter);
	timer->stop(timer);

	if (dumpit) {
		printf ("\n\nTime elapsed:");
	}
	printf ("  %fsec, total number of entities: %lu; tagged: %lu\n",
		timer->elapsed(timer),
		(unsigned long)num_entities, (unsigned long)num_entities_tagged);
	return H5_SUCCESS;
}

static h5_err_t
traverse_tets (
	h5t_mesh_t* const m,
	h5t_tagset_t* tagset,
	int dumpit,
	Timer* timer
	) {
	printf ( "  %-32s", "Traversing tetrahedra: ");
	timer->start(timer);
	/* get iterator for co-dim 0 entities, i.e. tetrahedra */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 0);
	/* iterate */
	h5_loc_id_t local_id;
	h5_size_t num_entities = 0;
	h5_size_t num_entities_tagged = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		num_entities++;
		if (dumpit) {
			printf ("\n");
			char v[256];
			char k[256];
			h5_loc_id_t local_vids[4];
			snprintf ( k, sizeof(k), "=%llx=", (long long)local_id );
			H5FedGetVertexIndicesOfEntity ( m, local_id, local_vids );
			snprintf ( v, sizeof(v), "=[%lld,%lld,%lld,%lld]=",
				   (long long)local_vids[0],
				   (long long)local_vids[1],
				   (long long)local_vids[2],
				   (long long)local_vids[3] );
			printf ( "| %-18s | %-18s |", k, v );
		}
		h5_size_t size = 3;
		h5_int64_t tval[3];
		if ((local_id = H5FedGetTag (tagset, local_id, &size, tval)) < 0) 
			continue;	// not tagged
		if (dumpit) {
			printf (" (%llx, %llx, %llx) |",
				(long long)tval[0], (long long)tval[1], (long long)tval[2]);
		}
		num_entities_tagged++;

		if (tval[0] != local_id || tval[1] != local_id+1 || tval[2] != local_id+2) {
			fprintf (stderr,
				 "!!! Wrong tag values for entity %lld\n",
				 (long long)local_id);
			exit (1);
		}
	}

	/* done */
	H5FedEndTraverseEntities (iter);
	timer->stop(timer);

	if (dumpit) {
		printf ("\n\nTime elapsed:");
	}
	printf ("  %fsec, total number of entities: %lu; tagged: %lu\n",
		timer->elapsed(timer),
		(unsigned long)num_entities, (unsigned long)num_entities_tagged);
	return H5_SUCCESS;
}

static h5_err_t
traverse_level (
	h5t_mesh_t* const m,
	h5t_tagset_t* tagset,
	const h5_lvl_idx_t level_id,
	int dumpit,
	Timer* timer
	) {
	printf ("Setting level to %d\n", level_id);
	H5FedSetLevel (m, level_id);
	traverse_vertices (m, tagset, dumpit, timer);
	traverse_edges (m, tagset, dumpit, timer);
	traverse_triangles (m, tagset, dumpit, timer);
	traverse_tets (m, tagset, dumpit, timer);
	return H5_SUCCESS;
}

static h5_err_t
traverse_mesh (
	const h5_file_t f,
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

	h5t_tagset_t* tagset;
	H5FedOpenMTagset (mesh, "testtag", &tagset);

	h5_size_t num_levels = H5FedGetNumLevels (mesh);
	printf ("    Number of levels in mesh: %lld\n", (long long)num_levels);

	/* loop over all levels */
	h5_lvl_idx_t level_id;
	for (level_id = 0; level_id < num_levels; level_id++) {
		traverse_level (mesh, tagset, level_id, dumpit, timer);
	}
	/* done */
	H5FedCloseMTagset (tagset);
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
	H5SetVerbosityLevel (h5_verbosity);
	H5SetDebugMask (h5_debugmsk);
	
	/* open file and get number of meshes */
	h5_file_t f = H5OpenFile (FNAME, H5_O_RDONLY, 0);
	h5_size_t num_meshes = H5FedGetNumTetrahedralMeshes (f);
	printf ("    Number of meshes: %lld\n", (long long)num_meshes);

	/* loop over all meshes */
	h5_id_t mesh_id;
	int dumpit = 1;
	for (mesh_id = 0; mesh_id < num_meshes; mesh_id++) {
		traverse_mesh (f, mesh_id, dumpit);
	}

	/* done */
	H5CloseFile (f);
	return 0;
}
