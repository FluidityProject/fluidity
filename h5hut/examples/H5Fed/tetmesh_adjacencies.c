#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "H5hut.h"

const char* FNAME = "simple_tet.h5";

typedef struct timer {
	clock_t _start;
	clock_t _end;
	clock_t _elapsed;
	struct timer* (*new)(void);
	void (*delete)(struct timer*);
	void (*start)(struct timer*);
	void (*stop)(struct timer*);
	void (*reset)(struct timer*);
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
	this->_elapsed += this->_end - this->_start;
	this->_start = this->_end;
}

static void
reset (Timer* this) {
	this->_start = this->_end = this->_elapsed = 0;
}

static double
elapsed (Timer* this) {
	return (double)(this->_elapsed)/CLOCKS_PER_SEC;
}

Timer Timer_ = {
	0,
	0,
	0,
	new,
	delete,
	start,
	stop,
	reset,
	elapsed
};

static h5_err_t
print_adjacencies_of_vertex (
	h5t_mesh_t* const m,
	h5_loc_id_t local_id,
	int dumpit,
	Timer* timer
	) {
	h5_loc_idlist_t* uadj_edges;
	h5_loc_idlist_t* uadj_triangles;
	h5_loc_idlist_t* uadj_tets;
	timer->reset(timer);
	timer->start(timer);
	H5FedGetAdjacencies (m, local_id, 1, &uadj_edges);
	H5FedGetAdjacencies (m, local_id, 2, &uadj_triangles);
	H5FedGetAdjacencies (m, local_id, 3, &uadj_tets);
	timer->stop(timer);
	int n = uadj_tets->num_items;
	if (uadj_triangles->num_items > n) n = uadj_triangles->num_items;
	if (uadj_edges->num_items > n) n = uadj_edges->num_items;
	if (dumpit) {
		for (int i = 0; i < n; i++) {
			char v[256];
			char k[256];
			char d[256];
			char t[256];
			h5_loc_id_t local_vids[4];
			if (i == 0) {
				H5FedGetVertexIndicesOfEntity (
					m, local_id, local_vids);
				snprintf (v, sizeof(v), "=[%lld]=",
					  (long long)local_vids[0]);
			} else {
				*v = '\0';
			}
			if (i < uadj_edges->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, uadj_edges->items[i], local_vids);
				snprintf (k, sizeof(k), "=[%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1]);
			} else {
				*k = '\0';
			}
			if (i < uadj_triangles->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, uadj_triangles->items[i], local_vids);
				snprintf (d, sizeof(d), "=[%lld,%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1],
					  (long long)local_vids[2]);
			} else {
				*d = '\0';
			}
			if (i < uadj_tets->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, uadj_tets->items[i], local_vids);
				snprintf (t, sizeof(t), "=[%lld,%lld,%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1],
					  (long long)local_vids[2],
					  (long long)local_vids[3]);
			} else {
				*t = '\0';
			}
			printf ("| %-18s | %-18s | %-18s | %-18s |\n", v, k, d, t);
		}
	}
	H5FedReleaseListOfAdjacencies (m, &uadj_edges);
	H5FedReleaseListOfAdjacencies (m, &uadj_triangles);
	H5FedReleaseListOfAdjacencies (m, &uadj_tets);
	return H5_SUCCESS;
}

static h5_err_t
print_adjacencies_of_edge (
	h5t_mesh_t* const m,
	h5_loc_id_t local_id,
	int dumpit,
	Timer* timer
	) {
	h5_loc_idlist_t* dadj_vertices;
	h5_loc_idlist_t* uadj_triangles;
	h5_loc_idlist_t* uadj_tets;
	timer->reset(timer);
	timer->start(timer);
	H5FedGetAdjacencies (m, local_id, 0, &dadj_vertices);
	H5FedGetAdjacencies (m, local_id, 2, &uadj_triangles);
	H5FedGetAdjacencies (m, local_id, 3, &uadj_tets);
	timer->stop(timer);
	int n = dadj_vertices->num_items;
	if (uadj_triangles->num_items > n) n = uadj_triangles->num_items;
	if (uadj_tets->num_items > n) n = uadj_tets->num_items;

	if (dumpit) {
		for (int i = 0; i < n; i++) {
			char v[256];
			char k[256];
			char d[256];
			char t[256];
			h5_loc_id_t local_vids[4];
			if (i < dadj_vertices->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, dadj_vertices->items[i], local_vids);
				snprintf (v, sizeof(v), "=[%lld]=",
					  (long long)local_vids[0]);
			} else {
				*v = '\0';
			}
			if (i == 0) {
				H5FedGetVertexIndicesOfEntity (
					m, local_id, local_vids);
				snprintf (k, sizeof(k), "=[%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1]);
			} else {
				*k = '\0';
			}
			if (i < uadj_triangles->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, uadj_triangles->items[i], local_vids);
				snprintf (d, sizeof(d), "=[%lld,%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1],
					  (long long)local_vids[2]);
			} else {
				*d = '\0';
			}
			if (i < uadj_tets->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, uadj_tets->items[i], local_vids);
				snprintf (t, sizeof(t), "=[%lld,%lld,%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1],
					  (long long)local_vids[2],
					  (long long)local_vids[3]);
			} else {
				*t = '\0';
			}
			printf ("| %-18s | %-18s | %-18s | %-18s |\n", v, k, d, t);
		}
	}
	H5FedReleaseListOfAdjacencies (m, &dadj_vertices);
	H5FedReleaseListOfAdjacencies (m, &uadj_triangles);
	H5FedReleaseListOfAdjacencies (m, &uadj_tets);
	return H5_SUCCESS;
}

static h5_err_t
print_adjacencies_of_triangle (
	h5t_mesh_t* const m,
	h5_loc_id_t local_id,
	int dumpit,
	Timer* timer
	) {
	h5_loc_idlist_t* dadj_vertices;
	h5_loc_idlist_t* dadj_edges;
	h5_loc_idlist_t* uadj_tets;
	timer->reset(timer);
	timer->start(timer);
	H5FedGetAdjacencies (m, local_id, 0, &dadj_vertices);
	H5FedGetAdjacencies (m, local_id, 1, &dadj_edges);
	H5FedGetAdjacencies (m, local_id, 3, &uadj_tets);
	timer->stop(timer);
	int n = dadj_vertices->num_items;
	if (dadj_edges->num_items > n) n = dadj_edges->num_items;
	if (uadj_tets->num_items > n) n = uadj_tets->num_items;
	if (dumpit) {
		for (int i = 0; i < n; i++) {
			char v[256];
			char k[256];
			char d[256];
			char t[256];
			h5_loc_id_t local_vids[4];
			if (i < dadj_vertices->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, dadj_vertices->items[i], local_vids);
				snprintf (v, sizeof(v), "=[%lld]=",
					  (long long)local_vids[0]);
			} else {
				*v = '\0';
			}
			if (i < dadj_edges->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, dadj_edges->items[i], local_vids);
				snprintf (k, sizeof(k), "=[%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1]);
			} else {
				*k = '\0';
			}
			if (i == 0) {
				H5FedGetVertexIndicesOfEntity (
					m, local_id, local_vids);
				snprintf (d, sizeof(k), "=[%lld,%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1],
					  (long long)local_vids[2]);
			} else {
				*d = '\0';
			}
			if (i < uadj_tets->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, uadj_tets->items[i], local_vids);
				snprintf (t, sizeof(t), "=[%lld,%lld,%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1],
					  (long long)local_vids[2],
					  (long long)local_vids[3]);
			} else {
				*t = '\0';
			}
			printf ( "| %-18s | %-18s | %-18s | %-18s |\n", v, k, d, t );
		}
	}
	H5FedReleaseListOfAdjacencies (m, &dadj_vertices);
	H5FedReleaseListOfAdjacencies (m, &dadj_edges);
	H5FedReleaseListOfAdjacencies (m, &uadj_tets);
	return H5_SUCCESS;
}

static h5_err_t
print_adjacencies_of_tet (
	h5t_mesh_t* const m,
	h5_loc_id_t local_id,
	int dumpit,
	Timer* timer
	) {
	h5_loc_idlist_t* dadj_vertices;
	h5_loc_idlist_t* dadj_edges;
	h5_loc_idlist_t* dadj_triangles;
	timer->reset(timer);
	timer->start(timer);
	H5FedGetAdjacencies (m, local_id, 0, &dadj_vertices);
	H5FedGetAdjacencies (m, local_id, 1, &dadj_edges);
	H5FedGetAdjacencies (m, local_id, 2, &dadj_triangles);
	timer->stop(timer);
	int n = dadj_vertices->num_items;
	if (dadj_edges->num_items > n) n = dadj_edges->num_items;
	if (dadj_triangles->num_items > n) n = dadj_triangles->num_items;
	if (dumpit) {
		for (int i = 0; i < n; i++) {
			char v[256];
			char k[256];
			char d[256];
			char t[256];
			h5_loc_id_t local_vids[4];
			if (i < dadj_vertices->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, dadj_vertices->items[i], local_vids);
				snprintf (v, sizeof(v), "=[%lld]=",
					  (long long)local_vids[0]);
			} else {
				*v = '\0';
			}
			if (i < dadj_edges->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, dadj_edges->items[i], local_vids);
				snprintf (k, sizeof(k), "=[%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1]);
			} else {
				*k = '\0';
			}
			if (i < dadj_triangles->num_items) {
				H5FedGetVertexIndicesOfEntity (
					m, dadj_triangles->items[i], local_vids);
				snprintf (d, sizeof(d), "=[%lld,%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1],
					  (long long)local_vids[2]);
			} else {
				*d = '\0';
			}
			if (i == 0) {
				H5FedGetVertexIndicesOfEntity (
					m, local_id, local_vids);
				snprintf (d, sizeof(k), "=[%lld,%lld,%lld,%lld]=",
					  (long long)local_vids[0],
					  (long long)local_vids[1],
					  (long long)local_vids[2],
					  (long long)local_vids[2]);
			} else {
				*t = '\0';
			}
			printf ("| %-18s | %-18s | %-18s | %-18s |\n", v, k, d, t);
		}
	}
	H5FedReleaseListOfAdjacencies (m, &dadj_vertices);
	H5FedReleaseListOfAdjacencies (m, &dadj_edges);
	H5FedReleaseListOfAdjacencies (m, &dadj_triangles);
	return H5_SUCCESS;
}

static h5_err_t
traverse_vertices (
	h5t_mesh_t* const m,
	int dumpit,
	Timer* timer
	) {
	h5_loc_id_t local_id;
	h5_loc_idx_t num = 0;
	double t_total = 0.0;
	double t_min = CLOCKS_PER_SEC;
	double t_max = 0.0;
	double t = 0.0;

	printf ("\nAdjacencies to vertices\n");
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 3);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		print_adjacencies_of_vertex (m, local_id, dumpit, timer);
		num++;
		t = timer->elapsed(timer);
		t_total += t;
		if (t < t_min) t_min = t;
		if (t > t_max) t_max = t;
	}
	fprintf (
		stderr,
		"%lld\ttotal: %f\tmin: %f\tavg: %f\tmax: %f\n",
		(long long)num,
		t_total,
		t_min,
		t_total / (double)num,
		t_max);
	return H5FedEndTraverseEntities (iter);
}

static h5_err_t
traverse_edges (
	h5t_mesh_t* const m,
	int dumpit,
	Timer* timer
	) {
	h5_loc_id_t local_id;
	h5_loc_idx_t num = 0;
	double t_total = 0.0;
	double t_min = CLOCKS_PER_SEC;
	double t_max = 0.0;
	double t = 0.0;
	printf ("\nAdjacencies to edges\n");
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 2);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		print_adjacencies_of_edge (m, local_id, dumpit, timer);
		num++;
		t = timer->elapsed(timer);
		t_total += t;
		if (t < t_min) t_min = t;
		if (t > t_max) t_max = t;
	}
	fprintf (
		stderr,
		"%lld\ttotal: %f\tmin: %f\tavg: %f\tmax: %f\n",
		(long long)num,
		t_total,
		t_min,
		t_total / (double)num,
		t_max);
	return H5FedEndTraverseEntities (iter);
}

static h5_err_t
traverse_triangles (
	h5t_mesh_t* const m,
	int dumpit,
	Timer* timer
	) {
	h5_loc_id_t local_id;
	h5_loc_idx_t num = 0;
	double t_min = CLOCKS_PER_SEC;
	double t_max = 0.0;
	double t_total = 0.0;
	double t = 0.0;
	printf ("\nAdjacencies to triangle\n");
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 1);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		print_adjacencies_of_triangle (m, local_id, dumpit, timer);
		num++;
		t = timer->elapsed(timer);
		t_total += t;
		if (t < t_min) t_min = t;
		if (t > t_max) t_max = t;
	}
	fprintf (
		stderr,
		"%lld\ttotal: %f\tmin: %f\tavg: %f\tmax: %f\n",
		(long long)num,
		t_total,
		t_min,
		t_total / (double)num,
		t_max);
	return H5FedEndTraverseEntities (iter);
}

static h5_err_t
traverse_elems (
	h5t_mesh_t* const m,
	int dumpit,
	Timer* timer
	) {
	h5_loc_id_t local_id;
	h5_loc_idx_t num = 0;
	double t_total = 0.0;
	double t_min = CLOCKS_PER_SEC;
	double t_max = 0.0;
	double t = 0.0;
	printf ("\nAdjacencies to tetrahedra\n");
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 0);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		print_adjacencies_of_tet (m, local_id, dumpit, timer);
		num++;
		t = timer->elapsed(timer);
		t_total += t;
		if (t < t_min) t_min = t;
		if (t > t_max) t_max = t;
	}
	fprintf (
		stderr,
		"%lld\ttotal: %f\tmin: %f\tavg: %f\tmax: %f\n",
		(long long)num,
		t_total,
		t_min,
		t_total / (double)num,
		t_max);
	return H5FedEndTraverseEntities (iter);
}

static h5_err_t
traverse_level (
	h5t_mesh_t* const m,
	const h5_lvl_idx_t level_id,
	int dumpit,
	Timer* timer
	) {
	printf ("    Setting level to %d\n", level_id);
	H5FedSetLevel (m, level_id);
	traverse_vertices (m, dumpit, timer);
	traverse_edges (m, dumpit, timer);
	traverse_triangles (m, dumpit, timer);
	traverse_elems (m, dumpit, timer);
	return H5_SUCCESS;
}

static h5_err_t
traverse_mesh (
	h5_file_t const f,
	const h5_id_t mesh_id
	) {
	h5t_mesh_t* mesh;
	int dumpit = 0;
	Timer* timer = Timer_.new();

	/* open mesh and get number of levels */
	printf ("    Opening mesh with id %lld\n", (long long)mesh_id);
	H5FedOpenTetrahedralMeshByIndex (f, mesh_id, &mesh);
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
	//H5SetVerbosityLevel (H5_DEBUG_ALL);
	H5SetVerbosityLevel (0);
	
	/* open file and get number of meshes */
	h5_file_t f = H5OpenFile (FNAME, H5_O_RDONLY, 0);
	h5_size_t num_meshes = H5FedGetNumTetrahedralMeshes (f);
	printf ("    Number of meshes: %lld\n", (long long)num_meshes);

	/* loop over all meshes */
	h5_id_t mesh_id;
	for (mesh_id = 0; mesh_id < num_meshes; mesh_id++) {
		traverse_mesh (f, mesh_id);
	}

	/* done */
	H5CloseFile (f);
	return 0;
}
