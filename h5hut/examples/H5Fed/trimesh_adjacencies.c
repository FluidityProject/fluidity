#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "H5hut.h"

#define PRINT_UPADJACENCIES 1
#define PRINT_DOWNADJACENCIES 0

const char* FNAME = "simple_triangle.h5";

static h5_err_t
print_adjacencies_of_vertex (
	h5t_mesh_t* const m,
	h5_loc_id_t local_id,
	clock_t* time_used
	) {
	h5_loc_idlist_t* uadj_edges;
	h5_loc_idlist_t* uadj_triangles;
	clock_t time_used_before = clock();
	H5FedGetAdjacencies (m, local_id, 1, &uadj_edges);
	H5FedGetAdjacencies (m, local_id, 2, &uadj_triangles);
	*time_used = clock() - time_used_before;
	int n = 0;
	if (uadj_triangles->num_items > n) n = uadj_triangles->num_items;
	if (uadj_edges->num_items > n) n = uadj_edges->num_items;
	int i;
	for (i = 0; i < n; i++) {
		char v[256];
		char k[256];
		char d[256];
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
		printf ("| %-18s | %-18s | %-18s |\n", v, k, d);
	}
	H5FedReleaseListOfAdjacencies (m, &uadj_edges);
	H5FedReleaseListOfAdjacencies (m, &uadj_triangles);
	return H5_SUCCESS;
}

static h5_err_t
print_adjacencies_of_edge (
	h5t_mesh_t* const m,
	h5_loc_id_t local_id,
	clock_t* time_used
	) {
	h5_loc_idlist_t* dadj_vertices;
	h5_loc_idlist_t* uadj_triangles;
	clock_t time_used_before = clock();
	H5FedGetAdjacencies (m, local_id, 0, &dadj_vertices);
	H5FedGetAdjacencies (m, local_id, 2, &uadj_triangles);
	*time_used = clock() - time_used_before;
	int n = dadj_vertices->num_items;
	if (uadj_triangles->num_items > n) n = uadj_triangles->num_items;

	int i;
	for (i = 0; i < n; i++) {
		char v[256];
		char k[256];
		char d[256];
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
		printf ("| %-18s | %-18s | %-18s |\n", v, k, d);
	}
	H5FedReleaseListOfAdjacencies (m, &dadj_vertices);
	H5FedReleaseListOfAdjacencies (m, &uadj_triangles);
	return H5_SUCCESS;
}

static h5_err_t
print_adjacencies_of_elem (
	h5t_mesh_t* const m,
	h5_loc_id_t local_id,
	clock_t* time_used
	) {
	h5_loc_idlist_t* dadj_vertices;
	h5_loc_idlist_t* dadj_edges;
	clock_t time_used_before = clock();
	H5FedGetAdjacencies (m, local_id, 0, &dadj_vertices);
	H5FedGetAdjacencies (m, local_id, 1, &dadj_edges);
	*time_used = clock() - time_used_before;
	int n = dadj_vertices->num_items;
	if (dadj_edges->num_items > n) n = dadj_edges->num_items;
	int i;
	for (i = 0; i < n; i++) {
		char v[256];
		char k[256];
		char d[256];
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
		printf ( "| %-18s | %-18s | %-18s |\n", v, k, d );
	}
	H5FedReleaseListOfAdjacencies (m, &dadj_vertices);
	H5FedReleaseListOfAdjacencies (m, &dadj_edges);
	return H5_SUCCESS;
}

static h5_err_t
traverse_vertices (
	h5t_mesh_t* const m
	) {
	h5_loc_id_t local_id;
	h5_loc_idx_t num = 0;
	clock_t t_total = 0;
	clock_t t_min = CLOCKS_PER_SEC;
	clock_t t_max = 0;
	clock_t t = 0;
	printf ("\nAdjacencies to vertices\n");
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 2);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		print_adjacencies_of_vertex (m, local_id, &t);
		num++;
		t_total += t;
		if (t < t_min) t_min = t;
		if (t > t_max) t_max = t;
	}
	fprintf (
		stderr,
		"%lld\ttotal: %f\tmin: %f\tavg: %f\tmax: %f\n",
		(long long)num,
		(double)t_total / (double)CLOCKS_PER_SEC,
		(double)t_min / (double)CLOCKS_PER_SEC,
		(double)t_total / (double)CLOCKS_PER_SEC / (double)num,
		(double)t_max / (double)CLOCKS_PER_SEC );
	return H5FedEndTraverseEntities (iter);
}

static h5_err_t
traverse_edges (
	h5t_mesh_t* const m
	) {
	h5_loc_id_t local_id;
	h5_loc_idx_t num = 0;
	clock_t t_total = 0;
	clock_t t_min = CLOCKS_PER_SEC;
	clock_t t_max = 0;
	clock_t t = 0;
	printf ("\nAdjacencies to edges\n");
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 1);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		print_adjacencies_of_edge (m, local_id, &t);
		num++;
		t_total += t;
		if (t < t_min) t_min = t;
		if (t > t_max) t_max = t;
	}
	fprintf (
		stderr,
		"%lld\ttotal: %f\tmin: %f\tavg: %f\tmax: %f\n",
		(long long)num,
		(double)t_total / (double)CLOCKS_PER_SEC,
		(double)t_min / (double)CLOCKS_PER_SEC,
		(double)t_total / (double)CLOCKS_PER_SEC / (double)num,
		(double)t_max / (double)CLOCKS_PER_SEC );
	return H5FedEndTraverseEntities (iter);
}

static h5_err_t
traverse_elems (
	h5t_mesh_t* const m
	) {
	h5_loc_id_t local_id;
	h5_loc_idx_t num = 0;
	clock_t t_total = 0;
	clock_t t_min = CLOCKS_PER_SEC;
	clock_t t_max = 0;
	clock_t t = 0;
	printf ("\nAdjacencies to triangles\n");
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 0);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		print_adjacencies_of_elem (m, local_id, &t);
		num++;
		t_total += t;
		if (t < t_min) t_min = t;
		if (t > t_max) t_max = t;
	}
	fprintf (
		stderr,
		"%lld\ttotal: %f\tmin: %f\tavg: %f\tmax: %f\n",
		(long long)num,
		(double)t_total / (double)CLOCKS_PER_SEC,
		(double)t_min / (double)CLOCKS_PER_SEC,
		(double)t_total / (double)CLOCKS_PER_SEC / (double)num,
		(double)t_max / (double)CLOCKS_PER_SEC);
	return H5FedEndTraverseEntities (iter);
}

static h5_err_t
traverse_level (
	h5t_mesh_t* const m,
	const h5_lvl_idx_t level_id
	) {
	printf ("    Setting level to %d\n", level_id);
	H5FedSetLevel (m, level_id);
	traverse_vertices (m);
	traverse_edges (m);
	traverse_elems (m);
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
	H5SetVerbosityLevel (2);

	/* open file and get number of meshes */
	h5_file_t f = H5OpenFile (FNAME, H5_O_RDONLY, 0);
	h5_size_t num_meshes = H5FedGetNumTriangleMeshes (f);
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
