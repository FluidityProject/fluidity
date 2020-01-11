#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

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
set_vertex_tags (
	h5t_mesh_t* m,
	h5t_tagset_t* const tagset,
	int verify,
	Timer* timer
	) {
	h5_loc_id_t local_id;
	h5_int64_t val[3];
	h5_info ("Tagging all vertices ...");
	timer->start(timer);
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 3);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		val[0] = local_id;
		val[1] = local_id+1;
		val[2] = local_id+2;
		H5FedSetTag (tagset, local_id, 3, val);
		if (verify) {
			h5_int64_t retval[3];
			h5_size_t dim = 3;
			H5FedGetTag (tagset, local_id, &dim, retval);
			if (memcmp (val, retval, sizeof(val))) {
				h5_warn ("Oops on entity %llx!", (long long)local_id);
			}
		}
		h5_debug ("Tagging %llx", (long long)local_id);
	}
	timer->stop(timer);
	h5_info ("  Time to tag to all vertices: %fsec", timer->elapsed(timer)); 
	return H5FedEndTraverseEntities (iter);
}

static h5_err_t
set_edge_tags (
	h5t_mesh_t* const m,
	h5t_tagset_t* tagset,
	int dumpit,
	Timer* timer
	) {
	h5_loc_id_t local_id;
	h5_int64_t val[3];
	h5_info ("Tagging all edges ...");
	timer->start(timer);
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 2);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		val[0] = local_id;
		val[1] = local_id+1;
		val[2] = local_id+2;
		H5FedSetTag (tagset, local_id, 3, val);
		h5_int64_t retval[3];
		h5_size_t dims;
		H5FedGetTag (tagset, local_id, &dims, retval);
		if (memcmp ( val, retval, sizeof(val))) {
			h5_warn ("Oops on entity %llx!", (long long)local_id);
		}
	}
	timer->stop(timer);
	h5_info ("  Time to tag all edges: %fsec", timer->elapsed(timer)); 
	return H5FedEndTraverseEntities (iter);
}

static h5_err_t
set_tri_tags (
	h5t_mesh_t* const m,
	h5t_tagset_t* tagset,
	int dumpit,
	Timer* timer
	) {
	h5_loc_id_t local_id;
	h5_int64_t val[3];
	h5_info ("Tagging all triangles ...");
	timer->start(timer);
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 1);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		val[0] = local_id;
		val[1] = local_id+1;
		val[2] = local_id+2;
		H5FedSetTag (tagset, local_id, 3, val);
		h5_int64_t retval[3];
		h5_size_t dims;
		H5FedGetTag (tagset, local_id, &dims, retval);
		if (memcmp (val, retval, sizeof(val))) {
			h5_warn ("Oops on entity %llx!", (long long)local_id);
		}
	}
	timer->stop(timer);
	h5_info ("  Time to tag all triangles: %fsec", timer->elapsed(timer)); 
	return H5FedEndTraverseEntities (iter);
}

static h5_err_t
set_tet_tags (
	h5t_mesh_t* const m,
	h5t_tagset_t* tagset,
	int dumpit,
	Timer* timer
	) {
	h5_loc_id_t local_id;
	h5_int64_t val[3];
	h5_info ("Tagging all tetrahedra ...");
	timer->start(timer);
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 0);
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		val[0] = local_id;
		val[1] = local_id+1;
		val[2] = local_id+2;
		H5FedSetTag (tagset, local_id, 3, val);
		h5_int64_t retval[3];
		h5_size_t dims;
		H5FedGetTag (tagset, local_id, &dims, retval);
		if (memcmp (val, retval, sizeof(val))) {
			h5_warn ("Oops on entity %llx!", (long long)local_id);
		}
	}
	timer->stop(timer);
	h5_info ("  Time to tag to all tetrahedra: %fsec", timer->elapsed(timer)); 
	return H5FedEndTraverseEntities (iter);
}

int
main (
	int argc,
	char* argv[]
	) {
	int verify = 1;
	/* abort program on error, so we don't have to handle them */
	H5SetErrorHandler (H5AbortErrorhandler);
	H5SetVerbosityLevel (3);

	/* open file and get number of meshes */
	h5_file_t f = H5OpenFile (FNAME, H5_O_RDWR, 0);
	h5t_mesh_t* mesh;
	Timer* timer = Timer_.new();
	timer->start(timer);
	H5FedOpenTetrahedralMesh (f, "0", &mesh);
	timer->stop(timer);
	h5_info ("  Time to open mesh %fsec", timer->elapsed(timer));

	/* open last level */
	//h5_size_t num_levels = H5FedGetNumLevels (mesh);
	//H5FedSetLevel (mesh, num_levels-1);
	H5FedSetLevel (mesh, 0);

	/* add new tagset and write some data to it */
	h5t_tagset_t* tagset = NULL;
	H5FedAddMTagset (mesh, "testtag", H5_INT64_T, &tagset);
	set_vertex_tags (mesh, tagset, verify, timer);
	set_edge_tags (mesh, tagset, verify, timer);
	set_tri_tags (mesh, tagset, verify, timer);
	set_tet_tags (mesh, tagset, verify, timer);

	// close tagset
	timer->start(timer);
	H5FedCloseMTagset (tagset);
	timer->stop(timer);
	h5_info ("  Time to write tagset to disk %fsec", timer->elapsed(timer));

	// close mesh
	timer->start(timer);
	H5FedCloseMesh (mesh);
	timer->stop(timer);
	h5_info ("  Time to close mesh %fsec", timer->elapsed(timer));

	// close file
	timer->start(timer);
	H5CloseFile (f);
	timer->stop(timer);
	h5_info ("  Time to close file %fsec", timer->elapsed(timer));
	return 0;
}
