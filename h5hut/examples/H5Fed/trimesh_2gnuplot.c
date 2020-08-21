#include <stdio.h>
#include <stdlib.h>

#include "H5hut.h"
#if defined (H5_HAVE_PARALLEL)
#include <mpi.h>
#endif

#define H5FedGetNumLeafElementsTotal H5FedGetNumElementsTotal
/*
  Traverse elements and output coordinates
 */
static h5_err_t
traverse_elems (
	h5t_mesh_t* const m
	) {
	/* get number of elements we have to expect */
	h5_size_t num_elems_expect = H5FedGetNumLeafElementsTotal (m);

	/* get iterator for co-dim 0 */
	h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 0);

	/* iterate over all co-dim 0 entities, i.e. elements */
	h5_loc_id_t local_id;
	h5_size_t num_elems = 0;
	while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
		h5_loc_id_t local_vids[4];
		H5FedGetVertexIndicesOfEntity (m, local_id, local_vids);

		printf ("# cell %llx\n", (long long unsigned)local_id );
		h5_float64_t P[3];
		for (int i = 0; i < 3; i++) {
			H5FedGetVertexCoordsByIndex (m, local_vids[i], P);
			printf (" %8.6f %8.6f\n", P[0], P[1]);
		}
		H5FedGetVertexCoordsByIndex (m, local_vids[0], P);
		printf (" %8.6f %8.6f\n", P[0], P[1]);
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

int
main (
	int argc,
	char* argv[]
	) {

	if (argc < 2 || argc > 3) {
		fprintf (stderr, "Usage: %s FILE [LEVEL]\n", argv[0]);
		exit (42);
	}

	/* abort program on error, so we don't have to handle them */
	H5SetErrorHandler (H5AbortErrorhandler);
	H5SetVerbosityLevel (0);

	/* open file and get number of meshes */
	h5_file_t f = H5OpenFile (argv[1], H5_O_RDONLY, H5_PROP_DEFAULT);
	h5t_mesh_t*  m;
	H5FedOpenTriangleMeshByIndex (f, 0, &m);
	int num_levels = H5FedGetNumLevels (m);
	int level = num_levels-1;
	if (argc >= 3) {
		level = atoi (argv[2]);
	}
	if (level >= num_levels || level < 0) {
		fprintf (stderr, "level %d out of range\n", level);
		goto done;
	}
	H5FedSetLevel (m, level);
	traverse_elems (m);


done:
	H5FedCloseMesh (m);
	H5CloseFile (f);

#if defined (H5_HAVE_PARALLEL)
	MPI_Finalize ();
#endif

	return 0;
}
