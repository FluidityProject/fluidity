#include <stdlib.h>
#include <H5hut.h>

#include "testframe.h"
#include "params.h"

/* from write.c */
void h5u_test_write1(void);
void h5u_test_write2(void);
#if H5_VERSION_LE(1,8,12)
void h5u_test_write3(void);
#endif
void h5u_test_write4(void);

/* from read.c */
void h5u_test_read1(void);
void h5u_test_read2(void);
void h5u_test_read3(void);
void h5u_test_read4(void);

int main(int argc, char **argv)
{
#ifdef H5_HAVE_PARALLEL
	MPI_Init(&argc, &argv);
	int procs;
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	if (procs > MAX_MPI_TASKS) {
		fprintf(stderr,
		        "ERROR: please use <= %d MPI tasks for the test.\n",
		        MAX_MPI_TASKS);
		exit(EXIT_FAILURE);
	}
#endif

	/* Initialize testing framework */
	TestInit(argv[0], NULL, NULL);

	/* Tests are generally arranged from least to most complexity... */
	AddTest("write1", h5u_test_write1, NULL, "Write 32-bit data", NULL);
	AddTest("read1", h5u_test_read1, NULL, "Read 32-bit data", NULL);
	AddTest("write2", h5u_test_write2, NULL, "Write 32-bit strided data", NULL);
	AddTest("read2", h5u_test_read2, NULL, "Read 32-bit strided data", NULL);
#if H5_VERSION_LE(1,8,12)
	AddTest("write3", h5u_test_write3, NULL, "Write 64-bit strided data", NULL);
	AddTest("read3", h5u_test_read3, NULL, "Read 64-bit strided data", NULL);
#endif
	AddTest("write4", h5u_test_write4, NULL, "Write 64-bit data", NULL);
	AddTest("read4", h5u_test_read4, NULL, "Read 64-bit data", NULL);

	/* Display testing information */
	TestInfo(argv[0]);

	/* Parse command line arguments */
	TestParseCmdLine(argc, argv);

	H5SetVerbosityLevel(GetTestVerbosity());

	/* Perform requested testing */
	PerformTests();

	/* Display test summary, if requested */
	if (GetTestSummary())
		TestSummary();

	/* Clean up test files, if allowed */
	//if (GetTestCleanup() && !getenv("HDF5_NOCLEANUP"))
	//    TestCleanup();

#ifdef H5_HAVE_PARALLEL
	TestPrintf ("reached end\n");
	fflush(stdout);
	MPI_Finalize();
#endif
	return GetTestNumErrs();
}

