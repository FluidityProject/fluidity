#include <stdlib.h>
#include "testframe.h"
#include "params.h"

static void
test_write_file_attribs(h5_file_t file, int position)
{
	h5_err_t status;
	char name[ATTR_NAME_SIZE];

	TEST("Writing file attributes");

	get_attr_name(name, "str", position);
	status = H5WriteFileAttribString(file, name, ATTR_STR_VAL);
	RETURN(status, H5_SUCCESS, "H5WriteFileAttribString");

	get_attr_name(name, "i32", position);
	h5_int32_t i32 = ATTR_INT32_VAL;
	status = H5WriteFileAttribInt32(file, name, &i32, 1);
	RETURN(status, H5_SUCCESS, "H5WriteFileAttribInt32");

	get_attr_name(name, "i64", position);
	h5_int64_t i64 = ATTR_INT64_VAL;
	status = H5WriteFileAttribInt64(file, name, &i64, 1);
	RETURN(status, H5_SUCCESS, "H5WriteFileAttribInt64");

	get_attr_name(name, "f32", position);
	h5_float32_t f32 = ATTR_FLOAT_VAL;
	status = H5WriteFileAttribFloat32(file, name, &f32, 1);
	RETURN(status, H5_SUCCESS, "H5WriteFileAttribFloat32");

	get_attr_name(name, "f64", position);
	h5_float64_t f64 = ATTR_FLOAT_VAL;
	status = H5WriteFileAttribFloat64(file, name, &f64, 1);
	RETURN(status, H5_SUCCESS, "H5WriteFileAttribFloat64");
}

static void
test_write_step_attribs(h5_file_t file, int position)
{
	h5_err_t status;
	char name[ATTR_NAME_SIZE];

	TEST("Writing step attributes");

	get_attr_name(name, "str", position);
	status = H5WriteStepAttribString(file, name, ATTR_STR_VAL);
	RETURN(status, H5_SUCCESS, "H5WriteStepAttribString");

	get_attr_name(name, "i32", position);
	h5_int32_t i32 = ATTR_INT32_VAL;
	status = H5WriteStepAttribInt32(file, name, &i32, 1);
	RETURN(status, H5_SUCCESS, "H5WriteStepAttribInt32");

	get_attr_name(name, "i64", position);
	h5_int64_t i64 = ATTR_INT64_VAL;
	status = H5WriteStepAttribInt64(file, name, &i64, 1);
	RETURN(status, H5_SUCCESS, "H5WriteStepAttribInt64");

	get_attr_name(name, "f32", position);
	h5_float32_t f32 = ATTR_FLOAT_VAL;
	status = H5WriteStepAttribFloat32(file, name, &f32, 1);
	RETURN(status, H5_SUCCESS, "H5WriteStepAttribFloat32");

	get_attr_name(name, "f64", position);
	h5_float64_t f64 = ATTR_FLOAT_VAL;
	status = H5WriteStepAttribFloat64(file, name, &f64, 1);
	RETURN(status, H5_SUCCESS, "H5WriteStepAttribFloat64");
}

static void
test_write_data64(h5_file_t file, int nparticles, int step)
{
	int i,t;
	h5_int64_t status, val;

	double *x,*y,*z;
	double *px,*py,*pz;
	h5_int64_t *id;

	x=(double*)malloc(nparticles*sizeof(double));
	y=(double*)malloc(nparticles*sizeof(double));
	z=(double*)malloc(nparticles*sizeof(double));
	px=(double*)malloc(nparticles*sizeof(double));
	py=(double*)malloc(nparticles*sizeof(double));
	pz=(double*)malloc(nparticles*sizeof(double));
	id=(h5_int64_t*)malloc(nparticles*sizeof(h5_int64_t));

	TEST("Writing 64-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		for (i=0; i<nparticles; i++)
		{
			x[i]  = 0.0 + (double)(i+nparticles*t);
			y[i]  = 0.1 + (double)(i+nparticles*t);
			z[i]  = 0.2 + (double)(i+nparticles*t);
			px[i] = 0.3 + (double)(i+nparticles*t);
			py[i] = 0.4 + (double)(i+nparticles*t);
			pz[i] = 0.5 + (double)(i+nparticles*t);
			id[i] = i + nparticles*t;
		}

		val = H5HasStep(file, t);

		status = H5SetStep(file, t);
		RETURN(status, H5_SUCCESS, "H5SetStep");

		if (val == 0) test_write_step_attribs(file, t);

		status = H5PartSetNumParticles(file, nparticles);
		RETURN(status, H5_SUCCESS, "H5PartSetNumParticles");

		status = H5PartWriteDataFloat64(file, "x", x);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "y", y);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "z", z);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "px", px);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "py", py);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "pz", pz);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataInt64(file, "id", id);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataInt64");
	}
}

static void
test_write_strided_data64(h5_file_t file, int nparticles, int step)
{
	int i,t;
	h5_int64_t status;

	double *data;

	data=(double*)malloc(6*nparticles*sizeof(double));

	status = H5PartSetNumParticlesStrided(file, nparticles, 6);
	RETURN(status, H5_SUCCESS, "H5PartSetNumParticlesStrided");

	TEST("Writing 64-bit strided data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		for (i=0; i<nparticles; i++)
		{
			data[6*i]   = 0.0 + (double)(i+nparticles*t);
			data[6*i+1] = 0.1 + (double)(i+nparticles*t);
			data[6*i+2] = 0.2 + (double)(i+nparticles*t);
			data[6*i+3] = 0.3 + (double)(i+nparticles*t);
			data[6*i+4] = 0.4 + (double)(i+nparticles*t);
			data[6*i+5] = 0.5 + (double)(i+nparticles*t);
		}

		status = H5SetStep(file, t);
		RETURN(status, H5_SUCCESS, "H5SetStep");

		status = H5PartSetNumParticlesStrided(file, nparticles, 6);
		RETURN(status, H5_SUCCESS, "H5PartSetNumParticlesStrided");

		status = H5PartWriteDataFloat64(file, "x", data+0);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "y", data+1);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "z", data+2);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "px", data+3);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "py", data+4);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "pz", data+5);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat64");

		test_write_step_attribs(file, t);
	}
}

static void
test_write_data32(h5_file_t file, int nparticles, int step)
{
	int i,t;
	h5_err_t status;
	h5_size_t val;

	float *x,*y,*z;
	float *px,*py,*pz;
	int *id;

	x=(float*)malloc(nparticles*sizeof(float));
	y=(float*)malloc(nparticles*sizeof(float));
	z=(float*)malloc(nparticles*sizeof(float));
	px=(float*)malloc(nparticles*sizeof(float));
	py=(float*)malloc(nparticles*sizeof(float));
	pz=(float*)malloc(nparticles*sizeof(float));
	id=(int*)malloc(nparticles*sizeof(int));

	status = H5PartSetNumParticles(file, nparticles);
	RETURN(status, H5_SUCCESS, "H5PartSetNumParticles");

#if defined(H5_HAVE_PARALLEL)
	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#else
	int rank = 0;
#endif

	TEST("Writing 32-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		for (i=0; i<nparticles; i++)
		{
			x[i]  = 0.0F + (float)(i+nparticles*t);
			y[i]  = 0.1F + (float)(i+nparticles*t);
			z[i]  = 0.2F + (float)(i+nparticles*t);
			px[i] = 0.3F + (float)(i+nparticles*t);
			py[i] = 0.4F + (float)(i+nparticles*t);
			pz[i] = 0.5F + (float)(i+nparticles*t);
			id[i] = i + nparticles*t;
		}

		val = H5HasStep(file, t);
		if (val == 0) {
			status = H5SetStep(file, t);
			RETURN(status, H5_SUCCESS, "H5SetStep");
		}

		/* test a two-part write using views */
		status = H5PartSetView(file,
		                       rank*nparticles,
		                       rank*nparticles + 31);
		RETURN(status, H5_SUCCESS, "H5PartSetView");

		status = H5PartWriteDataFloat32(file, "x", x);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		test_write_step_attribs(file, t);

		status = H5PartWriteDataFloat32(file, "y", y);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "z", z);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "px", px);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "py", py);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "pz", pz);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataInt32(file, "id", id);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataInt32");

		/* the second write phase... */
		status = H5PartSetView(file,
		                       rank*nparticles + 32,
		                       rank*nparticles + nparticles - 1);
		RETURN(status, H5_SUCCESS, "H5PartSetView");
		/* offset the input arrays */
		status = H5PartWriteDataFloat32(file, "x", x+32);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "y", y+32);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "z", z+32);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "px", px+32);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "py", py+32);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "pz", pz+32);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataInt32(file, "id", id+32);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataInt32");
	}
}

static void
test_write_strided_data32(h5_file_t file, int nparticles, int step)
{
	int i,t;
	h5_int64_t status;

	float *data;

	data=(float*)malloc(6*nparticles*sizeof(float));

	TEST("Writing 32-bit strided data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		for (i=0; i<nparticles; i++)
		{
			data[6*i]   = 0.0F + (float)(i+nparticles*t);
			data[6*i+1] = 0.1F + (float)(i+nparticles*t);
			data[6*i+2] = 0.2F + (float)(i+nparticles*t);
			data[6*i+3] = 0.3F + (float)(i+nparticles*t);
			data[6*i+4] = 0.4F + (float)(i+nparticles*t);
			data[6*i+5] = 0.5F + (float)(i+nparticles*t);
		}

		status = H5SetStep(file, t);
		RETURN(status, H5_SUCCESS, "H5SetStep");

		test_write_step_attribs(file, t);

		status = H5PartSetNumParticlesStrided(file, nparticles, 6);
		RETURN(status, H5_SUCCESS, "H5PartSetNumParticlesStrided");

		status = H5PartWriteDataFloat32(file, "x", data+0);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "y", data+1);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "z", data+2);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "px", data+3);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "py", data+4);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "pz", data+5);
		RETURN(status, H5_SUCCESS, "H5PartWriteDataFloat32");
	}
}

void h5u_test_write1(void)
{
	h5_file_t file1;
	h5_err_t status;

	TEST("Opening file once, write-truncate");
        h5_prop_t props = H5CreateFileProp ();

#if defined(H5_HAVE_PARALLEL)
        MPI_Comm comm = MPI_COMM_WORLD;
        status = H5SetPropFileMPIOCollective (props, &comm);
	RETURN(status, H5_SUCCESS, "H5SetPropFileMPIOCollective");
#endif
	file1 = H5OpenFile(FILENAME, H5_O_WRONLY, props);
	status = H5CheckFile(file1);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

        status = H5CloseProp (props);
	RETURN(status, H5_SUCCESS, "H5CloseProp");

	test_write_data32(file1, NPARTICLES, 1);
	test_write_file_attribs(file1, 0);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}

void h5u_test_write2(void)
{
	h5_file_t file1;
	h5_file_t file2;
	h5_err_t status;

	TEST("Opening file twice, write-append + read-only");
        h5_prop_t props = H5CreateFileProp ();

#if defined(H5_HAVE_PARALLEL)
        MPI_Comm comm = MPI_COMM_WORLD;
        status = H5SetPropFileMPIOCollective (props, &comm);
	RETURN(status, H5_SUCCESS, "H5SetPropFileMPIOCollective");
#endif
	file1 = H5OpenFile(FILENAME, H5_O_APPENDONLY, props);
	status = H5CheckFile(file1);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

	file2 = H5OpenFile(FILENAME, H5_O_RDONLY, props);
	status = H5CheckFile(file2);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

        status = H5CloseProp (props);
	RETURN(status, H5_SUCCESS, "H5CloseProp");

	test_write_strided_data32(file1, NPARTICLES, NTIMESTEPS+1);
	test_write_file_attribs(file1, 1);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
	status = H5CloseFile(file2);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}

void h5u_test_write3(void)
{
	h5_file_t file1;
	h5_err_t status;

#if H5_VERSION_LE(1,8,12)
	TEST("Opening file once, write-truncate, MPI-POSIX VFD");
#else
	TEST("Opening file once, write-truncate, MPI-IO Independent VFD");
#endif
	
        h5_prop_t props = H5CreateFileProp ();

#if defined(H5_HAVE_PARALLEL)
        MPI_Comm comm = MPI_COMM_WORLD;

#if H5_VERSION_LE(1,8,12)
        status = H5SetPropFileMPIOPosix (props, &comm);
	RETURN(status, H5_SUCCESS, "H5SetPropFileMPIOPosix");
#else
        status = H5SetPropFileMPIOIndependent (props, &comm);
	RETURN(status, H5_SUCCESS, "H5SetPropFileMPIOIndependent");
#endif

#endif // H5_HAVE_PARALLEL

	file1 = H5OpenFile(FILENAME, H5_O_WRONLY, props);

	status = H5CheckFile(file1);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

        status = H5CloseProp (props);
	RETURN(status, H5_SUCCESS, "H5CloseProp");

	TEST("Redefining step name");
	status = H5SetStepNameFormat(file1, "data", 16);
	RETURN(status, H5_SUCCESS, "H5SetStepNameFormat");

	test_write_strided_data64(file1, NPARTICLES, 0);
	test_write_file_attribs(file1, 0);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}

void h5u_test_write4(void)
{
	h5_file_t file1;
	h5_file_t file2;

	h5_err_t status;
        h5_prop_t props;

	TEST("Opening file twice, write-append + read-only, MPI-IO Independent VFD");

        props = H5CreateFileProp ();
#if defined(H5_HAVE_PARALLEL)
        MPI_Comm comm = MPI_COMM_WORLD;
        status = H5SetPropFileMPIOIndependent (props, &comm);
	RETURN(status, H5_SUCCESS, "H5SetPropFileMPIOIndependent");
#endif
	file1 = H5OpenFile(FILENAME, H5_O_APPENDONLY, props);
	status = H5CheckFile(file1);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

	file2 = H5OpenFile(FILENAME, H5_O_RDONLY, props);
	status = H5CheckFile(file2);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

        status = H5CloseProp (props);
	RETURN(status, H5_SUCCESS, "H5CloseProp");

	TEST("Redefining step name");
	status = H5SetStepNameFormat(file1, "data", 16);
	RETURN(status, H5_SUCCESS, "H5SetStepNameFormat");

	status = H5SetStepNameFormat(file2, "data", 16);
	RETURN(status, H5_SUCCESS, "H5SetStepNameFormat");

	status = H5PartSetChunkSize(file1, NPARTICLES);
	RETURN(status, H5_SUCCESS, "H5PartSetChunk");

	test_write_data64(file1, NPARTICLES, NTIMESTEPS-2);
	test_write_file_attribs(file1, 2);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
	status = H5CloseFile(file2);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}

