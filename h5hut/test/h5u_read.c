#include <stdlib.h>
#include <string.h>
#include "testframe.h"
#include "params.h"

static void
test_read_file_attribs(h5_file_t file, int position)
{
	h5_err_t status;
	char name[ATTR_NAME_SIZE];
	char str[ATTR_NAME_SIZE];
	h5_int32_t i32;
	h5_int64_t i64;
	h5_float32_t f32;
	h5_float64_t f64;

	TEST("Reading file attributes");

	i64 = H5GetNumFileAttribs(file);
	VALUE(i64 % 5, 0, "file attribute count");

	get_attr_name(name, "str", position);
	status = H5ReadFileAttribString(file, name, str);
	RETURN(status, H5_SUCCESS, "H5ReadFileAttribString");
	SVALUE(str, ATTR_STR_VAL, "string attribute");

	get_attr_name(name, "i32", position);
	status = H5ReadFileAttribInt32(file, name, &i32);
	RETURN(status, H5_SUCCESS, "H5ReadFileAttribInt32");
	IVALUE(i32, ATTR_INT32_VAL, "int32 attribute");

	get_attr_name(name, "i64", position);
	status = H5ReadFileAttribInt64(file, name, &i64);
	RETURN(status, H5_SUCCESS, "H5ReadFileAttribInt64");
	IVALUE(i64, ATTR_INT64_VAL, "int64 attribute");

	get_attr_name(name, "f32", position);
	status = H5ReadFileAttribFloat32(file, name, &f32);
	RETURN(status, H5_SUCCESS, "H5ReadFileAttribFloat32");
	FVALUE(f32, ATTR_FLOAT_VAL, "float32 attribute");

	get_attr_name(name, "f64", position);
	status = H5ReadFileAttribFloat64(file, name, &f64);
	RETURN(status, H5_SUCCESS, "H5ReadFileAttribFloat64");
	FVALUE(f64, ATTR_FLOAT_VAL, "float64 attribute");
}

static void
test_read_step_attribs(h5_file_t file, int position)
{
	h5_err_t status;
	char name[ATTR_NAME_SIZE];
	char str[ATTR_NAME_SIZE];
	h5_int32_t i32;
	h5_int64_t i64;
	h5_float32_t f32;
	h5_float64_t f64;

	TEST("Reading file attributes");

	i64 = H5GetNumStepAttribs(file);
	VALUE(i64 % 5, 0, "step attribute count");

	get_attr_name(name, "str", position);
	status = H5ReadStepAttribString(file, name, str);
	RETURN(status, H5_SUCCESS, "H5ReadStepAttribString");
	SVALUE(str, ATTR_STR_VAL, "string attribute");

	get_attr_name(name, "i32", position);
	status = H5ReadStepAttribInt32(file, name, &i32);
	RETURN(status, H5_SUCCESS, "H5ReadStepAttribInt32");
	IVALUE(i32, ATTR_INT32_VAL, "int32 attribute");

	get_attr_name(name, "i64", position);
	status = H5ReadStepAttribInt64(file, name, &i64);
	RETURN(status, H5_SUCCESS, "H5ReadStepAttribInt64");
	IVALUE(i64, ATTR_INT64_VAL, "int64 attribute");

	get_attr_name(name, "f32", position);
	status = H5ReadStepAttribFloat32(file, name, &f32);
	RETURN(status, H5_SUCCESS, "H5ReadStepAttribFloat32");
	FVALUE(f32, ATTR_FLOAT_VAL, "float32 attribute");

	get_attr_name(name, "f64", position);
	status = H5ReadStepAttribFloat64(file, name, &f64);
	RETURN(status, H5_SUCCESS, "H5ReadStepAttribFloat64");
	FVALUE(f64, ATTR_FLOAT_VAL, "float64 attribute");
}

static void
test_read_data64(h5_file_t file, int nparticles, int step)
{
	int i,t;
	int rank, nprocs;
	h5_err_t status;
	h5_int64_t val, start, end, type;
	char name1[4];
	char name2[8];
	h5_size_t indices[8];
	h5_size_t size;

	double *x,*y,*z;
	double *px,*py,*pz;
	h5_int64_t *id;


	TEST("Verifying dataset info");

#if defined(H5_HAVE_PARALLEL)
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#else
	nprocs = 1;
	rank = 2;
#endif

	x=(double*)malloc(nprocs*nparticles*sizeof(double));
	y=(double*)malloc(nprocs*nparticles*sizeof(double));
	z=(double*)malloc(nprocs*nparticles*sizeof(double));
	px=(double*)malloc(nprocs*nparticles*sizeof(double));
	py=(double*)malloc(nprocs*nparticles*sizeof(double));
	pz=(double*)malloc(nprocs*nparticles*sizeof(double));
	id=(h5_int64_t*)malloc(nprocs*nparticles*sizeof(h5_int64_t));

	val = H5PartGetNumParticles(file);
	IVALUE(val, nprocs*nparticles, "particle count");

	val = H5PartGetNumDatasets(file);
	IVALUE(val, 7, "dataset count");

	for (i=0; i<7; i++) {
		status = H5PartGetDatasetName(file, i, name1, 2);
		RETURN(status, H5_SUCCESS, "H5PartGetDatasetName");

		status = H5PartGetDatasetInfo(
		        file, i, name2, 4, &type, &size);
		RETURN(status, H5_SUCCESS, "H5PartGetDatasetInfo");
		CVALUE(name1[0], name2[0], "dataset name");

		status = H5PartGetDatasetName(file, i, name1, 4);
		RETURN(status, H5_SUCCESS, "H5PartGetDatasetName");
		CVALUE(name1[1], name2[1], "dataset name");

		IVALUE(size, nprocs*nparticles, "dataset size");
		if (name1[0] == 'i') IVALUE(type, H5_INT64_T, "dataset type");
		else IVALUE(type, H5_FLOAT64_T, "dataset type");
	}

	TEST("Reading 64-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		val = H5HasStep(file, t);
		IVALUE(val, 1, "has step");

		status = H5SetStep(file, t);
		RETURN(status, H5_SUCCESS, "H5SetStep");

		test_read_step_attribs(file, t);

		status = H5PartSetNumParticles(file, nparticles);
		RETURN(status, H5_SUCCESS, "H5PartSetNumParticles");

		status = H5PartResetView(file);
		RETURN(status, H5_SUCCESS, "H5PartResetView");

		start = rank;
		end = -1;

		status = H5PartSetView(file, start, end);
		RETURN(status, H5_SUCCESS, "H5PartSetView");

		val = H5PartGetView(file, &start, &end);
		IVALUE(val, nprocs*nparticles-start, "particle count");
		IVALUE(start, rank, "view start");
		IVALUE(end, nprocs*nparticles-1, "view end");

		status = H5PartSetView(file, -1, -1);
		RETURN(status, H5_SUCCESS, "H5PartSetView");

		status = H5PartReadDataFloat64(file, "x", x);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");
		IVALUE(x[rank], (double)(rank+nparticles*t), "x data");

		indices[0] = rank*2 + 0;
		indices[1] = rank*2 + 3;
		indices[2] = rank*2 + 9;
		indices[3] = rank*2 + 7;

		status = H5PartSetViewIndices(file, indices, 0);
		RETURN(status, H5_SUCCESS, "H5PartSetViewIndices");

		status = H5PartReadDataFloat64(file, "x", x);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");
		FVALUE(x[2*rank], (double)(2*rank+nparticles*t), "x data");

		status = H5PartResetView(file);
		RETURN(status, H5_SUCCESS, "H5PartResetView");

		status = H5PartSetViewIndices(file, indices, 4);
		RETURN(status, H5_SUCCESS, "H5PartSetViewIndices");

		val = H5PartGetNumParticles(file);
		IVALUE(val, 4, "particle count");

		double x2[4];
		status = H5PartReadDataFloat64(file, "x", x2);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");
		FVALUE(x2[0], (double)(2*rank+0+nparticles*t), "x data");
		FVALUE(x2[1], (double)(2*rank+3+nparticles*t), "x data");
		FVALUE(x2[2], (double)(2*rank+9+nparticles*t), "x data");
		FVALUE(x2[3], (double)(2*rank+7+nparticles*t), "x data");

		val = H5PartGetNumParticles(file);
		IVALUE(val, 4, "particle count");

		status = H5PartSetViewIndices(file, NULL, 4);
		RETURN(status, H5_SUCCESS, "H5PartSetViewIndices");

		status = H5PartReadDataFloat64(file, "x", x);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartSetCanonicalView(file);
		RETURN(status, H5_SUCCESS, "H5PartSetCanonicalView");

		status = H5PartReadDataFloat64(file, "x", x);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "y", y);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "z", z);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "px", px);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "py", py);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "pz", pz);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataInt64(file, "id", id);
		RETURN(status, H5_SUCCESS, "H5PartReadDataInt64");

		for (i=0; i<nparticles; i++)
		{
			FVALUE(x[i], 0.0 + (double)(i+nparticles*t), " x data");
			FVALUE(y[i], 0.1 + (double)(i+nparticles*t), " y data");
			FVALUE(z[i], 0.2 + (double)(i+nparticles*t), " z data");
			FVALUE(px[i], 0.3 + (double)(i+nparticles*t), " px data");
			FVALUE(py[i], 0.4 + (double)(i+nparticles*t), " py data");
			FVALUE(pz[i], 0.5 + (double)(i+nparticles*t), " pz data");
			IVALUE(id[i],               (i+nparticles*t), " id data");
		}
	}
}

static void
test_read_strided_data64(h5_file_t file, int nparticles, int step)
{
	int i,t;
	h5_int64_t status;

	double *data;

	data=(double*)malloc(6*nparticles*sizeof(double));

	TEST("Reading 64-bit strided data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		status = H5SetStep(file, t);
		RETURN(status, H5_SUCCESS, "H5SetStep");

		status = H5PartSetNumParticlesStrided(file, nparticles, 6);
		RETURN(status, H5_SUCCESS, "H5PartSetNumParticlesStrided");

		status = H5PartReadDataFloat64(file, "x", data+0);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "y", data+1);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "z", data+2);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "px", data+3);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		test_read_step_attribs(file, t);

		status = H5PartReadDataFloat64(file, "py", data+4);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "pz", data+5);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat64");

		for (i=0; i<nparticles; i++)
		{
			FVALUE(data[6*i], 0.0 + (double)(i+nparticles*t), "x data");
			FVALUE(data[6*i+1], 0.1 + (double)(i+nparticles*t), "y data");
			FVALUE(data[6*i+2], 0.2 + (double)(i+nparticles*t), "z data");
			FVALUE(data[6*i+3], 0.3 + (double)(i+nparticles*t), "px data");
			FVALUE(data[6*i+4], 0.4 + (double)(i+nparticles*t), "py data");
			FVALUE(data[6*i+5], 0.5 + (double)(i+nparticles*t), "pz data");
		}
	}
}

static void
test_read_data32(h5_file_t file, int nparticles, int step)
{
	int i,t;
	h5_int64_t status, val;

	float *x,*y,*z;
	float *px,*py,*pz;
	int *id;

	x=(float*)malloc(nparticles*sizeof(float));
	y=(float*)malloc(nparticles*sizeof(float));
	z=(float*)malloc(nparticles*sizeof(float));
	px=(float*)malloc(nparticles*sizeof(float));
	py=(float*)malloc(nparticles*sizeof(float));
	pz=(float*)malloc(nparticles*sizeof(float));
	id=(h5_int32_t*)malloc(nparticles*sizeof(h5_int32_t));

	TEST("Reading 32-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		val = H5HasStep(file, t);
		IVALUE(val, 1, "has step");

		status = H5SetStep(file, t);
		RETURN(status, H5_SUCCESS, "H5SetStep");

		status = H5PartSetCanonicalView(file);
		RETURN(status, H5_SUCCESS, "H5PartSetCanonicalView");

		test_read_step_attribs(file, t);
		return;
		status = H5PartReadDataFloat32(file, "x", x);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "y", y);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "z", z);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "px", px);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "py", py);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "pz", pz);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataInt32(file, "id", id);
		RETURN(status, H5_SUCCESS, "H5PartReadDataInt32");

		for (i=0; i<nparticles; i++)
		{
			FVALUE(x[i], 0.0F + (float)(i+nparticles*t), " x data");
			FVALUE(y[i], 0.1F + (float)(i+nparticles*t), " y data");
			FVALUE(z[i], 0.2F + (float)(i+nparticles*t), " z data");
			FVALUE(px[i], 0.3F + (float)(i+nparticles*t), " px data");
			FVALUE(py[i], 0.4F + (float)(i+nparticles*t), " py data");
			FVALUE(pz[i], 0.5F + (float)(i+nparticles*t), " pz data");
			IVALUE(id[i],               (i+nparticles*t), " id data");
		}
	}
}

static void
test_read_strided_data32(h5_file_t file, int nparticles, int step)
{
	int i,t;
	h5_int64_t status;

	float *data;

	data=(float*)malloc(6*nparticles*sizeof(float));

	TEST("Reading 32-bit strided data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		status = H5SetStep(file, t);
		RETURN(status, H5_SUCCESS, "H5SetStep");

		status = H5PartSetNumParticlesStrided(file, nparticles, 6);
		RETURN(status, H5_SUCCESS, "H5PartSetNumParticlesStrided");

		status = H5PartReadDataFloat32(file, "x", data+0);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "y", data+1);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "z", data+2);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "px", data+3);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "py", data+4);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "pz", data+5);
		RETURN(status, H5_SUCCESS, "H5PartReadDataFloat32");

		for (i=0; i<nparticles; i++)
		{
			FVALUE(data[6*i], 0.0F + (float)(i+nparticles*t), "x data");
			FVALUE(data[6*i+1], 0.1F + (float)(i+nparticles*t), "y data");
			FVALUE(data[6*i+2], 0.2F + (float)(i+nparticles*t), "z data");
			FVALUE(data[6*i+3], 0.3F + (float)(i+nparticles*t), "px data");
			FVALUE(data[6*i+4], 0.4F + (float)(i+nparticles*t), "py data");
			FVALUE(data[6*i+5], 0.5F + (float)(i+nparticles*t), "pz data");
		}

		test_read_step_attribs(file, t);
	}
}

void h5u_test_read1(void)
{
	h5_file_t file1;

	h5_int64_t status;

	TEST("Opening file once, read-only");
	file1 = H5OpenFile(FILENAME, H5_O_RDONLY, H5_PROP_DEFAULT);
	status = H5CheckFile(file1);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

	test_read_file_attribs(file1, 0);
	test_read_data32(file1, NPARTICLES, 1);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}

void h5u_test_read2(void)
{
	h5_file_t file1;
	h5_file_t file2;

	h5_int64_t status;

	TEST("Opening file twice, read-only");
	file1 = H5OpenFile(FILENAME, H5_O_RDONLY, H5_PROP_DEFAULT);
	status = H5CheckFile(file1);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

	file2 = H5OpenFile(FILENAME, H5_O_RDONLY, H5_PROP_DEFAULT);
	status = H5CheckFile(file2);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

	test_read_strided_data32(file1, NPARTICLES, NTIMESTEPS+1);
	test_read_file_attribs(file2, 1);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
	status = H5CloseFile(file2);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}

void h5u_test_read3(void)
{
	h5_file_t file1;
	h5_int64_t status;

#if H5_VERSION_LE(1,8,12)
	TEST("Opening file once, read-only, MPI-POSIX VFD");
#else
	TEST("Opening file once, read-only, MPI-IO Independent VFD");
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
	file1 = H5OpenFile(FILENAME, H5_O_RDONLY, props);
	status = H5CheckFile(file1);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

	TEST("Redefining step name");
	status = H5SetStepNameFormat(file1, "data", 16);
	RETURN(status, H5_SUCCESS, "H5SetStepNameFormat");

	test_read_strided_data64(file1, NPARTICLES, 0);
	test_read_file_attribs(file1, 0);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}

void h5u_test_read4(void)
{
	h5_file_t file1;
	h5_file_t file2;
	h5_err_t status;

	TEST("Opening file twice, read-only, MPI-IO Independent VFD");
        h5_prop_t props = H5CreateFileProp ();
#if defined(H5_HAVE_PARALLEL)
        MPI_Comm comm = MPI_COMM_WORLD;
        status = H5SetPropFileMPIOIndependent (props, &comm);
	RETURN(status, H5_SUCCESS, "H5SetPropFileMPIOIndependent");
#endif
	file1 = H5OpenFile(FILENAME, H5_O_RDONLY, props);
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

	test_read_file_attribs(file1, 1);

	status = H5SetStep(file2, NTIMESTEPS);
	RETURN(status, H5_SUCCESS, "H5SetStep");

	test_read_data64(file2, NPARTICLES, NTIMESTEPS-2);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
	status = H5CloseFile(file2);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}

