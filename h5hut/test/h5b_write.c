#include <stdlib.h>
#include "testframe.h"
#include "params.h"

static void
test_write_field_attribs(
	h5_file_t file,
	const char *field_name,
	int position)
{
	h5_err_t status;
	char name[ATTR_NAME_SIZE];

	TEST("Writing field attributes");

	get_attr_name(name, "str", position);
	status = H5BlockWriteFieldAttribString(
	        file, field_name, name, ATTR_STR_VAL);
	RETURN(status, H5_SUCCESS, "H5BlockWriteFieldAttribString");

	get_attr_name(name, "i32", position);
	h5_int32_t i32 = ATTR_INT32_VAL;
	status = H5BlockWriteFieldAttribInt32(
	        file, field_name, name, &i32, 1);
	RETURN(status, H5_SUCCESS, "H5BlockWriteFieldAttribInt32");

	get_attr_name(name, "i64", position);
	h5_int64_t i64 = ATTR_INT64_VAL;
	status = H5BlockWriteFieldAttribInt64(
	        file, field_name, name, &i64, 1);
	RETURN(status, H5_SUCCESS, "H5BlockWriteFieldAttribInt64");

	get_attr_name(name, "f32", position);
	h5_float32_t f32 = ATTR_FLOAT_VAL;
	status = H5BlockWriteFieldAttribFloat32(
	        file, field_name, name, &f32, 1);
	RETURN(status, H5_SUCCESS, "H5BlockWriteFieldAttribFloat32");

	get_attr_name(name, "f64", position);
	h5_float64_t f64 = ATTR_FLOAT_VAL;
	status = H5BlockWriteFieldAttribFloat64(
	        file, field_name, name, &f64, 1);
	RETURN(status, H5_SUCCESS, "H5BlockWriteFieldAttribFloat64");
}

static void
test_write_data64(h5_file_t file, int step)
{
	extern h5_size_t layout[6];

	int i,t;
	h5_int64_t status, val;

	double *e;
	double *ex,*ey,*ez;
	h5_int64_t *id;

	const size_t nelems =
	        (layout[1] - layout[0] + 1) *
	        (layout[3] - layout[2] + 1) *
	        (layout[5] - layout[4] + 1);

	e=(double*)malloc(nelems*sizeof(double));
	ex=(double*)malloc(nelems*sizeof(double));
	ey=(double*)malloc(nelems*sizeof(double));
	ez=(double*)malloc(nelems*sizeof(double));
	id=(h5_int64_t*)malloc(nelems*sizeof(h5_int64_t));

	TEST("Writing 64-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		for (i=0; i<nelems; i++)
		{
			e[i]   = 0.0 + (double)(i+nelems*t);
			ex[i]  = 0.1 + (double)(i+nelems*t);
			ey[i]  = 0.2 + (double)(i+nelems*t);
			ez[i]  = 0.3 + (double)(i+nelems*t);
			id[i] = i + nelems*t;
		}

		val = H5HasStep(file, t);

		status = H5SetStep(file, t);
		RETURN(status, H5_SUCCESS, "H5SetStep");

		if (val == 0) test_write_field_attribs(file, "e", t);

		status = H5Block3dSetView(file,
		                          layout[0], layout[1],
		                          layout[2], layout[3],
		                          layout[4], layout[5]);
		RETURN(status, H5_SUCCESS, "H5Block3dSetView");

		status = H5Block3dWriteScalarFieldFloat64(file, "e", e);
		RETURN(status, H5_SUCCESS, "H5Block3dWriteScalarFieldFloat64");

		status = H5Block3dWriteVector3dFieldFloat64(file,
		                                            "E", ex, ey, ez);
		RETURN(status, H5_SUCCESS,
		       "H5Block3dWriteVector3dFieldFloat64");

		status = H5Block3dWriteScalarFieldInt64(file, "id", id);
		RETURN(status, H5_SUCCESS, "H5Block3dWriteScalarFieldInt64");
	}
}

static void
test_write_data32(h5_file_t file, int step)
{
	int i,j,k,t;
	h5_int64_t status, val;

	float *e;
	float *ex,*ey,*ez;
	int *id;

	size_t nelems = NBLOCKX * (NBLOCKY+2) * (NBLOCKZ+4);

	e=(float*)malloc(nelems*sizeof(double));
	ex=(float*)malloc(nelems*sizeof(double));
	ey=(float*)malloc(nelems*sizeof(double));
	ez=(float*)malloc(nelems*sizeof(double));
	id=(int*)malloc(nelems*sizeof(int));

	nelems = NBLOCKX * NBLOCKY * NBLOCKZ;

	TEST("Writing 32-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		for (k=0; k<NBLOCKZ+4; k++)
			for (j=0; j<NBLOCKY+2; j++)
				for (i=0; i<NBLOCKX; i++)
				{
					int idx = i + j*NBLOCKX + k*NBLOCKX*(NBLOCKY+2);
					int n = i + (j-1)*NBLOCKX + (k-2)*NBLOCKX*NBLOCKY;
					e[idx]   = 0.0f + (float)(n+nelems*t);
					ex[idx]  = 0.1f + (float)(n+nelems*t);
					ey[idx]  = 0.2f + (float)(n+nelems*t);
					ez[idx]  = 0.3f + (float)(n+nelems*t);
					id[idx] = n + nelems*t;
				}

		val = H5HasStep(file, t);

		status = H5SetStep(file, t);
		RETURN(status, H5_SUCCESS, "H5SetStep");

		if (val == 0) test_write_field_attribs(file, "e", t);

#if defined(H5_HAVE_PARALLEL)
		extern h5_size_t grid[3];

		status = H5Block3dSetGrid(file, grid[0], grid[1], grid[2]);
		RETURN(status, H5_SUCCESS, "H5Block3dSetGrid");

		status = H5Block3dSetDims(file, NBLOCKX, NBLOCKY, NBLOCKZ);
		RETURN(status, H5_SUCCESS, "H5Block3dSetDims");
#endif
		status = H5Block3dSetHalo(file, 0, 1, 2);
		RETURN(status, H5_SUCCESS, "H5Block3dSetHalo");

		status = H5Block3dWriteScalarFieldFloat32(file, "e", e);
		RETURN(status, H5_SUCCESS, "H5Block3dWriteScalarFieldFloat32");

		status = H5Block3dWriteVector3dFieldFloat32(file,
		                                            "E", ex, ey, ez);
		RETURN(status, H5_SUCCESS,
		       "H5Block3dWriteVector3dFieldFloat32");

		status = H5Block3dWriteScalarFieldInt32(file, "id", id);
		RETURN(status, H5_SUCCESS, "H5Block3dWriteScalarFieldInt32");
	}
}

void h5b_test_write1(void)
{
	h5_file_t file1;
	h5_err_t status;

	TEST("Opening file once, write-truncate");
        h5_prop_t props = H5CreateFileProp ();
#if defined (H5_HAVE_PARALLEL)
        MPI_Comm comm = MPI_COMM_WORLD;
        status = H5SetPropFileMPIOCollective (props, &comm);
	RETURN(status, H5_SUCCESS, "H5SetPropFileMPIOCollective");
        status = H5SetPropFileThrottle (props, 2);
	RETURN(status, H5_SUCCESS, "H5SetPropFileThrottle");
#endif
	file1 = H5OpenFile(FILENAME, H5_O_WRONLY, props);
	status = H5CheckFile(file1);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

        status = H5CloseProp (props);
	RETURN(status, H5_SUCCESS, "H5CloseProp");

	test_write_data64(file1, 1);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}

void h5b_test_write2(void)
{
	h5_file_t file1;
	h5_file_t file2;

	h5_err_t status;

	TEST("Opening file twice, write-append + read-only");
	file1 = H5OpenFile(FILENAME, H5_O_APPENDONLY, H5_PROP_DEFAULT);

	status = H5CheckFile(file1);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

	file2 = H5OpenFile(FILENAME, H5_O_RDONLY, H5_PROP_DEFAULT);

	status = H5CheckFile(file2);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

	test_write_data32(file1, NTIMESTEPS+1);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
	status = H5CloseFile(file2);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}
