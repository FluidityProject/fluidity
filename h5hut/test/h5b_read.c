#include <stdlib.h>
#include <string.h>
#include "testframe.h"
#include "params.h"

static void
test_read_field_attribs(
	h5_file_t file,
	const char *field_name,
	int position)
{
	h5_err_t status;
	char name[ATTR_NAME_SIZE];
	char str[ATTR_NAME_SIZE];
	h5_int32_t i32;
	h5_int64_t i64;
	h5_float32_t f32;
	h5_float64_t f64;

	TEST("Reading field attributes");

	i64 = H5BlockGetNumFieldAttribs(file, field_name);
	VALUE(i64 % 5, 0, "file attribute count");

	get_attr_name(name, "str", position);
	status = H5BlockReadFieldAttribString(
	        file, field_name, name, str);
	RETURN(status, H5_SUCCESS, "H5BlockReadFieldAttribString");
	SVALUE(str, ATTR_STR_VAL, "string attribute");

	get_attr_name(name, "i32", position);
	status = H5BlockReadFieldAttribInt32(
	        file, field_name, name, &i32);
	RETURN(status, H5_SUCCESS, "H5BlockReadFieldAttribInt32");
	IVALUE(i32, ATTR_INT32_VAL, "int32 attribute");

	get_attr_name(name, "i64", position);
	status = H5BlockReadFieldAttribInt64(
	        file, field_name, name, &i64);
	RETURN(status, H5_SUCCESS, "H5BlockReadFieldAttribInt64");
	IVALUE(i64, ATTR_INT64_VAL, "int64 attribute");

	get_attr_name(name, "f32", position);
	status = H5BlockReadFieldAttribFloat32(
	        file, field_name, name, &f32);
	RETURN(status, H5_SUCCESS, "H5BlockReadFieldAttribFloat32");
	FVALUE(f32, ATTR_FLOAT_VAL, "float32 attribute");

	get_attr_name(name, "f64", position);
	status = H5BlockReadFieldAttribFloat64(
	        file, field_name, name, &f64);
	RETURN(status, H5_SUCCESS, "H5BlockReadFieldAttribFloat64");
	FVALUE(f64, ATTR_FLOAT_VAL, "float64 attribute");
}

static void
test_read_data64(h5_file_t file, int step)
{
	extern h5_size_t grid[3];
	extern h5_size_t layout[6];

	int i,t;
	h5_err_t status;
	h5_int64_t val, type[2];
	char name[4];
	h5_size_t field_rank[2], field_dims[6], elem_rank[2];

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

	TEST("Verifying dataset info");

#if defined(H5_HAVE_PARALLEL)
	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif

	status = H5SetStep(file, step);
	RETURN(status, H5_SUCCESS, "H5SetStep");

	val = H5BlockGetNumFields(file);
	IVALUE(val, 3, "field count");

	for (i=0; i<3; i++) {
		status = H5BlockGetFieldInfo(
		        file, i, name, 4,
		        field_rank, field_dims, elem_rank, type);
		RETURN(status, H5_SUCCESS, "H5BlockGetFieldInfo");

		status = H5BlockGetFieldInfoByName(
		        file, name,
		        field_rank+1, field_dims+3, elem_rank+1, type+1);
		RETURN(status, H5_SUCCESS, "H5BlockGetFieldInfoByName");
		IVALUE(field_rank[0], field_rank[1], "field rank");
		IVALUE(field_dims[0], field_dims[3], "field dims x");
		IVALUE(field_dims[1], field_dims[4], "field dims y");
		IVALUE(field_dims[2], field_dims[5], "field dims z");
		IVALUE(elem_rank[0], elem_rank[1], "elem rank");
		IVALUE(type[0], type[1], "field type");

		IVALUE(field_rank[0], 3, "field rank");
		IVALUE(field_dims[0], grid[0]*NBLOCKX, "field dims x");
		IVALUE(field_dims[1], grid[1]*NBLOCKY, "field dims y");
		IVALUE(field_dims[2], grid[2]*NBLOCKZ, "field dims z");
		if (i==1) {
			CVALUE(name[0], 'e', "field name");
			IVALUE(elem_rank[0], 1, "elem rank");
			IVALUE(type[0], H5_FLOAT64_T, "field type");
		} else if (i==0) {
			CVALUE(name[0], 'E', "field name");
			IVALUE(elem_rank[0], 3, "elem rank");
			IVALUE(type[1], H5_FLOAT64_T, "field type");
		} else if (i==2) {
			CVALUE(name[0], 'i', "field name");
			IVALUE(elem_rank[0], 1, "elem rank");
			IVALUE(type[1], H5_INT64_T, "field type");
		}
	}

	TEST("Reading 64-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		val = H5HasStep(file, t);
		IVALUE(val, 1, "has step");

		status = H5SetStep(file, t);
		RETURN(status, H5_SUCCESS, "H5SetStep");

		test_read_field_attribs(file, "e", t);

		status = H5Block3dSetView(file,
		                          layout[0], layout[1],
		                          layout[2], layout[3],
		                          layout[4], layout[5]);
		RETURN(status, H5_SUCCESS, "H5Block3dSetView");

		status = H5Block3dReadScalarFieldFloat64(file, "e", e);
		RETURN(status, H5_SUCCESS, "H5Block3dReadScalarFieldFloat64");

		status = H5Block3dReadVector3dFieldFloat64(file,
		                                           "E", ex, ey, ez);
		RETURN(status, H5_SUCCESS,
		       "H5Block3dReadVector3dFieldFloat64");

		status = H5Block3dReadScalarFieldInt64(file, "id", id);
		RETURN(status, H5_SUCCESS, "H5Block3dReadScalarFieldInt64");

		for (i=0; i<nelems; i++)
		{
			FVALUE(e[i], 0.0 + (double)(i+nelems*t), " e data");
			FVALUE(ex[i], 0.1 + (double)(i+nelems*t), " ex data");
			FVALUE(ey[i], 0.2 + (double)(i+nelems*t), " ey data");
			FVALUE(ez[i], 0.3 + (double)(i+nelems*t), " ez data");
			IVALUE(id[i],               (i+nelems*t), " id data");
		}
	}
}

static void
test_read_data32(h5_file_t file, int step)
{
	extern h5_size_t grid[3];

	int t;
	h5_err_t status;

	float *e;
	float *ex,*ey,*ez;
	int *id;

	const size_t nelems = NBLOCKX * NBLOCKY * NBLOCKZ;

	e=(float*)malloc(nelems*sizeof(double));
	ex=(float*)malloc(nelems*sizeof(double));
	ey=(float*)malloc(nelems*sizeof(double));
	ez=(float*)malloc(nelems*sizeof(double));
	id=(int*)malloc(nelems*sizeof(int));

	TEST("Reading 32-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		h5_int64_t val;
		val = H5HasStep(file, t);
		IVALUE(val, 1, "has step");

		status = H5SetStep(file, t);
		RETURN(status, H5_SUCCESS, "H5SetStep");

		test_read_field_attribs(file, "e", t);
#if defined(H5_HAVE_PARALLEL)
		status = H5Block3dSetGrid(file, grid[0], grid[1], grid[2]);
		RETURN(status, H5_SUCCESS, "H5Block3dSetGrid");

		status = H5Block3dSetDims(file, NBLOCKX, NBLOCKY, NBLOCKZ);
		RETURN(status, H5_SUCCESS, "H5Block3dSetDims");
#endif
		status = H5Block3dReadScalarFieldFloat32(file, "e", e);
		RETURN(status, H5_SUCCESS, "H5Block3dReadScalarFieldFloat32");

		status = H5Block3dReadVector3dFieldFloat32(file,
		                                           "E", ex, ey, ez);
		RETURN(status, H5_SUCCESS,
		       "H5Block3dReadVector3dFieldFloat32");

		status = H5Block3dReadScalarFieldInt32(file, "id", id);
		RETURN(status, H5_SUCCESS, "H5Block3dReadScalarFieldInt32");

		int i;
		for (i=0; i<nelems; i++)
		{
			FVALUE(e[i], 0.0f + (float)(i+nelems*t), " e data");
			FVALUE(ex[i], 0.1f + (float)(i+nelems*t), " ex data");
			FVALUE(ey[i], 0.2f + (float)(i+nelems*t), " ey data");
			FVALUE(ez[i], 0.3f + (float)(i+nelems*t), " ez data");
			IVALUE(id[i],              (i+nelems*t), " id data");
		}
	}
}

void h5b_test_read1(void)
{
	h5_file_t file1;
	h5_err_t status;

	TEST("Opening file once, read-only");
        h5_prop_t props = H5CreateFileProp ();
#if defined(H5_HAVE_PARALLEL)
        MPI_Comm comm = MPI_COMM_WORLD;
        status = H5SetPropFileMPIOCollective (props, &comm);
	RETURN(status, H5_SUCCESS, "H5SetPropFileMPIOCollective");
        status = H5SetPropFileThrottle (props, 2);
	RETURN(status, H5_SUCCESS, "H5SetPropFileThrottle");
#endif
	file1 = H5OpenFile(FILENAME, H5_O_RDONLY, props);
	status = H5CheckFile(file1);
	RETURN(status, H5_SUCCESS, "H5CheckFile");

        status = H5CloseProp (props);
	RETURN(status, H5_SUCCESS, "H5CloseProp");

	test_read_data64(file1, 1);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}

void h5b_test_read2(void)
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

	test_read_data32(file1, NTIMESTEPS+1);

	status = H5CloseFile(file1);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
	status = H5CloseFile(file2);
	RETURN(status, H5_SUCCESS, "H5CloseFile");
}
