#ifndef _H5HUT_TEST_PARAMS_H_
#define _H5HUT_TEST_PARAMS_H_

#define FILENAME "test.h5"
#define LONGNAME "thisisaverylongnamethatshouldexceedthelimitof64charcausingawarningtoprint"
#define LONGNAME2 "thisisaverylongnamethatshouldexceedthelimitof64charcausingawarni"
#define NTIMESTEPS 10

/* do not decrease this value below 99, or it will break assumptions
 * made in the read tests! */
#define NPARTICLES 99

#define NBLOCKX 16
#define NBLOCKY 12
#define NBLOCKZ 18

/* do not increase this value past 32! */
#define MAX_MPI_TASKS 32

#define ATTR_NAME_SIZE 16
#define ATTR_STR_VAL   "test"
#define ATTR_INT32_VAL -2147483648
#define ATTR_INT64_VAL 2147483648
#define ATTR_FLOAT_VAL 3.14159265F

#endif

