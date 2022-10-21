/* Test framework borrowed from HDF5 1.8.3:
 * test/testframe.c
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Copyright by The HDF Group.                                               *
* Copyright by the Board of Trustees of the University of Illinois.         *
* All rights reserved.                                                      *
*                                                                           *
* This file is part of HDF5.  The full HDF5 copyright notice, including     *
* terms governing use, modification, and redistribution, is contained in    *
* the files COPYING and Copyright.html.  COPYING can be found at the root   *
* of the source code distribution tree; Copyright.html can be found at the  *
* root level of an installed copy of the electronic HDF5 document set and   *
* is linked from the top-level documents page.  It can also be found at     *
* http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
* access to either file, you may request a copy from help@hdfgroup.org.     *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Programmer:  Quincey Koziol <koziol@ncsa.uiuc.edu>
 *              Tuesday, January  6, 2004
 *
 * Purpose:	Provides support functions for the testing framework.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "testframe.h"

/*
 * Definitions for the testing structure.
 */
#define MAXNUMOFTESTS   64
#define MAXTESTNAME     16
#define MAXTESTDESC     64

typedef struct TestStruct {
	int NumErrors;
	char Description[MAXTESTDESC];
	int SkipFlag;
	char Name[MAXTESTNAME];
	void (*Call)(void);
	void (*Cleanup)(void);
	const void *Parameters;
} TestStruct;


/*
 * Variables used by testing framework.
 */
static int num_errs = 0;        /* Total number of errors during testing */
static int Verbosity = VERBO_DEF;       /* Default Verbosity is Low */
static int Summary = 0;         /* Show test summary. Default is no. */
static int CleanUp = 1;         /* Do cleanup or not. Default is yes. */
static int TestExpress = -1;    /* Do TestExpress or not. -1 means not set yet. */
static TestStruct Test[MAXNUMOFTESTS];
static int Index = 0;
static const void *Test_parameters = NULL;
static const char *TestProgName = NULL;
static void (*TestPrivateUsage)(void) = NULL;
static int (*TestPrivateParser)(int ac, char *av[]) = NULL;


/*
 * Setup a test function and add it to the list of tests.
 *      It must have no parameters and returns void.
 * TheName--short test name.
 *    If the name starts with '-', do not run it by default.
 * TheCall--the test routine.
 * Cleanup--the cleanup routine for the test.
 * TheDescr--Long description of the test.
 * Parameters--pointer to extra parameters. Use NULL if none used.
 *    Since only the pointer is copied, the contents should not change.
 */
void
AddTest(const char *TheName, void (*TheCall)(void), void (*Cleanup)(void), const char *TheDescr, const void *Parameters)
{
	/* Sanity checking */
	if (Index >= MAXNUMOFTESTS) {
		printf("Too many tests added, increase MAXNUMOFTEST(%d).\n",
		       MAXNUMOFTESTS);
		exit(-1);
	}                       /* end if */
	if (strlen(TheDescr) >= MAXTESTDESC) {
		printf("Test description too long, increase MAXTESTDESC(%d).\n",
		       MAXTESTDESC);
		exit(-1);
	} /* end if */
	if (strlen(TheName) >= MAXTESTNAME) {
		printf("Test name too long, increase MAXTESTNAME(%d).\n",
		       MAXTESTNAME);
		exit(-1);
	} /* end if */

	/* Set up test function */
	strcpy(Test[Index].Description, TheDescr);
	if (*TheName != '-') {
		strcpy(Test[Index].Name, TheName);
		Test[Index].SkipFlag = 0;
	}
	else {  /* skip test by default */
		strcpy(Test[Index].Name, TheName+1);
		Test[Index].SkipFlag = 1;
	}
	Test[Index].Call = TheCall;
	Test[Index].Cleanup = Cleanup;
	Test[Index].NumErrors = -1;
	Test[Index].Parameters = Parameters;

	/* Increment test count */
	Index++;
}


/*
 * Initialize testing framework
 *
 * ProgName: Name of test program.
 * private_usage: Optional routine provided by test program to print the
 *      private portion of usage page.  Default to NULL which means none is
 *      provided.
 * private_parser: Optional routine provided by test program to parse the
 *      private options.  Default to NULL which means none is provided.
 *
 * Modifications:
 *      Albert Cheng 2004/08/17
 *      Added the ProgName, private_usage and private_parser arguments.
 */
void TestInit(const char *ProgName, void (*private_usage)(void), int (*private_parser)(int ac, char *av[]))
{
#if !(defined MAC || defined SYMANTEC_C)
	/* Un-buffer the stdout and stderr */
	setbuf(stderr, NULL);
	setbuf(stdout, NULL);
#endif

	/*
	 * Turn off automatic error reporting since we do it ourselves.  Besides,
	 * half the functions this test calls are private, so automatic error
	 * reporting wouldn't do much good since it's triggered at the API layer.
	 */
	//H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

	/*
	 * Record the program name and private routines if provided.
	 */
	TestProgName = ProgName;
	if (NULL != private_usage)
		TestPrivateUsage = private_usage;
	if (NULL != private_parser)
		TestPrivateParser = private_parser;
}


/*
 * Print test usage.
 *	First print the common test options, then the extra options if provided.
 *
 * Modification:
 *      2004/08/18 Albert Cheng.  Add TestPrivateUsage feature.
 */
void TestUsage(void)
{
	int i;

	printf("Usage: %s [-v[erbose] (l[ow]|m[edium]|h[igh]|0-9)] %s\n",
	       TestProgName, (TestPrivateUsage ? "<extra options>" : ""));
	printf("              [-[e]x[clude] name+] \n");
	printf("              [-o[nly] name+] \n");
	printf("              [-b[egin] name] \n");
	printf("              [-s[ummary]]  \n");
	printf("              [-c[leanoff]]  \n");
	printf("              [-h[elp]]  \n");
	printf("\n\n");
	printf("verbose   controls the amount of information displayed\n");
	printf("exclude   to exclude tests by name\n");
	printf("only      to name tests which should be run\n");
	printf("begin     start at the name of the test givin\n");
	printf("summary   prints a summary of test results at the end\n");
	printf("cleanoff  does not delete *.hdf files after execution of tests\n");
	printf("help      print out this information\n");
	if (TestPrivateUsage) {
		printf("\nExtra options\n");
		TestPrivateUsage();
	}
	printf("\n\n");
	printf("This program currently tests the following: \n\n");
	printf("%16s %s\n", "Name", "Description");
	printf("%16s %s\n", "----", "-----------");

	for (i = 0; i < Index; i++)
		printf("%16s %s\n", Test[i].Name, Test[i].Description);

	printf("\n\n");
}


/*
 * Print test info.
 */
void TestInfo(const char *ProgName)
{
	unsigned major, minor, release;

	H5get_libversion(&major, &minor, &release);

	printf("\nFor help use: %s -help\n",ProgName);
	printf("Linked with HDF5 version %u.%u release %u\n", major, minor, release);
	printf("Linked with H5hut version %s\n", H5_VER_STRING);
}


/*
 * Parse command line information.
 *      argc, argv: the usual command line argument count and strings
 *
 * Modification:
 *      2004/08/18 Albert Cheng.  Add extra_parse feature.
 */
void TestParseCmdLine(int argc, char *argv[])
{
	int ret_code;

	while (argv++, --argc > 0) {
		if ((strcmp(*argv, "-verbose") == 0) ||
		    (strcmp(*argv, "-v") == 0)) {
			if (argc > 0) {
				--argc; ++argv;
				ParseTestVerbosity(*argv);
			}else{
				TestUsage();
				exit(1);
			}
		}
		else if (((strcmp(*argv, "-exclude") == 0) ||
		          (strcmp(*argv, "-x") == 0))) {
			if (argc > 0) {
				--argc; ++argv;
				SetTest(*argv, SKIPTEST);
			}else{
				TestUsage();
				exit(1);
			}
		}
		else if (((strcmp(*argv, "-begin") == 0) ||
		          (strcmp(*argv, "-b") == 0))) {
			if (argc > 0) {
				--argc; ++argv;
				SetTest(*argv, BEGINTEST);
			}else{
				TestUsage();
				exit(1);
			}
		}
		else if (((strcmp(*argv, "-only") == 0) ||
		          (strcmp(*argv, "-o") == 0))) {
			if (argc > 0) {
				int Loop;
				--argc; ++argv;
				/* Skip all tests, then activate only one. */
				for (Loop = 0; Loop < Index; Loop++)
					Test[Loop].SkipFlag = 1;
				SetTest(*argv, ONLYTEST);
			}else{
				TestUsage();
				exit(1);
			}
		}
		else if ((strcmp(*argv, "-summary") == 0) || (strcmp(*argv, "-s") == 0))
			Summary = 1;
		else if ((strcmp(*argv, "-help") == 0) || (strcmp(*argv, "-h") == 0)) {
			TestUsage();
			exit(0);
		}
		else if ((strcmp(*argv, "-cleanoff") == 0) || (strcmp(*argv, "-c") == 0))
			SetTestNoCleanup();
		else {
			/* non-standard option.  Break out. */
			break;
		}

	}

	/* Call extra parsing function if provided. */
	if (NULL != TestPrivateParser) {
		ret_code=TestPrivateParser(argc+1, argv-1);
		if (ret_code != 0)
			exit(-1);
	}
}


/*
 * Perform Tests.
 */
void PerformTests(void)
{
	int Loop;

	for (Loop = 0; Loop < Index; Loop++)
		if (Test[Loop].SkipFlag) {
			MESSAGE(2, ("Skipping -- %s (%s) \n", Test[Loop].Description, Test[Loop].Name));
		} else {
			MESSAGE(2, ("Testing  -- %s (%s) \n", Test[Loop].Description, Test[Loop].Name));
			MESSAGE(5, ("===============================================\n"));
			Test[Loop].NumErrors = num_errs;
			Test_parameters = Test[Loop].Parameters;
			//ALARM_ON;
			Test[Loop].Call();
			//ALARM_OFF;
			Test[Loop].NumErrors = num_errs - Test[Loop].NumErrors;
			MESSAGE(5, ("===============================================\n"));
			MESSAGE(5, ("There were %d errors detected.\n\n", (int)Test[Loop].NumErrors));
		}

	Test_parameters = NULL; /* clear it. */
	MESSAGE(2, ("\n\n"));

	if (num_errs)
		printf("!!! %d Error(s) were detected !!!\n\n", (int) num_errs);
	else
		printf("All tests were successful. \n\n");
}


/*
 * Display test summary.
 */
void TestSummary(void)
{
	int Loop;

	printf("Summary of Test Results:\n");
	printf("Name of Test     Errors Description of Test\n");
	printf("---------------- ------ --------------------------------------\n");

	for (Loop = 0; Loop < Index; Loop++) {
		if (Test[Loop].NumErrors == -1)
			printf("%16s %6s %s\n", Test[Loop].Name, "N/A", Test[Loop].Description);
		else
			printf("%16s %6d %s\n", Test[Loop].Name, (int)Test[Loop].NumErrors, Test[Loop].Description);
	}

	printf("\n\n");
}


/*
 * Cleanup files from testing
 */
void TestCleanup(void)
{
	int Loop;

	MESSAGE(2, ("\nCleaning Up temp files...\n\n"));

	/* call individual cleanup routines in each source module */
	for (Loop = 0; Loop < Index; Loop++)
		if (!Test[Loop].SkipFlag && Test[Loop].Cleanup!=NULL)
			Test[Loop].Cleanup();
}


/*
 * Retrieve the verbosity level for the testing framework
 */
int GetTestVerbosity(void)
{
	return(Verbosity);
}

/*
 * Set the verbosity level for the testing framework.
 * Return previous verbosity level.
 */
int SetTestVerbosity(int newval)
{
	int oldval;

	oldval = Verbosity;
	Verbosity = newval;
	return(oldval);
}

/*
 * Retrieve the TestExpress mode for the testing framework
   Values:
   0: Exhaustive run
    Tests should take as long as necessary
   1: Full run.  Default if HDF5TestExpress is not defined
    Tests should take no more than 30 minutes
   2: Quick run
    Tests should take no more than 10 minutes
   3: Smoke test.  Default if HDF5TestExpress is set to a value other than 0-3
    Tests should take less than 1 minute

   Design:
   If the environment variable $HDF5TestExpress is defined,
   then test programs should skip some tests so that they
   complete sooner.

   Terms:
   A "test" is a single executable, even if it contains multiple
   sub-tests.
   The standard system for test times is a Linux machine running in
   NFS space (to catch tests that involve a great deal of disk I/O).

   Implementation:
   I think this can be easily implemented in the test library (libh5test.a)
   so that all tests can just call it to check the status of $HDF5TestExpress.
 */
int GetTestExpress(void)
{
	char * env_val;

	/* set it here for now.  Should be done in something like h5test_init(). */
	if(TestExpress==-1)
	{
		env_val = getenv("HDF5TestExpress");

		if(env_val == NULL)
			SetTestExpress(1);
		else if(strcmp(env_val, "0") == 0)
			SetTestExpress(0);
		else if(strcmp(env_val, "1") == 0)
			SetTestExpress(1);
		else if(strcmp(env_val, "2") == 0)
			SetTestExpress(2);
		else
			SetTestExpress(3);
	}

	return(TestExpress);
}

/*
 * Set the TestExpress mode for the testing framework.
 * Return previous TestExpress mode.
 * Values: non-zero means TestExpress mode is on, 0 means off.
 */
int SetTestExpress(int newval)
{
	int oldval;

	oldval = TestExpress;
	TestExpress = newval;
	return(oldval);
}

/*
 * Retrieve Summary request value.
 *     0 means no summary, 1 means yes.
 */
int GetTestSummary(void)
{
	return(Summary);
}

/*
 * Retrieve Cleanup request value.
 *     0 means no Cleanup, 1 means yes.
 */
int GetTestCleanup(void)
{
	return(CleanUp);
}

/*
 * Set cleanup to no.
 * Return previous cleanup value.
 */
int SetTestNoCleanup(void)
{
	int oldval;

	oldval = CleanUp;
	CleanUp = 0;
	return(oldval);
}

/*
 * Parse an argument string for verbosity level and set it.
 */
void ParseTestVerbosity(char *argv)
{
	if (*argv == 'l')
		SetTestVerbosity(VERBO_LO);
	else if (*argv == 'm')
		SetTestVerbosity(VERBO_MED);
	else if (*argv == 'h')
		SetTestVerbosity(VERBO_HI);
	else
		SetTestVerbosity(atoi(argv));
}


/*
 * Retrieve the number of testing errors for the testing framework
 */
int GetTestNumErrs(void)
{
	return(num_errs);
}


/*
 * Increment the number of testing errors
 */
void IncTestNumErrs(void)
{
	num_errs++;
}


/*
 * Retrieve the current Test Parameters pointer.
 */
const void *GetTestParameters(void)
{
	return(Test_parameters);
}

int
TestPrintf(const char *format, ...)
{
	va_list arglist;
	int ret_value = -1;

#if defined(H5_HAVE_PARALLEL)
	int nproc;
	MPI_Comm_rank(MPI_COMM_WORLD, &nproc);
	if ( nproc == 0 || VERBOSE_HI ) {
		char *format2 = malloc(strlen(format)+8);
		sprintf(format2, "[%d] %s", nproc, format);
		va_start(arglist, format);
		ret_value = vprintf(format2, arglist);
		va_end(arglist);
	}
#else
	va_start(arglist, format);
	ret_value = vprintf(format, arglist);
	va_end(arglist);
#endif

	/* Return the length of the string produced (like printf() does) */
	return ret_value;
}


/*
 * This routine is designed to provide equivalent functionality to 'printf'
 * and also increment the error count for the testing framework.
 */
int
TestErrPrintf(const char *format, ...)
{
	va_list arglist;
	int ret_value = -1;

	/* Increment the error count */
	num_errs++;

#if H5_HAVE_PARALLEL
	int nproc;
	MPI_Comm_rank(MPI_COMM_WORLD, &nproc);
	if ( nproc == 0 || VERBOSE_HI ) {
		char *format2 = malloc(strlen(format)+8);
		sprintf(format2, "[%d] %s", nproc, format);
		va_start(arglist, format);
		ret_value = vfprintf(stderr, format2, arglist);
		va_end(arglist);
	}
#else
	va_start(arglist, format);
	ret_value = vfprintf(stderr, format, arglist);
	va_end(arglist);
#endif

	/* Return the length of the string produced (like printf() does) */
	return ret_value;
}


/*
 * Set (control) which test will be tested.
 * SKIPTEST: skip this test
 * ONLYTEST: do only this test
 * BEGINETEST: skip all tests before this test
 *
 */
void SetTest(const char *testname, int action)
{
	int Loop;
	switch (action) {
	case SKIPTEST:
		for (Loop = 0; Loop < Index; Loop++)
			if (strcmp(testname, Test[Loop].Name) == 0) {
				Test[Loop].SkipFlag = 1;
				break;
			}
		break;
	case BEGINTEST:
		for (Loop = 0; Loop < Index; Loop++) {
			if (strcmp(testname, Test[Loop].Name) != 0)
				Test[Loop].SkipFlag = 1;
			else{
				/* Found it. Set it to run.  Done. */
				Test[Loop].SkipFlag = 0;
				break;
			}
		}
		break;
	case ONLYTEST:
		for (Loop = 0; Loop < Index; Loop++) {
			if (strcmp(testname, Test[Loop].Name) != 0)
				Test[Loop].SkipFlag = 1;
			else {
				/* Found it. Set it to run. Break to skip the rest. */
				Test[Loop].SkipFlag = 0;
				break;
			}
		}
		/* skip the rest */
		while (++Loop < Index)
			Test[Loop].SkipFlag = 1;
		break;
	default:
		/* error */
		printf("*** ERROR: Unknown action (%d) for SetTest\n", action);
		break;
	}
}

void
get_attr_name(char *name, char *tag, int id)
{
	sprintf(name, "Attr%d%s", id, tag);
}

void
test_open_objects(h5_file_t file, int max_objects)
{
	hid_t hfile = h5_get_hdf5_file(file);
	ssize_t nopen = H5Fget_obj_count(hfile, H5F_OBJ_ALL);
	if (nopen > max_objects)
	{
		TestErrPrintf(  "*** TOO MANY OBJECTS OPEN: %d > %d "
		                "at line %4d in %s\n", nopen, max_objects,
		                (int)__LINE__, __FILE__ );

		hid_t *list = malloc(sizeof(hid_t)*nopen);
		H5Fget_obj_ids(hfile, H5F_OBJ_ALL, nopen, list);

		H5O_info_t info;
		int i;
		for (i=0; i<nopen; i++) {
			H5Oget_info(list[i], &info);
			switch (info.type) {
			case H5O_TYPE_GROUP:
				TestErrPrintf("obj%d has type GROUP\n", i);
				break;
			case H5O_TYPE_DATASET:
				TestErrPrintf("obj%d has type DATASET\n", i);
				break;
			case H5O_TYPE_NAMED_DATATYPE:
				TestErrPrintf("obj%d has type NAMED_DATATYPE\n", i);
				break;
			default:
				TestErrPrintf("obj%d has unknown type\n", i);
			}
		}
		free(list);
	}
}

