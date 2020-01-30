/* Test framework borrowed from HDF5 1.8.3:
 * test/h5test.h
 * test/testhdf5.h
 */

#ifndef _H5HUT_TESTFRAME_H_
#define _H5HUT_TESTFRAME_H_

#include <H5hut.h>

/*
 * Predefined test verbosity levels.
 *
 * Convention:
 *
 * The higher the verbosity value, the more information printed.
 * So, output for higher verbosity also include output of all lower
 * verbosity.
 *
 *  Value     Description
 *  0         None:   No informational message.
 *  1                 "All tests passed"
 *  2                 Header of overall test
 *  3         Default: header and results of individual test
 *  4
 *  5         Low:    Major category of tests.
 *  6
 *  7         Medium: Minor category of tests such as functions called.
 *  8
 *  9         High:   Highest level.  All information.
 */
#define VERBO_DEF  0     /* Default */
#define VERBO_LO   2     /* Low     */
#define VERBO_MED  3     /* Medium  */
#define VERBO_HI   65535     /* High    */

/*
 * Verbose queries
 * Only None needs an exact match.  The rest are at least as much.
 */
#define VERBOSE_DEF     (GetTestVerbosity()>=VERBO_DEF)
#define VERBOSE_LO      (GetTestVerbosity()>=VERBO_LO)
#define VERBOSE_MED     (GetTestVerbosity()>=VERBO_MED)
#define VERBOSE_HI      (GetTestVerbosity()>=VERBO_HI)

#define SKIPTEST  1
#define ONLYTEST  2
#define BEGINTEST 3

/*
 * Print the current location on the standard output stream.
 */
#define AT()  TestPrintf ("	 at %s:%d in %s()...\n", \
                          __FILE__, __LINE__, __FUNCTION__);

/*
 * The name of the test is printed by saying TESTING("something") which will
 * result in the string `Testing something' being flushed to standard output.
 * If a test passes, fails, or is skipped then the PASSED(), H5_FAILED(), or
 * SKIPPED() macro should be called.  After H5_FAILED() or SKIPPED() the caller
 * should print additional information to stdout indented by at least four
 * spaces.  If the h5_errors() is used for automatic error handling then
 * the H5_FAILED() macro is invoked automatically when an API function fails.
 */
#define TESTING(WHAT)   {TestPrintf("Testing %-62s",WHAT); fflush(stdout); }
#define TESTING_2(WHAT) {TestPrintf(" Testing %-62s",WHAT); fflush(stdout); }
#define PASSED()        {TestPrintf(" PASSED"); fflush(stdout); }
#define H5_FAILED()     {TestPrintf("*FAILED*"); fflush(stdout); }
#define H5_WARNING()    {TestPrintf("*WARNING*"); fflush(stdout); }
#define SKIPPED()       {TestPrintf(" -SKIP-"); fflush(stdout); }
#define TEST_ERROR      {H5_FAILED(); AT(); goto error; }
#define FAIL_PUTS_ERROR(s) {H5_FAILED(); AT(); TestPrintf(s); goto error; }

/* Use %ld to print the value because long should cover most cases. */
/* Used to make certain a return value _is_not_ a value */
#define RETURN(ret, val, where) do {                                          \
		if (VERBOSE_HI) TestPrintf( "   Call to routine %15s at line %4d "        \
			                    "in %s returned %ld \n",                      \
			                    where, (int)__LINE__, __FILE__,               \
			                    (long)(ret));                                 \
		if ((ret) != (val)) {                                                     \
			TestErrPrintf("*** UNEXPECTED RETURN from %s is %ld at line %4d "     \
			              "in %s\n", where, (long)(ret), (int)__LINE__, __FILE__);   \
		}                                                                         \
} while(0)

#define VALUE(val, expected, what) do {                                 \
		if ((val) != (expected)) {                                      \
			TestErrPrintf(  "*** INCORRECT VALUE of %s at line "    \
			                "%d in %s\n", what, (int)__LINE__,      \
			                __FILE__);                              \
		}                                                               \
} while(0)

#define IVALUE(val, expected, what) do {                                \
		if (VERBOSE_HI) TestPrintf( "   Value of int    %15s at line %4d "   \
			                    "in %s is %ld =? %ld\n", what,          \
			                    (int)__LINE__, __FILE__,                \
			                    (long)(val), (long)(expected));         \
		VALUE(val, expected, what);                                         \
} while(0)

#define FVALUE(val, expected, what) do {                                \
		if (VERBOSE_HI) TestPrintf( "   Value of float  %15s at line %4d "  \
			                    "in %s is %g =? %g\n", what,            \
			                    (int)__LINE__, __FILE__,                \
			                    (val), (expected));                     \
		VALUE(val, expected, what);                                         \
} while(0)

#define CVALUE(val, expected, what) do {                                \
		if (VERBOSE_HI) TestPrintf( "   Value of char   %15s at line %4d "  \
			                    "in %s is %c =? %c\n", what,            \
			                    (int)__LINE__, __FILE__,                \
			                    (val), (expected));                     \
		VALUE(val, expected, what);                                         \
} while(0)


#define SVALUE(val, expected, what) do {                                \
		if (VERBOSE_HI) TestPrintf( "   Value of string %15s at line %4d "  \
			                    "in %s is %s =? %s\n", what,            \
			                    (int)__LINE__, __FILE__,                \
			                    (val), (expected));                     \
		if (strcmp(val,expected) != 0) {                                \
			TestErrPrintf(  "*** INCORRECT VALUE of %d at line "    \
			                "%4d in %s\n", what, (int)__LINE__,     \
			                __FILE__);                              \
		}                                                               \
} while(0)

/* Used to document process through a test */
#define MESSAGE(V,A) do {       \
		if (V) TestPrintf A;    \
} while(0)

#define TEST(WHAT)      MESSAGE(VERBOSE_DEF,(WHAT "\n"))

/* definitions for command strings */
#define VERBOSITY_STR   "Verbosity"
#define SKIP_STR        "Skip"
#define TEST_STR        "Test"
#define CLEAN_STR       "Cleanup"

void TestUsage(void);
void AddTest(const char *TheName, void (*TheCall)(void),
             void (*Cleanup)(void), const char *TheDescr,
             const void *Parameters);
void TestInfo(const char *ProgName);
void TestParseCmdLine(int argc, char *argv[]);
void PerformTests(void);
void TestSummary(void);
void TestCleanup(void);
void TestInit(const char *ProgName, void (*private_usage)(void), int (*private_parser)(int ac, char *av[]));
int  GetTestVerbosity(void);
int  SetTestVerbosity(int newval);
int  GetTestSummary(void);
int  GetTestCleanup(void);
int  SetTestNoCleanup(void);
int  GetTestExpress(void);
int  SetTestExpress(int newval);
void ParseTestVerbosity(char *argv);
int  GetTestNumErrs(void);
void  IncTestNumErrs(void);
const void *GetTestParameters(void);
int  TestPrintf(const char *format, ...);
int  TestErrPrintf(const char *format, ...);
void SetTest(const char *testname, int action);

void
get_attr_name(char *name, char *tag, int id);

void
test_open_objects(h5_file_t file, int max_objects);

#endif

