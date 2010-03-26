#include "fmangle.h"
#include "../FluxesReader.h"

using namespace std;

extern "C" {
#define test_FluxesReader_fc F77_FUNC(test_FluxesReader, TEST_FLUXESREADER)  
  void test_FluxesReader_fc();
}

extern void report_test(const string& title, const bool& fail, const bool& warn, const string& msg);

void test_FluxesReader_fc() {

#ifdef HAVE_LIBUDUNITS
  FluxesReader data;
  FluxesReader data2;
  data.VerboseOff();
  const int nFields = 1;
  double values[nFields];
  double value;
  double correct, correct_t, correct_d;
  int err = -1;
  bool fail = true;
  bool warn = false;
  char errorMessage[256];

  // lets first test the global file
  data.RegisterDataFile("../../tests/data/global_fluxes.nc");
  data.AddFieldOfInterest(string("_2t"));
  data.SetSimulationTimeUnits("seconds since 1960-01-01 06:00:0.0");
  data.SetTimeSeconds(0);

  // pick a few random points and grab the temperature
  // these should equal those given by ncview
  // point 1
  correct = 299.237;
  value = data.GetScalar("_2t",210.0,10.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: single point 1]",fail,warn,errorMessage);
  fail = true;
  // point 1a - useful for debugging
  correct = 298.730;
  value = data.GetScalar("_2t",212.5,10.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: single point 1a]",fail,warn,errorMessage);
  fail = true;
  // point 2
  correct = 299.120;
  value = data.GetScalar("_2t",357.5,0.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: single point 2]",fail,warn,errorMessage);
  fail = true;
  // point 3
  correct = 244.405;
  value = data.GetScalar("_2t",0.0,90.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: single point 3]",fail,warn,errorMessage);
  fail = true;


  // now we move onto interpolation
  // These answers are based using scipy to read in netcdf and interpolate 
  // The script is also included in this directory for completeness
  //
  // interpolate in longitude. 
  // Test near point 1
  correct = 299.217;
  value = data.GetScalar("_2t",210.1,10.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: longitude interp 1]",fail,warn,errorMessage);
  fail = true;
  // Test near point 2
  correct = 298.832;
  value = data.GetScalar("_2t",212.0,10.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: longitude interp 2]",fail,warn,errorMessage);
  fail = true;
  // Test half-way between
  correct = 298.984;
  value = data.GetScalar("_2t",211.25,10.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: longitude interp 3]",fail,warn,errorMessage);
  fail = true;  
  // Test the wrap-over across 0
  correct = 287.891;
  value = data.GetScalar("_2t",359,10.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: longitude interp over 0 deg]",fail,warn,errorMessage);
  fail = true; 

  // Do the same with latitude
  // Test near point 1
  correct = 299.070;
  value = data.GetScalar("_2t",210.,10.5);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: latitude interp 1]",fail,warn,errorMessage);
  fail = true; 
  // Test near point 2
  correct = 298.569;
  value = data.GetScalar("_2t",210.,12.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: latitude interp 2]",fail,warn,errorMessage);
  fail = true; 
  // Test half-way between
  correct = 298.819;
  value = data.GetScalar("_2t",210.,11.25);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: latitude interp 3]",fail,warn,errorMessage);
  fail = true; 

  // Now interpolate in lat and long
  // Near one corner
  correct = 298.988;
  value = data.GetScalar("_2t",210.5,10.5);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: lat-long interp 1]",fail,warn,errorMessage);
  fail = true; 
  // in the middle of the grid
  correct = 298.688;
  value = data.GetScalar("_2t",211.25,11.25);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: lat-long interp 2]",fail,warn,errorMessage);
  fail = true; 
  // in the middle, but to one side
  correct = 298.475;
  value = data.GetScalar("_2t",212.,12.);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: lat-long interp 3]",fail,warn,errorMessage);
  fail = true; 
  // Trying the wrap over at 0 again, but not on a latitude data point
  correct = 289.33535;
  value = data.GetScalar("_2t",359.,12.);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: lat-long interp over 0]",fail,warn,errorMessage);
  fail = true; 

  // Time for some temporal interpolation now
  // Near time 1
  data.SetTimeSeconds(3600);  // 1 hour in
  correct = 299.209;
  value = data.GetScalar("_2t",210.0,10.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: Temporal interp 1]",fail,warn,errorMessage);
  fail = true;
  // Near time 2
  data.SetTimeSeconds(18000);  // 5 hours in
  correct = 299.095;
  value = data.GetScalar("_2t",210.0,10.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: Temporal interp 2]",fail,warn,errorMessage);
  fail = true;  
  // Half-way between
  data.SetTimeSeconds(10800);  // 3 hours in
  correct = 299.152;
  value = data.GetScalar("_2t",210.0,10.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: Temporal interp 3]",fail,warn,errorMessage);
  fail = true; 

  // If the above work, then these should...
  // Lat, long, temporal interp
  // A corner, near time 1
  data.SetTimeSeconds(3600);  // 1 hour in
  correct = 298.959;
  value = data.GetScalar("_2t",210.5,10.5);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: All dims interp 1]",fail,warn,errorMessage);
  fail = true;  
  // Same corner, near time 2
  data.SetTimeSeconds(18000);  // 5 hours in
  correct = 298.843;
  value = data.GetScalar("_2t",210.5,10.5);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: All dims interp 2]",fail,warn,errorMessage);
  fail = true; 
  // in the middle, temporally and spatially
  data.SetTimeSeconds(10800);  // 3 hours in
  correct = 298.587;
  value = data.GetScalar("_2t",211.25,11.25);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader global: All dims interp 3]",fail,warn,errorMessage);
  fail = true; 

  
  // now we concentrate on getting arrays of values - only need to do one of each of these tests
  // but there are two values to check
  data.AddFieldOfInterest(string("_2d"));
  data.SetTimeSeconds(0);

  // point 1
  correct_t = 299.237;
  correct_d = 296.316;
  err = data.GetScalars(210.0,10.0,values);
  if (err == 0) fail = false;
  report_test("[test_FluxesReader global: return value OK from GetScalars]",fail,warn,"Incorrect return value for GetScalars");
  fail = true;
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader global: GetScalars single point: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader global: GetScalars single point: d]",fail,warn,errorMessage);
  fail = true;

  // point 2 - interp with lat
  correct_t = 298.819;
  correct_d = 295.855;
  err = data.GetScalars(210.0,11.25,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader global: GetScalars interp lat: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader global: GetScalars interp lat: d]",fail,warn,errorMessage);
  fail = true;

  // point 3 - interp with long
  correct_t = 298.984;
  correct_d = 296.061;
  err = data.GetScalars(211.25,10.0,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader global: GetScalars interp long: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader global: GetScalars interp long: d]",fail,warn,errorMessage);
  fail = true;

  // point 3 - interp with long/lat
  correct_t = 298.688;
  correct_d = 295.632;
  err = data.GetScalars(211.25,11.25,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader global: GetScalars interp long-lat: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader global: GetScalars interp long-lat: d]",fail,warn,errorMessage);
  fail = true;

  // point 3 - interp temporally
  data.SetTimeSeconds(10800);  // 3 hours in
  correct_t = 299.152;
  correct_d = 296.191;
  err = data.GetScalars(210.0,10.0,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader global: GetScalars interp time: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[0]); 
  report_test("[test_FluxesReader global: GetScalars interp time: d]",fail,warn,errorMessage);
  fail = true;

  // point n - interp in all 3 dims
  correct_t = 298.587;
  correct_d = 295.614;
  err = data.GetScalars(211.25,11.25,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader global: GetScalars interp all: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader global: GetScalars interp all: d]",fail,warn,errorMessage);
  fail = true;

  // that's it for the ones that should be correct, now move onto the edge cases
  // First, ask for a time that isn't in the file
  //err = data.SetTimeSeconds(2592000);  // 30 days in and the test file is only 2 days long...
  if (err == 0) fail = true; // this *should* return error, so fail if no error returned
  //report_test("[test_FluxesReader: test time error]",fail,warn,"Managed to get a time > than maximum in file");
  fail = true;

  // now ask for an unrealistic lat long combo
  data.SetTimeSeconds(0);
  correct = 299.24;
  //value = data.GetScalar("_2t",575.0,100.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  //report_test("[test_FluxesReader global: dodgy lat/long]",fail,warn,errorMessage);
  fail = true;


  /****************************************
   *  Now testing a subset of data - this
   *  does not contained the whole globe
   *  and exercises different parts of the
   *  code.
   *  This data is Station Papa location -
   *  or close enough
   ****************************************/
  // This tests switching files too
  //err = data.RegisterDataFile("../../tests/data/subset_fluxes.nc"); // should fail
  if (err == -1) fail = false; 
  //report_test("[test_FluxesReader: switch files]",fail,warn,"Failed to error when switching files");
  fail = true;

  err = data2.RegisterDataFile("../../tests/data/subset_fluxes.nc");
  data2.AddFieldOfInterest(string("_2t"));
  data2.SetSimulationTimeUnits("seconds since 1960-01-01 06:00:0.0");
  data2.SetTimeSeconds(0);

  // point 1
  correct = 279.30975;
  value = data2.GetScalar("_2t",215,50.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader subset: single point 1]",fail,warn,errorMessage);
  fail = true;

  // now we move onto interpolation
  // These answers are based using scipy to read in netcdf and interpolate 
  // The script is also included in this directory for completeness
  //
  // interpolate in longitude. 
  correct = 279.5136;
  value = data2.GetScalar("_2t",215.5,50.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader subset: longitude interp]",fail,warn,errorMessage);
  fail = true;
  // Do the same with latitude
  correct = 279.139;
  value = data2.GetScalar("_2t",215.,50.5);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader subset: latitude interp]",fail,warn,errorMessage);
  fail = true; 
  // Now interpolate in lat and long
  correct = 279.340;
  value = data2.GetScalar("_2t",215.5,50.5);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader subset: lat-long interp]",fail,warn,errorMessage);
  fail = true; 
  // Time for some temporal interpolation now
  // Near time 1
  data2.SetTimeSeconds(3600);  // 1 hour in
  correct = 279.298;
  value = data2.GetScalar("_2t",215.0,50.0);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader subset: Temporal interp]",fail,warn,errorMessage);
  fail = true;
  // in the middle, temporally and spatially
  data2.SetTimeSeconds(10800);  // 3 hours in
  correct = 279.213;
  value = data2.GetScalar("_2t",215.5,50.5);
  if (abs(value-correct) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct,value); 
  report_test("[test_FluxesReader subset: Interp all dims]",fail,warn,errorMessage);
  fail = true; 

  
  // now we concentrate on getting arrays of values - only need to do one of each of these tests
  // but there are two values to check
  data2.AddFieldOfInterest(string("_2d"));
  data2.SetTimeSeconds(0);

  // point 1
  correct_t = 279.30975;
  correct_d = 278.867;
  err = data2.GetScalars(215.0,50.0,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader subset: GetScalars single point: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader subset: GetScalars single point: d]",fail,warn,errorMessage);
  fail = true;
  // point 1
  correct_t = 278.458;
  correct_d = 277.210;
  err = data2.GetScalars(215.0,52.5,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader subset: GetScalars single point: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader subset: GetScalars single point: d]",fail,warn,errorMessage);
  fail = true;

  // point 2 - interp with lat
  correct_t = 279.139;
  correct_d = 278.535;
  err = data2.GetScalars(215.0,50.5,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader subset: GetScalars interp lat: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader subset: GetScalars interp lat: d]",fail,warn,errorMessage);
  fail = true;

  // point 3 - interp with long
  correct_t = 279.514;
  correct_d = 279.081;
  err = data2.GetScalars(215.5,50.0,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader subset: GetScalars interp long: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader subset: GetScalars interp long: d]",fail,warn,errorMessage);
  fail = true;

  // point 3 - interp with long/lat
  correct_t = 279.340;
  correct_d = 278.773;
  err = data2.GetScalars(215.5,50.5,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader subset: GetScalars interp long-lat: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader subset: GetScalars interp long-lat: d]",fail,warn,errorMessage);
  fail = true;

  // point 3 - interp temporally
  data2.SetTimeSeconds(10800);  // 3 hours in
  correct_t = 279.275;
  correct_d = 278.804;
  err = data2.GetScalars(215.0,50.0,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader subset: GetScalars interp time: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader subset: GetScalars interp time: d]",fail,warn,errorMessage);
  fail = true;

  // point n - interp in all 3 dims
  correct_t = 279.213;
  correct_d = 278.659;
  err = data2.GetScalars(215.5,50.5,values);
  if (abs(values[0]-correct_t) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_t,values[0]); 
  report_test("[test_FluxesReader subset: GetScalars interp all: t]",fail,warn,errorMessage);
  fail = true;
  if (abs(values[1]-correct_d) < 0.001 ) fail = false;
  sprintf(errorMessage, "Expected %f, got %f",correct_d,values[1]); 
  report_test("[test_FluxesReader subset: GetScalars interp all: d]",fail,warn,errorMessage);
  fail = true;

  // some edge cases
  // ask for lat/long not in the data set
  // Too big long
  //err = data2.GetScalars(255.0,50.,values);
  if (err == 0) fail = true;
  //report_test("[test_FluxesReader subset: Location not in file - lon]",fail,warn,"Should have got error");
  fail = true;
  // too big lat
  //err = data2.GetScalars(215.0,0.,values);
  if (err == 0) fail = true;
  //report_test("[test_FluxesReader subset: Location not in file - lat]",fail,warn,"Should have got error");
  fail = true;

  // handle missing data
  data2.SetTimeSeconds(10800);  // 30 hours in - missing data here
  //err = data2.GetScalars(215.0,50.0,values);
  if ( err == -1 ) fail = false;
  //report_test("[test_FluxesReader subset: Identify missing values]",fail,warn,"Failed to spot missing values");
  fail = true;
#else
  report_test("[dummy]",false,false,"Dummy");
#endif

}
