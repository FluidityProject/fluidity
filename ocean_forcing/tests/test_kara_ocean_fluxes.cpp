#include "fmangle.h"
#include "../FluxesReader.h"
#include "../BulkForcing.h"
#include <fstream>

using namespace std;


extern "C" {
#define test_kara_ocean_fluxes_fc F77_FUNC(test_kara_ocean_fluxes, TEST_KARA_OCEAN_FLUXES)  
    void test_kara_ocean_fluxes_fc();
}

extern int projections(int nPoints, double *x, double *y, double *z, string current_coord, string output_coord);
extern void report_test(const string& title, const bool& fail, const bool& warn, const string& msg);

void test_kara_ocean_fluxes_fc() {
  
   
#ifdef HAVE_LIBUDUNITS
    FluxesReader FluxesReader_ERAdata;
    int err = -1;
    bool fail = true;
    bool warn = false;
    char errorMessage[256];
    int NNodes = 1;
    double error;
    double time = 86400;
    vector<double> T(NNodes, 0.0), S(NNodes, 35.0), Vx(NNodes, 0.0), Vy(NNodes, 0.0), Vz(NNodes, 0.0),
    F_as(NNodes, 0.0), Q_as(NNodes, 0.0), tau_u(NNodes, 0.0), tau_v(NNodes, 0.0),
    X(NNodes, 0.0), Y(NNodes, 0.0), Z(NNodes, 0.0), Q_solar(NNodes, 0.0), e(NNodes,0.0),
    q_l(NNodes,0.0), q_h(NNodes,0.0), q_e(NNodes, 0.0), q_p(NNodes, 0.0), x(NNodes,0.0), 
    y(NNodes,0.0), z(NNodes,0.0);

    // set up the normal input file (ERA40)
    FluxesReader_global.RegisterDataFile("../../tests/data/stationPapa_1970.nc");
    FluxesReader_global.AddFieldOfInterest("u10");  //  0   | 10 metre U wind component
    FluxesReader_global.AddFieldOfInterest("v10");  //  1   | 10 metre V wind component
    FluxesReader_global.AddFieldOfInterest("ssrd"); //  2   | Surface solar radiation downwards
    FluxesReader_global.AddFieldOfInterest("strd"); //  3   | Surface thermal radiation downwards
    FluxesReader_global.AddFieldOfInterest("ro");   //  4   | Runoff
    FluxesReader_global.AddFieldOfInterest("tp");   //  5   | Total precipitation
    FluxesReader_global.AddFieldOfInterest("d2");   //  6   | Dewpoint temp at 2m 
    FluxesReader_global.AddFieldOfInterest("t2");   //  7   | Air temp at 2m 
    FluxesReader_global.AddFieldOfInterest("msl");  //  8   | Mean sea level pressure 
    FluxesReader_global.SetSimulationTimeUnits("seconds since 1970-01-02 00:00:0.0");
    FluxesReader_global.SetTimeSeconds(time);

    // set up the file containing the forcings from ERA40 - this is equivalent to our output
    // from the bulk formulae
    FluxesReader_ERAdata.RegisterDataFile("../../tests/data/stationPapa_1970_output.nc");
    FluxesReader_ERAdata.SetSimulationTimeUnits("seconds since 1970-01-02 00:00:0.0");
    FluxesReader_ERAdata.SetTimeSeconds(time);
    FluxesReader_ERAdata.AddFieldOfInterest("sshf"); // 0 | surface sensible heat flux
    FluxesReader_ERAdata.AddFieldOfInterest("slhf"); // 1 | surface latent heat flux
    FluxesReader_ERAdata.AddFieldOfInterest("ssr");  // 2 | surface solar radiation
    FluxesReader_ERAdata.AddFieldOfInterest("str");  // 3 | surface thermal radiation
    FluxesReader_ERAdata.AddFieldOfInterest("ewss"); // 4 | east-west surface stress 
    FluxesReader_ERAdata.AddFieldOfInterest("nsss"); // 5 | north-south surface stress 
    FluxesReader_ERAdata.AddFieldOfInterest("e");    // 6 | evaporation    

    int const nFields_in=9, nFields_out=7;
    double values_in[nFields_in];
    double values_out[nFields_out];
    int ret;
    double accumulated_correction = 6.0*60.0*60.0;

    int n = 1;
    // Ocean Station Papa (50, -145)
    // -3.35835e+06,-2.35154e+06,4.88594e+06
    X[0] = -3.35835e+06;
    Y[0] = -2.35154e+06;
    Z[0] = 4.88594e+06;
    y[0] = Y[0] , x[0] = X[0] , z[0] = Z[0];
    ret = projections(n, &x[0], &y[0], &z[0], "cart", "spherical");
    if (ret != 0) {
        cerr<<"Error converting coord system"<<endl;
    }
    
    FluxesReader_ERAdata.GetScalars(x[0], y[0], values_out);
    FluxesReader_global.GetScalars(x[0], y[0], values_in);
    double temp =  (values_in[3] - values_out[3]);
    temp = temp / accumulated_correction; // correcting for accumulated values
    temp = temp / 5.67e-8;                // SB constant
    temp = pow(temp, 1.0/4.0);
    T[0] = temp;
    double ocean_density = 1027.0;
    double ref_salinity = 35.0;

    double sumQ = 0, sumF = 0, sumTauX = 0, sumTauY = 0;
    double sumQ_era40 = 0, sumF_era40 = 0, sumTauX_era40 = 0, sumTauY_era40 = 0;
    double sumQ2 = 0, sumF2 = 0, sumTauX2 = 0, sumTauY2 = 0;
    double sumQ_era402 = 0, sumF_era402 = 0, sumTauX_era402 = 0, sumTauY_era402 = 0;
    int nPoints =  int(25920000-86400/21600);
    vector<double> Q_flux(nPoints,0), F_flux(nPoints,0), TauX_flux(nPoints,0), TauY_flux(nPoints,0);
    vector<double> Q_era40(nPoints,0), F_era40(nPoints,0), TauX_era40(nPoints,0), TauY_era40(nPoints,0);
    nPoints = 0;
    int formula;

#ifdef GRAPHIC
    ofstream test;
    test.open("test_data_kara.csv");
    test.clear();
    test <<"Time,ERA40_Q,Flux_Q,ERA40_Tau_X,Flux_Tau_X,ERA40_Tau_Y,Flux_Tau_Y,ERA40_F,Flux_F" << endl; 
#endif
    
    // Do the fluxes for several days then compare
    for (int t = 86400 ; t < 25920000; t += 21600) {
    
        FluxesReader_ERAdata.SetTimeSeconds(t);
        FluxesReader_ERAdata.GetScalars(x[0], y[0], values_out);
        FluxesReader_global.SetTimeSeconds(t);
        FluxesReader_global.GetScalars(x[0], y[0], values_in);
        // now let's work out the fluxes for our test location
        double Tau_u_out = values_out[4]/(accumulated_correction*ocean_density);
        double Tau_v_out = values_out[5]/(accumulated_correction*ocean_density);
        double F_out = values_in[5] + values_out[6] + values_in[4];
        F_out = -F_out / (accumulated_correction) * ref_salinity;
        double Q_p = -3.337e5 * values_in[5];
        double Q_out = values_out[2] + values_out[3] + values_out[0] + values_out[1] + Q_p;
        // 4000 = ocean heat capacity
        Q_out = Q_out / (accumulated_correction*ocean_density*4000.0);

        time = t;

        // get fluxes for test location
        formula = 2;
        get_era40_fluxes_fc(&time, &X[0], &Y[0], &Z[0], &T[0], &Vx[0], &Vy[0], &Vz[0], &S[0], &F_as[0],
                            &Q_as[0], &tau_u[0], &tau_v[0], &Q_solar[0], &n, false, &formula);

        
        sumQ2 += Q_as[0]*Q_as[0];
        sumF2 += F_as[0]*F_as[0];
        sumTauX2 += tau_u[0]*tau_u[0];
        sumTauY2 += tau_v[0]*tau_v[0];
        sumQ_era402 += Q_out*Q_out;
        sumF_era402 += F_out*F_out;
        sumTauX_era402 += Tau_u_out*Tau_u_out;
        sumTauY_era402 += Tau_v_out*Tau_v_out;        
        sumQ += Q_as[0];
        sumF += F_as[0];
        sumTauX += tau_u[0];
        sumTauY += tau_v[0];
        sumQ_era40 += Q_out;
        sumF_era40 += F_out;
        sumTauX_era40 += Tau_u_out;
        sumTauY_era40 += Tau_v_out;
        Q_flux[nPoints] = Q_as[0];
        F_flux[nPoints] = F_as[0];
        TauX_flux[nPoints] = tau_u[0];
        TauY_flux[nPoints] = tau_v[0];
        Q_era40[nPoints] = Q_out;
        F_era40[nPoints] = F_out;
        TauX_era40[nPoints] = Tau_u_out;
        TauY_era40[nPoints] = Tau_v_out;

/**        cout<<values_out[0]/accumulated_correction<<","<<q_h[0]<<endl;
        cout <<"Flux                 |    ERA40 value     |  Our value   "<<endl;
        cout <<"---------------------+--------------------+---------------"<<endl;
        printf("U stress (m2/s2)     | %-18g | %-18g\n",Tau_u_out,tau_u[0]);
        printf("V stress (m2/s2)     | %-18g | %-18g\n",Tau_v_out,tau_v[0]);
        printf("heat flux (K/s)      | %-18g | %-18g\n",Q_out,Q_as[0]);
        printf("Salinity flux        | %-18g | %-18g\n",F_out,F_as[0]);
        printf("Solar rad. (W/m2)    | %-18g | %-18g\n",values_out[2]/accumulated_correction,values_in[2]/accumulated_correction*(1-0.066));
        cout<<endl<<endl;
**/
        nPoints++;
#ifdef GRAPHIC
    test <<nPoints<< "," << Q_out << "," << Q_as[0]<< "," << Tau_u_out << "," <<tau_u[0]<<","<<Tau_v_out<<","<<tau_v[0]<<","<<F_out<<","<<F_as[0] << endl; 
#endif
    /*printf("%i,%-5g,%-5g,%-5g,%-5g,%-5g,%-5g,%-5g,%-5g\n", nPoints,
                                                 Q_out,Q_as[0],
                                                 Tau_u_out,tau_u[0],
                                                 Tau_v_out,tau_v[0],
                                                 F_out,F_as[0]);
*/
    }

#ifdef GRAPHIC
    test.close();
#endif

    sumQ /= nPoints;
    sumF /= nPoints;
    sumTauX /= nPoints;
    sumTauY /= nPoints;
    sumQ_era40 /= nPoints;
    sumF_era40 /= nPoints;
    sumTauX_era40 /= nPoints;
    sumTauY_era40 /= nPoints;

    // work out some correlation coefficient
    double xy = 0;
    for (int i=0; i<nPoints; i++) {
        xy += Q_flux[i] * Q_era40[i];
    }
    double r_q = (nPoints * xy - sumQ*sumQ_era40) /
                 (sqrt(nPoints*sumQ2-(sumQ*sumQ))*sqrt(nPoints*sumQ_era402-(sumQ_era40*sumQ_era40)));


    xy = 0;
    for (int i=0; i<nPoints; i++) {
        xy += F_flux[i] * F_era40[i];
    }
    double r_f = (nPoints * xy - sumF*sumF_era40) /
                 (sqrt(nPoints*sumF2-(sumF*sumF))*sqrt(nPoints*sumF_era402-(sumF_era40*sumF_era40)));

    xy = 0;
    for (int i=0; i<nPoints; i++) {
        xy += TauX_flux[i] * TauX_era40[i];
    }
    double r_tauX = (nPoints * xy - sumTauX*sumTauX_era40) /
                 (sqrt(nPoints*sumTauX2-(sumTauX*sumTauX))*sqrt(nPoints*sumTauX_era402-(sumTauX_era40*sumTauX_era40)));

    xy = 0;
    for (int i=0; i<nPoints; i++) {
        xy += TauY_flux[i] * TauY_era40[i];
    }
    double r_tauY = (nPoints * xy - sumTauY*sumTauY_era40) /
                 (sqrt(nPoints*sumTauY2-(sumTauY*sumTauY))*sqrt(nPoints*sumTauY_era402-(sumTauY_era40*sumTauY_era40)));


    // Check that we have roughly the right amount of forcing in total
    fail = true;
    warn = false;
    error = (sumQ_era40-sumQ)/sumQ;
    if (abs(error) < .75) fail = false;
    sprintf(errorMessage, "Error > .75 Is %g",error); 
    report_test("[test_kara_ocean_forcing: heat flux error to ERA40]",fail,warn,errorMessage);

    fail = true;
    error = (sumF-sumF_era40)/sumF;
    if (abs(error) < .75) fail = false;
    error = (sumF-sumF_era40)/sumF;
    sprintf(errorMessage, "Error > .75 Is %g",error);
    report_test("[test_kara_ocean_forcing: salinity flux error to ERA40]",fail,warn,errorMessage);

    fail = true;
    error = (sumTauX-sumTauX_era40)/sumTauX;
    if (abs(error) < .75) fail = false;
    sprintf(errorMessage, "Error > .75 Is %g",error);
    report_test("[test_kara_ocean_forcing: tau x momentum error to ERA40]",fail,warn,errorMessage);

    fail = true;
    error = (sumTauY-sumTauY_era40)/sumTauY;
    if (abs(error) < .75) fail = false;
    sprintf(errorMessage, "Error > .75 Is %g",error);
    report_test("[test_kara_ocean_forcing: tau y momentum error to ERA40]",fail,warn,errorMessage);

    // check we have the right direciton of forcing
    if (r_q > .88) fail = false;
    sprintf(errorMessage, "Correlation coefficient < .88: Is %g",r_q); 
    report_test("[test_kara_ocean_forcing: heat flux correlation]",fail,warn,errorMessage);

    fail = true;
    if (r_f > .88) fail = false;
    sprintf(errorMessage, "Correlation coefficient < .88: Is %g",r_f); 
    report_test("[test_kara_ocean_forcing: freshwater flux correlation]",fail,warn,errorMessage);

    fail = true;
    if (r_tauX > .88) fail = false;
    sprintf(errorMessage, "Correlation coefficient < .88: Is %g",r_tauX); 
    report_test("[test_kara_ocean_forcing: momentum (X) flux correlation]",fail,warn,errorMessage);

    fail = true;
    if (r_tauY > .88) fail = false;
    sprintf(errorMessage, "Correlation coefficient < .88: Is %g",r_tauY); 
    report_test("[test_kara_ocean_forcing: momentum (Y) flux correlation]",fail,warn,errorMessage);

#endif
}
