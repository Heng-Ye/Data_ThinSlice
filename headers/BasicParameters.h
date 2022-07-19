//#include <string.h>

const double m_proton=0.938272046; //GeV/c2
const int binsize=5; // cm

double dedx_min=30.;

//x-y cut ---------------------------//
double mean_x=-26.59; //prod4a data
double mean_y=423.5; //prod4a data
double dev_x=1.5*3.744; //prod4a data
double dev_y=1.5*4.364; //prod4a data

double mean_x_mc=-29.25; //prod4a mc
double mean_y_mc=422.1; //prod4a mc
double dev_x_mc=1.5*4.4; //prod4a mc
double dev_y_mc=1.5*3.883; //prod4a mc
//x-y cut ---------------------------//

//z0 offset (before SCE corr. ------------------------------------------//
double z0_before_sce=32.4672;
double sigma0_before_sce=0.872971;

//beam quality cut -----------------------------------------------------------------------------------------------------------------//

//for run5387 only
//double mean_StartZ=3.75897e+00; //prod4a-reco2[run5387]
//double sigma_StartZ=1.03771e+00; //prod4a-reco2[run5387]
//double mean_StartY=4.23949e+02; //prod4a-reco2[run5387]
//double sigma_StartY=4.61753e+00; //prod4a-reco2[run5387]
//double mean_StartX=-2.83979e+01; //prod4a-reco2[run5387]
//double sigma_StartX=3.88357e+00; //prod4a-reco2[run5387]

//Global
double mean_StartZ=3.14342; //prod4a-reco2[global]
double sigma_StartZ=1.25748; //prod4a-reco2[global]
double mean_StartY=424.221; //prod4a-reco2[global]
double sigma_StartY=4.63213; //prod4a-reco2[global]
double mean_StartX=-28.7395; //prod4a-reco2[global]
double sigma_StartX=3.89379; //prod4a-reco2[global]

double min1_z=mean_StartZ-3.*sigma_StartZ; //prod4a-reco2
double min2_z=mean_StartZ+3.*sigma_StartZ; //prod4a-reco2
double min1_y=mean_StartY-3.*sigma_StartY; //prod4a-reco2
double min2_y=mean_StartY+3.*sigma_StartY; //prod4a-reco2
double min1_x=mean_StartX-3.*sigma_StartX; //prod4a-reco2
double min2_x=mean_StartX+3.*sigma_StartX; //prod4a-reco2

double dx_min=-3.; double dx_max=3.;
double dy_min=-3.; double dy_max=3.;
double dz_min=-3.; double dz_max=3.;
double dxy_min=-1.; double dxy_max=3.;
double costh_min = 0.96; double costh_max = 2;

//double cosine_beam_primtrk_min=(9.92194e-01)-4.*(3.96921e-03); //new p4 [spec]

//stopping proton cut
double mean_norm_trklen_csda=8.77176e-01; //prod4a-reco2
double sigma_norm_trklen_csda=6.98794e-02; //prod4a-reco2
double min_norm_trklen_csda=mean_norm_trklen_csda-2.*sigma_norm_trklen_csda;
double max_norm_trklen_csda=mean_norm_trklen_csda+3.*sigma_norm_trklen_csda;

//new inel cut using chi^2 info
double pid_1=7.5;
double pid_2=10.;

//mu of beam ---------------------------------------------------------------//
/*
   double mu_data=4.45451e+02;
   double err_mu_data=7.21974e-01;
   double sigma_data=5.09380e+01;
   double err_sigma_data=5.22759e-01;

   double mu_stop_data=4.05686e+02; //range-based calc
   double err_mu_stop_data=1.23011e+00;
   double sigma_stop_data=5.25443e+01;
   double err_sigma_stop_data=1.20796e+00;

   double mu_calo_stop_data=3.84030e+02; //calo-based calc
   double err_mu_calo_stop_data=1.15843e+00;
   double sigma_calo_stop_data=5.65971e+01;
   double err_sigma_calo_stop_data=9.02569e-01;

//double de_upstream_data=mu_data-mu_stop_data;
//double err_de_upstream_data=sqrt(pow(err_mu_data,2)+pow(err_mu_stop_data,2));
//double de_upstream_data=mu_data-mu_calo_stop_data;
double de_upstream_data=mu_data-mu_stop_data;
//double err_de_upstream_data=sqrt(pow(err_mu_data,2)+pow(err_mu_calo_stop_data,2));
*/
//--------------------------------------------------------------------------//


//beam xy cut ---------------//
//double meanX_data=-31.3139; //run5387
//double rmsX_data=3.79366; //run5387
//double meanY_data=422.116; //run5387
//double rmsY_data=3.48005; //run5387

double meanX_data=-31.1041; //global
double rmsX_data=3.90508; //global
double meanY_data=422.142; //global
double rmsY_data=3.55531; //global


double meanX_mc=-29.1637;
double rmsX_mc=4.50311;
double meanY_mc=421.76;
double rmsY_mc=3.83908;

