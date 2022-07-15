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
const int n_RUN=26;
int RUN[n_RUN];
double mu_Z[n_RUN]; double sigma_Z[n_RUN];
double mu_Y[n_RUN]; double sigma_Y[n_RUN];
double mu_X[n_RUN]; double sigma_X[n_RUN];
RUN[0]=5219; mu_Z[0]=2.39249; sigma_Z[0]=0.950201; mu_Y[0]=423.994; sigma_Y[0]=4.61662; mu_X[0]=-28.5979; sigma_X[0]=3.60422;
RUN[1]=5225; mu_Z[1]=2.81494; sigma_Z[1]=0.769699; mu_Y[1]=423.824; sigma_Y[1]=4.37359; mu_X[1]=-28.165; sigma_X[1]=3.89314;
RUN[2]=5235; mu_Z[2]=2.00167; sigma_Z[2]=0.934767; mu_Y[2]=424.435; sigma_Y[2]=4.5774; mu_X[2]=-28.898; sigma_X[2]=3.87618;
RUN[3]=5240; mu_Z[3]=2.02353; sigma_Z[3]=0.730619; mu_Y[3]=423.644; sigma_Y[3]=4.34303; mu_X[3]=-29.0079; sigma_X[3]=3.96048;
RUN[4]=5244; mu_Z[4]=1.99276; sigma_Z[4]=0.904941; mu_Y[4]=423.822; sigma_Y[4]=4.69817; mu_X[4]=-29.0475; sigma_X[4]=3.75481;
RUN[5]=5308; mu_Z[5]=3.79232; sigma_Z[5]=0.871755; mu_Y[5]=424.259; sigma_Y[5]=4.6139; mu_X[5]=-28.6683; sigma_X[5]=3.87115;
RUN[6]=5311; mu_Z[6]=3.48456; sigma_Z[6]=1.16812; mu_Y[6]=424.047; sigma_Y[6]=4.59582; mu_X[6]=-28.6148; sigma_X[6]=3.866;
RUN[7]=5315; mu_Z[7]=4.30158; sigma_Z[7]=1.11292; mu_Y[7]=424.212; sigma_Y[7]=4.6911; mu_X[7]=-28.1816; sigma_X[7]=3.91249;
RUN[8]=5338; mu_Z[8]=4.11004; sigma_Z[8]=0.925806; mu_Y[8]=423.877; sigma_Y[8]=4.52686; mu_X[8]=-29.0855; sigma_X[8]=3.91689;
RUN[9]=5387; mu_Z[9]=3.75067; sigma_Z[9]=1.04507; mu_Y[9]=423.962; sigma_Y[9]=4.58983; mu_X[9]=-28.3621; sigma_X[9]=3.85745;
RUN[10]=5423; mu_Z[10]=3.92263; sigma_Z[10]=1.40972; mu_Y[10]=423.742; sigma_Y[10]=4.41215; mu_X[10]=-28.6561; sigma_X[10]=3.59162;
RUN[11]=5424; mu_Z[11]=4.16598; sigma_Z[11]=0.977282; mu_Y[11]=424.263; sigma_Y[11]=4.57237; mu_X[11]=-28.292; sigma_X[11]=3.8009;
RUN[12]=5426; mu_Z[12]=3.86441; sigma_Z[12]=0.995988; mu_Y[12]=424.462; sigma_Y[12]=4.52192; mu_X[12]=-28.3404; sigma_X[12]=3.7884;
RUN[13]=5455; mu_Z[13]=3.24509; sigma_Z[13]=0.977819; mu_Y[13]=424.167; sigma_Y[13]=4.75984; mu_X[13]=-28.9719; sigma_X[13]=3.95321;
RUN[14]=5456; mu_Z[14]=3.09068; sigma_Z[14]=0.826214; mu_Y[14]=424.425; sigma_Y[14]=4.60255; mu_X[14]=-29.692; sigma_X[14]=3.89421;
RUN[15]=5457; mu_Z[15]=3.81861; sigma_Z[15]=0.683909; mu_Y[15]=423.763; sigma_Y[15]=5.25214; mu_X[15]=-28.4411; sigma_X[15]=3.75078;
RUN[16]=5458; mu_Z[16]=3.11472; sigma_Z[16]=0.971477; mu_Y[16]=424.415; sigma_Y[16]=4.64395; mu_X[16]=-29.1837; sigma_X[16]=3.70776;
RUN[17]=5460; mu_Z[17]=3.31202; sigma_Z[17]=0.991936; mu_Y[17]=424.534; sigma_Y[17]=4.68872; mu_X[17]=-28.9989; sigma_X[17]=3.77217;
RUN[18]=5809; mu_Z[18]=2.38939; sigma_Z[18]=1.03706; mu_Y[18]=423.564; sigma_Y[18]=4.35458; mu_X[18]=-28.8828; sigma_X[18]=3.7523;
RUN[19]=5810; mu_Z[19]=1.85908; sigma_Z[19]=0.890347; mu_Y[19]=423.631; sigma_Y[19]=4.58065; mu_X[19]=-28.583; sigma_X[19]=3.80641;
RUN[20]=5814; mu_Z[20]=2.33636; sigma_Z[20]=0.926856; mu_Y[20]=423.926; sigma_Y[20]=4.47281; mu_X[20]=-28.7546; sigma_X[20]=3.92464;
RUN[21]=5816; mu_Z[21]=2.87698; sigma_Z[21]=1.07482; mu_Y[21]=424.269; sigma_Y[21]=4.49433; mu_X[21]=-28.9165; sigma_X[21]=3.81153;
RUN[22]=5817; mu_Z[22]=3.25487; sigma_Z[22]=0.957823; mu_Y[22]=424.434; sigma_Y[22]=4.59399; mu_X[22]=-28.9055; sigma_X[22]=3.89867;
RUN[23]=5842; mu_Z[23]=2.63629; sigma_Z[23]=0.965875; mu_Y[23]=424.508; sigma_Y[23]=4.70966; mu_X[23]=-28.9506; sigma_X[23]=3.92468;
RUN[24]=5843; mu_Z[24]=1.82563; sigma_Z[24]=1.03683; mu_Y[24]=424.612; sigma_Y[24]=4.62996; mu_X[24]=-29.3553; sigma_X[24]=3.60679;
RUN[25]=5844; mu_Z[25]=1.79371; sigma_Z[25]=0.940908; mu_Y[25]=424.479; sigma_Y[25]=4.83961; mu_X[25]=-29.0593; sigma_X[25]=3.88118;

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


//beam xy cut
double meanX_data=-31.3139;
double rmsX_data=3.79366;
double meanY_data=422.116;
double rmsY_data=3.48005;

double meanX_mc=-29.1637;
double rmsX_mc=4.50311;
double meanY_mc=421.76;
double rmsY_mc=3.83908;

