#include "math.h"

Double_t tanh_fit(Double_t *x, Double_t *par) {
 float result;
 result=0.5+0.5*(tanh((x[0]-par[0])/par[1]));
 return result;
}



Double_t eff_fit(Double_t *x, Double_t *par) {
 double result;
 double amp=par[0];
 double A=TMath::Exp(-(x[0]-par[1])/2);
 double B=par[2];
 result=amp*pow(1.-A, B);
 return result;
}

