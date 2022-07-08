#include "TMinuit.h"
#include "TH1D.h"
#include <iostream>

class TH1D;

class TemplateSimuFitter {
  
 public:
  
  TemplateSimuFitter();
  
  void SetHistograms(TH1D *hist0, TH1D *hist1, TH1D *hist2, TH1D *hist3);
  static void fcn(int &npar, double *gin, double &f, double *par, int iflag);
  void SetFitRange(int imin, int imax);

  void Fit();
  
  vector<double> GetPar();
  vector<double> GetParError();

  static TH1D *h0, *h1, *h2, *h3;
  static int i0, i1;

  TMinuit *gMinuit;

 private:

  bool fitsuccess;
  
};

TH1D* TemplateSimuFitter::h0;
TH1D* TemplateSimuFitter::h1;
TH1D* TemplateSimuFitter::h2;
TH1D* TemplateSimuFitter::h3;

int TemplateSimuFitter::i0;
int TemplateSimuFitter::i1;


TemplateSimuFitter::TemplateSimuFitter(){
  gMinuit = new TMinuit(2);
  gMinuit->SetFCN(fcn);
  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
}

void TemplateSimuFitter::SetHistograms(TH1D *hist0, TH1D *hist1, TH1D *hist2, TH1D *hist3){
  h0 = hist0;
  h1 = hist1;
  h2 = hist2;
  h3 = hist3;

  i0 = 1;
  i1 = h0->GetNbinsX();

}

void TemplateSimuFitter::SetFitRange(int imin, int imax){
  if (imin>=0) i0 = imin;
  if (imax>=0) i1 = imax;
}

void TemplateSimuFitter::fcn(int &npar, double *gin, double &f, double *par, int iflag){
  
  double chisq = 0;

  for (int i = i0; i<=i1; ++i){
    double x = h0->GetBinContent(i);
    double y = h1->GetBinContent(i);
    double z = h2->GetBinContent(i);
    double t = h3->GetBinContent(i);

    
    double ex = h0->GetBinError(i);
    double ey = h1->GetBinError(i);
    double ez = h2->GetBinError(i);
    double et = h3->GetBinError(i);


    if (!ex) ex = 1;
    chisq += pow(x-y-z*par[0]-t*par[1],2)/(pow(ex,2)+pow(ey,2)+pow(ez*par[0],2)+pow(et*par[1],2));
  }

  f = chisq;

}

void TemplateSimuFitter::Fit(){

  double arglist[10];
  int ierflg = 0;

  gMinuit->mncler();
  double vstart = 1;
  double step = 0.01;
  gMinuit->mnparm(0,"corr_fact0",vstart,step,0,10,ierflg);
  gMinuit->mnparm(1,"corr_fact1",vstart,step,0,10,ierflg);

  fitsuccess = false;

  if (h0->Integral(i0,i1) && h2->Integral(i0,i1) && h3->Integral(i0,i1)){
    arglist[0] = 500;
    arglist[1] = 1;
    gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
    
    double par, epar;
    double par2, epar2;
    gMinuit->GetParameter(0,par,epar);
    gMinuit->GetParameter(1,par2,epar2);

    TString test =  gMinuit->fCstatu.Data(); 
    if (test.EqualTo("CONVERGED ")){
      std::cout<<"Best fit1 = "<<par<<" error = "<<epar<<std::endl;
      std::cout<<"Best fit2 = "<<par2<<" error = "<<epar2<<std::endl;
      fitsuccess = true;
    }
  }
  else{
    std::cout<<"No fit was done because data and/or template are empty."<<std::endl;
  }
}

vector<double> TemplateSimuFitter::GetPar(){
  if (fitsuccess){
    double par, epar;
    double par2, epar2;
    gMinuit->GetParameter(0,par,epar);
    gMinuit->GetParameter(1,par2,epar2);
    vector<double> PAR;
    PAR.push_back(par);
    PAR.push_back(par2);
    return PAR;
  }
  else{
    vector<double> FPAR;
    FPAR.push_back(1);
    FPAR.push_back(1);
    return FPAR;
  }
}

vector<double> TemplateSimuFitter::GetParError(){
  if (fitsuccess){
    double par, epar;
    double par2, epar2;
    gMinuit->GetParameter(0,par,epar);
    gMinuit->GetParameter(1,par2,epar2);
    vector<double> EPAR;
    EPAR.push_back(epar);
    EPAR.push_back(epar2);
    return EPAR;
  }
  else{
    vector<double> FEPAR;
    FEPAR.push_back(2);
    FEPAR.push_back(2);

    return FEPAR;
  }
}


