import ROOT as RT
import sys
import math

from argparse import ArgumentParser as ap
import numpy as np
from array import array
from math import sqrt
from math import exp


def fitg(x, p):
  m=p[0]
  s=p[1]
  n=p[2]
  g=n*math.exp(-(x[0]-m)*(x[0]-m)/(2*s*s))
  return g


def VNFit(h, pre_mean, n_sigma):
  pre_max=h.GetMaximum()
  pre_rms=h.GetRMS()
  print('pre_max:',pre_max)
  print('pre_rms:',pre_rms)
  print('pre_max:',pre_max)
  
  #pre-fit
  gg=RT.TF1("gg", fitg, pre_mean-n_sigma*pre_rms, pre_mean+n_sigma*pre_rms, 3)
  gg.SetParameter(0, pre_mean)
  gg.SetParameter(1, pre_rms)
  gg.SetParameter(2, pre_max)
  h.Fit("gg","remn")

  #pos-fit
  g=RT.TF1("g", fitg, gg.GetParameter(0)-n_sigma*gg.GetParameter(1), gg.GetParameter(0)+n_sigma*gg.GetParameter(1), 3)
  g.SetParameter(0, gg.GetParameter(0))
  g.SetParameter(1, gg.GetParameter(1))
  g.SetParameter(2, gg.GetParameter(2))

  g.SetLineStyle(2)
  g.SetLineWidth(2)

  h.Fit("g","remn")
  return g

#Data File ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
file_data='/dune/data2/users/hyliao/protonana/v09_39_01/KEHY_ff/data_kehy_runAll.root'
file_data_after_corr='/dune/data2/users/hyliao/protonana/v09_39_01/KEHY_KEBEAMFF/proton_beamxy_beammom_runAll.root'
out_path='./keff_study'


#Plt Style -------------------------------------------------------------------------------
parser = ap()

RT.gROOT.LoadMacro("~/protoDUNEStyle.C")
RT.gROOT.SetBatch(); #PyROOT does NOT display any graphics in batch mode
RT.gStyle.SetOptStat(00000)
#RT.gStyle.SetErrorX(1.e-4)
RT.gStyle.SetTitleAlign(23)
RT.gStyle.SetTitleX(.5)
RT.gStyle.SetLineWidth(1)
tt = RT.TLatex();
tt.SetNDC();

#read files -------------------------------------------------------------
f_data=RT.TFile(file_data, "OPEN")

h2d_KEffbeam_KEhy_stop=f_data.Get("h2d_KEffbeam_KEhy_stop")
h2d_KEffbeam_KEhy_inel=f_data.Get("h2d_KEffbeam_KEhy_inel")
h2d_KEffbeam_KEhy_el=f_data.Get("h2d_KEffbeam_KEhy_el")
h1d_ratio_KEffbeam_KEhy_stop=f_data.Get("h1d_ratio_KEffbeam_KEhy_stop")
h1d_ratio_KEffbeam_KEhy_inel=f_data.Get("h1d_ratio_KEffbeam_KEhy_inel")
h1d_ratio_KEffbeam_KEhy_el=f_data.Get("h1d_ratio_KEffbeam_KEhy_el")

f_data_after_corr=RT.TFile(file_data_after_corr, "OPEN")
h1d_keffbeam_stop=f_data_after_corr.Get("h1d_keffbeam_stop")
h1d_kehy_stop=f_data_after_corr.Get("h1d_kehy_stop")


#[1] --------------------------------------------------------------------------
c0_ff_data=RT.TCanvas("c0_ff_data","",1200,900)
c0_ff_data.Divide(1,1)
c0_ff_data.cd(1)

f2d_data=RT.TH2D("f2d_data","", 400,200,600, 400,200,600)
f2d_data.SetTitle("Stopping Protons; KE_{beam}-#DeltaE [MeV]; KE(fit) [MeV]")
f2d_data.Draw("")

h2d_KEffbeam_KEhy_stop.Draw("colz same")
ll=RT.TLine(200,200,600,600)
ll.SetLineColor(2)
ll.SetLineStyle(2)
ll.Draw()
c0_ff_data.Print(out_path+'/keffbeam_kehy_stop.eps')

#[2]Ratio of (KEHY/(KEbeam-dE)) -----------------------------------------------------------------
c0_r_data=RT.TCanvas("c0_r_data","",1200,900)
c0_r_data.Divide(1,1)
c0_r_data.cd(1)
f2d_ratio_KEffbeam_KEhy_stop=RT.TH2D("f2d_ratio_KEffbeam_KEhy_stop","", 10,0.5,1.5, 100,0,3000)
f2d_ratio_KEffbeam_KEhy_stop.SetTitle("Stopping Protons; KE(fit)/(KE_{beam}-#DeltaE); ")
f2d_ratio_KEffbeam_KEhy_stop.Draw()
h1d_ratio_KEffbeam_KEhy_stop.Draw("same")

fit_ratio_KEffbeam_KEhy_stop=VNFit(h1d_ratio_KEffbeam_KEhy_stop, 1, 3)
fit_ratio_KEffbeam_KEhy_stop.SetLineColor(2)
fit_ratio_KEffbeam_KEhy_stop.SetMarkerStyle(2)
fit_ratio_KEffbeam_KEhy_stop.Draw("same")

mu=fit_ratio_KEffbeam_KEhy_stop.GetParameter(0)
er_mu=fit_ratio_KEffbeam_KEhy_stop.GetParError(0)
sigma=fit_ratio_KEffbeam_KEhy_stop.GetParameter(1)
er_sigma=fit_ratio_KEffbeam_KEhy_stop.GetParError(1)

print('mu:',mu)
print('sigma:',sigma)

leg0=RT.TLegend(0.1,0.7,.92,0.95)
leg0.SetFillStyle(0)
txt0=[]
txt0.append("Data: #mu={:.4f}#pm{:.4f} MeV, #sigma={:.4f}#pm{:.4f} MeV".format(mu,er_mu,sigma,er_sigma))
leg0.AddEntry(h1d_ratio_KEffbeam_KEhy_stop, txt0[0], "l")
#leg0.SetNColumns(2);
leg0.Draw()

c0_r_data.Print(out_path+'/ratio_kehy_keffbeam_stop.eps')


#[3] --------------------------------------------------------------------------
#After corr
c1_ff_data=RT.TCanvas("c1_ff_data","",1200,900)
c1_ff_data.Divide(1,1)
c1_ff_data.cd(1)

f2dx_data=RT.TH2D("f2d_data","", 400, 200, 600, 600, 0, 600)
f2dx_data.SetTitle("Stopping Protons; Proton KE at TPC FF [MeV]; ")
f2dx_data.Draw("")
h1d_kehy_stop.SetLineColor(4)
h1d_kehy_stop.SetMarkerColor(4)
h1d_keffbeam_stop.SetLineColor(1)
h1d_keffbeam_stop.SetMarkerColor(1)

h1d_keffbeam_stop.Draw("ep same")
h1d_kehy_stop.Draw("ep same")

#fit Gaussians
fit_h1d_keffbeam_stop=VNFit(h1d_keffbeam_stop, 400, 3)
fit_h1d_keffbeam_stop.SetLineColor(1)
fit_h1d_keffbeam_stop.SetMarkerStyle(2)

fit_h1d_kehy_stop=VNFit(h1d_kehy_stop, 400, 3)
fit_h1d_kehy_stop.SetLineColor(4)
fit_h1d_kehy_stop.SetMarkerStyle(2)

fit_h1d_keffbeam_stop.Draw("same")
fit_h1d_kehy_stop.Draw("same")

n_data=2

mu=[]
er_mu=[]
sigma=[]
er_sigma=[]

for i in range(n_data):
  fit=fit_h1d_keffbeam_stop
  if i==0: 
    fit=fit_h1d_keffbeam_stop
  if i==1: 
    fit=fit_h1d_kehy_stop
  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)
  mu.append(m)
  er_mu.append(er_m)
  sigma.append(s)
  er_sigma.append(er_s)
  print("i=",i," m=",m, "s=",s,"\n")

leg0=RT.TLegend(0.1,0.7,.92,0.85)
leg0.SetFillStyle(0)
txt0=[]
txt0.append("(KE(beam)-#DeltaE)*R: #mu={:.4f}#pm{:.4f} MeV, #sigma={:.4f}#pm{:.4f} MeV".format(mu[0],er_mu[0],sigma[0],er_sigma[0]))
txt0.append("KE(Fit): #mu={:.4f}#pm{:.4f} MeV, #sigma={:.4f}#pm{:.4f} MeV".format(mu[1],er_mu[1],sigma[1],er_sigma[1]))
leg0.AddEntry(h1d_keffbeam_stop, txt0[0], "l")
leg0.AddEntry(h1d_kehy_stop, txt0[1], "l")
leg0.Draw()


c1_ff_data.Print(out_path+'/keffbeamCorr_kehy_stop.eps')


