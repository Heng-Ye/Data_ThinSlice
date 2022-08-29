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
c0_ff_mc.Print(out_path+'/keffbeam_kehy_stop.eps')


