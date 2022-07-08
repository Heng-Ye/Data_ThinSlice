#!/bin/bash

#file_bmrw_data="../data_kebkg.root"
#file_bmrw_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_kecalobkg_bmrw_new.root"

file_bmrw_data="../data_kebkg_beamxy.root"
file_bmrw_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_kecalobkg_bmrw_beamxy_new.root"

#fout_path="./plots_bmrw"
#fout_path="./plots_beamxy_bmrw"
#fout_path="./plots_beamxy_beammom_bmrw"

fout_path="./plotKE_BKG_beamXY_new/"
#fout_path="./plotKE_BKG_new/"

obs="keff_reco_RecoEl"
obs_end="ke_reco_RecoEl"

fout_file=${fout_path}${obs}"_bmrw.eps"

xmin=220
xmax=570
ymax=600 

exe0_str="root -b -q 'plot_protonKE_BKG_Data_MC.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"", \""$obs\"", $xmin, $xmax, $ymax)'"
exe1_str="root -b -q 'plot_protonKE_BKG_Data_MC.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"", \""$obs_end\"", -200, 550, $ymax)'"
exe2_str="root -b -q 'plot_protonKE_BKG_Data_MC.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"", \""$obs_end\"", -200, 550, 1200)'"

#echo $exe0_str" ......"
eval $exe0_str
eval $exe1_str
eval $exe2_str


