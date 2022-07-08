#!/bin/bash

#file_bmrw_data="../data_proton_bmrw.root"
#file_bmrw_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_proton_bmrw.root"

file_bmrw_data="../data_proton_beamxy_bmrw.root"
file_bmrw_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_proton_beamxy_bmrw.root"

#file_bmrw_data="../data_proton_beamxy_beammom_bmrw.root"
#file_bmrw_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_proton_beamxy_beammom_bmrw.root"


file_babmrw_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_proton_after_bmrw.root"
file_babmrw_data="../data_proton_bmrw.root"

file_tightxy_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_proton_tightxy_after_bmrw.root"
file_tightxy_data="../data_proton_tightxy_bmrw.root"

file_hd_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_proton_after_bmrw_HD.root"
file_hd_data="../data_proton_bmrw_HD.root"

#file_bmrw_data="../data_proton_bmrw2_usedefault_range_calc.root"
#file_bmrw_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_proton_bmrw2_usedefault_range_calc.root"


#fout_path="./plots_bmrw"
#fout_path="./plots_beamxy_bmrw"
fout_path="./plots_beamxy_bmrw_newwithdataGaussian"
#fout_path="./plots_beamxy_beammom_bmrw"

fout_bkg_path="./plots_bkgrich_bmrw"

exe0_str="root -b -q 'plotXY.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"")'"
exe1_str="root -b -q 'plotBeamDataMC.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"")'"
exe2_str="root -b -q 'plotBMRW.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"")'"
exe3_str1="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_BQ\", \"trklen_BQ_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_path/trklen_data_mc.eps\"")'"
exe3_str2="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_RecoInel\", \"trklen_RecoInel_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_path/trklen_data_mc_recoinel.eps\"")'"

exe3_str5="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_CaloSz\", \"trklen_CaloSz_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_path/trklen_data_mc_calosz.eps\"")'"
exe3_str6="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_Pos\", \"trklen_Pos_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_path/trklen_data_mc_pos.eps\"")'"


exe3_str7="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_BQ\", \"trklen_BQ_bmrw\",\"Proton Track Length [cm]\",\""$file_tightxy_data\"", \""$file_tightxy_mc\"", \""$fout_path/trklen_data_mc_tightxy.eps\"")'"
exe3_str8="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_RecoInel\", \"trklen_RecoInel_bmrw\",\"Proton Track Length [cm]\",\""$file_tightxy_data\"", \""$file_tightxy_mc\"", \""$fout_path/trklen_data_mc_recoinel_tightxy.eps\"")'"


exe3_str9="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_BQ\", \"trklen_BQ_bmrw\",\"Proton Track Length [cm]\",\""$file_hd_data\"", \""$file_hd_mc\"", \""$fout_path/trklen_data_mc_HD_zoom.eps\"")'"
exe3_str10="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_RecoInel\", \"trklen_RecoInel_bmrw\",\"Proton Track Length [cm]\",\""$file_hd_data\"", \""$file_hd_mc\"", \""$fout_path/trklen_data_mc_recoinel_HD_zoom.eps\"")'"
exe3_str11="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_CaloSz\", \"trklen_CaloSz_bmrw\",\"Proton Track Length [cm]\",\""$file_hd_data\"", \""$file_hd_mc\"", \""$fout_path/trklen_data_mc_calosz_HD_zoom.eps\"")'"
exe3_str12="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_Pos\", \"trklen_Pos_bmrw\",\"Proton Track Length [cm]\",\""$file_hd_data\"", \""$file_hd_mc\"", \""$fout_path/trklen_data_mc_pos_HD_zoom.eps\"")'"




exe3_str13="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_MidpRich\", \"trklen_MidpRich_bmrw\",\"Proton Track Length [cm]\",\""$file_hd_data\"", \""$file_hd_mc\"", \""$fout_bkg_path/trklen_data_mc_midprich_HD.eps\"")'"
exe3_str14="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_ElRich\", \"trklen_ElRich_bmrw\",\"Proton Track Length [cm]\",\""$file_hd_data\"", \""$file_hd_mc\"", \""$fout_bkg_path/trklen_data_mc_elrich_HD.eps\"")'"
exe3_str15="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_ElRich\", \"trklen_ElRich_bmrw\",\"Proton Track Length [cm]\",\""$file_hd_data\"", \""$file_hd_mc\"", \""$fout_bkg_path/trklen_data_mc_elrich_HD_zoom.eps\"")'"
exe3_str16="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_CaloSz_ElRich\", \"trklen_CaloSz_ElRich_bmrw\",\"Proton Track Length [cm]\",\""$file_hd_data\"", \""$file_hd_mc\"", \""$fout_bkg_path/trklen_CaloSz_ElRich_HD.eps\"")'"
exe3_str17="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_CaloSz_ElRich\", \"trklen_CaloSz_ElRich_bmrw\",\"Proton Track Length [cm]\",\""$file_hd_data\"", \""$file_hd_mc\"", \""$fout_bkg_path/trklen_CaloSz_ElRich_HD_zoom.eps\"")'"
exe3_str18="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_Pos_ElRich\", \"trklen_Pos_ElRich_bmrw\",\"Proton Track Length [cm]\",\""$file_hd_data\"", \""$file_hd_mc\"", \""$fout_bkg_path/trklen_Pos_ElRich_HD_zoom.eps\"")'"
exe3_str19="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_Pos_ElRich\", \"trklen_Pos_ElRich_bmrw\",\"Proton Track Length [cm]\",\""$file_hd_data\"", \""$file_hd_mc\"", \""$fout_bkg_path/trklen_Pos_ElRich_HD.eps\"")'"


exe3_str20="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_MidpRich\", \"trklen_MidpRich_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_bkg_path/trklen_data_mc_midprich.eps\"")'"
exe3_str21="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_ElRich\", \"trklen_ElRich_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_bkg_path/trklen_data_mc_elrich.eps\"")'"


exe4_str1="root -b -q 'plotMisIDPBkgFitAfterBMRW.C(\"trklen_MidpRich\", \"trklen_MidpRich_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_bkg_path/bkgfit_trklen_data_mc_midprich.eps\"")'"
exe4_str2="root -b -q 'plotElBkgFitAfterBMRW.C(\"trklen_ElRich\", \"trklen_ElRich_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_bkg_path/bkgfit_trklen_data_mc_elrich.eps\"")'"
exe4_str3="root -b -q 'plotElBkgFitAfterBMRW_AfterMisIDPCorr.C(\"trklen_ElRich\", \"trklen_ElRich_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_bkg_path/bkgfit_trklen_data_mc_elrich_corr.eps\"")'"
exe4_str4="root -b -q 'plotSGBKGFitAfterBMRW.C(\"trklen_RecoInel\", \"trklen_RecoInel_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_bkg_path/bkgfit_trklen_data_mc_inelrich.eps\"")'"
exe4_str5="root -b -q 'plotSGBKGFitAfterBMRW_AfterMisIDPCorr.C(\"trklen_RecoInel\", \"trklen_RecoInel_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_bkg_path/bkgfit_trklen_data_mc_inelrich_corr.eps\"")'"

exe4_str6="root -b -q 'plotSGBKGFitAfterBMRW_AfterElMisIDPCorr.C(\"trklen_RecoInel\", \"trklen_RecoInel_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_bkg_path/bkgfit_trklen_data_mc_inelrich_corr2.eps\"")'"

exe4_str7="root -b -q 'plotInelMisIDPFitAfterBMRW.C(\"trklen_MidpRich\", \"trklen_MidpRich_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_bkg_path/inelmisidpfit_trklen_data_mc_misidprich.eps\"")'"

exe5_str10="root -b -q 'plotBeforeAfterBMRWAfterEffCorr.C(\"trklen_RecoInel\", \"trklen_RecoInel_bmrw\",\"Proton Track Length [cm]\",\""$file_babmrw_data\"", \""$file_babmrw_mc\"", \""$fout_path/effcorr_recoinel.eps\"")'"

#exe3_str2="root -b -q 'plotBeforeAfterBMRW.C(\"trklen_XY\", \"trklen_bmrw_XY\",\"Proton Track Length [cm]\",\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path/trklen_XY_data_mc.eps\"")'"
#exe3_str3="root -b -q 'plotBeforeAfterBMRW.C(\"zend\", \"zend_bmrw\",\"Proton EndZ [cm]\",\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path/endz_data_mc.eps\"")'"
#exe3_str4="root -b -q 'plotBeforeAfterBMRW.C(\"zend_XY\", \"zend_bmrw_XY\",\"Proton EndZ [cm]\",\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path/endz_XY_data_mc.eps\"")'"


#echo $exe0_str" ......"
#eval $exe0_str

#echo $exe1_str" ......"
#eval $exe1_str

echo $exe2_str" ......"
eval $exe2_str

#echo $exe3_str1" ......"
#eval $exe3_str1

#echo $exe3_str5" ......"
#eval $exe3_str2
#eval $exe3_str1
#eval $exe3_str5
#eval $exe3_str6
#eval $exe3_str7
#eval $exe3_str8

#eval $exe3_str9
#eval $exe3_str10
#eval $exe3_str11
#eval $exe3_str12

#eval $exe3_str13
#eval $exe3_str14
#eval $exe3_str15
#eval $exe3_str16
#eval $exe3_str17
#eval $exe3_str18
#eval $exe3_str19

#eval $exe3_str20
#eval $exe3_str21

#eval $exe4_str1
#eval $exe4_str2
#eval $exe4_str3
#eval $exe4_str4
#eval $exe4_str5
#eval $exe4_str6
#eval $exe4_str7



#eval $exe5_str10




#echo $exe3_str2" ......"
#eval $exe3_str2

#echo $exe3_str3" ......"
#eval $exe3_str3

#echo $exe3_str4" ......"
#eval $exe3_str4

