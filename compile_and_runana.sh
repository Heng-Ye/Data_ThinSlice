#!/bin/bash

#run=5225
#run=5460
run=5387
#run=5303
#run=5308

#class_name="ProtonACAmap_run"
#class_name="ProtonReweight_AjibACA_run"
#class_name="ProtonReweight_run"

#class_name="ProtonThinSlice_run"

class_name="ProtonMomentumReweight_run"
#class_name="ProtonBackgroundFit_run"
#class_name="ProtonEvtDisplay_run"
#class_name="ProtonMassProduction_run"
#class_name="ProtonRecombination_run"
#class_name="ProtonKE_run"

#class_name="ProtonESlice_run"

#class_name="ProtonBetheBlochKE_run"
#class_name="ProtonCaloKE_run"
#class_name="ProtonDataDrivenBKGMeas_BetheBloch_run"

#class_name="ProtonDataDrivenBKGMeas_run"

#class_name="ProtonSelector_run"
class_namex=$class_name"X"

selector_name=$class_name$run
ana_name="makeproton_ana_"$class_name$run

header_name=$class_name$run".h"
tmp_header_name=$class_name${run}"_tmp.h"

class_code=$class_name${run}".C"
tmp_class_code=$class_name$run"_tmp.C"



echo $class_namex

#[1]Generate the file list for analysis
#root -b -q 'HealthCheck_Files.C("/pnfs/dune/scratch/users/hyliao/v07_11_00/ana/protodune_1gev_180kv_run'${run}'")'
#root -b -q 'HealthCheck_Files.C("/pnfs/dune/scratch/users/hyliao/v07_11_00/ana/protodune_1gev_180kv_run'${run}'_calo_nohitxyz")'
#root -b -q 'HealthCheck_Files.C("/pnfs/dune/scratch/users/hyliao/v07_11_00/ana/protodune_1gev_180kv_run'${run}'_calo")'
#root -b -q 'HealthCheck_Files.C("/pnfs/dune/scratch/users/hyliao/v08_07_00/ana/protodune_1gev_180kv_run'${run}'_calo_conv")'
#root -b -q 'HealthCheck_Files.C("/pnfs/dune/scratch/users/hyliao/v08_07_00/ana/protodune_1gev_180kv_run'${run}'")'
#root -b -q 'HealthCheck_Files.C("/pnfs/dune/scratch/users/hyliao/protodunedata/Beam.root")'
#mv goodfile_list.txt goodfile_list_run$run.txt

#[2]Generate the ana module [to get the dat structure of selected trees]
#g++ makeproton_ana.cc `root-config --libs --cflags` -o makeproton_ana_run$run
g++ makeproton_ana.cc `root-config --libs --cflags` -o $ana_name


#[3]Run the ana module (input can be changable if needed but still need compile to loop over the selected files)
#./$ana_name tmp.txt $selector_name
#./$ana_name list_data_prod4_run$run'.txt' $selector_name
#./$ana_name file_run5387_reco1_old.txt $selector_name
#./$ana_name file_run5387_reco2.txt $selector_name
#./$ana_name file_run5387_reco2_new.txt $selector_name
#./$ana_name file_run5387_reco2_new2.txt $selector_name
#./$ana_name file_run5387_prod4reco2.txt $selector_name
./$ana_name file_run5387_reco2_new3.txt $selector_name

#Testing Ajib's ACA map
#ACA without interpolation
#./$ana_name ajib_acamap_sceoff.txt $selector_name
#ACA with interpolation
#./$ana_name ajib_acamap.txt $selector_name


#[4]Fix bugs in the generated makeclass module
#[4.1]Insert one line to make MakeClass work
#sed '/Init(tree)\;/i if (tree-\>InheritsFrom(\"TChain\")) ((TChain\*)tree)-\>LoadTree(0);' ProtonSelector_run${run}.h > ProtonSelector_0_run${run}.h
#mv ProtonSelector_0_run${run}.h ProtonSelector_run${run}.h

sed '/Init(tree)\;/i if (tree-\>InheritsFrom(\"TChain\")) ((TChain\*)tree)-\>LoadTree(0);' $header_name > $tmp_header_name
mv $tmp_header_name $header_name


#[4.2]Fix bug of GetEntry function (not necessary process if we a template already, which is the case for the moment)
#sed 's/GetEntriesFast/GetEntries/g' ProtonSelector_run${run}.C > ProtonSelector_0_run${run}.C

#[4.3]copy an existing code and replace the string to the selected run
#cp -prv ProtonSelector_runX.C  ProtonSelector_0_run${run}.C
#sed 's/5387/'${run}'/g' ProtonSelector_0_run${run}.C > ProtonSelector_run${run}.C
cp -prv $class_namex".C"  $tmp_class_code
sed 's/5387/'${run}'/g' $tmp_class_code > $class_code
rm -f $tmp_class_code

#cp -prv ProtonSelector_run5303.C  ProtonSelector_0_run${run}.C
#sed 's/5303/'${run}'/g' ProtonSelector_0_run${run}.C > ProtonSelector_run${run}.C

#cp -prv ProtonSelector_run5308.C  ProtonSelector_0_run${run}.C
#sed 's/5308/'${run}'/g' ProtonSelector_0_run${run}.C > ProtonSelector_run${run}.C

#[5]Run analysis code
#root -b -q 'RunAna.C('${run}')'

root_exe_str="root -b -q 'RunAna.C(\""$selector_name\"")'"

echo $root_exe_str" ......"
eval $root_exe_str

