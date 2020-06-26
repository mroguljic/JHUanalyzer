#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
xrdcp root://cmseos.fnal.gov//store/user/lcorcodi/10XwithNano.tgz ./
export SCRAM_ARCH=slc6_amd64_gcc700
scramv1 project CMSSW CMSSW_10_2_13
tar xzf 10XwithNano.tgz
rm 10XwithNano.tgz
cd CMSSW_10_2_13/src
eval `scramv1 runtime -sh`
git clone https://github.com/cmantill/JHUanalyzer.git -b cmantillHbb
cd ../..

mkdir tardir; cp tarball.tgz tardir/; cd tardir
tar xzvf tarball.tgz
echo "I AM HERE"
pwd
cp -r * ../CMSSW_10_2_13/src/JHUanalyzer/
cd ../CMSSW_10_2_13/src/JHUanalyzer/
eval `scramv1 runtime -sh`

echo make_preselection_HbbX.py $*
echo $*
python make_preselection_HbbX.py $*
ls -lrth
cp Hbbpreselection*.root ../../../

