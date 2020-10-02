#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
xrdcp root://cmseos.fnal.gov//store/user/lcorcodi/10XwithNano.tgz ./
export SCRAM_ARCH=slc7_amd64_gcc700
scramv1 project CMSSW CMSSW_10_6_14
tar xzf 10XwithNano.tgz
rm 10XwithNano.tgz
cd CMSSW_10_6_14/src
eval `scramv1 runtime -sh`
echo "git clone https://github.com/mroguljic/JHUanalyzer.git -b matejHbb"
git clone https://github.com/mroguljic/JHUanalyzer.git -b matejHbb
cd ../..

mkdir tardir; cp tarball.tgz tardir/; cd tardir
tar xzvf tarball.tgz
echo "I AM HERE"
pwd
cp -r * ../CMSSW_10_6_14/src/JHUanalyzer/
cd ../CMSSW_10_6_14/src/JHUanalyzer/
eval `scramv1 runtime -sh`

echo make_preselection_HbbSM.py $*
echo $*
python make_preselection_HbbSM.py $*
ls -lrth

