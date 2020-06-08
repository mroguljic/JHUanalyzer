# JHUanalyzer

Analyzes NanoAOD using a flexible framework to easily change cuts

## For ZbbNutple production (in cmslpc)
The condor workflow is optimized for cmslpc - modify accordingly for lxplus.
Nano skim (and post-processed) ntuples are here:
`/store/group/lpctlbsm/dbrehm/HH4bv5/`

To setup:
```
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src/
cmsenv
git clone https://github.com/cmantill/JHUanalyzer.git -b cmantillHbb
cd JHUanalyzer/
```

Hbb signal region pre-selection script: `make_preselection_HbbX.py`

Hbb muon CR region pre-selection script: `make_preselection_HbbSM.py`

Pre-selection functions and definition of WPs: `Presel_Functions.py`

To submit jobs (for deepAK8) -- the last arguments are the scripts transferred to condor:

For signalRegion:
``` python CondorHelper.py -r condor/run_scripts/run_presel_Zbb.sh -a condor/args/hh_presel_Zbb_args.txt -i "make_preselection_HbbX.py Presel_Functions.py GenParticleChecker.py trigger/ pileup/ weights/"
 python CondorHelper.py -r condor/run_scripts/run_presel_Zbb.sh -a condor/args/hh_presel_Zbb_args_vars.txt -i "make_preselection_HbbX.py Presel_Functions.py GenParticleChecker.py trigger/ pileup/ weights/"
 ```
For muonCR:
```
python CondorHelper.py -r condor/run_scripts/run_presel_ZbbSM.sh -a condor/args/hh_presel_ZbbSM_dAK8_args.txt -i "make_preselection_HbbSM.py Presel_Functions.py GenParticleChecker.py trigger/ pileup/ weights/"
python CondorHelper.py -r condor/run_scripts/run_presel_ZbbSM.sh -a condor/args/hh_presel_ZbbSM_dAK8_args_vars.txt -i "make_preselection_HbbSM.py Presel_Functions.py GenParticleChecker.py trigger/ pileup/ weights/"
```
Arguments of each script can be found in `condor/args/`. The (`*_vars`) files contain variations e.g. JERUp, JERDown, etc.
