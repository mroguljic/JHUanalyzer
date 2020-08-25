#################################################################
# make_preselection_HbbX.py - For Zbb selection #
# -----------------------------------------------------------   #
# Reads the jetNano trees on EOS, builds the 2D distributions   #
#################################################################

import ROOT
from ROOT import *

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object,Event
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import InputTree
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetSysColl import JetSysColl, JetSysObj
from PhysicsTools.NanoAODTools.postprocessing.framework.preskimming import preSkim

import pickle
from optparse import OptionParser
import copy
import math
from math import sqrt
import sys
import time

import Presel_Functions
from Presel_Functions import *

if __name__ == "__main__":
    
    parser = OptionParser()

    parser.add_option('-s', '--set', metavar='F', type='string', action='store',
                    default   =   'data',
                    dest      =   'set',
                    help      =   'dataset (ie data,ttbar etc)')
    parser.add_option('-r', '--region', metavar='F', type='string', action='store',
                    default   =   'default',
                    dest      =   'region',
                    help      =   'default, sideband, ttbar')
    parser.add_option('-y', '--year', metavar='FILE', type='string', action='store',
                    default   =   '',
                    dest      =   'year',
                    help      =   'Year (16,17,18)')
    parser.add_option('-J', '--JES', metavar='F', type='string', action='store',
                    default   =   'nominal',
                    dest      =   'JES',
                    help      =   'nominal, up, or down')
    parser.add_option('-R', '--JER', metavar='F', type='string', action='store',
                    default   =   'nominal',
                    dest      =   'JER',
                    help      =   'nominal, up, or down')
    parser.add_option('-a', '--JMS', metavar='F', type='string', action='store',
                    default   =   'nominal',
                    dest      =   'JMS',
                    help      =   'nominal, up, or down')
    parser.add_option('-b', '--JMR', metavar='F', type='string', action='store',
                    default   =   'nominal',
                    dest      =   'JMR',
                    help      =   'nominal, up, or down')
    parser.add_option('-j', '--job', metavar='F', type='string', action='store',
                    default   =   'all',
                    dest      =   'job',
                    help      =   'job number')
    parser.add_option('-n', '--njobs', metavar='F', type='string', action='store',
                    default   =   '1',
                    dest      =   'njobs',
                    help      =   'number of jobs')
    parser.add_option('-d', '--doublebtagger', metavar='F', type='string', action='store',
                    default   =   'btagHbb',
                    dest      =   'doublebtagger',
                    help      =   'Variable name in NanoAOD for double b tagger to be used. btagHbb (default), deepTagMD_HbbvsQCD, deepTagMD_ZHbbvsQCD, btagDDBvL')
    
    print("IMPORTING ROOT FROM: ")
    print(ROOT.__file__)

    (options, args) = parser.parse_args()
    print('run with options')

    ################################
    # Setup double b tagger option #
    ################################
    doubleB_titles = {'btagHbb':'Double b',
                      'deepTagMD_HbbvsQCD': 'DeepAK8 MD Hbb',
                      'deepTagMD_ZHbbvsQCD': 'DeepAK8 MD ZHbb',
                      'btagDDBvL': 'Deep Double b'}
    doubleB_abreviations = {'btagHbb':'doubleB',
                      'deepTagMD_HbbvsQCD': 'dak8MDHbb',
                      'deepTagMD_ZHbbvsQCD': 'dak8MDZHbb',
                      'btagDDBvL': 'DeepDB'}
    doubleB_name = options.doublebtagger
    doubleB_title = doubleB_titles[doubleB_name]
    doubleB_short = doubleB_abreviations[doubleB_name]

    ################################
    # data #
    ################################
    dataList = ['data','SinglePhoton']
    if any(x in options.set for x in dataList):
        isData = True
    else:
        isData = False

    ######################################
    # Make strings for final file naming #
    ######################################
    # JECs
    runOthers = True
    if not isData:
        mass_name = ''
    else:
        mass_name = '_nom'
    mod = ''
    if options.JES!='nominal':
        mod = '_JES' + '_' + options.JES
        mass_name = '_jesTotal'+options.JES.capitalize()
        runOthers = False
    if options.JER!='nominal':
        mod = '_JER' + '_' + options.JER
        mass_name = '_jer'+options.JER.capitalize()
        runOthers = False
    if options.JMS!='nominal':
        mod = '_JMS' + '_' + options.JMS
        mass_name = '_jms'+options.JMS.capitalize()
        runOthers = False
    if options.JMR!='nominal':
        mod = '_JMR' + '_' + options.JMR
        mass_name = '_jmr'+options.JMR.capitalize()
        runOthers = False

    #######################
    # Setup job splitting #
    #######################
    if int(options.job) > int(options.njobs):
        sys.stdout.write('Trying to run job '+options.job+' out of '+options.njobs)
        sys.stdout.write(options.set)
        print('ERROR: Trying to run job '+options.job+' out of '+options.njobs,options.set)
        raise RuntimeError('ERROR: Trying to run job '+options.job+' out of '+options.njobs)
    jobs=int(options.njobs)
    if jobs != 1:
        num=int(options.job)
        jobs=int(options.njobs)
        print "Running over " +str(jobs)+ " jobs"
        print "This will process job " +str(num)
    else:
        print "Running over all events"

    #################################
    # Load cut values and constants #
    #################################
    Cons = LoadConstants(options.year)
    lumi = Cons['lumi']

    Cuts = LoadCuts(options.region,options.year)

    #################################
    # Load N2ddt distribution
    #################################
    if options.year=='16':
        f_h2ddt = ROOT.TFile.Open("weights/GridOutput_v13_WP026.root", "read")
    elif options.year=='17':
        f_h2ddt = ROOT.TFile.Open("weights/Output_smooth_2017MC.root", "read")
    elif options.year=='2018':
        f_h2ddt = ROOT.TFile.Open("weights/n2ddtmap_2018bits_GaussianSmoothing1Sigma_CorrectVersion.root", "read")
    else:
        f_h2ddt = ROOT.TFile.Open("weights/GridOutput_v13_WP026.root", "read")
    trans_h2ddt = f_h2ddt.Get("Rho2D")
    trans_h2ddt.SetDirectory(0)
    f_h2ddt.Close()

    f_h2ddt = ROOT.TFile.Open("weights/n2ddt_hh%s.root"%options.year,"read")
    trans_h2ddt_new = f_h2ddt.Get("Rho2D")
    trans_h2ddt_new.SetDirectory(0)
    f_h2ddt.Close()

    ################################# 
    # Load mSD corrections
    ################################# 
    puppisd_corrGEN = ROOT.TF1("corrGEN", "[0]+[1]*pow(x*[2],-[3])", 200, 3500)
    puppisd_corrGEN.SetParameter(0,1.00626)
    puppisd_corrGEN.SetParameter(1, -1.06161)
    puppisd_corrGEN.SetParameter(2,0.0799900)
    puppisd_corrGEN.SetParameter(3,1.20454)
    puppisd_corrRECO_cen = ROOT.TF1("corrRECO_cen", "[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",200, 3500)
    puppisd_corrRECO_cen.SetParameter(0,1.09302)
    puppisd_corrRECO_cen.SetParameter(1,-0.000150068)
    puppisd_corrRECO_cen.SetParameter(2,3.44866e-07)
    puppisd_corrRECO_cen.SetParameter(3,-2.68100e-10)
    puppisd_corrRECO_cen.SetParameter(4,8.67440e-14)
    puppisd_corrRECO_cen.SetParameter(5,-1.00114e-17)
    puppisd_corrRECO_for = ROOT.TF1("corrRECO_for", "[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",200, 3500)
    puppisd_corrRECO_for.SetParameter(0,1.27212)
    puppisd_corrRECO_for.SetParameter(1,-0.000571640)
    puppisd_corrRECO_for.SetParameter(2,8.37289e-07)
    puppisd_corrRECO_for.SetParameter(3,-5.20433e-10)
    puppisd_corrRECO_for.SetParameter(4,1.45375e-13)
    puppisd_corrRECO_for.SetParameter(5,-1.50389e-17)
    
    ##########################################################
    # Load weights if not data #
    ##########################################################
    if not isData:
        # need to load singlephoton triggs eff too!
        if options.year=='16':
            effFile = 'trigger/RUNTriggerEfficiencies_SingleMuon_Run2016_V2p1_v03.root'
        elif options.year=='17':
            effFile = 'trigger/TrigEff_2017BtoF_noPS_Feb21.root'
        elif options.year=='18':
            effFile = 'trigger/TrigEff_2018_Feb21.root'
        else: 
            effFile = 'trigger/RUNTriggerEfficiencies_SingleMuon_Run2016_V2p1_v03.root'
        print "Using triggerEff file = ",effFile
        TrigFile = TFile.Open(effFile,"read")
        if options.year=='17' or options.year=='18':
            trig_denom = TrigFile.Get("h_denom"); trig_denom.SetDirectory(0)
            trig_numer = TrigFile.Get("h_numer"); trig_numer.SetDirectory(0)
        if options.year=='16':
            trig_denom = TrigFile.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtDenom_cutJet"); trig_denom.SetDirectory(0); trig_denom.RebinX(2); trig_denom.RebinY(5);
            trig_numer = TrigFile.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtPassing_cutJet"); trig_numer.SetDirectory(0); trig_numer.RebinX(2); trig_numer.RebinY(5);
        TrigPlot = TEfficiency()
        if (TEfficiency.CheckConsistency(trig_numer, trig_denom)):
            TrigPlot = TEfficiency(trig_numer, trig_denom)
            TrigPlot.SetDirectory(0)
        TrigFile.Close()

        # NLO k-factors
        f_kfactors = TFile.Open("weights/kfactors.root","read")
        hQCD_Z = f_kfactors.Get('ZJets_012j_NLO/nominal').Clone('QCD_Z')
        hQCD_W = f_kfactors.Get('WJets_012j_NLO/nominal').Clone('QCD_W')
        hLO_Z = f_kfactors.Get('ZJets_LO/inv_pt').Clone('LO_Z')
        hLO_W = f_kfactors.Get('WJets_LO/inv_pt').Clone('LO_W')
        hEWK_Z = f_kfactors.Get('EWKcorr/Z').Clone('EWK_Z')
        hEWK_W = f_kfactors.Get('EWKcorr/W').Clone('EWK_W')
        hQCD_Z.SetDirectory(0)
        hQCD_W.SetDirectory(0)
        hLO_Z.SetDirectory(0)
        hLO_W.SetDirectory(0)
        hEWK_Z.SetDirectory(0)
        hEWK_W.SetDirectory(0)
        f_kfactors.Close()
        hEWK_Z.Divide(hQCD_Z);
        hEWK_W.Divide(hQCD_W);
        hQCD_Z.Divide(hLO_Z);
        hQCD_W.Divide(hLO_W);

        f_ZNLO = ROOT.TFile.Open(os.path.expandvars("weights/ZJets_QCD_NLO.root"), "read")
        if options.year=='17' or options.year=='18':
            znlo = f_ZNLO.Get("Z_NLO_QCD_2017")
        elif options.year=='16':
            znlo = f_ZNLO.Get("Z_NLO_QCD_2016")
        znlo.SetDirectory(0)
        f_ZNLO.Close()

        f_WNLO = ROOT.TFile.Open(os.path.expandvars("weights/WJets_QCD_NLO.root"), "read")
        if options.year=='17' or options.year=='18':
            wnlo = f_WNLO.Get("W_NLO_QCD_2017")
        elif options.year=='16':
            wnlo = f_WNLO.Get("W_NLO_QCD_2016")
        wnlo.SetDirectory(0)
        f_WNLO.Close()

        histsk = {}
        histsk['hEWK_Z'] = hEWK_Z.Clone();
        histsk['hEWK_W'] = hEWK_W.Clone();
        histsk['znlo'] = znlo.Clone();
        histsk['wnlo'] = wnlo.Clone();
        for key in histsk.keys():
            histsk[key].SetDirectory(0);


    names = {
        'QCDHT700': 'qcd',
        'QCDHT1000': 'qcd',
        'QCDHT1500': 'qcd',
        'QCDHT2000': 'qcd',
        'QCD': 'qcd',
        'Zqq': 'zqq',
        'Wqq': 'wqq',
        'ttbar': 'tqq',
        'ttbar-semilep': 'tqq', # this is missing semilep!?
        'singletop': 'stqq',
        'singleantitop': 'stqq',
        'singletWtop': 'stqq',
        'singletWantitop': 'stqq',
        'ZZdiboson': 'vvqq',
        'WWdiboson': 'vvqq',
    }

    #############################
    # Make new file for storage #
    #############################
    if jobs!=1:
        f = TFile( "/afs/cern.ch/user/m/mrogulji/store/matej/wp_0.8/Hbbpreselection"+options.year+"_"+options.set+"_job"+options.job+"of"+options.njobs+mod+'_'+doubleB_short+'_'+options.region+".root", "recreate" )
    else:
        f = TFile( "/afs/cern.ch/user/m/mrogulji/store/matej/wp_0.8/Hbbpreselection"+options.year+"_"+options.set+mod+'_'+doubleB_short+'_'+options.region+".root", "recreate" )
    f.cd()

    ###################
    # Book histograms #
    ###################
    if not isData:
        hname = names[options.set.replace('ext','')]
    else:
        hname = 'data_obs'
    Hbb_cutflow = TH1D('%s_cutflow'%hname, 'Hbb_cutflow', 12, 0.5, 12.5)
    Hbb_cutflow.GetXaxis().SetBinLabel(1, "no cuts")
    Hbb_cutflow.GetXaxis().SetBinLabel(2, "p_{T}")
    Hbb_cutflow.GetXaxis().SetBinLabel(3, "rho")
    Hbb_cutflow.GetXaxis().SetBinLabel(4, "Lepton Veto")
    Hbb_cutflow.GetXaxis().SetBinLabel(5, "ttbar cut")
    Hbb_cutflow.GetXaxis().SetBinLabel(6, "MET<140")
    Hbb_cutflow.GetXaxis().SetBinLabel(7, "n2ddt<0")
    Hbb_cutflow.GetXaxis().SetBinLabel(8, "n2ddt_new<0")
    Hbb_cutflow.GetXaxis().SetBinLabel(9, "tau21ddt<0.43")
    Hbb_cutflow.GetXaxis().SetBinLabel(10, "L - "+doubleB_title)
    Hbb_cutflow.GetXaxis().SetBinLabel(11, "M - "+doubleB_title)
    Hbb_cutflow.GetXaxis().SetBinLabel(12, "T - "+doubleB_title)

    N2_map = TH3F('%s_n2ddt_map'%hname,'%s_n2ddt_map'%hname,
                  30, -6, -1.5,
                  100, 400, 1200,
                  200, 0, 0.5)
    N2_map.Sumw2()
    Hbb_doubleB = TH1F('%s_doubleB'%hname,''+doubleB_title+' tag',60,0,1)
    Hbb_doubleB.Sumw2()
    Hbb_deepAK8 = TH1F('%s_deepAK8'%hname,'deepAK8', 60,0,1)
    Hbb_deepAK8.Sumw2()
    Hbb_pT = TH1F('%s_pT'%hname,'hbb pT',400,200,1000)
    Hbb_pT.Sumw2()
    Hbb_pT_450 = TH1F('%s_pT_450'%hname,'hbb pT 450',100,400,1000)
    Hbb_pT_450.Sumw2()
    Hbb_mSD = TH1F('%s_mSD'%hname,'hbb mSD',60,40,220)
    Hbb_mSD.Sumw2()
    Hbb_mSD_450 = TH1F('%s_mSD_450'%hname,'hbb mSD 450',60,40,220)
    Hbb_mSD_450.Sumw2()
    Hbb_mSD_uncorr =  TH1F('%s_mSD_uncorr'%hname,'hbb mSD uncorr',60,40,220)
    Hbb_mSD.Sumw2()
    Hbb_eta = TH1F('%s_eta'%hname,'hbb eta',60,-3,3)
    Hbb_eta.Sumw2()
    Hbb_rho = TH1F('%s_rho'%hname,'hbb rho',60,-7,-1)
    Hbb_rho.Sumw2()
    Hbb_n2ddt = TH1F('%s_n2ddt'%hname,'hbb n2ddt',60,-0.5,0.5)
    Hbb_n2ddt.Sumw2()
    Hbb_n2ddt_new = TH1F('%s_n2ddt_new'%hname,'hbb n2ddt new',60,-0.5,0.5)
    Hbb_n2ddt_new.Sumw2()
    Hbb_n2 = TH1F('%s_n2'%hname,'hbb n2',60,0.,0.8)
    Hbb_n2.Sumw2()
    Hbb_tau21ddt = TH1F('%s_tau21ddt'%hname,'hbb tau21ddt',60,0,1.2)
    Hbb_tau21ddt.Sumw2()
    Hbb_tau21 = TH1F('%s_tau21'%hname,'hbb tau21',60,0,1)
    Hbb_tau21.Sumw2()
    MET = TH1F('%s_MET'%hname,'MET',60,0,500)
    MET.Sumw2()
    ttbarvar = TH1F('%s_ttbarvar'%hname,'Max Opp. Hem. AK4 DeepCSV b',60,0,1)
    ttbarvar.Sumw2()

    nev = TH1F("nev",   "nev",      1, 0, 1 )

    Mass_binBoundaries = []
    for i0 in range(0, 24):
        Mass_binBoundaries.append(40. + i0 * 7)
    Pt_binBoundaries= [450, 500, 550, 600, 675, 800, 1200]

    passfail = ['pass','fail']
    WPs = ['M','T','L']
    hists = {}
    for wp in WPs:
        hists[wp] = {}
        for cat in passfail:
            histname = "%s_%s%s"%(hname,wp,cat)
            histtitle = "mass of Hbb vs pT - %s %s"%(wp,cat)
            if options.JES != 'nominal': 
                histname += '_JES%s'%options.JES.capitalize()
                histtitle += '_JES%s'%options.JES.capitalize()
            if options.JER != 'nominal': 
                histname += '_JER%s'%options.JER.capitalize()
                histtitle += '_JER%s'%options.JER.capitalize()
            hists[wp][cat] = TH2F(histname, histtitle, len(Mass_binBoundaries)-1, array('d', Mass_binBoundaries), len(Pt_binBoundaries)-1, array('d', Pt_binBoundaries) )

    if runOthers == True:
        if not isData:
            # Pu, Trigger and matched/unmatched
            systs = ['triggerUp','triggerDown','PuUp','PuDown','matched','unmatched','semimatched']
            for wp in WPs:
                for sys in systs:
                    for cat in passfail:
                        hists[wp]['%s_%s'%(cat,sys)] = TH2F("%s_%s%s_%s"%(hname,wp,cat,sys), "mass of Hbb vs pT - %s %s %s"%(cat,sys,wp), len(Mass_binBoundaries)-1, array('d', Mass_binBoundaries), len(Pt_binBoundaries)-1, array('d', Pt_binBoundaries) )
                        hists[wp]['%s_%s'%(cat,sys)].Sumw2()  

    print("Histograms booked")
    ###############################
    # Grab root file that we want #
    ###############################
    file_string = Load_jetNano(options.set,options.year)
    file = TFile.Open(file_string)
    print("root file"+file_string+" loaded")

    ################################
    # Grab event tree from nanoAOD #
    ################################
    inTree = file.Get("Events")
    elist,jsonFilter = preSkim(inTree,None,'')
    inTree = InputTree(inTree,elist)
    treeEntries = inTree.entries
    print("event tree loaded")

    #############################
    # Get process normalization #
    #############################
    norm_weight = 1
    if not isData:
        runs_tree = file.Get("Runs")
        nevents_gen = 0
        
        for i in runs_tree:
            nevents_gen+=i.genEventCount

        xsec = Cons[options.set.replace('ext','')+'_xsec']
        norm_weight = lumi*xsec/float(nevents_gen)

    #####################################
    # Design the splitting if necessary #
    #####################################
    if jobs != 1:
        evInJob = int(treeEntries/jobs)
        lowBinEdge = evInJob*(num-1)
        highBinEdge = evInJob*num
        if num == jobs:
            highBinEdge = treeEntries
    else:
        lowBinEdge = 0
        highBinEdge = treeEntries

    #highBinEdge = 10000 # !!! warning for testing
    print "Range of events: (" + str(lowBinEdge) + ", " + str(highBinEdge) + ")"

    count = 0
    
    ##############
    # Begin Loop #
    ##############
    start = time.time()
    last_event_time = start
    event_time_sum = 0
    for entry in range(lowBinEdge,highBinEdge):        
        count   =   count + 1
        if count % 10000 == 0 :
            print  '--------- Processing Event ' + str(count) +'   -- percent complete ' + str(100*count/(highBinEdge-lowBinEdge)) + '% -- '

        # Grab the event
        event = Event(inTree, entry)

        # grab collections
        ak8JetsColl = Collection(event, "FatJet")
        ak4JetsColl = Collection(event, "Jet")
        electronColl = Collection(event, "Electron")
        muonColl = Collection(event, "Muon")
        tauColl = Collection(event, "Tau")
        ak8subJetsColl = Collection(event,"SubJet")
        if not isData:
            GenParticles = Collection(event,'GenPart')

        # Apply triggers first
        if isData: 
            if 'data' in options.set:
                if options.year=='16':
                    triggers = ['HLT_PFHT800','HLT_PFHT900','HLT_AK8PFJet360_TrimMass30','HLT_AK8PFHT700_TrimR0p1PT0p03Mass50','HLT_PFHT650_WideJetMJJ950DEtaJJ1p5','HLT_PFHT650_WideJetMJJ900DEtaJJ1p5','HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20','HLT_PFJet450']
                if options.year=='17':
                    triggers = ['HLT_AK8PFJet330_PFAK8BTagCSV_p17','HLT_PFHT1050','HLT_AK8PFJet400_TrimMass30','HLT_AK8PFJet420_TrimMass30','HLT_AK8PFHT800_TrimMass50','HLT_PFJet500','HLT_AK8PFJet500']
                if options.year=='18':
                    triggers = ['HLT_AK8PFJet400_TrimMass30','HLT_AK8PFJet420_TrimMass30','HLT_AK8PFHT800_TrimMass50','HLT_PFHT1050','HLT_PFJet500','HLT_AK8PFJet500','HLT_AK8PFJet330_PFAK8BTagCSV_p17','HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4']
            elif 'SinglePhoton' in options.set:
                triggers = ['HLT_Photon175','HLT_Photon165_HE10']
            else: 
                triggers = []
    
            isTriggered = False
            for t in triggers:
                if hasattr(event,t): 
                    isTriggered = True
            if not isTriggered:
                continue

        # Filters
        filters = [event.Flag_goodVertices,
                   event.Flag_HBHENoiseFilter,
                   event.Flag_HBHENoiseIsoFilter,
                   event.Flag_globalTightHalo2016Filter,
                   event.Flag_EcalDeadCellTriggerPrimitiveFilter,
                   event.Flag_eeBadScFilter]
        filterFails = 0
        for thisFilter in filters:
            if thisFilter == 0:
                filterFails += 1
        if filterFails > 0:
            continue

        # Selections
        Hbbsel = {}

        # Now jetID which (in binary #s) is stored with bit1 as loose, bit2 as tight, and filters (after grabbing jet coLections)
        if len(ak8JetsColl) >= 1 and (ak8JetsColl[0].jetId & 2 == 2):
            Hbbsel['jetIds'] = True
        else:
            continue

        # check if we have enough jets
        if len(ak8JetsColl) >= 1:
            Hbbsel['nFatJet'] = True
        elif len(ak8JetsColl) < 1:
            continue

        # lepton veto
        Hbbsel['leptonExists'] = False
        if len(electronColl) > 0:
            for e in electronColl:
                if e.cutBased >= 2:
                    Hbbsel['leptonExists'] = True
        if len(muonColl) > 0:
            for m in muonColl:
                if m.looseId:
                    Hbbsel['leptonExists'] = True
        if len(tauColl) > 0:
            for t in tauColl:
                if t.idMVAnewDM2017v2 >= 4:
                    Hbbsel['leptonExists'] = True

        # Start preselections
        # pt cut
        Hbbsel['pT'] = ak8JetsColl[0].pt_nom > 450
        Hbbpt = ak8JetsColl[0].pt_nom

        # mass selection
        if isData:
            Hbbmsd = ak8JetsColl[0].msoftdrop_nom
        else:
            Hbbmsd = CorrectMSD(ak8JetsColl[0],ak8subJetsColl,puppisd_corrGEN,puppisd_corrRECO_cen,puppisd_corrRECO_for)
        if Hbbmsd < 20: continue
        Hbbsel['msd'] = Hbbmsd > 40

        # eta - this is ok as a cut
        if abs(ak8JetsColl[0].eta) > 2.4: continue

        # preselection that is not inmediately applied
        # rho selection
        if ak8JetsColl[0].pt_nom >= 200:
            Hbbrho = 2*log(Hbbmsd/ak8JetsColl[0].pt_nom)
            Hbbsel['rho'] = -6 < Hbbrho < -2.1
        else:
            Hbbrho = -1
            Hbbsel['rho'] = False

        # check ttbar veto based on topology of ak4 jets
        # no medium DeepCSV b-tagged jet in the opposite hemisphere of the AK8 jet among the four leading-pT AK4 jets
        candidateAK4s = Hemispherize(ak8JetsColl,ak4JetsColl) # get opposite hemisphere ak4 jets
        Hbbsel['TTbarCut'] = True
        maxdeepcsv = 0
        if len(candidateAK4s)>0:
            for i in range(min(4,len(candidateAK4s))):
                if candidateAK4s[i].btagDeepB > maxdeepcsv:
                    maxdeepcsv = candidateAK4s[i].btagDeepB
            #if Cuts['deepbtag'][0] < maxdeepcsv < Cuts['deepbtag'][1]: # check if any of the leading 4 with deepbtag within that range
            #0.2219 0.6324 0.8958
            if 0.2219 < maxdeepcsv < 1:  #tight?
            #if 0.6324 < maxdeepcsv < 1: #medium
                Hbbsel['TTbarCut'] = False

        # MET cut
        Hbbsel['MET'] = inTree.readBranch('MET_pt') < 140

        # N2ddt cut
        n2ddt = getN2dtt(trans_h2ddt,ak8JetsColl[0].n2b1,Hbbrho,Hbbpt)
        n2ddt_new = getN2dtt(trans_h2ddt_new,ak8JetsColl[0].n2b1,Hbbrho,Hbbpt)
        Hbbsel['n2ddt'] = n2ddt < 0
        Hbbsel['n2ddt_new'] = n2ddt_new < 0

        # tau21ddt cut
        # old t21ddt: t21+0.063*math.log(msd**2/pt)
        # jmar t21ddt: t21+0.080*math.log(msd**2/pt)
        tau21 = (ak8JetsColl[0].tau2/ak8JetsColl[0].tau1)
        tau21ddt = tau21 + 0.080*math.log((Hbbmsd**2)/Hbbpt)
        Hbbsel['tau21ddt'] = tau21ddt < Cuts['tau21ddt'][1]

        jet0 = TLorentzVector()
        jet0.SetPtEtaPhiM(ak8JetsColl[0].pt_nom, ak8JetsColl[0].eta, ak8JetsColl[0].phi, Hbbmsd)        
        # jet0.SetPt(ak8JetsColl[0].pt_nom)
        # jet0.SetEta(ak8JetsColl[0].eta)
        # jet0.SetPhi(ak8JetsColl[0].phi)
        # jet0.SetM(Hbbmsd)
        # match jet
        Hbbsel['unmatched'] = False
        Hbbsel['semimatched'] = False
        Hbbsel['matched'] = False
        if not isData:
            if 'zqq' in names[options.set.replace('ext','')]: vid = 23
            if 'wqq' in names[options.set.replace('ext','')]: vid = 24
            if 'tqq' in names[options.set.replace('ext','')]:
                matched = TopJetMatching(jet0, GenParticles)
                genVPt = 0
                if matched==3:
                    Hbbsel['matched'] = True
                elif matched==2:
                    Hbbsel['semimatched'] = True
                else:
                    Hbbsel['unmatched'] = True
            else:
                matched,genVPt = VJetMatching(jet0,GenParticles,vid)
                if matched:
                    Hbbsel['matched'] = True
                else:
                    Hbbsel['unmatched'] = True

            '''
            dphi = 9999; dpt= 9999; dmass = 9999;
            if ak8GenJetsColl[0].pt > 0 and ak8GenJetsColl[0].mass > 0:
                dphi = abs(deltaPhi(ak8GenJetsColl[0].phi , jet0.Phi()))
                dpt = math.fabs(ak8GenJetsColl[0].pt - jet0.Pt()) / ak8GenJetsColl[0].pt
                dmass = math.fabs(ak8GenJetsColl[0].mass - jet0.M()) / ak8GenJetsColl[0].mass
            if dphi < 0.8 and dpt < 0.5 and dmass < 0.3:
                Hbbsel['matched'] = True
            else:
                Hbbsel['unmatched'] = True
            '''

        # bb selection
        Hbbsel['DoubleB_lead_tight'] = (Cuts['doublebtagTight'][0] < getattr(ak8JetsColl[0],doubleB_name) < Cuts['doublebtagTight'][1])
        Hbbsel['DoubleB_lead_medium'] = (Cuts['doublebtag'][0] < getattr(ak8JetsColl[0],doubleB_name) < Cuts['doublebtag'][1])
        Hbbsel['DoubleB_lead_loose'] = (Cuts['doublebtagLoose'][0] < getattr(ak8JetsColl[0],doubleB_name) < Cuts['doublebtagLoose'][1])
        
        # pass
        Hbbsel['pass_T'] = Hbbsel['DoubleB_lead_tight']
        Hbbsel['pass_M'] = Hbbsel['DoubleB_lead_medium']
        Hbbsel['pass_L'] = Hbbsel['DoubleB_lead_loose']

        # fail
        Hbbsel['fail_T'] = not Hbbsel['DoubleB_lead_tight']
        Hbbsel['fail_M'] = not Hbbsel['DoubleB_lead_medium']
        Hbbsel['fail_L'] = not Hbbsel['DoubleB_lead_loose']

        # define preselection
        #preselection = Hbbsel['pT'] and Hbbsel['msd'] and Hbbsel['rho'] and Hbbsel['jetIds'] and not Hbbsel['leptonExists'] and Hbbsel['TTbarCut'] and Hbbsel['MET'] and Hbbsel['tau21ddt'] 
        preselection = Hbbsel['pT'] and Hbbsel['msd'] and Hbbsel['rho'] and Hbbsel['jetIds'] and not Hbbsel['leptonExists'] and Hbbsel['TTbarCut'] and Hbbsel['MET'] and Hbbsel['n2ddt_new'] 
        if not isData:
            Hbb_cutflow.Fill(1)
            if Hbbsel['pT']:
                Hbb_cutflow.Fill(2)
                if Hbbsel['rho']:
                    Hbb_cutflow.Fill(3)
                    if not Hbbsel['leptonExists']:
                        Hbb_cutflow.Fill(4)
                        if Hbbsel['TTbarCut']:
                            Hbb_cutflow.Fill(5)
                            if Hbbsel['MET']:
                                Hbb_cutflow.Fill(6)
                                if Hbbsel['n2ddt']: # this is not applied
                                    Hbb_cutflow.Fill(7)
                                if Hbbsel['tau21ddt']: # this goes in the middle
                                    Hbb_cutflow.Fill(9)
                                if Hbbsel['n2ddt_new']: # this is not applied                                                                                                     
                                    Hbb_cutflow.Fill(8)
                                    if Hbbsel['pass_L']: # this is for loose W/Z templates
                                        Hbb_cutflow.Fill(10)
                                        if Hbbsel['pass_M']:
                                            Hbb_cutflow.Fill(11) # this is the current working point
                                            if Hbbsel['pass_T']:
                                                Hbb_cutflow.Fill(12)
        
        #########################################
        # Weights
        #########################################
        weights = { 'Pu':{},
                    'trigger':{},
                    'kFactor':{},
                }

        if not isData:
            # Pileup reweighting applied         
            weights['Pu']['nom'] = inTree.readBranch('puWeight')
            weights['Pu']['up'] = inTree.readBranch('puWeightUp')
            weights['Pu']['down'] = inTree.readBranch('puWeightDown')

            # Trigger weight
            massForTrig = min(max(Hbbmsd,0), 300.)
            ptForTrig = max(200., min(ak8JetsColl[0].pt_nom, 1000.))
            weights['trigger']['nom'] = Trigger_Lookup( massForTrig, ptForTrig, TrigPlot )[0]
            weights['trigger']['up'] = Trigger_Lookup( massForTrig, ptForTrig, TrigPlot )[1]
            weights['trigger']['down'] = Trigger_Lookup( massForTrig, ptForTrig, TrigPlot )[2]

            # k-factor weight
            weights['kFactor']['nom'] = kFactor_Lookup( genVPt, histsk, names[options.set.replace('ext','')]+'_'+options.year )[0]
            weights['kFactor']['up'] = weights['kFactor']['nom']
            weights['kFactor']['down'] = weights['kFactor']['nom']

        #########################################
        # Fill control histograms #
        #########################################
        if not isData:
            if Hbbmsd > 20 and ak8JetsColl[0].pt_nom > 200 and not Hbbsel['leptonExists']:
                N2_map.Fill(2*log(Hbbmsd/ak8JetsColl[0].pt_nom),ak8JetsColl[0].pt_nom,ak8JetsColl[0].n2b1,norm_weight*Weightify(weights,'nominal'))
        if preselection:
            Hbb_rho.Fill(Hbbrho,norm_weight*Weightify(weights,'nominal'))
            Hbb_pT.Fill(jet0.Pt(),norm_weight*Weightify(weights,'nominal'))
            Hbb_mSD_uncorr.Fill(ak8JetsColl[0].msoftdrop,norm_weight*Weightify(weights,'nominal'))
            Hbb_mSD.Fill(jet0.M(),norm_weight*Weightify(weights,'nominal'))
            Hbb_eta.Fill(jet0.Eta(),norm_weight*Weightify(weights,'nominal'))
            MET.Fill(inTree.readBranch('MET_pt'),norm_weight*Weightify(weights,'nominal'))
            ttbarvar.Fill(maxdeepcsv,norm_weight*Weightify(weights,'nominal'))
            Hbb_doubleB.Fill(getattr(ak8JetsColl[0],doubleB_name),norm_weight*Weightify(weights,'nominal'))
            Hbb_deepAK8.Fill(getattr(ak8JetsColl[0],'deepTagMD_ZHbbvsQCD'),norm_weight*Weightify(weights,'nominal'))
            Hbb_n2.Fill(ak8JetsColl[0].n2b1,norm_weight*Weightify(weights,'nominal'))
            Hbb_n2ddt.Fill(n2ddt,norm_weight*Weightify(weights,'nominal'))
            Hbb_n2ddt_new.Fill(n2ddt_new,norm_weight*Weightify(weights,'nominal'))
            Hbb_tau21.Fill(tau21,norm_weight*Weightify(weights,'nominal'))
            Hbb_tau21ddt.Fill(tau21ddt,norm_weight*Weightify(weights,'nominal'))
            if ak8JetsColl[0].pt_nom > 450:
                Hbb_pT_450.Fill(jet0.Pt(),norm_weight*Weightify(weights,'nominal'))
                Hbb_mSD_450.Fill(jet0.M(),norm_weight*Weightify(weights,'nominal'))

        ######################################### 
        # Check preselection #
        ######################################### 
        if preselection:
            for wp in WPs:
                for cat in ['pass','fail']:
                    if Hbbsel['%s_%s'%(cat,wp)]:
                        hists[wp][cat].Fill(jet0.M(),jet0.Pt(),norm_weight*Weightify(weights,'nominal'))

            # Fill variations
            if runOthers and not isData:
                for wp in WPs:
                    for cat in ['pass','fail']:
                        for sys in ['matched','unmatched','semimatched']:
                            if Hbbsel['%s_%s'%(cat,wp)] and Hbbsel[sys]:
                                mass_corr = jet0.M()
                                if sys=='matched': mass_corr = jet0.M()*Cons['shift_SF']
                                hists[wp]['%s_%s'%(cat,sys)].Fill(mass_corr,jet0.Pt(),norm_weight*Weightify(weights,'nominal'))
                            
                        for sys in ['Pu','trigger']:
                            if Hbbsel['%s_%s'%(cat,wp)]:
                                hists[wp]['%s_%sUp'%(cat,sys)].Fill(jet0.M(),jet0.Pt(),norm_weight*Weightify(weights,'%s_up'%sys))
                                hists[wp]['%s_%sDown'%(cat,sys)].Fill(jet0.M(),jet0.Pt(),norm_weight*Weightify(weights,'%s_down'%sys))

    end = time.time()
    print '\n'
    print str((end-start)/60.) + ' min'
    f.cd()
    f.Write()
    f.Close()
