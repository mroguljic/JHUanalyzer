#################################################################
# make_preselection_HbbSM.py - For Zbb SM control region selection #
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
    dataList = ['data','SinglePhoton','singlemuon','EGamma']
    if any(x in options.set for x in dataList):
        isData = True
    else:
        isData = False

    ######################################
    # Make strings for final file naming #
    ######################################
    # JECs
    runOthers = True
    mod = ''
    if options.JES!='nominal':
        mod = '_JES' + '_' + options.JES
        runOthers = False
    if options.JER!='nominal':
        mod = '_JER' + '_' + options.JER
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
        # for 2016 only for now

        # Muon trigger efficiency
        lumi_GH = 16.146
        lumi_BCDEF = 19.721
        lumi_total = lumi_GH + lumi_BCDEF
        f_mutrig_BCDEF = ROOT.TFile.Open("trigger/EfficienciesAndSF_RunBtoF.root", "read")
        mutrig_eff_BCDEF = f_mutrig_BCDEF.Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA")
        mutrig_eff_BCDEF.Sumw2()
        mutrig_eff_BCDEF.Scale(lumi_BCDEF/lumi_total)
        
        f_mutrig_GH    = ROOT.TFile.Open("trigger/EfficienciesAndSF_Period4.root", "read") 
        mutrig_eff_GH = f_mutrig_GH.Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA")
        mutrig_eff_GH.Sumw2()
        mutrig_eff_GH.Scale(lumi_GH / lumi_total)
        
        mutrig_eff = mutrig_eff_BCDEF.Clone("pt_abseta_DATA_mutrig_ave")
        mutrig_eff.Add(mutrig_eff_GH)
        mutrig_eff.SetDirectory(0)
        
        # Muon isolation efficiency
        f_muiso_GH = ROOT.TFile.Open("trigger/EfficienciesAndSF_ISO_GH.root", "read")
        muiso_eff_GH = f_muiso_GH.Get("LooseISO_LooseID_pt_eta/efficienciesDATA/pt_abseta_DATA")
        muiso_eff_GH.Sumw2()
        muiso_eff_GH.SetDirectory(0)
        f_muiso_GH.Close()
        
        f_muiso_BCDEF = ROOT.TFile.Open("trigger/EfficienciesAndSF_ISO_BCDEF.root", "read")
        muiso_eff_BCDEF = f_muiso_BCDEF.Get("LooseISO_LooseID_pt_eta/efficienciesDATA/pt_abseta_DATA")
        muiso_eff_BCDEF.Sumw2()
        muiso_eff_BCDEF.SetDirectory(0)
        f_muiso_BCDEF.Close()
        
        muiso_eff = muiso_eff_GH.Clone('pt_abseta_DATA_muiso_ave')
        muiso_eff.Scale(lumi_GH / lumi_total)
        muiso_eff.Add(muiso_eff_BCDEF, lumi_BCDEF / lumi_total)

        #2016 muon ID
        f_muid_GH = ROOT.TFile.Open("trigger/EfficienciesAndSF_GH.root", "read")
        muid_eff_GH = f_muid_GH.Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/pt_abseta_DATA")
        muid_eff_GH.Sumw2()
        muid_eff_GH.SetDirectory(0)
        f_muid_GH.Close()
        f_muid_BCDEF = ROOT.TFile.Open("trigger/EfficienciesAndSF_BCDEF.root", "read")
        muid_eff_BCDEF = f_muid_BCDEF.Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/pt_abseta_DATA")
        muid_eff_BCDEF.Sumw2()
        muid_eff_BCDEF.SetDirectory(0)
        f_muid_BCDEF.Close()
        
        muid_eff = muid_eff_GH.Clone('pt_abseta_DATA_muid_ave')
        muid_eff.Scale(lumi_GH / lumi_total)
        muid_eff.Add(muid_eff_BCDEF, lumi_BCDEF / lumi_total)

    names = {
        'QCDHT700': 'qcd',
        'QCDHT1000': 'qcd',
        'QCDHT1500': 'qcd',
        'QCDHT2000': 'qcd',
        'QCD': 'qcd',
        'Zqq': 'zqq',
        'Wqq': 'wqq',
        'ttbar': 'tqq',
        'singletop': 'stqq',
        'singleantitop': 'stqq',
        'singletWtop': 'stqq',
        'singletWantitop': 'stqq',
        'ZZdiboson': 'vvqq',
        'WWdiboson': 'vvqq',
        'wjetstolnu': 'wlnu',
        'dyjetstoll1jet': 'zll',
    }

    #############################
    # Make new file for storage #
    #############################
    if jobs!=1:
        f = TFile( "/afs/cern.ch/user/m/mrogulji/2020/Zbb/CristinaInstructions/CMSSW_10_2_13/src/JHUanalyzer/condorResults/HbbSM"+options.year+"_"+options.set+"_job"+options.job+"of"+options.njobs+mod+'_'+doubleB_short+'_'+options.region+".root", "recreate" )
    else:
        f = TFile( "/afs/cern.ch/user/m/mrogulji/2020/Zbb/CristinaInstructions/CMSSW_10_2_13/src/JHUanalyzer/condorResults/HbbSM"+options.year+"_"+options.set+mod+'_'+doubleB_short+'_'+options.region+".root", "recreate" )
    f.cd()

    ###################
    # Book histograms #
    ###################
    if not isData:
        hname = names[options.set.replace('ext','')]
    else:
        hname = 'data_obs'
    Hbb_cutflow = TH1D('%s_cutflow'%hname, 'Hbb_cutflow', 10, 0.5, 10.5)
    Hbb_cutflow.GetXaxis().SetBinLabel(1, "no cuts")
    Hbb_cutflow.GetXaxis().SetBinLabel(2, "p_{T}")
    Hbb_cutflow.GetXaxis().SetBinLabel(3, "mass")
    Hbb_cutflow.GetXaxis().SetBinLabel(4, "muon kin")
    Hbb_cutflow.GetXaxis().SetBinLabel(5, "muon dphi")
    Hbb_cutflow.GetXaxis().SetBinLabel(6, "btag")
    Hbb_cutflow.GetXaxis().SetBinLabel(7, "n2ddt_new<0")
    Hbb_cutflow.GetXaxis().SetBinLabel(8, "L - "+doubleB_title)
    Hbb_cutflow.GetXaxis().SetBinLabel(9, "M - "+doubleB_title)
    Hbb_cutflow.GetXaxis().SetBinLabel(10, "T - "+doubleB_title)

    N2_map = TH3F('%s_n2ddt_map'%hname,'%s_n2ddt_map'%hname,
                  30, -6, -1.5,
                  100, 400, 1200,
                  200, 0, 0.5)
    N2_map.Sumw2()
    Hbb_doubleB = TH1F('%s_doubleB'%hname,''+doubleB_title+' tag',40,0,1)
    Hbb_doubleB.Sumw2()
    Hbb_deepAK8 = TH1F('%s_deepAK8'%hname,'deepAK8', 40,0,1)
    Hbb_deepAK8.Sumw2()
    Hbb_pT = TH1F('%s_pT'%hname,'hbb pT',40,200,1000)
    Hbb_pT.Sumw2()
    Hbb_pT_450 = TH1F('%s_pT_450'%hname,'hbb pT 450',40,400,1000)
    Hbb_pT_450.Sumw2()
    Hbb_mSD = TH1F('%s_mSD'%hname,'hbb mSD',20,40,220)
    Hbb_mSD.Sumw2()
    Hbb_mSD_450 = TH1F('%s_mSD_450'%hname,'hbb mSD 450',20,40,220)
    Hbb_mSD_450.Sumw2()
    Hbb_mSD_uncorr =  TH1F('%s_mSD_uncorr'%hname,'hbb mSD uncorr',20,40,220)
    Hbb_mSD.Sumw2()
    Hbb_mSD_passM_n2ddt = TH1F('%s_mSD_passM_n2ddt'%hname,'hbb mSD',20,40,220)
    Hbb_mSD.Sumw2()
    Hbb_mSD_failM_n2ddt = TH1F('%s_mSD_failM_n2ddt'%hname,'hbb mSD',20,40,220)
    Hbb_mSD.Sumw2()
    Hbb_eta = TH1F('%s_eta'%hname,'hbb eta',40,-3,3)
    Hbb_eta.Sumw2()
    Hbb_rho = TH1F('%s_rho'%hname,'hbb rho',40,-7,-1)
    Hbb_rho.Sumw2()
    Hbb_n2ddt_new = TH1F('%s_n2ddt_new'%hname,'hbb n2ddt new',40,-0.5,0.5)
    Hbb_n2ddt_new.Sumw2()
    Hbb_n2 = TH1F('%s_n2'%hname,'hbb n2',40,0.,0.8)
    Hbb_n2.Sumw2()
    Muon_pT = TH1F('%s_muonpT'%hname,'muon pT',40,0,500)
    Muon_pT.Sumw2()
    Muon_eta = TH1F('%s_muoneta'%hname,'muon eta',40,-2.5,2.5)
    Muon_eta.Sumw2()
    Muon_phi = TH1F('%s_muonphi'%hname,'muon phi',40,-4,4)
    Muon_phi.Sumw2()
    Muon_dphi = TH1F('%s_muondphi'%hname,'muon dphi',40,-4,4)
    Muon_dphi.Sumw2()
    Muon_pfRelIso04_all = TH1F('%s_muonpfRelIso04'%hname,'muon reliso 0.4',40,0,1)
    Muon_pfRelIso04_all.Sumw2()
    MET = TH1F('%s_MET'%hname,'MET',60,0,500)
    MET.Sumw2()

    nev = TH1F("nev",   "nev",      1, 0, 1 )

    Mass_binBoundaries = []
    for i0 in range(0, 24):
        Mass_binBoundaries.append(40. + i0 * 7)

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
            hists[wp][cat] = TH1F(histname, histtitle, len(Mass_binBoundaries)-1, array('d', Mass_binBoundaries))

    if runOthers == True:
        if not isData:
            systs = ['mutriggerUp','mutriggerDown','PuUp','PuDown','muisoUp','muisoDown','muidUp','muidDown','matched','unmatched','semimatched']
            for wp in WPs:
                for sys in systs:
                    for cat in passfail:
                        hists[wp]['%s_%s'%(cat,sys)] = TH1F("%s_%s%s_%s"%(hname,wp,cat,sys), "mass of Hbb vs pT - %s %s %s"%(cat,sys,wp), len(Mass_binBoundaries)-1,  array('d', Mass_binBoundaries))
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
    print "Range of events: (" + str(lowBinEdge) + ", " + str(highBinEdge) + ")"

    count = 0
    
    ##############
    # Begin Loop #
    ##############
    start = time.time()
    last_event_time = start
    event_time_sum = 0
    nevents = {'total': 0,
               'pT': 0,
               'msd': 0,
               'jetIds': 0,
               'muonkin': 0,
               'muon': 0,
               'btag': 0,
               'n2ddt_new': 0
           }

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
            triggers = ['HLT_Mu50','HLT_Mu55']
    
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

        # lepton veto except muons
        Hbbsel['leptonExists'] = False
        if len(electronColl) > 0:
            for e in electronColl:
                if e.cutBased >= 2:
                    Hbbsel['leptonExists'] = True
        if len(tauColl) > 0:
            for t in tauColl:
                if t.idMVAnewDM2017v2 >= 4:
                    Hbbsel['leptonExists'] = True

        # check muon
        Hbbsel['muon'] = False
        Hbbsel['muonkin'] = False
        muon_pt = 0
        muon_eta = 0
        muon_phi = 0
        muon_dphi = 0
        muon_pfRelIso04_all = 0
        if len(muonColl) > 0:
            muon_pfRelIso04_all = muonColl[0].pfRelIso04_all
            #if( muonColl[0].looseId and muonColl[0].pfIsoId == 1):
            if( muonColl[0].looseId and muon_pfRelIso04_all < 0.40):
            #if( muonColl[0].looseId and muon_pfRelIso04_all < 0.25):
            #if( muonColl[0].looseId and muon_pfRelIso04_all < 0.15):
                muon_pt = muonColl[0].pt 
                muon_eta = muonColl[0].eta
                muon_phi = muonColl[0].phi
                muon_dphi = abs(deltaPhi(muon_phi,ak8JetsColl[0].phi))
                if(muon_pt > 55) and (abs(muon_eta) < 2.1):
                    Hbbsel['muonkin'] = True
                    if(muon_dphi > 2.0*TMath.Pi()/3.0):
                       Hbbsel['muon'] = True

        # check Jet
        # argh why not nom here?
        # pt cut
        Hbbsel['pT'] = ak8JetsColl[0].pt > 400
        Hbbpt = ak8JetsColl[0].pt

        # mass selection
        #Hbbmsd = ak8JetsColl[0].msoftdrop
        Hbbmsd = CorrectMSD(ak8JetsColl[0],ak8subJetsColl,puppisd_corrGEN,puppisd_corrRECO_cen,puppisd_corrRECO_for)
        if Hbbmsd < 20: continue
        Hbbsel['msd'] = Hbbmsd > 40

        # eta - this is ok as a cut
        if abs(ak8JetsColl[0].eta) > 2.4: continue

        #jet0 = TLorentzVector(); jet0.SetPtEtaPhiM(ak8JetsColl[0].pt_nom, ak8JetsColl[0].eta, ak8JetsColl[0].phi, Hbbmsd)
        jet0 = TLorentzVector(); jet0.SetPtEtaPhiM(ak8JetsColl[0].pt, ak8JetsColl[0].eta, ak8JetsColl[0].phi, Hbbmsd)

        # check AK4 jet
        Hbbsel['btag'] = False
        if len(ak4JetsColl)>0:
            for i in range(0,len(ak4JetsColl)):
                if ak4JetsColl[i].btagDeepB > 0.6324 and ak4JetsColl[i].pt > 50 and abs(ak4JetsColl[i].eta) < 2.5:
                    ak4vec =  TLorentzVector(); ak4vec.SetPtEtaPhiM(ak4JetsColl[i].pt,ak4JetsColl[i].eta,ak4JetsColl[i].phi,ak4JetsColl[i].mass);
                    if(ak4vec.DeltaR(jet0) > 0.8):
                        Hbbsel['btag'] = True

        # rho
        if ak8JetsColl[0].pt >= 200:
            Hbbrho = 2*log(jet0.M()/jet0.Pt())
            Hbbsel['rho'] = -6 < Hbbrho < -2.1
        else:
            Hbbrho = -1
            Hbbsel['rho'] = False

        # N2ddt cut
        n2ddt_new = getN2dtt(trans_h2ddt_new,ak8JetsColl[0].n2b1,Hbbrho,Hbbpt)
        Hbbsel['n2ddt_new'] = n2ddt_new < 0

        # match jet                                                                                                                                                             
        Hbbsel['unmatched'] = False
        Hbbsel['semimatched'] = False
        Hbbsel['matched'] = False
        if not isData:
            if 'zqq' in names[options.set.replace('ext','')]: vid = 23
            if 'wqq' in names[options.set.replace('ext','')]: vid = 24
            if 'tqq' in names[options.set.replace('ext','')]:
                matched = TopJetMatching(jet0, GenParticles)
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
        preselection = (Hbbsel['pT']) and (Hbbsel['msd']) and (Hbbsel['jetIds']) and (Hbbsel['muonkin']) and (Hbbsel['muon']) and (Hbbsel['btag']) #and (Hbbsel['n2ddt_new'])

        Hbb_cutflow.Fill(1)
        if Hbbsel['pT']:
            Hbb_cutflow.Fill(2)
            nevents['pT'] +=1
            if Hbbsel['msd']:
                Hbb_cutflow.Fill(3)
                nevents['msd'] +=1
                if Hbbsel['jetIds']:
                    nevents['jetIds'] +=1
                    if Hbbsel['muonkin']:
                        Hbb_cutflow.Fill(4)   
                        nevents['muonkin']+=1
                        if Hbbsel['muon']:
                            Hbb_cutflow.Fill(5)
                            nevents['muon']+=1
                            if Hbbsel['btag']:
                                Hbb_cutflow.Fill(6)
                                nevents['btag']+=1
                                if Hbbsel['n2ddt_new']:
                                    Hbb_cutflow.Fill(7)
                                    nevents['n2ddt_new']+=1
                                    if Hbbsel['pass_L']:
                                        Hbb_cutflow.Fill(8)
                                        #print('loose')
                                        if Hbbsel['pass_M']:
                                            Hbb_cutflow.Fill(9) # this is the current working point
                                            nevents['total']+=1
                                            #print('medium')
                                            if Hbbsel['pass_T']:
                                                Hbb_cutflow.Fill(10)
        
        #########################################
        # Weights
        #########################################
        weights = { 'Pu':{},
                    'mutrigger':{},
                    'muid':{},
                    'muiso':{},
                }

        if not isData:
            # Pileup reweighting applied         
            weights['Pu']['nom'] = inTree.readBranch('puWeight')
            weights['Pu']['up'] = inTree.readBranch('puWeightUp')
            weights['Pu']['down'] = inTree.readBranch('puWeightDown')

            # Trigger weight
            if muon_pt>0 and abs(muon_eta)>0:
                muPtForTrig = muon_pt
                muEtaForTrig = abs(muon_eta)
                muTrigs = Hist2d_Lookup(muPtForTrig,muEtaForTrig,mutrig_eff)

                weights['mutrigger']['nom'] = muTrigs[0]
                weights['mutrigger']['up'] = muTrigs[1]
                weights['mutrigger']['down'] = muTrigs[2]

                muPtForId = max(20., min(muon_pt, 100.))
                muEtaForId = min(abs(muon_eta), 2.3)
                muIds = Hist2d_Lookup(muPtForId,muEtaForId,muid_eff)
                weights['muid']['nom'] = muIds[0]
                weights['muid']['up'] = muIds[1]
                weights['muid']['down'] = muIds[2]

                muIsos = Hist2d_Lookup(muPtForId,muEtaForId,muiso_eff)
                weights['muiso']['nom'] = muIsos[0]
                weights['muiso']['up'] = muIsos[1]
                weights['muiso']['down'] = muIsos[2]
            else:
                weights['mutrigger']['nom'] = 1
                weights['mutrigger']['up'] = 1
                weights['mutrigger']['down'] = 1
                weights['muid']['nom'] = 1
                weights['muid']['up'] = 1
                weights['muid']['down'] = 1
                weights['muiso']['nom'] = 1
                weights['muiso']['up'] = 1
                weights['muiso']['down'] = 1

        #########################################
        # Fill control histograms #
        #########################################
        nominalweight = Weightify(weights,'nominal',['Pu','mutrigger','muid','muiso'])
        if not isData:
            if ak8JetsColl[0].msoftdrop > 20 and ak8JetsColl[0].pt > 200 and not Hbbsel['leptonExists']:
                N2_map.Fill(2*log(ak8JetsColl[0].msoftdrop/ak8JetsColl[0].pt),ak8JetsColl[0].pt,ak8JetsColl[0].n2b1,norm_weight*nominalweight)
        if preselection:
            Hbb_rho.Fill(Hbbrho,norm_weight*nominalweight)
            Hbb_pT.Fill(jet0.Pt(),norm_weight*nominalweight)
            Hbb_mSD_uncorr.Fill(ak8JetsColl[0].msoftdrop,norm_weight*nominalweight)
            Hbb_mSD.Fill(jet0.M(),norm_weight*nominalweight)
            Hbb_eta.Fill(jet0.Eta(),norm_weight*nominalweight)
            MET.Fill(inTree.readBranch('MET_pt'),norm_weight*nominalweight)
            Hbb_doubleB.Fill(getattr(ak8JetsColl[0],doubleB_name),norm_weight*nominalweight)
            Hbb_deepAK8.Fill(getattr(ak8JetsColl[0],'deepTagMD_ZHbbvsQCD'),norm_weight*nominalweight)
            Hbb_n2.Fill(ak8JetsColl[0].n2b1,norm_weight*nominalweight)
            Hbb_n2ddt_new.Fill(n2ddt_new,norm_weight*nominalweight)
            Muon_pT.Fill(muon_pt,norm_weight*nominalweight)
            Muon_eta.Fill(muon_eta,norm_weight*nominalweight)
            Muon_phi.Fill(muon_phi,norm_weight*nominalweight)
            Muon_dphi.Fill(muon_dphi,norm_weight*nominalweight)
            Muon_pfRelIso04_all.Fill(muon_pfRelIso04_all,norm_weight*nominalweight)
            if ak8JetsColl[0].pt > 450:
                Hbb_pT_450.Fill(jet0.Pt(),norm_weight*nominalweight)
                Hbb_mSD_450.Fill(jet0.M(),norm_weight*nominalweight)
            if Hbbsel['n2ddt_new']:
                if Hbbsel['pass_M']:
                    Hbb_mSD_passM_n2ddt.Fill(jet0.M(),norm_weight*nominalweight)
                else:
                    Hbb_mSD_failM_n2ddt.Fill(jet0.M(),norm_weight*nominalweight)

        ######################################### 
        # Check preselection #
        ######################################### 

        if(getattr(ak8JetsColl[0],doubleB_name)<0.2):
            continue

        if preselection:
            #print('filling')
            for wp in WPs:
                for cat in ['pass','fail']:
                    if Hbbsel['%s_%s'%(cat,wp)]:
                        hists[wp][cat].Fill(jet0.M(),norm_weight*nominalweight)

            # Fill variations
            if runOthers and not isData:
                for wp in WPs:
                    for cat in ['pass','fail']:
                        for sys in ['matched','unmatched','semimatched']:
                            if Hbbsel['%s_%s'%(cat,wp)] and Hbbsel[sys]:
                                mass_corr = jet0.M()
                                hists[wp]['%s_%s'%(cat,sys)].Fill(mass_corr,norm_weight*nominalweight)
                        for sys in ['Pu','mutrigger','muiso','muid']:
                            if Hbbsel['%s_%s'%(cat,wp)]:
                                hists[wp]['%s_%sUp'%(cat,sys)].Fill(jet0.M(),norm_weight*Weightify(weights,'%s_up'%sys,['Pu','mutrigger','muid','muiso']))
                                hists[wp]['%s_%sDown'%(cat,sys)].Fill(jet0.M(),norm_weight*Weightify(weights,'%s_down'%sys,['Pu','mutrigger','muid','muiso']))

    end = time.time()
    print '\n'
    print str((end-start)/60.) + ' min'
    print(nevents)
    #print('nevts',nevents/(highBinEdge-lowBinEdge))
    f.cd()
    f.Write()
    f.Close()
