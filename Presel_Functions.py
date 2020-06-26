##################################################################
##                                                              ##
## Name: Bstar_Functions.py                                     ##
## Author: Kevin Nash                                           ##
## Date: 5/13/2015                                              ##
## Purpose: This contains all functions used by the             ##
##      analysis.  A method is generally placed here if         ##
##      it is called more than once in reproducing all          ##
##      analysis results.  The functions contained here         ##
##      Are capable of tuning the analysis - such as changing   ##
##      cross sections, updating lumi, changing file            ##
##      locations, etc. with all changes propegating            ##
##      to all relevant files automatically.                    ##
##                                                              ##
##################################################################

import os
import array
import glob
import math
import itertools
from random import random
from math import sqrt, exp, log
import ROOT
import sys
import time
import subprocess
import cppyy
import pickle
from array import *
from ROOT import *
from PhysicsTools.NanoAODTools.postprocessing.tools import *

#This is the most impostant Function.  Correct information here is essential to obtaining valid results.
#In order we have Luminosity, top tagging scale factor, cross sections for wprime right,left,mixed,ttbar,qcd, and singletop and their corresponding event numbers
#If I wanted to access the left handed W' cross section at 1900 GeV I could do Xsecl1900 = LoadConstants()['xsec_wpl']['1900']
#FYI Lumi and xsecs are in femtobarns here
def LoadConstants(year):
    constants = {
        'QCDHT700_xsec':6802000,
        'QCDHT1000_xsec':1206000,
        'QCDHT1500_xsec':120400,
        'QCDHT2000_xsec':25250,
        'RadNar-1000_xsec':1,
        'RadNar-1500_xsec':1,
        'RadNar-2000_xsec':1,
        'RadNar-2500_xsec':1,
        'RadNar-3000_xsec':1,
        'GravNar-1000_xsec':1,
        'GravNar-1500_xsec':1,
        'GravNar-2000_xsec':1,
        'GravNar-2500_xsec':1,
        'GravNar-3000_xsec':1
    }
    if year == '16':
        constants['ttbar_xsec'] = 831760
	constants['Wqq_xsec'] = 2788000
	constants['Zqq_xsec'] = 1187000
        constants['lumi'] = 35.872301001
        constants['GJetsHT100_xsec']= 9238000
        constants['GJetsHT200_xsec']= 2305000
	constants['GJetsHT400_xsec']= 275200
        constants['GJetsHT600_xsec']= 93190
	constants['ZGToJJG_xsec'] = 570
        constants['WGToJJG_xsec'] = 1220
        constants['singletop_xsec'] = 136.02 # t?
        constants['singleantitop_xsec'] = 81.42
        constants['singletWtop_xsec'] = 35.85
        constants['singletWantitop_xsec'] = 35.85
        constants['ZZdiboson_xsec'] = 16.523
        constants['WWdiboson_xsec'] = 63.21*1.878
        constants['wjetstolnu_xsec'] = 50380.0
        constants['dyjetstoll1jet_xsec'] = 859.589402
        constants['shift_SF'] = 1.001 # n2ddt 26%
        #constants['shift_SF'] = 1.014 # tau21ddt 0.43

    elif year == '17':
        constants['lumi'] = 41.518865298 #35851.0,
        constants['ttbar_xsec'] = 377960 #uncertainty +4.8%-6.1%
        constants['ttbar-semilep_xsec'] = 365340
        constants['WjetsHT400_xsec'] = 315600
        constants['WjetsHT600_xsec'] = 68570
        constants['WjetsHT800_xsec'] = 34900
        constants['ZjetsHT400_xsec'] = 145400
        constants['ZjetsHT600_xsec'] = 68570
        constants['ZjetsHT800_xsec'] = 34900
        constants['shift_SF'] = 1

    elif year == '18':
        constants['lumi'] = 59.660725495
        constants['ttbar_xsec'] = 377960 #uncertainty +4.8%-6.1%
        constants['ttbar-semilep_xsec'] = 365340
        constants['WjetsHT400_xsec'] = 315600
        constants['WjetsHT600_xsec'] = 68570
        constants['WjetsHT800_xsec'] = 34900
        constants['ZjetsHT400_xsec'] = 145400
        constants['ZjetsHT600_xsec'] = 68570
        constants['ZjetsHT800_xsec'] = 34900
        constants['shift_SF'] = 1

    return constants
    
def LoadCuts(region,year):
    cuts = {
        'hpt':[300.0,float("inf")],
        'bpt':[30.0,float("inf")],
        'hmass':[105.0,135.0],
        'bbmass':[90.,140.],
	'btag':[0.6321,1.0],
        'deepbtag':[0.8953,1.0], #medium 2016
        'deepbtag_medium':[0.6321,1.0], 
        'eta':[0.0,2.4],
        'tau21ddt':[0,0.43], # loose WP
    }

    # double-check these values
    if region=='deepTagMD_ZHbbvsQCD':
        cuts['doublebtag'] = [0.9,1.0]
        cuts['doublebtagTight'] = [0.95,1.0]
        cuts['doublebtagLoose'] = [0.8,1.0]

    if region=='btagDDBvL':
        cuts['doublebtag'] = [0.89,1.0] #medium
        cuts['doublebtagTight'] = [0.9,1.0]
        cuts['doublebtagLoose'] = [0.8,1.0] # use this for loose W/Z templates

    return cuts

def CorrectMSD(jet,subJets,puppisd_corrGEN,puppisd_corrRECO_cen,puppisd_corrRECO_for):
    jets_msoftdrop_raw = 0.0
    jets_msoftdrop_nom = 0.0
    if jet.subJetIdx1 >= 0 and jet.subJetIdx2 >= 0 :
        j1 = subJets[ jet.subJetIdx1 ].p4()
        j1_uncorr = j1
        #j1_uncorr.SetPtEtaPhiM(j1.Pt()*(1-jet.rawFactor), j1.Eta(), j1.Phi(), j1.M()*(1-jet.rawFactor))
        j1_uncorr.SetPt(j1.Pt()*(1-jet.rawFactor))
        j1_uncorr.SetEta(j1.Eta())
        j1_uncorr.SetPhi(j1.Phi())
        j1_uncorr.SetM(j1.M()*(1-jet.rawFactor))
        j2 = subJets[ jet.subJetIdx2 ].p4()
        j2_uncorr = j2
        #j2_uncorr.SetPtEtaPhiM(j2.Pt()*(1-jet.rawFactor), j2.Eta(), j2.Phi(), j2.M()*(1-jet.rawFactor))
        j2_uncorr.SetPt(j2.Pt()*(1-jet.rawFactor))
        j2_uncorr.SetEta(j2.Eta())
        j2_uncorr.SetPhi(j2.Phi())
        j2_uncorr.SetM(j2.M()*(1-jet.rawFactor))        
        groomedP4 = j1_uncorr+ j2_uncorr
    else :
        groomedP4 = None

    # CMS: uncorrect groomed subjets
    if groomedP4 != None:
        #groomedP4.SetPtEtaPhiM(groomedP4.Pt(), groomedP4.Eta(), groomedP4.Phi(), groomedP4.M())
        groomedP4.SetPt(groomedP4.Pt())
        groomedP4.SetPhi(groomedP4.Phi())
        groomedP4.SetEta(groomedP4.Eta())
        groomedP4.SetM(groomedP4.M())        
        jets_msoftdrop_raw = groomedP4.M()
    
    # LC: Apply PUPPI SD mass correction https://github.com/cms-jet/PuppiSoftdropMassCorr/
    puppisd_genCorr = puppisd_corrGEN.Eval(jet.pt)
    if abs(jet.eta) <= 1.3:
        puppisd_recoCorr = puppisd_corrRECO_cen.Eval(jet.pt)
    else:
        puppisd_recoCorr = puppisd_corrRECO_for.Eval(jet.pt)
        
    puppisd_total = puppisd_genCorr * puppisd_recoCorr
    
    if groomedP4 != None:
       #groomedP4.SetPtEtaPhiM(groomedP4.Perp(), groomedP4.Eta(), groomedP4.Phi(), groomedP4.M()*puppisd_total)
        groomedP4.SetPt(groomedP4.Perp())
        groomedP4.SetPhi(groomedP4.Phi())
        groomedP4.SetEta(groomedP4.Eta())
        groomedP4.SetM(groomedP4.M()*puppisd_total)
        jets_msoftdrop_nom = groomedP4.M()

    jets_groomed_corr_PUPPI = puppisd_total
    return jets_msoftdrop_nom

#This function loads up Ntuples based on what type of set you want to analyze.  
#This needs to be updated whenever new Ntuples are produced (unless the file locations are the same).
def Load_jetNano(string,year):
    print 'running on ' + string 
    #return 'root://cmseos.fnal.gov//store/user/dbrehm/data18andTTbarSignalMC/rootfiles/'+string+'_hh'+year+'.root'
    return 'root://cmseos.fnal.gov//store/group/lpctlbsm/dbrehm/HH4bv5/'+string+'_hh'+year+'.root'

def PDF_Lookup(pdfs , pdfOP ):
    # Computes the variance of the pdf weights to estimate the up and down uncertainty for each set (event)
    ilimweight = 0.0

    limitedpdf = []
    for ipdf in range(pdfs.GetSize()):
        curpdf = pdfs[ipdf]
        if abs(curpdf)<1000.0:
            limitedpdf.append(curpdf)

    if len(limitedpdf) == 0:
        return 1

    limave =  limitedpdf
    limave =  reduce(lambda x, y: x + y, limitedpdf) / len(limitedpdf)
    #print ave
    for limpdf in limitedpdf :
        ilimweight = ilimweight + (limpdf-limave)*(limpdf-limave)

    if pdfOP == "up" :
        return min(13.0,1.0+sqrt((ilimweight) / (len(limitedpdf))))
    else :
        return max(-12.0,1.0-sqrt((ilimweight) / (len(limitedpdf))))

def Trigger_Lookup( massForTrig, ptForTrig , TRP ):
    Weight = TRP.GetEfficiency(TRP.FindFixBin(massForTrig, ptForTrig))
    WeightUp = Weight + TRP.GetEfficiencyErrorUp(TRP.FindFixBin(massForTrig, ptForTrig))
    WeightDown = Weight - TRP.GetEfficiencyErrorLow(TRP.FindFixBin(massForTrig, ptForTrig))
    if Weight <= 0 or WeightDown <= 0 or WeightUp <= 0:
        Weight = 1.0
        WeightUp = 1.0
        WeightDown = 1.0
    return [Weight,WeightUp,WeightDown]

def Hist2d_Lookup(var1, var2, Hist):
    Weight = Hist.GetBinContent(Hist.FindBin(var1,var2))
    WeightUp = Weight + Hist.GetBinError(Hist.FindBin(var1,var2))
    WeightDown = Weight - Hist.GetBinError(Hist.FindBin(var1,var2))
    if Weight <= 0 or WeightDown <= 0 or WeightUp <= 0:
        Weight = 1.0
        WeightUp = 1.0
        WeightDown = 1.0
    return [Weight,WeightUp,WeightDown]

def kFactor_Lookup( genpt, histsk, name):
    vjetsKF = 1.
    if 'zqq_16' in name: #for 2016legacy 
        ptForNLO = max(250., min(genpt, 1000.))
        iEWKKF = histsk['hEWK_Z'].GetBinContent(histsk['hEWK_Z'].FindBin(ptForNLO));
        iQCDKF = histsk['znlo'].GetBinContent(histsk['znlo'].FindBin(ptForNLO))        # New QCD KF for 2017                                                                                                  
        vjetsKF = iQCDKF*iEWKKF;
    elif 'zqq_' in name:
        ptForNLO = max(250., min(genpt, 1000.))
        iEWKKF = histsk['hEWK_Z'].GetBinContent(histsk['hEWK_Z'].FindBin(ptForNLO));  # same EWK as 2016                                                                                                          
        iQCDKF = histsk['znlo'].GetBinContent(histsk['znlo'].FindBin(ptForNLO))
        vjetsKF = iQCDKF*iEWKKF;
    if 'wqq_16' in name: #for 2016legacy 
        ptForNLO = max(250., min(genpt, 1000.))
        iQCDKF = histsk['wnlo'].GetBinContent(histsk['wnlo'].FindBin(ptForNLO))        # New QCD KF for bacon 13+                                                                                               
        iEWKKF  = histsk['hEWK_W'].GetBinContent(histsk['hEWK_W'].FindBin(ptForNLO));
        vjetsKF = iQCDKF*iEWKKF;
    elif 'wqq_' in name:
        ptForNLO = max(250., min(genpt, 1000.))
        iEWKKF = histsk['hEWK_W'].GetBinContent(histsk['hEWK_W'].FindBin(ptForNLO));  # same EWK as 2016                                                                                                         
        iQCDKF = histsk['wnlo'].GetBinContent(histsk['wnlo'].FindBin(ptForNLO))
        vjetsKF = iQCDKF*iEWKKF;
    return [vjetsKF]

class myGenParticle:
    def __init__ (self, index, genpart):
        self.idx = index
        self.genpart = genpart
        self.status = genpart.status
        self.pdgId = genpart.pdgId
        self.vect = TLorentzVector()
        #self.vect.SetPtEtaPhiM(genpart.pt,genpart.eta,genpart.phi,genpart.mass)
        self.vect.SetPt(genpart.pt)
        self.vect.SetEta(genpart.eta)
        self.vect.SetPhi(genpart.phi)
        self.vect.SetM(genpart.mass)
        self.motherIdx = genpart.genPartIdxMother

# Compares ak4 jets against leading ak8 and looks for any in opposite hemisphere 
def Hemispherize(fatjetCollection,jetCollection):
    # First find the highest pt ak8 jet with mass > 40 geV                                                                                                                      
    candidateFatJetIndex = -1
    for fjet in range(0,len(fatjetCollection)):
        if 40 < fatjetCollection[fjet].msoftdrop < 220:
            candidateFatJetIndex = fjet
            break
            
    if candidateFatJetIndex == -1:
        return []
        
    leadFatJet = fatjetCollection[candidateFatJetIndex]

    # Maintain same indexing for these throughout next bit                                                                                                                    
    candidateJetIndices = []
    
    # Check the AK4s against the AK8                                                                                                                                              
    for ijet in range(0,len(jetCollection)):
        if abs(deltaPhi(leadFatJet.phi,jetCollection[ijet].phi))>TMath.Pi()/2.0:
            candidateJetIndices.append(ijet)

    if len(candidateJetIndices) > 0:
        return [jetCollection[j] for j in candidateJetIndices]
    else:
        return []

def Weightify(wd,outname,corrections=['Pu','trigger','kFactor']):
    final_w = 1.0

    if outname == 'nominal':
        for c in corrections:
            if 'nom' in wd[c].keys():
                final_w = final_w*wd[c]['nom']

    elif outname.split('_')[0] in corrections:
        if outname.split('_')[1] in wd[outname.split('_')[0]].keys():
            final_w = wd[outname.split('_')[0]][outname.split('_')[1]]
        for c in corrections:
            if (c != outname.split('_')[0]) and ('nom' in wd[c].keys()):
                final_w = final_w*wd[c]['nom']

    else:
        final_w = wd[outname.split('_')[0]][outname.split('_')[1]]
        for c in corrections:
            if 'nom' in wd[c].keys():
                final_w = final_w*wd[c]['nom']
    
    return final_w

#This is just a quick function to automatically make a tree
#This is used right now to automatically output branches used to validate the cuts used in a run
def Make_Trees(Floats,name="Tree"):
    t = TTree(name, name);
    print "Booking trees"
    for F in Floats.keys():
        t.Branch(F, Floats[F], F+"/D")
    return t

#Makes the blue pull plots
def Make_Pull_plot( DATA,BKG,BKGUP,BKGDOWN ):
    pull = DATA.Clone("pull")
    pull.Add(BKG,-1)
    sigma = 0.0
    FScont = 0.0
    BKGcont = 0.0
    for ibin in range(1,pull.GetNbinsX()+1):
        FScont = DATA.GetBinContent(ibin)
        BKGcont = BKG.GetBinContent(ibin)
        if FScont>=BKGcont:
            FSerr = DATA.GetBinErrorLow(ibin)
            BKGerr = abs(BKGUP.GetBinContent(ibin)-BKG.GetBinContent(ibin))
        if FScont<BKGcont:
            FSerr = DATA.GetBinErrorUp(ibin)
            BKGerr = abs(BKGDOWN.GetBinContent(ibin)-BKG.GetBinContent(ibin))
        sigma = sqrt(FSerr*FSerr + BKGerr*BKGerr)
        if FScont == 0.0:
            pull.SetBinContent(ibin, 0.0 )  
        else:
            if sigma != 0 :
                pullcont = (pull.GetBinContent(ibin))/sigma
                pull.SetBinContent(ibin, pullcont)
            else :
                pull.SetBinContent(ibin, 0.0 )
    return pull

# v jet matching  by looking up the deltaR separation 
# of the daughter particle from the jet axis. If passes, return 1.
# this is needed for Ws/Zs and Hs
def VJetMatching(vjetVect,genparticles,vid=24,iMass=80.4):
    import GenParticleChecker
    from GenParticleChecker import GenParticleTree,GenParticleObj
    
    # Build the tree of gen particles that we care about
    particle_tree = GenParticleTree()
    Vs = []
    quarks = []
    prongs = [] # Final particles we'll check
    for i,p in enumerate(genparticles):
        # Internal class info
        this_gen_part = GenParticleObj(i,p)
        this_gen_part.SetStatusFlags()
        this_gen_part.SetPDGName(abs(this_gen_part.pdgId))
        
        # Add particles to tree and keep track of them in external lists
        if abs(this_gen_part.pdgId) == vid:# and this_gen_part.status == 22: # 22 means intermediate part of hardest subprocess, only other to appear is 52 (outgoing copy of recoiler, with changed momentum)
            particle_tree.AddParticle(this_gen_part)
            Vs.append(this_gen_part)

        elif abs(this_gen_part.pdgId) >= 1 and abs(this_gen_part.pdgId) <= 5:
            particle_tree.AddParticle(this_gen_part)
            quarks.append(this_gen_part)

    for q in quarks:
        # if parent is a v and inside jet cone 
        if particle_tree.GetParent(q) and abs(particle_tree.GetParent(q).pdgId) == vid and particle_tree.GetParent(q).vect.DeltaR(vjetVect) < 0.8 and q.vect.DeltaR(vjetVect) < 0.8:
            #if particle_tree.GetParent(q) and abs(particle_tree.GetParent(q).pdgId) == vid and particle_tree.GetParent(q).vect.DeltaR(vjetVect) < 0.6 and q.vect.DeltaR(vjetVect) < 0.6:
            prongs.append(q)
    
    matched = False
    if len(Vs)>0:
        genVPt = Vs[0].pt # taking leading as genVPt
        if len(prongs) == 2: 
            pPhi = math.fabs(Vs[0].phi - vjetVect.Phi())
            pPt = math.fabs(Vs[0].pt - vjetVect.Pt())/Vs[0].pt
            pMass = math.fabs(iMass - vjetVect.M())/iMass
            if pPhi < 0.8 and pPt < 0.5 and pMass < 0.3:
                matched = True
    else:
        genVPt = -1

    return [matched,genVPt]

def getN2dtt(trans_h2ddt,n2,rho,pt):
    cur_rho_index = trans_h2ddt.GetXaxis().FindBin(rho)
    cur_pt_index = trans_h2ddt.GetYaxis().FindBin(pt)
    if rho > trans_h2ddt.GetXaxis().GetBinUpEdge(trans_h2ddt.GetXaxis().GetNbins()): cur_rho_index = trans_h2ddt.GetXaxis().GetNbins()
    if rho < trans_h2ddt.GetXaxis().GetBinLowEdge(1): cur_rho_index = 1
    if pt > trans_h2ddt.GetYaxis().GetBinUpEdge(trans_h2ddt.GetYaxis().GetNbins()): cur_pt_index = trans_h2ddt.GetYaxis().GetNbins()
    if pt < trans_h2ddt.GetYaxis().GetBinLowEdge(1): cur_pt_index = 1
    n2ddt = n2 - trans_h2ddt.GetBinContent(cur_rho_index, cur_pt_index)
    return n2ddt

def TopJetMatching(jet, genparticles):
    from GenParticleChecker import GenParticleTree,GenParticleObj
    # Build the tree of gen particles that we care about
    particle_tree = GenParticleTree()
    tops = []
    Ws = []
    quarks = []
    prongs = [] # Final particles we'll check
    for i,p in enumerate(genparticles):
        # Internal class info
        this_gen_part = GenParticleObj(i,p)
        this_gen_part.SetStatusFlags()
        this_gen_part.SetPDGName(abs(this_gen_part.pdgId))
        
        # Add particles to tree and keep track of them in external lists
        if abs(this_gen_part.pdgId) == 6 and this_gen_part.DeltaR(jet)<0.8:# and this_gen_part.status == 62: # 22 means intermediate part of hardest subprocess, only other to appear is 62 (outgoing subprocess particle with primordial kT included)
            particle_tree.AddParticle(this_gen_part)
            tops.append(this_gen_part)

        elif abs(this_gen_part.pdgId) == 24:# and this_gen_part.status == 22: # 22 means intermediate part of hardest subprocess, only other to appear is 52 (outgoing copy of recoiler, with changed momentum)
            particle_tree.AddParticle(this_gen_part)
            Ws.append(this_gen_part)

        elif abs(this_gen_part.pdgId) >= 1 and abs(this_gen_part.pdgId) <= 5 and this_gen_part.status == 23:
            particle_tree.AddParticle(this_gen_part)
            quarks.append(this_gen_part)

        elif this_gen_part.DeltaR(jet)<0.8:
            particle_tree.AddParticle(this_gen_part)

    for W in Ws:
        # If parent is a top that matches with the jet
        wParent = particle_tree.GetParent(W)
        if wParent != False:
            if abs(wParent.pdgId) == 6 and wParent.DeltaR(jet) < 0.8:
                # Loop over the daughters of the W
                this_W = W
                # Skip down chain of W's to last one
                if len(particle_tree.GetChildren(this_W)) == 1 and particle_tree.GetChildren(this_W)[0].pdgId == W.pdgId:
                    this_W = particle_tree.GetChildren(this_W)[0]
                
                for c in particle_tree.GetChildren(this_W):
                    if abs(c.pdgId) >= 1 and abs(c.pdgId) <= 5:
                        # Append daughter quarks to prongs
                        prongs.append(c)

    for q in quarks:
        # if bottom      and     has a parent              and   parent is a top                          and    the top matches to the jet
        if abs(q.pdgId) == 5:
            bottomParent = particle_tree.GetParent(q)
            # if parent exists
            if bottomParent != False:
                # if parent is a top matched to the jet
                if abs(bottomParent.pdgId) == 6 and bottomParent.DeltaR(jet) < 0.8:
                    prongs.append(q)

    # Now that we have the prongs, check how many are merged
    merged_particles = 0
    if len(prongs) < 3: # you've either tagged a QCD jet and can't match to a gen top or the W decayed leptonically so don't apply a SF in either case
        # if len(prongs) == 1: 
        #     particle_tree.PrintTree(ievent,['idx','status'],wp,jet)
        #     raw_input(prongs[0].idx)
        return [1,1,1],-len(prongs)-1


    for p in prongs:
        if p.DeltaR(jet) < 0.8:
            merged_particles += 1

    #if merged_particles == 3:
    #    print('mergedTop')
    #elif merged_particles == 2:
    #    print('semimerged')
    #elif merged_particles == 1:
    #    print('nomerged')
    #else:
    #return [1,1,1],-len(prongs)-1

    del particle_tree

    return merged_particles
