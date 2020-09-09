import os
import glob
import math
import array
import sys
import time
from optparse import OptionParser
import ROOT
from ROOT import *
from math import sqrt

def makeDDTmap(iH,iWP=0.26):
    lDDT = iH.Project3D("yx")
    for i0 in range(iH.GetXaxis().GetNbins()):
        for i1 in range(iH.GetYaxis().GetNbins()):
            pProj = iH.ProjectionZ(iH.GetName()+str(i0)+str(i1)+"ddt",i0,i0,i1,i1)
            if pProj.Integral() == 0: continue
            lp = array.array('d', [iWP])
            lq = array.array('d', [0.0]*len(lp))
            pProj.GetQuantiles(len(lp), lq, lp)
            lDDT.SetBinContent( i0, i1, lq[0] )
    return lDDT

def getRatio(hist, reference):
    ratio = hist.Clone("%s_ratio"%hist.GetName())
    ratio.SetDirectory(0)
    ratio.SetLineColor(hist.GetLineColor())
    ratio.SetMarkerSize(hist.GetMarkerSize())
    for xbin in xrange(1,reference.GetNbinsX()+1):
        ref = reference.GetBinContent(xbin)
        val = hist.GetBinContent(xbin)
        refE = reference.GetBinError(xbin)
        valE = hist.GetBinError(xbin)
        try:
            ratio.SetBinContent(xbin, val/ref)
            ratio.SetBinError(xbin, math.sqrt( (val*refE/(ref**2))**2 + (valE/ref)**2 ))
        except ZeroDivisionError:
            ratio.SetBinError(xbin, 0.0)
    return ratio

def TOTerror(hmc, ratio ):
    hmc.Sumw2()
    den1 = hmc.Clone ("den1");
    nvar = hmc.GetNbinsX();
    for km in range(1,nvar+1):
        delta = hmc.GetBinError(km)
        den1.SetBinError(km,0)
    ratiop = hmc.Clone("ratiop");
    ratiop.Divide(den1);
    return ratiop;

def makeCanvasComparisonStackWData(hd,hb,legname,color,style,outname,lumi=30,normalize=None,ratio=True):
    hd.SetMarkerColor(1)
    hd.SetMarkerSize(10)
    hd.SetLineStyle(1)
    hd.SetLineColor(1)

    for name, h in sorted(hb.iteritems(),key=lambda (k,v): v.Integral()):
        error = array.array('d',[0.0])
        integral = h.IntegralAndError(1,h.GetNbinsX(),error)
        print(name,integral)
    error = array.array('d',[0.0])
    integral = hd.IntegralAndError(1,hd.GetNbinsX(),error)
    dataInt = integral
    dataErr = error[0]

    maxval = -99
    hstack = ROOT.THStack("hstack","hstack")
    for name, h in sorted(hb.iteritems(),key=lambda (k,v): v.Integral()):
        hstack.Add(h)
        h.SetFillColor(color[name])
        h.SetLineColor(1)
        h.SetLineStyle(1)
        h.SetLineWidth(1)
        h.SetFillStyle(1001)
        if h.GetMaximum() > maxval: maxval = h.GetMaximum()
    allMC=hstack.GetStack().Last().Clone()
    maxval = max(hd.GetMaximum(),maxval)
    fullmc = hstack.GetStack().Last();
    scalefactor = hd.Integral()/fullmc.Integral();
    print "data/mc scale factor = ", scalefactor
    # normalizing 
    if normalize is not None:
        for name, h in sorted(hb.iteritems(),key=lambda (k,v): v.Integral()):
            if 'qcd' in name and normalize=='qcd':
                print('scaling ',h,'by ',scalefactor)
                h.Scale( scalefactor );
            if 'tqq' in name and normalize=='tqq':
                print('scaling ',h,'by ',scalefactor)
                h.Scale( scalefactor );

    hstack2 = ROOT.THStack("hstack2","hstack2");
    tempFile = ROOT.TFile("test.root","UPDATE")
    tempFile.cd()
    hd.Write()
    for name, h in sorted(hb.iteritems(),key=lambda (k,v): v.Integral()):
        hstack2.Add(h);
        h.Write()
        h.SetFillColor(color[name])
        h.SetLineColor(1)
        h.SetLineStyle(1)
        h.SetLineWidth(1)
        h.SetFillStyle(1001)
    tempFile.Close()
    leg_y = 0.88 - (2+int(len(hb)/3))*0.03
    leg = ROOT.TLegend(0.2,leg_y,0.5,0.88) 
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.035)
    leg.SetTextFont(42)
    leg2 = ROOT.TLegend(0.5,leg_y,0.78,0.88,)
    leg2.SetFillStyle(0)
    leg2.SetBorderSize(0)
    leg2.SetTextSize(0.035)
    leg2.SetTextFont(42)
    leg3 = ROOT.TLegend(0.65,leg_y,0.98,0.88,)
    leg3.SetFillStyle(0)
    leg3.SetBorderSize(0)
    leg3.SetTextSize(0.035)
    leg3.SetTextFont(42)
    count=1
    for name, h in sorted(hb.iteritems(),key=lambda (k,v): -v.Integral()):
        if count <4:
            #if name in 'qcd' : leg.AddEntry(h,legname[name]+" (k-factor %.2f)"%scalefactor,"f")
            leg.AddEntry(h,legname[name],"f")
            #else : leg.AddEntry(h,legname[name],"f")
        elif count >3 and count<7 : leg2.AddEntry(h,legname[name],"f")
        elif count >6 : leg3.AddEntry(h,legname[name],"f")
        count = count+1
    leg3.AddEntry(hd,'Data',"pe");

    c = ROOT.TCanvas("c"+outname,"c"+outname,1000,800)
    c.SetFillStyle(4000)
    c.SetFrameFillStyle(1000)
    c.SetFrameFillColor(0)
    if ratio:
        oben = ROOT.TPad('oben','oben',0,0.3 ,1.0,1.0)
        unten = ROOT.TPad('unten','unten',0,0.0,1.0,0.3)
        oben.SetBottomMargin(0)
        unten.SetTopMargin(0.)
        unten.SetBottomMargin(0.35)
    else:
        oben = ROOT.TPad('oben','oben',0,0.03 ,1.0,1.0)
        unten = ROOT.TPad('unten','unten',0,0.0,1.0,0.0)
    oben.SetFillStyle(4000)
    oben.SetFrameFillStyle(1000)
    oben.SetFrameFillColor(0)
    unten.SetFillStyle(4000)
    unten.SetFrameFillStyle(1000)
    unten.SetFrameFillColor(0)
    oben.Draw()
    unten.Draw()
    oben.cd()

    hd.Draw('PE')
    hstack2.Draw('hist same')
    hstack2.SetMinimum(1.)
    hstack2.SetMaximum(1.2*maxval)
    hstack2.GetYaxis().SetRangeUser(1.,1.2*maxval)
    hstack2.GetYaxis().SetTitle('Events')
    hstack2.GetYaxis().SetTitleOffset(1.0)
    hstack2.GetXaxis().SetTitle(allMC.GetXaxis().GetTitle())
    hstack2.GetXaxis().SetTitleOffset(1.3)
    hstack2.GetXaxis().SetLabelSize(0.04)
    hstack2.GetXaxis().SetTitleSize(0.045)
    hstack2.Draw('hist')
    leg.Draw()
    leg2.Draw()
    leg3.Draw()
    hstack2.SetMinimum(1)
    hd.Draw('PE same');
    allMC2=hstack2.GetStack().Last().Clone()
    '''
    for name, h in sorted(hb.iteritems(),key=lambda (k,v): -v.Integral()):
        if name in 'qcd' :
        #if name in 'tqq' :
            herr = h.Clone('herr')
            herr2 = h.Clone('herr2')
    for name, h in sorted(hb.iteritems(),key=lambda (k,v): -v.Integral()):
        #if name in 'qcd' : continue
        if name in 'tqq' : continue
        for ibin in range(1,h.GetNbinsX()+1):
           valA  = herr.GetBinContent(ibin);
           evalA = herr.GetBinError(ibin);
           valB  = h.GetBinContent(ibin);
           evalB = h.GetBinError(ibin);
           herr.SetBinContent(ibin,(valA+valB));
           herr.SetBinError(ibin,sqrt(evalA*evalA+evalB*evalB));
           if(valA+valB >0): herr2.SetBinContent(ibin,(valA+valB+sqrt(evalA*evalA+evalB*evalB))/(valA+valB));
           else : herr2.SetBinContent(ibin,1);

    theErrorGraph = ROOT.TGraphErrors(herr)
    theErrorGraph.SetFillColor(ROOT.kGray+2)
    theErrorGraph.SetFillStyle(3002)
    herr.SetFillColor(ROOT.kGray+2)
    herr.SetFillStyle(3002)
    herr.SetMarkerColor(1111);
    leg3.AddEntry(herr,"MC uncert. (stat.)","fl")
    hd.Draw('PEsame');
    theErrorGraph.Draw('SAME E')
    '''

    tag1 = ROOT.TLatex(0.67,0.92,"%.1f fb^{-1} (13 TeV)"%lumi)
    tag1.SetNDC(); tag1.SetTextFont(42)
    tag1.SetTextSize(0.045)
    tag2 = ROOT.TLatex(0.17,0.92,"CMS")
    tag2.SetNDC()
    tag2.SetTextFont(62)
    tag3 = ROOT.TLatex(0.27,0.92,"Preliminary")
    tag3.SetNDC()
    tag3.SetTextFont(52)
    tag2.SetTextSize(0.055)
    tag1.Draw()
    tag2.Draw()
    tag3.Draw()

    if ratio:
        unten.cd()
        ratio = getRatio(hd,allMC2)
        herr3= TOTerror(allMC2,ratio);
        toterree = ROOT.TGraphErrors(herr3)
        ksScore = hd.KolmogorovTest( allMC2 )
        chiScore = hd.Chi2Test( allMC2 , "UWCHI2/NDF")
        ratio.SetStats(0)
        ratio.GetYaxis().SetRangeUser(0,2)
        ratio.GetYaxis().SetNdivisions(504)
        ratio.GetYaxis().SetTitle("Data/Simulation")
        ratio.GetXaxis().SetTitle(allMC.GetXaxis().GetTitle())
        ratio.GetXaxis().SetTitleSize(0.14)
        ratio.GetXaxis().SetTitleOffset(1.0)
        ratio.GetYaxis().SetTitleOffset(0.5)
        ratio.GetYaxis().SetLabelSize(0.12)
        ratio.GetYaxis().SetTitleSize(0.11)
        ratio.GetXaxis().SetLabelSize(0.11)

        line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1.0,
                          ratio.GetXaxis().GetXmax(), 1.0)
        line.SetLineColor(ROOT.kGray)
        line.SetLineStyle(2)
        line.Draw()
        tKsChi = ROOT.TLatex()
        tKsChi.SetNDC()
        tKsChi.SetTextFont(42)
        tKsChi.SetTextSize(0.09)
        ratio.Draw("P E ")

        toterree.SetFillColor(ROOT.kGray+2);
        toterree.SetLineColor(ROOT.kGray+2);
        toterree.SetFillStyle(3002);
        toterree.Draw("2 same");
        line.Draw("same")

        leg4 = ROOT.TLegend(0.7,0.89,0.5,0.8)
        leg4.SetFillStyle(0)
        leg4.SetBorderSize(0)
        leg4.SetTextSize(0.05)
        leg4.SetTextFont(42)
        leg4.AddEntry(toterree,"MC uncert. (stat.)","fl")
        leg4.Draw()

    c.SaveAs(outname+".png")
    oben.SetLogy()
    c.SaveAs(outname+"_log.png")
    return c

def main(options,args):
    variations = []
    if options.muonCR:
        plots = [
            #'h_Lpass',
            'h_Mpass',
            #'h_Tpass',
            #'h_Lfail',
            'h_Mfail',
            #'h_Tfail',
            #'h_cutflow',
            'h_doubleB',
            'h_deepAK8',
            #'h_pT',
            #'h_mSD',
            #'h_muonpT',
            #'h_muoneta',
            #'h_muonphi',
            #'h_muondphi',
            #'h_muonpfRelIso04',
        ]
        bkgSamples = ['qcd','wlnu','tqq','stqq','vvqq'] 
        #bkgSamples = ['qcd','tqq']
        variations = ['mutriggerUp','mutriggerDown','PuUp','PuDown','muisoUp','muisoDown','muidUp','muidDown','JESUp','JESDown','JERUp','JERDown']  
    else:
        plots = [
            'h_Lpass',
            'h_Mpass',
            'h_Tpass',
            'h_Lfail',
            'h_Mfail',
            'h_Tfail',
            #'h_cutflow',
            #'h_doubleB',
            'h_deepAK8',
            #'h_pT',
            #'h_pT_450',
            #'h_mSD',
            #'h_mSD_450',
            #'h_mSD_uncorr',
            #'h_rho',
            #'h_eta',
            #'h_n2ddt',
            #'h_n2ddt_new',
            #'h_n2',
            #'h_tau21ddt',
            #'h_tau21',
            #'h_MET',
            #'h_ttbarvar',
        ]

        bkgSamples = ['qcd','zqq','wqq','tqq'] #,'stqq','vvqq']
        variations =  ['triggerUp','triggerDown','PuUp','PuDown','matched','unmatched','JESUp','JESDown','JERUp','JERDown']

    dataSamples = ['data_obs']

    ifile = options.ifile
    lumi = options.lumi

    color = {'tqq':  ROOT.kMagenta+2,
             'zqq':  ROOT.kRed,
             'wqq':  ROOT.kGreen+2,
             'qcd': ROOT.kYellow-4,
             'stqq': ROOT.kTeal-2,
             'vvqq': ROOT.kAzure-4,
             'data_obs': ROOT.kBlack,
             'wlnu': ROOT.kGreen,
             'zll': ROOT.kPink,
             'tqq_matched': ROOT.kMagenta-1,
             'tqq_unmatched': ROOT.kMagenta-6,
             'tqq_semimatched': ROOT.kMagenta+6,
         }
    legname = {'tqq': 't#bar{t}+jets',
               'qcd': 'QCD',
               'data_obs': 'Data',
               'wqq': 'W+jets',
               'zqq': 'Z+jets',
               'stqq': 'Single t',
               'vvqq': 'VV',
               'wlnu': 'W(l#nu)',
               'zll': 'DY(ll)',
               'tqq_matched': 't#bar{t} merged',
               'tqq_unmatched': 't#bar{t} un-merged',
               'tqq_semimatched': 't#bar{t} semi-merged',
               }
    style = {'tqq': 1,
             'qcd': 1,
             'data_obs': 1,
             'wqq': 1,
             'zqq': 1,
             'stqq': 1,
             'vvqq': 1,
             'wlnu': 1,
             'zll': 1,
             'tqq_matched': 1,
             'tqq_unmatched': 1,
             'tqq_semimatched': 1,
             }

    canvases = []
    
    ofile = ROOT.TFile.Open(ifile,'read')
    hnew = {}
    hnew['L'] = {}
    hnew['M'] = {}
    hnew['T'] = {}
    for plot in plots:
        print(plot)
        hb = {}; hb2d = {}
        for process in bkgSamples:
            try:
                hb[process] = ofile.Get(plot.replace('h_',process+'_')).ProjectionX().Clone()
            except:
                print(plot.replace('h_',process+'_'))
                hb[process] = ofile.Get(plot.replace('h_',process+'_')).Clone()
            hb[process].SetDirectory(0)
            hb2d[process] = ofile.Get(plot.replace('h_',process+'_')).Clone()
            hb2d[process].SetDirectory(0)

        try:
            hd = ofile.Get(plot.replace('h_','data_obs_')).ProjectionX().Clone()
        except:
            hd = ofile.Get(plot.replace('h_','data_obs_')).Clone()
        hd.SetDirectory(0)
        hd2d = ofile.Get(plot.replace('h_','data_obs_')).Clone()
        hd2d.SetDirectory(0)

        if options.muonCR:
            new = hb.copy()
            if '_Mpass' in plot or '_Mfail' in plot:
                new['tqq_matched'] = ofile.Get(plot.replace('h_','tqq_')+'_matched').Clone()
                new['tqq_unmatched'] = ofile.Get(plot.replace('h_','tqq_')+'_unmatched').Clone()
                new['tqq_semimatched'] = ofile.Get(plot.replace('h_','tqq_')+'_semimatched').Clone()
                # remove tqq
                n = new.pop("tqq", None) 
            # scale QCD for muon CR 
            new['qcd'].Scale(0.8467)
            hb['qcd'].Scale(0.8467)
            #ratio = makeCanvasComparisonStackWData(hd,hb,legname,color,style,plot.replace('h_','stack_'),lumi,normalize='tqq')
            ratio = makeCanvasComparisonStackWData(hd,new,legname,color,style,plot.replace('h_','stack_'),lumi)
        else:
            ratio = makeCanvasComparisonStackWData(hd,hb,legname,color,style,plot.replace('h_','stack_'),lumi,normalize='qcd')

        print(hb.keys())
        for wp in ['L','M','T']: # ['L','M','T']:
            print(wp)
            if plot=='h_%spass'%wp or plot=='h_%sfail'%wp:
                for key in hb.keys():
                    hnew[wp][key+plot] = hb2d[key].Clone(hb2d[key].GetName().replace(wp,''))
                    hnew[wp][key+plot].SetDirectory(0)
                #if wp=='M':
                #if 'pass' in plot:
                #    hnew[wp]['qcd'+plot].Scale(0.501264267183)
                #else:
                #    hnew[wp]['qcd'+plot].Scale(0.410469905994)
                hnew[wp]['data_obs'+plot] = hd2d.Clone(hd2d.GetName().replace(wp,''))
                hnew[wp]['data_obs'+plot].SetDirectory(0)
                for var in variations:
                    for process in bkgSamples:
                        #if process in ['qcd','stqq','vvqq']: continue
                        if process in ['qcd']: continue
                        print(plot.replace('h_',process+'_')+'_'+var)
                        hnew[wp][var+plot+process] = ofile.Get(plot.replace('h_',process+'_')+'_'+var).Clone(plot.replace('h_',process+'_').replace(wp,'')+'_'+var)
                        hnew[wp][var+plot+process].SetDirectory(0)

    hqcd =ofile.Get("qcd_n2ddt_map").Clone()
    hqcd.SetDirectory(0)
    ofile.Close()

    #print(hnew['L'])
    print(hnew['T'])
    for wp in ['M','T']:
        newfile = ROOT.TFile.Open(ifile.replace('.root','_%s.root'%wp),'RECREATE')
        for key,h in hnew[wp].iteritems():
            h.Write()
        newfile.Close()
        # write loose templates for wqq and zqq
        if wp=='M' and not 'hbbSM' in ifile:
            newfile = ROOT.TFile.Open(ifile.replace('.root','_%s_looserWZ.root'%wp),'RECREATE')
            for process in ['wqq','zqq']:
                for cat in ['pass','fail']:
                    name = '%sh_L%s'%(process,cat)
                    h = hnew['L'][name].Clone()
                    h.Write()
                    for var in variations:
                        name = '%sh_L%s%s'%(var,cat,process)
                        h = hnew['L'][name].Clone()
                        h.Write()
            newfile.Close()

    #if not options.muonCR:
    # make ddt map
    #newfile = ROOT.TFile.Open(ifile.replace('.root','_qcdddtmap.root'),'RECREATE')
    #ddt = makeDDTmap(hqcd,0.26)
    #ddt.SetDirectory(0)
    #ddt.SetName("Rho2D")
    #    ddt.Write()
    #    newfile.Close()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option("--lumi", dest="lumi", type=float, default = 35.9,help="luminosity", metavar="lumi")
    parser.add_option('-i','--ifile', dest='ifile', default = 'input.root',help='directory with data', metavar='idir')
    parser.add_option("--muonCR", action='store_true', default =False, help="muonCR")
    (options, args) = parser.parse_args()

    import tdrstyle
    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.10)
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetPaintTextFormat("1.1f")
    ROOT.gStyle.SetOptFit(0000)
    ROOT.gROOT.SetBatch()

    main(options,args)
