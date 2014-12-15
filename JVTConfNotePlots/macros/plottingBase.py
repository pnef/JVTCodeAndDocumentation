#!/usr/bin/python
##################################################################
# pascal nef                             April 23rd, 2013        #
##################################################################
import sys
import os
import re
import math
import time
import numpy as n
from optparse import OptionParser
from sys import argv,exit
from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut, AddressOf
from ROOT import TMVA
import ROOT
import array
from multiprocessing import Process, Lock, Manager

##################################################################

class Plotter:
    def __init__(self, MultiProcess=False):
        self.manager_          = Manager()
        self.mp_               = MultiProcess
        self.logy_             = False
        self.textfont_         = 132
        self.logofont_         =  32
        self.atlasinfo_        = "Internal"
        self.normalized_       = False
        self.normalizedToData_ = False
        self.normalizeToLumi_  = False
        self.lumi_             = 1
        self.scale_            = 1
        self.drawRatio_        = False
        self.drawRatioUncert_  = ""
        self.atlasStyle_       = False
        self.saveHistos_       = False
        self.lock              = Lock()
        self.verbose_          = False
        self.ROC_no_1_1_point_ = False
        self.slacinfo_         = False
        self.customlabel_      = ""
        self.plotsavename_     = None
        self.ROC_rejection_    = False
        if(self.mp_):
            self.H_=self.manager_.dict()
            self.HlistRatio_=self.manager_.list()
            print "Plotter:__init__: enabling multiprocessing! " 
        else:
            self.H_={}
            self.HlistRatio_       = []

        ROOT.gROOT.ProcessLine(".L ~/bin/root_macros/MT2Helpers.hh");
        ROOT.gROOT.ProcessLine(".x ~/bin/root_macros/SetStyle_MT2Ratio.C");
        ROOT.gROOT.ProcessLine(".x ~/bin/root_macros/SetStyle_MT2.C");
        ROOT.gROOT.ForceStyle()

    def SetLogY(self, yesno):
        self.logy_ = yesno

    def SetATLASStyle(self):
        self.atlasStyle_ = True
        self.textfont_ = 42
        self.logofont_ = 72
        ROOT.gStyle.SetLegendFont(  42)
        ROOT.gStyle.SetLabelFont  ( 42   ,"X")
        ROOT.gStyle.SetLabelFont  ( 42   ,"Y")
        ROOT.gStyle.SetLabelFont  ( 42   ,"Z")
        ROOT.gStyle.SetTitleFont  ( 42   ,"X")
        ROOT.gStyle.SetTitleFont  ( 42   ,"Y")
        ROOT.gStyle.SetTitleFont  ( 42   ,"Z")
        ROOT.gStyle.SetTitleFont  ( 42)
        ROOT.gROOT.ForceStyle()

    def SetAtlasInfo(self, info):
        self.atlasinfo_=info

    def SetSLACInfo(self):
        self.textfont_ = 42
        self.logofont_ = 72
        ROOT.gStyle.SetLegendFont(  42)
        ROOT.gStyle.SetLabelFont  ( 42   ,"X")
        ROOT.gStyle.SetLabelFont  ( 42   ,"Y")
        ROOT.gStyle.SetLabelFont  ( 42   ,"Z")
        ROOT.gStyle.SetTitleFont  ( 42   ,"X")
        ROOT.gStyle.SetTitleFont  ( 42   ,"Y")
        ROOT.gStyle.SetTitleFont  ( 42   ,"Z")
        ROOT.gStyle.SetTitleFont  ( 42)
        ROOT.gROOT.ForceStyle()
        self.slacinfo_ = True

    def Set2DCanvas(self, YesNo, specialPalette=False):
        if YesNo:
            if specialPalette:
                ROOT.gROOT.ProcessLine(".x ~/bin/root_macros/SetStyle_MT2_2D_newPalette.C");
                ROOT.gROOT.ForceStyle()
            else:
                ROOT.gROOT.ProcessLine(".x ~/bin/root_macros/SetStyle_MT2_2D.C");
                ROOT.gROOT.ForceStyle()
        else:
            ROOT.gROOT.ProcessLine(".L ~/bin/root_macros/MT2Helpers.hh");
            ROOT.gROOT.ProcessLine(".x ~/bin/root_macros/SetStyle_MT2Ratio.C");
            ROOT.gROOT.ProcessLine(".x ~/bin/root_macros/SetStyle_MT2.C");
            ROOT.gROOT.ForceStyle()



    def DrawRatio(self, IDnum, IDden, IDratio, new=False):
        if new:
            self.Clone(ID, IDratio)
            for ibin in range(1,self.H_[IDratio].h_.GetNbinsX()+1):
                cont = self.H_[IDratio].h_.GetBinContent(ibin)
                err  = self.H_[IDratio].h_.GetBinError(ibin)
                self.H_[IDratio].h_.SetBinContent(ibin, 1)
                self.H_[IDratio].h_.SetBinError(ibin, math.fabs(err/cont))
        else:
            self.Divide(IDnum, IDden, IDratio, "", True)

        self.drawRatio_ = True
        self.HlistRatio_.extend([IDratio])
        ROOT.gROOT.SetStyle("MT2-ratio-Style");
        if (self.atlasStyle_): self.SetATLASStyle()
        ROOT.gROOT.UseCurrentStyle();
        ROOT.gROOT.ForceStyle();
        self.H_[IDratio].h_.UseCurrentStyle();
        ROOT.gROOT.SetStyle("MT2-Style");
        ROOT.gROOT.UseCurrentStyle();
        ROOT.gROOT.ForceStyle();
        self.H_[IDratio].setStyle()
    
    def DrawRatioMCUncert(self, ID, title, new=False):
        if new:
            self.Clone(ID, ID+"uncert")
            for ibin in range(1,self.H_[ID+"uncert"].h_.GetNbinsX()+1):
                cont = self.H_[ID+"uncert"].h_.GetBinContent(ibin)
                err  = self.H_[ID+"uncert"].h_.GetBinError(ibin)
                self.H_[ID+"uncert"].h_.SetBinContent(ibin, 1)
                self.H_[ID+"uncert"].h_.SetBinError(ibin, math.fabs(err/cont))
            self.drawRatio_ = True
            self.HlistRatio_.extend([ID+"uncert"])
            ROOT.gROOT.SetStyle("MT2-ratio-Style");
            if (self.atlasStyle_): self.SetATLASStyle()
            ROOT.gROOT.UseCurrentStyle();
            ROOT.gROOT.ForceStyle();
            self.H_[ID+"uncert"].h_.UseCurrentStyle();
            ROOT.gROOT.SetStyle("MT2-Style");
            ROOT.gROOT.UseCurrentStyle();
            ROOT.gROOT.ForceStyle();
            self.H_[ID+"uncert"].setStyle()
        else:
            self.DrawRatio(ID, ID, ID+"uncert")
        self.H_[ID+"uncert"].setStyle(1,921,"E2","", title ,3001)
        self.drawRatioUncert_ = ID+"uncert"

    def NormalizeToData(self,TrueOrFalse):
        self.normalizedToData_ = TrueOrFalse
    
    def NormalizeToLumi(self, TrueOrFalse):
        self.normalizeToLumi_ =TrueOrFalse

    def SetLumi(self, lumi):
        self.lumi_ = lumi

    def SetMCScale(self, scale=1, Hlist=[]):
        self.scale_ = scale
        for ID in Hlist:
            self.H_[ID].setScale(scale)

    def GetHisto(self, histname, identifier, tfile):
        H=Hist(histname, tfile)
        H.gethist()
        self.H_[identifier]=H

    def GetHistoFromTree(self, ID, tfile, tree, expression, selection, nbins, xmin, xmax, nbiny=0, ymin=-99, ymax=-99, nbinz=0, zmin=-99, zmax=-99):
        H=Hist(ID, tfile)
        H.gethistfromtree(tree, expression, selection, nbins, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax)
        self.H_[ID] = H
    
    def AddHistos(self, hlist, IDnew):
        if isinstance(hlist, str):
            raise Exception("AddHistos: argument must be a list...")
        else:
            htmp = self.H_[hlist[0]].h_.Clone(IDnew)
            for hname in hlist[1:]:
                htmp.Add(self.H_[hname].h_)
            H=Hist("", "")
            H.sethist(htmp)
            self.H_[IDnew]=H

    def Profile(self, ID, axis, IDnew):
        htmp=0
        if(axis is "x"):
            htmp=self.H_[ID].h_.ProfileX().Clone(IDnew)
        elif (axis is "y"):
            htmp=self.H_[ID].h_.ProfileY().Clone(IDnew)
        else:
            print "you suck"
        H=Hist("", "")
        H.sethist(htmp)
        self.H_[IDnew]=H

    def Projector(self, identifier, proj, newID):
        if(identifier is newID):
            if self.verbose_:
                print "Plotter::Projector: this is dangerous..."

        htmp=0
        try:
            htmp=self.H_[identifier].h_.Project3D(proj).Clone(newID)
        except AttributeError:
            if(proj is "y"):
                htmp=self.H_[identifier].h_.ProjectionY().Clone(newID)
            else:
                htmp=self.H_[identifier].h_.ProjectionX().Clone(newID)

        H=Hist("", "")
        H.sethist(htmp)
        self.H_[newID]=H
    
    def ProjectX(self, identifier, newID):
        htmp=self.H_[identifier].h_.ProjectionX().Clone(newID)
        H=Hist("", "")
        H.sethist(htmp)
        self.H_[newID]=H
    
    def ProjectZ(self, identifier, newID):
        htmp=self.H_[identifier].h_.ProjectionZ().Clone(newID)
        H=Hist("", "")
        H.sethist(htmp)
        self.H_[newID]=H

    def SetRange(self, ID, axis, cut1, cut2, newID):
        htmp=self.H_[ID].h_.Clone(newID)
        if(axis is "x"):
            bin1=1; bin2 =  htmp.GetXaxis().GetNbins()+1
            if(cut1 != -99): bin1 = htmp.GetXaxis().FindBin(cut1)
            if(cut2 != -99): bin2 = htmp.GetXaxis().FindBin(cut2)
            htmp.GetXaxis().SetRange(bin1, bin2-1)
        elif(axis is "y"):
            bin1=1; bin2 =  htmp.GetYaxis().GetNbins()+1
            if(cut1 != -99): bin1 = htmp.GetYaxis().FindBin(cut1)
            if(cut2 != -99): bin2 = htmp.GetYaxis().FindBin(cut2)
            htmp.GetYaxis().SetRange(bin1, bin2-1)
        elif(axis is "z"):
            bin1=1; bin2 =  htmp.GetZaxis().GetNbins()+1
            if(cut1 != -99): bin1 = htmp.GetZaxis().FindBin(cut1)
            if(cut2 != -99): bin2 = htmp.GetZaxis().FindBin(cut2)
            htmp.GetZaxis().SetRange(bin1, bin2-1)
        H=Hist("", "")
        H.sethist(htmp)
        self.H_[newID]=H

    def Clone(self, ID, IDnew):
        htmp = self.H_[ID].h_.Clone(IDnew)
        H=Hist("", "")
        H.sethist(htmp)
        self.H_[IDnew]=H

    def Divide(self, ID1, ID2, newID, binomial="", ratiostyle=False):
        htmp = self.H_[ID1].h_.Clone(newID)
        if binomial is "B":
            htmp.Divide(htmp, self.H_[ID2].h_,1,1,"B")
        elif binomial is "BayesDivide":
            htmp = ROOT.TGraphAsymmErrors()
            htmp.Divide(self.H_[ID1].h_, self.H_[ID2].h_,"cl=0.683 b(1,1) mode")
        else:
            htmp.Divide(self.H_[ID2].h_)
        H=Hist("", "")
        H.sethist(htmp)
        self.H_[newID]=H

    
    def Rebin(self, ID, axis, rebin, IDnew):
        htmp = self.H_[ID].h_.Clone(IDnew)
        try:
            if(axis is "x"): htmp=htmp.RebinX(rebin, "tmp")
            if(axis is "y"): htmp=htmp.RebinY(rebin, "tmp")
            if(axis is "z"): htmp=htmp.RebinZ(rebin, "tmp")
        except AttributeError:
            htmp=htmp.Rebin(rebin)
        H=Hist("", "")
        H.sethist(htmp)
        self.H_[IDnew]=H


    def MCSum(self, IDs, IDnew):
        tmp = self.H_[IDs[0]].h_.Clone(IDnew)
        for ID in IDs[1:]:
            tmp.Add(self.H_[ID].h_)
        H=Hist("", "")
        H.sethist(tmp)
        self.H_[IDnew]=H


    def JetMultiplicity1(self, hist, refhist,  xmin, xmax, ymin, ymax, zmin, zmax, xrebin, newID):
        tmpID1=hist+"x"+str(xmin)+"to"+str(xmax)
        tmpID2=tmpID1+"y"+str(ymin)+"to"+str(ymax)
        tmpID3=tmpID2+"z"+str(zmin)+"to"+str(zmax)
        self.SetRange (hist,   "x", xmin, xmax, tmpID1)
        self.SetRange (tmpID1, "y", ymin, ymax, tmpID2)
        self.SetRange (tmpID2, "z", zmin, zmax, tmpID3)

        self.Projector(tmpID3, "x", tmpID3)
        self.Rebin(tmpID3, "x", xrebin, tmpID3)
        self.Rebin(refhist,"x", xrebin, refhist+"rebin")
        self.Divide(tmpID3, refhist+"rebin", newID)

    def JetMultiplicity(self, hist, xmin, xmax, ymin, ymax, zmin, zmax, refhist, yminref, ymaxref, xrebin, newID):
        tmp   =hist   +"x"+str(xmin)+"to"+str(xmax)+"y"+str(ymin)+"to"+str(ymax)+"z"+str(zmin)+"to"+str(zmax)
        tmpref=refhist+"y"+str(yminref)+"to"+str(ymaxref)

        self.SetRange (hist,   "x", xmin, xmax, tmp)
        self.SetRange (tmp,    "y", ymin, ymax, tmp)
        self.SetRange (tmp,    "z", zmin, zmax, tmp)

        self.SetRange (refhist,"y", yminref, ymaxref, tmpref)

        self.Projector(tmp,    "x", tmp)
        self.ProjectX(tmpref,  tmpref)

        self.Rebin(tmp,   "x", xrebin, tmp)
        self.Rebin(tmpref,"x", xrebin, tmpref)
        self.Divide(tmp, tmpref, newID)

    def GetBinForEff(self, ID, axis, eff, HSonTheRight=True):
        self.Projector(ID, axis, "htmp")
        self.H_["htmp"].h_.Scale(1./self.H_["htmp"].h_.Integral())
        nbins = self.H_["htmp"].h_.GetNbinsX()
#        self.H_["htmp"].h_.DrawNormalized()
#        wait()
        bin=0;
        for i in range(nbins):
            tmpbin = i+1
            if(HSonTheRight     and self.H_["htmp"].h_.Integral(tmpbin,nbins) > eff):
#                print bin, self.H_["htmp"].h_.Integral(tmpbin,nbins), self.H_["htmp"].h_.GetBinLowEdge(bin)
                bin = tmpbin
            if(not HSonTheRight and 1-self.H_["htmp"].h_.Integral(tmpbin,nbins) > eff):
                bin = tmpbin
                break
        if HSonTheRight:
            eff1 = self.H_["htmp"].h_.Integral(bin  ,nbins) ; print eff1
            eff2 = self.H_["htmp"].h_.Integral(bin+1,nbins) ; print eff2
            if(math.fabs(eff1-eff) < math.fabs(eff2-eff)):
                return {"bin":bin,"cut":round(self.H_["htmp"].h_.GetBinLowEdge(bin),20)}
            else:
                return {"bin":bin+1,"cut":round(self.H_["htmp"].h_.GetBinLowEdge(bin+1),20)}
        else: 
            eff1 = 1-self.H_["htmp"].h_.Integral(bin  ,nbins) ; print eff1
            eff2 = 1-self.H_["htmp"].h_.Integral(bin-1,nbins) ; print eff2
            if(math.fabs(eff1-eff) < math.fabs(eff2-eff)):
                return {"bin":bin,"cut":round(self.H_["htmp"].h_.GetBinLowEdge(bin),20)}
            else:
                return {"bin":bin-1,"cut":round(self.H_["htmp"].h_.GetBinLowEdge(bin-1),20)}


    def EfficiencyVsNPV(self, ID, xrebin, ymin, ymax, zmin, zmax, newID):
        self.Rebin(ID, "x", xrebin, newID)
        self.Rebin(ID, "x", xrebin, ID+"tmp")
        self.SetRange(newID,    "y", ymin, ymax, newID)
        self.SetRange(ID+"tmp", "y", ymin, ymax, ID+"tmp")

        self.SetRange(newID, "z", zmin, zmax, newID) # do not apply cut to reference histo

        self.Projector(newID, "x", newID+"projx")
        self.Projector(ID+"tmp", "x", ID+"tmp")

        self.Divide(newID+"projx", ID+"tmp", newID, "B")
#        self.H_[newID].h_.Draw("same")
#        wait()
#        for bin in range(self.H_[newID+"projx"].h_.GetNbinsX()):
#            print newID+"projx", bin+1, "center", self.H_[newID+"projx"].h_.GetBinCenter(bin+1), "content num", self.H_[newID+"projx"].h_.GetBinContent(bin+1),\
#                  "content den", self.H_[ID+"tmp"].h_.GetBinContent(bin+1), "eff", self.H_[newID].h_.GetBinContent(bin+1), "err", self.H_[newID].h_.GetBinError(bin+1), \
#                  "new err", math.sqrt(self.H_[newID].h_.GetBinContent(bin+1)*(1-self.H_[newID].h_.GetBinContent(bin+1))/self.H_[ID+"tmp"].h_.GetBinContent(bin+1)  )

    def ROC(self, IDhs, IDpu, xmin, xmax, ymin, ymax, zmin, zmax, newID, HSonTheRight=True):
        self.SetRange(IDhs         , "x", xmin, xmax, IDhs+"forROC")
        self.SetRange(IDpu         , "x", xmin, xmax, IDpu+"forROC")
        self.SetRange(IDhs+"forROC", "y", ymin, ymax, IDhs+"forROC")
        self.SetRange(IDpu+"forROC", "y", ymin, ymax, IDpu+"forROC")
        self.SetRange(IDhs+"forROC", "z", zmin, zmax, IDhs+"forROC")
        self.SetRange(IDpu+"forROC", "z", zmin, zmax, IDpu+"forROC")
        self.Projector(IDpu+"forROC","z", IDpu+"forROC")
        self.Projector(IDhs+"forROC","z", IDhs+"forROC")
        self.H_[IDpu+"forROC"].h_.Scale(1./self.H_[IDpu+"forROC"].h_.Integral())
        self.H_[IDhs+"forROC"].h_.Scale(1./self.H_[IDhs+"forROC"].h_.Integral())
        nbins     = self.H_[IDpu+"forROC"].h_.GetNbinsX()
        eff       = n.zeros(nbins, dtype=float)
        fake      = n.zeros(nbins, dtype=float)
        reject    = n.zeros(nbins, dtype=float)
        efferr    = n.zeros(nbins, dtype=float)
        fakeerr   = n.zeros(nbins, dtype=float)
        rejecterr = n.zeros(nbins, dtype=float)
        bincounter=0
        for i in range(nbins):
            ibin     = i+1
            efferror, fakeerror, failerr = ROOT.Double(0), ROOT.Double(0), ROOT.Double(0)
            if HSonTheRight:
                tmpeff    = self.H_[IDhs+"forROC"].h_.IntegralAndError(ibin, nbins, efferror)
                tmpfake   = self.H_[IDpu+"forROC"].h_.IntegralAndError(ibin, nbins, fakeerror) 
                if tmpeff >0.00000001 and not (tmpeff > 0.99999 and self.ROC_no_1_1_point_ == True) and not (tmpfake==0 and self.ROC_rejection_ == True) :
                    eff       [bincounter]=self.H_[IDhs+"forROC"].h_.IntegralAndError(ibin, nbins, efferror)
                    fail                  =self.H_[IDhs+"forROC"].h_.IntegralAndError(1,    ibin-1,failerr)
                    fake      [bincounter]=self.H_[IDpu+"forROC"].h_.IntegralAndError(ibin, nbins, fakeerror)
                    reject    [bincounter]=1/fake[bincounter]
                    efferr    [bincounter]=math.sqrt(math.pow(efferror*(1. - eff[bincounter]), 2) + math.pow(failerr*eff[bincounter], 2))
                    fakeerr   [bincounter]=fakeerror
                    rejecterr [bincounter]=1/(fake[bincounter]*fake[bincounter])*fakeerror 
#                    print  bincounter,eff [bincounter], fake[bincounter],  reject    [bincounter]
                    if (eff[bincounter] >0.99999  and fake [bincounter] > 0.99999): efferr [bincounter]=0
                    bincounter+=1
#                    print "bin", ibin, "cut", self.H_[IDpu+"forROC"].h_.GetBinLowEdge(ibin), "eff", eff[bincounter-1],  "fake", fake   [bincounter-1], "fakeerr", fakeerr [bincounter-1] 
            else:
                tmpeff    = 1-self.H_[IDhs+"forROC"].h_.IntegralAndError(ibin, nbins, efferror)
                if tmpeff >0.00000001 and not ( tmpeff > 0.99999 and self.ROC_no_1_1_point_ == True):
                    fail               =self.H_[IDhs+"forROC"].h_.IntegralAndError(ibin, nbins, failerr)
                    eff    [bincounter]=self.H_[IDhs+"forROC"].h_.IntegralAndError(1,    ibin-1,efferror)
                    fake   [bincounter]=self.H_[IDpu+"forROC"].h_.IntegralAndError(1,    ibin-1,fakeerror)
                    reject [bincounter]=1/fake[bincounter]
                    efferr [bincounter]=math.sqrt(math.pow(efferror*(1. - eff[bincounter]), 2) + math.pow(failerr*eff[bincounter], 2))
                    fakeerr[bincounter]=fakeerror
                    rejecterr [bincounter]=1/(fake[bincounter]*fake[bincounter])*fakeerror 
                    if (eff[bincounter] >0.99999  and fake [bincounter] > 0.99999): efferr [bincounter]=0
                    bincounter+=1
#                print "bin", ibin, "cut", self.H_[IDpu+"forROC"].h_.GetBinLowEdge(ibin), "eff", eff[i], "fake", fake   [i]

        if (self.ROC_rejection_ == False):
            tgraph = ROOT.TGraphErrors(bincounter, eff, fake, efferr, fakeerr)
        else:
            tgraph = ROOT.TGraphErrors(bincounter, eff, reject, efferr, rejecterr)
#        print eff, fake, fakeerr
        H=Hist("", "")
        H.sethist(tgraph)
        self.H_[newID]=H
#        self.H_[IDpu+"forROC"].h_.Draw("")
#        self.H_[IDhs+"forROC"].h_.Draw("same")
#        tgraph.Draw("ACP")
#        wait()
    
    def ROCforFixedCut(self, IDhs, IDpu, xmin, xmax, ymin, ymax, zmin, zmax, newID, HSonTheRight, varcut):
        self.SetRange(IDhs         , "x", xmin, xmax, IDhs+"forROC")
        self.SetRange(IDpu         , "x", xmin, xmax, IDpu+"forROC")
        self.SetRange(IDhs+"forROC", "y", ymin, ymax, IDhs+"forROC")
        self.SetRange(IDpu+"forROC", "y", ymin, ymax, IDpu+"forROC")
        self.SetRange(IDhs+"forROC", "z", zmin, zmax, IDhs+"forROC")
        self.SetRange(IDpu+"forROC", "z", zmin, zmax, IDpu+"forROC")
        nbins = self.H_[IDpu].h_.GetNbinsX()
        eff    = n.zeros(nbins+2, dtype=float)
        fake   = n.zeros(nbins+2, dtype=float)
        efferr = n.zeros(nbins+2, dtype=float)
        fakeerr= n.zeros(nbins+2, dtype=float)
        bincounter=0
        for i in range(1,nbins+1):
            self.H_[IDhs+"forROC"].h_.GetXaxis().SetRange(i, i+1)
            self.H_[IDpu+"forROC"].h_.GetXaxis().SetRange(i, i+1)
            self.Projector(IDpu+"forROC","z", IDpu+"forROC2")
            self.Projector(IDhs+"forROC","z", IDhs+"forROC2")
            self.H_[IDpu+"forROC2"].h_.Scale(1./self.H_[IDpu+"forROC2"].h_.Integral())
            self.H_[IDhs+"forROC2"].h_.Scale(1./self.H_[IDhs+"forROC2"].h_.Integral())
            cutbin   = self.H_[IDhs+"forROC2"].h_.FindBin(varcut)
            nCutbins = self.H_[IDhs+"forROC2"].h_.GetNbinsX()
            efferror, fakeerror, failerr = ROOT.Double(0), ROOT.Double(0), ROOT.Double(0)
            if HSonTheRight:
                tmpeff    = self.H_[IDhs+"forROC2"].h_.IntegralAndError(cutbin, nCutbins, efferror)
                if tmpeff >0.00000001:
                    eff    [bincounter]=self.H_[IDhs+"forROC2"].h_.IntegralAndError(cutbin, nCutbins, efferror)
                    fail               =self.H_[IDhs+"forROC2"].h_.IntegralAndError(1,      cutbin-1,failerr)
                    fake   [bincounter]=self.H_[IDpu+"forROC2"].h_.IntegralAndError(cutbin, nCutbins, fakeerror)
                    efferr [bincounter]=math.sqrt(math.pow(efferror*(1. - eff[i]), 2) + math.pow(failerr*eff[i], 2))
                    fakeerr[bincounter]=fakeerror
                    if (eff[bincounter] >0.99999  and fake [bincounter] > 0.99999): efferr [bincounter]=0
                    bincounter+=1
            else:
                tmpeff    = 1-self.H_[IDhs+"forROC2"].h_.IntegralAndError(cutbin, nCutbins, efferror)
                if tmpeff >0.00000001:
                    fail               =self.H_[IDhs+"forROC2"].h_.IntegralAndError(cutbin, nCutbins, failerr)
                    eff    [bincounter]=self.H_[IDhs+"forROC2"].h_.IntegralAndError(1,    cutbin-1,efferror)
                    fake   [bincounter]=self.H_[IDpu+"forROC2"].h_.IntegralAndError(1,    cutbin-1,fakeerror)
                    efferr [bincounter]=math.sqrt(math.pow(efferror*(1. - eff[i]), 2) + math.pow(failerr*eff[i], 2))
                    fakeerr[bincounter]=fakeerror
                    if (eff[bincounter] >0.99999  and fake [bincounter] > 0.99999): efferr [bincounter]=0
                    bincounter+=1
            print "for x in", self.H_[IDhs+"forROC"].h_.GetXaxis().GetBinLowEdge(i), "and", self.H_[IDhs+"forROC"].h_.GetXaxis().GetBinLowEdge(i+1), "efficiency", eff    [bincounter-1] , "fake", fake    [bincounter-1]
        # set dummy values for fix plotting range    
        eff    [nbins]=0;    eff    [nbins+1]=1.1
        fake   [nbins]=0.1;  fake   [nbins+1]=0.1
        efferr [nbins]=0;    efferr [nbins+1]=0
        fakeerr[nbins]=0;    fakeerr[nbins+1]=0

        tgraph = ROOT.TGraphErrors(bincounter+2, eff, fake, efferr, fakeerr)
        H=Hist("", "")
        H.sethist(tgraph)
        self.H_[newID]=H
#        self.H_[IDpu+"forROC2"].h_.Draw("")
#        self.H_[IDhs+"forROC2"].h_.Draw("same")
#        tgraph.Draw("ACP")
#        wait()

    def NJetsVsPtThreshold(self, hID, hIDref, xmin, xmax, ymin, ymax, zmin, zmax, hIDnew):
    # compute average njets as a function of minimal pt threshold, 
    # -> integrated over all x ranges (NPV)
    # -> for interval of z ranges (mu)
        self.SetRange(hID,          "x", xmin, xmax, hIDnew)
        self.SetRange(hIDref,       "x", xmin, xmax, hIDnew+"ref")
        self.SetRange(hIDnew,       "y", ymin, ymax, hIDnew)
        self.SetRange(hIDnew,       "z", zmin, zmax, hIDnew)
        self.SetRange(hIDnew+"ref", "y", zmin, zmax, hIDnew+"ref") # !!! href is 2D X vs Y: apply zmin, zmax cut to y-axis
        self.Projector(hIDnew+"ref","x", hIDnew+"ref")
        self.Projector(hIDnew,"yx", hIDnew)
        nptbins  = self.H_[hIDnew].h_.GetYaxis().GetNbins()
        njets    = n.zeros(nptbins, dtype=float)
        ptthres  = n.zeros(nptbins, dtype=float)
        njetserr = n.zeros(nptbins, dtype=float)
        pterr    = n.zeros(nptbins, dtype=float)
        for iBin in range(nptbins):
#            print self.H_[hIDnew].h_.GetYaxis().GetBinLowEdge(iBin+1)
            ptthres [iBin] = self.H_[hIDnew].h_.GetYaxis().GetBinLowEdge(iBin+1)
            self.SetRange(hIDnew, "y", ptthres [iBin], ymax, "htmp")
            self.Projector("htmp",  "x", "htmp")
            err1, err2 =ROOT.Double(0), ROOT.Double(0)
            int1       = self.H_["htmp"].h_.IntegralAndError(1,-1,err1) 
            int2       = self.H_[hIDnew+"ref"].h_.IntegralAndError(1,-1,err2)
            njets   [iBin] = int1/int2
            try:
                njetserr[iBin] = njets   [iBin] * math.sqrt(math.pow(err1/int1,2)+math.pow(err2/int2,2))
            except ZeroDivisionError:
                njetserr[iBin] = njets   [iBin] * math.sqrt(math.pow(err2/int2,2))
#            print "pt", ptthres [iBin] ,  "njets", njets    [iBin] 
        tgraph = ROOT.TGraphErrors(nptbins, ptthres, njets, pterr, njetserr )
        H=Hist("", "")
        H.sethist(tgraph)
        self.H_[hIDnew]=H
#        tgraph.Draw("ACP")

    def GetSlope(self, ID, IDnew):
        nbins  = self.H_[ID].h_.GetNbinsX()
        xbin   = n.zeros(nbins, dtype=float)
        xbinerr= n.zeros(nbins, dtype=float)
        slope  = n.zeros(nbins, dtype=float)
        err    = n.zeros(nbins, dtype=float)
        for iBin in range(0,nbins-1):
            y1     = self.H_[ID].h_.GetBinContent(iBin+1)
            y1err  = self.H_[ID].h_.GetBinError(iBin+1)
            y2     = self.H_[ID].h_.GetBinContent(iBin+2)
            y2err  = self.H_[ID].h_.GetBinError(iBin+2)
            DeltaX = self.H_[ID].h_.GetBinCenter(iBin+2)-self.H_[ID].h_.GetBinCenter(iBin+1)
            slope[iBin]=(y2-y1)/DeltaX
            err  [iBin]=1/DeltaX*math.sqrt(y1err*y1err + y2err*y2err)
            xbin [iBin]=(self.H_[ID].h_.GetBinCenter(iBin+2) + self.H_[ID].h_.GetBinCenter(iBin+1))/2
            xbinerr[iBin]=0
            print iBin, y1, y2, DeltaX, slope[iBin]

        tgraph = ROOT.TGraphErrors(nbins-1, xbin, slope,xbinerr, err)
        H=Hist("", "")
        H.sethist(tgraph)
        self.H_[IDnew]=H

    def DrawRMSvsVariable(self, ID, rmsAxis, variableAxis, variableRebin, IDnew, otherAxis="", OtherMin=-99, OtherMax=-99):
        self.Rebin(ID, variableAxis, variableRebin, ID+"tmp")
        if(otherAxis is not ""):
            self.SetRange(ID+"tmp", otherAxis, OtherMin, OtherMax, ID+"tmp")

        axis = 0
        if (variableAxis=="x"):
            axis = self.H_[ID+"tmp"].h_.GetXaxis()
        elif (variableAxis=="y"):
            axis = self.H_[ID+"tmp"].h_.GetYaxis()
        elif (variableAxis=="z"):
            axis = self.H_[ID+"tmp"].h_.GetZaxis()

        nbins  = axis.GetNbins()
        xbin   = n.zeros(nbins, dtype=float)
        xbinerr= n.zeros(nbins, dtype=float)
        rms    = n.zeros(nbins, dtype=float)
        rmserr = n.zeros(nbins, dtype=float)
        for ibin in range(1, nbins+1):
            self.SetRange(ID+"tmp", variableAxis, axis.GetBinLowEdge(ibin), axis.GetBinLowEdge(ibin+1), ID+"tmp2")
            self.Projector(ID+"tmp2", rmsAxis, ID+"tmp3")
            rms    [ibin-1] = self.H_[ID+"tmp3"].h_.GetRMS() 
            rmserr [ibin-1] = self.H_[ID+"tmp3"].h_.GetRMSError() 
            xbin   [ibin-1] = axis.GetBinCenter(ibin)
            xbinerr[ibin-1]= 0
            print ibin,ibin+1, xbin[ibin-1], rms[ibin-1]
        print xbin, rms

        tgraph = ROOT.TGraphErrors(nbins, xbin, rms,xbinerr, rmserr)
        H=Hist("", "")
        H.sethist(tgraph)
        self.H_[IDnew]=H

    def DrawMeanvsVariable(self, ID, meanAxis, variableAxis, variableRebin, IDnew, otherAxis="", OtherMin=-99, OtherMax=-99):
        self.Rebin(ID, variableAxis, variableRebin, ID+"tmp")
        if(otherAxis is not ""):
            self.SetRange(ID+"tmp", otherAxis, OtherMin, OtherMax, ID+"tmp")

        axis = 0
        if (variableAxis=="x"):
            axis = self.H_[ID+"tmp"].h_.GetXaxis()
        elif (variableAxis=="y"):
            axis = self.H_[ID+"tmp"].h_.GetYaxis()
        elif (variableAxis=="z"):
            axis = self.H_[ID+"tmp"].h_.GetZaxis()

        nbins  = axis.GetNbins()
        xbin   = n.zeros(nbins, dtype=float)
        xbinerr= n.zeros(nbins, dtype=float)
        mean    = n.zeros(nbins, dtype=float)
        meanerr = n.zeros(nbins, dtype=float)
        for ibin in range(1, nbins+1):
            self.SetRange(ID+"tmp", variableAxis, axis.GetBinLowEdge(ibin), axis.GetBinLowEdge(ibin+1), ID+"tmp2")
            self.Projector(ID+"tmp2", meanAxis, ID+"tmp3")
            mean    [ibin-1] = self.H_[ID+"tmp3"].h_.GetMean() 
            meanerr [ibin-1] = self.H_[ID+"tmp3"].h_.GetMeanError() 
            xbin    [ibin-1] = axis.GetBinCenter(ibin)
            xbinerr [ibin-1]= 0
            print ibin,ibin+1, xbin[ibin-1], mean[ibin-1]
        print xbin, mean

        tgraph = ROOT.TGraphErrors(nbins, xbin, mean,xbinerr, meanerr)
        H=Hist("", "")
        H.sethist(tgraph)
        self.H_[IDnew]=H

    def GetResolutionVsPtTrueFromResponse(self, ID, responseAxis, ptTrueAxis, IDnew):
        self.SetRange(ID, ptTrueAxis, 20, 30, ID+"tmp")
        self.Projector(ID+"tmp", responseAxis, IDnew)



        
    def GetPtThresholdVsEtaForFakeRate(self, ID_EtaVsPT_PU, ID_EtaVsPt_Incl, fakerate, IDnew):
        axis       = self.H_[ID_EtaVsPT_PU].h_.GetXaxis()
        nbins      = axis.GetNbins()
        xbin       = n.zeros(nbins, dtype=float)
        xbinerr    = n.zeros(nbins, dtype=float)
        ptthres    = n.zeros(nbins, dtype=float)
        ptthreserr = n.zeros(nbins, dtype=float)
        for ibin in range(1, nbins+1):
            xmin = self.H_[ID_EtaVsPT_PU].h_.GetXaxis().GetBinLowEdge(ibin)
            xmax = self.H_[ID_EtaVsPT_PU].h_.GetXaxis().GetBinLowEdge(ibin+1)
            print "------------------- eta ", xmin, xmax, "----------------------------"
            self.SetRange(ID_EtaVsPt_Incl, "x", xmin, xmax, "htmp_Incl")
            self.SetRange(ID_EtaVsPT_PU,   "x", xmin, xmax, "htmp_PU")
            self.Projector("htmp_Incl", "y", "htmp_Incl_y")
            self.Projector("htmp_PU",   "y", "htmp_PU_y")
            thresholdForFakeRate = 0
            efffakeRate          = 0
            for ptbin in range(1, self.H_["htmp_Incl_y"].h_.GetNbinsX()+1):
                if(self.H_["htmp_Incl_y"].h_.Integral(ptbin, -1) >0):
                    fakeRate = self.H_["htmp_PU_y"].h_.Integral(ptbin, -1)/self.H_["htmp_Incl_y"].h_.Integral(ptbin, -1)
                    print ptbin, self.H_["htmp_Incl_y"].h_.GetBinLowEdge(ptbin), fakeRate
                    if(fakeRate <=fakerate):
                        thresholdForFakeRate    = self.H_["htmp_Incl_y"].h_.GetBinLowEdge(ptbin)
                        efffakeRate             = fakeRate
                        break
            print "fakeRate", efffakeRate, "for ptthres", thresholdForFakeRate, "and eta " , (xmin+xmax)/2
            xbin     [ibin-1] =  (xmin+xmax)/2
            ptthres  [ibin-1] =  thresholdForFakeRate
        tgraph = ROOT.TGraphErrors(nbins, xbin, ptthres, xbinerr, ptthreserr)
        H=Hist("", "")
        H.sethist(tgraph)
        self.H_[IDnew]=H
                    
    def GetJVFSystHist2(self, ID, IDnew):
        nbins = self.H_[ID].h_.GetXaxis().GetNbins()
        result= {'bin':-999, 'cut':-999}
        hNew  = self.H_[ID].h_.Clone("hNew"); hNew.Sumw2()
        for ibin in range(1, nbins+1):
            if(self.H_[ID]  .h_.Integral()>0):
                hNew.SetBinContent(ibin, self.H_[ID]  .h_.Integral(ibin, -1) / self.H_[ID]  .h_.Integral())
                hNew.SetBinError  (ibin, 0)
            else:
                hNew.SetBinContent(ibin, 0)
                hNew.SetBinError  (ibin, 0)
        H=Hist("", ""); H.sethist(hNew); self.H_[IDnew]=H



    def GetJVFSyst(self, IDdata, IDMC, efficiency):
        # get nominal cut value for MC to reach target efficiency
        resultMC     = {}
        resultData   = {}
        resultMC     = self.GetCutAndBinForEfficiency(IDMC, efficiency)

        # get nominal eff in data with cut from above
        resultData['eff']      = self.H_[IDdata]  .h_.Integral(resultMC['bin'], -1) / self.H_[IDdata]  .h_.Integral()
        resultMC  ['deltaEff'] = math.fabs(resultMC['eff'] - resultData['eff'])


        # get cutUpDown for MC to get eff = efficiency +- DeltaEff
        resultMCup   = self.GetCutAndBinForEfficiency(IDMC, efficiency+resultMC['deltaEff'], "up")
        resultMCdown = self.GetCutAndBinForEfficiency(IDMC, efficiency-resultMC['deltaEff'], "down")

        resultMC['effUp']   = resultMCup['eff']
        resultMC['binUp']   = resultMCup['bin']
        resultMC['cutUp']   = resultMCup['cut']
        resultMC['effDown'] = resultMCdown['eff']
        resultMC['binDown'] = resultMCdown['bin']
        resultMC['cutDown'] = resultMCdown['cut']
        return   [resultMC, resultData]

#        print [resultMC, resultData]


    def GetCutAndBinForEfficiency(self, ID, eff, updown=None):    
        nbins = self.H_[ID].h_.GetXaxis().GetNbins()
        result= {'bin':-999, 'cut':-999}
        for ibin in range(1, nbins+1):
            tmpeff   = self.H_[ID]  .h_.Integral(ibin, -1) / self.H_[ID]  .h_.Integral()
            if(tmpeff   < eff): 
                tmpeff2 = self.H_[ID]  .h_.Integral(ibin-1, -1) / self.H_[ID]  .h_.Integral()
                if(math.fabs(tmpeff2 - eff) < math.fabs(tmpeff - eff)):
                    result['bin']=ibin-1
                    result['eff']=tmpeff2
                    result['cut']=self.H_[ID]  .h_.GetBinLowEdge(ibin-1)
                else:
                    result['bin']=ibin
                    result['eff']=tmpeff
                    result['cut']=self.H_[ID]  .h_.GetBinLowEdge(ibin)
                if(updown == "up"):
                    result['bin']=ibin-1
                    result['eff']=tmpeff2
                    result['cut']=self.H_[ID]  .h_.GetBinLowEdge(ibin-1)
                if(updown == "down"):
                    result['bin']=ibin
                    result['eff']=tmpeff
                    result['cut']=self.H_[ID]  .h_.GetBinLowEdge(ibin)
                break
        return result



    def GetCutForEffValue(self, ID, axisforbinning):
        if(axisforbinning=="x"):
            nbins = self.H_[ID].h_.GetXaxis().GetNbins()
        elif (axisforbinning=="y"):
            nbins = self.H_[ID].h_.GetYaxis().GetNbins()

        quant       = n.zeros(1, dtype=float)
        xvals       = n.zeros(1, dtype=float)
        quant[0]=0.9

        for ibin in range(1, nbins):
            self.H_[ID].h_.ProjectionY("Proj",ibin, ibin+1).GetQuantiles(1,xvals,quant);
            print xvals, quant
    
            

    def Draw(self, Hlist, xtitle="x", ytitle="y", normalized=False, ATLASlogo=False, labels=[], legends=False):

        # check if Hlist is a list of strings of a list of Hist objects
        hlist = []
        dataNEvents = 0
        for h in Hlist:
            try:
                self.H_.has_key(h)
                hlist.append(self.H_[h])
            except KeyError:
                hlist.append(h)

        for h in hlist:
            if not h in self.H_.values() and not self.mp_:
                print "values", self.H_.values()
                print "ERROR: " , h, "is not in list of histos",  self.H_
                exit(1)
            if h.isData_:
                dataNEvents=h.h_.Integral()

        #--------

        canv=0
        if(self.drawRatio_):
            canv = ROOT.MT2RatioCanvas()
            canv.UseMainPad()
        else:
            canv = ROOT.TCanvas()

        if self.logy_: ROOT.gPad.SetLogy()
        legnentries=0
        for hist in hlist: # n entries for legend
            if(hist.addToLeg_): legnentries+=1
        leg = ROOT.TLegend(0.5, 0.85-legnentries*0.04, 0.8, 0.85);
#        leg = ROOT.TLegend(0.2, 0.68-legnentries*0.04, 0.5, 0.68);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.SetTextFont(self.textfont_);
        solidcount=0; opencount =0; filledcounter=0; legcount  =0; linecount =0; datacount=0; filledwithlinecount=0
            
        for index, hist in zip(range(len(hlist)), hlist):
            if hist.isData_:
                hist.h_.SetLineColor(1)
#                hist.h_.SetFillColor(1)
                hist.h_.SetMarkerColor(1)
                hist.h_.SetMarkerStyle(self.MARKERSData[datacount])
                datacount+=1
            elif not hist.styleset_:
                if hist.filledWithLine_:
                    hist.h_.SetLineColor  (self.COLORSLine[filledwithlinecount])
                    hist.h_.SetFillColor  (self.COLORSFilled[filledwithlinecount])
                    hist.h_.SetFillStyle  (1001)
                    filledwithlinecount+=1
                elif hist.filled_:
                    hist.h_.SetLineColor  (self.COLORSFilled[filledcounter])
                    hist.h_.SetFillColor  (self.COLORSFilled[filledcounter])
                    hist.h_.SetFillStyle  (self.STYLESFilled[filledcounter])
                    filledcounter+=1
                elif hist.solid_:
                    hist.h_.SetLineColor  (self.COLORSSolid[solidcount])
                    hist.h_.SetMarkerColor(self.COLORSSolid[solidcount])
                    hist.h_.SetMarkerStyle(self.MARKERSSolid[solidcount])
                    solidcount+=1
                elif hist.line_:
                    hist.h_.SetLineColor  (self.COLORSLine[linecount])
                    hist.h_.SetMarkerColor(self.COLORSLine[linecount])
                    linecount+=1
                else:
                    hist.h_.SetLineColor  (self.COLORSOpen[opencount])
                    hist.h_.SetMarkerColor(self.COLORSOpen[opencount])
                    hist.h_.SetMarkerStyle(self.MARKERSOpen[opencount])
                    opencount+=1
            if(legends and hist.addToLeg_): 
                if hist.legstyle_ != "":
                    leg.AddEntry( hist.h_, legends[legcount], hist.legstyle_)
                elif hist.filled_ or hist.drawstyle_=="box":
                    leg.AddEntry( hist.h_, legends[legcount], "f")
                elif hist.line_:
                    leg.AddEntry( hist.h_, legends[legcount], "l")
                else:
                    leg.AddEntry( hist.h_, legends[legcount], "lp")
                legcount+=1
            if (index == 0):
                hist.h_.GetXaxis().SetTitle(xtitle)
                hist.h_.GetYaxis().SetTitle(ytitle)
                if(self.normalizedToData_):
                    if not hist.isData_:
                        hist.h_.Scale(dataNEvents / hist.h_.Integral() )
                    hist.h_.Draw(hist.drawstyle_)
                elif (self.normalizeToLumi_):
                    if not hist.isData_:
                        hist.h_.Scale(self.lumi_)
                        print "scale", self.lumi_
                    hist.h_.Draw(hist.drawstyle_)
                else:
                    if(normalized or self.normalized_):
                        hist.h_.DrawNormalized(hist.drawstyle_)
                    else:
                        hist.h_.Draw(hist.drawstyle_)
            else:
                if(self.normalizedToData_):
                    if not hist.isData_:
                        hist.h_.Scale(dataNEvents / hist.h_.Integral() )
                    hist.h_.Draw(hist.drawstyle_+"same")
                elif (self.normalizeToLumi_):
                    if not hist.isData_:
                        hist.h_.Scale(self.lumi_)
                        print "scale", self.lumi_
                    hist.h_.Draw(hist.drawstyle_+"same")
                else:
                    if(normalized or self.normalized_):
                        hist.h_.DrawNormalized(hist.drawstyle_+"same")
                    else:
                        hist.h_.Draw(hist.drawstyle_+"same")
        leg.Draw()
        if self.customlabel_ == "":
            if not self.slacinfo_:
                if(ATLASlogo): ATLASInternal(0.2, 0.87, 1, 0.04, self.logofont_, self.textfont_, self.atlasinfo_)
            else:
                SLACJetWorkingGroup(0.2, 0.87,  46, 0.04, self.logofont_, self.textfont_, )
        else:
            CustomLabel(0.2, 0.87, 1, 0.05, self.logofont_, self.textfont_, self.customlabel_)

        for i,label in zip(range(len(labels)), labels):
		    myText(0.2, 0.82-i*0.04, label, 0.03, 1, self.textfont_);

            
        if(self.drawRatio_):
            canv.UseRatioPad()
            if not self.drawRatioUncert_ == "":
                self.H_[self.drawRatioUncert_].h_.Draw(self.H_[self.drawRatioUncert_].drawstyle_)
                l=ROOT.TLine(self.H_[self.drawRatioUncert_].h_.GetXaxis().GetBinLowEdge(self.H_[self.drawRatioUncert_].h_.GetXaxis().GetFirst()), 
                                 1.00, 
                                 self.H_[self.drawRatioUncert_].h_.GetXaxis().GetBinUpEdge(self.H_[self.drawRatioUncert_].h_.GetXaxis().GetLast()), 
                                 1.00)
                l.SetLineWidth(2);
                l.SetLineStyle(7);
                l.Draw()
            for index, ID in zip(range(len(self.HlistRatio_)),self.HlistRatio_):
                if index==0 and self.drawRatioUncert_ == "":
                    self.H_[ID].h_.Draw(self.H_[ID].drawstyle_)

                    l=ROOT.TLine(self.H_[ID].h_.GetXaxis().GetBinLowEdge(self.H_[ID].h_.GetXaxis().GetFirst()), 
                                 1.00, 
                                 self.H_[ID].h_.GetXaxis().GetBinUpEdge(self.H_[ID].h_.GetXaxis().GetLast()), 
                                 1.00)
                    l.SetLineWidth(2);
                    l.SetLineStyle(7);
                    l.Draw()
                else :
                    self.H_[ID].h_.Draw(self.H_[ID].drawstyle_+"same")

            self.drawRatio_  = False # prevent next histo to have a ratio plot per default
            self.HlistRatio_ = []


        save(canv, self.plotsavename_)
        if(self.saveHistos_):
            for hist in hlist:
                SaveHist(hist)

        return [canv, leg]

            

    def MakeDiscriminantFromXML(self, xmlfile, TMVAMethodName, LikelihoodHistoName, varlistnames, spectlistnames, LikelihoodHistoNbins, LikelihoodHistoAxisMins, LikelihoodHistoAxisMaxs, IDnew):
        self.Set2DCanvas(True,True)
        if(self.atlasStyle_):
            self.SetATLASStyle()

        reader = TMVA.Reader( "!Color:!Silent" )
        # Vars and Spectators 
        varlist   =[]
        spectlist =[]
        for var, varindex in zip(varlistnames,range(len(varlistnames))):
            varlist.append(array.array('f',[0]))
            reader.AddVariable(var, varlist[varindex])

        for spec, specindex in zip(spectlistnames,range(len(spectlistnames))):
            spectlist.append(array.array('f',[0]))
            reader.AddSpectator(spec, spectlist[specindex])

        reader.BookMVA(TMVAMethodName, xmlfile)

        histoLikelihood = ROOT.TH2F(LikelihoodHistoName,"",LikelihoodHistoNbins[0], LikelihoodHistoAxisMins[0], LikelihoodHistoAxisMaxs[0], LikelihoodHistoNbins[1], LikelihoodHistoAxisMins[1], LikelihoodHistoAxisMaxs[1])
        for i in range(1,histoLikelihood.GetNbinsX() + 1):
            for j in range(1,histoLikelihood.GetNbinsY() + 1):
                # find the bin center coordinates
                (varlist   [0])[0] = histoLikelihood.GetXaxis().GetBinCenter(i)
                (varlist   [1])[0] = histoLikelihood.GetYaxis().GetBinCenter(j)
                if   ( (varlist[0])[0] == 0):     Likelihood = 0
                else:                             Likelihood = reader.EvaluateMVA(TMVAMethodName)
#                print varlistnames[0], (varlist[0])[0], varlistnames[1], (varlist[1])[0], TMVAMethodName, Likelihood
                histoLikelihood.SetBinContent(i,j,Likelihood)

        H=Hist("", "")
        H.sethist(histoLikelihood)
        self.H_[IDnew]=H
        self.H_[IDnew].setdrawstyle("colz")
        self.H_[IDnew].AddToLeg(False)

    
    def AddDiscriminantToBranches(self, fFileIn, tTree, hFile, hName, xvar, yvar, newVarname):

        self.GetHisto(hName, newVarname, hFile)


        Tfile = ROOT.TFile(fFileIn, "update")
        Tfile.cd()
        tree=ROOT.gDirectory.Get(tTree)
        if tree == None:
            raise NameError('cannot find TTree '+tTree)

        newVar = array.array('f',[0])
        Xvar   = array.array('f',[0])
        Yvar   = array.array('f',[0])
        newBranch = tree.Branch(newVarname, newVar, newVarname+"/F")
        tree.SetBranchAddress(xvar,Xvar);
        tree.SetBranchAddress(yvar,Yvar);

        nentries = tree.GetEntries();
        for entry in range(nentries):
            tree.GetEntry(entry)

            binnr      = self.H_[newVarname].h_.FindBin(Xvar[0], Yvar[0])
            newVar[0]  = self.H_[newVarname].h_.GetBinContent(binnr)
#            print xvar, Xvar[0], yvar, Yvar[0], newVar[0]  
            newBranch.Fill()
        tree.Print()
        tree.Write()



    
    def GetFractionWithinBoundsGraph(self, histos, leg,  newID, lower, upper):
        nbins = len(histos)
        fraction  = n.zeros(nbins, dtype=float)
        fractione = n.zeros(nbins, dtype=float)
        trimcut   = n.zeros(nbins, dtype=float)
        trimcute  = n.zeros(nbins, dtype=float)
        for index, h in zip(range(len(histos)),histos):
            low = self.H_[h].h_.FindBin(lower)
            up  = self.H_[h].h_.FindBin(upper)
            print index
            fraction[index] = self.H_[h].h_.Integral(low, up)/self.H_[h].h_.Integral()
            trimcut[index]  = leg[index]

        tgraph = ROOT.TGraphErrors(nbins, trimcut, fraction, trimcute, fractione)
        H=Hist("", "")
        H.sethist(tgraph)
        self.H_[newID]=H

    def GetSigmaOverMeanGraph(self, histos, leg,  newID, lower, upper):
        nbins = len(histos)
        sigmaovermean  = n.zeros(nbins, dtype=float)
        sigmaovermeane = n.zeros(nbins, dtype=float)
        trimcut        = n.zeros(nbins, dtype=float)
        trimcute       = n.zeros(nbins, dtype=float)
        for index, h in zip(range(len(histos)),histos):
            ibin = self.H_[h].h_.GetMaximumBin()
            f1 = ROOT.TF1('myf1',Gauss(),0, 200 ,3)
            f1.SetParameter(0,100);
            f1.SetParameter(1,90);
            f1.SetParameter(2,10);
            self.H_[h].h_.Fit("myf1","R","",self.H_[h].h_.GetBinCenter(ibin)-20, self.H_[h].h_.GetBinCenter(ibin)+20)
#            self.H_[h].h_.Draw()
#            f1.Draw("same")
            par = f1.GetParameters()
            print par[0], par[1], par[2]
#            wait()
            sigmaovermean[index] = par[2]/par[1]
            trimcut[index]  = leg[index]

        tgraph = ROOT.TGraphErrors(nbins, trimcut, sigmaovermean, trimcute, sigmaovermeane)
        H=Hist("", "")
        H.sethist(tgraph)
        self.H_[newID]=H
        
    def GetPeakPlusMinusRangeGraph(self, histos, leg,  newID, plusminusrange):
        nbins = len(histos)
        fraction  = n.zeros(nbins, dtype=float)
        fractione = n.zeros(nbins, dtype=float)
        trimcut   = n.zeros(nbins, dtype=float)
        trimcute  = n.zeros(nbins, dtype=float)
        for index, h in zip(range(len(histos)),histos):
            ibin = self.H_[h].h_.GetMaximumBin()
            center = self.H_[h].h_.GetBinCenter(ibin)
            low  = self.H_[h].h_.FindBin(center-plusminusrange)
            up   = self.H_[h].h_.FindBin(center+plusminusrange)
            print index
            fraction[index] = self.H_[h].h_.Integral(low, up)/self.H_[h].h_.Integral()
            trimcut[index]  = leg[index]

        tgraph = ROOT.TGraphErrors(nbins, trimcut, fraction, trimcute, fractione)
        H=Hist("", "")
        H.sethist(tgraph)
        self.H_[newID]=H

    # COLORSSolid  = [ROOT.kViolet+9,ROOT.kTeal-5,ROOT.kMagenta+3, ROOT.kBlack]
    # COLORSOpen   = [ROOT.kViolet+9,ROOT.kTeal-5,ROOT.kMagenta+3, ROOT.kBlack]
    COLORSSolid  = [ROOT.kMagenta+3, ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    COLORSOpen   = [ROOT.kMagenta+3, ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    COLORSFilled = [ROOT.kTeal-5,ROOT.kViolet+9, ROOT.kMagenta+3 ]
    COLORSLine   = [ROOT.kGray+3, ROOT.kRed+1, ROOT.kMagenta+3,ROOT.kTeal-5, ROOT.kBlack,  ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1,]
    MARKERSSolid = [20, 21, 22, 23, 33, 34, 28, 29] 
    MARKERSOpen  = [24, 25, 26, 32, 27, 28, 30]
    STYLESFilled = [3002, 3004 , 3002]
    MARKERSData  = [20, 4] 


class Linear:
    def __call__(self, x, par ):
        return par[0] + x[0]*par[1]

class Const:
    def __call__(self, x, par ):
        return par[0]
    
class Gauss:
    def __call__(self, x, par):
        return par[0]*ROOT.TMath.Exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]))

class Hist:
    def __init__(self, histname, histfile):
        self.histname_       = histname
        self.histfile_       = histfile
        self.tree_           = 0
        self.tfile_          = 0
        self.h_              = 0
        self.drawstyle_      = "EX"
        self.addToLeg_       = True
        self.solid_          = True
        self.filled_         = False
        self.filledWithLine_ = False
        self.line_           = False
        self.isData_         = False
        self.normalized_     = False
        self.scale_          = 1
        self.styleset_       = False
        self.legstyle_       = ""

    def gethist(self):
        self.tfile_ = ROOT.TFile(self.histfile_)
        self.tfile_.cd()
        self.h_       = ROOT.gDirectory.Get(self.histname_).Clone(self.histname_)

    def gethistfromtree(self, tree, expression, selection, nbins, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax):
        self.tfile_ = ROOT.TFile(self.histfile_)
        self.tfile_.cd()
        self.tree_  = ROOT.gDirectory.Get(tree)
        if self.tree_ == None:
            raise NameError('cannot find TTree '+tree)

        print ">>>> plottingBase: getting histo with sel:", selection
        name = expression.replace("(", "").replace(")", "").replace("/","").replace(":","").replace("[","").replace("]","").replace("+","")
        tree = tree.replace("/","_")
        if(nbinsz !=0):
            self.h_ = ROOT.TH3D("h_%s_%s"%(tree,name) ,"",nbins, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax)
        elif(nbinsy !=0):
            self.h_ = ROOT.TH2D("h_%s_%s"%(tree,name) ,"",nbins, xmin, xmax, nbinsy, ymin, ymax)
        else:
            self.h_ = ROOT.TH1D("h_%s_%s"%(tree,name) ,"",nbins, xmin, xmax)
        self.h_.Sumw2()
        self.tree_.Draw(expression+">>h_%s_%s"%(tree,name),selection,"goff")

    def sethist(self, hist):
        self.h_       = hist
    def setopen(self):
        self.solid_   = False
    def setfilled(self, yesno=True):
        self.filled_  = yesno
    def setfilledwithline(self):
        self.filledWithLine_ = True
        self.filled_         = True
        self.drawstyle_      = "hist"
    def setline(self):
        self.line_    = True
        self.setfilled(False)
        self.solid_   = False
    def setdrawstyle(self, style):
        self.drawstyle_=style
        if(self.drawstyle_ == "hist"):
            self.setfilled()
    def AddToLeg(self, yesno):
        self.addToLeg_=yesno
    def isData(self):
        self.isData_    = True
        self.drawstyle_ ="EX"
        self.solid_     = True
    def setScale(self, scale):
        self.scale_ = scale
        self.h_.Scale(scale)
    def setNormalized(self):
        self.normalized_ = True
    def setStyle(self, marker=20, color=1, style="hist", xtitle="", ytitle="", fillstyle=0):
        self.h_.SetMarkerStyle(marker)
        self.h_.SetMarkerColor(color)
        self.h_.SetLineColor  (color)
        self.h_.SetFillColor  (color)
        self.h_.SetFillStyle  (fillstyle)
        self.h_.SetXTitle     (xtitle)
        self.h_.SetYTitle     (ytitle)
        self.setdrawstyle     (style)


        


def wait():
    var = raw_input("click 'y' to save ")
    if var == 'y':
        savename =  raw_input("enter name to save plot. ")
        SaveCan(canv, savename)

def save(canv, name=None):
    if name == None:
        name = os.path.basename(sys.argv[0]).replace(".py", "")
    var = raw_input("click 'y' to save")
    if var == 'y':
        yes = raw_input("under name:  plots/" + name+"? (y/n)")
        if yes == 'y':
            SaveCan(canv, "plots/"+name)
        else:
            savename =  raw_input("enter name to save plot. ")
            SaveCan(canv, savename)


def SaveCan(can, name):
    try:
        print ">>> saving Canvas",can.GetName(),"as",name," .pdf/eps/pnd/root/C <<<"
        can.SaveAs(name+".pdf")
        can.SaveAs(name+".eps")
        can.SaveAs(name+".root")
        can.SaveAs(name+".C")
    except AttributeError:
        print ">>> saving Canvas", can.canv.GetName(),"as",name," .pdf/eps/pnd/root/C <<<"
        can.canv.SaveAs(name+".pdf")
        can.canv.SaveAs(name+".eps")
        can.canv.SaveAs(name+".root")
        can.canv.SaveAs(name+".C")

def SaveHist(h, name="hist.root"):
    file = ROOT.TFile(name.replace(".root","")+".root", "UPDATE")
    print file
    h.h_.Write()
    file.Close()
    

def ATLASInternal(x, y, color=1, textSize=0.04, logofont=32, font=132, label="Internal"):
	l=ROOT.TLatex()
 	l.SetNDC()
 	l.SetTextFont(logofont)
 	l.SetTextColor(color)
 	l.SetTextSize(textSize)
	l.SetLineWidth(1)
 	l.DrawLatex(x,y,"ATLAS")
 	myText(x+0.13, y, label, textSize, 1, font)

def SLACJetWorkingGroup(x, y, color=1, textSize=0.04, logofont=32, font=132, label="Jet Working Group"):
	l=ROOT.TLatex()
 	l.SetNDC()
 	l.SetTextFont(62)
 	l.SetTextColor(color)
 	l.SetTextSize(textSize)
	l.SetLineWidth(1)
 	l.DrawLatex(x,y,"SLAC")
 	myText(x+0.11, y, label, textSize, 46, 62)

def CustomLabel(x, y, color=1, textSize=0.04, logofont=32, font=132, label="mylabel"):
 	myText(x, y, label, textSize, color, font)

def myText(x,y,text, tsize=0.04,color=1, font=132):
  	l=ROOT.TLatex()
  	l.SetTextSize(tsize)
  	l.SetNDC()
 	l.SetTextFont(font)
  	l.SetTextColor(color)
  	l.DrawLatex(x,y,text)


def round_to_value(number,roundto):
    return (round(number / roundto) * roundto)



