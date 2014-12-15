#!/usr/bin/python
##################################################################
# pascal nef                             April 23rd, 2013        #
##################################################################
import ROOT as R
import array
import sys
import os
import numpy as n
from plottingBase import Plotter
from plottingBase import Hist
##################################################################

class DataToMCPlottingBase(Plotter):
    def __init__(self, ):
        Plotter.__init__(self, MultiProcess=False) # MultiProcess=True does not yet work
        self.doNormMCtoData_   = True
        self.MCscaleFact_      = 1
        self.drawMCcomponents_ = True
        self.doRatio_          = False
        self.printInfo_        = True
        self.drawMCsum_        = True
        self.exclude_          = [] # (hist, name) tuples
        self.linestyle_        = []

        self.COLORSSolid  = [R.kViolet+9,R.kMagenta+3, R.kTeal-5,   R.kBlack, R.kGray+3, R.kGray+2, R.kGray+1, R.kOrange+10]
        self.COLORSOpen   = [R.kViolet+9,R.kMagenta+3, R.kTeal-5,   R.kBlack, R.kGray+3, R.kGray+2, R.kGray+1, R.kOrange+10]
        self.COLORSFilled = [R.kViolet+9,  R.kTeal-5,  R.kViolet+9, R.kMagenta+3 ]
        self.STYLESFilled = [3004, 3002, 3004 , 3002]

        # for fitting
        self.fit_ = {}

    def MLfit(self, IDsig, IDbkg, IDdata, varranges, fixHSscale=0):
        self.fit_ = {}
        self.fit_['var']   = R.RooRealVar("var", "var", varranges[0], varranges[1])

        self.fit_['DataH'] = R.RooDataHist ("data"   ,"data"  , R.RooArgList(self.fit_['var']), self.H_[IDdata].h_)
        self.fit_['SigH']  = R.RooDataHist ("Sig"    ,"Sig"   , R.RooArgList(self.fit_['var']), self.H_[IDsig].h_) 
        self.fit_['BkgH']  = R.RooDataHist ("Bkg"    ,"Bkg"   , R.RooArgList(self.fit_['var']), self.H_[IDbkg].h_) 

        self.fit_['Sig_pdf']= R.RooHistPdf     ("Sig_pdf", "Sig_pdf", R.RooArgSet(self.fit_['var']), self.fit_['SigH']) 
        self.fit_['Bkg_pdf']= R.RooHistPdf     ("Bkg_pdf" ,"Bkg_pdf", R.RooArgSet(self.fit_['var']), self.fit_['BkgH']); 
	    
        if fixHSscale >0:
            self.fit_['nBkg']   = R.RooRealVar  ("nBkg"   ,"number of PU jets ",    self.H_[IDdata].h_.Integral(), 0 ,self.H_[IDdata].h_.Integral())
            self.fit_['nSig']   = R.RooRealVar  ("nSig"   ,"number of HS jets ",    self.H_[IDsig].h_.Integral()*fixHSscale);
            self.fit_['nSig'].setConstant(R.kTRUE)
        else:
            self.fit_['nSig']   = R.RooRealVar  ("nSig"   ,"number of HS jets ",    self.H_[IDdata].h_.Integral(), 0 ,self.H_[IDdata].h_.Integral());
            self.fit_['nBkg']   = R.RooRealVar  ("nBkg"   ,"number of PU jets ",    self.H_[IDdata].h_.Integral(), 0 ,self.H_[IDdata].h_.Integral())

        model = R.RooAddPdf ("model","model", R.RooArgList(self.fit_['Sig_pdf'],self.fit_['Bkg_pdf']), R.RooArgList(self.fit_['nSig'], self.fit_['nBkg']))
        model.fitTo(self.fit_['DataH'],R.RooFit.SumW2Error(R.kFALSE),R.RooFit.Extended()); 
        
        print ""
        print ""
        print ""
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        print "---- scale factors ---- "
        print "nBkg", self.fit_['nBkg'].getVal(), "pm", self.fit_['nBkg'].getError(), "bkg scale factor", self.fit_['nBkg'].getVal()/self.H_[IDbkg].h_.Integral()
        print "nSig", self.fit_['nSig'].getVal(), "pm", self.fit_['nBkg'].getError(), "sig scale factor", self.fit_['nSig'].getVal()/self.H_[IDsig].h_.Integral()
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        print ""
        print ""
        print ""
        self.fit_['sigScaleFact'] = self.fit_['nSig'].getVal()/self.H_[IDsig].h_.Integral()
        self.fit_['bkgScaleFact'] = self.fit_['nBkg'].getVal()/self.H_[IDbkg].h_.Integral()


        


    def DataToMC(self, HMClist, HData, xtitle="x", ytitle="y", labels=[], legends=[]):
        print "                             "
        print "                             "
        print "                             "
        print "                             "
        print ">>> DataToMCPlottingBase: new plot <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" 
        print ">>> DataToMCPlottingBase: Data", HData
        print ">>> DataToMCPlottingBase: MC  ", HMClist  


        self.AddHistos(HMClist,"MCsum")
        self.H_[HData].isData()

        if self.doNormMCtoData_:
            scale = self.H_[HData].h_.Integral()/self.H_["MCsum"].h_.Integral()
            print ">>> DataToMCPlottingBase: Normalizing MC to Data: Scale Factor =", scale  
            self.SetMCScale(scale,HMClist+["MCsum"])

        # fix styles ---------------------------
        self.H_["MCsum"].setdrawstyle("hist")


        # decide which histos to draw
        if self.drawMCcomponents_:
            HList  = ["MCsum"]+HMClist+[HData]
            HNames = ["MC incl."]+legends+["Data"]
        else:
            HList  = ["MCsum"]+[HData]
            HNames = ["MC incl.","Data"]
        if not self.drawMCsum_:
            HList .remove("MCsum")
            HNames.remove("MC incl.")


        for (name, hist) in self.exclude_:
            try:
                HNames.remove(name)
            except ValueError:
                print name, "not in", HNames
            try:
                HList.remove(hist)
            except ValueError:
                print name, "not in", HNames

        if self.doRatio_:
            self.DrawRatioMCUncert("MCsum", "Data / MC")
            self.DrawRatio(HData, "MCsum", "MCsumRatio")
            self.H_["MCsumRatio"].setStyle(20,1,"EX","", "Data / MC")

        # info
        if self.printInfo_:
            for hmc,name in zip(HList,HNames):
                print '{0:25}: Nevents {1:10}: {2:10.1f} -> Fraction wrt data: {3}'.format(">>> DataToMCPlottingBase",name,self.H_[hmc].h_.Integral(),repr(self.H_[hmc].h_.Integral()/self.H_[HData].h_.Integral()))

        # draw 
        print ">>> DataToMCPlottingBase: draw plot"; print; print
        self.Draw(HList, xtitle, ytitle, False, True, labels, HNames )



        


    def DataToMCSyst(self,  HMCSyst, HMC, HData, xtitle="x", ytitle="y", labels=[], legends=[], AddHlist=[]):
        print "                             "
        print "                             "
        print "                             "
        print "                             "
        print ">>> DataToMCPlottingBase: new plot <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" 
        print ">>> DataToMCPlottingBase: Data", HData
        print ">>> DataToMCPlottingBase: MC  ", HMC
        print ">>> DataToMCPlottingBase: MC  ", HMCSyst  


        if self.doRatio_:
            self.DrawRatioMCUncert(HMCSyst, "Data / MC", True)
            self.DrawRatio(HData, HMCSyst, "MCsumRatio")
            self.H_["MCsumRatio"].setStyle(20,1,"EX","", "Data / MC")


        # draw 
        print ">>> DataToMCPlottingBase: draw plot"; print; print
        self.H_[HMCSyst].setStyle(1,920,"E2","", "" ,1001)
        self.H_[HMCSyst].h_.SetFillColor(920)
        self.H_[HMCSyst].h_.SetLineColor(920)
        self.H_[HMCSyst].h_.SetMarkerColor(920)
        self.H_[HMCSyst].styleset_ = True
        self.H_[HMCSyst].legstyle_ = "fl"


        self.Draw([HMCSyst,HMC,HData]+AddHlist, xtitle, ytitle, False, True, labels, legends )



