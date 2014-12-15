#!/usr/bin/env python
##################################################################
# pascal nef                             April 23rd, 2013        #
##################################################################
import ROOT as R
import sys
import os
import re
import math
import time
from optparse import OptionParser
from sys import argv,exit
from DataToMCPlottingBase import *

##################################################################

def testEventSelection():
    # dPhi distribution
    plotter.GetHistoFromTree("JVFData",      tfileData,  tree ,"fabs(J0ZDPhi)" , "Weight*(1                                        &&"+JetSelection+"&&ZPt>20&&ZPt<40"+")", 30, 0, 3.15)
    plotter.GetHistoFromTree("JVFMC_HS",     tfileMC  ,  tree ,"fabs(J0ZDPhi)" , "Weight*(JisHS&&Jtruthpt>10                       &&"+JetSelection+"&&ZPt>20&&ZPt<40"+")", 30, 0, 3.15)
    plotter.GetHistoFromTree("JVFMC_nHS",    tfileMC  ,  tree ,"fabs(J0ZDPhi)" , "Weight*(!(JisHS&&Jtruthpt>10)                    &&"+JetSelection+"&&ZPt>20&&ZPt<40"+")", 30, 0, 3.15)
    plotter.GetHistoFromTree("JVFMC_PU",     tfileMC  ,  tree ,"fabs(J0ZDPhi)" , "Weight*(JdRTruthPt4>0.4                          &&"+JetSelection+"&&ZPt>20&&ZPt<40"+")", 30, 0, 3.15)
    plotter.GetHistoFromTree("JVFMC_other",  tfileMC  ,  tree ,"fabs(J0ZDPhi)" , "Weight*(!(JdRTruthPt4>0.4)&&!(JisHS&&Jtruthpt>10)&&"+JetSelection+"&&ZPt>20&&ZPt<40"+")", 30, 0, 3.15)
    plotter.SetMCScale(ScaleFact,    ["JVFMC_HS","JVFMC_nHS","JVFMC_PU","JVFMC_other"])
    plotter.SetAtlasInfo("Internal")
    plotter.DataToMC(["JVFMC_HS","JVFMC_PU", "JVFMC_other" ], "JVFData" , "|#Delta#phi(Z,jet)|", "Events", [lumiLabel,MClabel,jetlabel,"20<p_{T}<50 GeV, |#eta|<2.4", "30 < Z p_{T} < 40 GeV" ], [ "HS", "PU" , "other"])
    
    # dPhi distribution JVT > 0.4
    plotter.GetHistoFromTree("JVFData",      tfileData,  tree ,"fabs(J0ZDPhi)" , "Weight*(1                                        &&"+JetSelection+"&&ZPt>30&&ZPt<40&&kNN100trim_pt20to50_Likelihood>0.4"+")", 30, 0, 3.15)
    plotter.GetHistoFromTree("JVFMC_HS",     tfileMC  ,  tree ,"fabs(J0ZDPhi)" , "Weight*(JisHS&&Jtruthpt>10                       &&"+JetSelection+"&&ZPt>30&&ZPt<40&&kNN100trim_pt20to50_Likelihood>0.4"+")", 30, 0, 3.15)
    plotter.GetHistoFromTree("JVFMC_nHS",    tfileMC  ,  tree ,"fabs(J0ZDPhi)" , "Weight*(!(JisHS&&Jtruthpt>10)                    &&"+JetSelection+"&&ZPt>30&&ZPt<40&&kNN100trim_pt20to50_Likelihood>0.4"+")", 30, 0, 3.15)
    plotter.GetHistoFromTree("JVFMC_PU",     tfileMC  ,  tree ,"fabs(J0ZDPhi)" , "Weight*(JdRTruthPt4>0.4                          &&"+JetSelection+"&&ZPt>30&&ZPt<40&&kNN100trim_pt20to50_Likelihood>0.4"+")", 30, 0, 3.15)
    plotter.GetHistoFromTree("JVFMC_other",  tfileMC  ,  tree ,"fabs(J0ZDPhi)" , "Weight*(!(JdRTruthPt4>0.4)&&!(JisHS&&Jtruthpt>10)&&"+JetSelection+"&&ZPt>30&&ZPt<40&&kNN100trim_pt20to50_Likelihood>0.4"+")", 30, 0, 3.15)
    plotter.SetMCScale(ScaleFact,    ["JVFMC_HS","JVFMC_nHS","JVFMC_PU","JVFMC_other"])
    plotter.SetAtlasInfo("Internal")
    plotter.DataToMC(["JVFMC_HS","JVFMC_PU", "JVFMC_other" ], "JVFData" , "|#Delta#phi(Z,jet)|", "Events", [lumiLabel,MClabel,jetlabel,"20<p_{T}<50 GeV, |#eta|<2.4", "30 < Z p_{T} < 40 GeV", "JVT>0.4" ], [ "HS", "PU" , "other"])

    
    # JVT distribution
    plotter.GetHistoFromTree("JVFData",      tfileData,  tree ,"kNN100trim_pt20to50_Likelihood" , "Weight*(1                                        &&"+JetSelection+" && "+dPhiVar+">"+str(Signal_dPhiMin)+")", 30, -0.2,1.1)
    plotter.GetHistoFromTree("JVFMC_HS",     tfileMC  ,  tree ,"kNN100trim_pt20to50_Likelihood" , "Weight*(JisHS&&Jtruthpt>10                       &&"+JetSelection+" && "+dPhiVar+">"+str(Signal_dPhiMin)+")", 30, -0.2,1.1)
    plotter.GetHistoFromTree("JVFMC_nHS",    tfileMC  ,  tree ,"kNN100trim_pt20to50_Likelihood" , "Weight*(!(JisHS&&Jtruthpt>10 )                   &&"+JetSelection+" && "+dPhiVar+">"+str(Signal_dPhiMin)+")", 30, -0.2,1.1)
    plotter.GetHistoFromTree("JVFMC_PU",     tfileMC  ,  tree ,"kNN100trim_pt20to50_Likelihood" , "Weight*(JdRTruthPt4>0.4                          &&"+JetSelection+" && "+dPhiVar+">"+str(Signal_dPhiMin)+")", 30, -0.2,1.1)
    plotter.GetHistoFromTree("JVFMC_other",  tfileMC  ,  tree ,"kNN100trim_pt20to50_Likelihood" , "Weight*(!(JdRTruthPt4>0.4)&&!(JisHS&&Jtruthpt>10)&&"+JetSelection+" && "+dPhiVar+">"+str(Signal_dPhiMin)+")", 30, -0.2,1.1)

    plotter.SetMCScale(ScaleFact,    ["JVFMC_HS","JVFMC_nHS","JVFMC_PU","JVFMC_other"])

    print ">>> Scale Factor Ok?? <<< "
    plotter.DataToMC(["JVFMC_HS","JVFMC_PU", "JVFMC_other" ], "JVFData" , "JVT", "Events", [lumiLabel,MClabel,jetlabel,"20<p_{T}<50 GeV, |#eta|<2.4"], [ "HS", "PU" , "other"])




def PrintWorkingPoint():
    plotter.GetHistoFromTree("JVTMC_HS",     tfileMC,tree ,"kNN100trim_pt20to50_Likelihood" ,"Weight*("+HSDef              +" && "+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+")", 1000, -0.,1.1)
    plotter.GetHistoFromTree("JVTMC_PU",     tfileMC,tree ,"kNN100trim_pt20to50_Likelihood" ,"Weight*("+PUDef              +" && "+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+")", 1000, -0.,1.1)
    plotter.GetHistoFromTree("JVTMC_HStight",tfileMC,tree ,"kNN100trim_pt20to50_Likelihood" ,"Weight*(JisHS && Jtruthpt>10    && "+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+")", 1000, -0.,1.1)

    try:
        cut = CUTVal[0]
    except:
        cut = CUTVal
    mybin = plotter.H_["JVTMC_HS"].h_.FindBin(cut)
    print ">>>>>>>>>>>>>>>>>>>>> Efficiencies for WP <<<<<<<<<<<<<<<<<<<<<<<<<<<"
    print ">>> PU def  ", "fake rate", plotter.H_["JVTMC_PU"].h_.Integral(mybin, -1)     /plotter.H_["JVTMC_PU"].h_.Integral()
    print ">>> HS def  ", "eff      ", plotter.H_["JVTMC_HS"].h_.Integral(mybin, -1)     /plotter.H_["JVTMC_HS"].h_.Integral()
    print ">>> HS tight", "eff      ", plotter.H_["JVTMC_HStight"].h_.Integral(mybin, -1)/plotter.H_["JVTMC_HStight"].h_.Integral()
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<"



def PUJetDefFlatness():
    from plottingBase import Const
    from plottingBase import wait
    fPU1 = R.TF1('myf_PU1',Const(),0, 3.1416 ,1)
    fPU2 = R.TF1('myf_PU2',Const(),0, 3.1416 ,1)
    

    plotter.GetHistoFromTree("JVFMC_PU1",     tfileMC_pow  ,  tree ,"fabs(J0ZDPhi)" , "Weight*(JdRTruthPt4>0.4                          &&"+JetSelection+"&&ZPt>30 && ZPt<40)", 10, 0, 3.15)
    plotter.GetHistoFromTree("JVFMC_PU2",     tfileMC_pow  ,  tree ,"fabs(J0ZDPhi)" , "Weight*(JisPU                                    &&"+JetSelection+"&&ZPt>30 && ZPt<40)", 10, 0, 3.15)

    plotter.SetAtlasInfo("Internal")
    plotter.SetMCScale(ScaleFact,    ["JVFMC_PU1","JVFMC_PU2"])
    plotter.H_["JVFMC_PU1"].h_.Fit("myf_PU1","R","",0,2)
    plotter.H_["JVFMC_PU2"].h_.Fit("myf_PU2","R","",0,2)
    plotter.Draw(["JVFMC_PU1","JVFMC_PU2"], "|#Delta#phi(Z,jet)|", "Entries", False, True, [lumiLabel,MClabel,jetlabel,"20<p_{T}<50 GeV, |#eta|<2.4", "20 < Z p_{T} < 30 GeV"],["PU (#DeltaR>0.4)", "PU (#DeltaR>0.6)"])
    fPU1.SetRange(0,3.1415)
    fPU1.Draw()
    fPU2.Draw()
    wait()


# Main ------------------------------------------------------------------------
if __name__ == "__main__":
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # parse arguments --------------------------------------------
    usage = """
    usage: %prog [options] InputFile.root
    	
    """
    parser = OptionParser(usage)
    
    parser.add_option("--verbose"            ,dest="verbose"           ,    help="set verbose"          , action="store_true")
    parser.add_option("--WP"                 ,dest="WP"                ,    help="set verbose"          , action="store_true")
    parser.add_option("--EventSel"           ,dest="EventSel"          ,    help="set verbose"          , action="store_true")
    parser.add_option("--flatness"           ,dest="flatness"          ,    help="set verbose"          , action="store_true")
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    global tfileData, tfileMC, tree, jetlabel,ptlabel,lumiLabel, MClabel, ScaleFact, JetSelection, plotter

    path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/"
    tfileData    = path+"2014021x.PileUpStudies_DiLeptonSelection_rev351830.MuonsSMWZPeriodABCDEGHIJL.jetmet2012pileupcustom.root"
    tfileMC_pow  = path+"20140214.22.39_PileUpStudies_DiLeptonSelection_rev351830.Zmumu_PowhegPythia8_MC12_COMMON.jetmet2012pileupcustom.root"
    tfileMC      = path+"20140305.08.12_PileUpStudies_DiLeptonSelection_351830.Zmumu_MC12_Sherpa_COMMON.jetmet2012pileupcustom.root"
    tree         = "JetTree0"

    # labels 
    jetlabel = "Anti-k_{t} LCW+JES R=0.4";
    ptlabel  = "20 < p_{T} < 50 GeV, |#eta| < 2.4";
    lumiLabel= "#sqrt{s} = 8 TeV, L = 20.3 fb^{-1}"
    MClabel = "Sherpa Z#rightarrow#mu#mu"
    Zsel     = "20 < p_{T}^{ref} < 30 GeV"

    # parameters ---------------------------
    ScaleFact        = 20.3*1.2*1000*1.2 # additional 20% scale for Sherpa
    JetSelection     = "ZPt>20 && Jpt>20 && Jpt<50 && fabs(Jeta)<2.4 && Jindex==0 &&ZPt>20 && ZPt<30" 
    Variable         = "Jeta"
    VariableSyst     = "fabs(Jeta)"
    VarName          = "jet #eta"
    VariableNBins    = 10
    VariableLow      = -2.4
    VariableHigh     = 2.4
    VariableNBinsSyst= 5
    VariableLowSyst  = 0
    VariableHighSyst = 2.4
    CUT              = "kNN100trim_pt20to50_Likelihood >"
    CUTVal           = [0.4, 0.7, 0.2, 0.4,0.7]
    dPhiVar          = "fabs(J0ZDPhi)"
    dPhiNBins        = 300
    dPhiLow          = 0
    dPhiHigh         = 3.15
    PUDef            = "JdRTruthPt4>0.4"
    HSDef            = "!(JdRTruthPt4>0.4)"
    TrueHSDef        = "Jtruthpt>10 && JisHS"
    #---
    PUControl_dPhiMax = 1.2
    Signal_dPhiMin    = 2.8
    Signal_dPhiMax    = R.TMath.Pi()
    HSuncert          = 0.3  # uncertainty on HS contribution
    doDataErr         = False
    NSmoothSyst       = 1 # 1 for pt 40 to 50 JVT > 0.7
    # parameters ---------------------------





    #++++++++++++++++ BASIC SETUP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = DataToMCPlottingBase()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Preliminary")
    plotter.doRatio_ = True
    plotter.printInfo_ = True
    plotter.doNormMCtoData_ = False
    plotter.drawMCcomponents_ = True
    plotter.exclude_ = [("other","JVFMC_other"), ("EffRecoCorrSyst", "EffRecoCorrSyst")]


    # -----------------------------------------------------------------------------------------
    if(options.flatness):
        PUJetDefFlatness()
    if options.WP:
        PrintWorkingPoint()
    if options.EventSel:
        testEventSelection()

       
    #--------------------------------------------------------------------------------------------
    # -------------------- PU bkg estimate ---------------------------------------------------
    plotter.GetHistoFromTree("dPhi_nHS",    tfileMC  ,  tree ,dPhiVar+":"+Variable , "Weight*("+PUDef+"&&"+JetSelection+")", VariableNBins, VariableLow, VariableHigh, dPhiNBins,  dPhiLow, dPhiHigh)
    plotter.GetHistoFromTree("dPhi_HS",     tfileMC  ,  tree ,dPhiVar+":"+Variable , "Weight*("+HSDef+"&&"+JetSelection+")", VariableNBins, VariableLow, VariableHigh, dPhiNBins,  dPhiLow, dPhiHigh)
    plotter.GetHistoFromTree("dPhi_Inc",    tfileMC  ,  tree ,dPhiVar+":"+Variable , "Weight*(   1     &&"+JetSelection+")", VariableNBins, VariableLow, VariableHigh, dPhiNBins,  dPhiLow, dPhiHigh)
    plotter.SetMCScale(ScaleFact,    ["dPhi_nHS","dPhi_HS","dPhi_Inc"])
    
    nbins        = plotter.H_["dPhi_nHS"].h_.GetNbinsX()
    eta          = n.zeros(nbins, dtype=float)
    etaerr       = n.zeros(nbins, dtype=float)
    nBKG         = n.zeros(nbins, dtype=float)
    nBKGerr      = n.zeros(nbins, dtype=float)
    nBKGtruth    = n.zeros(nbins, dtype=float)
    nBKGerrtruth = n.zeros(nbins, dtype=float)
    tmp          = R.Double(0)

    for ibin in range(1,plotter.H_["dPhi_nHS"].h_.GetNbinsX()+1):
        binlow  = plotter.H_["dPhi_nHS"].h_.GetXaxis().GetBinLowEdge(ibin)
        binhigh = plotter.H_["dPhi_nHS"].h_.GetXaxis().GetBinLowEdge(ibin)+plotter.H_["dPhi_nHS"].h_.GetXaxis().GetBinWidth(ibin)
        plotter.SetRange("dPhi_nHS", "x", binlow, binhigh, "htmp")
        plotter.SetRange("dPhi_HS",  "x", binlow, binhigh, "htmpHS")
        plotter.SetRange("dPhi_Inc",  "x", binlow, binhigh, "htmpInc")
        plotter.Projector("htmp", "y", "htmp")
        plotter.Projector("htmpHS", "y", "htmpHS")
        plotter.Projector("htmpInc", "y", "htmpInc")
        bindPhiLimitBkg    = plotter.H_["htmp"].h_.FindBin(PUControl_dPhiMax)
        bindPhiLimitSig    = plotter.H_["htmp"].h_.FindBin(Signal_dPhiMin)
        bindPhiLimitSigMax = plotter.H_["htmp"].h_.FindBin(Signal_dPhiMax)
        eta [ibin-1]         = plotter.H_["dPhi_nHS"].h_.GetXaxis().GetBinCenter(ibin) 
        nBKG[ibin-1]         = (plotter.H_["htmpInc"].h_.IntegralAndError(0, bindPhiLimitBkg,tmp) - plotter.H_["htmpHS"].h_.Integral(0, bindPhiLimitBkg) ) * (Signal_dPhiMax - Signal_dPhiMin)/PUControl_dPhiMax
        nBKGerr[ibin-1]      = (Signal_dPhiMax- Signal_dPhiMin)/PUControl_dPhiMax* math.sqrt(tmp*tmp + (HSuncert*plotter.H_["htmpHS"].h_.Integral(0, bindPhiLimitBkg))*(HSuncert*plotter.H_["htmpHS"].h_.Integral(0, bindPhiLimitBkg)))
        nBKGtruth[ibin-1]        = plotter.H_["htmp"].h_.IntegralAndError(bindPhiLimitSig,bindPhiLimitSigMax,tmp) 
        nBKGerrtruth[ibin-1]     = tmp

    tgraph  = R.TGraphErrors(nbins, eta, nBKG,      etaerr , nBKGerr )
    tgraph2 = R.TGraphErrors(nbins, eta, nBKGtruth, etaerr , nBKGerrtruth )
    H=Hist("", "")
    H.sethist(tgraph)
    plotter.H_["PUcontamination"]=H

    H=Hist("", "")
    H.sethist(tgraph2)
    plotter.H_["PUcontaminationTruth"]=H

    plotter.H_["PUcontamination"]      .setdrawstyle("AP")
    plotter.H_["PUcontaminationTruth"].setdrawstyle("P")
    plotter.SetAtlasInfo("Simulation Internal")
    plotter.Draw(["PUcontamination", "PUcontaminationTruth"], VarName, "N_{j}^{PU} in signal region", False, True, [lumiLabel, MClabel, jetlabel, ptlabel,Zsel], ["MC prediction", "MC"])

    #-------------------------------------------------------------------------------------------------
    # -------------------- Data PU bkg estimate ---------------------------------------------------
    plotter.GetHistoFromTree("dPhi_data",   tfileData,  tree ,dPhiVar+":"+Variable , "Weight*(1        &&"+JetSelection+")", VariableNBins, VariableLow, VariableHigh, dPhiNBins,  dPhiLow, dPhiHigh)

    D_nbins        = plotter.H_["dPhi_data"].h_.GetNbinsX()
    D_eta          = n.zeros(nbins, dtype=float)
    D_etaerr       = n.zeros(nbins, dtype=float)
    D_nBKG         = n.zeros(nbins, dtype=float)
    D_nBKGerr      = n.zeros(nbins, dtype=float)
    D_tmp          = R.Double(0)
    


    for ibin in range(1,plotter.H_["dPhi_data"].h_.GetNbinsX()+1):
        binlow  = plotter.H_["dPhi_data"].h_.GetXaxis().GetBinLowEdge(ibin)
        binhigh = plotter.H_["dPhi_data"].h_.GetXaxis().GetBinLowEdge(ibin)+plotter.H_["dPhi_data"].h_.GetXaxis().GetBinWidth(ibin)
        plotter.SetRange("dPhi_data", "x", binlow, binhigh, "htmp")
        plotter.SetRange("dPhi_HS",  "x", binlow, binhigh, "htmpHS")
        plotter.Projector("htmp",   "y", "htmp")
        plotter.Projector("htmpHS", "y", "htmpHS")
        bindPhiLimitBkg = plotter.H_["htmp"].h_.FindBin(PUControl_dPhiMax)
        D_eta [ibin-1]         = plotter.H_["dPhi_data"].h_.GetXaxis().GetBinCenter(ibin) 
        D_nBKG[ibin-1]         = (plotter.H_["htmp"].h_.IntegralAndError(0, bindPhiLimitBkg,D_tmp) - plotter.H_["htmpHS"].h_.Integral(0, bindPhiLimitBkg) ) * (Signal_dPhiMax - Signal_dPhiMin)/PUControl_dPhiMax
        D_nBKGerr[ibin-1]      = (Signal_dPhiMax - Signal_dPhiMin)/PUControl_dPhiMax * math.sqrt(D_tmp*D_tmp + (HSuncert*plotter.H_["htmpHS"].h_.Integral(0, bindPhiLimitBkg))*(HSuncert*plotter.H_["htmpHS"].h_.Integral(0, bindPhiLimitBkg)))
        print "data", plotter.H_["dPhi_data"].h_.GetXaxis().GetBinCenter(ibin) , 
        print (plotter.H_["htmp"].h_.IntegralAndError(0, bindPhiLimitBkg,D_tmp)  ) * (Signal_dPhiMax - Signal_dPhiMin)/PUControl_dPhiMax,
        print (plotter.H_["htmp"].h_.IntegralAndError(0, bindPhiLimitBkg,D_tmp) - plotter.H_["htmpHS"].h_.Integral(0, bindPhiLimitBkg) ) * (Signal_dPhiMax - Signal_dPhiMin)/PUControl_dPhiMax

    D_tgraph  = R.TGraphErrors(nbins, D_eta, D_nBKG, D_etaerr , D_nBKGerr )

    H=Hist("", "")
    H.sethist(D_tgraph)
    plotter.H_["PUcontaminationData"]=H
    plotter.H_["PUcontaminationData"].isData()
    plotter.H_["PUcontaminationData"].setdrawstyle("P")
    plotter.SetAtlasInfo("Internal")
    plotter.Draw(["PUcontamination", "PUcontaminationTruth","PUcontaminationData"], 
                 VarName, "N_{j}^{PU} in signal region", False, True, [lumiLabel, MClabel, jetlabel, ptlabel,Zsel], ["MC prediction", "MC","Data"])
       
    

    #----------------------------------------------------------------------------------------------
    # Different JVT cuts

    for JVTcut in CUTVal:
        # -------------------- efficiency truth --------------------------------------------------------
        plotter.GetHistoFromTree("EtaHS"   , tfileMC, tree, Variable, "Weight*("+HSDef+" && "+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+")"                      ,VariableNBins, VariableLow, VariableHigh)
        plotter.GetHistoFromTree("EtaHSJVT", tfileMC, tree, Variable, "Weight*("+HSDef+" && "+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+"&&"+CUT+str(JVTcut)+")" ,VariableNBins, VariableLow, VariableHigh)
        plotter.Divide("EtaHSJVT", "EtaHS", "EffTruth", "B")
        
        # Sherpa
        plotter.GetHistoFromTree("EtaHS_2"   , tfileMC_pow, tree, Variable, "Weight*("+HSDef+" && "+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+")"                      ,VariableNBins, VariableLow, VariableHigh)
        plotter.GetHistoFromTree("EtaHSJVT_2", tfileMC_pow, tree, Variable, "Weight*("+HSDef+" && "+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+"&&"+CUT+str(JVTcut)+")" ,VariableNBins, VariableLow, VariableHigh)
        plotter.Divide("EtaHSJVT_2", "EtaHS_2", "EffTruth_2", "B")

        plotter.Draw(["EffTruth", "EffTruth_2"], VarName, "HS efficiency", False, True, [lumiLabel, "Z#rightarrow#mu#mu", jetlabel, ptlabel,Zsel, "JVT > "+str(JVTcut)], ["MC, Sherpa", "MC, PowHeg"])

        # ------------ truth efficiency for systematics --------------------------------
        # Serpa
        plotter.GetHistoFromTree("EtaHSSyst"   ,tfileMC,tree,VariableSyst,"Weight*("+HSDef+"&&"+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+")"                      ,VariableNBinsSyst, VariableLowSyst, VariableHighSyst)
        plotter.GetHistoFromTree("EtaHSJVTSyst",tfileMC,tree,VariableSyst,"Weight*("+HSDef+"&&"+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+"&&"+CUT+str(JVTcut)+")" ,VariableNBinsSyst, VariableLowSyst, VariableHighSyst)
        plotter.Divide("EtaHSJVTSyst", "EtaHSSyst", "EffTruthSyst", "B")
        
        # Powheg
        plotter.GetHistoFromTree("EtaHSSyst2"   ,tfileMC_pow,tree,VariableSyst,"Weight*("+HSDef+"&&"+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+")"                     ,VariableNBinsSyst,VariableLowSyst,VariableHighSyst)
        plotter.GetHistoFromTree("EtaHSJVTSyst2",tfileMC_pow,tree,VariableSyst,"Weight*("+HSDef+"&&"+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+"&&"+CUT+str(JVTcut)+")",VariableNBinsSyst,VariableLowSyst,VariableHighSyst)
        plotter.Divide("EtaHSJVTSyst2", "EtaHSSyst2", "EffTruthSyst2", "B")
        
        plotter.Clone("EffTruthSyst","EffTruthSyst_diff")
        plotter.H_["EffTruthSyst_diff"].h_.Add(plotter.H_["EffTruthSyst2"].h_, -1)
        plotter.Divide("EffTruthSyst_diff","EffTruthSyst", "EffTruthSyst_diff_frac")

        plotter.Draw(["EffTruthSyst_diff"],      VarName, "#Delta #epsilon",            False, True, [lumiLabel, "Z#rightarrow#mu#mu", jetlabel, ptlabel, Zsel,"JVT > "+str(JVTcut)], ["Sherpa - Powheg"])
        plotter.Draw(["EffTruthSyst_diff_frac"], VarName, "#Delta #epsilon / #epsilon", False, True, [lumiLabel, "Z#rightarrow#mu#mu", jetlabel, ptlabel, Zsel,"JVT > "+str(JVTcut)], ["(Sherpa - Powheg) / Sherpa" ])

        plotter.Clone("EffTruthSyst_diff_frac","EffTruthSyst_diff_frac_smooth")
        plotter.H_["EffTruthSyst_diff_frac_smooth"].h_.Smooth(NSmoothSyst)
        plotter.Draw(["EffTruthSyst_diff_frac","EffTruthSyst_diff_frac_smooth"], VarName, "#Delta #epsilon / #epsilon", False, True, [lumiLabel, "Z#rightarrow#mu#mu", jetlabel, ptlabel,Zsel, "JVT > "+str(JVTcut)], ["nominal", "smoothed" ])


        # ---------------------- efficiency reco ---------------------------------------------------------
        plotter.GetHistoFromTree("EtaReco"   , tfileMC, tree, Variable, "Weight*("+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+")"         ,             VariableNBins, VariableLow, VariableHigh)
        plotter.GetHistoFromTree("EtaRecoJVT", tfileMC, tree, Variable, "Weight*("+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+"&&"+CUT+str(JVTcut)+")" ,VariableNBins, VariableLow, VariableHigh)
        plotter.SetMCScale(ScaleFact,    ["EtaReco","EtaRecoJVT"])
        plotter.Divide("EtaRecoJVT", "EtaReco", "EffReco", "B")

        # correct for PU contamination
        plotter.Clone("EtaReco", "EtaRecoCorr")
        for ibin in range(1, plotter.H_["EtaRecoCorr"].h_.GetNbinsX()+1):
            uncorr    = plotter.H_["EtaRecoCorr"].h_.GetBinContent(ibin)
            uncorrErr = plotter.H_["EtaRecoCorr"].h_.GetBinError(ibin)
            bincenter = plotter.H_["EtaRecoCorr"].h_.GetBinCenter(ibin)
            correction= plotter.H_["PUcontamination"].h_.Eval(bincenter)
            corrErr   = plotter.H_["PUcontamination"].h_.GetErrorY(ibin-1)
            x = R.Double(0); y = R.Double(0)
            plotter.H_["PUcontamination"].h_.GetPoint(ibin-1, x, y)
            print "center", bincenter, "x", x, "y", y, "err", corrErr 
            plotter.H_["EtaRecoCorr"].h_.SetBinContent(ibin, uncorr - correction)
            plotter.H_["EtaRecoCorr"].h_.SetBinError  (ibin, math.sqrt(uncorrErr*uncorrErr + corrErr*corrErr))


        plotter.Divide("EtaRecoJVT", "EtaRecoCorr", "EffRecoCorr", "B")
        plotter.Draw(["EffTruth","EffReco","EffRecoCorr"],VarName, "HS efficiency", False, True, [lumiLabel, MClabel, jetlabel, ptlabel,Zsel], ["MC", "MC reco uncorr", "MC reco corr"])


        # ---------------------- efficiency data ---------------------------------------------------------
        plotter.GetHistoFromTree("EtaData"   , tfileData, tree, Variable, "("+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+")"          ,            VariableNBins, VariableLow, VariableHigh)
        plotter.GetHistoFromTree("EtaDataJVT", tfileData, tree, Variable, "("+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+"&&"+CUT+str(JVTcut)+")" ,VariableNBins, VariableLow, VariableHigh)
        plotter.Divide("EtaDataJVT", "EtaData", "EffData", "B")
        plotter.H_["EffData"].isData()
        plotter.Draw(["EffTruth","EffReco","EffRecoCorr", "EffData"], VarName , "HS efficiency", False, True, [lumiLabel, MClabel, jetlabel, ptlabel,Zsel], ["MC truth", "MC reco uncorr", "MC reco corr", "data uncorr"])


        # correct for PU contamination
        plotter.Clone("EtaData", "EtaDataCorr")
        for ibin in range(1, plotter.H_["EtaDataCorr"].h_.GetNbinsX()+1):
            uncorr    = plotter.H_["EtaDataCorr"].h_.GetBinContent(ibin)
            uncorrErr = plotter.H_["EtaDataCorr"].h_.GetBinError(ibin)
            bincenter = plotter.H_["EtaDataCorr"].h_.GetBinCenter(ibin)
            correction= plotter.H_["PUcontaminationData"].h_.Eval(bincenter)
            corrErr   = plotter.H_["PUcontaminationData"].h_.GetErrorY(ibin-1)
            x = R.Double(0); y = R.Double(0)
            plotter.H_["PUcontaminationData"].h_.GetPoint(ibin-1, x, y)
            print "center", bincenter, "x", x, "y", y, "err", corrErr 
            plotter.H_["EtaDataCorr"].h_.SetBinContent(ibin, uncorr - correction)
            if(doDataErr):
                plotter.H_["EtaDataCorr"].h_.SetBinError  (ibin, math.sqrt(uncorrErr*uncorrErr + corrErr*corrErr))
            else:
                plotter.H_["EtaDataCorr"].h_.SetBinError  (ibin, uncorrErr)

        plotter.Divide("EtaDataJVT", "EtaDataCorr", "EffDataCorr", "B")
        plotter.H_["EffDataCorr"].isData()
        plotter.Draw(["EffTruth","EffReco","EffRecoCorr", "EffData","EffDataCorr"], VarName , "HS efficiency", False, True, 
                     [lumiLabel, MClabel, jetlabel, ptlabel,Zsel], ["MC truth", "MC reco uncorr", "MC reco corr", "data uncorr", "data"])

        # final efficiency plot 
        plotter.Draw(["EffTruth", "EffTruth_2", "EffRecoCorr", "EffDataCorr"], VarName , "Efficiency", False, True, 
                     [lumiLabel, MClabel, jetlabel, ptlabel,Zsel, "JVT > "+str(JVTcut)], ["MC", "MC Powheg", "MC pred",  "Data"])
        
        # final efficiency plot 
        plotter.SetAtlasInfo("Simulation Internal")
        plotter.Draw(["EffTruth","EffRecoCorr"],  VarName , "Efficiency", False, True, 
                     [lumiLabel, MClabel, jetlabel, ptlabel,Zsel,"JVT > "+str(JVTcut)], ["MC", "MC pred"  ])
        
        # final efficiency plot 
        plotter.SetAtlasInfo("Preliminary")
        plotter.Draw(["EffRecoCorr", "EffDataCorr"], VarName , "Efficiency", False, True, 
                     [lumiLabel, MClabel, jetlabel, ptlabel,Zsel, "JVT > "+str(JVTcut)], ["MC",  "Data"])
        
        plotter.drawMCsum_=False
        plotter.SetAtlasInfo("Preliminary")
        plotter.DataToMC(["EffRecoCorr"], "EffDataCorr" , VarName, "Efficiency", [lumiLabel, MClabel, jetlabel, ptlabel,Zsel, "JVT > "+str(JVTcut)], ["MC"])

        # final plot with systmatics
        plotter.Clone("EffRecoCorr", "EffRecoCorrSyst")
        plotter.Clone("EffRecoCorr", "EffRecoCorrSystSmooth")
        for ibin in range(1, plotter.H_["EffRecoCorrSyst"].h_.GetNbinsX()+1):
            print "bin ", ibin, "error", plotter.H_["EffRecoCorrSyst"].h_.GetBinError(ibin)
            bincenter   =  plotter.H_["EffRecoCorrSyst"].h_.GetBinCenter(ibin)
            stat        = plotter.H_["EffRecoCorrSyst"].h_.GetBinError(ibin)
            syst        = math.fabs(plotter.H_["EffTruth"].h_.GetBinContent(ibin) - plotter.H_["EffTruth_2"].h_.GetBinContent(ibin))/plotter.H_["EffTruth"].h_.GetBinContent(ibin) * plotter.H_["EffRecoCorrSyst"].h_.GetBinContent(ibin)

            # find bin for truth syst hist
            systbin     = plotter.H_["EffTruthSyst_diff_frac_smooth"].h_.FindFixBin(math.fabs(bincenter))
            print       ibin, "bin for syst " , systbin 
            systsmooth  = math.fabs(plotter.H_["EffTruthSyst_diff_frac_smooth"].h_.GetBinContent(systbin))* plotter.H_["EffRecoCorrSyst"].h_.GetBinContent(ibin)

            total       = math.sqrt(stat*stat + syst*syst)
            totalsmooth = math.sqrt(stat*stat + systsmooth*systsmooth)

            totalFrac   = total / plotter.H_["EffRecoCorrSyst"].h_.GetBinContent(ibin)
            totalFracSm = totalsmooth / plotter.H_["EffRecoCorrSyst"].h_.GetBinContent(ibin)
            print "bin", ibin, "eff", plotter.H_["EffRecoCorrSyst"].h_.GetBinContent(ibin), "stat", stat, "syst", syst, "syst smooth", systsmooth, "tot", total, "total smooth", totalsmooth, "totalFrac", totalFrac, "totalFracSmooth", totalFracSm
            plotter.H_["EffRecoCorrSyst"].h_.SetBinError(ibin, total) 
            plotter.H_["EffRecoCorrSystSmooth"].h_.SetBinError(ibin, totalsmooth) 

        
        plotter.DataToMCSyst("EffRecoCorrSyst","EffRecoCorr","EffDataCorr" ,       VarName, "Efficiency", [lumiLabel, MClabel, jetlabel, ptlabel,Zsel, "JVT > "+str(JVTcut)], ["Total Uncertainty old", "MC", "Data"])
        plotter.DataToMCSyst("EffRecoCorrSystSmooth","EffRecoCorr","EffDataCorr" , VarName, "Efficiency", [lumiLabel, MClabel, jetlabel, ptlabel,Zsel, "JVT > "+str(JVTcut)], ["Total Uncertainty", "MC", "Data"])


        # true HS efficiency
        plotter.GetHistoFromTree("EtaHStrue"   , tfileMC, tree, Variable, "Weight*("+TrueHSDef+" && "+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+")"                      ,VariableNBins, VariableLow, VariableHigh)
        plotter.GetHistoFromTree("EtaHStrueJVT", tfileMC, tree, Variable, "Weight*("+TrueHSDef+" && "+dPhiVar+">"+str(Signal_dPhiMin)+"&&"+dPhiVar+"<"+str(Signal_dPhiMax)+"&&"+JetSelection+"&&"+CUT+str(JVTcut)+")" ,VariableNBins, VariableLow, VariableHigh)
        plotter.Divide("EtaHStrueJVT", "EtaHStrue", "EffVeryTruth", "B")
        
        plotter.DataToMCSyst("EffRecoCorrSyst","EffRecoCorr","EffDataCorr",  VarName, "Efficiency", [lumiLabel, MClabel, jetlabel, ptlabel,Zsel,"JVT > "+str(JVTcut)], ["Total Uncertainty", "MC", "Data", "MC HS"], ["EffVeryTruth"])



