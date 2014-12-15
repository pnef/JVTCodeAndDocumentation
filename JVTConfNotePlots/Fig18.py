#!/usr/bin/env python
##################################################################
# pascal nef                             April 23rd, 2013        #
##################################################################
import ROOT
import sys
import os
import re
import math
import time
from optparse import OptionParser
from sys import argv,exit
from plottingBase import *

##################################################################


# Main ------------------------------------------------------------------------
if __name__ == "__main__":

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # parse arguments --------------------------------------------
    usage = """
    usage: %prog [options] InputFile.root
    	
    """
    parser = OptionParser(usage)
    
    parser.add_option("--verbose"            ,dest="verbose"           ,    help="set verbose"          , action="store_true")
    parser.add_option("--doEffFakeVsNPV"     ,dest="doEffFakeVsNPV"    ,    help="do Eff vs Fake vs NPV", action="store_true")
    parser.add_option("--doJetMult"          ,dest="doJetMult"         ,    help="do jet multiplicity", action="store_true")
    parser.add_option("--doROC"              ,dest="doROC"             ,    help="do ROC curve", action="store_true")
    parser.add_option("--drawFromTree"       ,dest="drawFromTree"      ,    help="draw from Tree", action="store_true")
    parser.add_option("--doValidationPlots"  ,dest="doValidationPlots" ,    help="validation plots", action="store_true")
    parser.add_option("--doPUValidationPlots",dest="doPUValidationPlots" ,  help="PU validation plots", action="store_true")
    parser.add_option("--dataToMC"           ,dest="dataToMC"          ,    help="", action="store_true")
    parser.add_option("--PUselectionPlots"   ,dest="PUselectionPlots"  ,    help="", action="store_true")
    parser.add_option("--doRhoPU"            ,dest="doRhoPU"           ,    help="calculate PU rho from tracks", action="store_true")
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/"
    tfileData = path+"2014021x.PileUpStudies_DiLeptonSelection_rev351830.MuonsSMWZPeriodABCDEGHIJL.jetmet2012pileupcustom.root"
    tfileMC   = path+"20140305.08.12_PileUpStudies_DiLeptonSelection_351830.Zmumu_MC12_Sherpa_COMMON.jetmet2012pileupcustom.root"

    # labels 
    LumiSqrtS= "#sqrt{s}=8 TeV"
    jetlabel = "Anti-k_{t} LCW+JES R=0.4";
    etalabel = "|#eta| < 2.4";
    gaplabel = " ";
    lumiLabel= "#sqrt{s}=8 TeV, L = 20.3 fb^{-1}"
    normalizationLabel = "MC normalized to data"
    MClabel = "Sherpa Z#rightarrow#mu#mu"

    # Histos from Tree
    tree               = "JetTree0"
    
    
    #++++++++++++++++ BASIC SETUP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = Plotter()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Preliminary")

    plotter.COLORSSolid  = [ROOT.kViolet+9, ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.COLORSOpen   = [ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.MARKERSSolid = [21, 22, 23, 33, 34, 28, 29] 
    plotter.MARKERSOpen  = [25, 26, 32, 27, 28, 30]
    
    ## Jet Mult vs NPV from Tree
    plotter.GetHistoFromTree("NPVVsNjets20D",      tfileData,  tree , "Mu", "Weight*(Jindex!=-999 && fabs(Jeta)<2.4 && Jpt>20                                         && ZPt>20)", 8, 4, 41)
    plotter.GetHistoFromTree("NPVVsNjets20DJVF",   tfileData,  tree , "Mu", "Weight*(Jindex!=-999 && fabs(Jeta)<2.4 && Jpt>20 &&kNN100trim_pt20to50_Likelihood>0.7    && ZPt>20)", 8, 4, 41)
    plotter.GetHistoFromTree("NPVD",               tfileData,  tree , "Mu", "Weight*(Jindex==-999                                                                     && ZPt>20)", 8, 4, 41)
    plotter.GetHistoFromTree("NPVVsNjets20M",      tfileMC,    tree , "Mu", "Weight*(Jindex!=-999 && fabs(Jeta)<2.4 && Jpt>20                                         && ZPt>20)", 8, 4, 41)
    plotter.GetHistoFromTree("NPVVsNjets20MJVF",   tfileMC,    tree , "Mu", "Weight*(Jindex!=-999 && fabs(Jeta)<2.4 && Jpt>20 &&kNN100trim_pt20to50_Likelihood>0.7    && ZPt>20)", 8, 4, 41)
    plotter.GetHistoFromTree("NPVM",               tfileMC,    tree , "Mu", "Weight*(Jindex==-999                                                                     && ZPt>20)", 8, 4, 41)
    plotter.Divide("NPVVsNjets20D",    "NPVD", "NJetsVsNPVD"   ,"B")
    plotter.Divide("NPVVsNjets20M",    "NPVM", "NJetsVsNPVM"   ,"B")
    plotter.Divide("NPVVsNjets20DJVF", "NPVD", "NJetsVsNPVDJVF","B")
    plotter.Divide("NPVVsNjets20MJVF", "NPVM", "NJetsVsNPVMJVF","B")
    plotter.H_["NJetsVsNPVD"].isData()
    plotter.H_["NJetsVsNPVDJVF"].isData()
    plotter.H_["NJetsVsNPVM"].setdrawstyle("EX0")
    plotter.H_["NJetsVsNPVD"].setdrawstyle("EX0")
    plotter.H_["NJetsVsNPVMJVF"].setdrawstyle("EX0")
    plotter.H_["NJetsVsNPVDJVF"].setdrawstyle("EX0")

    plotter.Draw(["NJetsVsNPVM","NJetsVsNPVD","NJetsVsNPVMJVF","NJetsVsNPVDJVF"], "#mu","#LT N_{j} #GT", False, True, 
                [lumiLabel, MClabel, jetlabel,"p_{T} > 20 GeV, |#eta|<2.4","Z boson p_{T} > 20 GeV"], ["MC", "Data", "MC, JVT>0.7", "Data, JVT>0.7"])


    # dNj / dNPV
    plotter.GetSlope("NJetsVsNPVM",    "NJetsVsNPVMSlope")
    plotter.GetSlope("NJetsVsNPVMJVF", "NJetsVsNPVMJVFSlope")
    plotter.GetSlope("NJetsVsNPVD",    "NJetsVsNPVDSlope")
    plotter.GetSlope("NJetsVsNPVDJVF", "NJetsVsNPVDJVFSlope")
    plotter.H_["NJetsVsNPVDSlope"].isData()
    plotter.H_["NJetsVsNPVDJVFSlope"].isData()
    plotter.H_["NJetsVsNPVMSlope"].setdrawstyle("AP")
    plotter.H_["NJetsVsNPVMJVFSlope"].setdrawstyle("P")
    plotter.H_["NJetsVsNPVDSlope"].setdrawstyle("P")
    plotter.H_["NJetsVsNPVDJVFSlope"].setdrawstyle("P")
    canv, leg = plotter.Draw(["NJetsVsNPVMSlope","NJetsVsNPVDSlope","NJetsVsNPVMJVFSlope","NJetsVsNPVDJVFSlope"], "#mu", "#partial #LT N_{j} #GT / #partial #mu", False, True,  
                 [lumiLabel, MClabel, jetlabel, "p_{T} > 20 GeV, |#eta|<2.4","Z boson p_{T} > 20 GeV"], ["MC", "Data","MC, JVT>0.7", "Data, JVT>0.7"]) 


    class Linear:
        def __call__( self, x, par ):
            return par[0] + x[0]*par[1]

    fdataNoCut = ROOT.TF1('myf_dataNoCut',Linear(),7, 34.,2)
    fMCNoCut   = ROOT.TF1('myf_MCNoCut',Linear(),7, 34.,2)
    fdataCut   = ROOT.TF1('myf_dataCut',Linear(),7, 34.,2)
    fMCCut     = ROOT.TF1('myf_MCCut',Linear(),7, 34.,2)

    fdataNoCut.SetLineColor(1)
    fMCNoCut  .SetLineColor(ROOT.kMagenta+3)
    fdataCut  .SetLineColor(ROOT.kViolet+9)
    fMCCut    .SetLineColor(1)

    
    plotter.H_["NJetsVsNPVMSlope"]   .h_.Fit("myf_MCNoCut","","",7,34)
    plotter.H_["NJetsVsNPVMJVFSlope"].h_.Fit("myf_MCCut"  ,"","",7,34)
    plotter.H_["NJetsVsNPVDSlope"]   .h_.Fit("myf_dataNoCut","","",7,34)
    plotter.H_["NJetsVsNPVDJVFSlope"].h_.Fit("myf_dataCut","","",7,34)
    canv.Draw()

    save(canv)
        
    print plotter.H_.keys()
    

