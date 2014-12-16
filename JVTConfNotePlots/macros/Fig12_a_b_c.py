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


# Main ------------------------------------------------------------------------
if __name__ == "__main__":

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # parse arguments --------------------------------------------
    usage = """
    usage: %prog [options] InputFile.root
    	
    """
    parser = OptionParser(usage)
    
    parser.add_option("--verbose"            ,dest="verbose"           ,    help="set verbose"          , action="store_true")
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/"
    path = "/atlas/output/pnef/"
    tfileData = path+"2014021x.PileUpStudies_DiLeptonSelection_rev351830.MuonsSMWZPeriodABCDEGHIJL.jetmet2012pileupcustom.root"
    tfileMC   = path+"20140305.08.12_PileUpStudies_DiLeptonSelection_351830.Zmumu_MC12_Sherpa_COMMON.jetmet2012pileupcustom.root"

    # labels 
    LumiSqrtS= "#sqrt{s}=8 TeV"
    jetlabel = "Anti-k_{t} LCW+JES R=0.4";
    etalabel = "|#eta| < 2.4";
    ptlabel  = "p_{T} > 20 GeV, |#eta| < 2.4";
    gaplabel = " ";
    lumiLabel= "#sqrt{s} = 8 TeV, L = 20.3 fb^{-1}"
    normalizationLabel = "MC normalized to data"
    MClabel  = ""
    MClabel = "Sherpa Z#rightarrow#mu#mu"

    # Histos from Tree
    tree               = "JetTree0"
    
    # Signal norm
    SigScaleFact = 20.3*1.2*1000*1.15
    
    #++++++++++++++++ BASIC SETUP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = DataToMCPlottingBase()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Preliminary")
    plotter.doRatio_ = True
    plotter.printInfo_ = True
    plotter.doNormMCtoData_ = False
    plotter.drawMCcomponents_ = True
    plotter.exclude_ = [("other","JVFMC_other")]

    
    



    # ----------------- SIGNAL ----------------------------------------------------
    signalSelection = "ZPt>30 && Jpt>20 && Jpt<50 && fabs(Jeta)<2.4 && Jindex==0"

    # JVT leading Jet -----------------------------
    low = -0.2; high = 1.1
    plotter.GetHistoFromTree("JVFData",      tfileData,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*( 1                            &&"+signalSelection+")", 35, low, high)
    plotter.GetHistoFromTree("JVFMC_HS",     tfileMC  ,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(   JisHS&&Jtruthpt>10         &&"+signalSelection+")", 35, low, high)
    plotter.GetHistoFromTree("JVFMC_nHS",    tfileMC  ,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*( !(JisHS&&Jtruthpt>10)        &&"+signalSelection+")", 35, low, high)
    plotter.GetHistoFromTree("JVFMC_PU",     tfileMC  ,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*( JisPU                        &&"+signalSelection+")", 35, low, high)
    plotter.GetHistoFromTree("JVFMC_other",  tfileMC  ,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*( !JisPU&&!(JisHS&&Jtruthpt>10)&&"+signalSelection+")", 35, low, high)

    plotter.SetMCScale(SigScaleFact, ["JVFMC_HS"])
    plotter.SetMCScale(SigScaleFact,    ["JVFMC_PU","JVFMC_other"])

    plotter.DataToMC(["JVFMC_HS", "JVFMC_PU","JVFMC_other"], "JVFData" , "JVT", "Events", [lumiLabel,MClabel,jetlabel,"20<p_{T}<50 GeV, |#eta|<2.4", "Z p_{T}>30 GeV", "inclusive"], [  "HS", "PU", "other"])

    # corrJVF leading Jet -----------------------------
    low = -1; high = 1.1
    plotter.GetHistoFromTree("JVFData",      tfileData,  tree , "JnPUTrkCorrJVF", "Weight*( 1                            &&"+signalSelection+")", 30, low, high)
    plotter.GetHistoFromTree("JVFMC_HS",     tfileMC  ,  tree , "JnPUTrkCorrJVF", "Weight*(   JisHS&&Jtruthpt>10         &&"+signalSelection+")", 30, low, high)
    plotter.GetHistoFromTree("JVFMC_nHS",    tfileMC  ,  tree , "JnPUTrkCorrJVF", "Weight*( !(JisHS&&Jtruthpt>10)        &&"+signalSelection+")", 30, low, high)
    plotter.GetHistoFromTree("JVFMC_PU",     tfileMC  ,  tree , "JnPUTrkCorrJVF", "Weight*( JisPU                        &&"+signalSelection+")", 30, low, high)
    plotter.GetHistoFromTree("JVFMC_other",  tfileMC  ,  tree , "JnPUTrkCorrJVF", "Weight*( !JisPU&&!(JisHS&&Jtruthpt>10)&&"+signalSelection+")", 30, low, high)

    plotter.SetMCScale(SigScaleFact, ["JVFMC_HS"])
    plotter.SetMCScale(SigScaleFact,    ["JVFMC_PU","JVFMC_other"])

    plotter.DataToMC(["JVFMC_HS", "JVFMC_PU","JVFMC_other"], "JVFData" , "corrJVF", "Events", [lumiLabel,MClabel,jetlabel,"20<p_{T}<50 GeV, |#eta|<2.4", "Z p_{T}>30 GeV", "inclusive"], [ "HS", "PU", "other"])



    # R_{pT} leading Jet -----------------------------
    low = 0; high = 1.5
    plotter.GetHistoFromTree("JVFData",      tfileData,  tree , "JHSPVtrkSumOverPt", "Weight*( 1                            &&"+signalSelection+")", 30, low, high)
    plotter.GetHistoFromTree("JVFMC_HS",     tfileMC  ,  tree , "JHSPVtrkSumOverPt", "Weight*(   JisHS&&Jtruthpt>10         &&"+signalSelection+")", 30, low, high)
    plotter.GetHistoFromTree("JVFMC_nHS",    tfileMC  ,  tree , "JHSPVtrkSumOverPt", "Weight*( !(JisHS&&Jtruthpt>10)        &&"+signalSelection+")", 30, low, high)
    plotter.GetHistoFromTree("JVFMC_PU",     tfileMC  ,  tree , "JHSPVtrkSumOverPt", "Weight*( JisPU                        &&"+signalSelection+")", 30, low, high)
    plotter.GetHistoFromTree("JVFMC_other",  tfileMC  ,  tree , "JHSPVtrkSumOverPt", "Weight*( !JisPU&&!(JisHS&&Jtruthpt>10)&&"+signalSelection+")", 30, low, high)

    plotter.SetMCScale(SigScaleFact, ["JVFMC_HS"])
    plotter.SetMCScale(SigScaleFact,    ["JVFMC_PU","JVFMC_other"])

    plotter.DataToMC(["JVFMC_HS", "JVFMC_PU","JVFMC_other"], "JVFData" , "R_{pT}", "Events", [lumiLabel,MClabel,jetlabel,"20<p_{T}<50 GeV, |#eta|<2.4", "Z p_{T}>30 GeV", "inclusive"], [ "HS", "PU", "other"])





    print plotter.H_.keys()
    

