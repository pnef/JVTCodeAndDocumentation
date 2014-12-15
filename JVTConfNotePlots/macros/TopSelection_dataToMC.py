#!/usr/bin/python
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
sys.path.append(os.path.abspath("/Users/pascal/bin/python/plotting"))
from DataToMCPlottingBase import *
from multiprocessing import Process, Lock, Queue

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
    parser.add_option("--validation"         ,dest="validation"        ,    help="set validation"       , action="store_true")
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/TopSelection/"
    tfile           = path+"../skimmed/skimmedPt20to50Eta2p4.20140213.00.12_PileUpStudies_rev351681.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN20140213.root"
    tfileData       = path+"20140215_PileUpStudies_TopSelection_rev351830.MuonsSMWZPeriodABCDEGHIJL.jetmet2012pileupcustom.root"
    tfileTTbar      = path+"20140215.15.13_PileUpStudies_TopSelection_rev351830.PowhegPythia_ttbar_LeptonFilterSMWZ.jetmet2012pileupcustom.root"
    tfileSingleTop  = path+"20140215.17.25_PileUpStudies_TopSelection_rev351830.SingleTopFullSMWZ.jetmet2012pileupcustom.root"
    tfileWJets      = path+"20140215.15.57_PileUpStudies_TopSelection_rev351830.WJetsFullSMWZ.jetmet2012pileupcustom.root"

    # labels 
    jetlabel = "Anti-k_{t} LCW+JES R=0.4";
    lumiLabel= "#sqrt{s}=8 TeV, L = 20.3 fb^{-1}"
    normalizationLabel = "MC normalized to data"
    MClabel = "Top Selection"

    # Histos from Tree
    tree               = "JetTree0"

    # Signal norm
    SigScaleFact = 20.3*1.*1000
    
    #++++++++++++++++ BASIC SETUP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = DataToMCPlottingBase()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Preliminary")
    plotter.doNormMCtoData_ = True
    plotter.doRatio_ = True
    plotter.exclude_ = [("other","MC_other")]


    
    ##########################################################################################################################
    if options.validation:
        plotter.drawMCcomponents_ = True
        plotter.doRatio_          = True

        signalSelection = " JisBtagged && Jpt>20 && Jpt<50 && fabs(Jeta)<2.4 "
        plotter.GetHistoFromTree("Data",       tfileData,                tree ,"Jpt", "Weight*(1   &&"+signalSelection+")", 15, 0, 60)
        plotter.GetHistoFromTree("TTbar",      tfileTTbar,               tree ,"Jpt", "Weight*(1   &&"+signalSelection+")", 15, 0, 60)
        plotter.GetHistoFromTree("SingleTop",  tfileSingleTop,           tree ,"Jpt", "Weight*(1   &&"+signalSelection+")", 15, 0, 60)
        plotter.GetHistoFromTree("WJets",      tfileWJets,               tree ,"Jpt", "Weight*(1   &&"+signalSelection+")", 15, 0, 60)
        plotter.SetMCScale(SigScaleFact, ["TTbar","SingleTop","WJets"])
        plotter.DataToMC(["TTbar","SingleTop","WJets"],"Data",  "b-j p_{T} [GeV]", "Normalized Events", ["Top Selection"], ["TTbar","SingleTop","WJets"])

        plotter.GetHistoFromTree("Data",       tfileData,                tree ,"NPV", "Weight*(Jindex==-999)", 15, 0, 30)
        plotter.GetHistoFromTree("TTbar",      tfileTTbar,               tree ,"NPV", "Weight*(Jindex==-999)", 15, 0, 30)
        plotter.GetHistoFromTree("SingleTop",  tfileSingleTop,           tree ,"NPV", "Weight*(Jindex==-999)", 15, 0, 30)
        plotter.GetHistoFromTree("WJets",      tfileWJets,               tree ,"NPV", "Weight*(Jindex==-999)", 15, 0, 30)
        plotter.SetMCScale(SigScaleFact, ["TTbar","SingleTop","WJets"])
        plotter.DataToMC(["TTbar","SingleTop","WJets"],"Data", "NPV", "Normalized Events", ["Top Selection"], ["TTbar","SingleTop","WJets"])






    # ----------------- SIGNAL ----------------------------------------------------
    signalSelection = " JisBtagged && Jpt>20 && Jpt<50 && fabs(Jeta)<2.4 "
    
    # JVT leading Jet -----------------------------
    low = -0.2; high = 1.1
    plotter.GetHistoFromTree("Data",             tfileData,   tree , "kNN100trim_pt20to50_Likelihood", "Weight*(1              &&"+signalSelection+")",               20, low, high)

    plotter.GetHistoFromTree("ttMC_HS",          tfileTTbar,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(JisHS&&Jtruthpt>10           &&"+signalSelection+")", 20, low, high)
    plotter.GetHistoFromTree("ttMC_PU",          tfileTTbar,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(JisPU                        &&"+signalSelection+")", 20, low, high)
    plotter.GetHistoFromTree("ttMC_o" ,          tfileTTbar,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(!(JisHS&&Jtruthpt>10)&&!JisPU&&"+signalSelection+")" ,20, low, high)

    plotter.GetHistoFromTree("WMC_HS",           tfileWJets,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(JisHS&&Jtruthpt>10           &&"+signalSelection+")", 20, low, high)
    plotter.GetHistoFromTree("WMC_PU",           tfileWJets,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(JisPU                        &&"+signalSelection+")", 20, low, high)
    plotter.GetHistoFromTree("WMC_o" ,           tfileWJets,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(!(JisHS&&Jtruthpt>10)&&!JisPU&&"+signalSelection+")" ,20, low, high)
    
    plotter.GetHistoFromTree("tMC_HS",       tfileSingleTop,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(JisHS&&Jtruthpt>10           &&"+signalSelection+")", 20, low, high)
    plotter.GetHistoFromTree("tMC_PU",       tfileSingleTop,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(JisPU                        &&"+signalSelection+")", 20, low, high)
    plotter.GetHistoFromTree("tMC_o" ,       tfileSingleTop,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(!(JisHS&&Jtruthpt>10)&&!JisPU&&"+signalSelection+")", 20, low, high)

    plotter.AddHistos(["ttMC_HS" ,"WMC_HS" ,"tMC_HS"],  "MC_HS")
    plotter.AddHistos(["ttMC_PU" ,"WMC_PU" ,"tMC_PU"],  "MC_PU")
    plotter.AddHistos(["ttMC_o"  ,"WMC_o"  ,"tMC_o"],   "MC_other")

    plotter.SetMCScale(SigScaleFact,    ["MC_HS","MC_PU","MC_other"])

    plotter.DataToMC(["MC_HS", "MC_PU","MC_other"], "Data" , "JVT", "Events", [lumiLabel,MClabel,jetlabel,"20<p_{T}<50 GeV, |#eta|<2.4","b-tagged"], [ "HS", "PU", "other"])


    # corrJVF leading Jet -----------------------------
    low = -1; high = 1.1
    plotter.GetHistoFromTree("Data",             tfileData,   tree , "JnPUTrkCorrJVF", "Weight*(1                            &&"+signalSelection+")", 15, low, high)

    plotter.GetHistoFromTree("ttMC_HS",          tfileTTbar,  tree , "JnPUTrkCorrJVF", "Weight*(JisHS&&Jtruthpt>10           &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("ttMC_PU",          tfileTTbar,  tree , "JnPUTrkCorrJVF", "Weight*(JisPU                        &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("ttMC_o" ,          tfileTTbar,  tree , "JnPUTrkCorrJVF", "Weight*(!(JisHS&&Jtruthpt>10)&&!JisPU&&"+signalSelection+")" ,15, low, high)

    plotter.GetHistoFromTree("WMC_HS",           tfileWJets,  tree , "JnPUTrkCorrJVF", "Weight*(JisHS&&Jtruthpt>10           &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("WMC_PU",           tfileWJets,  tree , "JnPUTrkCorrJVF", "Weight*(JisPU                        &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("WMC_o" ,           tfileWJets,  tree , "JnPUTrkCorrJVF", "Weight*(!(JisHS&&Jtruthpt>10)&&!JisPU&&"+signalSelection+")" , 15, low, high)
    
    plotter.GetHistoFromTree("tMC_HS",       tfileSingleTop,  tree , "JnPUTrkCorrJVF", "Weight*(JisHS&&Jtruthpt>10           &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("tMC_PU",       tfileSingleTop,  tree , "JnPUTrkCorrJVF", "Weight*(JisPU                        &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("tMC_o" ,       tfileSingleTop,  tree , "JnPUTrkCorrJVF", "Weight*(!(JisHS&&Jtruthpt>10)&&!JisPU&&"+signalSelection+")",  15, low, high)

    plotter.AddHistos(["ttMC_HS" ,"WMC_HS" ,"tMC_HS"],  "MC_HS")
    plotter.AddHistos(["ttMC_PU" ,"WMC_PU" ,"tMC_PU"],  "MC_PU")
    plotter.AddHistos(["ttMC_o"  ,"WMC_o"  ,"tMC_o"],   "MC_other")

    plotter.SetMCScale(SigScaleFact,    ["MC_HS","MC_PU","MC_other"])
    plotter.DataToMC(["MC_HS", "MC_PU","MC_other"], "Data" , "corrJVF", "Events", [lumiLabel,MClabel,jetlabel,"20<p_{T}<50 GeV, |#eta|<2.4","b-tagged"], [ "HS", "PU", "other"])

    # RpT leading Jet -----------------------------
    low = 0; high = 1.5
    plotter.GetHistoFromTree("Data",             tfileData,   tree , "JHSPVtrkSumOverPt", "Weight*(1              &&"+signalSelection+")", 15, low, high)

    plotter.GetHistoFromTree("ttMC_HS",          tfileTTbar,  tree , "JHSPVtrkSumOverPt", "Weight*(JisHS&&Jtruthpt>10           &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("ttMC_PU",          tfileTTbar,  tree , "JHSPVtrkSumOverPt", "Weight*(JisPU                        &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("ttMC_o" ,          tfileTTbar,  tree , "JHSPVtrkSumOverPt", "Weight*(!(JisHS&&Jtruthpt>10)&&!JisPU&&"+signalSelection+")" , 15, low, high)

    plotter.GetHistoFromTree("WMC_HS",           tfileWJets,  tree , "JHSPVtrkSumOverPt", "Weight*(JisHS&&Jtruthpt>10           &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("WMC_PU",           tfileWJets,  tree , "JHSPVtrkSumOverPt", "Weight*(JisPU                        &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("WMC_o" ,           tfileWJets,  tree , "JHSPVtrkSumOverPt", "Weight*(!(JisHS&&Jtruthpt>10)&&!JisPU&&"+signalSelection+")" , 15, low, high)
    
    plotter.GetHistoFromTree("tMC_HS",       tfileSingleTop,  tree , "JHSPVtrkSumOverPt", "Weight*(JisHS&&Jtruthpt>10           &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("tMC_PU",       tfileSingleTop,  tree , "JHSPVtrkSumOverPt", "Weight*(JisPU                        &&"+signalSelection+")", 15, low, high)
    plotter.GetHistoFromTree("tMC_o" ,       tfileSingleTop,  tree , "JHSPVtrkSumOverPt", "Weight*(!(JisHS&&Jtruthpt>10)&&!JisPU&&"+signalSelection+")",  15, low, high)

    plotter.AddHistos(["ttMC_HS" ,"WMC_HS" ,"tMC_HS"],  "MC_HS")
    plotter.AddHistos(["ttMC_PU" ,"WMC_PU" ,"tMC_PU"],  "MC_PU")
    plotter.AddHistos(["ttMC_o"  ,"WMC_o"  ,"tMC_o"],   "MC_other")

    plotter.SetMCScale(SigScaleFact,    ["MC_HS","MC_PU","MC_other"])

    plotter.DataToMC(["MC_HS", "MC_PU","MC_other"], "Data" , "R_{pT}", "Events", [lumiLabel,MClabel,jetlabel,"20<p_{T}<50 GeV, |#eta|<2.4","b-tagged"], [ "HS", "PU", "other"])




















    print plotter.H_.keys()
    

