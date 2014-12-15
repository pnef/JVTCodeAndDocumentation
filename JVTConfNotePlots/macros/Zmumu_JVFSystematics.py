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
sys.path.append(os.path.abspath("/Users/pascal/Work/CERN/Projects/ATLASAnalysis/PileUpStudies/Code/Scripts/Plotting/"))
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
    parser.add_option("--JVF",                dest="JVF" ,                  help="", action="store_true")
    parser.add_option("--nPUtrkCorrJVF",      dest="nPUtrkCorrJVF" ,                 help="", action="store_true")
    parser.add_option("--NPVcorrJVF",         dest="NPVcorrJVF" ,           help="", action="store_true")
    parser.add_option("--PtRatio",            dest="PtRatio" ,              help="", action="store_true")
    parser.add_option("--PUPtRatio",          dest="PUPtRatio" ,            help="", action="store_true")
    parser.add_option("--ANN",                dest="ANN" ,                  help="", action="store_true")
    parser.add_option("--MC"                 ,dest="ZmumuMC"           ,    help="Zmumu MC version")
    parser.add_option("--nPUTracks"          ,dest="nPUTracks"           ,  help="Zmumu MC version", action="store_true")
    parser.add_option("--BDT"                ,dest="BDT"                 ,  help="Zmumu MC version", action="store_true")
    parser.add_option("--doSyst"             ,dest="doSyst"                 ,  help="Zmumu MC version", action="store_true")

    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/"
    tfileData = path+"2014021x.PileUpStudies_DiLeptonSelection_rev351830.MuonsSMWZPeriodABCDEGHIJL.jetmet2012pileupcustom.root"
    tfileMC= path+"20140214.22.39_PileUpStudies_DiLeptonSelection_rev351830.Zmumu_PowhegPythia8_MC12_COMMON.jetmet2012pileupcustom.root"

    # labels 
    LumiSqrtS= "#sqrt{s}=8 TeV"
    jetlabel = "Anti-k_{t} LCW+JES R=0.4";
    etalabel = "|#eta| < 2.4";
    ptlabel  = "20<p_{T}<50 GeV, |#eta| < 2.4";
    gaplabel = " ";
    lumiLabel= "#sqrt{s}=8 TeV, L = 20.3 fb^{-1}"
    normalizationLabel = "MC normalized to data"
    MClabel = "Powheg+Pythia8 Z#rightarrow #mu#mu"

    # Histos from Tree
    tree               = "JetTree0"
    
    
    #++++++++++++++++ BASIC SETUP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = Plotter()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Preliminary")

    
        
    ##########################################################################################################################
    #### PU Validation Plots ##################################################################################################

    if(options.doSyst):
        # ----------------- SIGNAL ----------------------------------------------------
        signalSelection = "fabs(J0ZDPhi) >2.6 && ZPt>30 && Jpt>20 && Jpt<50 && fabs(Jeta)<2.4 && Jindex==0 && kNN100trim_pt20to50_Likelihood>=0"
        

        plotter.GetHistoFromTree("JVFDataSyst",  tfileData,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(1    &&"+signalSelection+" )", 100, 0, 1.1)
        plotter.GetHistoFromTree("JVFMCSyst"  ,  tfileMC  ,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(1    &&"+signalSelection+" )", 100, 0, 1.1)
        plotter.GetJVFSystHist2("JVFDataSyst", "JVFDataSyst_eff")
        plotter.GetJVFSystHist2("JVFMCSyst",   "JVFMCSyst_eff"  )
        plotter.SetRange("JVFDataSyst_eff", "x", 0, 1., "JVFDataSyst_eff")

        plotter.DrawRatio("JVFDataSyst_eff", "JVFMCSyst_eff", "EffRatio")
        plotter.H_["JVFDataSyst_eff"].isData()
        plotter.H_["JVFDataSyst_eff"].setdrawstyle("P")
        plotter.H_["JVFMCSyst_eff"]  .setopen()
        plotter.H_["JVFMCSyst_eff"]  .setdrawstyle("P")        
        plotter.H_["EffRatio"]       .setStyle(20,1,"P","", "ratio")
        legends = ["data", "MC", "data / MC"]

        plotter.Draw(["JVFDataSyst_eff","JVFMCSyst_eff"], "JVT cut value", "efficiency", 
                     False, True, [jetlabel,MClabel, ptlabel, lumiLabel, "Z p_{T} > 30 GeV"], legends)


        # BDT leading Jet
        plotter.GetHistoFromTree("JVFData",  tfileData,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(1    &&"+signalSelection+" )", 30, -0.11, 1.2)
        plotter.GetHistoFromTree("JVFMC_PU", tfileMC  ,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(JisPU&&"+signalSelection+" )", 30, -0.11, 1.2)
        plotter.GetHistoFromTree("JVFMC_HS", tfileMC  ,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*(JisHS&&"+signalSelection+" )", 30, -0.11, 1.2)
        plotter.GetHistoFromTree("JVFMC_Inc",tfileMC  ,  tree , "kNN100trim_pt20to50_Likelihood", "Weight*( 1   &&"+signalSelection+" )", 30, -0.11, 1.2)
        scale = plotter.H_["JVFData"].h_.Integral() / plotter.H_["JVFMC_Inc"].h_.Integral()
        print scale
        plotter.SetMCScale(scale,["JVFMC_PU","JVFMC_HS","JVFMC_Inc"])
        plotter.DrawRatioMCUncert("JVFMC_Inc", "Data / MC")
        plotter.DrawRatio("JVFData", "JVFMC_Inc", "JVFRatio")
        plotter.H_["JVFRatio"].setStyle(20,1,"EX","", "Data / MC")
        plotter.H_["JVFData"].isData()
        plotter.H_["JVFMC_Inc"].setdrawstyle("hist")
        plotter.H_["JVFMC_PU"].setdrawstyle("hist");  plotter.H_["JVFMC_PU"].setline()
        plotter.H_["JVFMC_HS"].setdrawstyle("hist");  plotter.H_["JVFMC_HS"].setline()
        legends = ["MC inc","MC, PU", "MC, HS","Data"]
        plotter.Draw(["JVFMC_Inc","JVFMC_PU","JVFMC_HS","JVFData"], "leading jet BDT", "events", 
                     False, True, [jetlabel,MClabel, ptlabel, lumiLabel, "Z p_{T} > 40 GeV"], legends)

        

    print plotter.H_.keys()
    

