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
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/TopSelection/"
    tfileData = path+"20140215_PileUpStudies_TopSelection_rev351830.MuonsSMWZPeriodABCDEGHIJL.jetmet2012pileupcustom.root"
    tfileMC   = path+"20140215.15.13_PileUpStudies_TopSelection_rev351830.PowhegPythia_ttbar_LeptonFilterSMWZ.jetmet2012pileupcustom.root"

    # labels 
    LumiSqrtS= "#sqrt{s}=8 TeV"
    jetlabel = "Anti-k_{t} LCW+JES R=0.4";
    etalabel = "|#eta| < 2.4";
    gaplabel = " ";
    lumiLabel= "#sqrt{s}=8 TeV, L = 5.8 fb^{-1}"
    normalizationLabel = "MC normalized to data"

    # Histos from Tree
    tree               = "JetTree0"
    
    
    #++++++++++++++++ BASIC SETUP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = Plotter()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Simulation Preliminary")

    
    ## Jet Mult vs NPV from Tree
    plotter.GetHistoFromTree("NPVVsNjets20D",      tfileData,  tree , "Mu", "(Jindex!=-999 && fabs(Jeta)<2.4 && Jpt>20                                         )", 9, 4, 41)
    plotter.GetHistoFromTree("NPVVsNjets20DJVF",   tfileData,  tree , "Mu", "(Jindex!=-999 && fabs(Jeta)<2.4 && Jpt>20 &&kNN100trim_pt20to50_Likelihood>0.9    )", 9, 4, 41)
    plotter.GetHistoFromTree("NPVD",               tfileData,  tree , "Mu", "(Jindex==-999                                                                     )", 9, 4, 41)
    plotter.GetHistoFromTree("NPVVsNjets20M",      tfileMC,    tree , "Mu", "(Jindex!=-999 && fabs(Jeta)<2.4 && Jpt>20                                         )", 9, 4, 41)
    plotter.GetHistoFromTree("NPVVsNjets20MJVF",   tfileMC,    tree , "Mu", "(Jindex!=-999 && fabs(Jeta)<2.4 && Jpt>20 &&kNN100trim_pt20to50_Likelihood>0.9    )", 9, 4, 41)
    plotter.GetHistoFromTree("NPVM",               tfileMC,    tree , "Mu", "(Jindex==-999                                                                     )", 9, 4, 41)
    plotter.Divide("NPVVsNjets20D",    "NPVD", "NJetsVsNPVD"   ,"B")
    plotter.Divide("NPVVsNjets20M",    "NPVM", "NJetsVsNPVM"   ,"B")
    plotter.Divide("NPVVsNjets20DJVF", "NPVD", "NJetsVsNPVDJVF","B")
    plotter.Divide("NPVVsNjets20MJVF", "NPVM", "NJetsVsNPVMJVF","B")
    plotter.H_["NJetsVsNPVD"].isData()
    plotter.H_["NJetsVsNPVDJVF"].isData()
    plotter.Draw(["NJetsVsNPVM","NJetsVsNPVD","NJetsVsNPVMJVF","NJetsVsNPVDJVF"], "#mu","#LT Jet multiplicity#GT", False, True, 
                [lumiLabel, jetlabel,"p_{T} > 20 GeV, |#eta|<2.4","Z boson p_{T} > 20 GeV"], ["MC", "Data", "MC, JVT>0.6", "Data, JVT>0.6"])


    # dNj / dNPV
    plotter.GetSlope("NJetsVsNPVM",    "NJetsVsNPVMSlope")
    plotter.GetSlope("NJetsVsNPVMJVF", "NJetsVsNPVMJVFSlope")
    plotter.GetSlope("NJetsVsNPVD",    "NJetsVsNPVDSlope")
    plotter.GetSlope("NJetsVsNPVDJVF", "NJetsVsNPVDJVFSlope")
    plotter.H_["NJetsVsNPVDSlope"].isData()
    plotter.H_["NJetsVsNPVDJVFSlope"].isData()
    plotter.H_["NJetsVsNPVMSlope"].setdrawstyle("ALP")
    plotter.H_["NJetsVsNPVMJVFSlope"].setdrawstyle("LP")
    plotter.H_["NJetsVsNPVDSlope"].setdrawstyle("LP")
    plotter.H_["NJetsVsNPVDJVFSlope"].setdrawstyle("LP")
    plotter.Draw(["NJetsVsNPVMSlope","NJetsVsNPVMJVFSlope","NJetsVsNPVDSlope","NJetsVsNPVDJVFSlope"], "#mu", "#partial #LT N_{j} #GT / #partial #mu", False, True,  
                 [lumiLabel, jetlabel, "p_{T} > 20 GeV, |#eta|<2.1","Z boson p_{T} > 20 GeV"], ["MC", "MC, JVT>0.6","Data", "Data, JVT>0.6"]) 
        
    print plotter.H_.keys()
    

