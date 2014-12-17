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
from PerformancePlotter import PerformancePlotter

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
    path = "/atlas/local/pnef/Pileup/Data/"
    tfileData = path+"20140211.16.32_PileUpStudies_DiLeptonSelection_rev351611.VBFH125_ZZ4lep_PowhegPythia8MC12_COMMON.jetmet2012pileupcustom.root"

    # labels 
    LumiSqrtS= "#sqrt{s}=8 TeV"
    jetlabel = "Anti-k_{t} LCW+JES R=0.4";
    etalabel = "|#eta| < 2.4";
    ptlabel  = "p_{T} > 20 GeV, |#eta| < 2.4";
    gaplabel = " ";
    lumiLabel= "#sqrt{s} = 8 TeV, L = 20.3 fb^{-1}"
    normalizationLabel = "MC normalized to data"
    MClabel = "qq' #rightarrow Hqq', H#rightarrowZZ#rightarrow lllll"

    # Histos from Tree
    tree               = "EventTree"
    
    
    #++++++++++++++++ BASIC SETUP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = PerformancePlotter()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Simulation Preliminary")

    plotter.COLORSSolid  = [R.kMagenta+3, R.kViolet+9, R.kTeal-5, R.kBlack, R.kGray+3, R.kGray+2, R.kGray+1, R.kOrange+10]
    plotter.COLORSOpen   = [R.kMagenta+3, R.kViolet+9, R.kTeal-5, R.kBlack, R.kGray+3, R.kGray+2, R.kGray+1, R.kOrange+10]
    plotter.MARKERSData  = [20, 4, 21] 

    
    

    ##########################################################################################################################

    EventSelection = "&& fabs(Jeta[0] - Jeta[1])>3 && Jpt[0]>20 && Jpt[1]>20 && JisHS[0] && JisHS[1]"
    graphsall   = plotter.JetVetoEfficiencyVsNPV(tfileData, tree, [20,30], [0,5,10,15,20,25,30], " && fabs(Jeta)<2.4 && Jpt!=Jpt[0] && Jpt!=Jpt[1]"                                        ,EventSelection, False, "JetVetoEffall")
    graphscut   = plotter.JetVetoEfficiencyVsNPV(tfileData, tree, [20,30], [0,5,10,15,20,25,30], " && fabs(Jeta)<2.4 && Jpt!=Jpt[0] && Jpt!=Jpt[1] &&kNN100trim_pt20to50_Likelihood>0.7",   EventSelection, False ,"JetVetoEffcut")
    graphstruth = plotter.JetVetoEfficiencyVsNPV(tfileData, tree, [20,30], [0,5,10,15,20,25,30], " && fabs(Jeta)<2.4 && Jpt!=Jpt[0] && Jpt!=Jpt[1] &&Jtruthpt>10",                          EventSelection, True  ,"JetVetoEfftruth")

    PTLeg = ["p_{T} > 20 GeV", "p_{T} > 30 GeV"]
    print plotter.H_.keys()
    for graph in graphscut+graphstruth:
        plotter.H_[graph].setopen()
        plotter.H_[graph].AddToLeg(False)

#    plotter.H_[graphstruth[0]].isData()
#    plotter.H_[graphstruth[1]].isData()
#    plotter.H_[graphstruth[2]].isData()
    graphs= graphsall+graphscut
#    graphs= graphsall+graphscut+graphstruth
    for graph in graphs:
        plotter.H_[graph].setdrawstyle("LP")
    plotter.H_[graphs[0]].setdrawstyle("ALP")
    
    print graphs
    plotter.Draw(graphs, "N_{Vtx}", "jet veto efficiency", False, True, [MClabel,jetlabel,"solid markers: inclusive", "open markers: JVT>0.7"],PTLeg)

    

