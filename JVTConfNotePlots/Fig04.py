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
    parser.add_option("--doRhoPU"            ,dest="doRhoPU"           ,    help="calculate PU rho from tracks", action="store_true")
    parser.add_option("--prelim"             ,dest="prelim"            ,    help="Add ATLAS Preliminary Label", action="store_true")
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/skimmed/"
    fT   =path+"skimmedPt20to50Eta2p4.20140213.00.12_PileUpStudies_rev351681.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN20140213.root"

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    treeSelection_PU      = " NPV>=5 && JisPU &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1  "
    treeSelection_HS      = " NPV>=5 && JisHS &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 && Jtruthpt >10  "

    # labels 
    MClabel            = "Pythia8 dijets"
    jetlabel           = "Anti-k_{t} LCW+JES R=0.4";
    etalabel           = "|#eta| < 2.4";
    ptLegend           = ", 20 < p_{T} < 50 GeV"



    plotter = Plotter()
    plotter.SetAtlasInfo("Simulation Preliminary")
    if(options.prelim):
        plotter.SetAtlasInfo("Simulation Preliminary")

    plotter.Set2DCanvas(True, False)
    plotter.SetATLASStyle()

    # CorrJVF vs RpT
    plotter.GetHistoFromTree("CorrJVFvsRpT", fT,  "JetTree0", "JHSPVtrkSumOverPt:JnPUTrkCorrJVF" , treeSelection_HS+"&&JHSPVtrkSumOverPt<=1.5" , 30, 0, 1.1, 40, 0, 2.0 )
    plotter.H_["CorrJVFvsRpT"].setdrawstyle("colz")
    plotter.H_["CorrJVFvsRpT"].h_.Scale(1./plotter.H_["CorrJVFvsRpT"].h_.Integral())
    plotter.H_["CorrJVFvsRpT"].h_.SetMinimum(0.0000001)
    plotter.H_["CorrJVFvsRpT"].h_.SetMaximum(1)
    plotter.Draw(["CorrJVFvsRpT"], "corrJVF", "R_{pT}", True, True, [MClabel,jetlabel,etalabel+ptLegend,"HS jets"])

    plotter.GetHistoFromTree("CorrJVFvsRpT", fT,  "JetTree0", "JHSPVtrkSumOverPt:JnPUTrkCorrJVF" , treeSelection_PU+"&&JHSPVtrkSumOverPt<=1.5" , 30, 0, 1.1, 40, 0, 2. )
    plotter.H_["CorrJVFvsRpT"].setdrawstyle("colz")
    plotter.H_["CorrJVFvsRpT"].h_.Scale(1./plotter.H_["CorrJVFvsRpT"].h_.Integral())
    plotter.H_["CorrJVFvsRpT"].h_.SetMinimum(0.0000001)
    plotter.H_["CorrJVFvsRpT"].h_.SetMaximum(1)
    plotter.Draw(["CorrJVFvsRpT"], "corrJVF", "R_{pT}", True, True, [MClabel,jetlabel,etalabel+ptLegend,"PU jets"])




    print plotter.H_.keys()
    

