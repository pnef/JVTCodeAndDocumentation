#!/usr/bin/python
##################################################################
# pascal nef                             April 23rd, 2013        #
##################################################################
import sys
import os
import re
import math
import time
from optparse import OptionParser
from sys import argv,exit
from plottingBase import *
import ROOT
from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut
from ROOT import TMVA
import numpy as np
import array

##################################################################
# Main ------------------------------------------------------------------------
if __name__ == "__main__":

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # parse arguments --------------------------------------------
    usage = """
    usage: %prog options
    	
    """
    parser = OptionParser(usage)
    
    parser.add_option("--verbose"            ,dest="verbose"           ,    help="set verbose"          , action="store_true")
    parser.add_option("--makeLikelihood"     ,dest="makeLikelihood"    ,    help=""                    , action="store_true")
    parser.add_option("--evaluate"           ,dest="evaluate"          ,    help=""                    , action="store_true")
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    
    #++++++++++++++++ END INPUT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    xmlfilepath            = "/atlas/local/pnef/Pileup/TMVA-xml/"
    xmlfile                = xmlfilepath+"skimmedPt20to50Eta2p4.20140213.00.12_PileUpStudies_rev351681.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN100.nPUtrkCorrJVF_RpT.JHStrkPtSumGE0_KNN100trim.weights.xml"
    LikelihoodHistoFile    = "kNN100trim_pt20to50_Likelihood_Histo.root" 
    LikelihoodHistoName    = "kNN100trim_pt20to50_Likelihood" 
    LikelihoodVarX         = "JnPUTrkCorrJVF"
    LikelihoodVarY         = "JHSPVtrkSumOverPt"
    EvalVarName            = "kNN100trim_pt20to50_Likelihood"
    TMVAMethodName         = "KNN100trim"
    VarListNames           = ["JnPUTrkCorrJVF", "JHSPVtrkSumOverPt"]
    SpectListNames         = ["NPV","Mu", "Jpt", "Jeta", "JisPU", "JisHS", "Weight", "VtxDzTruth", "Jtruthpt", "JJVF", "JCorrPUtrkPtSumOverPt", "nPUTracks", "JHStrkPtSum" ]
    LikelihoodHistoNbins   = [101,       101]
    LikelihoodHistoAxisMins= [-0.005, -0.005]
    LikelihoodHistoAxisMaxs= [1.005,   1.005]

    # TTree where new discriminator is stored
    path                   = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/skimmed/"
    tF                     = path+"skimmedPt20to50Eta2p4.20140213.00.12_PileUpStudies_rev351681.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN20140213.root"
    treeName               = "JetTree0"

    #---------------------------------------

    # labels 
    MClabel            = "Pythia8 dijets (2#rightarrow2)"
    funkylabel         = ""
    jetlabel           = "Anti-k_{t} LCW+JES R=0.4";
    etalabel           = "|#eta| < 2.4";
    ptLegend           = ", 20 < p_{T} < 50 GeV"

    #++++++++++++++++++++++++++++++
    plotter = Plotter()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Internal Simulation")


    if options.makeLikelihood:
        plotter.MakeDiscriminantFromXML(xmlfile, TMVAMethodName, LikelihoodHistoName, VarListNames, SpectListNames, LikelihoodHistoNbins, LikelihoodHistoAxisMins, LikelihoodHistoAxisMaxs,EvalVarName  )
        plotter.H_[EvalVarName].h_.SetZTitle("JVT Likelihood")
        plotter.Draw([EvalVarName], "corrJVF", "R_{pT}", False, True, [MClabel, jetlabel, etalabel+ptLegend], [])
        SaveHist(plotter.H_[EvalVarName], LikelihoodHistoFile)

    if options.evaluate:
        # read likelihood file
        plotter.GetHisto(LikelihoodHistoName, LikelihoodHistoName, LikelihoodHistoFile)

        # read Ttree
        Tfile = ROOT.TFile(tF, "update")
        Tfile.cd()
        tree=ROOT.gDirectory.Get(treeName)
        if tree == None:
            raise NameError('cannot find TTree '+treeName)

        # declare vars
        newVar              = array.array('f',[0])
        JHSPVtrkSumOverPt   = array.array('f',[0])
        JnPUTrkCorrJVF      = array.array('f',[0])

        #make new branch
        newBranch = tree.Branch(EvalVarName, newVar, EvalVarName+"/F")

        # set branch address
        tree.SetBranchAddress("JHSPVtrkSumOverPt",JHSPVtrkSumOverPt);
        tree.SetBranchAddress("JnPUTrkCorrJVF"   ,JnPUTrkCorrJVF);

        # loop over entries
        nentries = tree.GetEntries();
        for entry in range(nentries):
            tree.GetEntry(entry)
            # RpT>1
            R_pT = 0
            if(JHSPVtrkSumOverPt[0]>1): 
                R_pT=1 
            else:
                R_pT=JHSPVtrkSumOverPt[0]

            # no tracks:
            if(JnPUTrkCorrJVF[0]==-1):
                newVar[0] = -0.1
            # no HS tracks:
            elif(JnPUTrkCorrJVF[0]==0):
                newVar[0] = 0.0
            else:
                binnr      = plotter.H_[LikelihoodHistoName].h_.FindBin(JnPUTrkCorrJVF[0], R_pT)
                newVar[0]  = plotter.H_[LikelihoodHistoName].h_.GetBinContent(binnr)
#                print JnPUTrkCorrJVF[0], R_pT, newVar[0], binnr


            newBranch.Fill()

        # write tree    
        tree.Print()
        tree.Write()

    print plotter.H_.keys()
    

