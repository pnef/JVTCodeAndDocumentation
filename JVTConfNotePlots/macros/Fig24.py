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

""" 

    !!!!!!!!!!!!!!!!!!!!!
    First run JetMult_VsPtThreshold_vsEta_GetWorkingPoint.py to decide cut values
    !!!!!!!!!!!!!!!!!!!!!

"""

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
    tree      = "JetTree0"

    # labels 
    pythialabel = "Sherpa Z#rightarrow#mu#mu"
    jetlabel    = "Anti-k_{t} LCW+JES R=0.4, #sqrt{s}=8 TeV";
    npvtruth    = "#LT N_{Vtx}^{truth} #GT = 23"

    # var cut values
    CutVarName       = "#mu"
    MuCutsLow        = 0
    MuCutsHigh       = 50
    NPVtruthCutsLow  = 0
    NPVtruthCutsHigh = 60

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # get histos
    plotter = Plotter()

    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Simulation Preliminary")

    plotter.COLORSSolid  = [ROOT.kMagenta+3, ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.COLORSOpen   = [ROOT.kMagenta+3, ROOT.kMagenta+3, ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.COLORSFilled = [ROOT.kTeal-5,ROOT.kViolet+9, ROOT.kMagenta+3 ]
    plotter.COLORSLine   = [ROOT.kGray+3, ROOT.kRed+1, ROOT.kMagenta+3,ROOT.kTeal-5, ROOT.kBlack,  ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1,]
    plotter.MARKERSSolid = [20, 21, 22, 23, 33, 34, 28, 29] 
    plotter.MARKERSOpen  = [24, 3, 25, 26, 32, 27, 28, 30]
    plotter.SetLogY(True)

    # PU Jet Mult vs Mu  --------------------------------------------
    # no mu weight since drawn as a function of mu
    plotter.SetLogY(False)
    plotter.GetHistoFromTree("NPVM_1",      tfileMC,    tree , "Mu", "(Jindex!=-999 && fabs(Jeta)<2.4 && Jpt>20 && JisPU                       )", 20,0, 40)
    plotter.GetHistoFromTree("NPVM_2",      tfileMC,    tree , "Mu", "(Jindex!=-999 && fabs(Jeta)>2.4 && fabs(Jeta)<4.5 && Jpt>20 && JisPU     )", 20,0, 40)
    plotter.GetHistoFromTree("NPVM",        tfileMC,    tree , "Mu", "(Jindex==-999                                                            )", 20,0, 40)
    plotter.Divide("NPVM_1",    "NPVM", "NJetsNPVM_1"   ,"B")
    plotter.Divide("NPVM_2",    "NPVM", "NJetsNPVM_2"   ,"B")
    plotter.Draw(["NJetsNPVM_1","NJetsNPVM_2"], "#mu","#LT Pileup jet multiplicity #GT", False, True, 
                [pythialabel,jetlabel,"p_{T} > 20 GeV"], ["|#eta|<2.4", "2.4<|#eta|<4.5"])


    #-----------------------------------------------------------------------------------------------
    ## Jet Mult vs NPV from Tree
    # per unit eta
    # weights important since mu inclusive
    plotter.GetHistoFromTree("PU_NPVvsPTvsMu_1",    tfileMC, "JetTree0", "Mu:Jpt:NPVtruth", "1/(2*2.4)* Weight*(Jindex>=0    && fabs(Jeta)<2.4 && JisPU)",                                         50,0,50, 100, 0, 200, 50,0,50)
    plotter.GetHistoFromTree("PU_NPVvsPTvsMu_1JVT", tfileMC, "JetTree0", "Mu:Jpt:NPVtruth", "1/(2*2.4)* Weight*(Jindex>=0    && fabs(Jeta)<2.4 && JisPU && kNN100trim_pt20to50_Likelihood>0.6)",   50,0,50, 100, 0, 200, 50,0,50)
    plotter.GetHistoFromTree("PU_NPVvsPTvsMu_1JVF", tfileMC, "JetTree0", "Mu:Jpt:NPVtruth", "1/(2*2.4)* Weight*(Jindex>=0    && fabs(Jeta)<2.4 && JisPU && JstdJVF>0.5)"                       ,   50,0,50, 100, 0, 200, 50,0,50)
    plotter.GetHistoFromTree("PU_NPVvsPTvsMu_2",    tfileMC, "JetTree0", "Mu:Jpt:NPVtruth", "1/(2*0.8)* Weight*(Jindex>=0    && fabs(Jeta)>2.4 && fabs(Jeta)<3.2 && JisPU)",                       50,0,50, 100, 0, 200, 50,0,50)
    plotter.GetHistoFromTree("PU_NPVvsPTvsMu_3",    tfileMC, "JetTree0", "Mu:Jpt:NPVtruth", "1/(2*1.3)* Weight*(Jindex>=0    && fabs(Jeta)>3.2 && fabs(Jeta)<4.5 && JisPU)",                       50,0,50, 100, 0, 200, 50,0,50)
    plotter.GetHistoFromTree("NPVvsMu",             tfileMC, "JetTree0", "Mu:NPVtruth",     "Weight*(Jindex==-999 )"                                            ,                       50,0,50, 50,0,50)
    plotter.NJetsVsPtThreshold("PU_NPVvsPTvsMu_1"    ,"NPVvsMu",  NPVtruthCutsLow, NPVtruthCutsHigh, 20, 200, MuCutsLow, MuCutsHigh, "PU_NPVvsPTvsMu_NjetsVsPt1")
    plotter.NJetsVsPtThreshold("PU_NPVvsPTvsMu_1JVT" ,"NPVvsMu",  NPVtruthCutsLow, NPVtruthCutsHigh, 20, 200, MuCutsLow, MuCutsHigh, "PU_NPVvsPTvsMu_NjetsVsPt1JVT")
    plotter.NJetsVsPtThreshold("PU_NPVvsPTvsMu_1JVF" ,"NPVvsMu",  NPVtruthCutsLow, NPVtruthCutsHigh, 20, 200, MuCutsLow, MuCutsHigh, "PU_NPVvsPTvsMu_NjetsVsPt1JVF")
    plotter.NJetsVsPtThreshold("PU_NPVvsPTvsMu_2"    ,"NPVvsMu",  NPVtruthCutsLow, NPVtruthCutsHigh, 20, 200, MuCutsLow, MuCutsHigh, "PU_NPVvsPTvsMu_NjetsVsPt2")
    plotter.NJetsVsPtThreshold("PU_NPVvsPTvsMu_3"    ,"NPVvsMu",  NPVtruthCutsLow, NPVtruthCutsHigh, 20, 200, MuCutsLow, MuCutsHigh, "PU_NPVvsPTvsMu_NjetsVsPt3")
    plotter.SetRange("PU_NPVvsPTvsMu_NjetsVsPt1", "x", 20, 60, "PU_NPVvsPTvsMu_NjetsVsPt1")
    plotter.H_["PU_NPVvsPTvsMu_NjetsVsPt1"]   .setdrawstyle("ALP")
    plotter.H_["PU_NPVvsPTvsMu_NjetsVsPt2"]   .setdrawstyle("LP")
    plotter.H_["PU_NPVvsPTvsMu_NjetsVsPt3"]   .setdrawstyle("LP")
    plotter.H_["PU_NPVvsPTvsMu_NjetsVsPt1JVT"].setopen()
    plotter.H_["PU_NPVvsPTvsMu_NjetsVsPt1JVF"].setopen()
    plotter.H_["PU_NPVvsPTvsMu_NjetsVsPt1JVT"].setdrawstyle("LP")
    plotter.H_["PU_NPVvsPTvsMu_NjetsVsPt1JVF"].setdrawstyle("LP")
    plotter.H_["PU_NPVvsPTvsMu_NjetsVsPt1JVF"].h_.SetMinimum(0.000001); plotter.H_["PU_NPVvsPTvsMu_NjetsVsPt1JVF"].h_.SetMaximum(10)
    plotter.Draw(["PU_NPVvsPTvsMu_NjetsVsPt1","PU_NPVvsPTvsMu_NjetsVsPt1JVT","PU_NPVvsPTvsMu_NjetsVsPt1JVF","PU_NPVvsPTvsMu_NjetsVsPt2","PU_NPVvsPTvsMu_NjetsVsPt3"], 
            "p_{T} threshold [GeV]", "#LT Pileup jet multiplicity per unit #eta #GT", False, True, [pythialabel,jetlabel,npvtruth], ["|#eta|<2.4", "|#eta|<2.4, JVT > 0.6", "|#eta|<2.4, JVF > 0.5", "2.4<|#eta|<3.2", "3.2<|#eta|<4.5"])




    print plotter.H_.keys()
    

