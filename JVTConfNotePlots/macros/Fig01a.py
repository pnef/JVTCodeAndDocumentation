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
    parser.add_option("--doRhoPU"            ,dest="doRhoPU"           ,    help="calculate PU rho from tracks", action="store_true")
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/skimmed/"
    path = "/atlas/local/pnef/Pileup/Data/"
    fT   =path+"skimmedPt20to50Eta2p4.20140213.00.12_PileUpStudies_rev351681.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN20140213.root"
    
    # Draw Variable ------------------------------------------------------------------------
    putree                = "JetTree0"
    hstree                = "JetTree0"
    treeSelection_PU      = " JisPU &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1  "
    treeSelection_HS      = " JisHS &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 && Jtruthpt >10  "
    selection             = ""

    treevariable          = [  "JstdJVF"                      ]  
    VarName               = [  "JVF"                          ] 
    nBins                 = [   50                            ] 
    binMin                = [  -1                             ] 
    binMax                = [  1.1                            ] 
    HSonTheRight          = [  True                           ]
    tfiles                = [  fT                             ]
    additionalCutPU       = [  ""                             ]
    additionalCutHS       = [  ""                             ]

    ApplyWeights          = True

    HStargetEff        = 0.90
    EffFakeVsVar       = "NPV"
    ptBinsLow          = [20]
    ptBinsHigh         = [30]
    npvBinsLow         = [0  ,15, 30]
    npvBinsHigh        = [15, 30, 45]

    # labels 
    MClabel            = "Pythia8 dijets"
    funkylabel         = ""
    jetlabel           = "Anti-k_{t} LCW+JES R=0.4";
    etalabel           = "|#eta| < 2.4";
    ptLegend           = ", 20 < p_{T} < 30 GeV"
    efflabel           = "Target signal efficiency ="+str(HStargetEff)
    #-------------------------------------------


    
    #++++++++++++++++ END INPUT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = Plotter()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Simulation Preliminary")

    plotter.COLORSSolid  = [ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.COLORSOpen   = [ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.COLORSFilled = [ROOT.kViolet+9, ROOT.kTeal-5,ROOT.kViolet+9, ROOT.kMagenta+3 ]
    plotter.COLORSLine   = [ROOT.kGray+3, ROOT.kRed+1, ROOT.kMagenta+3,ROOT.kTeal-5, ROOT.kBlack,  ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1,]

    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # histo from tree
    for varIndex,var in zip(range(len(treevariable)),treevariable):
        for ptbin,ptlow,pthigh in zip(range(len(ptBinsLow)),ptBinsLow,ptBinsHigh):
            print ptbin,ptlow,pthigh 
            histos =[]
            legends=[]
            for npvbin,npvlow,npvhigh in zip(range(len(npvBinsLow)),npvBinsLow,npvBinsHigh):
                hisths = "HSJVarptbin"+str(ptbin)+"npvbin"+str(npvbin)
                selection_HS = treeSelection_HS + "&& (Jpt >"+str(ptlow)+"&& Jpt<"+str(pthigh)+"&& NPVtruth >"+str(npvlow)+"&& NPVtruth<"+str(npvhigh)+")"
                print  selection_HS
                plotter.GetHistoFromTree(hisths, tfiles[varIndex],  hstree , var, selection_HS, nBins[varIndex], binMin[varIndex], binMax[varIndex])

                histHS = "HSJVarptbin"+str(ptbin)+"npvbin"+str(npvbin)
                if(npvbin==0):
                    plotter.H_[histHS].setdrawstyle("hist")
                if(npvbin==2):
                    plotter.H_[histHS].setopen()

                    
                histos.append(histHS)
                if(len(npvBinsLow)>1):
                    if(npvbin==0):
                        legends.append(str(npvlow)+"   #leq N_{Vtx}^{truth} < "+str(npvhigh))
                    else:
                        legends.append(str(npvlow)+" #leq N_{Vtx}^{truth} < "+str(npvhigh))
                else:
                    legends.append("HS jets")

            ptLegend   = ", "+str(ptlow)+" < p_{T} < "+str(pthigh)+" GeV"
            npvLegend = "0 #leq N_{Vtx}^{truth} #leq 30"
            plotter.Draw(histos, VarName[varIndex], "Normalized Entries", True, True, 
                          [MClabel, jetlabel, etalabel+ptLegend, funkylabel], legends )



    print plotter.H_.keys()
    

