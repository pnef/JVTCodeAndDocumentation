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
    #path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/skimmed/"
    path = "/atlas/output/pnef/"
    fT   =path+"skimmedPt20to50Eta2p4.20140213.00.12_PileUpStudies_rev351681.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN20140213.root"
    
    # Draw Variable ------------------------------------------------------------------------
    putree                = "JetTree0"
    hstree                = "JetTree0"
    treeSelection_PU      = " NPV>=5 && JisPU &&  fabs(Jeta)<2.4 && Jpt<30 && Jpt>20 && VtxDzTruth<0.1  "
    treeSelection_HS      = " NPV>=5 && JisHS &&  fabs(Jeta)<2.4 && Jpt<30 && Jpt>20 && VtxDzTruth<0.1 && Jtruthpt >10  "
    selection             = ""
    treevariable          = [    "JnPUTrkCorrJVF" ]  
    VarName               = [    "corrJVF"                        ] 
    nBins                 = [    50                          ] 
    binMin                = [    -1                         ] 
    binMax                = [    1.1                          ] 
    HSonTheRight          = [    True                         ]
    tfiles                = [    fT                           ]
    additionalCutPU       = [    ""                           ]
    additionalCutHS       = [    ""                           ]

    ApplyWeights          = True


    HStargetEff        = 0.9
    EffFakeVsVar       = "NPV"
    ptBinsLow          = [20]
    ptBinsHigh         = [30]
    npvBinsLow         = [5   ]
    npvBinsHigh        = [30]
    etaBinsLow         = [0  ,1]
    etaBinsHigh        = [1, 2.5]

    # labels 
    MClabel            = "Pythia8 dijets"
#    MClabel            = "Powheg Z#rightarrow #mu#mu"
#    funkylabel         = "Z p_{T} > 30 GeV"
#    funkylabel         = "#DeltaR(#mu,jet)>1"
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
    if(options.prelim):
        plotter.SetAtlasInfo("Simulation Preliminary")
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # histo from tree
    if(True):
        for varIndex,var in zip(range(len(treevariable)),treevariable):
            for ptbin,ptlow,pthigh in zip(range(len(ptBinsLow)),ptBinsLow,ptBinsHigh):
                print ptbin,ptlow,pthigh 
                histos =[]
                legends=[]
                for npvbin,npvlow,npvhigh in zip(range(len(npvBinsLow)),npvBinsLow,npvBinsHigh):
                    histpu = "PUJVarptbin"+str(ptbin)+"npvbin"+str(npvbin)
                    hisths = "HSJVarptbin"+str(ptbin)+"npvbin"+str(npvbin)
                    selection_PU = treeSelection_PU + "&& (Jpt >"+str(ptlow)+"&& Jpt<"+str(pthigh)+"&& NPV >"+str(npvlow)+"&& NPV<"+str(npvhigh)+")"
                    selection_HS = treeSelection_HS + "&& (Jpt >"+str(ptlow)+"&& Jpt<"+str(pthigh)+"&& NPV >"+str(npvlow)+"&& NPV<"+str(npvhigh)+")"
                    print selection_PU, selection_HS
                    plotter.GetHistoFromTree(histpu, tfiles[varIndex],  putree , var, selection_PU, nBins[varIndex], binMin[varIndex], binMax[varIndex])
                    plotter.GetHistoFromTree(hisths, tfiles[varIndex],  hstree , var, selection_HS, nBins[varIndex], binMin[varIndex], binMax[varIndex])

                    histHS = "HSJVarptbin"+str(ptbin)+"npvbin"+str(npvbin)
                    histPU = "PUJVarptbin"+str(ptbin)+"npvbin"+str(npvbin)
                    if(npvbin==0):
                        plotter.H_[histHS].setdrawstyle("hist")
                        plotter.H_[histPU].setdrawstyle("hist")

                        
                    histos.append(histPU)
                    histos.append(histHS)
                    if(len(npvBinsLow)>1):
                        legends.append("PU jets, "+str(npvlow)+"<N_{PV}<"+str(npvhigh))
                        legends.append("HS jets, "+str(npvlow)+"<N_{PV}<"+str(npvhigh))
                    else:
                        legends.append("PU jets")
                        legends.append("HS jets")

                    plotter.H_[histHS].h_.Scale(1./plotter.H_[histHS].h_.Integral())
                    plotter.H_[histPU].h_.Scale(1./plotter.H_[histPU].h_.Integral())
                    plotter.H_[histHS].h_.SetMinimum(0.0001)
                    plotter.H_[histHS].h_.SetMaximum(100)
                    plotter.H_[histPU].h_.SetMinimum(0.0001)
                    plotter.H_[histPU].h_.SetMaximum(100)

                ptLegend   = ", "+str(ptlow)+" < p_{T} < "+str(pthigh)+" GeV"
                npvLegend = "0 #leq N_{Vtx} #leq 30"
                plotter.SetLogY(True)
                plotter.Draw(histos, VarName[varIndex], "Normalized Entries", True, True, 
                              [MClabel, jetlabel, etalabel+ptLegend, npvLegend, funkylabel], legends )



    print plotter.H_.keys()
    

