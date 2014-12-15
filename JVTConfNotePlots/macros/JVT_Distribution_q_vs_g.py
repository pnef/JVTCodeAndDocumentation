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
    parser.add_option("--doRhoPU"            ,dest="doRhoPU"           ,    help="calculate PU rho from tracks", action="store_true")
    parser.add_option("--prelim"             ,dest="prelim"            ,    help="Add ATLAS Preliminary Label", action="store_true")
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/skimmed/"
    fT   =path+"skimmedPt20to50Eta2p4.20140213.00.12_PileUpStudies_rev351681.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN20140213.root"
    
    # Draw Variable ------------------------------------------------------------------------
    putree                = "JetTree0"
    hstree                = "JetTree0"
    treeSelection_HS      = " NPV>=5 && JisHS &&  fabs(Jeta)<2.4 && Jpt<30 && Jpt>20 && VtxDzTruth<0.1 && Jtruthpt >10  "
    selection             = ""
    treevariable          = [    "kNN100trim_pt20to50_Likelihood" ]  
    VarName               = [    "JVT"                            ] 
    nBins                 = [    30                               ] 
    binMin                = [    -0.2                             ] 
    binMax                = [    1.2                              ] 
    HSonTheRight          = [    True                             ]
    tfiles                = [    fT                               ]
    additionalCutPU       = [    ""                               ]
    additionalCutHS       = [    ""                               ]

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

    
    plotter.COLORSSolid  = [ROOT.kViolet+9, ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.COLORSFilled = [ROOT.kViolet+9, ROOT.kTeal-5,   ROOT.kViolet+9, ROOT.kMagenta+3 ]
    plotter.COLORSLine   = [ ROOT.kGray+3,  ROOT.kGray+3,  ROOT.kGray+3]
    plotter.STYLESFilled = [3004 , 3002]

    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # histo from tree
    if(True):
        for varIndex,var in zip(range(len(treevariable)),treevariable):
            for ptbin,ptlow,pthigh in zip(range(len(ptBinsLow)),ptBinsLow,ptBinsHigh):
                print ptbin,ptlow,pthigh 
                histos =[]
                legends=[]
                for npvbin,npvlow,npvhigh in zip(range(len(npvBinsLow)),npvBinsLow,npvBinsHigh):
                    histHS_q = "HS_q_JVarptbin"+str(ptbin)+"npvbin"+str(npvbin)
                    histHS_g = "HS_g_JVarptbin"+str(ptbin)+"npvbin"+str(npvbin)
                    selection_HS_q = "("+treeSelection_HS + "&& (Jpt >"+str(ptlow)+"&& Jpt<"+str(pthigh)+"&& NPV >"+str(npvlow)+"&& NPV<"+str(npvhigh)+" && JisParton_uds))"
                    selection_HS_g = "("+treeSelection_HS + "&& (Jpt >"+str(ptlow)+"&& Jpt<"+str(pthigh)+"&& NPV >"+str(npvlow)+"&& NPV<"+str(npvhigh)+" && JisParton_g))"
                    print selection_HS_q 
                    print selection_HS_g
                    plotter.GetHistoFromTree(histHS_q, tfiles[varIndex],  putree , var, selection_HS_q, nBins[varIndex], binMin[varIndex], binMax[varIndex])
                    plotter.GetHistoFromTree(histHS_g, tfiles[varIndex],  hstree , var, selection_HS_g, nBins[varIndex], binMin[varIndex], binMax[varIndex])

                    if(varIndex==0):
                        plotter.H_[histHS_q].setdrawstyle("hist")

                        
                    histos.append(histHS_q)
                    histos.append(histHS_g)
                    if(len(npvBinsLow)>1):
                        legends.append("q jets, "+str(npvlow)+"<N_{PV}<"+str(npvhigh))
                        legends.append("g jets, "+str(npvlow)+"<N_{PV}<"+str(npvhigh))
                    else:
                        legends.append("q jets")
                        legends.append("g jets")

                    plotter.H_[histHS_q].h_.Scale(1./plotter.H_[histHS_q].h_.Integral())
                    plotter.H_[histHS_g].h_.Scale(1./plotter.H_[histHS_g].h_.Integral())
                    plotter.H_[histHS_q].h_.SetMinimum(0.001)
                    plotter.H_[histHS_q].h_.SetMaximum(10)
                    plotter.H_[histHS_g].h_.SetMinimum(0.001)
                    plotter.H_[histHS_g].h_.SetMaximum(10)

                ptLegend   = ", "+str(ptlow)+" < p_{T} < "+str(pthigh)+" GeV"
                npvLegend = "0 #leq N_{Vtx} #leq 30"
                plotter.SetLogY(True)
                plotter.Draw(histos, VarName[varIndex], "Normalized Entries", True, True, 
                              [MClabel, jetlabel, etalabel+ptLegend, npvLegend, funkylabel], legends )



    print plotter.H_.keys()
    

