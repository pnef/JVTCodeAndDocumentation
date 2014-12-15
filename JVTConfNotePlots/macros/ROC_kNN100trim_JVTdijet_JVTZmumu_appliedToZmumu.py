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
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fTnew         ="/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/20140214.22.39_PileUpStudies_DiLeptonSelection_rev351830.Zmumu_PowhegPythia8_MC12_COMMON.jetmet2012pileupcustom.root"
    
    # Draw Variable ------------------------------------------------------------------------
    treeSelectionROC_PU   = " NPV>=5 && JisPU &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1   "
    treeSelectionROC_HS   = " NPV>=5 && JisHS &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 &&  Jtruthpt >10"
    selection             = ""

    treevariable          = [     "kNN100trim_pt20to50_Likelihood"              ,"kNN100trim_pt20to50_Likelihood_Zmumu"                   ]
    VarName               = [     "JVT(dijets #rightarrow Z(#mu#mu))"           ,"JVT(Z(#mu#mu) #rightarrow Z(#mu#mu))" ] 
    nBins                 = [     500                                           ,500                                                      ] 
    binMin                = [     -1                                            ,-1                                                       ] 
    binMax                = [     1.1                                           ,1.1                                                      ] 
    HSonTheRight          = [     True                                          ,True                                                     ]
    tfiles                = [     fTnew                                         ,fTnew                                                    ]
    additionalCutPU       = [     ""                                            ,""                                                       ]
    additionalCutHS       = [     ""                                            ,""                                                       ]
    TreeName              = [     "JetTree0"                                    ,"JetTree0"                                                ]


    ApplyWeights          = True

    HStargetEff        = 0.9
    EffFakeVsVar       = "NPV"
    ptBinsLow          = [20,20,30,40]
    ptBinsHigh         = [50,30,40,50]
    npvBinsLow         = [5  ,15]
    npvBinsHigh        = [15, 30]
    etaBinsLow         = [0  ,1]
    etaBinsHigh        = [1, 2.5]

    # labels 
    MClabel            = "Pythia8 dijets (2#rightarrow2)"
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

    plotter.COLORSSolid  = [ROOT.kMagenta+3, ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.COLORSOpen   = [ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # histo from tree

    #-----------------------------------------------------------------------------------------------
    # ROC curve 
    if options.doROC:



        # Get3D histo from tree
        for ptlow, pthigh in zip(ptBinsLow,ptBinsHigh):
            listROC=[]
            for VarIndex,Var in zip(range(len(treevariable)),treevariable):
                SelectionROC_PU = treeSelectionROC_PU+additionalCutPU[VarIndex]+"&&Jpt>"+str(ptlow)+"&&Jpt<"+str(pthigh)
                SelectionROC_HS = treeSelectionROC_HS+additionalCutHS[VarIndex]+"&&Jpt>"+str(ptlow)+"&&Jpt<"+str(pthigh)
                if ApplyWeights :
                    SelectionROC_PU = "Weight*(1 &&"+SelectionROC_PU+")"
                    SelectionROC_HS = "Weight*(1 &&"+SelectionROC_HS+")"


                    
                print ">>>plotter: Getting 3D Histos for ROC, variable", Var
                print ">>>plotter: selection HS", SelectionROC_HS
                print ">>>plotter: selection PU", SelectionROC_PU

                plotter.GetHistoFromTree("PU_Var_ROC_FromTree_"+VarName[VarIndex], tfiles[VarIndex],  TreeName[VarIndex] , Var+":Jpt:NPV", SelectionROC_PU , 51,-0.5,50.5,100,0,500,nBins[VarIndex],binMin[VarIndex],binMax[VarIndex])
                plotter.GetHistoFromTree("HS_Var_ROC_FromTree_"+VarName[VarIndex], tfiles[VarIndex],  TreeName[VarIndex] , Var+":Jpt:NPV", SelectionROC_HS , 51,-0.5,50.5,100,0,500,nBins[VarIndex],binMin[VarIndex],binMax[VarIndex])
                
                plotter.ROC("HS_Var_ROC_FromTree_"+VarName[VarIndex], "PU_Var_ROC_FromTree_"+VarName[VarIndex], -99, -99, -99, -99, -99, -99, "Var_ROC_fromTree_"+VarName[VarIndex], HSonTheRight[VarIndex])

                print "PU_Var_ROC_FromTree_"+VarName[VarIndex], plotter.H_["PU_Var_ROC_FromTree_"+VarName[VarIndex]].h_.GetEntries()
                print "HS_Var_ROC_FromTree_"+VarName[VarIndex], plotter.H_["HS_Var_ROC_FromTree_"+VarName[VarIndex]].h_.GetEntries()

                if VarIndex==0:
                    plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].setdrawstyle("ALP")
                else:
                    plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].setdrawstyle("LP")
                listROC.append("Var_ROC_fromTree_"+VarName[VarIndex])
                if VarIndex%2==1:
                     plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].setopen()


            # draw
            ptlabel  = str(ptlow)+" < p_{T} < "+str(pthigh)+" GeV"
            plotter.Draw(listROC,"Efficiency", "Fake Rate", False, True, [jetlabel, etalabel, ptlabel ], VarName)

    print plotter.H_.keys()
    

