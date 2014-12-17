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
    path = "/atlas/local/pnef/Pileup/Data/"
    fT   =path+"skimmedPt20to50Eta2p4.20140213.00.12_PileUpStudies_rev351681.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN20140213.root"
    
    # Draw Variable ------------------------------------------------------------------------
    putree                = "JetTree0"
    hstree                = "JetTree0"
    treeSelection_PU      = "  NPV>=5 && JisPU &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1  "
    treeSelection_HS      = "  NPV>=5 && JisHS &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 && Jtruthpt >10  "
    treeSelectionROC_PU   = "  NPV>=5 && JisPU &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 "
    treeSelectionROC_HS   = "  NPV>=5 && JisHS &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 && Jtruthpt >10"
    treeSelectionNPVdep_PU= "  NPV>=5 && JisPU &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 "
    treeSelectionNPVdep_HS= "  NPV>=5 && JisHS &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 && Jtruthpt >10"
    selection             = ""
    treevariable          = [ "kNN100trim_pt20to50_Likelihood", "kNN100trim_pt20to50_Likelihood","kNN100trim_pt20to50_Likelihood",]  
    VarName               = [ "JVT(q)"                        , "JVT(g)"                        ,"JVT(b)"                        ,] 
    nBins                 = [ 100                             , 100                             ,100                             ,] 
    binMin                = [ -0.2                            , -0.2                            ,-0.2                            ,] 
    binMax                = [ 1.2                             , 1.2                             ,1.2                             ,] 
    HSonTheRight          = [ True                            , True                            ,True                            ,]
    tfiles                = [ fT                              , fT                              ,fT                              ,]
    additionalCutPU       = [ "                "              , "              "                ,"             "                 ,]
    additionalCutHS       = [ "&& JisParton_uds"              , "&& JisParton_g"                ,"&&JisParton_b"                 ,]

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
    ptLegend           = ", 20 < p_{T} < 50 GeV"
    efflabel           = "Target signal efficiency ="+str(HStargetEff)
    #-------------------------------------------


    
    #++++++++++++++++ END INPUT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = Plotter()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Simulation Preliminary")
    if(options.prelim):
        plotter.SetAtlasInfo("Simulation Preliminary")
    
    plotter.COLORSSolid  = [ROOT.kMagenta+3, ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.COLORSOpen   = [ROOT.kGray+3, ROOT.kMagenta+3, ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.COLORSFilled = [ROOT.kTeal-5,ROOT.kViolet+9, ROOT.kMagenta+3 ]
    plotter.COLORSLine   = [ROOT.kGray+3, ROOT.kRed+1, ROOT.kMagenta+3,ROOT.kTeal-5, ROOT.kBlack,  ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1,]


    #-----------------------------------------------------------------------------------------------
    # ROC curve 
    if True:



        listROC=[]
        # Get3D histo from tree
        for VarIndex,Var in zip(range(len(treevariable)),treevariable):
            SelectionROC_PU = treeSelectionROC_PU+additionalCutPU[VarIndex]
            SelectionROC_HS = treeSelectionROC_HS+additionalCutHS[VarIndex]
            if ApplyWeights :
                SelectionROC_PU = "Weight*(1 &&"+SelectionROC_PU+")"
                SelectionROC_HS = "Weight*(1 &&"+SelectionROC_HS+")"


                
            print ">>>plotter: Getting 3D Histos for ROC, variable", Var
            print ">>>plotter: selection HS", SelectionROC_HS
            print ">>>plotter: selection PU", SelectionROC_PU

            plotter.GetHistoFromTree("PU_Var_ROC_FromTree_"+VarName[VarIndex], tfiles[VarIndex],  putree , Var+":Jpt:NPV", SelectionROC_PU , 51,-0.5,50.5,100,0,500,nBins[VarIndex],binMin[VarIndex],binMax[VarIndex])
            plotter.GetHistoFromTree("HS_Var_ROC_FromTree_"+VarName[VarIndex], tfiles[VarIndex],  hstree , Var+":Jpt:NPV", SelectionROC_HS , 51,-0.5,50.5,100,0,500,nBins[VarIndex],binMin[VarIndex],binMax[VarIndex])
            
            plotter.ROC("HS_Var_ROC_FromTree_"+VarName[VarIndex], "PU_Var_ROC_FromTree_"+VarName[VarIndex], -99, -99, -99, -99, -99, -99, "Var_ROC_fromTree_"+VarName[VarIndex], HSonTheRight[VarIndex])

            print "PU_Var_ROC_FromTree_"+VarName[VarIndex], plotter.H_["PU_Var_ROC_FromTree_"+VarName[VarIndex]].h_.GetEntries()
            print "HS_Var_ROC_FromTree_"+VarName[VarIndex], plotter.H_["HS_Var_ROC_FromTree_"+VarName[VarIndex]].h_.GetEntries()

            if VarIndex==0:
                plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].setdrawstyle("ALP")
            else:
                plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].setdrawstyle("LP")
            if VarName[VarIndex] == "JVF":
                plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].setdrawstyle("LP")
                plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].setopen()

            plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].h_.SetMinimum(0.003)
            plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].h_.SetMaximum(0.3)
            plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].h_.GetXaxis().SetRangeUser(0.78,1)
            listROC.append("Var_ROC_fromTree_"+VarName[VarIndex])

        if plotter.H_[listROC[0]].drawstyle_.find("A")==-1:
            plotter.H_[listROC[0]].drawstyle_ +="A"

        # draw
        ptlabel  = "20 < p_{T} < 50 GeV"
        plotter.SetLogY(True)
        plotter.Draw(listROC,"Efficiency", "Fake Rate", False, True, [MClabel, jetlabel, etalabel, ptlabel ], VarName)

    print plotter.H_.keys()
    

