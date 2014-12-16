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
    parser.add_option("--doQvsG"             ,dest="doQvsG"            ,    action="store_true")
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/PUID/skimmed/"
    path = "/atlas/output/pnef/"
    fT_z0      =path+"skimmedPt20to50Eta2p4.20140204.19.26_PileUpStudies_rev349447.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.root" # trk sel0 (JetTree0), dzsinTheta 1mm requirement, no btrack recovery
    fT_vtx     =path+"skimmedPt20to50Eta2p4.20140204.15.11_PileUpStudies_rev349439.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.root" # trk sel0 (JetTree0), TrkFromVertex, btrack recovery 3mm 
    fT_vtx_nob = path+"skimmed_Jpt20to50Eta2p4.20140128.16.47_PileUpStudies_rev348911.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.root"  # trk sel (JetTree), TrkFromVertex
    
    # Draw Variable ------------------------------------------------------------------------
    treeSelectionROC_PU   = "  NPV>=5 && JisPU &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 "
    treeSelectionROC_HS   = "  NPV>=5 && JisHS &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 && Jtruthpt >10"

    treevariable          = [ "JHSPVtrkSumOverPt>0?JHSPVtrkSumOverPt:0"     ,"JHSPVtrkSumOverPt>0?JHSPVtrkSumOverPt:0"      ,"JHSPVtrkSumOverPt>0?JHSPVtrkSumOverPt:0"   ]
    VarName               = [ "R_{pT}, Trk-to-Vtx 1"                        ,"R_{pT}, Trk-to-Vtx 2"                         ,"R_{pT}, Trk-to-Vtx 3"  ] 
    nBins                 = [  200                                          , 200                                           , 200               ] 
    binMin                = [ -1                                            ,-1                                             ,-1                 ] 
    binMax                = [ 1.5                                           ,1.5                                            ,1.5                ] 
    HSonTheRight          = [ True                                          ,True                                           ,True               ] 
    tfiles                = [ fT_z0                                         ,fT_vtx_nob                                     ,fT_vtx             ] 
    additionalCutPU       = [ ""                                            ,""                                             ,""                 ] 
    additionalCutHS       = [ "             "                               ,"              "                               ,"               "  ]
    TreeName              = [  "JetTree0"                                   ,"JetTree"                                      ,"JetTree0"         ]
    ApplyWeights          = True
    
    HStargetEff        = 0.9
    EffFakeVsVar       = "NPV"
    ptBinsLow          = [20,30,40,20]
    ptBinsHigh         = [50,40,50,30]
    npvBinsLow         = [5  ,15]
    npvBinsHigh        = [15, 30]
    etaBinsLow         = [0  ,1]
    etaBinsHigh        = [1, 2.5]

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
    plotter.COLORSSolid  = [ROOT.kMagenta+3,ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.COLORSOpen   = [ROOT.kMagenta+3,ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
    plotter.COLORSFilled = [ROOT.kViolet+9,ROOT.kTeal-5,ROOT.kViolet+9, ROOT.kMagenta+3 ]
    plotter.ROC_no_1_1_point_ = True


    #-----------------------------------------------------------------------------------------------
    # ROC curve 
    if True:
        for ptlow, pthigh in zip(ptBinsLow,ptBinsHigh):
            print "pT bin" , ptlow ,pthigh 
            listROC=[]
            # Get3D histo from tree
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
                plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].h_.SetMinimum(0.001)
                plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].h_.SetMaximum(10)
                plotter.H_["Var_ROC_fromTree_"+VarName[VarIndex]].h_.GetXaxis().SetRangeUser(0.8,1.01)

            # draw
            ptlabel  = str(ptlow)+" < p_{T} < "+str(pthigh)+" GeV"
            plotter.SetLogY(True)
            plotter.Draw(listROC,"Efficiency", "Fake Rate", False, True, [MClabel, jetlabel, etalabel, ptlabel ], VarName)

    print plotter.H_.keys()
    

