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
    treeSelectionNPVdep_PU= " NPV>=5 && JisPU &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 "
    treeSelectionNPVdep_HS= " NPV>=5 && JisHS &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 && Jtruthpt >10"
    selection             = ""

    treevariable          = [  "JstdJVF"                , "JstdJVF"                , "JstdJVF"                 ]  
    VarName               = [  "20<p_{T}<30 GeV"        , "30<p_{T}<40 GeV"        , "40<p_{T}<50 GeV"         ] 
    nBins                 = [   500                     ,  500                     ,  500                      ] 
    binMin                = [  -1                       , -1                       , -1                        ] 
    binMax                = [  1.1                      , 1.1                      , 1.1                       ] 
    HSonTheRight          = [  True                     , True                     , True                      ]
    tfiles                = [  fT                       , fT                       , fT                        ]
    additionalCutPU       = [  "&&Jpt<30 && Jpt>20"     , "&&Jpt<40 && Jpt>30"     , "&&Jpt<50 && Jpt>40"        ]
    additionalCutHS       = [  "&&Jpt<30 && Jpt>20"     , "&&Jpt<40 && Jpt>30"     , "&&Jpt<50 && Jpt>40"        ]

    ApplyWeights          = False

    EffFakeVsVar       = "NPVtruth"

    # labels 
    MClabel            = "Pythia8 dijets"
#    MClabel            = "Powheg Z#rightarrow #mu#mu"
#    funkylabel         = "Z p_{T} > 30 GeV"
#    funkylabel         = "#DeltaR(#mu,jet)>1"
    funkylabel         = ""
    jetlabel           = "Anti-k_{t} LCW+JES R=0.4";
    etalabel           = "|#eta| < 2.4";
    ptLegend           = ", 20 < p_{T} < 50 GeV"
    efflabel           = "JVF > 0.5"
    #-------------------------------------------


    
    #++++++++++++++++ END INPUT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = Plotter()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Simulation Preliminary")
    if(options.prelim):
        plotter.SetAtlasInfo("Simulation Preliminary")


    


    #-----------------------------------------------------------------------------------------------
    # Efficiency And Fake Rate vs NPV --------------------------------------------------------------

    # Get3D histo from tree
    print ">>>plotter: Getting 3D Histos for doEffFakeVsNPV" 

    H_HS_Var_EffFake=[]
    H_PU_Var_EffFake=[]
    Cuts_Var_EffFake=[]
    for VarIndex,Var in zip(range(len(treevariable)),treevariable):
        SelectionNPVdep_PU = treeSelectionNPVdep_PU+additionalCutPU[VarIndex]
        SelectionNPVdep_HS = treeSelectionNPVdep_HS+additionalCutHS[VarIndex]
        if ApplyWeights :
            SelectionNPVdep_PU = "Weight*(1 &&"+SelectionNPVdep_PU+")"
            SelectionNPVdep_HS = "Weight*(1 &&"+SelectionNPVdep_HS+")"

        print ">>>plotter: Getting 3D Histos for doEffFakeVsNPV"
        print ">>>", SelectionNPVdep_PU
        print ">>>", SelectionNPVdep_HS
        plotter.GetHistoFromTree("PU_Var_EffFake_FromTree"+str(VarIndex), tfiles[VarIndex],  putree , Var+":Jpt:"+EffFakeVsVar, SelectionNPVdep_PU , 50,5,55,50,0,100,nBins[VarIndex],binMin[VarIndex],binMax[VarIndex])
        plotter.GetHistoFromTree("HS_Var_EffFake_FromTree"+str(VarIndex), tfiles[VarIndex],  hstree , Var+":Jpt:"+EffFakeVsVar, SelectionNPVdep_HS , 50,5,55,50,0,100,nBins[VarIndex],binMin[VarIndex],binMax[VarIndex])

        # make efficiency plots
        if(HSonTheRight[VarIndex]):
            plotter.EfficiencyVsNPV("HS_Var_EffFake_FromTree"+str(VarIndex), 2, -99, -99, 0.5,  -99, "HS_Var_EffFake_effVal"+str(VarIndex))
            plotter.EfficiencyVsNPV("PU_Var_EffFake_FromTree"+str(VarIndex), 4, -99, -99, 0.5,  -99, "PU_Var_EffFake_effVal"+str(VarIndex))
        else:
            plotter.EfficiencyVsNPV("HS_Var_EffFake_FromTree"+str(VarIndex), 2, -99, -99, -99, 0.5,  "HS_Var_EffFake_effVal"+str(VarIndex))
            plotter.EfficiencyVsNPV("PU_Var_EffFake_FromTree"+str(VarIndex), 4, -99, -99, -99, 0.5,  "PU_Var_EffFake_effVal"+str(VarIndex))

        # format 
        plotter.SetRange("HS_Var_EffFake_effVal"+str(VarIndex), "x", 0, 52,"HS_Var_EffFake_effVal"+str(VarIndex))
        plotter.SetRange("PU_Var_EffFake_effVal"+str(VarIndex), "x", 0, 52,"PU_Var_EffFake_effVal"+str(VarIndex))
        plotter.H_["HS_Var_EffFake_effVal"+str(VarIndex)].h_.SetMinimum(0.6)
        plotter.H_["HS_Var_EffFake_effVal"+str(VarIndex)].h_.SetMaximum(1.4)

        # append
        H_HS_Var_EffFake.append("HS_Var_EffFake_effVal"+str(VarIndex))
        H_PU_Var_EffFake.append("PU_Var_EffFake_effVal"+str(VarIndex))
        Cuts_Var_EffFake.append(VarName[VarIndex] )


    # format and draw
    plotter.Draw(H_HS_Var_EffFake, 
            "Truth vertex multiplicity N_{Vtx}^{truth}", "Efficiency", False, True, [MClabel, jetlabel, etalabel+ptLegend, efflabel], 
            Cuts_Var_EffFake)
    plotter.Draw(H_PU_Var_EffFake, 
            EffFakeVsVar, "Fake rate",  False, True, [MClabel, jetlabel, etalabel+ptLegend, efflabel], 
            Cuts_Var_EffFake)


    print plotter.H_.keys()
    

