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
    treeSelectionNPVdep_PU= "   NPV>=5 &&JisPU &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 "
    treeSelectionNPVdep_HS= "   NPV>=5 &&JisHS &&  fabs(Jeta)<2.4 && Jpt<50 && Jpt>20 && VtxDzTruth<0.1 && Jtruthpt >10"
    selection             = ""
    treevariable          = [  "JJVF"                  ,"JJVF"                    , "kNN100trim_pt20to50_Likelihood", "kNN100trim_pt20to50_Likelihood"]  
    VarName               = [  "JVF"                   ,"JVF"                     , "JVT"                           , "JVT"                        ] 
    nBins                 = [   500                    , 500                      ,  1000                           ,  1000           ] 
    binMin                = [  -1                      ,-1                        , -0.1                            , -0.1              ] 
    binMax                = [  1.1                     ,1.1                       , 1.1                             , 1.1             ] 
    HSonTheRight          = [  True                    ,True                      , True                            , True            ]
    tfiles                = [  fT                      ,fT                        , fT                              , fT              ]
    additionalCutPU       = [  "&&Jpt<30 && Jpt>20"    ,"&&Jpt<40 && Jpt>30"      , "&&Jpt<30 && Jpt>20"            , "&&Jpt<40 && Jpt>30"              ]
    additionalCutHS       = [  "&&Jpt<30 && Jpt>20"    ,"&&Jpt<40 && Jpt>30"      , "&&Jpt<30 && Jpt>20"            , "&&Jpt<40 && Jpt>30"              ]

    ApplyWeights          = False


    HStargetEff        = 0.95
    EffFakeVsVar       = "Mu"
    ptBinsLow          = [20]
    ptBinsHigh         = [30]
    npvBinsLow         = [5  ,15]
    npvBinsHigh        = [15, 30]
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
    efflabel           = "Target signal efficiency = "+str(HStargetEff)
    #-------------------------------------------


    
    #++++++++++++++++ END INPUT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = Plotter()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Simulation Preliminary")
    if(options.prelim):
        plotter.SetAtlasInfo("Simulation Preliminary")

    #-----------------------------------------------------------------------------------------------
    # Efficiency And Fake Rate vs NPV --------------------------------------------------------------

    if True:
        # Get3D histo from tree
        print ">>>plotter: Getting 3D Histos for doEffFakeVsNPV" 


        H_HS_Var_EffFake=[]
        H_PU_Var_EffFake=[]
        Cuts_Var_EffFake=[]
        for VarIndex,Var in zip(range(len(treevariable)),treevariable):
            print ">>>plotter: Getting 3D Histos for doEffFakeVsNPV" 
            if ApplyWeights :
                CUTS_PU= "Weight*(1 &&"+treeSelectionNPVdep_PU+additionalCutPU[VarIndex]+")"
                CUTS_HS= "Weight*(1 &&"+treeSelectionNPVdep_HS+additionalCutPU[VarIndex]+")"
            else:
                CUTS_PU = treeSelectionNPVdep_PU+additionalCutPU[VarIndex]
                CUTS_HS = treeSelectionNPVdep_HS+additionalCutPU[VarIndex]

            print CUTS_PU
            print CUTS_HS

            plotter.GetHistoFromTree("PU_Var_EffFake_FromTree"+str(VarIndex), tfiles[VarIndex],  putree , Var+":Jpt:"+EffFakeVsVar, CUTS_PU , 45,5,50,50,0,100,nBins[VarIndex],binMin[VarIndex],binMax[VarIndex])
            plotter.GetHistoFromTree("HS_Var_EffFake_FromTree"+str(VarIndex), tfiles[VarIndex],  hstree , Var+":Jpt:"+EffFakeVsVar, CUTS_HS , 45,5,50,50,0,100,nBins[VarIndex],binMin[VarIndex],binMax[VarIndex])

            # determine cut value for specified efficiency
            BinAndCutForEff = plotter.GetBinForEff("HS_Var_EffFake_FromTree"+str(VarIndex), "z",HStargetEff, HSonTheRight[VarIndex])
            print VarName[VarIndex], HStargetEff, "HS eff reached for bin", BinAndCutForEff["bin"], "with cut", BinAndCutForEff["cut"]
            
            # make efficiency plots
            if(HSonTheRight[VarIndex]):
                plotter.EfficiencyVsNPV("HS_Var_EffFake_FromTree"+str(VarIndex), 2, -99, -99, BinAndCutForEff["cut"],  -99, "HS_Var_EffFake_effVal"+str(VarIndex))
                plotter.EfficiencyVsNPV("PU_Var_EffFake_FromTree"+str(VarIndex), 4, -99, -99, BinAndCutForEff["cut"],  -99, "PU_Var_EffFake_effVal"+str(VarIndex))
            else:
                plotter.EfficiencyVsNPV("HS_Var_EffFake_FromTree"+str(VarIndex), 2, -99, -99, -99, BinAndCutForEff["cut"],  "HS_Var_EffFake_effVal"+str(VarIndex))
                plotter.EfficiencyVsNPV("PU_Var_EffFake_FromTree"+str(VarIndex), 4, -99, -99, -99, BinAndCutForEff["cut"],  "PU_Var_EffFake_effVal"+str(VarIndex))

            # format 
            plotter.SetRange("HS_Var_EffFake_effVal"+str(VarIndex), "x", 5, 40,"HS_Var_EffFake_effVal"+str(VarIndex))
            plotter.SetRange("PU_Var_EffFake_effVal"+str(VarIndex), "x", 5, 40,"PU_Var_EffFake_effVal"+str(VarIndex))
            plotter.H_["HS_Var_EffFake_effVal"+str(VarIndex)].h_.SetMinimum(0.6)
            plotter.H_["HS_Var_EffFake_effVal"+str(VarIndex)].h_.SetMaximum(1.3)
            if(VarIndex%2==1):
                plotter.H_["HS_Var_EffFake_effVal"+str(VarIndex)].setopen()

            # append
            H_HS_Var_EffFake.append("HS_Var_EffFake_effVal"+str(VarIndex))
            H_PU_Var_EffFake.append("PU_Var_EffFake_effVal"+str(VarIndex))
            cut = round(round_to_value(BinAndCutForEff["cut"],0.05),2)
            Cuts_Var_EffFake.append(VarName[VarIndex]+">"+str(cut))
            


        # format and draw
        plotter.Draw(H_HS_Var_EffFake, 
                "#mu", "Efficiency", False, True, [MClabel, jetlabel, etalabel, efflabel,"solid markers: 20 < p_{T} < 30 GeV","open markers: 30 < p_{T} < 40 GeV"], 
                Cuts_Var_EffFake)
        plotter.Draw(H_PU_Var_EffFake, 
                EffFakeVsVar, "Fake rate",  False, True, [MClabel, jetlabel, etalabel, efflabel], 
                Cuts_Var_EffFake)

    

    print plotter.H_.keys()
    

