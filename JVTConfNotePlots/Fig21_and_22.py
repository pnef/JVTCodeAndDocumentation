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
    parser.add_option("--drawFromTree"       ,dest="drawFromTree"      ,    help="draw from Tree", action="store_true")
    parser.add_option("--approval1"          ,dest="approval1"       , action="store_true")
    parser.add_option("--approval2"          ,dest="approval2"       , action="store_true")
    parser.add_option("--ariel"              ,dest="ariel"       , action="store_true")
    parser.add_option("--fCutOpt"            ,dest="fCutOpt"           , action="store_true")
    parser.add_option("--TwoDplot"             ,dest="TwoDplot"      , action="store_true")
    parser.add_option("--prelim"             ,dest="prelim"            ,    help="Add ATLAS Preliminary Label", action="store_true")
    
    global options, args, files, RootFile, TDir
    (options, args) = parser.parse_args()
	

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++ USER INPUT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    path = "/Users/pascal/Data/ATLAS/ProofAnaOutput/PileUpStudies/LargeR/"
    fT   =path+"20140409.15.49_LargeR_vtx_fixedCleansingII.Pythia8_Wprime_WZqqqq_m1000_COMMON.jetmet2012pileupcustom.root"
    
    # Draw Variable ------------------------------------------------------------------------
    
    if (options.approval1):
        tree                  = "LargeRJetTree"
        treeSelection         = "JZdR<0.5 && JZdR>0 && Jpt>300 && fabs(Jeta)<1.5&& JZc1dR<1 && JZc2dR<1"
        treevariable          = [ "Jm"         ,"jet_0p00TrimCorrJVF_m"     ,"jet_0p04TrimCorrJVF_m"     ,"jet_0p10TrimCorrJVF_m"]  
        VarName               = [ "Mass [GeV]" ,"Mass [GeV]"                ,"Mass [GeV]"                ,"Mass [GeV]"           ] 
        VarLeg                = [ "ungroomed"  ,"corrJVF>0.6"               ,"corrJVF>0.6,f_{cut}=0.04"  ,"corrJVF>0.6,f_{cut}=0.10" ] 
        VarCut                = [ 0            ,0.01                        ,0.04                        ,0.10                   ] 
        nBins                 = [  100          ,100                        ,100                         ,100                    ] 
        binMin                = [   0          ,0                           ,0                           ,0                      ] 
        binMax                = [ 200          ,200                         ,200                         ,200                    ] 
        tfiles                = [ fT           ,fT                          ,fT                          ,fT                     ] 
        
    
    elif (options.approval2):
        tree                  = "LargeRJetTree"
        treeSelection         = "JZdR<0.5 && JZdR>0 && Jpt>300 && fabs(Jeta)<1.5&& JZc1dR<1 && JZc2dR<1"
        treevariable          = [ "Jm"         , "jet_0p05Trim_m"   , "jet_0p04TrimCorrJVF_m"   ,"jet_LinCleansed0p55Trim0p04_m"  ,"jet_JVFCleansedTrim0p04_m"             ]  
        VarName               = [ "Mass [GeV]" , "Mass [GeV]"       , "Mass [GeV]"              ,"Mass [GeV]"                     ,"Mass [GeV]"                            ] 
        VarLeg                = [ "ungroomed"  , "f_{cut}=0.05"     , "f_{cut}=0.04,corrJVF>0.6","f_{cut}=0.04,linCleansing"      ,"f_{cut}=0.04,JVFCleansing"             ] 
        nBins                 = [  100         , 100                , 100                       ,100                              ,100                                     ] 
        binMin                = [   0          , 0                  , 0                         ,0                                ,0                                       ] 
        binMax                = [ 200          , 200                , 200                       ,200                              ,200                                     ] 
        tfiles                = [ fT           , fT                 , fT                        ,fT                               ,fT                                      ] 
    
    elif (options.ariel):
        tree                  = "LargeRJetTree"
        treeSelection         = "JZdR<0.5 && JZdR>0 && Jpt>300 && fabs(Jeta)<1.5&& JZc1dR<1 && JZc2dR<1"
        treevariable          = [ "Jm"         , "jet_0p00TrimCorrJVF_m"      , "jet_0p05TrimCorrJVF_m"      , "jet_0p05Trim_m"       ]  
        VarName               = [ "Mass [GeV]" , "Mass [GeV]"                 , "Mass [GeV]"                 , "Mass [GeV]"           ] 
        VarLeg                = [ "ungroomed"  , "CorrJVF>0.6"                , "CorrJVF>0.6,f_{Cut}=0.05"   , "f_{Cut}=0.05"         ] 
        nBins                 = [  100         , 100                          , 100                          , 100                    ] 
        binMin                = [   0          , 0                            , 0                            , 0                      ] 
        binMax                = [ 200          , 200                          , 200                          , 200                    ] 
        tfiles                = [ fT           , fT                           , fT                           , fT                     ] 
    
    elif (options.fCutOpt):
        tree                  = "LargeRJetTree"
        treeSelection         = "JZdR<0.5 && JZdR>0 && Jpt>300 && fabs(Jeta)<1.5&& JZc1dR<1 && JZc2dR<1"
        treevariable=[]; VarName=[]; VarLeg=[]; VarCut=[]; nBins=[]; binMin=[]; binMax=[]; tfiles=[]


        linCleans                  ="jet_LinCleansed0p55Trim0p"
        treevariable.append(       [ linCleans+"00_m" ,linCleans+"01_m" ,linCleans+"02_m" ,linCleans+"03_m" ,linCleans+"04_m" ,linCleans+"05_m",linCleans+"06_m" ,linCleans+"07_m",linCleans+"08_m"   ]) 
        VarName     .append(       [ "Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"    ,"Mass [GeV]"     ,"Mass [GeV]"    ,"Mass [GeV]"       ]) 
        VarCut      .append(       [ 0.00             ,0.01             ,0.02             ,0.03             ,0.04             ,0.05            ,0.06             ,0.07            ,0.08               ]) 
        tfiles      .append(       [ fT               ,fT               ,fT               ,fT               ,fT               ,fT              ,fT               ,fT              ,fT                 ]) 
        
        JVFCleans                  ="jet_JVFCleansedTrim0p"
        treevariable.append(       [ JVFCleans+"00_m" ,JVFCleans+"01_m" ,JVFCleans+"02_m" ,JVFCleans+"03_m" ,JVFCleans+"04_m" ,JVFCleans+"05_m",JVFCleans+"06_m" ,JVFCleans+"07_m"  ,JVFCleans+"08_m"  ]) 
        VarName     .append(       [ "Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"    ,"Mass [GeV]"     ,"Mass [GeV]"      ,"Mass [GeV]"     ]) 
        VarCut      .append(       [ 0.00             ,0.01             ,0.02             ,0.03             ,0.04             ,0.05            ,0.06             ,0.07              ,0.08             ]) 
        binMax      .append(       [ 140              ,140              ,140              ,140              ,140              ,140             ,140              ,140               ,140              ]) 
        tfiles      .append(       [ fT               ,fT               ,fT               ,fT               ,fT               ,fT              ,fT               ,fT                ,fT               ]) 

        CorrJVF                   ="TrimCorrJVF_m"
        treevariable.append(      ["jet_0p00"+CorrJVF,"jet_0p01"+CorrJVF,"jet_0p02"+CorrJVF,"jet_0p03"+CorrJVF,"jet_0p04"+CorrJVF,"jet_0p05"+CorrJVF,"jet_0p06"+CorrJVF,"jet_0p07"+CorrJVF,"jet_0p08"+CorrJVF ])
        VarName     .append(      ["Mass [GeV]"      ,"Mass [GeV]"      ,"Mass [GeV]"      ,"Mass [GeV]"      ,"Mass [GeV]"      ,"Mass [GeV]"      ,"Mass [GeV]"      ,"Mass [GeV]"      ,"Mass [GeV]"       ]) 
        VarCut      .append(      [0.00              ,0.01              ,0.02              ,0.03              ,0.04              ,0.05              ,0.06              ,0.07              , 0.08              ]) 
        tfiles      .append(      [fT                ,fT                ,fT                ,fT                ,fT                ,fT                ,fT                ,fT                , fT                ]) 
        
        treevariable.append(       [ "jet_0p01Trim_m" ,"jet_0p02Trim_m" ,"jet_0p03Trim_m" ,"jet_0p04Trim_m" ,"jet_0p05Trim_m","jet_0p06Trim_m" ,"jet_0p07Trim_m" ,"jet_0p08Trim_m"     ])
        VarName     .append(       [ "Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"    ,"Mass [GeV]"     ,"Mass [GeV]"     ,"Mass [GeV]"          ]) 
        VarCut      .append(       [ 0.01             ,0.02             ,0.03             ,0.04             ,0.05            ,0.06             ,0.07             ,0.08                  ]) 
        tfiles      .append(       [ fT               ,fT               ,fT               ,fT               ,fT              ,fT               ,fT               ,fT                    ]) 

    # labels 
    MClabel            = "Pythia8 (W'#rightarrow WZ#rightarrow qqqq)"
    masslabel          = "M_{W'} = 1 TeV"
    jetlabel           = "Anti-k_{t} LCW R=1.0";
    #-------------------------------------------


    
    #++++++++++++++++ END INPUT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    plotter = Plotter()
    plotter.SetATLASStyle()
    plotter.SetAtlasInfo("Simulation Preliminary")
    if(options.prelim):
        plotter.SetAtlasInfo("Simulation Preliminary")

    
    plotter.COLORSFilled = [ROOT.kViolet+9, ROOT.kTeal-5, ROOT.kMagenta+3 ]
    plotter.STYLESFilled = [3004, 3002, 3002]
    if(options.approval2):
        plotter.COLORSSolid  = [ROOT.kMagenta+3,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
        plotter.COLORSOpen   = [ROOT.kViolet+9,ROOT.kTeal-5, ROOT.kBlack, ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1, ROOT.kOrange+10]
        plotter.COLORSLine   = [ROOT.kGray+3, ROOT.kRed+1, ROOT.kMagenta+3,ROOT.kTeal-5, ROOT.kBlack,  ROOT.kGray+3, ROOT.kGray+2, ROOT.kGray+1,]
        plotter.MARKERSSolid = [20, 22, 23, 33, 34, 28, 29] 
        plotter.MARKERSOpen  = [25, 26, 32, 27, 28, 30]

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # histo from tree
    if(options.drawFromTree):
        histos =[]
        legends=[]
        for varIndex,var in zip(range(len(treevariable)),treevariable):
            hist = "Var"+str(varIndex)
            plotter.GetHistoFromTree(hist, tfiles[varIndex],  tree , var, treeSelection, nBins[varIndex], binMin[varIndex], binMax[varIndex])
                    
            if(varIndex==0 ) :
                plotter.H_[hist].setdrawstyle("hist")
            else:
                pass

            if (options.approval1 or options.approval2):
                plotter.H_[hist].h_.Scale(1./plotter.H_[hist].h_.Integral())
                plotter.H_[hist].h_.SetMinimum(0)
                plotter.H_[hist].h_.SetMaximum(0.12)

            if (options.approval2 and varIndex==2):
                plotter.H_[hist].setopen()

            legends.append(VarLeg[varIndex])
            histos.append(hist)

        plotter.Draw(histos, VarName[varIndex], "Normalized Entries", True, True, 
                [MClabel, masslabel, jetlabel,"0 #leq #mu #leq 40"], legends )



    if(options.fCutOpt):
        for index, varlist in zip(range(len(treevariable)),treevariable):
            histos =[]
            cuts   =[]
            for varIndex,var in zip(range(len(varlist)),varlist):
                hist = "Var"+str(varIndex)
                plotter.GetHistoFromTree(hist, tfiles[index][varIndex],  tree , var, treeSelection, 1000, 0, 500)
                cuts   .append(VarCut[index][varIndex])
                histos.append(hist)
#            plotter.GetFractionWithinBoundsGraph(histos, cuts, "graph"+str(index), 60, 120)
#            plotter.GetSigmaOverMeanGraph(histos, cuts, "graph"+str(index), 60, 120)
            plotter.GetPeakPlusMinusRangeGraph(histos, cuts, "graph"+str(index), 10)
            plotter.H_["graph"+str(index)].setdrawstyle("ALP")
#            plotter.Draw(["graph"+str(index)])
        
        plotter.H_["graph0"].setdrawstyle("ALP")
        plotter.H_["graph1"].setdrawstyle("LP")
        plotter.H_["graph2"].setdrawstyle("LP")
        plotter.H_["graph3"].setdrawstyle("LP")
        plotter.Draw(["graph0","graph1","graph2","graph3"], "f_{Cut}", "fraction of jets under peak pm 10 GeV", False, True, [MClabel,masslabel,jetlabel], 
                     ["linear cleansing", "JVF cleansing", "trimming & corrJVF>0.6","trimming"])


    if(options.TwoDplot):
        plotter.Set2DCanvas(True,False)
        plotter.SetATLASStyle()
        plotter.SetAtlasInfo("Simulation Preliminary")
        plotter.GetHistoFromTree("2Dplot", fT, "LargeRJetTree", "TMath::Log10(Jsubpt/Jpt):JsubCorrJVF", "JZdR<0.5 && JZdR>0 && Jpt>300 && JZc1dR<1 && JZc2dR<1 && abs(Jeta)<1.5 &&JsubCorrJVF>=0", 50,0,1.1,50,-2.5,1)
        plotter.H_["2Dplot"].setdrawstyle("colz")
        plotter.H_["2Dplot"].h_.SetZTitle("Normalized Entries")
        plotter.Draw(["2Dplot"],"corrJVF", "log_{10}(p_{T}^{subj}/p_{T}^{ungroomed})",True,True, [MClabel,masslabel,jetlabel+" jet","k_{t} LCW R=0.3 subjets"])
#
