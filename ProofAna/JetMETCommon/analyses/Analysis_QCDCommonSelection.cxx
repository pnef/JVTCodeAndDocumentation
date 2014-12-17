/**************************************************************************
 **
 **   File:         Analysis_QCDCommonSelection.h
 **
 **   Description:  See header
 **                 
 **   Authors:      M. Swiatlowski, B. Nachman
 **
 **************************************************************************/

#define Analysis_QCDCommonSelection_cxx

#include "Analysis_QCDCommonSelection.h"
#include "AnaConfig.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "TKey.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "PileupReweighting/TPileupReweighting.h"
#include "TileTripReader/TTileTripReader.h" 



///=========================================
/// WorkerBegin: setup binnaing, etc
///=========================================
void Analysis_QCDCommonSelection::WorkerBegin()
{
  if (Debug()) cout << "Analysis_QCDCommonSelection: DEBUG In WorkerBegin()" << endl;

  Analysis_JetMET_Base::WorkerBegin();

  //SetMaps(&m_h1d, &m_h2d);

  m_treader=new Root::TTileTripReader("myTripReader");


  if (Debug()) cout << "Analysis_QCDCommonSelection: DEBUG Finish WorkerBegin()" << endl;

} 

///=========================================
/// WorkerTerminate: Usually nothing
///=========================================
void Analysis_QCDCommonSelection::WorkerTerminate(){

}

///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_QCDCommonSelection::ProcessEvent()
{
  if (Debug()) cout << "Analysis_QCDCommonSelection: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;


  OutputDir()->cd();
  //MakeTruthPlots(); 


  ///////////////////////////////////
  ////      TRIGGER SETUP        ////
  ///////////////////////////////////

 TString mPRWVEC = "PRWVEC";
 const static MomKey Passed("Passed"), j("jets"), LumiFrac("LumiFrac"), LumiFracFull("LumiFracFull");
 const static MomKey AntiKt4LCTopo("AntiKt4LCTopo"), AntiKt4TopoEM("AntiKt4TopoEM"),
  AntiKt6LCTopo("AntiKt6LCTopo"), AntiKt6TopoEM("AntiKt6TopoEM");


  if(Debug()){
    ChainCfg()->Show();
  }
  if(ChainCfg()->Bool("DOJETTRIGGERS")){
      for(Int_t i = 0; i<ChainCfg()->Objs(mPRWVEC); ++i) {
        TString trig = ((TObjString*)ChainCfg()->Obj(mPRWVEC,i))->GetString();

        bool passedAntiKt4LCTopo = false;

        float lumi = ChainCfg()->Float("LUMI");
        float lumifrac = 0.;
        float lumifracfull = 0.;

        if(Bool(trig+Passed)){
          if(Exists(j+AntiKt4LCTopo) && jets(AntiKt4LCTopo) > 0){
            if(trig.Contains("_j15_") && jet(0, AntiKt4LCTopo).p.Perp() > 21.){
              passedAntiKt4LCTopo = true;
              lumifrac = 0.00728262 / lumi;
              lumifracfull = 0.014734 / lumi;
            } else if(trig.Contains("_j25_") && jet(0, AntiKt4LCTopo).p.Perp() > 34.){
              passedAntiKt4LCTopo = true;
              lumifrac = 0.0373692 / lumi;
              lumifracfull = 0.0784406 / lumi;
            } else if(trig.Contains("_j35_") && jet(0, AntiKt4LCTopo).p.Perp() > 45.){
              passedAntiKt4LCTopo = true;
              lumifrac = 0.255253 / lumi;
              lumifracfull = 0.452736 / lumi;
            } else if(trig.Contains("_j45_a4tchad_L2FS_L1J15") && jet(0, AntiKt4LCTopo).p.Perp() > 57.){
              passedAntiKt4LCTopo = true;
              lumifrac = 0.0378231 / lumi;
              lumifracfull = 0.137455 / lumi;
            } else if(trig.Contains("_j45_") && jet(0, AntiKt4LCTopo).p.Perp() > 57.){
              passedAntiKt4LCTopo = true;
              lumifrac = 0.427999 / lumi;
              lumifracfull = 0.426189 / lumi;
            } else if(trig.Contains("_j55_") && jet(0, AntiKt4LCTopo).p.Perp() > 68.){
              passedAntiKt4LCTopo = true;
              lumifrac = 0.210446 / lumi;
              lumifracfull = 0.210446 / lumi;
            } else if(trig.Contains("_j80_") && jet(0, AntiKt4LCTopo).p.Perp() > 97.){
              passedAntiKt4LCTopo = true;
              lumifrac = 1.10322 / lumi;
              lumifracfull = 2.3164 / lumi;
            } else if(trig.Contains("_j110_") && jet(0, AntiKt4LCTopo).p.Perp() > 134.){
              passedAntiKt4LCTopo = true;
              lumifrac = 4.74964 / lumi;
              lumifracfull = 9.81141 / lumi;
            } else if(trig.Contains("_j145_") && jet(0, AntiKt4LCTopo).p.Perp() > 177.){
              passedAntiKt4LCTopo = true;
              lumifrac = 17.8494 / lumi;
              lumifracfull = 36.2647 / lumi;
            } else if(trig.Contains("_j180_") && jet(0, AntiKt4LCTopo).p.Perp() > 207.){
              passedAntiKt4LCTopo = true;
              lumifrac = 30.1438 / lumi;
              lumifracfull = 78.7553 / lumi;
            } else if(trig.Contains("_j220_") && jet(0, AntiKt4LCTopo).p.Perp() > 247.){
              passedAntiKt4LCTopo = true;
              lumifrac = 126.352 / lumi;
              lumifracfull = 261.379 / lumi;
            } else if(trig.Contains("_j360_") && jet(0, AntiKt4LCTopo).p.Perp() > 410.){
              passedAntiKt4LCTopo = true;
              lumifrac = 5905.84 / lumi;
              lumifracfull = 20100. / lumi;
            } else if(trig.Contains("_j460_") && jet(0, AntiKt4LCTopo).p.Perp() > 523.){
              passedAntiKt4LCTopo = true;
              lumifrac = 5905.84 / lumi;
              lumifracfull = 20100. / lumi;
            }
          } // end check on jets
        } // end check if trigger passed
        Set(trig+Passed+AntiKt4LCTopo, passedAntiKt4LCTopo);
        Set(trig+Passed+AntiKt4LCTopo+LumiFrac, lumifrac);
        Set(trig+Passed+AntiKt4LCTopo+LumiFracFull, lumifracfull);

        //delete type;
      }
  }

  if(Debug()){
    cout << " At end of QCDCommonSelection: " << endl;
    Show();
  }

  ////////////////////////////////////////////
  ////             CUTFLOW                ////
  ////////////////////////////////////////////

  bool all = true;

  // -1: GRL
  if(!isMC() && !Bool("grl")){
    Set("QCDSelection_GRL", false);
    all = false;
  } else {
    Set("QCDSelection_GRL", true);
    if(all) BookCutflowNoPrefix("QCDSelection_GRL");
  }

  // 0. larerror
 if(Int("larerror") > 1 || Int("tileerror") == 2 || (Int("coreflags")&0x40000) !=0 || m_treader->checkEvent(RunNumber(),Float("LBN"),EventNumber()) == 0) {
    Set("QCDSelection_LarError", false);
    all = false;
  } else {
    Set("QCDSelection_LarError", true);
    if(all) BookCutflowNoPrefix("QCDSelection_LarError");
  }

  // 1. Good vertex
  if(vtxs() < 1){
    Set("QCDSelection_Vtxs", false);
    all = false;
  } else {
    Set("QCDSelection_Vtxs", true);
    if(all) BookCutflowNoPrefix("QCDSelection_Vtxs");
  }

  // 2. No badLoose jets in AntiKt4LCTopo
  bool hasBadLoose = false;
  for(int iJ = 0; iJ < jets("AntiKt4LCTopo"); iJ++){
    if(jet(iJ, "AntiKt4LCTopo").Bool("isBadLoose")){
      hasBadLoose = true;
      break;
    }
  }
  if(hasBadLoose){
    all = false;
    Set("QCDSelection_NoBadLooseJets", false);
  } else {
    Set("QCDSelection_NoBadLooseJets", true);
    if(all) BookCutflowNoPrefix("QCDSelection_NoBadLooseJets");
  }

  // All done
  Set("QCDSelection_All", all);

  return kTRUE; 

}

