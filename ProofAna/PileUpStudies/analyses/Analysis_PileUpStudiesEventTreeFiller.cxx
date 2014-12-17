/**************************************************************************
 **
 **   File:         Analysis_PileUpStudiesEventTreeFiller.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      P. Nef
 **
 **************************************************************************/

#define Analysis_PileUpStudiesEventTreeFiller_cxx

#include "Analysis_PileUpStudiesEventTreeFiller.h"
#include "AnaConfig.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>      // std::setprecision
#include <cstdlib>
#include <sstream>
#include "TKey.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"


///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_PileUpStudiesEventTreeFiller::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_PileUpStudiesEventTreeFiller: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();
  
  fTrkSel   = "tracksGoodSel0";
  Cfg()->Get("TrackSel", fTrkSel);


  // trees -------------------
  OutputDir()->cd();
  fEventTree = new TTree("EventTree", "Tree with event-by-event variables");
  AddBranches(fEventTree);
} 


///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_PileUpStudiesEventTreeFiller::ProcessEvent()
{

  if (Debug()) cout << "Analysis_PileUpStudiesEventTreeFiller: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();
  

  if(! Exists("PileUpStudiesSelection")) return true;
  if(! Bool("PileUpStudiesSelection")  ) return true;

  // Trk definition ------------------------------------------------------------------------
  MomKey myTracksGhost(fTrkSel+"Ghost");
  MomKey myJVFLinks("JVFLinks_"+fTrkSel+"Ghost");

  // Fill Tree-----------------------------------------------------------------------------
  FillTree("AntiKt4LCTopo_ptgt20", myJVFLinks, myTracksGhost);


  // end----------------------------------------------------------------------
  if (Debug()) cout << "Analysis_PileUpStudiesEventTreeFiller: DEBUG End ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;
  return true;
}



///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_PileUpStudiesEventTreeFiller::WorkerTerminate()
{
  gDirectory->cd("/");
    fEventTree->Write();


  // Nothing more

}

///===========================================
/// Fill Tree
///========================================
void Analysis_PileUpStudiesEventTreeFiller::FillTree(const MomKey JetKey, const MomKey JVFKey, const MomKey TrkKey){
  if(Debug()) cout <<"Analysis_PileUpStudiesEventTreeFiller::FillTree Start" << endl;

  ResetBranches(fEventTree);
  FillEventVars(fEventTree, JetKey, JVFKey, TrkKey);
  fEventTree->Fill();
    
  if(Debug()) cout <<"Analysis_PileUpStudiesEventTreeFiller::FillTree End" << endl;
  return;
}

///=============================================
/// Add Branches To Tree
///=============================================
void Analysis_PileUpStudiesEventTreeFiller::AddBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_PileUpStudiesEventTreeFiller::AddBranches Start" << endl;
  
    // Event Info
    tree->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
    tree->Branch("RunNumber",                 &fTRunNumber,              "RunNumber/I");
    tree->Branch("Weight" ,                   &fTWeight,                 "Weight/F");
    tree->Branch("Mu" ,                       &fTMu,                     "Mu/F");
    tree->Branch("NPVtruth" ,                 &fTNPVtruth,               "NPVtruth/I");
    tree->Branch("NPV" ,                      &fTNPV,                    "NPV/I");
  
    // Jet vars --------------------------------------------------------------------------------------
    tree->Branch("NJets",                     &fTNJets,                  "NJets/I");
    tree->Branch("NJetsFilled",               &fTNJetsFilled,            "NJetsFilled/I");
    tree->Branch("Jpt",                       &fTJPt,                    "Jpt[NJetsFilled]/F");
    tree->Branch("Jeta",                      &fTJEta,                   "Jeta[NJetsFilled]/F");
    tree->Branch("Jphi",                      &fTJPhi,                   "Jphi[NJetsFilled]/F");
    tree->Branch("JisPU",                     &fTJisPU,                  "JisPU[NJetsFilled]/I");
    tree->Branch("JisHS",                     &fTJisHS,                  "JisHS[NJetsFilled]/I");
    tree->Branch("JisBtagged",                &fTJisBtagged,             "JisBtagged[NJetsFilled]/I");
    tree->Branch("JJVF",                      &fTJJVF,                   "JJVF[NJetsFilled]/F");
    tree->Branch("JRpT",                      &fTJRpT,                   "JRpT[NJetsFilled]/F");
    tree->Branch("JnPUTrkCorrJVF",            &fTJnPUTrkCorrJVF,         "JnPUTrkCorrJVF[NJetsFilled]/F");
    tree->Branch("Jtruthpt",                  &fTJtruthpt,               "Jtruthpt[NJetsFilled]/F");

  if(Debug()) cout <<"Analysis_PileUpStudiesEventTreeFiller::AddBranches End" << endl;
    return;
}

void Analysis_PileUpStudiesEventTreeFiller::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_PileUpStudiesEventTreeFiller::ResetBranches Start" << endl;
  
    // Event Info
    fTEventNumber           = -999;
    fTRunNumber             = -999;
    fTWeight                = -999;
    fTNPVtruth              = -999;
    fTNPV                   = -999;
    fTMu                    = -999;
    
    fTNJets                 = 0;
    fTNJetsFilled           = 0;
    for(int i=0;i<MaxNJets; ++i){
        fTJPt             [i] = -999;
        fTJEta            [i] = -999;
        fTJPhi            [i] = -999;
        fTJM              [i] = -999;
        fTJisPU           [i] = -999; 
        fTJisHS           [i] = -999;
        fTJisBtagged      [i] = -999;
        fTJJVF            [i] = -999;
        fTJRpT            [i] = -999;
        fTJnPUTrkCorrJVF  [i] = -999;
        fTJtruthpt        [i] = -999;
    }

  if(Debug()) cout <<"Analysis_PileUpStudiesEventTreeFiller::ResetBranches End" << endl;
    return;
}


void Analysis_PileUpStudiesEventTreeFiller::FillEventVars(TTree *tree, const MomKey JetKey, const MomKey JVFKey , const MomKey TrkKey){
  if(Debug()) cout <<"Analysis_PileUpStudiesEventTreeFiller::FillEventVars Begin" << endl;
  
    // Event Info -----------------------------------
    fTEventNumber                 = Int("EventNumber");
    fTRunNumber                   = Int("RunNumber");
    fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
    fTNPV                         = Exists("NPV")?      Int("NPV"):-1;
    fTWeight                      = Float("EventWeight");
    fTMu                          = Float("averageIntPerXing");                  

    // fill jets ----------------------
    for(int iJ=0; iJ<jets(JetKey); ++iJ){
        Particle  *myjet   = &(jet(iJ, JetKey));

        fTNJets++;
        fTNJetsFilled++;
        if(iJ>=MaxNJets)       continue;

        fTJPt           [iJ] = myjet->p.Pt();
        fTJEta          [iJ] = myjet->p.Eta();
        fTJPhi          [iJ] = myjet->p.Phi();
        fTJM            [iJ] = myjet->p.M();
        fTJisPU         [iJ] = myjet->Bool("isPUJet");
        fTJisHS         [iJ] = myjet->Bool("isHSJet");
        fTJisBtagged    [iJ] = myjet->Exists("btag")? myjet->Float("btag"): -1;
        fTJJVF          [iJ] = myjet->Float(JVFKey+"_JVF");
        fTJRpT          [iJ] = myjet->Float(JVFKey+"_HSPVtrkSumOverPt")>0 ? myjet->Float(JVFKey+"_HSPVtrkSumOverPt"):0;
        fTJnPUTrkCorrJVF[iJ] = myjet->Float(JVFKey+"_nPUTrkCorrJVF");

        if(myjet->Bool("isHSJet") && myjet->Exists("AntiKt4Truth_match")){
            Particle *truth = (Particle*) myjet->Obj("AntiKt4Truth_match");
            fTJtruthpt[iJ]= truth->p.Pt();                                                                                         
        }
    }

  if(Debug()) cout <<"Analysis_PileUpStudiesEventTreeFiller::FillEventVars End" << endl;
  return;
}

