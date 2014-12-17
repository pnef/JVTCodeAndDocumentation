/**************************************************************************
 **
 **   File:         Analysis_EventTreeFiller.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      P. Nef
 **
 **************************************************************************/

#define Analysis_EventTreeFiller_cxx

#include "Analysis_EventTreeFiller.h"
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
 void Analysis_EventTreeFiller::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_EventTreeFiller: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();


  // trees -------------------
  fEventTree = new TTree("EventTree", "Tree with event-by-event variables");
  AddBranches(fEventTree);
} 


///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_EventTreeFiller::ProcessEvent()
{

  if (Debug()) cout << "Analysis_EventTreeFiller: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		            << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();
  
  // Event Selection goes here... ---------------------------
  // ...

  // Fill Tree-----------------------------------------------
  FillTree("AntiKt4LCTopo_ptgt20");


  // end----------------------------------------------------------------------
  if (Debug()) cout << "Analysis_EventTreeFiller: DEBUG End ProcessEvent(): RunNumber = " << RunNumber() 
		            << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;
  return true;
}



///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_EventTreeFiller::WorkerTerminate()
{
    fEventTree->Write();


  // Nothing more

}

///===========================================
/// Fill Tree
///========================================
void Analysis_EventTreeFiller::FillTree(const MomKey JetKey){
  if(Debug()) cout <<"Analysis_EventTreeFiller::FillTree Start" << endl;

  ResetBranches(fEventTree);
  FillEventVars(fEventTree, JetKey);
  fEventTree->Fill();
    
  if(Debug()) cout <<"Analysis_EventTreeFiller::FillTree End" << endl;
  return;
}

///=============================================
/// Add Branches To Tree
///=============================================
void Analysis_EventTreeFiller::AddBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_EventTreeFiller::AddBranches Start" << endl;
  
    // Event Info
    tree->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
    tree->Branch("RunNumber",                 &fTRunNumber,              "RunNumber/I");
    tree->Branch("Weight" ,                   &fTWeight,                 "Weight/F");
    tree->Branch("Mu" ,                       &fTMu,                     "Mu/F");
    tree->Branch("NPVtruth" ,                 &fTNPVtruth,               "NPVtruth/I");
    tree->Branch("NPV" ,                      &fTNPVtruth,               "NPV/I");
  
    // Jet vars --------------------------------------------------------------------------------------
    tree->Branch("NJets",                     &fTNJets,                  "NJets/I");
    tree->Branch("NJetsFilled",               &fTNJetsFilled,            "NJetsFilled/I");
    tree->Branch("Jpt",                       &fTJPt,                    "Jpt[NJetsFilled]");
    tree->Branch("Jeta",                      &fTJEta,                   "Jeta[NJetsFilled]");
    tree->Branch("Jphi",                      &fTJPhi,                   "Jphi[NJetsFilled]");
    tree->Branch("Jm",                        &fTJM,                     "Jm[NJetsFilled]");

  if(Debug()) cout <<"Analysis_EventTreeFiller::AddBranches End" << endl;
    return;
}

void Analysis_EventTreeFiller::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_EventTreeFiller::ResetBranches Start" << endl;
  
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
        fTJPt             [i] = -999.99;
        fTJEta            [i] = -999.99;
        fTJPhi            [i] = -999.99;
        fTJM              [i] = -999.99;
    }

  if(Debug()) cout <<"Analysis_EventTreeFiller::ResetBranches End" << endl;
    return;
}


void Analysis_EventTreeFiller::FillEventVars(TTree *tree, const MomKey JetKey){
  if(Debug()) cout <<"Analysis_EventTreeFiller::FillEventVars Begin" << endl;
  
    // Event Info -----------------------------------
    fTEventNumber                 = Int("EventNumber");
    fTRunNumber                   = Int("RunNumber");
    fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
    fTNPV                         = Exists("NPV")?      Int("NPV"):-1;
    fTWeight                      = DefaultWeight();
    fTMu                          = Float("averageIntPerXing");                  

    // fill jets ----------------------
    for(int iJ=0; iJ<jets(JetKey); ++iJ){
        Particle  *myjet   = &(jet(iJ, JetKey));

        fTNJets++;
        fTNJetsFilled++;
        if(iJ>=MaxNJets)       continue;

        fTJPt  [iJ]       = myjet->p.Pt();
        fTJEta [iJ]       = myjet->p.Eta();
        fTJPhi [iJ]       = myjet->p.Phi();
        fTJM   [iJ]       = myjet->p.M();
    }

  if(Debug()) cout <<"Analysis_EventTreeFiller::FillEventVars End" << endl;
  return;
}

