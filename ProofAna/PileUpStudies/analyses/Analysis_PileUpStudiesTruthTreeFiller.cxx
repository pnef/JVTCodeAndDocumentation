/**************************************************************************
 **
 **   File:         Analysis_PileUpStudiesTruthTreeFiller.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      P. Nef
 **
 **************************************************************************/

#define Analysis_PileUpStudiesTruthTreeFiller_cxx

#include "Analysis_PileUpStudiesTruthTreeFiller.h"
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
 void Analysis_PileUpStudiesTruthTreeFiller::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_PileUpStudiesTruthTreeFiller: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();


  // config
  fisMuScan          = false;
  ChainCfg()->Get("isMuScan",           fisMuScan);

  // trees -------------------
  OutputDir()->cd();
  fJetTree = new TTree("TruthJetTree", "Tree with truthJets");
  AddBranches(fJetTree);
} 


///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_PileUpStudiesTruthTreeFiller::ProcessEvent()
{

  if (Debug()) cout << "Analysis_PileUpStudiesTruthTreeFiller: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();
  

  if(! Exists("PileUpStudiesSelection")) return true;
  if(! Bool("PileUpStudiesSelection")  ) return true;


  // Fill Tree-----------------------------------------------------------------------------
  FillTree("AntiKt4Truth");


  // end----------------------------------------------------------------------
  if (Debug()) cout << "Analysis_PileUpStudiesTruthTreeFiller: DEBUG End ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;
  return true;
}



///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_PileUpStudiesTruthTreeFiller::WorkerTerminate()
{
  gDirectory->cd("/");
    fJetTree->Write();


  // Nothing more

}

///===========================================
/// Fill Tree
///========================================
void Analysis_PileUpStudiesTruthTreeFiller::FillTree(const MomKey JetKey){
  if(Debug()) cout <<"Analysis_PileUpStudiesTruthTreeFiller::FillTree Start" << endl;
  ResetBranches(fJetTree);
  FillEventVars(fJetTree, JetKey);
  fJetTree->Fill();
  for(int iJet=0; iJet<jets(JetKey); ++iJet){
    if(jet(iJet,JetKey).p.Pt()<20) continue; // cut for muscan samples
    ResetBranches(fJetTree);
    FillEventVars(fJetTree, JetKey);
    FillJetVars(fJetTree, iJet, JetKey);

    fJetTree->Fill();
  }
    
  if(Debug()) cout <<"Analysis_PileUpStudiesTruthTreeFiller::FillTree End" << endl;
 return;

}

///=============================================
/// Add Branches To Tree
///=============================================
void Analysis_PileUpStudiesTruthTreeFiller::AddBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_PileUpStudiesTruthTreeFiller::AddBranches Start" << endl;
  
    // Event Info
    tree->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
    tree->Branch("RunNumber",                 &fTRunNumber,              "RunNumber/I");
    tree->Branch("Weight" ,                   &fTWeight,                 "Weight/F");
    tree->Branch("NPVtruth" ,                 &fTNPVtruth,               "NPVtruth/F");
    tree->Branch("Mu" ,                       &fTMu,                     "Mu/F");
  
    // Jet vars --------------------------------------------------------------------------------------
    tree->Branch("Jindex",                    &fTJindex,                 "Jindex/I");
    tree->Branch("Jpt"  ,                     &fTJpt,                    "Jpt/F");
    tree->Branch("Jeta" ,                     &fTJeta,                   "Jeta/F");
    tree->Branch("Jphi" ,                     &fTJphi,                   "Jphi/F");
    tree->Branch("JmatchConstscalePt"  ,      &fTJmatchConstitPt,        "JmatchConstscalePt/F");
    tree->Branch("JmatchAreacorrPt"  ,        &fTJmatchAreacorrPt,       "JmatchAreacorrPt/F");
    tree->Branch("JmatchPt"  ,                &fTJmatchPt,               "JmatchPt/F");
    tree->Branch("JmatchJVF" ,                &fTJmatchJVF,              "JmatchJVF/F");
    tree->Branch("JmatchcorrJVF" ,            &fTJmatchcorrJVF,          "JmatchcorrJVF/F");
    tree->Branch("JmatchRpT" ,                &fTJmatchRpT,              "JmatchRpT/F");

  if(Debug()) cout <<"Analysis_PileUpStudiesTruthTreeFiller::AddBranches End" << endl;
    return;
}

void Analysis_PileUpStudiesTruthTreeFiller::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_PileUpStudiesTruthTreeFiller::ResetBranches Start" << endl;
  
    // Event Info
    fTEventNumber           = -9999;
    fTRunNumber             = -9999;
    fTWeight                = -999.99;
    fTNPVtruth              = -99;
    fTMu                    = -99;
    
    // jet vars
    fTJindex                 = -999;
    fTJpt                    = -999.99;         
    fTJeta                   = -999.99;
    fTJphi                   = -999.99;
    fTJmatchConstitPt        = -999.99;
    fTJmatchAreacorrPt       = -999.99;
    fTJmatchPt               = -999.99;         
    fTJmatchJVF              = -999.99;         
    fTJmatchcorrJVF          = -999.99;         
    fTJmatchRpT              = -999.99;         

  if(Debug()) cout <<"Analysis_PileUpStudiesTruthTreeFiller::ResetBranches End" << endl;
    return;
}


void Analysis_PileUpStudiesTruthTreeFiller::FillEventVars(TTree *tree, const MomKey JetKey){
  if(Debug()) cout <<"Analysis_PileUpStudiesTruthTreeFiller::FillEventVars Begin" << endl;
  

    // Event Info
    fTEventNumber                 = Int("EventNumber");
    fTRunNumber                   = Int("RunNumber");
    fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
    fTWeight                      = DefaultWeight();
    fTMu                          = Float("averageIntPerXing");                  

  if(Debug()) cout <<"Analysis_PileUpStudiesTruthTreeFiller::FillEventVars End" << endl;
  return;
}

///============================================================
/// Fill jet vars 
///============================================================
void Analysis_PileUpStudiesTruthTreeFiller::FillJetVars(TTree* tree, int jindex, const MomKey JetKey){
  if(Debug()) cout <<"Analysis_PileUpStudiesTruthTreeFiller::FillJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex, JetKey));

  fTJindex  = jindex;
  fTJpt     = myjet->p.Pt();
  fTJeta    = myjet->p.Eta();
  fTJphi    = myjet->p.Phi();
  if(myjet->Exists("AntiKt4LCTopo_match")){
    Particle *matchedJet = (Particle*) myjet->Obj("AntiKt4LCTopo_match");
    fTJmatchConstitPt    = matchedJet->Float("constscale_pt");
    fTJmatchAreacorrPt   = matchedJet->Float("areacorr_pt");
    fTJmatchPt           = matchedJet->p.Pt();
    fTJmatchJVF          = matchedJet->Exists("JVFLinks_tracksGoodBenGhost_JVF")             ?matchedJet->Float("JVFLinks_tracksGoodBenGhost_JVF")          :-1;
    fTJmatchcorrJVF      = matchedJet->Exists("JVFLinks_tracksGoodBenGhost_nPUTrkCorrJVF")   ?matchedJet->Float("JVFLinks_tracksGoodBenGhost_nPUTrkCorrJVF"):-1;
    fTJmatchRpT          = matchedJet->Exists("JVFLinks_tracksGoodBenGhost_HSPVtrkSumOverPt")?matchedJet->Float("JVFLinks_tracksGoodBenGhost_HSPVtrkSumOverPt"):-1; }


  if(Debug()) cout <<"Analysis_PileUpStudiesTruthTreeFiller::FillJetVars End" << endl;
  return;
}
