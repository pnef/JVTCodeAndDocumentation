/**************************************************************************
 **
 **   File:         Analysis_LargeRTreeFiller.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      P. Nef
 **
 **************************************************************************/

#define Analysis_LargeRTreeFiller_cxx

#include "Analysis_LargeRTreeFiller.h"
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
 void Analysis_LargeRTreeFiller::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_LargeRTreeFiller: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();


  // config
  fisMuScan          = false;
  ChainCfg()->Get("isMuScan",           fisMuScan);

  TObjArray* arr = ChainCfg()->String("TRIMMEDJETNAMES").Tokenize(",");
  for(Int_t i = 0; i<arr->GetEntries(); ++i) { 
    TObjString* TypeObj = (TObjString*)arr->At(i);
    TString Type        = TypeObj->GetString();
    cout << Type << endl;
    trimmedJetNames.push_back(Type);
  }
  
  // Specific Config
  fTrkSel = "tracksGoodSel0";
  Cfg()->Get("TrackSel", fTrkSel);

  // trees -------------------
  OutputDir()->cd();
  fJetTree = new TTree("LargeRJetTree", "Tree with LargeR");
  AddBranches(fJetTree);
} 


///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_LargeRTreeFiller::ProcessEvent()
{

  if (Debug()) cout << "Analysis_LargeRTreeFiller: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();
  

  // Event Sel 
  if(! Exists("PileUpStudiesSelection")) return true;
  if(! Bool("PileUpStudiesSelection")  ) return true;



  // Fill Tree-----------------------------------------------------------------------------
  FillTree("Fat10AntiKt");


  // end----------------------------------------------------------------------
  if (Debug()) cout << "Analysis_LargeRTreeFiller: DEBUG End ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;
  return true;
}



///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_LargeRTreeFiller::WorkerTerminate()
{
  gDirectory->cd("/");
    fJetTree->Write();


  // Nothing more

}

///===========================================
/// Fill Tree
///========================================
void Analysis_LargeRTreeFiller::FillTree(const MomKey JetKey){
  if(Debug()) cout <<"Analysis_LargeRTreeFiller::FillTree Start" << endl;
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
    
  if(Debug()) cout <<"Analysis_LargeRTreeFiller::FillTree End" << endl;
 return;

}

///=============================================
/// Add Branches To Tree
///=============================================
void Analysis_LargeRTreeFiller::AddBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_LargeRTreeFiller::AddBranches Start" << endl;
  
    // Event Info
    tree->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
    tree->Branch("RunNumber",                 &fTRunNumber,              "RunNumber/I");
    tree->Branch("Weight" ,                   &fTWeight,                 "Weight/F");
    tree->Branch("NPVtruth" ,                 &fTNPVtruth,               "NPVtruth/F");
    tree->Branch("Mu" ,                       &fTMu,                     "Mu/F");
  
    // Jet vars --------------------------------------------------------------------------------------
    tree->Branch("Jindex",                    &fTJindex,                 "Jindex/I");
    tree->Branch("Jpt"  ,                     &fTJpt,                    "Jpt/F");
    tree->Branch("Jm"   ,                     &fTJm,                     "Jm/F");
    tree->Branch("Jeta" ,                     &fTJeta,                   "Jeta/F");
    tree->Branch("Jphi" ,                     &fTJphi,                   "Jphi/F");
    tree->Branch("JTopdR" ,                   &fTJTopdR,                 "JTopdR[2]/F");
    tree->Branch("JTopc1dR" ,                 &fTJTopc1dR,               "JTopc1dR[2]/F");
    tree->Branch("JTopc2dR" ,                 &fTJTopc2dR,               "JTopc2dR[2]/F");
    tree->Branch("JZdR" ,                     &fTJZdR,                   "JZdR/F");
    tree->Branch("JZc1dR" ,                   &fTJZc1dR,                 "JZc1dR/F");
    tree->Branch("JZc2dR" ,                   &fTJZc2dR,                 "JZc2dR/F");
    tree->Branch("JWdR" ,                     &fTJWdR,                   "JWdR/F");
    tree->Branch("JWc1dR" ,                   &fTJWc1dR,                 "JWc1dR/F");
    tree->Branch("JWc2dR" ,                   &fTJWc2dR,                 "JWc2dR/F");
    tree->Branch("JNSub"  ,                   &fTJNSub,                  "JNSub/I");
    tree->Branch("Jsubpt"  ,                  &fTJsubPt,                 "Jsubpt[JNSub]/F");
    tree->Branch("Jsubm"   ,                  &fTJsubM,                  "Jsubm[JNSub]/F");
    tree->Branch("Jsubeta" ,                  &fTJsubEta,                "Jsubeta[JNSub]/F");
    tree->Branch("Jsubphi" ,                  &fTJsubPhi,                "Jsubphi[JNSub]/F");
    tree->Branch("JsubJVF" ,                  &fTJsubJVF,                "JsubJVF[JNSub]/F");
    tree->Branch("JsubCorrJVF" ,              &fTJsubnPUtrkCorrJVF,      "JsubCorrJVF[JNSub]/F");

    for(int iN=0; iN<trimmedJetNames.size(); ++iN){
        tree->Branch(TString::Format("%s_pt", trimmedJetNames[iN].Data()).Data(), &(fTJTrimmedPt[iN]),   TString::Format("%s_pt/F", trimmedJetNames[iN].Data()).Data());
        tree->Branch(TString::Format("%s_m",  trimmedJetNames[iN].Data()).Data(), &(fTJTrimmedM [iN]),   TString::Format("%s_m/F" , trimmedJetNames[iN].Data()).Data());
    }

  if(Debug()) cout <<"Analysis_LargeRTreeFiller::AddBranches End" << endl;
    return;
}

void Analysis_LargeRTreeFiller::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_LargeRTreeFiller::ResetBranches Start" << endl;
  
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
    fTJm                     = -999.99;
    fTJZdR                   = -999;
    fTJZc1dR                 = -999;
    fTJZc2dR                 = -999;
    fTJWdR                   = -999;
    fTJWc1dR                 = -999;
    fTJWc2dR                 = -999;

    for (int i=0; i<2; ++i){
        fTJTopdR  [i]                 = -999;
        fTJTopc1dR[i]                 = -999;
        fTJTopc2dR[i]                 = -999;
    }

    fTJNSub                  = 0;
    // subjets 
    for(int i=0; i<300; ++i){
        fTJsubPt[i]           = -999;
        fTJsubEta[i]          = -999;
        fTJsubPhi[i]          = -999;
        fTJsubM[i]            = -999;
        fTJsubJVF[i]          = -999;
        fTJsubnPUtrkCorrJVF[i]= -999;
    }
    
    for(int i=0; i<300; ++i){
        fTJTrimmedPt[i]       = -999;
        fTJTrimmedM[i]        = -999;
    }

  if(Debug()) cout <<"Analysis_LargeRTreeFiller::ResetBranches End" << endl;
    return;
}


void Analysis_LargeRTreeFiller::FillEventVars(TTree *tree, const MomKey JetKey){
  if(Debug()) cout <<"Analysis_LargeRTreeFiller::FillEventVars Begin" << endl;
  

    // Event Info
    fTEventNumber                 = Int("EventNumber");
    fTRunNumber                   = Int("RunNumber");
    fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
    fTWeight                      = DefaultWeight();
    fTMu                          = Float("averageIntPerXing");                  

  if(Debug()) cout <<"Analysis_LargeRTreeFiller::FillEventVars End" << endl;
  return;
}

///============================================================
/// Fill jet vars 
///============================================================
void Analysis_LargeRTreeFiller::FillJetVars(TTree* tree, int jindex, const MomKey JetKey){
  if(Debug()) cout <<"Analysis_LargeRTreeFiller::FillJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex, JetKey));

  fTJindex       = jindex;
  fTJpt          = myjet->p.Pt();
  fTJeta         = myjet->p.Eta();
  fTJphi         = myjet->p.Phi();
  fTJm           = myjet->p.M();
  fTJTopdR[0]    = myjet->Exists("recosTop_dR")           ? myjet->Float("recosTop_dR"):-999;
  fTJTopc1dR[0]  = myjet->Exists("recosTop_child1_dR")    ? myjet->Float("recosTop_child1_dR"):-999;
  fTJTopc2dR[0]  = myjet->Exists("recosTop_child2_dR")    ? myjet->Float("recosTop_child2_dR"):-999;
  fTJTopdR[1]    = myjet->Exists("recosAntiTop_dR")       ? myjet->Float("recosAntiTop_dR"):-999;
  fTJTopc1dR[1]  = myjet->Exists("recosAntiTop_child1_dR")? myjet->Float("recosAntiTop_child1_dR"):-999;
  fTJTopc2dR[1]  = myjet->Exists("recosAntiTop_child2_dR")? myjet->Float("recosAntiTop_child2_dR"):-999;
  fTJZdR         = myjet->Exists("recosZ_dR")             ? myjet->Float("recosZ_dR"):-999;
  fTJZc1dR       = myjet->Exists("recosZ_child1_dR")      ? myjet->Float("recosZ_child1_dR"):-999;
  fTJZc2dR       = myjet->Exists("recosZ_child2_dR")      ? myjet->Float("recosZ_child2_dR"):-999;
  fTJWdR         = myjet->Exists("recosW_dR")             ? myjet->Float("recosW_dR"):-999;
  fTJWc1dR       = myjet->Exists("recosW_child1_dR")      ? myjet->Float("recosW_child1_dR"):-999;
  fTJWc2dR       = myjet->Exists("recosW_child2_dR")      ? myjet->Float("recosW_child2_dR"):-999;

  for(int iS=0; iS<myjet->Objs("Subjets"); ++iS){
    if(iS >=300) continue;
    fTJNSub++;
    fTJsubPt [iS]           =  ((Particle*) myjet->Obj("Subjets",iS))->p.Pt();
    fTJsubEta[iS]           =  ((Particle*) myjet->Obj("Subjets",iS))->p.Eta();
    fTJsubPhi[iS]           =  ((Particle*) myjet->Obj("Subjets",iS))->p.Phi();
    fTJsubM  [iS]           =  ((Particle*) myjet->Obj("Subjets",iS))->p.Pt();
    fTJsubJVF[iS]           =  ((Particle*) myjet->Obj("Subjets",iS))->Float("JVFLinks_"+fTrkSel+"Ghost_JVF");
    fTJsubnPUtrkCorrJVF[iS] =  ((Particle*) myjet->Obj("Subjets",iS))->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF");
  }
  
  for(int iN=0; iN<trimmedJetNames.size(); ++iN){
    if(iN >=300) continue;
    Particle *trimmedJ = NULL;
    if(myjet->Exists(trimmedJetNames[iN])) trimmedJ = (Particle*) myjet->Obj(trimmedJetNames[iN]);
    if(trimmedJ!=NULL){
        fTJTrimmedPt[iN] = trimmedJ->p.Pt();
        fTJTrimmedM [iN] = trimmedJ->p.M();
    }
  }

  if(Debug()) cout <<"Analysis_LargeRTreeFiller::FillJetVars End" << endl;
  return;
}
