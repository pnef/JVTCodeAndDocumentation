/**************************************************************************
 **
 **   File:         Analysis_EventTreeFiller_ClustersAndTruth.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      P. Nef
 **
 **************************************************************************/

#define Analysis_EventTreeFiller_ClustersAndTruth_cxx

#include "Analysis_EventTreeFiller_ClustersAndTruth.h"
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
 void Analysis_EventTreeFiller_ClustersAndTruth::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_EventTreeFiller_ClustersAndTruth: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();


  // trees -------------------
  fEventTree = new TTree("EventTree", "Tree with event-by-event variables");
  AddBranches(fEventTree);
} 


///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_EventTreeFiller_ClustersAndTruth::ProcessEvent()
{

  if (Debug()) cout << "Analysis_EventTreeFiller_ClustersAndTruth: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		            << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();
  
  // Event Selection goes here... ---------------------------
  // ...

  AddStableParticles();
  // Fill Tree-----------------------------------------------
  FillTree();


  MakeJets(fastjet::antikt_algorithm, 0.4, "clustersLCTopo", "test");
  MakeJets(fastjet::antikt_algorithm, 0.4, "truthsStable",   "test");
  cout << "EventNumber" <<  Int("EventNumber") << endl;
  for(int ij=0; ij<jets("AntiKt4LCTopotest"); ++ij){
//    cout << "reco jet " << ij << " pt " << jet(ij, "AntiKt4LCTopotest").p.Pt() << endl;
  }
  for(int ij=0; ij<jets("AntiKt4LCTopo"); ++ij){
//    cout << "reco nom jet " << ij << " pt " << jet(ij, "AntiKt4LCTopo").p.Pt() << endl;
  }
  for(int ij=0; ij<jets("AntiKt4Truthtest"); ++ij){
//    cout << "truth jet " << ij << " pt " << jet(ij, "AntiKt4Truthtest").p.Pt() << endl;
  }
  for(int ij=0; ij<jets("AntiKt4Truth"); ++ij){
//    cout << "truth nom jet " << ij << " pt " << jet(ij, "AntiKt4Truth").p.Pt() << endl;
  }

  // end----------------------------------------------------------------------
  if (Debug()) cout << "Analysis_EventTreeFiller_ClustersAndTruth: DEBUG End ProcessEvent(): RunNumber = " << RunNumber() 
		            << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;
  return true;
}



///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_EventTreeFiller_ClustersAndTruth::WorkerTerminate()
{
    fEventTree->Write();


  // Nothing more

}

///===========================================
/// Fill Tree
///========================================
void Analysis_EventTreeFiller_ClustersAndTruth::FillTree(){
  if(Debug()) cout <<"Analysis_EventTreeFiller_ClustersAndTruth::FillTree Start" << endl;

  ResetBranches(fEventTree);
  FillEventVars(fEventTree);
  fEventTree->Fill();
    
  if(Debug()) cout <<"Analysis_EventTreeFiller_ClustersAndTruth::FillTree End" << endl;
  return;
}

///=============================================
/// Add Branches To Tree
///=============================================
void Analysis_EventTreeFiller_ClustersAndTruth::AddBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_EventTreeFiller_ClustersAndTruth::AddBranches Start" << endl;
  
    // Event Info
    tree->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
    tree->Branch("RunNumber",                 &fTRunNumber,              "RunNumber/I");
    tree->Branch("Weight" ,                   &fTWeight,                 "Weight/F");
    tree->Branch("Mu" ,                       &fTMu,                     "Mu/F");
    tree->Branch("NPVtruth" ,                 &fTNPVtruth,               "NPVtruth/I");
    tree->Branch("NPV" ,                      &fTNPVtruth,               "NPV/I");
  
    // clusters --------------------------------------------------------------------------------------
    tree->Branch("cl_lc_n",                   &fTcl_lc_n,                "cl_lc_n/I");
    tree->Branch("cl_lc_px",                  &fTcl_lc_px,               "cl_lc_px[cl_lc_n]/F");
    tree->Branch("cl_lc_py",                  &fTcl_lc_py,               "cl_lc_py[cl_lc_n]/F");
    tree->Branch("cl_lc_pz",                  &fTcl_lc_pz,               "cl_lc_pz[cl_lc_n]/F");
    tree->Branch("cl_lc_E",                   &fTcl_lc_E,                "cl_lc_E[cl_lc_n]/F");
    tree->Branch("cl_em_n",                   &fTcl_em_n,                "cl_em_n/I");
    tree->Branch("cl_em_px",                  &fTcl_em_px,               "cl_em_px[cl_em_n]/F");
    tree->Branch("cl_em_py",                  &fTcl_em_py,               "cl_em_py[cl_em_n]/F");
    tree->Branch("cl_em_pz",                  &fTcl_em_pz,               "cl_em_pz[cl_em_n]/F");
    tree->Branch("cl_em_E",                   &fTcl_em_E,                "cl_em_E[cl_em_n]/F");

    // truths -----------------------------------------------------------------------
    tree->Branch("truth_n",                   &fTtruth_n,                "truth_n/I");
    tree->Branch("truth_px",                  &fTtruth_px,               "truth_px [truth_n]/F");
    tree->Branch("truth_py",                  &fTtruth_py,               "truth_py[truth_n]/F");
    tree->Branch("truth_pz",                  &fTtruth_pz,               "truth_pz[truth_n]/F");
    tree->Branch("truth_E",                   &fTtruth_E,                "truth_E  [truth_n]/F");
    tree->Branch("truth_id",                  &fTtruth_id,               "truth_id  [truth_n]/F");


  if(Debug()) cout <<"Analysis_EventTreeFiller_ClustersAndTruth::AddBranches End" << endl;
    return;
}

void Analysis_EventTreeFiller_ClustersAndTruth::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_EventTreeFiller_ClustersAndTruth::ResetBranches Start" << endl;
  
    // Event Info
    fTEventNumber           = -999;
    fTRunNumber             = -999;
    fTWeight                = -999;
    fTNPVtruth              = -999;
    fTNPV                   = -999;
    fTMu                    = -999;
    
    fTcl_lc_n= 0;
    fTcl_em_n= 0;
    fTtruth_n= 0;
    for(int i=0;i<MaxNClus; ++i){
        fTcl_lc_px   [i] = -999.99;
        fTcl_lc_py  [i] = -999.99;
        fTcl_lc_pz  [i] = -999.99;
        fTcl_lc_E   [i] = -999.99;
        fTcl_em_px   [i] = -999.99;
        fTcl_em_py  [i] = -999.99;
        fTcl_em_pz  [i] = -999.99;
        fTcl_em_E   [i] = -999.99;
    }

    for( int it=0; it<MaxNTruths; ++it){
        fTtruth_px [it] = -999.99;
        fTtruth_py [it] = -999.99;
        fTtruth_pz [it] = -999.99;
        fTtruth_E  [it] = -999.99;
        fTtruth_id [it] = -999.99;
    }

  if(Debug()) cout <<"Analysis_EventTreeFiller_ClustersAndTruth::ResetBranches End" << endl;
    return;
}


void Analysis_EventTreeFiller_ClustersAndTruth::FillEventVars(TTree *tree){
  if(Debug()) cout <<"Analysis_EventTreeFiller_ClustersAndTruth::FillEventVars Begin" << endl;
  
    // Event Info -----------------------------------
    fTEventNumber                 = Int("EventNumber");
    fTRunNumber                   = Int("RunNumber");
    fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
    fTNPV                         = Exists("NPV")?      Int("NPV"):-1;
    fTWeight                      = DefaultWeight();
    fTMu                          = Float("averageIntPerXing");                  

    // fill clusters 
    for(int iJ=0; iJ<clusters("LCTopo"); ++iJ){
        if(cluster(iJ, "LCTopo").p.E() <0) continue;
        if(iJ>=MaxNClus)    continue;
        fTcl_lc_px  [fTcl_lc_n]       = cluster(iJ, "LCTopo").p.Px();
        fTcl_lc_py  [fTcl_lc_n]       = cluster(iJ, "LCTopo").p.Py();
        fTcl_lc_pz  [fTcl_lc_n]       = cluster(iJ, "LCTopo").p.Pz();
        fTcl_lc_E   [fTcl_lc_n]       = cluster(iJ, "LCTopo").p.E();
        fTcl_lc_n++;
    }
    for(int iJ=0; iJ<clusters("EM"); ++iJ){
        if(cluster(iJ, "EM").p.E() <0) continue;
        if(iJ>=MaxNClus)    continue;
        fTcl_em_px  [fTcl_em_n]       = cluster(iJ, "EM").p.Px();
        fTcl_em_py  [fTcl_em_n]       = cluster(iJ, "EM").p.Py();
        fTcl_em_pz  [fTcl_em_n]       = cluster(iJ, "EM").p.Pz();
        fTcl_em_E   [fTcl_em_n]       = cluster(iJ, "EM").p.E();
        fTcl_em_n++;
    }

    //fill truth 
    for(int iT=0; iT<truths("Stable"); ++iT){
        if(iT>=MaxNTruths) continue;
        fTtruth_px  [fTtruth_n] = truth(iT, "Stable").p.Px();
        fTtruth_py  [fTtruth_n] = truth(iT, "Stable").p.Py();
        fTtruth_pz  [fTtruth_n] = truth(iT, "Stable").p.Pz();
        fTtruth_E   [fTtruth_n] = truth(iT, "Stable").p.E();
        fTtruth_id  [fTtruth_n] = truth(iT, "Stable").Int("pdgId");
        fTtruth_n++;
    }


  if(Debug()) cout <<"Analysis_EventTreeFiller_ClustersAndTruth::FillEventVars End" << endl;
  return;
}

