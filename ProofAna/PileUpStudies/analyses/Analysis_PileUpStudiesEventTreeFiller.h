/**************************************************************************
 **
 **   File:         Analysis_PileUpStudiesEventTreeFiller.h
 **
 **   Description:  Prototype PileUp Analysis
 **                 
 ** 
 **   Authors:      P. Nef
 **
 **   Created:      4/12/2013
 **
 **************************************************************************/

#ifndef Analysis_PileUpStudiesEventTreeFiller_h
#define Analysis_PileUpStudiesEventTreeFiller_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

 
using std::cout;
using std::endl;


class Analysis_PileUpStudiesEventTreeFiller : public Analysis_JetMET_Base {

 public :
  
  Analysis_PileUpStudiesEventTreeFiller(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_PileUpStudiesEventTreeFiller() { }
  
  ClassDef(Analysis_PileUpStudiesEventTreeFiller, 0);
  
  Bool_t  fDetail;
  Bool_t  fAddTacksToTree;
  TString fTrkSel;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();


  TTree *fEventTree;
  void AddBranches(TTree *tree);
  void FillTree(const MomKey JetKey, const MomKey JVFKey, const MomKey TrkKey);
  void FillEventVars(TTree *tree, const MomKey JetKey, const MomKey JVFKey, const MomKey TrkKey);
  void ResetBranches(TTree *tree);
  
  // per Event variables
  int   fTEventNumber;
  int   fTRunNumber;
  float fTWeight;
  float fTMu;
  int   fTNPVtruth;
  int   fTNPV;

  // jets ----------------
  static const int MaxNJets = 20;
  int   fTNJets;
  int   fTNJetsFilled;
  float fTJPt              [MaxNJets];
  float fTJEta             [MaxNJets];
  float fTJPhi             [MaxNJets];
  float fTJM               [MaxNJets];
  int   fTJisPU            [MaxNJets];
  int   fTJisHS            [MaxNJets];
  int   fTJisBtagged       [MaxNJets];
  float fTJJVF             [MaxNJets];
  float fTJRpT             [MaxNJets];
  float fTJnPUTrkCorrJVF   [MaxNJets];
  float fTJtruthpt         [MaxNJets];





};

#endif

