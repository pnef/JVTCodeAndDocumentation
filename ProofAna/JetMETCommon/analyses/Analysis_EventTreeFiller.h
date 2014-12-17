/**************************************************************************
 **
 **   File:         Analysis_EventTreeFiller.h
 **
 **   Description:  Template for event-by-event tree filler
 **                 
 ** 
 **   Authors:      P. Nef
 **
 **   Created:      2/7/2014
 **
 **************************************************************************/

#ifndef Analysis_EventTreeFiller_h
#define Analysis_EventTreeFiller_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

 
using std::cout;
using std::endl;


class Analysis_EventTreeFiller : public Analysis_JetMET_Base {

 public :
  
  Analysis_EventTreeFiller(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_EventTreeFiller() { }
  
  ClassDef(Analysis_EventTreeFiller, 0);
  
  Bool_t  fDetail;
  Bool_t  fAddTacksToTree;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();


  TTree *fEventTree;
  void AddBranches(TTree *tree);
  void FillTree(const MomKey JetKey);
  void FillEventVars(TTree *tree, const MomKey JetKey);
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
  float fTJPt [MaxNJets];
  float fTJEta[MaxNJets];
  float fTJPhi[MaxNJets];
  float fTJM  [MaxNJets];

};

#endif

