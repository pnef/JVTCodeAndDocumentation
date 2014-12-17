/**************************************************************************
 **
 **   File:         Analysis_LargeRTreeFiller.h
 **
 **   Description:  Prototype PileUp Analysis
 **                 
 ** 
 **   Authors:      P. Nef
 **
 **   Created:      4/12/2013
 **
 **************************************************************************/

#ifndef Analysis_LargeRTreeFiller_h
#define Analysis_LargeRTreeFiller_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

 
using std::cout;
using std::endl;


class Analysis_LargeRTreeFiller : public Analysis_JetMET_Base {

 public :
  
  Analysis_LargeRTreeFiller(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_LargeRTreeFiller() { }
  
  ClassDef(Analysis_LargeRTreeFiller, 0);
  
  Bool_t  fDetail;
  Bool_t  fAddTacksToTree;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();

  TString fTrkSel;


  TTree *fJetTree;
  TTree *fEventTree;
  void AddBranches(TTree *tree);
  void FillTree(const MomKey JetKey);
  void FillEventVars(TTree *tree, const MomKey JetKey);
  void FillJetVars(TTree *tree, int jindex, const MomKey JetKey);
  void ResetBranches(TTree *tree);
  bool fisMuScan;
  
  // per Event variables (but filled on jet-by-jet basis)
  int   fTEventNumber;
  int   fTRunNumber;
  float fTWeight;
  float fTMu;
  float fTNPVtruth;

  // jet by jet
  float fTJpt;
  float fTJm;
  float fTJeta;
  float fTJphi;
  float fTJWdR;
  float fTJWc1dR;
  float fTJWc2dR;
  float fTJZdR;
  float fTJZc1dR;
  float fTJZc2dR;
  float fTJTopdR[2];
  float fTJTopc1dR[2];
  float fTJTopc2dR[2];
  int   fTJindex;

  int   fTJNSub;
  float fTJsubPt [300];
  float fTJsubEta[300];
  float fTJsubPhi[300];
  float fTJsubM  [300];
  float fTJsubJVF [300];
  float fTJsubnPUtrkCorrJVF[300];

  float fTJTrimmedPt[300];
  float fTJTrimmedM [300];

  // 
  vector<TString> trimmedJetNames;

};

#endif

