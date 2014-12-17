/**************************************************************************
 **
 **   File:         Analysis_PileUpStudiesTruthTreeFiller.h
 **
 **   Description:  Prototype PileUp Analysis
 **                 
 ** 
 **   Authors:      P. Nef
 **
 **   Created:      4/12/2013
 **
 **************************************************************************/

#ifndef Analysis_PileUpStudiesTruthTreeFiller_h
#define Analysis_PileUpStudiesTruthTreeFiller_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

 
using std::cout;
using std::endl;


class Analysis_PileUpStudiesTruthTreeFiller : public Analysis_JetMET_Base {

 public :
  
  Analysis_PileUpStudiesTruthTreeFiller(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_PileUpStudiesTruthTreeFiller() { }
  
  ClassDef(Analysis_PileUpStudiesTruthTreeFiller, 0);
  
  Bool_t  fDetail;
  Bool_t  fAddTacksToTree;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();



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
  float fTJeta;
  float fTJphi;
  int   fTJindex;
  float fTJmatchPt;
  float fTJmatchConstitPt;
  float fTJmatchAreacorrPt;
  float fTJmatchJVF;
  float fTJmatchcorrJVF;
  float fTJmatchRpT;

};

#endif

