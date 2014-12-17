/**************************************************************************
 **
 **   File:         Analysis_EventTreeFiller_ClustersAndTruth.h
 **
 **   Description:  Template for event-by-event tree filler
 **                 
 ** 
 **   Authors:      P. Nef
 **
 **   Created:      2/7/2014
 **
 **************************************************************************/

#ifndef Analysis_EventTreeFiller_ClustersAndTruth_h
#define Analysis_EventTreeFiller_ClustersAndTruth_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

 
using std::cout;
using std::endl;


class Analysis_EventTreeFiller_ClustersAndTruth : public Analysis_JetMET_Base {

 public :
  
  Analysis_EventTreeFiller_ClustersAndTruth(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_EventTreeFiller_ClustersAndTruth() { }
  
  ClassDef(Analysis_EventTreeFiller_ClustersAndTruth, 0);
  
  Bool_t  fDetail;
  Bool_t  fAddTacksToTree;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();


  TTree *fEventTree;
  void AddBranches(TTree *tree);
  void FillTree();
  void FillEventVars(TTree *tree);
  void ResetBranches(TTree *tree);
  
  // per Event variables
  int   fTEventNumber;
  int   fTRunNumber;
  float fTWeight;
  float fTMu;
  int   fTNPVtruth;
  int   fTNPV;

  // clusters and truths 
  static const int MaxNClus = 5000;
  int   fTcl_lc_n;
  float fTcl_lc_px  [MaxNClus];
  float fTcl_lc_py  [MaxNClus];
  float fTcl_lc_pz  [MaxNClus];
  float fTcl_lc_E   [MaxNClus];

  int   fTcl_em_n;
  float fTcl_em_px  [MaxNClus];
  float fTcl_em_py  [MaxNClus];
  float fTcl_em_pz  [MaxNClus];
  float fTcl_em_E   [MaxNClus];

  static const int MaxNTruths = 1000;
  int   fTtruth_n;
  float fTtruth_px   [MaxNTruths];
  float fTtruth_py   [MaxNTruths];
  float fTtruth_pz   [MaxNTruths];
  float fTtruth_E    [MaxNTruths];
  float fTtruth_id   [MaxNTruths];


};

#endif

