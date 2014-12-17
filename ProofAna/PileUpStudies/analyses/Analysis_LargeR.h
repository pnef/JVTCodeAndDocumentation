/**************************************************************************
 **
 **   File:         Analysis_LargeR.h
 **
 **   Description:  Prototype PileUp Analysis
 **                 
 ** 
 **   Authors:      P. Nef
 **
 **   Created:      4/12/2013
 **
 **************************************************************************/

#ifndef Analysis_LargeR_h
#define Analysis_LargeR_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

 
using std::cout;
using std::endl;


class Analysis_LargeR : public Analysis_JetMET_Base {

 public :
  
  Analysis_LargeR(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_LargeR() { }
  
  ClassDef(Analysis_LargeR, 0);
  
  Bool_t  fDetail;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();

  AnaKey GetPtBin(const MomKey JetKey, const int iJet);
  AnaKey GetEtaBin(const MomKey JetKey, const int iJet);
  
  // helper functions to facilitate similar types of plots
  void MakeJetPlots(const MomKey JetKey);
  void MakeJetCollections(const MomKey JetKey);
  void CalculateJVF(MomKey JetType, MomKey TrackType, bool requireCuts=false);
  float GetNPUTrkCorrJVF(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, float scale);
  void AddTrimmedJetCollections(const MomKey JetKey);
  void CalculateDR(const MomKey JetKey, const MomKey reco);


  private :			  
    vector<double>  ptBounds;
    vector<double>  etaBounds;
    TString fTrkSel;


};

#endif

