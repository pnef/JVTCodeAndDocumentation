/**************************************************************************
 **
 **   File:         Analysis_PileUpStudies.h
 **
 **   Description:  Prototype PileUp Analysis
 **                 
 ** 
 **   Authors:      P. Nef
 **
 **   Created:      4/12/2013
 **
 **************************************************************************/

#ifndef Analysis_PileUpStudies_h
#define Analysis_PileUpStudies_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"
#include "JetVertexTagger/JetVertexTagger.h"

 
using std::cout;
using std::endl;


class Analysis_PileUpStudies : public Analysis_JetMET_Base {

 public :
  
  Analysis_PileUpStudies(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_PileUpStudies() { }
  
  ClassDef(Analysis_PileUpStudies, 0);
  
  Bool_t  fDetail;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();

  AnaKey GetPtBin(const MomKey JetKey, const int iJet);
  AnaKey GetEtaBin(const MomKey JetKey, const int iJet);
  
  bool   MakeJetMassCut(MomKey JetKey);
  bool   PassJPtCut(MomKey JetKey, int index, float ptcut, bool isGood );
  void PrintEvent();

  // Event and Trigger
  bool EventSelection();
  bool TriggerSelection();

  // helper functions to facilitate similar types of plots
  void MakeJetPlots(const MomKey JetKey);
  void MakeJetPlots(const MomKey JetKey, const MomKey JVFKey);
  void MakeTruthJetPlots(const MomKey JetKey);
  void MakeVertexPlots();
  void JetByJetHistFiller(Particle* jet, const MomKey JVFKey, const TString prefix);

  // JVFHelpers
  void          CalculateJVF(MomKey JetType, MomKey TrackType, bool requireCuts=false);
  float         GetValueOrderedJVF(Particle *jet, const MomKey JVFKey, unsigned int position);
  vector<float> GetValueOrderedJVF(Particle *jet, const MomKey JVFKey);
  vector<float> GetValueOrderedJVPt(Particle *jet, const MomKey JVFKey);
  float         GetJVPt(Particle *jet, const MomKey JVFKey, const unsigned int position);
  vector<float> GetJVF(Particle *jet, const MomKey JVFKey);
  float         GetMeanJVF(Particle *jet, const MomKey JVFKey);
  float         GetRMSJVF(Particle *jet, const MomKey JVFKey);
  float         GetNxJVF(Particle* jet, const MomKey JVFKey, float Xvalue);
  float         GetJVFEntropy(Particle* jet, const MomKey JVFKey);
  float         GetJetPVWidth(Particle* jet, const MomKey JVFKey, const int PV);
  float         GetJetCaloWidth(Particle* jet, bool squared);
  float         GetJetPVTrkSum(Particle* jet, const MomKey JVFKey, const int PV, const bool dRWeight, const bool useTrackAxis=false);
  float         GetJetPVTrkMeanDR(Particle* jet, const MomKey JVFKey, const int PV);
  float         GetJetAnnulusPtRatio(Particle *jet, float dRmin, float dRmax);
  float         GetAverageDeltaRSquared(Particle *jet);
  float         GetJetLeadingTrkDeltaZ(Particle *jet, const MomKey TrackKey, bool log, bool sqrt);
  float         GetJetAverageTrkDeltaZ(Particle *jet, const MomKey TrackKey, bool log, bool sqrt);
  float         GetJetBetaStar(Particle *jet, const MomKey TrackKey);
  float         GetMuCorrJVF(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, float scale);
  float         GetNPVCorrJVF(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, float scale);
  float         GetNPUTrkCorrJVF(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, float scale);
  float         GetNContribVtxCorrJVF(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, float scale);
  int           GetNContribVtxs(Particle *jet, const MomKey JVFKey);
  float         GetWeightedPUTrkSumPt(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, bool dRWeight, bool pTWeight);
  float         GetWeightedHSTrkSumPt(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, bool dRWeight, bool pTWeight);
  float         GetClosebyDistance(Particle *jet, const MomKey JetKey, bool kt_metric, int& index);
  float         GetConstitPtWeightedMeanOverDr(Particle* myjet);
  float         GetConstitPtWeightedStdDevOverDr(Particle* myjet);
  float         GetConstitPtWeightedSkewnessOverDr(Particle* myjet); 
  float         GetConstitPtWeightedKurtosisOverDr(Particle* myjet); 


  // Tracks
  vector<Particle*> GetJetPVTracks(Particle* jet, const MomKey JVFKey, const int PV);
  void        TrackSelector(MomKey TrkKey);

  // print interesting event
  void JetAnalysisDisplay();
  void JVFPrintHist(const MomKey JetKey, bool pileup);


  //make jet collections
  void MakeJetCollections(const MomKey JetKey);
  void CalculateJetVars(const MomKey JetKey, const MomKey JVFKey, const MomKey TrackKey);
  void CalculateJetVarsCalo(const MomKey JetKey);


  private :			  
    vector<double>  ptBounds;
    vector<double>  etaBounds;


    TString  fTrkSel;
    bool    fisMuScan;

    int pucounter;
    int hscounter;

    JetVertexTagger* jvt;

};

#endif

