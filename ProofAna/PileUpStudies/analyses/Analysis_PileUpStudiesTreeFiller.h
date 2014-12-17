/**************************************************************************
 **
 **   File:         Analysis_PileUpStudiesTreeFiller.h
 **
 **   Description:  Prototype PileUp Analysis
 **                 
 ** 
 **   Authors:      P. Nef
 **
 **   Created:      4/12/2013
 **
 **************************************************************************/

#ifndef Analysis_PileUpStudiesTreeFiller_h
#define Analysis_PileUpStudiesTreeFiller_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

 
using std::cout;
using std::endl;


class Analysis_PileUpStudiesTreeFiller : public Analysis_JetMET_Base {

 public :
  
  Analysis_PileUpStudiesTreeFiller(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_PileUpStudiesTreeFiller() { }
  
  ClassDef(Analysis_PileUpStudiesTreeFiller, 0);
  
  Bool_t  fDetail;
  Bool_t  fAddTacksToTree;
  Bool_t  fAddLeptonZInfoToTree;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();



  TTree *fJetTree;
  TTree *fEventTree;
  TString fTreeName;
  TString fTrkSel;

  bool fDoTrkFromVxt;
  void AddBranches(TTree *tree);
  void FillTree(const MomKey JetKey, const MomKey JVFKey,  const MomKey trackKey);
  void FillEventVars(TTree *tree, const MomKey JetKey, const MomKey trackKey);
  void FillJetVars(TTree *tree, int jindex, const MomKey JetKey, const MomKey JVFKey, const MomKey trackKey);
  void ResetBranches(TTree *tree);
  bool fisMuScan;
  
  // per Event variables (but filled on jet-by-jet basis)
  int   fTEventNumber;
  int   fTRunNumber;
  float fTWeight;
  float fTPUWeight;
  float fTLumiXSecWeight;
  float fTm1MuScaleFactor;
  float fTm1MuTrigScaleFactor;
  float fTNPV;
  float fTNPVtruth;
  float fTMu;
  int   fTNBJets;
  int   fTGRL;
  float fTRhoKt4LC;
  int   fTnPUTracks;
  int   fTnTracksNotFromPV;

  int   fTNJetspt20eta2p1;
  int   fTNJetspt30eta2p1;
  int   fTNJetspt40eta2p1;
  int   fTNJetspt50eta2p1;
  int   fTNJetspt20eta2p1HS;
  int   fTNJetspt30eta2p1HS;
  int   fTNJetspt40eta2p1HS;
  int   fTNJetspt50eta2p1HS;

  float fTZpt;
  float fTZmass;
  float fTZeta;
  float fTZphi;

  static const int MaxNMuons =2;
  int   fTNMuons;
  float fTMuonPt[MaxNMuons];
  float fTMuonEta[MaxNMuons];
  float fTMuonPhi[MaxNMuons];
  float fTMuonPtUnsmeared[MaxNMuons];
  float fTMuonEtaUnsmeared[MaxNMuons];
  float fTMuonPhiUnsmeared[MaxNMuons];


  float fTJ0Pt;
  float fTJ1Pt;
  float fTJ2Pt;
  float fTJ0Eta;
  float fTJ1Eta;
  float fTJ2Eta;
  float fTJ0Phi;
  float fTJ1Phi;
  float fTJ2Phi;
  float fTJ0ZDPhi;
  float fTJ1ZDPhi;
  float fTJ2ZDPhi;
  float fTJ0MuonLeadDR;
  float fTJ1MuonLeadDR;
  float fTJ2MuonLeadDR;
  float fTJ0MuonTrailDR;
  float fTJ1MuonTrailDR;
  float fTJ2MuonTrailDR;
  float fTJLeadingTruthPt;

  // jet by jet
  float fTJpt;
  float fTJeta;
  float fTJphi;
  float fTJconstitpt;
  float fTJareacorrpt;
  float fTJtruthpt;
  float fTJtrutheta;
  float fTJtruthphi;
  float fTJtruthMatchDR;
  float fTJtruthInTimept;
  float fTJtruthInTimeMatchDR;
  float fTJdRTruthPt4;
  float fTJdRTruthPt6;
  float fTJdRTruthPt8;
  float fTJdRTruthPt10;
  float fTClosebydR;
  float fTClosebykt;
  float fTClosebydRpT;
  int   fTJindex;
  int   fTJisPU;
  int   fTJisHS;
  int   fTJisParton_uds;
  int   fTJisParton_c;
  int   fTJisParton_b;
  int   fTJisParton_g;
  int   fTJisBtagged;
  float fTJArea;
  float fTJJVF;
  float fTJstdJVF;
  float fTJcorrJVF;
  float fTJJVFZeroOverMax;
  float fTJMaxJVF;
  float fTJHSPVtrkPtlead;
  float fTJHSPVtrkPttrail;
  float fTJLeadTrkPt;
  float fTJHSPVTrk12dR;
  float fTJHSPVtrkSumOverPt;
  float fTJHSPVtrkSumOverPtdRW;
  float fTJHSPVTrkWidth;
  float fTJHSPVTrkMeanDR;
  float fTJMVA1result;
  float fTJMVA2result;
  float fTJMVA3result;
  float fTJMVA4result;
  float fTJTiming;
  float fTJEMFrac;       
  float fTJPUtrkPtSumOverPt;
  float fTJCorrPUtrkPtSumOverPt;
  float fTJCorr2PUtrkPtSumOverPt;
  float fTJPUtrkPtSum;
  float fTJPUtrkPtSumSquared;
  float fTJPUtrkPtSumdRWeight;
  float fTRhoGoodPUTracks;
  float fTRhoPUTracks;
  float fTJTruthMF;       
  int   fTJisBadLoose;    
  float fTJMuonMinDR;
  float fTJZDR;
  float fTJZDPhi;
  float fTJMuonLeadDR;
  float fTJMuonTrailDR;
  float fTJCaloWidth;
  int   fTJNClusters;
  float fTJNumTowers;
  float fTJCaloWidth2;
  float fTJAnnulusPtRatio0to1;
  float fTJAnnulusPtRatio0to2;
  float fTJAnnulusPtRatio1to2;
  float fTJAnnulusPtRatio2to4;
  float fTJAverageDeltaRSquared;
  float fTJLeadTrkAllHSPVdZ;
  float fTJLeadTrkAllHSPVdZTrunc;
  float fTJLeadTrkAllHSPVdZlog;
  float fTJLeadTrkAllHSPVdZsqrt;
  float fTJAverageTrkAllHSPVdZ;
  float fTJAverageTrkAllHSPVdZlog;
  float fTJAverageTrkAllHSPVdZsqrt;
  float fTJBetaStar;
  float fTJHStrkPtSum;
  float fTJHStrkPtSumdRWeight;
  float fTJHStrkPtSumSquared;
  int   fTJNContribVtxs;
  float fTJNPVCorrJVF;
  float fTJNContribVtxsCorrJVF;
  float fTJMuCorrJVF;
  float fTJMuCorrJVF0p1;
  float fTJMuCorrJVF0p14;
  float fTJnPUTrkCorrJVF;
  float fTJJVT_RootCore     ;
  float fTJcorrJVF_RootCore ;
  float fTJRpT_RootCore     ;








  static const int MaxNTrks = 20;
  int   fTNTracks;
  int   fTNTracksHS;
  int   fTNTracksPU;
  int   fTNTracksFilled;
  float fTTrkPt [MaxNTrks];
  float fTTrkEta[MaxNTrks];
  float fTTrkPhi[MaxNTrks];
  float fTTrkd0[MaxNTrks];
  float fTTrkZ0[MaxNTrks];
  float fTTrkZ0errWrtPV[MaxNTrks];
  float fTTrkd0errWrtPV[MaxNTrks];
  int   fTTrkNPixHits[MaxNTrks];
  int   fTTrkNSCTHits[MaxNTrks];
  int   fTTrkNTRTHits[MaxNTrks];
  float fTTrkChiSq[MaxNTrks];
  float fTTrkNDof[MaxNTrks];
  int   fTTrkMCBarCode[MaxNTrks];
  int   fTTrkOrigin[MaxNTrks];

  float fTVtxDzTruth;
  float fTHSVtxTruthminDz;
  float fTHSPV_mindZPUVtx;
  int   fTNMCVertWithin1mmFromHS;

  // TMVA 
  TMVA::Reader *fTMVAreader;
  bool  doMVAEval;
  bool  runCluster;
  TString fMVA1name;
  TString fMVA2name;
  TString fMVA3name;
  TString fMVA4name;
  TString fMVA1file;
  TString fMVA2file;
  TString fMVA3file;
  TString fMVA4file;


};

#endif

