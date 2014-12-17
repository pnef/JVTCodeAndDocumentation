/**************************************************************************
 **
 **   File:         Analysis_JetMET_Base.h
 **
 **   Description:  Base Analysis class for JetMET analysis, has some nice functions
 **                                  
 **
 **   Authors:      M. Swiatlowski
 **
 **   Created:      2013-01-22
 **
 **************************************************************************/

#ifndef Analysis_JetMET_Base_h
#define Analysis_JetMET_Base_h

class LeptonTriggerSF;
namespace Analysis {
  class AnalysisMuonConfigurableScaleFactors;
  class CalibrationDataInterfaceROOT;
}
class DataPeriod;


#include "AnalysisBase.h"


//#include "TLVC.h"

#include <map>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
//#include "PlotHelper.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"


using std::cout;
using std::endl;


class PJ_Info: public fastjet::PseudoJet::UserInfoBase {
public:
    PJ_Info(bool GhostTrackIn, bool GhostParticleIn, Particle* PAPointer){
       GhostTrack = GhostTrackIn;
       GhostParticle = GhostParticleIn; 
       Pointer = PAPointer;
       Cluster = ! (GhostTrack || GhostParticle);
    }
    bool GhostTrack, GhostParticle, Cluster;
    Particle* Pointer;
};



class Analysis_JetMET_Base : public AnalysisBase {

public :

  Analysis_JetMET_Base(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  virtual ~Analysis_JetMET_Base() { }    
  
  ClassDef(Analysis_JetMET_Base, 0);


  Bool_t  fDetail;
  
  virtual bool    ProcessEvent(){/* To be overridden by derived analysis class */ 
    std::cout<<"%%%%% Analysis_JetMET_Base::ProcessEvent() %%%%%%%%%%%%%"<<std::endl;
    return kTRUE;};
  virtual void    WorkerBegin();
  virtual void    WorkerTerminate(){/* To be overridden by derived analysis class */;};  


  //making jets
  vector<fastjet::PseudoJet>    ObjsToPJ(const MomKey type, MomentObj* jet = 0);
  MomKey                        MakeJets(const fastjet::JetAlgorithm algo, const double jetR, const MomKey constType, const MomKey extra = "", const double minpt=10);
  MomKey                        MakePrunedJets(const MomKey JetType, const fastjet::JetAlgorithm algo, const double zcut, const double dcut);
  MomKey                        MakePrunedJetsName(const MomKey JetType, const fastjet::JetAlgorithm algo, const double zcut, const double dcut);
  void                          RemoveJetsPtCut(const MomKey JetType, float PtCut);

  // ghost assocation
  void                          AddGhostMatch(const MomKey JetKey, const MomKey MatchObject);
  void                          AddGhostMatch(const MomKey JetKey, const MomKey MatchObject, const MomKey ConstType, const fastjet::JetAlgorithm AlgoIn, const float Radius);
  void                          AddAsGhosts(const MomKey GhostObject, vector<fastjet::PseudoJet>* TheThings, bool IsTrack);

  // delta r assoc
  bool                          AddDRMatch(const MomKey JetKey, const MomKey MatchObject);

  // related: truth stable particles for truth jets
  void                          AddTruthStableCharged(const MomKey JetKey);



  //setting up truth stable particles for truth jets
  MomKey                        AddStableParticles();


  // setting up tracks for track jets
  MomKey                        AddGoodTracks();


  //making substructure
  void                          AddNsub(const MomKey constType, const MomKey jetCollection, const bool min = true, const float beta = 1);
  void                          CalcNsub(Particle* myjet, const MomKey constType, const bool min, const float beta, const float jetRad);
  MomKey                        AddSubjets(const MomKey JetType, const fastjet::JetAlgorithm algo, const double fcut = 0.04, const double ptcut = 40.);
  MomKey                        SetEventSubjettiness(const MomKey JetType, const MomKey nsubType, const int NJets = 4);
  MomKey                        SetEventMass(const MomKey JetType, int NJets = 4);
  MomKey                        SetHt(const MomKey JetType, int NJets = -1);
  MomKey                        SetEventSubjets(const MomKey JetType, const MomKey SubjetType, int NJets = 4);



  //utility functions
  void                          SetFatJetOverlap(const MomKey FatKey, const MomKey SmallKey, float dR);
  void                          MakeDRMatch(const MomKey JetKey1, const MomKey JetKey2);
  void                          CalculateResponseOffset(const MomKey JetKey, const MomKey MatchKey);
  void                          MakeDRResponseOffset(const MomKey JetKey, const MomKey TruthJetKey);
  void                          PrintParticle(Particle* part);
  bool                          SumMatchedCollection(Particle* theJet, const MomKey collection);
  bool                          SumMatchedCollectionJets(const MomKey JetKey, const MomKey collection);
  bool                          ManuallyAddConstituents(const MomKey JetKey);

  // muon scale factor calculation, for both trigger and nominal
  void                          DoMuonSF();
  MomKey                        GetMuonSFKey(const int SystType);
  MomKey                        GetMuonTrigSFKey(const int SystType);

  // b-tag scale factor calculation
  void                          DoBtaggingSF(bool sys, const MomKey JetType);
  MomKey                        GetBTagKey(const int SystType);

  // b-labelling functions
  bool                          MakeTruthBHadronCollection();
  bool                          IsBHadron(int mc_pdgId);
  bool                          IsWeaklyDecayingBHadron(int mc_pdgId);
  bool                          IsBHadronFromTop(int iTruth);
  bool                          GetBHadronDecayProducts(Particle* part, bool printFinalProducts, bool printTree);
  bool                          LabelBJetsDR(const MomKey JetKey);
  bool                          ReduceTruthBHadronCollection(const MomKey TruthKey = MomKey("BHadron"));


  ////////////////////////////
  ///2D histos
  ////////////////////////////
  TH2* Fill(AnaKey name, double XVal, double YVal, double weight, 
	      int NumBinsX, double XMin, double XMax, int NumBinsY, double YMin, double YMax);
  TH2* Fill(AnaKey name, double XVal, double YVal, double weight, 
	    int NumBinsX, const double* Xbins, int NumBinsY, double YMin, double YMax);
  TH2* Fill(AnaKey name, double XVal, double YVal, double weight, 
	    int NumBinsX, double XMin, double XMax, int NumBinsY, const double* Ybins);
  TH2* Fill(AnaKey name, double XVal, double YVal, double weight, 
	      int NumBinsX, const float* Xbins, int NumBinsY, const float* Ybins);
  TH2* Fill(AnaKey name, std::string xlabel, std::string ylabel, double weight,
	      std::vector<std::string> x_labels, std::vector<std::string> y_labels);

  ////////////////////////////
  ///3D histos
  ////////////////////////////
  TH3* Fill(AnaKey name, double XVal, double YVal, double ZVal, double weight, 
	      int NumBinsX, double XMin, double XMax, int NumBinsY, double YMin, double YMax, int NumBinsZ, double ZMin, double ZMax);

  
  ////////////////////////////
  ///1D histos
  ////////////////////////////
  TH1* Fill(AnaKey name, double XVal, double weight, int NumBinsX, double XMin, double XMax);   
  TH1* Fill(AnaKey name, double XVal, double weight, int NumBinsX, const float* Xbins);
  TH1* Fill(AnaKey name, std::string label, double weight, std::vector<std::string> x_labels);

  // Tools for muon scale factors
  LeptonTriggerSF*                                    MuonTriggerSF;
  Analysis::AnalysisMuonConfigurableScaleFactors*     MuonSF;
  DataPeriod*                                         myDataPeriod;

  // tools for b-tag scale factors

   Analysis::CalibrationDataInterfaceROOT* myBTagCalib;
   Analysis::CalibrationDataInterfaceROOT* myBTagCalib70;
   Analysis::CalibrationDataInterfaceROOT* myBTagCalib80;
   Analysis::CalibrationDataInterfaceROOT* myBTagCalib90;
   Analysis::CalibrationDataInterfaceROOT* myBTagCalibLead;
   Analysis::CalibrationDataInterfaceROOT* myBTagCalibLead3;
   Analysis::CalibrationDataInterfaceROOT* myJVFBTagCalib;
   Analysis::CalibrationDataInterfaceROOT* myJVFBTagCalib70;
   Analysis::CalibrationDataInterfaceROOT* myJVFBTagCalib80;
   Analysis::CalibrationDataInterfaceROOT* myJVFBTagCalib90;
   Analysis::CalibrationDataInterfaceROOT* myJVFBTagCalibLead;
   Analysis::CalibrationDataInterfaceROOT* myJVFBTagCalibLead3;


   vector<MomKey> keys;
   vector<float> opvals;
   vector<std::string> OPs;
   vector<Analysis::CalibrationDataInterfaceROOT*> BTagCalibs;
   vector<Analysis::CalibrationDataInterfaceROOT*> JVFBTagCalibs;
};

#endif
