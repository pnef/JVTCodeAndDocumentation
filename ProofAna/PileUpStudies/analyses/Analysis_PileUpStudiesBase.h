/**************************************************************************
 **
 **   File:         Analysis_PileUpStudiesBase.h
 **
 **   Description:  Prototype PileUp Analysis
 **                 
 ** 
 **   Authors:      P. Nef
 **
 **   Created:      4/12/2013
 **
 **************************************************************************/

#ifndef Analysis_PileUpStudiesBase_h
#define Analysis_PileUpStudiesBase_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

 
using std::cout;
using std::endl;


class Analysis_PileUpStudiesBase : public Analysis_JetMET_Base {

 public :
  
  Analysis_PileUpStudiesBase(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_PileUpStudiesBase() { }
  
  ClassDef(Analysis_PileUpStudiesBase, 0);
  
  Bool_t  fDetail;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();


  // Event and Trigger
  bool EventSelection();

  //make jet collections
  void MakeJetCollections(const MomKey JetKey);

  // truth info
  void AddTruthHZZDaugthers();

  // event weight
  void SetEventWeight();

  private :			  
    bool  fDoLeptonSelection;
    bool  fDoQCDSelection;
    bool  fDoHZZSelection;
    bool  fDoTopSelection;
    bool  fRequireTrueHSvertex;
    bool  fisMuScan;

    int    SystType;
    MomKey BKey;
    MomKey MuonSFKey;
    MomKey MuonTrigSFKey;

};

#endif

