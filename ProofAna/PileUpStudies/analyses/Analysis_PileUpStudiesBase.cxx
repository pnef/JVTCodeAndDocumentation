/**************************************************************************
 **
 **   File:         Analysis_PileUpStudiesBase.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      P. Nef
 **
 **************************************************************************/

#define Analysis_PileUpStudiesBase_cxx

#include "Analysis_PileUpStudiesBase.h"
#include "AnaConfig.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>      // std::setprecision
#include <cstdlib>
#include <sstream>
#include "TKey.h"
#include "TObjString.h"
#include "TObjArray.h"


///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_PileUpStudiesBase::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_PileUpStudiesBase: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();

  // config -------------
  fDoLeptonSelection = false;
  fDoQCDSelection    = false;
  fDoHZZSelection    = false;
  fDoTopSelection    = false;
  fRequireTrueHSvertex = false;
  fisMuScan            = false;
  SystType             = 0;

  ChainCfg()->Get("requireTrueHSvertex",fRequireTrueHSvertex); 
  ChainCfg()->Get("doLeptonSelection",  fDoLeptonSelection);
  ChainCfg()->Get("doQCDSelection",     fDoQCDSelection);
  ChainCfg()->Get("doHZZSelection",     fDoHZZSelection);
  ChainCfg()->Get("doTopSelection",     fDoTopSelection);
  ChainCfg()->Get("isMuScan",           fisMuScan);
  ChainCfg()->Get("SYSTTYPE"          , SystType     );

  BKey          = GetBTagKey(SystType);
  MuonSFKey     = GetMuonSFKey(SystType);
  MuonTrigSFKey = GetMuonTrigSFKey(SystType);

  if (Debug()) cout << "Analysis_PileUpStudiesBase: DEBUG Finish WorkerBegin()" << endl;  
} 


///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_PileUpStudiesBase::ProcessEvent()
{

  if (Debug()) cout << "Analysis_PileUpStudiesBase: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();
  
  // add truth H->ZZ->llll particles
  if(fDoHZZSelection){AddTruthHZZDaugthers();}
  
  // Jet Collections   --------------------------------------------------------
  static const MomKey JetKey4LC("AntiKt4LCTopo");
  MakeJetCollections(JetKey4LC);

  // SetEventWeight
  SetEventWeight();

  // event and trigger selection---------------------------------------------------------------
  if( EventSelection() )                          Set("PileUpStudiesSelection", true ); 
  else                                            Set("PileUpStudiesSelection", false);
  if(fDoHZZSelection &&  truths("ZDaugthers")==4) Set("PileUpStudiesSelection", true ); 

  if(! Bool("PileUpStudiesSelection"))        return true;

  return true;
}

///---------------------------------------------------
/// Add H->ZZ->llll e and mu daugthers
/// --------------------------------------------------
void Analysis_PileUpStudiesBase::AddTruthHZZDaugthers(){
      AddVec("truthsZDaugthers"); // electrons or muons from H->ZZ->llll
      for(int iT=0; iT<truths(); ++iT){
          if (abs(truth(iT).Int("pdgId"))!=13 && abs(truth(iT).Int("pdgId"))!=11) continue;
          if (truth(iT).Int("status") != 1)                                       continue;
          if (! (truth(iT).Objs("parents")>0))                                    continue;

          Particle *parent = (Particle*) truth(iT).Obj("parents",0);
          if (parent->Int("pdgId") != 23)                                         continue;
          Add("truthsZDaugthers", &truth(iT));
      }
      return;
}

///-----------------------------------------------------
/// SetEventWeight                                      
///---------------------------------------------------
void Analysis_PileUpStudiesBase::SetEventWeight(){
    if      (Exists("PeriodAB_lumi")) Set("PUWeight", Float("PeriodAB_lumi"));
    else if (Exists("Full_lumi")    ) Set("PUWeight", Float("Full_lumi"));
    else                              Set("PUWeight", 0);

    if(!isMC()){ // data
        Set("EventWeight", 1);
    }else if(isMC() && fDoTopSelection){
        if(!Exists(BKey))          Set(BKey,1.);
        if(!Exists(MuonSFKey))     Set(MuonSFKey,1.);
        if(!Exists(MuonTrigSFKey)) Set(MuonTrigSFKey,1.);
        if(!Exists("PUWeight"))    Set("PUWeight",1.); 
        Set("EventWeight", (DefaultWeight() * Float(MuonSFKey) * Float(MuonTrigSFKey) * Float(BKey) * Float("PUWeight")));
    }else if(isMC() && fDoLeptonSelection){
        Set("EventWeight", (DefaultWeight() *  Float("PUWeight")));
    }else if(isMC() && fDoHZZSelection){
        Set("EventWeight", (DefaultWeight() *  Float("PUWeight")));
    }else if(isMC() && fDoQCDSelection){
        Set("EventWeight", (DefaultWeight()));         // no PU reweighting
    }else{
        Set("EventWeight", 0);                // something's wrong
    }
}

/// -----------------------------------------------------
/// Event Selection 
/// -----------------------------------------------------
bool Analysis_PileUpStudiesBase::EventSelection(){
  if(Debug()) cout <<"Analysis_PileUpStudiesBase::EventSelection Start" << endl;

  if(isMC() && fRequireTrueHSvertex) {
      if(Debug()) cout << "Recognized MC...dzTruth is " << vtx(0).Float("dzTruth") << endl;
  
      if(fabs(vtx(0).Float("dzTruth")) > 0.1) {
      cout << "PV[0] is more than 1 mm away from the hard scatter!" << endl;
      return false;
      }
  }

  // Event Selection ---------------------------------------------------------
  if (fDoLeptonSelection){
      if(Debug()) cout << "Analysis_PileUpStudiesBase: doLeptonSelection"     << endl;
      if( ! (Objs("muonsgood")==2))                                                  {return false; } 
      if( ! Bool("TopSelection_GRL"))                                                {return false; }
      if( ! Bool("TopSelection_LarError"))                                           {return false; }
      if( ! Bool("TopSelection_NoBadLooseJets"))                                     {return false; }
      if( ! (muon(0, "good").Float("charge")*muon(1, "good").Float("charge") == -1)) {return false; }
      if( ! Bool("1muontrigger"))                                                    {return false;}
      // DO NOT CUT ON JETS!! THIS WILL BIAS THE JET HS MULTIPLICTY WITH NPV
//      if (! (jets("AntiKt4LCTopo")>0))                        {return false; } 
//      if (! (jet(0,"AntiKt4LCTopo").p.Pt()>20))               {return false; }
//      if (! (fabs(jet(0,"AntiKt4LCTopo").p.Eta())<2.1))       {return false; }

      // reconstruct Z
      AddVec("recosZCandidate");
      Particle* Z = new Particle();
      Z->p = muon(0, "good").p + muon(1, "good").p;
      Add("recosZCandidate", Z);
      if( ! ((reco(0, "ZCandidate").p.M() > 71) && 
             (reco(0, "ZCandidate").p.M()  < 111)))           return false;  
  }
  
  if (fDoHZZSelection){
      if(Debug()) cout << "Analysis_PileUpStudiesBase: fDoHZZSelection"     << endl;
      if(truths("ZDaugthers") != 4)                           return false;
  }

  if (fDoQCDSelection){
      if(Debug()) cout << "Analysis_PileUpStudiesBase: doQCDSelection"        << endl;
      if (! Bool("QCDSelection_All"))                         return false;
      if (! (jets("AntiKt4LCTopo")>0))                        return false;
      if (! (jet(0,"AntiKt4LCTopo").p.Pt()>20))               return false;
      if (! (fabs(jet(0,"AntiKt4LCTopo").p.Eta())<2.4))       return false;
      
  }

  if (fDoTopSelection){
      if(Debug()) cout << "Analysis_PileUpStudiesBase: doTopSelection"        << endl;
      if(! Bool("TopSelection_All"))                          return false;
  }

  if(Debug()) cout <<"Analysis_PileUpStudiesBase::EventSelection End" << endl;
  return true;

}

// --------------------------------------------------
// Make jet collections
// ----------------------------------------------
void Analysis_PileUpStudiesBase::MakeJetCollections(const MomKey JetKey){
  if(Debug()) cout <<"Analysis_PileUpStudiesBase::MakeJetCollections Start" << endl;
    AddVec("jets"+JetKey+"_ptgt20");
    AddVec("jets"+JetKey+"_ptgt20_etalt2p4");
    AddVec("jets"+JetKey+"_ptgt20_etagt2p4lt3p2");
    AddVec("jets"+JetKey+"_ptgt20_etagt3p2lt4p5");

    if(fDoHZZSelection){
        cout << "Warning: rejecting jets overlapping with truth muons and electrons" << endl;
    }
    for(int iJet = 0; iJet < jets(JetKey); iJet++){
        Particle  *myjet        = &(jet(iJet, JetKey));
        // calculate jet area and set area corr pt
        TLorentzVector areaV(0,0,0,0);
        areaV.SetPxPyPzE(myjet->Float("ActiveAreaPx"), myjet->Float("ActiveAreaPy") ,myjet->Float("ActiveAreaPz"), myjet->Float("ActiveAreaE"));
        myjet->Set("Area", areaV.Pt());
        float areacorrpt = (myjet->Float("constscale_pt") - Float("Eventshape_rhoKt4LC")/1000.*myjet->Float("Area") >0) ? myjet->Float("constscale_pt") - Float("Eventshape_rhoKt4LC")/1000.*myjet->Float("Area") : 0;
        myjet->Set("areacorr_pt", areacorrpt);

        // apply cuts
        if((myjet->p.Pt() < 20 || myjet->Int("isBadLoose")==true) && fisMuScan==false) continue; // cut on pt>20 and badloose in case of no mu-scan
        bool overlap(false);
        if(Exists("muonsgood")){
            for(int iMu=0; iMu < muons("good"); iMu++){ // Muon - Jet Overlap ---------
                float dR = muon(iMu,"good").p.DeltaR(myjet->p);
                if(dR<0.4) overlap=true;
            }
        }
        if(fDoHZZSelection){ // reject jets overlapping with truth muons
            for(int iT=0; iT<truths("ZDaugthers"); ++iT){
              float dR = truth(iT,"ZDaugthers").p.DeltaR(myjet->p);
              if(dR<0.4) overlap=true; 
            }
        }
        if(overlap) continue;                       // ---------------------------

        // make collections
        Add("jets"+JetKey+"_ptgt20", myjet); 
        if(fabs(myjet->p.Eta())<2.4 )                                   Add("jets"+JetKey+"_ptgt20_etalt2p4", myjet);
        else if(fabs(myjet->p.Eta())>2.4 && fabs(myjet->p.Eta())<3.2)   Add("jets"+JetKey+"_ptgt20_etagt2p4lt3p2", myjet); 
        else if(fabs(myjet->p.Eta())>3.2 && fabs(myjet->p.Eta())<4.5)   Add("jets"+JetKey+"_ptgt20_etagt3p2lt4p5", myjet); 
    }
  if(Debug()) cout <<"Analysis_PileUpStudiesBase::MakeJetCollections End" << endl;
}

///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_PileUpStudiesBase::WorkerTerminate()
{

  // Nothing more

}

