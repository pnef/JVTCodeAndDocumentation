/**************************************************************************
 **
 **   File:         EventBuilder_example.cxx
 **
 **   Description:  EventBuilder example class
 **                 
 **
 **   Author:       B. Butler
 **
 **   Created:      2-23-11
 **   Modified:
 **
 **************************************************************************/
// Use this file as a template for your EventBuilder source file (delete the 
// one produced by TTree->MakeSelector("EventBuilder_myFile")). Only the
// CopyEvent function must be changed.

#define EventBuilder_ntupbtag2012_cxx
/// Type used for internal caching

#include "EventBuilder_ntupbtag2012.h"
#include "AnaConfig.h"
#include "AnalysisBase.h"
#include "TMath.h"
#include "TVector3.h"
#include "PileupReweighting/TPileupReweighting.h"
#include "GoodRunsLists/TGoodRunsList.h"
#include "TKey.h"
#include <iostream>

EventBuilder_ntupbtag2012::EventBuilder_ntupbtag2012(TTree * /*tree*/)
{
  m_trigdefs = 0;
  myGRL = 0;
  calibAntiKt4TopoEM = 0;
  calibAntiKt6TopoEM = 0;
  calibAntiKt4LCTopo = 0;
  calibAntiKt6LCTopo = 0;

  mJETTYPEVEC    = "JETTYPEVEC";
  mBTAGVEC       = "BTAGVEC";
  // nothing to init yet
}

EventBuilder_ntupbtag2012::~EventBuilder_ntupbtag2012(){
  // nothing to delete
  //delete m_trigdefs;
  if(m_trigdefs) delete m_trigdefs;
  // if(calibAntiKt4TopoEM) delete calibAntiKt4TopoEM;
  // if(calibAntiKt6TopoEM) delete calibAntiKt6TopoEM;
  // if(calibAntiKt4LCTopo) delete calibAntiKt4LCTopo;
  // if(calibAntiKt6LCTopo) delete calibAntiKt6LCTopo;
  
}

void EventBuilder_ntupbtag2012::Initialize(){
  //myGRL = (Root::TGoodRunsList*)fInput->FindObject("myGRL");	
	//if(!myGRL) {
		//cout << "EventBuilder_ntupbtag2012: ERROR cannot retrieve GRL object from input list!" << endl;
		//exit(1); // whatever don't care for now
	//}

  // TString JES_config_file="CalibrationConfigs/Rel17_JES.config";
  // calibAntiKt4TopoEM = new JetCalibrationTool("AntiKt4TopoEM",JES_config_file);      
  // calibAntiKt4TopoEM->UseGeV(true);
  // calibAntiKt6TopoEM = new JetCalibrationTool("AntiKt6TopoEM",JES_config_file);      
  // calibAntiKt6TopoEM->UseGeV(true);
  // calibAntiKt4LCTopo = new JetCalibrationTool("AntiKt4LCTopo",JES_config_file);      
  // calibAntiKt4LCTopo->UseGeV(true);
  // calibAntiKt6LCTopo = new JetCalibrationTool("AntiKt6LCTopo",JES_config_file);      
  // calibAntiKt6LCTopo->UseGeV(true);

}

Bool_t EventBuilder_ntupbtag2012::CopyEvent(AnalysisBase* evt)
{

  // Copy the event from the TTree to the event class. Return false if you want to reject
  // event at the EventBuilder level using D3PD snippets, etc. 

  // ignore events with "future" mu profile
  //if(evt->isMC() && RunNumber == 185761) return kFALSE;


  doJet4 = false;
  doJet6 = false;
  doTruth = false;
  doTruthParentChild = false;
  doTruthJets = false;
  doLCJets = false;
  doEMJets = false;
  doTrack = false;
  doLCCluster = false;
  doEMCluster = false;
  doTruthHeavyQ = false;
  doTruthLightQ = false;
  doTruthStable = false;
  doVertex = false;
  doPhotons = false;
  doOffset = false;
  doMuons = false;
  muRef = 0.;
  npvRef = 0.;

  //evt->Cfg()->Get("DOJET4", doJet4);
  //evt->Cfg()->Get("DOJET6", doJet6);
  //evt->Cfg()->Get("DOTRUTHJETS", doTruthJets);
  //evt->Cfg()->Get("DOLCJETS", doLCJets);
  //evt->Cfg()->Get("DOEMJETS", doEMJets);
  evt->Cfg()->Get("DOTRACK", doTrack);
  evt->Cfg()->Get("DOPHOTON", doPhotons);
  evt->Cfg()->Get("DOLCCLUSTER", doLCCluster);
  evt->Cfg()->Get("DOEMCLUSTER", doEMCluster);
  evt->Cfg()->Get("DOTRUTH",doTruth);
  evt->Cfg()->Get("DOTRUTHPARENTCHILD",doTruthParentChild);
  evt->Cfg()->Get("DOTRUTHHEAVYQ",doTruthHeavyQ);
  evt->Cfg()->Get("DOTRUTHLIGHTQ",doTruthLightQ);
  evt->Cfg()->Get("DOTRUTHSTABLE",doTruthStable);
  evt->Cfg()->Get("DOVTX",doVertex);
  evt->Cfg()->Get("DOOFFSET",doOffset);
  evt->Cfg()->Get("DOMUONS",doMuons);

  evt->Cfg()->Get("MUREF",muRef);
  evt->Cfg()->Get("NPVREF",npvRef);
  
  ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  /// Add jet type list to event config
  ///=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if(!evt->Cfg()->Exists(mJETTYPEVEC)) {
    evt->Cfg()->AddVec(mJETTYPEVEC);
    //evt->AddVec(mJETTYPEVEC);
    TObjArray* arr = evt->Cfg()->String("JETTYPES").Tokenize(",");
    for(Int_t i = 0; i<arr->GetEntries(); ++i) {
      if (evt->Debug()) {
        TObjString* jetTypeObj = (TObjString*)arr->At(i);
    TString jetType        = jetTypeObj->GetString();
        cout << "EventBuilder_ntupbtag2012::CopyEvent(): DEBUG Adding jet type: " << jetType << endl; 
      }
      
      evt->Cfg()->Add(mJETTYPEVEC, (TObjString*)arr->At(i));
    }
    
    arr->SetOwner(kFALSE);  delete arr;
  }

  if(!evt->Cfg()->Exists(mBTAGVEC)) {
    evt->Cfg()->AddVec(mBTAGVEC);
    //evt->AddVec(mJETTYPEVEC);
    TObjArray* arr = evt->Cfg()->String("BTAGS").Tokenize(",");
    for(Int_t i = 0; i<arr->GetEntries(); ++i) {
      if (evt->Debug()) {
        TObjString* jetTypeObj = (TObjString*)arr->At(i);
    TString jetType        = jetTypeObj->GetString();
        cout << "EventBuilder_ntupbtag2012::CopyEvent(): DEBUG Adding btag type: " << jetType << endl; 
      }
      
      evt->Cfg()->Add(mBTAGVEC, (TObjString*)arr->At(i));
    }
    
    arr->SetOwner(kFALSE);  delete arr;
  }

  float temp = 0.;

  evt->Set("evt_weight",Get<vector<vector<double> > >("mcevt_weight").at(0).at(0));

  evt->Cfg()->Get("LUMI",temp);
  evt->Set("LUMI",temp);

  //doTruth = doTruthHeavyQ || doTruthLightQ || doTruthStable;

  //cout << "gonna try to ask for a trigger" << endl;
  //cout << m_trigdefs->IsPassed("EF_j20_a4tc_EFFS",TrigDefs::Physics) << endl;
 
  fEvt = evt;

  if(evt->Debug()) cout << " set up configs " << endl;

	//cout << "init global vars" << endl;
	evt->Set("RunNumber", (int) Get<UInt_t>("RunNumber"));
	evt->Set("EventNumber", (int)Get<UInt_t>("EventNumber"));
	evt->Set("EventWeight", 1.0); // default value
	evt->Set("larerror", Get<int>("larError")); 
  //evt->Set("isMC",(bool)(Get<UInt_t>("RunNumber")<152166));
	evt->Set("averageIntPerXing", Get<float>("averageIntPerXing"));
	evt->Set("LBN", (int)Get<UInt_t>("lbn"));
	evt->Set("BunchCrossingID", (int)Get<UInt_t>("bcid"));

  static const BranchKey bmc_channel_number("mc_channel_number");
	static const MomKey misMC("isMC"), mChannelNumber("ChannelNumber");
	if(BranchNames().find(bmc_channel_number)!=BranchNames().end()) {
		evt->Set(mChannelNumber,(int)Get<UInt_t>(bmc_channel_number));
		evt->Set(misMC,true);	
	}
	else evt->Set(misMC,false);	

	//GRL
	static const MomKey mgrl("grl");
	if(fEvt->isMC()) fEvt->Set(mgrl,true);
	else fEvt->Set(mgrl,myGRL->HasRunLumiBlock(Get<UInt_t>("RunNumber"),Get<UInt_t>("lbn")));


  if(evt->Debug()) cout << "set event quantities" << endl;
  
  if(evt->isMC()){
    float maxpt6 = 0.;
    if(Get<vector<float> >("jet_antikt6topoem_pt").size() > 0) maxpt6 = Get<vector<float> >("jet_antikt6topoem_pt").at(0)/1000;   

    int chan = evt->Int("ChannelNumber");
    if(((chan == 105009 || chan == 113204 || chan == 126127) && maxpt6 > 0) || 
        ((chan == 105010 || chan == 113205 || chan == 126128) && maxpt6 > 70) || 
        ((chan == 105011 || chan == 113206 || chan == 126129) && maxpt6 > 140) ||
        ((chan == 105012 || chan == 113207 || chan == 126130) && maxpt6 > 240) || // max : 210 instead of 240
        ((chan == 105013 || chan == 113208 || chan == 126131) && maxpt6 > 440) ||
        ((chan == 105014 || chan == 113209 || chan == 126132) && maxpt6 > 800) ||
        ((chan == 105015 || chan == 113210 || chan == 126133) && maxpt6 > 1300) ||
        ((chan == 105016 || chan == 126134) && maxpt6 > 2400) ){
      if(evt->Debug()) cout << "MC event with high-pt pileup jet!" << endl;
      return kFALSE;
    }

    //double pthat = Get<vector<double> >("mcevt_event_scale").at(0);
    //cout << "pthat = " << pthat << endl;
    double pthat = Get<vector<double> >("mcevt_pdf_scale")[0];

    if((chan == 108087 && (pthat > 45 || pthat < 0)) ||
        (chan == 108081 && (pthat > 85 || pthat < 45)) ||
        (chan == 108082 && (pthat > 150 || pthat < 85)) ||
        (chan == 108083 && (pthat > 310 || pthat < 150)) ||
        (chan == 108084 && (pthat > 1000000 || pthat < 310))){ // NEEDS WORK

        if(evt->Debug()) cout << "killed event with pthat " << pthat << " for channel " << evt->Int("ChannelNumber") << endl;
      return kFALSE;
    } 



    if(chan >= 147910 && chan <= 147917){
      bool status = true;
      if(Get<vector<float> >("jet_antikt4lctopo_pt").size() == 0 || Get<vector<float> >("jet_antikt4lctopo_pt").size() == 1 || Get<vector<float> >("jet_antikt4truth_pt").size() == 0){
        status = true;
      }
      else{
        float ptave = (Get<vector<float> >("jet_antikt4lctopo_pt").at(0)/1000. + Get<vector<float> >("jet_antikt4lctopo_pt").at(1)/1000.)/2.;
        float pttrue = (Get<vector<float> >("jet_antikt4truth_pt").at(0)/1000.);

        if(ptave/pttrue > 0.6 && ptave/pttrue < 1.4)
          status = true;
        else
          status = false;
      }
      if(!status){
        if(fEvt->Debug()) cout << "found and killed a bad jet!" << endl;

        return kFALSE;

      }
    }



  }
  if(evt->Debug()) cout << "killed high-pt pileup jets" << endl;

  bool fail = true;

  fail = fail && CopyVertices();
  fail = fail && CopyClusters();
  fail = fail && CopyTracks();
  fail = fail && CopyMuons();

  fail = fail && CopyTruthFull();
  fail = fail && AddTruthParentChild();

  fail = fail && CopyJets();
  fail = fail && AddBtags();
  fail = fail && CopyPhotons();
  evt->SortPtAll(); // above need to be sorted by pt, Truth done by E	
  //if(evt->isMC()) fail = fail && CopyTruth();
  
  if(fEvt->Debug()) cout << "abandon all hope, ye who enter here" << endl; 
  fail = fail && DoCustomSetup();
  
  if(fEvt->Debug()) cout << "set up full event, status is " << fail << endl;
	return fail;
}


///=========================================
/// CopyTruthFull
///=========================================
Bool_t EventBuilder_ntupbtag2012::CopyTruthFull()
{
  
  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012: DEBUG CopyTruthFull()  and doTruth = " << doTruth << endl;  
  
  /// If scheduled, set the truth particles
  if(fEvt->isMC() && doTruth) {
    
    fEvt->AddVec("truths");
    
    for (unsigned int i = 0; i< Get<UInt_t>("mcpart_n"); i++) {

      Particle * mc = new Particle();
      float pt  = Get<vector<float> >("mcpart_pt") .at(i);
      float eta = Get<vector<float> >("mcpart_eta").at(i);
      float phi = Get<vector<float> >("mcpart_phi").at(i);
      float m   = Get<vector<float> >("mcpart_m")  .at(i);
      mc->p.SetPtEtaPhiM(pt/1000., eta, phi, m/1000.);
      mc->Set("index",i);      
      mc->Set("type",  Get<vector<int> >("mcpart_type") .at(i)); 
      mc->Set("status", Get<vector<int> >("mcpart_status").at(i)); 
     // mc->Set("charge", Get<vector<float> >("mcpart_charge").at(i));
      mc->Set("barcode", Get<vector<int> >("mcpart_barcode").at(i));


      fEvt->Add("truths",mc);
      if(fEvt->Debug()){
        //cout << "Going to print particle " << i << endl;
        //mc->Show();
      }
    }
  }
  
  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::CopyTruthFull(): DEBUG Scheduled Truth Particles" << endl;
  return kTRUE;
}

///=========================================
/// AddTruthParentChild
///=========================================
Bool_t EventBuilder_ntupbtag2012::AddTruthParentChild()
{
  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::AddTruthParentChild(): DEBUG Starting" << endl;

  if(doTruthParentChild)
  {
    for(int iTr = 0; iTr < fEvt->truths(); iTr++)
    {
      //if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::AddTruthParentChild(): DEBUG On particle " << iTr << " of " << fEvt->truths() << endl;

      fEvt->truth(iTr).AddVec("children");
      fEvt->truth(iTr).AddVec("parents", true);
      //if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::AddTruthParentChild(): DEBUG Has children " << Get<vector<int> >("mcpart_child_n").at(iTr)  << endl;

      for(int iChild = 0; iChild < Get<vector<int> >("mcpart_child_n").at(iTr); iChild++)
      {
        int index = Get<vector<vector<int> > >("mcpart_child_index").at(iTr).at(iChild);
        if(index > 0)
          fEvt->truth(iTr).Add("children",&(fEvt->truth(index)));
      }

      //if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::AddTruthParentChild(): DEBUG Has parents " << Get<vector<int> >("mcpart_mother_n").at(iTr)  << endl;

      for(int iParent = 0; iParent < Get<vector<int> >("mcpart_mother_n").at(iTr); iParent++)
      {
        int index = Get<vector<vector<int> > >("mcpart_mother_index").at(iTr).at(iParent);
        //if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::AddTruthParentChild(): DEBUG Adding parent " << index << endl;
        if(index > 0)
          fEvt->truth(iTr).Add("parents",&(fEvt->truth(index)));
      }
    }
  } 
  return true;
  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::AddTruthParentChild(): DEBUG Ending" << endl;

}

Bool_t EventBuilder_ntupbtag2012::CopyTruth(){
    //if(evt->Debug()) cout << "size of turth4 = " << evt->jets("truth_akt4") << endl;
    if(doTruth){
      if(doTruthHeavyQ) fEvt->AddVec("truthsHeavyQ");
      if(doTruthLightQ) fEvt->AddVec("truthsLightQ");
      if(doTruthStable) fEvt->AddVec("truthsStable");

      for(unsigned int i = 0; i< Get<vector<int> >("mc_pdgId").size(); i++){

        bool stable = true;
        bool pHeavy = false; // init different from stable
        bool pLight = false;

        if(fabs(Get<vector<float > >("mc_eta").at(i)) > 5.5){
          continue; // kill the far far forward stuff
        }
        int mc_pdgId = Get<vector<int> >("mc_pdgId").at(i);
        int mc_status = Get<vector<int> >("mc_status").at(i);

        if(mc_status % 1000 != 1) stable = false;

        if(fabs(mc_pdgId) == 13 || 
            fabs(mc_pdgId) == 12 ||
            fabs(mc_pdgId) == 14 ||
            fabs(mc_pdgId) == 16)  stable = false;

        if(fabs(mc_pdgId) == 1 ||
            fabs(mc_pdgId) == 2 ||
            fabs(mc_pdgId) == 3 ||
            fabs(mc_pdgId) == 9 ||
            fabs(mc_pdgId) == 21) pLight = true;

        if(fabs(mc_pdgId) == 4 ||
            fabs(mc_pdgId) == 5 ||
            fabs(mc_pdgId) == 6 ) pHeavy = true;

        if(fEvt->Debug()){
          if(pHeavy && pLight && stable) cout << "Major problem: all three mc on!" << endl;
        }


        if(!(stable || pHeavy || pLight)) continue; 

        float pt = Get<vector<float> >("mc_pt").at(i)/1000.;
        float eta = Get<vector<float> >("mc_eta").at(i);
        float phi = Get<vector<float> >("mc_phi").at(i);
        float m = Get<vector<float> >("mc_m").at(i)/1000.;
        if(pt < 0.001) continue;

        Particle * mc = new Particle();
        mc->Set("mc_pdgId", mc_pdgId); 
        mc->Set("mc_status", mc_status); 

        mc->p.SetPtEtaPhiM(pt, eta, phi, m);
        //cout << "setting up mc " << mc.p.Perp() << endl;
        if(stable && doTruthStable) fEvt->Add("truthsStable", mc); 
        else if(pHeavy && doTruthHeavyQ) fEvt->Add("truthsHeavyQ", mc);
        else if(pLight && doTruthLightQ) fEvt->Add("truthsLightQ", mc);  
        else delete mc;
      } 

    if(doTruthLightQ)fEvt->SortE("truthsHeavyQ");
    if(doTruthHeavyQ)fEvt->SortE("truthsLightQ");
    if(doTruthStable)fEvt->SortE("truthsStable");
	}
  if(fEvt->Debug()) cout << "set up truth" << endl;
  return kTRUE;
}

/*Bool_t EventBuilder_ntupbtag2012::CopyJets(){

  

  if(doLCJets){
    if(doJet4){
      


      fEvt->AddVec("jetsakt4LCTopo");
      for(unsigned int iJet = 0; iJet < Get<vector<float> >("jet_AntiKt4LCTopo_phi").size(); iJet++){
        float pt = Get<vector<float> >("jet_AntiKt4LCTopo_pt").at(iJet);
        float eta = Get<vector<float> >("jet_AntiKt4LCTopo_eta").at(iJet);
        float phi = Get<vector<float> >("jet_AntiKt4LCTopo_phi").at(iJet);
        float e = Get<vector<float> >("jet_AntiKt4LCTopo_E").at(iJet);

        Particle* jet = new Particle();
        jet->p.SetPtEtaPhiE(pt/1000.,eta,phi,e/1000.);
        jet->Set("jvf",Get<vector<float> >("jet_AntiKt4LCTopo_jvtxf").at(iJet));	
        jet->Set("isUgly",Get<vector<UInt_t> >("jet_AntiKt4LCTopo_isUgly").at(iJet));
        jet->Set("isBadLoose",Get<vector<UInt_t> >("jet_AntiKt4LCTopo_isBadLoose").at(iJet));
        jet->Set("isBadLooser",Get<vector<UInt_t> >("jet_AntiKt4LCTopo_isBadLooser").at(iJet));
        jet->Set("ptConst",Get<vector<float> >("jet_AntiKt4LCTopo_constscale_pt").at(iJet));
        if(fEvt->isMC()) jet->Set("truth_label",Get<vector<int> >("jet_AntiKt4LCTopo_flavor_truth_label").at(iJet));
        fEvt->Add("jetsakt4LCTopo",jet);
      }
    }
    if(doJet6){ 
      fEvt->AddVec("jetsakt6LCTopo");
      for(unsigned int iJet = 0; iJet < Get<vector<float> >("jet_AntiKt6LCTopo_phi").size(); iJet++){
        float pt = Get<vector<float> >("jet_AntiKt6LCTopo_pt").at(iJet);
        float eta = Get<vector<float> >("jet_AntiKt6LCTopo_eta").at(iJet);
        float phi = Get<vector<float> >("jet_AntiKt6LCTopo_phi").at(iJet);
        float e = Get<vector<float> >("jet_AntiKt6LCTopo_E").at(iJet);

        Particle* jet = new Particle();
        jet->p.SetPtEtaPhiE(pt/1000.,eta,phi,e/1000.);
        jet->Set("jvf",Get<vector<float> >("jet_AntiKt6LCTopo_jvtxf").at(iJet));	
        jet->Set("isUgly",Get<vector<UInt_t> >("jet_AntiKt6LCTopo_isUgly").at(iJet));
        jet->Set("isBadLoose",Get<vector<UInt_t> >("jet_AntiKt6LCTopo_isBadLoose").at(iJet));
        jet->Set("isBadLooser",Get<vector<UInt_t> >("jet_AntiKt6LCTopo_isBadLooser").at(iJet));
        jet->Set("ptConst",Get<vector<float> >("jet_AntiKt6LCTopo_constscale_pt").at(iJet));
        if(fEvt->isMC()) jet->Set("truth_label",Get<vector<int> >("jet_AntiKt6LCTopo_flavor_truth_label").at(iJet));
        fEvt->Add("jetsakt6LCTopo",jet);
      }
    }
  }

  if(doEMJets){
    if(doJet4){ 
      fEvt->AddVec("jetsakt4TopoEM");
      for(unsigned int iJet = 0; iJet < Get<vector<float> >("jet_AntiKt4TopoEM_phi").size(); iJet++){
        Particle* jet = new Particle();
       
        if(!doOffset){ 
          float pt = Get<vector<float> >("jet_AntiKt4TopoEM_pt").at(iJet);
          float eta = Get<vector<float> >("jet_AntiKt4TopoEM_eta").at(iJet);
          float phi = Get<vector<float> >("jet_AntiKt4TopoEM_phi").at(iJet);
          float e = Get<vector<float> >("jet_AntiKt4TopoEM_E").at(iJet);

          jet->p.SetPtEtaPhiE(pt/1000.,eta,phi,e/1000.);
        }
        else{
          //float pt_const = Get<vector<float> >("jet_AntiKt4TopoEM_constscale_pt").at(iJet);
          float eta_const = Get<vector<float> >("jet_AntiKt4TopoEM_constscale_eta").at(iJet);
          //float phi_const = Get<vector<float> >("jet_AntiKt4TopoEM_constscale_phi").at(iJet);
          float e_const = Get<vector<float> >("jet_AntiKt4TopoEM_constscale_E").at(iJet)/1000.;

          float eta_origin = Get<vector<float> >("jet_AntiKt4TopoEM_EtaOrigin").at(iJet);
          float phi_origin= Get<vector<float> >("jet_AntiKt4TopoEM_PhiOrigin").at(iJet);
          float m_origin= Get<vector<float> >("jet_AntiKt4TopoEM_MOrigin").at(iJet)/1000.;

          jet->p = calibAntiKt4TopoEM->ApplyOffsetEtaJES(e_const, eta_const, eta_origin, phi_origin, m_origin, fEvt->Float("averageIntPerXing"), fEvt->Float("NPV")); 
        }

        
  
        jet->Set("jvf",Get<vector<float> >("jet_AntiKt4TopoEM_jvtxf").at(iJet));	
        jet->Set("isUgly",Get<vector<UInt_t> >("jet_AntiKt4TopoEM_isUgly").at(iJet));
        jet->Set("isBadLoose",Get<vector<UInt_t> >("jet_AntiKt4TopoEM_isBadLoose").at(iJet));
        jet->Set("ptConst",Get<vector<float> >("jet_AntiKt4TopoEM_constscale_pt").at(iJet));
        if(fEvt->isMC()) jet->Set("truth_label",Get<vector<int> >("jet_AntiKt4TopoEM_flavor_truth_label").at(iJet));
        fEvt->Add("jetsakt4TopoEM",jet);
      }
    }
    if(doJet6){ 
      fEvt->AddVec("jetsakt6TopoEM");
      for(unsigned int iJet = 0; iJet < Get<vector<float> >("jet_AntiKt6TopoEM_phi").size(); iJet++){
        Particle* jet = new Particle();
        if(!doOffset){ 
          float pt = Get<vector<float> >("jet_AntiKt6TopoEM_pt").at(iJet);
          float eta = Get<vector<float> >("jet_AntiKt6TopoEM_eta").at(iJet);
          float phi = Get<vector<float> >("jet_AntiKt6TopoEM_phi").at(iJet);
          float e = Get<vector<float> >("jet_AntiKt6TopoEM_E").at(iJet);

          jet->p.SetPtEtaPhiE(pt/1000.,eta,phi,e/1000.);
        }
        else{
          //float pt_const = Get<vector<float> >("jet_AntiKt6TopoEM_constscale_pt").at(iJet);
          float eta_const = Get<vector<float> >("jet_AntiKt6TopoEM_constscale_eta").at(iJet);
          //float phi_const = Get<vector<float> >("jet_AntiKt6TopoEM_constscale_phi").at(iJet);
          float e_const = Get<vector<float> >("jet_AntiKt6TopoEM_constscale_E").at(iJet)/1000.;

          float eta_origin = Get<vector<float> >("jet_AntiKt6TopoEM_EtaOrigin").at(iJet);
          float phi_origin= Get<vector<float> >("jet_AntiKt6TopoEM_PhiOrigin").at(iJet);
          float m_origin= Get<vector<float> >("jet_AntiKt6TopoEM_MOrigin").at(iJet)/1000.;

          jet->p = calibAntiKt6TopoEM->ApplyOffsetEtaJES(e_const, eta_const, eta_origin, phi_origin, m_origin, fEvt->Float("averageIntPerXing"), fEvt->Float("NPV")); 
        }



        jet->Set("jvf",Get<vector<float> >("jet_AntiKt6TopoEM_jvtxf").at(iJet));	
        jet->Set("isUgly",Get<vector<UInt_t> >("jet_AntiKt6TopoEM_isUgly").at(iJet));
        jet->Set("isBadLoose",Get<vector<UInt_t> >("jet_AntiKt6TopoEM_isBadLoose").at(iJet));
        jet->Set("ptConst",Get<vector<float> >("jet_AntiKt6TopoEM_constscale_pt").at(iJet));
        if(fEvt->isMC()) jet->Set("truth_label",Get<vector<int> >("jet_AntiKt6TopoEM_flavor_truth_label").at(iJet));
        fEvt->Add("jetsakt6TopoEM",jet);
      }
    }
  }

  if(doTruthJets && fEvt->isMC()){
    if(doJet4){ 
      fEvt->AddVec("jetsakt4Truth");
      for(unsigned int iJet = 0; iJet < Get<vector<float> >("AntiKt4TruthNew_phi").size(); iJet++){
        float pt = Get<vector<float> >("AntiKt4TruthNew_pt").at(iJet);
        float eta = Get<vector<float> >("AntiKt4TruthNew_eta").at(iJet);
        float phi = Get<vector<float> >("AntiKt4TruthNew_phi").at(iJet);
        float e = Get<vector<float> >("AntiKt4TruthNew_E").at(iJet);

        Particle* jet = new Particle();
        jet->p.SetPtEtaPhiE(pt/1000.,eta,phi,e/1000.);
        jet->Set("ptConst",Get<vector<float> >("AntiKt4TruthNew_pt").at(iJet));
        fEvt->Add("jetsakt4Truth",jet);
      }
    }
    if(doJet6){ 
      fEvt->AddVec("jetsakt6Truth");
      for(unsigned int iJet = 0; iJet < Get<vector<float> >("AntiKt6TruthNew_phi").size(); iJet++){
        float pt = Get<vector<float> >("AntiKt6TruthNew_pt").at(iJet);
        float eta = Get<vector<float> >("AntiKt6TruthNew_eta").at(iJet);
        float phi = Get<vector<float> >("AntiKt6TruthNew_phi").at(iJet);
        float e = Get<vector<float> >("AntiKt6TruthNew_E").at(iJet);

        Particle* jet = new Particle();
        jet->Set("ptConst",Get<vector<float> >("AntiKt6TruthNew_pt").at(iJet));
        jet->p.SetPtEtaPhiE(pt/1000.,eta,phi,e/1000.);
        fEvt->Add("jetsakt6Truth",jet);
      }
    }
  }
  if(fEvt->Debug()) cout << "set up jets" << endl;
  return kTRUE;
}*/

Bool_t EventBuilder_ntupbtag2012::CopyClusters(){
  if(doLCCluster){
    fEvt->AddVec("clustersLCTopo");
    for(unsigned int iCl = 0; iCl < Get<UInt_t>("cl_had_n"); iCl++){
      Particle* cl = new Particle();
      float pt = Get<vector<float> >("cl_had_pt").at(iCl);
      float eta = Get<vector<float> >("cl_had_eta").at(iCl);
      float phi = Get<vector<float> >("cl_had_phi").at(iCl);
      cl->p.SetPtEtaPhiM(pt/1000.,eta,phi,0);
      fEvt->Add("clustersLCTopo",cl);
    }
  }
  if(fEvt->Debug()) cout << "set up clusters" << endl;
  return kTRUE;
}

Bool_t EventBuilder_ntupbtag2012::CopyTracks(){
  if(doTrack){
    fEvt->AddVec("tracks");
    for(unsigned int iTr = 0; iTr < Get<UInt_t>("trk_n"); iTr++){
      Particle* tr = new Particle();
      float pt = Get<vector<float> >("trk_pt").at(iTr);
      float eta = Get<vector<float> >("trk_eta").at(iTr);
      float phi = Get<vector<float> >("trk_phi").at(iTr);
      tr->p.SetPtEtaPhiM(pt/1000.,eta,phi,0);
      tr->Set("nPixHits", Get<vector<int> >("trk_nPixHits").at(iTr));
      tr->Set("nPixHoles", Get<vector<int> >("trk_nPixHoles").at(iTr));
      tr->Set("nPixelDeadSensors", Get<vector<int> >("trk_nPixelDeadSensors").at(iTr));
      tr->Set("nSCTHits", Get<vector<int> >("trk_nSCTHits").at(iTr));
      tr->Set("nSCTHoles", Get<vector<int> >("trk_nSCTHoles").at(iTr));
      tr->Set("nSCTDeadSensors", Get<vector<int> >("trk_nSCTDeadSensors").at(iTr));
      tr->Set("d0", Get<vector<float> >("trk_d0_wrtPV").at(iTr));
      tr->Set("z0", Get<vector<float> >("trk_z0_wrtPV").at(iTr));
      tr->Set("chi2", Get<vector<float> >("trk_chi2").at(iTr));
      tr->Set("theta", Get<vector<float> >("trk_theta").at(iTr));
      tr->Set("ndof", Get<vector<int> >("trk_ndof").at(iTr));
      tr->Set("nTRTOutliers", Get<vector<int> >("trk_nTRTOutliers").at(iTr));
      tr->Set("nTRTHits", Get<vector<int> >("trk_nTRTHits").at(iTr));

      fEvt->Add("tracks",tr);
    }
  }
  if(fEvt->Debug()) cout << "set up tracks" << endl;
  return kTRUE;
}

Bool_t EventBuilder_ntupbtag2012::CopyVertices()
{
  if(doVertex){
    fEvt->AddVec("vtxs");

    //Vertices
    int NPV=0;
    for(unsigned pvi=0;pvi<Get<vector<int> >("primvx_nTracks").size();++pvi) {
      if(Get<vector<int> >("primvx_nTracks").at(pvi) >= 2){
        ++NPV;
        Point* vtx = new Point();
        vtx->x.SetXYZ(Get<vector<float> >("primvx_x").at(pvi),Get<vector<float> >("primvx_y").at(pvi),Get<vector<float> >("primvx_z").at(pvi));
        vtx->Set("ntrk",Get<vector<int> >("primvx_nTracks").at(pvi));

        fEvt->Add("vtxs",vtx);
      }
    }
    fEvt->Set("NPV",NPV);
  }
  if(fEvt->Debug()) cout << "set up vtx" << endl;
  return kTRUE;
}

Bool_t EventBuilder_ntupbtag2012::CopyPhotons(){
  if(doPhotons){
    fEvt->AddVec("photons");
    for(unsigned int iPh = 0; iPh < Get<vector<float> >("ph_pt").size(); iPh++){
      float pt = Get<vector<float> >("ph_pt").at(iPh);
      float eta = Get<vector<float> >("ph_eta").at(iPh);
      float phi = Get<vector<float> >("ph_phi").at(iPh);
      float e = Get<vector<float> >("ph_E").at(iPh);      

      Particle* ph = new Particle();
      ph->p.SetPtEtaPhiE(pt/1000.,eta,phi,e/1000.);

      ph->Set("Etcone40_corrected", Get<vector<float> >("ph_Etcone40_corrected").at(iPh)/1000.);
     // ph->Set("Etcone40_ED_corrected", Get<vector<float> >("ph_Etcone40_ED_corrected").at(iPh));
      //ph->Set("ptIsolationCone", Get<vector<float> >("ph_ptIsolationCone").at(iPh));
      //ph->Set("topoEtcone40", Get<vector<float> >("ph_topoEtcone40").at(iPh));
      //ph->Set("Etcone40_pt_corrected", Get<vector<float> >("ph_Etcone40_pt_corrected").at(iPh));
      
      ph->Set("tight", Get<vector<int> >("ph_tight").at(iPh));
      ph->Set("tightIso", Get<vector<int> >("ph_tightIso").at(iPh));
      ph->Set("loose", Get<vector<int> >("ph_loose").at(iPh));
      ph->Set("looseIso", Get<vector<int> >("ph_looseIso").at(iPh));

      ph->Set("isEM", Get<vector<UInt_t> >("ph_isEM").at(iPh));
      ph->Set("isEMLoose", Get<vector<UInt_t> >("ph_isEMLoose").at(iPh));
      ph->Set("isEMMedium", Get<vector<UInt_t> >("ph_isEMMedium").at(iPh));
      ph->Set("isEMTight", Get<vector<UInt_t> >("ph_isEMTight").at(iPh));

      ph->Set("OQ", Get<vector<UInt_t> >("ph_OQ").at(iPh));
      //ph->Set("OQRecalc", Get<vector<UInt_t> >("ph_OQRecalc").at(iPh));

      fEvt->Add("photons",ph);
    }
  }
  if(fEvt->Debug()) cout << "set up photons" << endl;
  return kTRUE;
}


void EventBuilder_ntupbtag2012::WriteAuxTrees(TDirectory* outfile)
{
	if(!m_isNewInit) return;
	m_isNewInit = false;

	TFile* file = 0;
	TChain* chain = dynamic_cast< TChain* >( Tree() );
	if( chain ) {
		// We are running locally...
		file = chain->GetFile();
	} else {
		// We are running on PROOF:
		file = Tree()->GetCurrentFile();
	}
	
	//Thie following doesn't work, as TMNUtil modifies the tree in some way to make it unclonable.
	//Instead, we have to do some annoying stuff to get the unadulterated disk copy
	//TTree* confTree = dynamic_cast< TTree* >( file->Get( "susyMeta/TrigConfTree" ) ); 
	
	TDirectory* origQCDMeta = (TDirectory*)file->Get("qcdMeta");
	TKey *key;
	TTree* confTree = 0;
   	TIter nextkey(origQCDMeta->GetListOfKeys());
   	while ((key = (TKey *) nextkey())) {
   		if (strcmp("TrigConfTree",key->GetName()) == 0) { //This should grab highest cycle, if there are more than one
			TDirectory::TContext ctxt(origQCDMeta);
			confTree = (TTree*)key->ReadObj();
			break;
		}
	}
	if(!confTree) Abort("EventBuilder_qcdMeta: ERROR cannot retrieve TrigConfTree from input file!");

	TDirectory* qcdMeta = (TDirectory*)outfile->Get("qcdMeta");
	if(!qcdMeta) {
		qcdMeta = outfile->mkdir("qcdMeta");
	}
	qcdMeta->cd();

	TTree* TrigConfTree = (TTree*)qcdMeta->Get("TrigConfTree"); 
	if(!TrigConfTree) {
		TrigConfTree = confTree->CloneTree();
	}
	else { //Merge trees
		TList *list = new TList;
		list->Add(TrigConfTree);
		list->Add(confTree);

		TTree *newtree = TTree::MergeTrees(list);

		TrigConfTree->SetDirectory(0);
		delete TrigConfTree; //make sure we load the new one next time
		
		newtree->SetName("TrigConfTree");
		
		delete list;
	}
	delete confTree;
}

///=========================================
/// AddBtags
///=========================================
Bool_t EventBuilder_ntupbtag2012::AddBtags() 
{

  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012: DEBUG In AddBtags() " << endl;   
  bool fail = true;
  /// Copy all jet types
  for(Int_t i = 0; i<fEvt->Cfg()->Objs(mBTAGVEC); ++i) {
    TString type = ((TObjString*)fEvt->Cfg()->Obj(mBTAGVEC,i))->GetString();
    if (type.Contains("Truth")) continue; 
    if(fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::AddBtags(): DEBUG Copying btag type " << type << " which is the " << i << "th jet collection" << endl;
    if(fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::AddBtagss(): DEBUG type in vec is " << ((TObjString*) fEvt->Cfg()->Obj(mBTAGVEC,i))->GetString() << endl;
    fail = fail && AddBtag(type,"");
    //delete type;
  }
  
  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::AddBtags(): DEBUG Scheduled All Jet B-tags" << endl;
  return fail;
}

///=========================================
/// AddBtag
///=========================================
Bool_t EventBuilder_ntupbtag2012::AddBtag(TString jetType, TString jetClass) 
{

  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012: DEBUG In AddBtag() for jetType = " << jetType << "; jetClass = " << jetClass << endl;  
  
  // Handle the case where this is called without an argument
  if (jetType == "") { 
    Abort("EventBuilder_ntupbtag2012::AddBtag(): ERROR No jet type selected. Aborting!");
  }
  
  // Now try to figure out what type of jet this is
  if (jetClass == "") {
    if (fEvt->Debug()) {
      cout << "EventBuilder_ntupbtag2012::AddBtag(): INFO The jetClass variable is empty. ";
      cout << "Assuming that the class (truth, track, topo, etc) can be determined from the name." << endl;
    } // END IF
    
    // Parse the jetType: truth, track, topo
    if (jetType.Contains("Truth")) {
     
      cout << "EventBuilder_ntupbtag2012::AddBtag(): ERROR Cannot add B-tagging for Truth jet collections!" << endl;
      return kFALSE;
    
    }
    else if (jetType.Contains("Track"))  jetClass = "track";
    else if (jetType.Contains("LCTopo")) jetClass = "topo";
    else if (jetType.Contains("TopoEM")) jetClass = "topo";
    else if (jetType.Contains("Tower"))  {
      cout << "EventBuilder_ntupbtag2012::AddBtag(): ERROR Tower jets are not currently supported." << endl;
      return kFALSE;
    } // END IF
    
  } // END IF JET CLASS 
  
  // Otherwise, try to figure out what class of jet this really is
  else if (jetClass != "topo" && jetClass != "track") {
    cout << "EventBuilder_ntupbtag2012::AddBtag(): ERROR Currently, the only supported classes of jets are: truth, track, topo." << endl;
    return kFALSE;
  }
  
  // Collect jet flavor information
  for(unsigned int iJet = 0; iJet < Get<vector<float> >("jet_" + jetType + "_phi").size(); iJet++){
  
  // get the pointer to the jet already created
  Particle* jet = &fEvt->jet(iJet,jetType);  
  /// Spit out some useful information about this jet
  if (fEvt->Debug()) {
    cout << "EventBuilder_ntupbtag2012::AddBtag(): DEBUG ";
    TString info = 
    TString::Format("Jet %d: (E, Et, m, eta, phi) = (%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)",
            iJet, 
            jet->p.E(),jet->p.Et(),jet->p.M(),jet->p.Eta(),jet->p.Phi()
          );
   
    cout << info << endl;
  }
    
    jet->Set("flavor_weight_Comb",                   Get<vector<float> >("jet_" + jetType + "_flavor_weight_Comb")                  .at(iJet));
    jet->Set("flavor_weight_IP2D",                   Get<vector<float> >("jet_" + jetType + "_flavor_weight_IP2D")                  .at(iJet));
    jet->Set("flavor_weight_IP3D",                   Get<vector<float> >("jet_" + jetType + "_flavor_weight_IP3D")                  .at(iJet));
    jet->Set("flavor_weight_SV0",                    Get<vector<float> >("jet_" + jetType + "_flavor_weight_SV0")                   .at(iJet));
    jet->Set("flavor_weight_SV1",                    Get<vector<float> >("jet_" + jetType + "_flavor_weight_SV1")                   .at(iJet));
    jet->Set("flavor_weight_SV2",                    Get<vector<float> >("jet_" + jetType + "_flavor_weight_SV2")                   .at(iJet));
    jet->Set("flavor_weight_JetFitterTagNN",         Get<vector<float> >("jet_" + jetType + "_flavor_weight_JetFitterTagNN")        .at(iJet));
    jet->Set("flavor_weight_JetFitterCOMBNN",        Get<vector<float> >("jet_" + jetType + "_flavor_weight_JetFitterCOMBNN")       .at(iJet));
    jet->Set("flavor_weight_MV1",                    Get<vector<float> >("jet_" + jetType + "_flavor_weight_MV1")                   .at(iJet));
    jet->Set("flavor_weight_MV2",                    Get<vector<float> >("jet_" + jetType + "_flavor_weight_MV2")                   .at(iJet));
    jet->Set("flavor_weight_SoftMuonTag",            Get<vector<float> >("jet_" + jetType + "_flavor_weight_SoftMuonTag")           .at(iJet));
    jet->Set("flavor_weight_SecondSoftMuonTag",      Get<vector<float> >("jet_" + jetType + "_flavor_weight_SecondSoftMuonTag")     .at(iJet)); 
    jet->Set("flavor_weight_SoftMuonTagChi2",        Get<vector<float> >("jet_" + jetType + "_flavor_weight_SoftMuonTagChi2")       .at(iJet));
    jet->Set("flavor_weight_SecondSoftMuonTagChi2",  Get<vector<float> >("jet_" + jetType + "_flavor_weight_SecondSoftMuonTagChi2") .at(iJet));

    jet->Set("svp_x",                                Get<vector<float> >("jet_" + jetType + "_flavor_component_svp_x")              .at(iJet));
    jet->Set("svp_y",                                Get<vector<float> >("jet_" + jetType + "_flavor_component_svp_y")              .at(iJet));
    jet->Set("svp_z",                                Get<vector<float> >("jet_" + jetType + "_flavor_component_svp_z")              .at(iJet));


    if(fEvt->isMC()){
      jet->Set("flavor_truth_label",                 Get<vector<int> >  ("jet_" + jetType + "_flavor_truth_label")                  .at(iJet));
    }

    if(doTrack){ // now copy over the track assignments to the SV
      jet->AddVec("svp_track"); //flavor_component_svp_trk_index
      vector<int> trackIndices = Get<vector<vector<int> > >("jet_" + jetType + "_flavor_component_svp_trk_index")                  .at(iJet);
      if(fEvt->Debug()) cout << "Adding svp_track information, there are " << trackIndices.size() << " of them!" << endl;

      for(unsigned int iTr = 0; iTr < trackIndices.size(); iTr++){
        if(fEvt->Debug()) cout << "Adding svp_track information, iTr = " << iTr << " and trackIndices = " << trackIndices[iTr] << endl;

        jet->Add("svp_track",&fEvt->track(trackIndices[iTr]));
      }
    }
  
  } // END FOR LOOP
  
  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::AddBtag(): DEBUG Scheduled b-tags for jet collection = " << jetType << " with " << fEvt->jets(jetType) << " jets" << endl;
  return kTRUE;  
  
} // END AddBtag


///=========================================
/// CopyJets
///=========================================
Bool_t EventBuilder_ntupbtag2012::CopyJets() 
{

  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012: DEBUG In CopyJets() " << endl;   

  bool fail       = true;
  TString jetType = "";
  
  /// Copy all jet types
  for(Int_t i = 0; i<fEvt->Cfg()->Objs(mJETTYPEVEC); ++i) {
    TString type = ((TObjString*)fEvt->Cfg()->Obj(mJETTYPEVEC,i))->GetString();
    if (type.Contains("Truth") && !fEvt->isMC()) continue; 
    if(fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::CopyJets(): DEBUG Copying jet type " << type << " which is the " << i << "th jet collection" << endl;
    if(fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::CopyJets(): DEBUG type in vec is " << ((TObjString*) fEvt->Cfg()->Obj(mJETTYPEVEC,i))->GetString() << endl;
    fail = fail && CopyJet(type,"");
    //delete type;
  }
  
  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::CopyJets(): DEBUG Scheduled All Jets" << endl;
  return fail;
}


  ///=========================================
/// CopyJet
///=========================================
Bool_t EventBuilder_ntupbtag2012::CopyJet(TString jetType, TString jetClass) 
{

  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012: DEBUG In CopyJet() for jetType = " << jetType << "; jetClass = " << jetClass << endl;  
  
  // Handle the case where this is called without an argument
  if (jetType == "") { 
    Abort("EventBuilder_ntupbtag2012::CopyJet(): ERROR No jet type selected. Aborting!");
  }
  
  // Now try to figure out what type of jet this is
  if (jetClass == "") {
    if (fEvt->Debug()) {
      cout << "EventBuilder_ntupbtag2012::CopyJet(): INFO The jetClass variable is empty. ";
      cout << "Assuming that the class (truth, track, topo, etc) can be determined from the name." << endl;
    } // END IF
    
    // Parse the jetType: truth, track, topo
    if (jetType.Contains("truth")) {
      if (fEvt->isMC()) jetClass = "truth";
      else { 
        cout << "EventBuilder_ntupbtag2012::CopyJet(): ERROR This is not MC and you have requested truth jets. Reconfigure analysis and run again!." << endl;
        return kFALSE;
      }
    }
    else if (jetType.Contains("track"))  jetClass = "track";
    else if (jetType.Contains("lctopo")) jetClass = "topo";
    else if (jetType.Contains("topoem")) jetClass = "topo";
    else if (jetType.Contains("tower"))  {
      cout << "EventBuilder_ntupbtag2012::CopyJet(): ERROR Tower jets are not currently supported." << endl;
      return kFALSE;
    } // END IF
    
  } // END IF JET CLASS 
  
  // Otherwise, try to figure out what class of jet this really is
  else if (jetClass != "truth" && jetClass != "topo" && jetClass != "track") {
    cout << "EventBuilder_ntupbtag2012::CopyJet(): ERROR Currently, the only supported classes of jets are: truth, track, topo." << endl;
    return kFALSE;
  }
  
  // Collect jet information
  fEvt->AddVec("jets" + jetType);
  for(unsigned int iJet = 0; iJet < Get<vector<float> >("jet_" + jetType + "_phi").size(); iJet++){
  
  /// Get Basic kinematics
  float pt  = Get<vector<float> >("jet_" + jetType + "_pt")  .at(iJet);
  float eta = Get<vector<float> >("jet_" + jetType + "_eta") .at(iJet);
  float phi = Get<vector<float> >("jet_" + jetType + "_phi") .at(iJet);
  float e   = Get<vector<float> >("jet_" + jetType + "_E")   .at(iJet);

    /// Set Basic kinematics
  Particle* jet = new Particle();
  jet->p.SetPtEtaPhiE(pt/1000.,eta,phi,e/1000.);
  
  /// Spit out some useful information about this jet
  if (fEvt->Debug()) {
    cout << "EventBuilder_ntupbtag2012::CopyJet(): DEBUG ";
    TString info = 
    TString::Format("Jet %d: (E, Et, m, eta, phi) = (%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)",
            iJet, 
            jet->p.E(),jet->p.Et(),jet->p.M(),jet->p.Eta(),jet->p.Phi()
          );
   
    cout << info << endl;
  }
    
    
  // Handle reconstructed jets
  if (jetClass == "topo" ) {
    jet->Set("ptEMScale", Get<vector<float> >("jet_" + jetType + "_emscale_pt").at(iJet));
    
    jet->Set("JVF",Get<vector<float> >("jet_" + jetType + "_jvtxf").at(iJet));  
        
    // Cleaning information is only present for AKT4 jets
    if (jetType.Contains("antikt4") || jetType.Contains("antikt6")) {
      jet->Set("isUgly",          Get<vector<UInt_t> >("jet_" + jetType + "_isUgly")         .at(iJet));
    jet->Set("isBadLoose",      Get<vector<UInt_t> >("jet_" + jetType + "_isBadLoose")     .at(iJet));
    jet->Set("isBadLooseMinus", Get<vector<UInt_t> >("jet_" + jetType + "_isBadLooseMinus").at(iJet));
    } // ENDIF AntiKt4
    
  } // ENDIF topo

  // Handle track jets
  else if (jetClass == "track" ) { /* Nothing to do */ }
  
    // Handle truth jets
  else if (jetClass == "truth" ) { /* Nothing to do */ }
  
  // Add the jet to the event
  fEvt->Add("jets" + jetType , jet);
  
  } // END FOR LOOP
  
  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::CopyJet(): DEBUG Scheduled jet collection = " << jetType << " with " << fEvt->jets(jetType) << " jets" << endl;
  return kTRUE;  
  
} // END CopyJet


///=========================================
/// CopyMuons
///=========================================
Bool_t EventBuilder_ntupbtag2012::CopyMuons() 
{
  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012: DEBUG CopyMuons() BEGIN " << endl; 

  if(doMuons) {
    fEvt->AddVec("muons");
    for(unsigned int iMu = 0; iMu < Get<vector<float> >("muoninjet_pt").size(); iMu++) {
      
      if(fEvt->Debug()) cout << "EventBuilder_ntupbtag2012: DEBUG CopyMuons(): copying muon number " << iMu << " out of " << Get<vector<float> >("muoninjet_pt").size() << endl;

      Particle* mu = new Particle();
      float pt  = Get<vector<float> >("muoninjet_pt") .at(iMu);
      float eta = Get<vector<float> >("muoninjet_eta").at(iMu);
      float phi = Get<vector<float> >("muoninjet_phi").at(iMu);
      float m   = Get<vector<float> >("muoninjet_m")  .at(iMu);
      mu->p.SetPtEtaPhiM(pt/1000.,eta,phi,m/1000.);
      mu->Set("author",                  Get<vector<int> >("muoninjet_author")                .at(iMu));

      
      // mu->Set("hastrack",                Get<vector<int> >("muoninjet_hastrack")               .at(iMu));

      // mu->Set("charge",                  Get<vector<float> >("muoninjet_charge")               .at(iMu));
      // mu->Set("loose",                    Get<vector<int> >("muoninjet_loose")                 .at(iMu));
      // mu->Set("medium",                   Get<vector<int> >("muoninjet_medium")                .at(iMu));
      // mu->Set("tight",                    Get<vector<int> >("muoninjet_tight")                 .at(iMu));
      mu->Set("eloss",                    Get<vector<float> >("muoninjet_energyLossPar")               .at(iMu)/1000.);
      mu->Set("nOutliersOnTrack",         Get<vector<int> >("muoninjet_nOutliersOnTrack")       .at(iMu));
      mu->Set("nMDTHits",                 Get<vector<int> >("muoninjet_nMDTHits")               .at(iMu));
      mu->Set("nMDTHoles",                 Get<vector<int> >("muoninjet_nMDTHoles")               .at(iMu));
      mu->Set("nCSCEtaHits",                 Get<vector<int> >("muoninjet_nCSCEtaHits")               .at(iMu));
      mu->Set("nCSCEtaHoles",                 Get<vector<int> >("muoninjet_nCSCEtaHoles")               .at(iMu));
      mu->Set("nCSCPhiHits",                 Get<vector<int> >("muoninjet_nCSCPhiHits")               .at(iMu));
      mu->Set("nCSCPhiHoles",                 Get<vector<int> >("muoninjet_nCSCPhiHoles")               .at(iMu));
      mu->Set("nRPCEtaHits",                 Get<vector<int> >("muoninjet_nRPCEtaHits")               .at(iMu));
      mu->Set("nRPCEtaHoles",                 Get<vector<int> >("muoninjet_nRPCEtaHoles")               .at(iMu));
      mu->Set("nRPCPhiHits",                 Get<vector<int> >("muoninjet_nRPCPhiHits")               .at(iMu));
      mu->Set("nRPCPhiHoles",                 Get<vector<int> >("muoninjet_nRPCPhiHoles")               .at(iMu));
      mu->Set("nTGCEtaHits",                 Get<vector<int> >("muoninjet_nTGCEtaHits")               .at(iMu));
      mu->Set("nTGCEtaHoles",                 Get<vector<int> >("muoninjet_nTGCEtaHoles")               .at(iMu));
      mu->Set("nTGCPhiHits",                 Get<vector<int> >("muoninjet_nTGCPhiHits")               .at(iMu));
      mu->Set("nTGCPhiHoles",                 Get<vector<int> >("muoninjet_nTGCPhiHoles")               .at(iMu));

      mu->Set("primtrk_ndof",                 Get<vector<int> >("muoninjet_primtrk_ndof")               .at(iMu));
      mu->Set("primtrk_chi2",                 Get<vector<float> >("muoninjet_primtrk_chi2")               .at(iMu));

      mu->Set("hastrack",                 Get<vector<int> >("muoninjet_me_hastrack")               .at(iMu));
      mu->Set("d0",                 Get<vector<float> >("muoninjet_me_d0")               .at(iMu));
      mu->Set("z0",                 Get<vector<float> >("muoninjet_me_z0")               .at(iMu));

      mu->AddVec("track");

      if(mu->Bool("hastrack") && Get<vector<int> >("muoninjet_trk_index").at(iMu) >= 0){
        mu->Add("track",&(fEvt->track(Get<vector<int> >("muoninjet_trk_index").at(iMu))));          
      }
      else
        mu->Set("hastrack",false);

      fEvt->Add("muons",mu);
    }
  }
  
  if (fEvt->Debug()) cout << "EventBuilder_ntupbtag2012::CopyPhotons(): scheduled photons" << endl;
  return kTRUE;
}
