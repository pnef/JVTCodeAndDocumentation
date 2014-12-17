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

#define EventBuilder_jetmet2012_cxx
/// Type used for internal caching

#include "EventBuilder_jetmet2012.h"
#include "AnaConfig.h"
#include "AnalysisBase.h"
#include "TMath.h"
#include "TVector3.h"
#include "PileupReweighting/TPileupReweighting.h"
#include "GoodRunsLists/TGoodRunsList.h"
#include "TKey.h"
#include <iostream>

 //Overlap removal w/z samples
#include "SUSYTools/HforToolD3PD.h"

// muon mc smearing
#include "MuonMomentumCorrections/SmearingClass.h"

 //calibration tools
#include "ApplyJetCalibration/ApplyJetCalibration.h"
#include "JetUncertainties/MultijetJESUncertaintyProvider.h"
#include "JetResolution/JERProvider.h"
#include "ApplyJetResolutionSmearing/ApplyJetSmearing.h"



#include "MissingETUtility/METUtility.h" 

struct BTagIndex {
  unsigned int b[2];
  unsigned int c[2];
  unsigned int t[2];
  unsigned int l[2];
};

EventBuilder_jetmet2012::EventBuilder_jetmet2012(TTree * /*tree*/)
{

  m_trigdefs = 0;
  myGRL = 0;
  myPRW = 0;

  data4LCCalibrator = 0;
  data4EMCalibrator = 0;

  data6LCCalibrator = 0;
  data6EMCalibrator = 0;

  mc4LCCalibrator = 0;
  mc4EMCalibrator = 0;

  mc6LCCalibrator = 0;
  mc6EMCalibrator = 0;


  afiimc4LCCalibrator = 0;
  afiimc4EMCalibrator = 0;

  afiimc6LCCalibrator = 0;
  afiimc6EMCalibrator = 0;


  myHforTool = new HforToolD3PD();
  myHforTool->initialize("angularbased", 0.4);


  mJETTYPEVEC    = "JETTYPEVEC";
  mBTAGVEC       = "BTAGVEC";
  mPRWVEC        = "PRWVEC";
  // nothing to init yet
}

EventBuilder_jetmet2012::~EventBuilder_jetmet2012(){
  // nothing to delete
  //delete m_trigdefs;
  //if(m_trigdefs) delete m_trigdefs;
  //if(myHforTool) delete myHforTool;
  //if(myGRL) delete myGRL;
  //if(myPRW) delete myPRW;
}

Bool_t EventBuilder_jetmet2012::SetupPRW(){
  if(!fEvt->Cfg()->Exists(mPRWVEC)) {
    fEvt->Cfg()->AddVec(mPRWVEC);
    //fEvt->AddVec(mJETTYPEVEC);
    TObjArray* arr = fEvt->Cfg()->String("PRWTYPES").Tokenize(",");
    for(Int_t i = 0; i<arr->GetEntries(); ++i) {

      TObjString* prwTypeObj = (TObjString*)arr->At(i);
      TString prwType        = prwTypeObj->GetString();
      if (fEvt->Debug()) {
        cout << "EventBuilder_jetmet2012::CopyEvent(): DEBUG Adding prw type: " << prwType << endl; 
      }
      Root::TPileupReweighting* prw = (Root::TPileupReweighting*) fInput->FindObject(prwType);
      prwMap.insert(pair<MomKey, Root::TPileupReweighting*>(MomKey(prwType), prw));
      
      fEvt->Cfg()->Add(mPRWVEC, (TObjString*)arr->At(i));
    }
    
    arr->SetOwner(kFALSE);  delete arr;
  }

  return kTRUE;
}

void EventBuilder_jetmet2012::Initialize(){
  myGRL = (Root::TGoodRunsList*)fInput->FindObject("myGRL");  
  if(!myGRL) {
    cout << "EventBuilder_jetmet2012: ERROR cannot retrieve GRL object from input list!" << endl;
    exit(1);
  }
  GetJetCalibrationTools();

  //MCPSmear = (MuonSmear::SmearingClass*) fInput->FindObject("MCPSmear");

  std::string maindir(gSystem->Getenv("PROOFANADIR"));
  if(maindir==".") maindir+="/libProofAna";
  cout << maindir << " is the maindir we work from! " << endl;
  cout << gSystem->Exec(TString("ls "+maindir).Data()) << endl;
  MCPSmear = new MuonSmear::SmearingClass("Data12","muid","q_pT","Rel17.2Repro", maindir+"/utils/MuonMomentumCorrections/share/");
  //MCPSmear->Initialize("Data12","muid","q_pT","Rel17.2Repro", "../utils/MuonMomentumCorrections/share/");

	m_util = new METUtility;
	 //vector<float> jet_AntiKt4LCTopo_recalibrated_pt ;
	
}

Bool_t EventBuilder_jetmet2012::CopyEvent(AnalysisBase* evt)
{

////////////////////////////////////////////////////////////
//// Initialize configurations for this event
////////////////////////////////////////////////////////////

  doPRW                     = false;
  doAFII                    = false;
  doJet4                    = false;
  doJet6                    = false;
  doTruthJets               = false;
  doOOTtruthJet4            = false;
  doInTimeTruthJet4         = false;
  doLCJets                  = false;
  doEMJets                  = false;
  doJetCalibrations         = false; 
  doSMWZ                    = false;
  doCOMMON                  = false;
  doTruth                   = false;
  doTruthParentChild        = false;
  doTrack                   = false;
  doLCCluster               = false;
  doEMCluster               = false;
  doVertex                  = false;
  doPhotons                 = false;
  doMuons                   = false;
  doMuonTriggers            = false;
  doElectronTriggers        = false;
  doMuonTriggersMatch       = true;
  doJetTriggers             = false;
  doHFOR                    = false;
  doVectorJVFLinks          = false;
  doTruthMFLinks            = false;
  doConstitLinks            = false;
  doMETRefFinal             = false;
  doBarcode                 = false;
  Rename                    = false;
  doTrkJetBTag              = false;
  doLeptonTruth             = false;
  doITKFixes                = false;
  
  muRef                     = 0.;
  npvRef                    = 0.;
  whichsyste                = JetMETSyst::NONE;
  doMuSmear                 = true;
  doMCGRL                   = false;

  evt->Cfg()->Get("DOHFOR",               doHFOR);
  evt->Cfg()->Get("DOJET4",               doJet4); 
  evt->Cfg()->Get("DOJET6",               doJet6); 
  evt->Cfg()->Get("DOTRUTHJETS",          doTruthJets);
  evt->Cfg()->Get("DOOOTTRUTHJET4",       doOOTtruthJet4);
  evt->Cfg()->Get("DOINTIMETRUTHJET4",    doInTimeTruthJet4);
  evt->Cfg()->Get("DOLCJETS",             doLCJets);
  evt->Cfg()->Get("DOEMJETS",             doEMJets);
  evt->Cfg()->Get("PILE",                 doPRW);
  evt->Cfg()->Get("DOVECTORJVFLINKS",     doVectorJVFLinks);
  evt->Cfg()->Get("DOTRUTHMFLINKS",       doTruthMFLinks);
  evt->Cfg()->Get("DOCONSTITLINKS",       doConstitLinks);
  evt->Cfg()->Get("DOMETREFFINAL",        doMETRefFinal);
  evt->Cfg()->Get("DOSMWZ",               doSMWZ);
  evt->Cfg()->Get("DOCOMMON",             doCOMMON);
  evt->Cfg()->Get("DOAFII",               doAFII);;
  evt->Cfg()->Get("DOTRACK",              doTrack);
  evt->Cfg()->Get("DOPHOTON",             doPhotons);
  evt->Cfg()->Get("DOLCCLUSTER",          doLCCluster);
  evt->Cfg()->Get("DOEMCLUSTER",          doEMCluster);
  evt->Cfg()->Get("DOTRUTH",              doTruth);
  evt->Cfg()->Get("DOPARENTCHILD",        doTruthParentChild);
  evt->Cfg()->Get("DOVTX",                doVertex);
  evt->Cfg()->Get("DOJETCALIBRATIONS",    doJetCalibrations);
  evt->Cfg()->Get("DOMUONS",              doMuons);
  evt->Cfg()->Get("DOMUONTRIGGERS",       doMuonTriggers);
  evt->Cfg()->Get("DOELECTRONTRIGGERS",   doElectronTriggers);
  evt->Cfg()->Get("DOMUONTRIGGERSMATCH",  doMuonTriggersMatch);
  evt->Cfg()->Get("DOJETTRIGGERS",        doJetTriggers);
  evt->Cfg()->Get("DOELECTRONS",          doElectrons);
  evt->Cfg()->Get("SYSTTYPE",             whichsyste);
  evt->Cfg()->Get("DOBARCODE",            doBarcode);
  evt->Cfg()->Get("RENAME",               Rename   );
  evt->Cfg()->Get("DOTRKJETBTAG",         doTrkJetBTag);
  evt->Cfg()->Get("DOLEPTONTRUTH",        doLeptonTruth);
  evt->Cfg()->Get("DOMUSMEAR",            doMuSmear);
  evt->Cfg()->Get("DOMCGRL",              doMCGRL);
  evt->Cfg()->Get("DOITKFIXES",           doITKFixes);

  //cout << "whichsyste is " << whichsyste << endl;

  evt->Cfg()->Get("MUREF",                muRef);
  evt->Cfg()->Get("NPVREF",               npvRef);
  
  fEvt = evt;

  if(evt->Debug()) cout << " set up configs " << endl;

////////////////////////////////////////////////////////////
//// Set some global variable
////////////////////////////////////////////////////////////
	
	jet_AntiKt4LCTopo_recalibrated_pt.clear();

	//cout << "init global vars" << endl;
	evt->Set("RunNumber", (int) Get<UInt_t>("RunNumber"));
	evt->Set("EventNumber", (int)Get<UInt_t>("EventNumber"));

	
	if(fEvt->Debug()) cout << fEvt->RunNumber() << " is run number  and " << fEvt->EventNumber() << "is enum" << endl;
	
	evt->Set("larerror", Get<unsigned int>("larError")); 
  evt->Set("tileerror", Get<unsigned int>("tileError")); 
  evt->Set("coreflags", Get<unsigned int>("coreFlags")); 
  //evt->Set("isMC",(bool)(Get<UInt_t>("RunNumber")<152166));
	evt->Set("averageIntPerXing", Get<float>("averageIntPerXing"));
  evt->Set("Eventshape_rhoKt4LC", Get<float>("Eventshape_rhoKt4LC"));
  evt->Set("Eventshape_rhoKt4EM", Get<float>("Eventshape_rhoKt4EM"));
	evt->Set("LBN", (int)Get<UInt_t>("lbn"));
	evt->Set("BunchCrossingID", (int)Get<UInt_t>("bcid"));

  static const BranchKey bmc_channel_number("mc_channel_number");
	static const MomKey misMC("isMC"), mChannelNumber("ChannelNumber");
	if(BranchNames().find(bmc_channel_number)!=BranchNames().end()) {
		evt->Set(mChannelNumber,(int)Get<UInt_t>(bmc_channel_number));
		evt->Set(misMC,true);	
	}
	else evt->Set(misMC,false);	


  static const MomKey mEventWeight("EventWeight");
  static const BranchKey bmcevt_weight("mcevt_weight");
  if(fEvt->isMC()){
    if(Get<vector<vector<double> > >(bmcevt_weight).size()) {
      if(Get<vector<vector<double> > >(bmcevt_weight).at(0).size()) { //This stuff necessary due to a bug it seems
        fEvt->Set(mEventWeight,(float)Get<vector<vector<double> > >(bmcevt_weight).at(0).at(0));
      }
    }
    else Abort("ERROR: mcevt_weight vector is empty!");
  }
  else{
    evt->Set("EventWeight", 1.0); // default value
  }

  if(evt->Debug()){
    cout << "showing at the event quantities" << endl;
    evt->Show();
  }


  if(fEvt->Debug()) cout << " basic event quantities done" << endl;
////////////////////////////////////////////////////////////
//// Do GRL Parsing, and HFOR/PU 
////////////////////////////////////////////////////////////

  if(doPRW){
    SetupPRW();
  }

  if(evt->Debug()){
    cout << "showing at the setup prw" << endl;
    evt->Show();
  }

  static const MomKey mgrl("grl"), LBN("LBN"), RunNumber("RunNumber");

  if (!(fEvt->isMC())){
    if (evt->Debug()){
      cout << " checking GRL " << endl;
      cout << LBN << " " << RunNumber << endl;
      if (myGRL){
	cout << "yes " << endl;
      }
      else{
	cout << "no" << endl;
      }
    }
  }

  if(fEvt->isMC() && !doMCGRL){ 
    fEvt->Set(mgrl,true);
  } else if(myGRL){
    fEvt->Set(mgrl,myGRL->HasRunLumiBlock(fEvt->Int(RunNumber),fEvt->Int(LBN)));
  } 


  if(evt->Debug()) cout << "set grl" << endl;
   

  if(fEvt->isMC()) {
    if(!doITKFixes && !KillPileupFluctuations())
      return kFALSE;
    if(!DoHFOverlap()) 
      return kFALSE;
  }

  if(evt->Debug()){
    cout << "showing at the hfor" << endl;
    evt->Show();
  }

  if(fEvt->Debug()) cout << " hfor done " << endl;
////////////////////////////////////////////////////////////
//// Copy over all the objects
////////////////////////////////////////////////////////////


  bool fail = true;

  // NOTE: I put the truth copy in the head for track-truth association.
  fail = fail && CopyTruthFull();

  if(evt->Debug()){
    cout << "showing at the end of truthfull" << endl;
    evt->Show();
  }

  fail = fail && AddTruthParentChild();

  fail = fail && ParseJetStrings();
  fail = fail && CopyTracks();
  fail = fail && CopyVertices();
  fail = fail && CopyClusters();

  if(evt->Debug()){
    cout << "showing at the end of tracks" << endl;
    evt->Show();
    cout << "HaHaHa" << endl;
  }

  fail = fail && CopyMuons();
  
  if(evt->Debug()){
    cout << "end of CopyMuons" << endl;
  }
  
  fail = fail && CopyMuonTriggers();

  if(evt->Debug()) cout << "end of CopyMuonTriggers" << endl;
  
  fail = fail && CopyElectronTriggers();
  
  if(evt->Debug()) cout << "end of copyElectronTriggers" << endl;

  if(evt->Debug()){
    cout << "showing at the of meuons" << endl;
    evt->Show();
  }

  if(evt->Debug()){
    cout << "showing at the end of truths" << endl;
    evt->Show();
  }
  fail = fail && CopyJets();
  fail = fail && AddBtags();
  fail = fail && CopyPhotons();
  fail = fail && CopyElectrons();
  fail = fail && AddPRWWeights();
  fail = fail && CopyJetTriggers();
  fail = fail && CopyMET();

  if(evt->Debug()){
    cout << "showing at the end of copies" << endl;
    evt->Show();
  }

  evt->SortPtAll();
  
  if(fEvt->Debug()) cout << "abandon all hope, ye who enter here" << endl; 
  fail = fail && DoCustomSetup();
  
  if(fEvt->Debug()) cout << "set up full event, status is " << fail << endl;
  
  if(evt->Debug()){
    cout << "showing at the end" << endl;
    evt->Show();
  }

	return fail;
}

///=========================================
/// CopyMET
///=========================================
Bool_t EventBuilder_jetmet2012::CopyMET(){

  if(doMETRefFinal){
	  if(doSMWZ || doCOMMON){
     m_util->reset();
     m_util->setJetParameters(&jet_AntiKt4LCTopo_recalibrated_pt, 
                              GetP<vector<float> >("jet_AntiKt4LCTopo_eta"), 
                              GetP<vector<float> >("jet_AntiKt4LCTopo_phi"), 
                              GetP<vector<float> >("jet_AntiKt4LCTopo_E"), 
                              GetP<vector<vector<float> > >("jet_AntiKt4LCTopo_MET_wpx"),
                              GetP<vector<vector<float> > >("jet_AntiKt4LCTopo_MET_wpy"),
                              GetP<vector<vector<float> > >("jet_AntiKt4LCTopo_MET_wet"),
                              GetP<vector<vector<unsigned int> > >("jet_AntiKt4LCTopo_MET_statusWord"));
     
     m_util->setMETTerm(METUtil::RefEle, Get<float>("MET_RefEle_et")*cos(Get<float>("MET_RefEle_phi")),Get<float>("MET_RefEle_et")*sin(Get<float>("MET_RefEle_phi")),Get<float>("MET_RefEle_sumet"));
     m_util->setMETTerm(METUtil::RefTau,Get<float>("MET_RefTau_et")*cos(Get<float>("MET_RefTau_phi")),Get<float>("MET_RefTau_et")*sin(Get<float>("MET_RefTau_phi")),Get<float>("MET_RefTau_sumet"));
     m_util->setMETTerm(METUtil::SoftTerms, Get<float>("MET_CellOut_Eflow_et")*cos(Get<float>("MET_CellOut_Eflow_phi")),Get<float>("MET_CellOut_Eflow_et")*sin(Get<float>("MET_CellOut_Eflow_phi")),Get<float>("MET_CellOut_Eflow_sumet"));
     m_util->setElectronParameters(GetP<vector<float> >("el_pt"),GetP<vector<float> >("el_eta"),GetP<vector<float> >("el_phi"),GetP<vector<vector<float> > >("el_MET_wet"),GetP<vector<vector<float> > >("el_MET_wpx"),GetP<vector<vector<float> > >("el_MET_wpy"),GetP<vector<vector<unsigned int> > >("el_MET_statusWord"));
     m_util->setMuonParameters(GetP<vector<float> >("mu_pt"),GetP<vector<float> >("mu_eta"),GetP<vector<float> >("mu_phi"),GetP<vector<vector<float> > >("mu_MET_wet"),GetP<vector<vector<float> > >("mu_MET_wpx"),GetP<vector<vector<float> > >("mu_MET_wpy"),GetP<vector<vector<unsigned int> > >("mu_MET_statusWord"));
     m_util->setExtraMuonParameters(GetP<vector<float> >("mu_ms_qoverp"),GetP<vector<float> >("mu_ms_theta"),GetP<vector<float> >("mu_ms_phi"),GetP<vector<float> >("mu_charge"));
     
	  //These are from the twiki, but HSG3 does deals with the leptons slightly differently; we do it the way Higgs way.
	  //m_util->setMETTerm(METUtil::RefGamma, Get<float>("MET_RefGamma_et")*cos(Get<float>("MET_RefGamma_phi")),Get<float>("MET_RefGamma_et")*sin(Get<float>("MET_RefGamma_phi")),Get<float>("MET_RefGamma_sumet"));
	  //m_util->setMETTerm(METUtil::MuonTotal, Get<float>("MET_MuonBoy_et")*cos(Get<float>("MET_MuonBoy_phi")),Get<float>("MET_MuonBoy_et")*sin(Get<float>("MET_MuonBoy_phi")),Get<float>("MET_MuonBoy_sumet"));
     
	  //to get METRefFinal
     METUtility::METObject refFinal;
  switch(whichsyste) {
    /// Soft Terms scale
    case JetMETSyst::METSCALESTUP: 
      refFinal = m_util->getMissingET(METUtil::RefFinal,METUtil::ScaleSoftTermsUp); 
      break;
    case JetMETSyst::METSCALESTDOWN: 
      refFinal = m_util->getMissingET(METUtil::RefFinal,METUtil::ScaleSoftTermsDown); 
      break;
    /// Soft Terms resolution
    case JetMETSyst::METRESOST: 
      refFinal = m_util->getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsUp); 
      break;
    default:
      refFinal = m_util->getMissingET(METUtil::RefFinal);
      break;
  }



     double refFinal_et = refFinal.et()/1000.;
     double refFinal_phi = refFinal.phi();
     
     fEvt->Set("NewMET",refFinal_et);	  
     fEvt->Set("NewMETphi",refFinal_phi);
     fEvt->Set("NewMETx",refFinal_et*cos(refFinal_phi));	  
     fEvt->Set("NewMETy",refFinal_et*sin(refFinal_phi));
	  }
    float MET_et = Get<float>("MET_RefFinal_et") / 1000.;
    float MET_phi = Get<float>("MET_RefFinal_phi");

    fEvt->Set("MET_RefFinal_et", MET_et);
    fEvt->Set("MET_RefFinal_phi", MET_phi);

    fEvt->Set("MET_RefFinal_x", MET_et*cos(MET_phi));
    fEvt->Set("MET_RefFinal_y", MET_et*sin(MET_phi));
    
    if(Rename){
      fEvt->Set("Metx",      Get<Float_t>("MET_RefFinal_etx")/1000.);
      fEvt->Set("Mety",      Get<Float_t>("MET_RefFinal_ety")/1000.);
      fEvt->Set("Met_et",    Get<Float_t>("MET_RefFinal_et")/1000.);
      fEvt->Set("Met_phi",   Get<Float_t>("MET_RefFinal_phi"));
      fEvt->Set("Met_sumet", Get<Float_t>("MET_RefFinal_sumet")/1000.);
    }  

    //cout << "refFinal_et " << refFinal_et << " old ref = " << MET_et << endl;
  }

  return kTRUE;
}


///=========================================
/// CopyTruthFull
///=========================================
Bool_t EventBuilder_jetmet2012::CopyTruthFull()
{
  
  static const MomKey mtruths("truths");


  static const BranchKey bmc_pt("mc_pt"), bmc_eta("mc_eta"), bmc_phi("mc_phi"), bmc_m("mc_m"),
    bmc_pdgId("mc_pdgId"), bmc_status("mc_status"), bmc_charge("mc_charge"), bmc_barcode("mc_barcode"), bmc_vx_barcode("mc_vx_barcode");

  static const MomKey index("index"), type("type"), pdgId("pdgId"), charge("charge"), barcode("barcode"), mpt("pt"), status("status"), vx_barcode("vx_barcode");


  barcode_link = 0;

  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012: DEBUG CopyTruthFull()  and doTruth = " << doTruth << endl;  
  
  /// If scheduled, set the truth particles
  if(fEvt->isMC() && doTruth) {
    
    fEvt->AddVec(mtruths);
    
    if(doBarcode){
      fEvt->AddVec("barcode_database");
      barcode_link = new MomentObj();
      fEvt->Add("barcode_database", barcode_link);
    }
    
    for (int i = 0; i< Get<Int_t>("mc_n"); i++) {

      Particle * mc = new Particle();
      float pt  = Get<vector<float> >(bmc_pt) .at(i);
      float eta = Get<vector<float> >(bmc_eta).at(i);
      float phi = Get<vector<float> >(bmc_phi).at(i);
      float m   = Get<vector<float> >(bmc_m)  .at(i);
      mc->Set("pt", pt);
      mc->p.SetPtEtaPhiM(pt/1000., eta, phi, m/1000.);
      mc->Set(index, i );      
      //mc->Set(mpt  , pt);

      mc->Set(type,      Get<vector<int> >  (bmc_pdgId)     .at(i)); 
      mc->Set(pdgId,     Get<vector<int> >  (bmc_pdgId)     .at(i)); 
      mc->Set(status,    Get<vector<int> >  (bmc_status)    .at(i)); 
      mc->Set(charge,    Get<vector<float> >(bmc_charge)    .at(i));
      mc->Set(barcode,   Get<vector<int> >  (bmc_barcode)   .at(i));
      mc->Set(vx_barcode,Get<vector<int> >  (bmc_vx_barcode).at(i));

      if(Rename){
        mc->Set("mc_type",    Get<vector<int> >(bmc_pdgId)    .at(i));
      	mc->Set("mc_pdgId",   Get<vector<int> >(bmc_pdgId)    .at(i));
	      mc->Set("mc_status",  Get<vector<int> >(bmc_status)   .at(i));
	      mc->Set("mc_charge",  Get<vector<float> >(bmc_charge) .at(i));
	      mc->Set("mc_barcode", Get<vector<int> >(bmc_barcode)  .at(i));
      }

      if(doBarcode){
        //MomentObj* barcode_link = fEvt->Obj<MomentObj>("barcode_database");
	barcode_link->AddVec(TString::Itoa(mc->Int("mc_barcode"), 10));
	barcode_link->Add(TString::Itoa(mc->Int("mc_barcode"), 10), mc);
      }

      fEvt->Add(mtruths,mc);
      if(fEvt->Debug()){
        //cout << "Going to print particle " << i << endl;
        //mc->Show();
      }
    } // end for loop over particles
  }
  
  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::CopyTruthFull(): DEBUG Scheduled Truth Particles" << endl;
  return kTRUE;
}

///=========================================
/// AddTruthParentChild
///=========================================
Bool_t EventBuilder_jetmet2012::AddTruthParentChild()
{
  TString prefix;
  if(Rename) prefix = "mc_";
  else       prefix = "";

  static const MomKey children(prefix + "children"), parents(prefix + "parents");

  TString child_branch, parent_branch;
  if(!doCOMMON){
    child_branch  = "mc_children";
    parent_branch = "mc_parents";
  }
  else{
    child_branch  = "mc_child_index";
    parent_branch = "mc_parent_index";
  }

  static const BranchKey bmc_children(child_branch), bmc_parents(parent_branch);

  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::AddTruthParentChild(): DEBUG Starting" << endl;

  if(fEvt->isMC() && doTruthParentChild)
  {
    for(int iTr = 0; iTr < fEvt->truths(); iTr++)
    {
      if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::AddTruthParentChild(): DEBUG On particle " << iTr << " of " << fEvt->truths() << endl;

      fEvt->truth(iTr).AddVec(children);
      fEvt->truth(iTr).AddVec(parents, true);
      if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::AddTruthParentChild(): DEBUG Has children " << Get<vector<vector<int> > >(bmc_children).at(iTr).size()  << endl;

      int nChildren = Get<vector<vector<int> > >(bmc_children).at(iTr).size();
      for(int iChild = 0; iChild < nChildren; iChild++)
      {
        int index = Get<vector<vector<int> > >(bmc_children).at(iTr).at(iChild);
        if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::AddTruthParentChild(): DEBUG Adding child " << index  << endl;
           
        if(!doCOMMON){     
          if((index > -1) && (index < fEvt->truths()) && (index != iTr))
            fEvt->truth(iTr).Add(children,&(fEvt->truth(index)));
          else{
            if(fEvt->Debug()) cout << "Ops ... illegal truth index" << endl;
          }
        }
        else{
          fEvt->truth(iTr).Add(children,&(fEvt->truth(index)));
        }
      }

      if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::AddTruthParentChild(): DEBUG Has parents " << Get<vector<vector<int> > >(bmc_parents).at(iTr).size()  << endl;

      int nParents = Get<vector<vector<int> > >(bmc_parents).at(iTr).size();
      for(int iParent = 0; iParent < nParents; iParent++)
      {
        int index = Get<vector<vector<int> > >(bmc_parents).at(iTr).at(iParent);
        if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::AddTruthParentChild(): DEBUG Adding parent " << index << endl;
        
        if(!doCOMMON){
          if((index > -1) && (index < fEvt->truths()) && (index != iTr))
            fEvt->truth(iTr).Add(parents,&(fEvt->truth(index)));
          else{
            if(fEvt->Debug()) cout << "Ops ... illegal truth index" << endl;
           }
        }
        else{
          fEvt->truth(iTr).Add(parents,&(fEvt->truth(index)));
        }
      }
      
      if(fEvt->Debug()){
        cout << parents.Data() << " Has member " << fEvt->truth(iTr).Objs(parents) << endl;
        cout << children.Data() << " Has memebr " << fEvt->truth(iTr).Objs(children) << endl;
      }
    }
  } 

  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::AddTruthParentChild(): DEBUG Ending" << endl;
  return kTRUE;

}


Bool_t EventBuilder_jetmet2012::CopyClusters(){
  if(doLCCluster){
    static const MomKey clustersKey("clustersLCTopo");
    fEvt->AddVec(clustersKey);
    BranchKey bn, bpt, beta, bphi;
    if(doSMWZ ){
      bn   = BranchKey("cl_n");
      bpt  = BranchKey("cl_pt");
      beta = BranchKey("cl_eta");
      bphi = BranchKey("cl_phi");
    } else if (doCOMMON){
      bn   = BranchKey("cl_lc_n");
      bpt  = BranchKey("cl_lc_pt");
      beta = BranchKey("cl_lc_eta");
      bphi = BranchKey("cl_lc_phi");
    } else  {
      bn   = BranchKey("cl_had_n");
      bpt  = BranchKey("cl_had_pt");
      beta = BranchKey("cl_had_eta");
      bphi = BranchKey("cl_had_phi");
    }
    for(int iCl = 0; iCl < Get<int>(bn); iCl++){
      Particle* cl = new Particle();
      float pt = Get<vector<float> >(bpt).at(iCl);
      float eta = Get<vector<float> >(beta).at(iCl);
      float phi = Get<vector<float> >(bphi).at(iCl);
      cl->p.SetPtEtaPhiE(pt/1000.,eta,phi,pt/1000.*cosh(eta));
      if(doCOMMON){
        cl->Set("cl_centerlambda", Get<vector<float> >("cl_centerlambda").at(iCl));
        cl->Set("cl_firstEdens",   Get<vector<float> >("cl_firstEdens")  .at(iCl));
      }
      fEvt->Add(clustersKey,cl);
    }
  }
  if(doEMCluster){
    static const MomKey clustersKey("clustersEM");
    fEvt->AddVec(clustersKey);
    BranchKey bn, bpt, beta, bphi;
    if(doSMWZ ){
      bn   = BranchKey("cl_n");
      bpt  = BranchKey("cl_pt");
      beta = BranchKey("cl_eta");
      bphi = BranchKey("cl_phi");
    } else if (doCOMMON){
      bn   = BranchKey("cl_em_n");
      bpt  = BranchKey("cl_em_pt");
      beta = BranchKey("cl_em_eta");
      bphi = BranchKey("cl_em_phi");
    } else  {
      bn   = BranchKey("cl_had_n");
      bpt  = BranchKey("cl_had_pt");
      beta = BranchKey("cl_had_eta");
      bphi = BranchKey("cl_had_phi");
    }
    for(int iCl = 0; iCl < Get<int>(bn); iCl++){
      Particle* cl = new Particle();
      float pt = Get<vector<float> >(bpt).at(iCl);
      float eta = Get<vector<float> >(beta).at(iCl);
      float phi = Get<vector<float> >(bphi).at(iCl);
      cl->p.SetPtEtaPhiE(pt/1000.,eta,phi,pt/1000.*cosh(eta));
      if(doCOMMON){
        cl->Set("cl_centerlambda", Get<vector<float> >("cl_centerlambda").at(iCl));
        cl->Set("cl_firstEdens",   Get<vector<float> >("cl_firstEdens")  .at(iCl));
      }
      fEvt->Add(clustersKey,cl);
    }
  }
  if(fEvt->Debug()) cout << "set up clusters" << endl;
  return kTRUE;
}


Bool_t EventBuilder_jetmet2012::CopyTracks(){

  static const MomKey mtracks("tracks"),
    nPixHits("nPixHits"), nPixHoles("nPixHoles"), nPixelDeadSensors("nPixelDeadSensors"), nSCTHits("nSCTHits"),
    nSCTHoles("nSCTHoles"), nSCTDeadSensors("nSCTDeadSensors"), chi2("chi2"), ndof("ndof"), nTRTOutliers("nTRTOutliers"),
    nTRTHits("nTRTHits"), qoverp("qoverp"), theta("theta"), d0_wrtPV("d0_wrtPV"), z0_wrtPV("z0_wrtPV"), z0("z0"), d0_wrtBL("d0_wrtBL"), d0_wrtBS("d0_wrtBS"),
    eta_atCalo("eta_atCalo"), phi_atCalo("phi_atCalo"), mc_barcode("mc_barcode"), index("index");



  static const BranchKey btrk_n("trk_n"), btrk_pt("trk_pt"), btrk_eta("trk_eta"), btrk_nPixHits("trk_nPixHits"),
    btrk_nPixHoles("trk_nPixHoles"), btrk_nPixelDeadSensors("trk_nPixelDeadSensors"), btrk_nSCTHits("trk_nSCTHits"), 
    btrk_nSCTHoles("trk_nSCTHoles"), btrk_nSCTDeadSensors("trk_nSCTDeadSensors"), btrk_d0_wrtPV("trk_d0_wrtPV"),
    btrk_z0_wrtPV("trk_z0_wrtPV"), btrk_z0("trk_z0"), btrk_d0_wrtBL("trk_d0_wrtBL"), btrk_d0_wrtBS("trk_d0_wrtBS"), btrk_chi2("trk_chi2"), btrk_ndof("trk_ndof"),
    btrk_nTRTOutliers("trk_nTRTOutliers"), btrk_nTRTHits("trk_nTRTHits"),
    btrk_eta_atCalo("trk_eta_atCalo"), btrk_phi_atCalo("trk_phi_atCalo"), btrk_mc_barcode("trk_mc_barcode");

  BranchKey fix, fix_phi;

  if((doSMWZ || doCOMMON) && !doITKFixes){
    fix_phi = BranchKey("_wrtPV");
    fix     = BranchKey("_wrtPV");
  }else if((doSMWZ || doCOMMON) && doITKFixes){
    fix_phi = BranchKey("_wrtPV");
    fix     = BranchKey("");
  } else {
    fix     = BranchKey("");
    fix_phi = BranchKey("");
  }

  const BranchKey btrk_phi("trk_phi"+fix_phi), btrk_qoverp("trk_qoverp"+fix), btrk_theta("trk_theta"+fix);

  if(doTrack){

    fEvt->AddVec(mtracks);
    
    if(doTrkJetBTag) fEvt->AddVec("tracks_trkZjet", true);
    
    for(int iTr = 0; iTr < Get<Int_t>(btrk_n); iTr++){
      Particle* tr = new Particle();
      
      if(!doTrkJetBTag){
        float pt =  Get<vector<float> >(btrk_pt) .at(iTr);
        float eta = Get<vector<float> >(btrk_eta).at(iTr);
        float phi = Get<vector<float> >(btrk_phi).at(iTr);
        tr->p.SetPtEtaPhiM(pt/1000.,eta,phi,0);

        tr->Set(index,             iTr); 
        tr->Set(qoverp,            Get<vector<float> >(btrk_qoverp)           .at(iTr));
        tr->Set(nPixHits,          Get<vector<int> >  (btrk_nPixHits)         .at(iTr));
        tr->Set(nPixHoles,         Get<vector<int> >  (btrk_nPixHoles)        .at(iTr));
        tr->Set(nPixelDeadSensors, Get<vector<int> >  (btrk_nPixelDeadSensors).at(iTr));
        tr->Set(nSCTHits,          Get<vector<int> >  (btrk_nSCTHits)         .at(iTr));
        tr->Set(nSCTHoles,         Get<vector<int> >  (btrk_nSCTHoles)        .at(iTr));
        tr->Set(nSCTDeadSensors,   Get<vector<int> >  (btrk_nSCTDeadSensors)  .at(iTr));
        tr->Set(d0_wrtPV,          Get<vector<Float_t> >(btrk_d0_wrtPV)         .at(iTr));
        tr->Set(d0_wrtBL,          Get<vector<float> >(btrk_d0_wrtBL)         .at(iTr));
        if(!doSMWZ && !doCOMMON)
          tr->Set(d0_wrtBS,          Get<vector<float> >(btrk_d0_wrtBS)         .at(iTr));
        tr->Set(z0_wrtPV,          Get<vector<float> >(btrk_z0_wrtPV)         .at(iTr));
        if(!doSMWZ && !doCOMMON)
          tr->Set(z0,                Get<vector<float> >(btrk_z0)               .at(iTr));
        tr->Set(chi2,              Get<vector<float> >(btrk_chi2)             .at(iTr));
        tr->Set(theta,             Get<vector<float> >(btrk_theta)            .at(iTr));
        tr->Set(ndof,              Get<vector<int> >  (btrk_ndof)             .at(iTr));
        //tr->Set(nTRTOutliers,      Get<vector<int> >  (btrk_nTRTOutliers)     .at(iTr));
        tr->Set(nTRTHits,          Get<vector<int> >  (btrk_nTRTHits)         .at(iTr));
        if(!doSMWZ)
        tr->Set(eta_atCalo,        Get<vector<float> > (btrk_eta_atCalo)      .at(iTr));
        if(!doSMWZ)
        tr->Set(phi_atCalo,        Get<vector<float> > (btrk_phi_atCalo)      .at(iTr));
        if(!doSMWZ && fEvt->isMC())
        tr->Set(mc_barcode,        Get<vector<int> > (btrk_mc_barcode)        .at(iTr));
        if(!doSMWZ){
        tr->Set("err_d0_wrtPV",      TMath::Sqrt(Get<vector<float> > ("trk_cov_d0_wrtPV")    .at(iTr)));
        tr->Set("err_z0_wrtPV",      TMath::Sqrt(Get<vector<float> > ("trk_cov_z0_wrtPV")    .at(iTr)));
        }
      }
      else{
        float pt  = Get<vector<float> >("trk_pt") .at(iTr);
        float eta = Get<vector<float> >("trk_eta").at(iTr);
        //float qoverp = Get<vector<float> >("trk_qoverp_wrtPV").at(iTr);
        //float theta = Get<vector<float> >("trk_theta_wrtPV").at(iTr);

        //float pt = TMath::Sin( theta ) / fabs( qoverp );
        //float eta = -1.0 * log( TMath::Tan( theta / 2.0 )  );
        float phi = Get<vector<float> >("trk_phi_wrtPV").at(iTr);
        float Mpi = 139.57018; // MeV
        tr->p.SetPtEtaPhiM(pt/1000.,eta,phi,Mpi/1000.);
      
        tr->Set("d0",              Get<vector<float> >("trk_d0_wrtPV")    .at(iTr));
        if(!doCOMMON)
          tr->Set("d0_wrtBS",        Get<vector<float> >("trk_d0_wrtBS")    .at(iTr));
        else
          tr->Set("d0_wrtBS", tr->Float("d0"));  // Because I don't want to change all following labels ... 
        tr->Set("z0",              Get<vector<float> >("trk_z0_wrtPV")    .at(iTr));
        tr->Set("z0SinTheta",      tr->Float("z0")*TMath::Sin(tr->p.Theta()));
        tr->Set("err_d0",          TMath::Sqrt(Get<vector<float> > ("trk_cov_d0_wrtPV")    .at(iTr)));
        tr->Set("err_z0",          TMath::Sqrt(Get<vector<float> > ("trk_cov_z0_wrtPV")    .at(iTr)));
        tr->Set("err_phi",         TMath::Sqrt(Get<vector<float> > ("trk_cov_phi_wrtPV")   .at(iTr)));
        tr->Set("err_theta",       TMath::Sqrt(Get<vector<float> > ("trk_cov_theta_wrtPV") .at(iTr)));
        tr->Set("err_qoverp",      TMath::Sqrt(Get<vector<float> > ("trk_cov_qoverp_wrtPV").at(iTr)));
     
        //NOTE: ignoring COV(z0, theta) term... since we don't have it
        double err_z0SinTheta = TMath::Sqrt( pow(TMath::Sin(tr->p.Theta()) * tr->Float("err_z0"), 2) + pow(TMath::Cos(tr->p.Theta()) * tr->Float("z0") * tr->Float("err_theta") , 2) );
        tr->Set("err_z0SinTheta", err_z0SinTheta);
        
        tr->Set("IPtSig", (tr->Float("d0"))/(tr->Float("err_d0")));
        tr->Set("IPzSig", (tr->Float("z0SinTheta"))/(tr->Float("err_z0SinTheta")));
 
        tr->Set("nPixHits",          Get<vector<int> >  ("trk_nPixHits")             .at(iTr));
        tr->Set("nBLHits",           Get<vector<int> >  ("trk_nBLHits")              .at(iTr));
        tr->Set("nSCTHits",          Get<vector<int> >  ("trk_nSCTHits")             .at(iTr));
        tr->Set("nPixHoles",         Get<vector<int> >  ("trk_nPixHoles")            .at(iTr));
        tr->Set("nSCTHoles",         Get<vector<int> >  ("trk_nSCTHoles")            .at(iTr));
        tr->Set("nPixDeadSensors",   Get<vector<int> >  ("trk_nPixelDeadSensors")    .at(iTr));
        tr->Set("nSCTDeadSensors",   Get<vector<int> >  ("trk_nSCTDeadSensors")      .at(iTr));

        tr->Set("chi2",              Get<vector<float> >("trk_chi2")                 .at(iTr));
        tr->Set("ndof",              Get<vector<int> >  ("trk_ndof")                 .at(iTr));

        if(!doCOMMON){
          tr->Set("nBLSharedHits",     Get<vector<int> >  ("trk_nBLSharedHits")        .at(iTr));
          tr->Set("nPixSharedHits",    Get<vector<int> >  ("trk_nPixSharedHits")       .at(iTr));
          tr->Set("nSCTSharedHits",    Get<vector<int> >  ("trk_nSCTSharedHits")       .at(iTr));
        } 

        tr->Set("mc_barcode",        Get<vector<int> >  ("trk_mc_barcode")           .at(iTr));
        tr->Set("mc_probability",    Get<vector<float> >("trk_mc_probability")       .at(iTr));
        if(doCOMMON){
          tr->Set("mc_index",          Get<vector<int> >("trk_mc_index")               .at(iTr));
          int mcindex = tr->Int("mc_index");

          // track-truth association (before sorting so that mc_index is useful)
          if(fEvt->Debug()) cout << "In CopyTracks: track-truth association for track " << iTr << " and truth " << mcindex << endl;
          if(mcindex != -1){
            tr->AddVec("MatchedTruth_D3PD", true);
            tr->Add("MatchedTruth_D3PD", &fEvt->truth(mcindex));
          }
        }
      }

      fEvt->Add(mtracks,tr);
      
      // TrackZ collection
      if(doTrkJetBTag){
        if( tr->p.Pt() > 0.5                                &&
	          fabs(tr->p.Eta()) < 2.5                         &&
	          tr->Int("nPixHits") > 0                         &&
            (tr->Int("nPixHits") + tr->Int("nBLHits")) > 0  &&
	          tr->Int("nSCTHits") > 5                         &&
         	  (tr->Int("nPixHits") + tr->Int("nSCTHits")) > 6 &&
	          fabs(tr->Float("d0")) < 1                       &&
	          //fabs(tr->Float("z0")*sin(tr->Float("theta"))) < 1.5 &&
	          fabs(tr->Float("z0")) < 1.5                     &&
	          //fabs(tr->Float("z0_wrtBS")) < 200 &&
	          tr->Float("chi2") /  tr->Float("ndof") < 3
          )
          	fEvt->Add("tracks_trkZjet",tr);
      }
      
    }
  }
  if(fEvt->Debug()) cout << "set up tracks" << endl;
  return kTRUE;
}

Bool_t EventBuilder_jetmet2012::CopyVertices()
{

  static const MomKey mvtxs("vtxs"), ntrk("ntrk"), mNPV("NPV"), sumPt("sumPt"), chisq("chisq"), ndof("ndof");

  static const BranchKey bvxp_nTracks("vxp_nTracks"), bvxp_x("vxp_x"), bvxp_y("vxp_y"), 
                         bvxp_z("vxp_z"), bvxp_chisq("vxp_chi2"), bvxp_ndof("vxp_ndof"), bvxp_sumpt("vxp_sumPt"), bvxp_trk_index("vxp_trk_index");

  static const MomKey mvtxsTruth("vtxsTruth"), mNPVTruth("NPVTruth");
  static const BranchKey bmcVx_x("mcVx_x"), bmcVx_y("mcVx_y"), bmcVx_z("mcVx_z");
  static const BranchKey bmcVxSMWZ_x("mc_vx_x"), bmcVxSMWZ_y("mc_vx_y"), bmcVxSMWZ_z("mc_vx_z");

  if(doVertex){
    if(fEvt->isMC() ) {
      fEvt->AddVec(mvtxsTruth);

      unsigned int npvTruth = doSMWZ ?Get<vector<float> >(bmcVxSMWZ_z).size():Get<vector<float> >(bmcVx_z).size();
      fEvt->Set(mNPVTruth, npvTruth);
      
      for(unsigned int pvi = 0; pvi < npvTruth; pvi++) {
	    Point* vtx = new Point();
        if(!doSMWZ){
	        vtx->x.SetXYZ(Get<vector<float> >(bmcVx_x).at(pvi),Get<vector<float> >(bmcVx_y).at(pvi),Get<vector<float> >(bmcVx_z).at(pvi));
        }else{
	        vtx->x.SetXYZ(Get<vector<float> >(bmcVxSMWZ_x).at(pvi),Get<vector<float> >(bmcVxSMWZ_y).at(pvi),Get<vector<float> >(bmcVxSMWZ_z).at(pvi));
        }
	    fEvt->Add(mvtxsTruth, vtx);
      }
    }

    fEvt->AddVec(mvtxs);

    //Vertices
    int NPV=0;
    for(unsigned pvi=0;pvi<Get<vector<int> >(bvxp_nTracks).size();++pvi) {
      if(Get<vector<int> >(bvxp_nTracks).at(pvi) >= 2){
        ++NPV;
        Point* vtx = new Point();
        vtx->x.SetXYZ(Get<vector<float> >(bvxp_x).at(pvi),Get<vector<float> >(bvxp_y).at(pvi),Get<vector<float> >(bvxp_z).at(pvi));
        vtx->Set(ntrk,Get<vector<int> >(bvxp_nTracks).at(pvi));
        vtx->Set(sumPt, Get<vector<float> >(bvxp_sumpt).at(pvi)/1000.);
        vtx->Set(chisq, Get<vector<float> >(bvxp_chisq).at(pvi));
        vtx->Set(ndof, Get<vector<int> >(bvxp_ndof).at(pvi));
        vtx->AddVec("vtxTracks", true);
        vector<int> trk_index = Get<vector<vector<int> > >(bvxp_trk_index).at(pvi);
        for(int iTrk=0; iTrk<trk_index.size(); ++iTrk){
            vtx->Add("vtxTracks", &(fEvt->track(trk_index[iTrk])), true);
        }

        fEvt->Add(mvtxs,vtx);
      }
    }
    fEvt->Set(mNPV,NPV);

  }
  if(fEvt->Debug()) cout << "set up vtx" << endl;
  return kTRUE;
}

Bool_t EventBuilder_jetmet2012::CopyPhotons(){
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


void EventBuilder_jetmet2012::WriteAuxTrees(TDirectory* outfile){
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
Bool_t EventBuilder_jetmet2012::AddBtags() 
{

  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012: DEBUG In AddBtags() " << endl;   
  bool fail = true;
  /// Copy all jet types
  for(Int_t i = 0; i<fEvt->Cfg()->Objs(mBTAGVEC); ++i) {
    TString type = ((TObjString*)fEvt->Cfg()->Obj(mBTAGVEC,i))->GetString();
    if (type.Contains("Truth")) continue; 
    if(fEvt->Debug()) cout << "EventBuilder_jetmet2012::AddBtags(): DEBUG Copying btag type " << type << " which is the " << i << "th jet collection" << endl;
    if(fEvt->Debug()) cout << "EventBuilder_jetmet2012::AddBtagss(): DEBUG type in vec is " << ((TObjString*) fEvt->Cfg()->Obj(mBTAGVEC,i))->GetString() << endl;
    fail = fail && AddBtag(type,"");
    //delete type;
  }
  
  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::AddBtags(): DEBUG Scheduled All Jet B-tags" << endl;
  return fail;
}

///=========================================
/// AddBtag
///=========================================
Bool_t EventBuilder_jetmet2012::AddBtag(TString jetType, TString jetClass) 
{

  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012: DEBUG In AddBtag() for jetType = " << jetType << "; jetClass = " << jetClass << endl;  
  
  // Handle the case where this is called without an argument
  if (jetType == "") { 
    Abort("EventBuilder_jetmet2012::AddBtag(): ERROR No jet type selected. Aborting!");
  }
  
  // Now try to figure out what type of jet this is
  if (jetClass == "") {
    if (fEvt->Debug()) {
      cout << "EventBuilder_jetmet2012::AddBtag(): INFO The jetClass variable is empty. ";
      cout << "Assuming that the class (truth, track, topo, etc) can be determined from the name." << endl;
    } // END IF
    
    // Parse the jetType: truth, track, topo
    if (jetType.Contains("Truth")) {
     
      cout << "EventBuilder_jetmet2012::AddBtag(): ERROR Cannot add B-tagging for Truth jet collections!" << endl;
      return kFALSE;
    
    }
    else if (jetType.Contains("Track"))  jetClass = "track";
    else if (jetType.Contains("LCTopo")) jetClass = "topo";
    else if (jetType.Contains("TopoEM")) jetClass = "topo";
    else if (jetType.Contains("Tower"))  {
      cout << "EventBuilder_jetmet2012::AddBtag(): ERROR Tower jets are not currently supported." << endl;
      return kFALSE;
    } // END IF
    
  } // END IF JET CLASS 
  
  // Otherwise, try to figure out what class of jet this really is
  else if (jetClass != "topo" && jetClass != "track") {
    cout << "EventBuilder_jetmet2012::AddBtag(): ERROR Currently, the only supported classes of jets are: truth, track, topo." << endl;
    return kFALSE;
  }
  

  // Collect jet flavor information
  for(unsigned int iJet = 0; iJet < Get<vector<float> >("jet_" + jetType + "_phi").size(); iJet++){
  
  // get the pointer to the jet already created
  Particle* jet = &fEvt->jet(iJet,jetType);  
  /// Spit out some useful information about this jet
  if (fEvt->Debug()) {
    cout << "EventBuilder_jetmet2012::AddBtag(): DEBUG ";
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
    if(!doCOMMON)  jet->Set("flavor_weight_MV2",                    Get<vector<float> >("jet_" + jetType + "_flavor_weight_MV2")                   .at(iJet));
    //jet->Set("flavor_weight_SoftMuonTag",            Get<vector<float> >( jetType + "_flavor_weight_SoftMuonTag")           .at(iJet));
    //jet->Set("flavor_weight_SecondSoftMuonTag",      Get<vector<float> >( jetType + "_flavor_weight_SecondSoftMuonTag")     .at(iJet));
    jet->Set("flavor_weight_SoftMuonTagChi2",        Get<vector<float> >("jet_" + jetType + "_flavor_weight_SoftMuonTagChi2")       .at(iJet));
    jet->Set("flavor_weight_SecondSoftMuonTagChi2",  Get<vector<float> >("jet_" + jetType + "_flavor_weight_SecondSoftMuonTagChi2") .at(iJet));

    jet->Set("svp_x",                                Get<vector<float> >("jet_" + jetType + "_flavor_component_svp_x")              .at(iJet));
    jet->Set("svp_y",                                Get<vector<float> >("jet_" + jetType + "_flavor_component_svp_y")              .at(iJet));
    jet->Set("svp_z",                                Get<vector<float> >("jet_" + jetType + "_flavor_component_svp_z")              .at(iJet));


    if(fEvt->isMC()){
      jet->Set("flavor_truth_label",                 Get<vector<int> >  ("jet_" + jetType + "_flavor_truth_label")                  .at(iJet));
    }

    // this does not work in ntup_jetmet

    // if(doTrack){ // now copy over the track assignments to the SV
    //   jet->AddVec("svp_track"); //flavor_component_svp_trk_index
    //   vector<int> trackIndices = Get<vector<vector<int> > >("jet_" + jetType + "_flavor_component_svp_trk_index")                  .at(iJet);
    //   if(fEvt->Debug()) cout << "Adding svp_track information, there are " << trackIndices.size() << " of them!" << endl;

    //   for(unsigned int iTr = 0; iTr < trackIndices.size(); iTr++){
    //     if(fEvt->Debug()) cout << "Adding svp_track information, iTr = " << iTr << " and trackIndices = " << trackIndices[iTr] << endl;

    //     jet->Add("svp_track",&fEvt->track(trackIndices[iTr]));
    //   }
    // }
  
  } // END FOR LOOP
  
  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::AddBtag(): DEBUG Scheduled b-tags for jet collection = " << jetType << " with " << fEvt->jets(jetType) << " jets" << endl;
  return kTRUE;  
  
} // END AddBtag


///=========================================
/// CopyJets
///=========================================
Bool_t EventBuilder_jetmet2012::CopyJets() 
{

  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012: DEBUG In CopyJets() " << endl;   
  if (fEvt->Debug()) cout << "doLCJets " << doLCJets << " doEMJets " << doEMJets << " doTruthJets " << doTruthJets << endl;  

  if(doLCJets){
    if(doJet4){
      CopyStandardJet(JET4LC);
    }
    if(doJet6){
      CopyStandardJet(JET6LC);
    }
  }
  if(doEMJets){
    if(doJet4){
      CopyStandardJet(JET4EM);
    }
    if(doJet6){
      CopyStandardJet(JET6EM);
    }
  }
  if(doTruthJets && fEvt->isMC()){
    if(doJet4){
      CopyStandardJet(JET4TRUTH);
    }
    if(doJet6){
      CopyStandardJet(JET6TRUTH);
    }
    if(doOOTtruthJet4){
      CopyStandardJet(JET4OOTTRUTH);
    }
    if(doInTimeTruthJet4){
      CopyStandardJet(JET4INTIMETRUTH);
    }
  }

  bool fail       = true;
  TString jetType = "";
  
  /// Copy all jet types
  for(Int_t i = 0; i<fEvt->Cfg()->Objs(mJETTYPEVEC); ++i) {
    TString type = ((TObjString*)fEvt->Cfg()->Obj(mJETTYPEVEC,i))->GetString();
    if (type.Contains("Truth") && !fEvt->isMC()) continue; 
    if(fEvt->Debug()) cout << "EventBuilder_jetmet2012::CopyJets(): DEBUG Copying jet type " << type << " which is the " << i << "th jet collection" << endl;
    if(fEvt->Debug()) cout << "EventBuilder_jetmet2012::CopyJets(): DEBUG type in vec is " << ((TObjString*) fEvt->Cfg()->Obj(mJETTYPEVEC,i))->GetString() << endl;
    if(type=="") continue;
    fail = fail && CopyJet(type,"");
    //delete type;
  }
  
  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::CopyJets(): DEBUG Scheduled All Jets" << endl;
  return fail;
}


///=========================================
/// CopyJet: Generic jet copier for non-standard jets. slower, but more versatile.
///=========================================
Bool_t EventBuilder_jetmet2012::CopyJet(TString jetType, TString jetClass) 
{

  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012: DEBUG In CopyJet() for jetType = " << jetType << "; jetClass = " << jetClass << endl;  
  
  // Handle the case where this is called without an argument
  if (jetType == "") { 
    Abort("EventBuilder_jetmet2012::CopyJet(): ERROR No jet type selected. Aborting!");
  }
  
  // Now try to figure out what type of jet this is
  if (jetClass == "") {
    if (fEvt->Debug()) {
      cout << "EventBuilder_jetmet2012::CopyJet(): INFO The jetClass variable is empty. ";
      cout << "Assuming that the class (truth, track, topo, etc) can be determined from the name." << endl;
    } // END IF
    
    // Parse the jetType: truth, track, topo
    if (jetType.Contains("Truth")) {
      if (fEvt->isMC()) jetClass = "truth";
      else { 
        cout << "EventBuilder_jetmet2012::CopyJet(): ERROR This is not MC and you have requested truth jets. Reconfigure analysis and run again!." << endl;
        return kFALSE;
      }
    }
    else if (jetType.Contains("Track"))  jetClass = "track";
    else if (jetType.Contains("LCTopo")) jetClass = "topo";
    else if (jetType.Contains("TopoEM")) jetClass = "topo";
    else if (jetType.Contains("tower"))  {
      cout << "EventBuilder_jetmet2012::CopyJet(): ERROR Tower jets are not currently supported." << endl;
      return kFALSE;
    } // END IF
    
  } // END IF JET CLASS 
  
  // Otherwise, try to figure out what class of jet this really is
  else if (jetClass != "truth" && jetClass != "topo" && jetClass != "track") {
    cout << "EventBuilder_jetmet2012::CopyJet(): ERROR Currently, the only supported classes of jets are: truth, track, topo." << endl;
    return kFALSE;
  }
 
  if(fEvt->Debug()){
    cout << " about to look for jetPre!" << endl;
  }

  // Collect jet information
  fEvt->AddVec("jets" + jetType);
  for(unsigned int iJet = 0; iJet < Get<vector<float> >( TString("jet_") + jetType + "_phi").size(); iJet++){
  
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
    cout << "EventBuilder_jetmet2012::CopyJet(): DEBUG ";
    TString info = 
    TString::Format("Jet %d: (E, Et, m, eta, phi) = (%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)",
            iJet, 
            jet->p.E(),jet->p.Et(),jet->p.M(),jet->p.Eta(),jet->p.Phi()
          );
   
    cout << info << endl;
  }
    
    
  // Handle reconstructed jets
  if (jetClass == "topo" ) {
    jet->Set("emscale_pt", Get<vector<float> >("jet_" + jetType + "_emscale_pt").at(iJet));
    jet->Set("constscale_pt", Get<vector<float> >("jet_" + jetType + "_constscale_pt").at(iJet));
    
    jet->Set("JVF",Get<vector<float> >("jet_" + jetType + "_jvtxf").at(iJet));  

  } // ENDIF topo

  // Handle track jets
  else if (jetClass == "track" ) {
    // add constits
    jet->AddVec("constit");
    vector<int> indices = Get<vector<vector<int> > >(TString("jet_") + jetType + "_constit_index").at(iJet);
    for(int iC = 0; iC < Get<vector<int> >(TString("jet_") + jetType + "_constit_n").at(iJet); iC++){
      Particle& constitTrack = fEvt->track(indices[iC]);
      jet->Add("constit", &constitTrack);
    }
  }
  
    // Handle truth jets
  else if (jetClass == "truth" ) {
    // add constits
    jet->AddVec("constit");
    vector<int> indices = Get<vector<vector<int> > >(TString("jet_") + jetType + "_constit_index").at(iJet);
    for(int iC = 0; iC < Get<vector<int> >(TString("jet_") + jetType + "_constit_n").at(iJet); iC++){
      Particle& constitTrack = fEvt->truth(indices[iC]);
      jet->Add("constit", &constitTrack);
    }
  }
  
  // Add the jet to the event
  fEvt->Add("jets" + jetType , jet);
  
  } // END FOR LOOP
  
  if (fEvt->Debug()) cout << "EventBuilder_jetmet2012::CopyJet(): DEBUG Scheduled jet collection = " << jetType << " with " << fEvt->jets(jetType) << " jets" << endl;
  return kTRUE;  
  
} // END CopyJet

///=========================================
/// CopyMuons
///=========================================
Bool_t EventBuilder_jetmet2012::CopyMuons()
{

  if(doMuons){
  static const MomKey mmuons("muons"), 
    mmuonsbase("muonsbase"), 
    mmuonscosmic("muonscosmic"), 
    mmuonsbad("muonsbad"),
    mTrigMatch("TrigMatch"),
    mTrigMatchNonIso("TrigMatchNonIso");

  static const BranchKey bmu_muid_n("mu_muid_n"), bmu_muid_pt("mu_muid_pt"), bmu_muid_eta("mu_muid_eta"), bmu_muid_phi("mu_muid_phi"), bmu_muid_E("mu_muid_E"),
    bmu_muid_me_qoverp_exPV("mu_muid_me_qoverp_exPV"), bmu_muid_id_qoverp_exPV("mu_muid_id_qoverp_exPV"),
    bmu_muid_me_theta_exPV("mu_muid_me_theta_exPV"), bmu_muid_id_theta_exPV("mu_muid_id_theta_exPV"), bmu_muid_id_theta("mu_muid_id_theta"),
    bmu_muid_charge("mu_muid_charge"),
    bmu_muid_isStandAloneMuon("mu_muid_isStandAloneMuon"),bmu_muid_isCombinedMuon("mu_muid_isCombinedMuon"), bmu_muid_isSegmentTaggedMuon("mu_muid_isSegmentTaggedMuon"),
    bmu_muid_loose("mu_muid_loose"), bmu_muid_medium("mu_muid_medium"), bmu_muid_tight("mu_muid_tight"), bmu_muid_expectBLayerHit("mu_muid_expectBLayerHit"), bmu_muid_nBLHits("mu_muid_nBLHits"), 
    bmu_muid_nPixHits("mu_muid_nPixHits"), bmu_muid_nPixelDeadSensors("mu_muid_nPixelDeadSensors"), bmu_muid_nPixHoles("mu_muid_nPixHoles"), 
    bmu_muid_nSCTHits("mu_muid_nSCTHits"),bmu_muid_nSCTDeadSensors("mu_muid_nSCTDeadSensors"), bmu_muid_nSCTHoles("mu_muid_nSCTHoles"),
    bmu_muid_nTRTHits("mu_muid_nTRTHits"), bmu_muid_nTRTOutliers("mu_muid_nTRTOutliers"), bmu_muid_ptcone20("mu_muid_ptcone20"), 
    bmu_muid_z0_exPV("mu_muid_z0_exPV"), bmu_muid_d0_exPV("mu_muid_d0_exPV"), bmu_muid_qoverp_exPV("mu_muid_qoverp_exPV"), 
    bmu_muid_cov_qoverp_exPV("mu_muid_cov_qoverp_exPV"), 
    bmu_muid_etcone30("mu_muid_etcone30"), bmu_muid_z0wrtPV("bmu_muid_z0wrtPV"), bmu_muid_author("mu_muid_author"),
    bmu_muid_id_phi("mu_muid_id_phi"), bmu_muid_ptcone30("mu_muid_ptcone30"), bmu_muid_etcone20("mu_muid_etcone20"),
    bmu_muid_trackd0pvunbiased("mu_muid_trackd0pvunbiased"), bmu_muid_trackz0pvunbiased("mu_muid_trackz0pvunbiased"),
    bmu_muid_isSiliconAssociatedForwardMuon("mu_muid_isSiliconAssociatedForwardMuon");

  static const MomKey me_qoverp_exPV("me_qoverp_exPV"), id_qoverp_exPV("id_qoverp_exPV"),
    me_theta_exPV("me_theta_exPV"), id_theta_exPV("id_theta_exPV"), id_theta("id_theta"),
    charge("charge"),
    isStandAloneMuon("isStandAloneMuon"), isCombinedMuon("isCombinedMuon"), isSegmentTaggedMuon("isSegmentTaggedMuon"),
    loose("loose"), medium("medium"), tight("tight"), expectBLayerHit("expectBLayerHit"), nBLHits("nBLHits"),
    nPixHits("nPixHits"), nPixelDeadSensors("nPixelDeadSensors"), nPixHoles("nPixHoles"), 
    nSCTHits("nSCTHits"),nSCTDeadSensors("nSCTDeadSensors"), nSCTHoles("nSCTHoles"),
    nTRTHits("nTRTHits"), nTRTOutliers("nTRTOutliers"), ptcone20("ptcone20"), 
    z0_exPV("z0_exPV"), d0_exPV("d0_exPV"), qoverp_exPV("qoverp_exPV"), 
    cov_qoverp_exPV("cov_qoverp_exPV"),  
    etcone30("Etcone30"), z0wrtPV("z0wrtPV"), author("author"),
    id_phi("id_phi"), ptcone30("Ptcone30"), etcone20("Etcone20"),
    d0pvunbiased("d0pvunbiased"), z0pvunbiased("z0pvunbiased"),
    isSiliconAssociatedForwardMuon("isSiliconAssociatedForwardMuon");

  static const BranchKey btrig_EF_trigmuonef_n("trig_EF_trigmuonef_n"), btrig_EF_trigmuonef_track_n("trig_EF_trigmuonef_track_n"),
    btrig_EF_trigmuonef_track_CB_pt("trig_EF_trigmuonef_track_CB_pt"), btrig_EF_trigmuonef_track_CB_eta("trig_EF_trigmuonef_track_CB_eta"),
    btrig_EF_trigmuonef_track_CB_phi("trig_EF_trigmuonef_track_CB_phi"), btrig_EF_trigmuonef_track_CB_hasCB("trig_EF_trigmuonef_track_CB_hasCB");

  static const BranchKey bEF_mu24i_tight("EF_mu24i_tight"), btrig_EF_trigmuonef_EF_mu24i_tight("trig_EF_trigmuonef_EF_mu24i_tight"), 
    bEF_mu36_tight("EF_mu36_tight"), btrig_EF_trigmuonef_EF_mu36_tight("trig_EF_trigmuonef_EF_mu36_tight"),
    bEF_mu24_tight("EF_mu24_tight"), btrig_EF_trigmuonef_EF_mu24_tight("trig_EF_trigmuonef_EF_mu24_tight");

  static vector<BranchKey> MuTrigSignatures;
  if(!MuTrigSignatures.size()) {
    MuTrigSignatures.push_back(btrig_EF_trigmuonef_EF_mu24i_tight);
    MuTrigSignatures.push_back(btrig_EF_trigmuonef_EF_mu36_tight);
  }
  
  static vector<BranchKey> MuTrigSignaturesNonIso;
  if(!MuTrigSignaturesNonIso.size()) {
    MuTrigSignaturesNonIso.push_back(btrig_EF_trigmuonef_EF_mu24_tight);
    MuTrigSignaturesNonIso.push_back(btrig_EF_trigmuonef_EF_mu36_tight);
  }
  
  //trigger
  bool onemuon = true;
  static const MomKey m1MUON("1MUON"), m2MUON("2MUON");
  fEvt->Cfg()->Get(m1MUON,onemuon);
  if(!onemuon) fEvt->Cfg()->Get(m2MUON,onemuon);

  fEvt->AddVec(mmuons);

  
  for(unsigned int iMu = 0; iMu<Get<vector<float> >(bmu_muid_pt).size(); ++iMu) {
    Particle* mu = new Particle();
    fEvt->Add(mmuons, mu);
    float pt  = Get<vector<float> >(bmu_muid_pt) .at(iMu);
    float eta = Get<vector<float> >(bmu_muid_eta).at(iMu);
    float phi = Get<vector<float> >(bmu_muid_phi).at(iMu);
    float e   = Get<vector<float> >(bmu_muid_E)  .at(iMu);




    mu->Set(me_qoverp_exPV,     Get<vector<float> >(bmu_muid_me_qoverp_exPV)     .at(iMu));
    mu->Set(id_qoverp_exPV,     Get<vector<float> >(bmu_muid_id_qoverp_exPV)     .at(iMu));
    mu->Set(me_theta_exPV,      Get<vector<float> >(bmu_muid_me_theta_exPV)      .at(iMu));
    mu->Set(id_theta_exPV,      Get<vector<float> >(bmu_muid_id_theta_exPV)      .at(iMu));  
    mu->Set(id_theta,           Get<vector<float> >(bmu_muid_id_theta)           .at(iMu));
    mu->Set(charge,             Get<vector<float> >(bmu_muid_charge)             .at(iMu));
    mu->Set(isStandAloneMuon,   Get<vector<int> >  (bmu_muid_isStandAloneMuon)   .at(iMu));
    mu->Set(isCombinedMuon,     Get<vector<int> >  (bmu_muid_isCombinedMuon)     .at(iMu)); 
    mu->Set(isSegmentTaggedMuon,Get<vector<int> >  (bmu_muid_isSegmentTaggedMuon).at(iMu));
    mu->Set(loose,              Get<vector<int> >  (bmu_muid_loose)              .at(iMu));
    mu->Set(medium,             Get<vector<int> >  (bmu_muid_medium)             .at(iMu));    
    mu->Set(tight,              Get<vector<int> >  (bmu_muid_tight)              .at(iMu)); 
    mu->Set(expectBLayerHit,    Get<vector<int> >  (bmu_muid_expectBLayerHit)    .at(iMu)); 
    mu->Set(nBLHits,            Get<vector<int> >  (bmu_muid_nBLHits)            .at(iMu)); 
    mu->Set(nPixHits,           Get<vector<int> >  (bmu_muid_nPixHits)           .at(iMu));
    mu->Set(nPixelDeadSensors,  Get<vector<int> >  (bmu_muid_nPixelDeadSensors)  .at(iMu)); 
    mu->Set(nPixHoles,          Get<vector<int> >  (bmu_muid_nPixHoles)          .at(iMu)); 
    mu->Set(nSCTHits,           Get<vector<int> >  (bmu_muid_nSCTHits)           .at(iMu));
    mu->Set(nSCTDeadSensors,    Get<vector<int> >  (bmu_muid_nSCTDeadSensors)    .at(iMu)); 
    mu->Set(nSCTHoles,          Get<vector<int> >  (bmu_muid_nSCTHoles)          .at(iMu)); 
    mu->Set(nTRTHits,           Get<vector<int> >  (bmu_muid_nTRTHits)           .at(iMu)); 
    mu->Set(nTRTOutliers,       Get<vector<int> >  (bmu_muid_nTRTOutliers)       .at(iMu));
    mu->Set(z0wrtPV,            Get<vector<float> >(bmu_muid_z0_exPV)            .at(iMu));
    mu->Set(etcone30,           Get<vector<float> >(bmu_muid_etcone30)           .at(iMu));
    mu->Set(etcone20,           Get<vector<float> >(bmu_muid_etcone20)           .at(iMu));
    mu->Set(ptcone20,           Get<vector<float> >(bmu_muid_ptcone20)           .at(iMu));
    mu->Set(ptcone30,           Get<vector<float> >(bmu_muid_ptcone30)           .at(iMu));
    mu->Set(author,             Get<vector<int  > >(bmu_muid_author)             .at(iMu));
    mu->Set(id_phi,             Get<vector<float> >(bmu_muid_id_phi)             .at(iMu));
    mu->Set(d0pvunbiased,       Get<vector<float> >(bmu_muid_trackd0pvunbiased)  .at(iMu));
    mu->Set(z0pvunbiased,       Get<vector<float> >(bmu_muid_trackz0pvunbiased)  .at(iMu));
    mu->Set(isSiliconAssociatedForwardMuon, Get<vector<int> >(bmu_muid_isSiliconAssociatedForwardMuon).at(iMu));


    const static MomKey ptUnsmeared("ptUnsmeared"), etaUnsmeared("etaUnsmeared"), phiUnsmeared("phiUnsmeared"), eUnsmeared("eUnsmeared");
    if(fEvt->isMC() && doMuSmear){ // need to do muon smearing

      float mypt = pt;

      float ptcb = pt;
      float ptid = (mu->Float(id_qoverp_exPV) != 0.) ? fabs(sin(mu->Float(id_theta_exPV))/mu->Float(id_qoverp_exPV)) : 0.;
      float ptms = (mu->Float(me_qoverp_exPV) != 0.) ? fabs(sin(mu->Float(me_theta_exPV))/mu->Float(me_qoverp_exPV)) : 0.;

      //if(ptms <= 0.)
        //cout << "starting pt " << mypt << " ptid = " << ptid << " ptms = " << ptms << " isCOmbined = " << mu->Bool(isCombinedMuon) << " isSegment " << mu->Bool(isSegmentTaggedMuon) << endl;
      int seed = int(fabs(phi*1.e+5));
      if(!seed) ++seed;
      MCPSmear->SetSeed(seed);

      std::string THESTRING = "";
      bool doSyst = false;
      switch(whichsyste) {
        case JetMETSyst::MMSDOWN:    THESTRING = "MSLOW";  doSyst = true; break;
        case JetMETSyst::MMSUP:      THESTRING = "MSUP";  doSyst = true;  break;
        case JetMETSyst::MIDDOWN:    THESTRING = "IDLOW"; doSyst = true; break;
        case JetMETSyst::MIDUP:      THESTRING = "IDUP"; doSyst = true;  break;
        case JetMETSyst::MSCALEDOWN: THESTRING = "SCALELOW"; doSyst = true; break;
        case JetMETSyst::MSCALEUP:   THESTRING = "SCALEUP"; doSyst = true;  break;
        default: break;
      }
      if (mu->Bool(isCombinedMuon) && ptms > 0.)
        MCPSmear->Event(ptms,ptid,ptcb,eta,mu->Float(charge));
      else if (mu->Bool(isSegmentTaggedMuon) || ptms <= 0.)
        MCPSmear->Event(ptid,eta,"ID",mu->Float(charge));
      else
        MCPSmear->Event(ptms,eta,"MS",mu->Float(charge));

      if (!doSyst) {
       if (mu->Bool(isCombinedMuon) && ptms > 0.)
          mypt = MCPSmear->pTCB();
        else if (mu->Bool(isSegmentTaggedMuon) || ptms <= 0.)
          mypt = MCPSmear->pTID();
        else
          mypt = MCPSmear->pTMS();
      } else {
        double pTMS_smeared = mypt;
        double pTID_smeared = mypt;
        double pTCB_smeared = mypt;

        // Valid values for "THESTRING":{"IDLOW", IDUP", "MSLOW", MSUP","SCALELOW", "SCALEUP"}
        MCPSmear->PTVar(pTMS_smeared, pTID_smeared, pTCB_smeared, THESTRING);

        if (mu->Bool(isCombinedMuon) && ptms > 0.) mypt = pTCB_smeared;
        else if (mu->Bool(isSegmentTaggedMuon) || ptms <= 0.) mypt = pTID_smeared;
        else mypt = pTMS_smeared;
      }
      //cout << " muon pt " << mypt << " started with " << pt << endl; 

      // if (mu->Bool(isCombinedMuon) && ptms > 0.){
      //   MCPSmear->Event(ptms,ptid,ptcb,eta,mu->Float(charge));
      //   mypt = MCPSmear->pTCB();
      // } else if (mu->Bool(isSegmentTaggedMuon) || ptms <= 0.){
      //   //cout << "segment tagged muon " << endl;
      //   MCPSmear->Event(ptid,eta,"ID",mu->Float(charge));
      //   mypt = MCPSmear->pTID();
      // } else {
      //   MCPSmear->Event(ptms,eta,"MS",mu->Float(charge));
      //   mypt = MCPSmear->pTMS();
      // }

      mu->p.SetPtEtaPhiE(mypt/1000.,eta,phi,e/1000.);

      // add tracks associated to muon (TLorentzVector obj)
      if(fEvt->Debug()) cout << "In CopyMuons: Adding tracks associated to muon" << endl;
      
      mu->AddVec("AssocTrack");
      float track_phi    = Get<vector<float> >("mu_muid_trackphi").at(iMu);
      float track_theta  = Get<vector<float> >("mu_muid_tracktheta").at(iMu);
      float track_qoverp = Get<vector<float> >("mu_muid_trackqoverp").at(iMu);
      float track_p = fabs(1./track_qoverp);
      TLorentzVector* track4V = new TLorentzVector();
      track4V->SetPtEtaPhiM(track_p*sin(track_theta)/1000., -1.0*log(tan(track_theta/2)), track_phi, 0.105658);
      mu->Add("AssocTrack", track4V);
   
      if(doLeptonTruth){
        // truth information
        if(fEvt->Debug()) cout << "In CopyMuons: Adding muon truth information" << endl;
      
        mu->Set("truth_matched",  Get<vector<int> >("mu_muid_truth_matched").at(iMu));
        mu->Set("truth_barcode", Get<vector<int> >("mu_muid_truth_barcode").at(iMu));
        //mu->Set("truth_index",    Get<vector<int> >("mu_muid_truth_index").at(iMu));
        bool mu_truth_matched = doBarcode ? barcode_link->Exists(TString::Itoa(mu->Int("truth_barcode"), 10)) : false;

        // truth obj has to be added here. index will be useless after everythign is sorted ...
        if(mu_truth_matched){
          mu->AddVec("TruthParticle", true);
          if(doBarcode){
            mu->Add("TruthParticle", barcode_link->Obj<Particle>(TString::Itoa(mu->Int("truth_barcode"), 10)));
          }
         
          //mu->AddVec("TruthParticle", true);
          //mu->Add("TruthParticle", &fEvt->truth(mu->Int("truth_index")));
        }

        if(fEvt->Debug()) cout << "In CopyMuons: End of muon truth information copy" << endl;
      }

      mu->Set(ptUnsmeared, pt/1000.);
      mu->Set(etaUnsmeared, eta);
      mu->Set(phiUnsmeared, phi);
      mu->Set(eUnsmeared, e/1000.);

      if(fEvt->Debug()){
        cout << "IN MUON SMEARING, STARTED WITH " << pt << " AND ENDED WITH " << mypt << endl;
      }

    } else { // data is easy
      mu->p.SetPtEtaPhiE(pt/1000.,eta,phi,e/1000.);

      mu->Set(ptUnsmeared, pt/1000.);
      mu->Set(etaUnsmeared, eta);
      mu->Set(phiUnsmeared, phi);
      mu->Set(eUnsmeared, e/1000.);
    }


    // trigger muon matching
    if(doMuonTriggers && doMuonTriggersMatch) {
      bool trigmatch = false;
      vector<BranchKey>::iterator iter = MuTrigSignatures.begin();
      int EFindex = -1;
      int EFtrackindex = -1;
      for(; iter!=MuTrigSignatures.end(); ++iter) {
        if(MuonHasTriggerMatch(
           Get<vector<float> >(bmu_muid_eta).at(iMu),
           Get<vector<float> >(bmu_muid_phi).at(iMu),
           GetP<vector<int> >(*iter),
           EFindex,EFtrackindex,Get<Int_t>(btrig_EF_trigmuonef_n),
           GetP<vector<vector<float> > >(btrig_EF_trigmuonef_track_CB_eta),
           GetP<vector<vector<float> > >(btrig_EF_trigmuonef_track_CB_phi),
           GetP<vector<vector<int> > >(btrig_EF_trigmuonef_track_CB_hasCB))) {
            trigmatch = true;
            break;              
          }
      }
      mu->Set(mTrigMatch,trigmatch);
    }
  }
  }
  return kTRUE;
}


///=========================================
/// CopyElectrons
///=========================================
Bool_t EventBuilder_jetmet2012::CopyElectrons() 
{
  if (fEvt->Debug()) cout << "EventBuilder_GRJETSv7_Hbb: DEBUG CopyElectrons() " << endl; 

  if(doElectrons) {
    fEvt->AddVec("electrons");
    for(int iEl = 0; iEl < Get<Int_t>("el_n"); iEl++) {
      
      Particle* el = new Particle();
 //      float pt  = Get<vector<float> >("el_pt") .at(iEl);
//       float eta = Get<vector<float> >("el_eta").at(iEl);
//       float phi = Get<vector<float> >("el_phi").at(iEl);
//       float e   = Get<vector<float> >("el_E")  .at(iEl);
//       el->p.SetPtEtaPhiE(pt/1000.,eta,phi,e/1000.);

      float cl_eta = Get<vector<float> >("el_cl_eta").at(iEl);
      float cl_phi = Get<vector<float> >("el_cl_phi").at(iEl);
      float cl_e   = Get<vector<float> >("el_cl_E")  .at(iEl);

      float trk_pt  = Get<vector<float> >("el_trackpt") .at(iEl);
      float trk_eta = Get<vector<float> >("el_tracketa").at(iEl);
      float trk_theta = Get<vector<float> >("el_tracktheta").at(iEl);
      float trk_phi = Get<vector<float> >("el_trackphi").at(iEl);

      float nPixHits  = Get<vector<int> >("el_nPixHits") .at(iEl);
      float nSCTHits  = Get<vector<int> >("el_nSCTHits") .at(iEl);


      float eta = (nPixHits+nSCTHits >= 4 ? trk_eta : cl_eta);
      float phi = (nPixHits+nSCTHits >= 4 ? trk_phi : cl_phi);
      float e   = cl_e;
//      float et  = e / cosh(eta);
      float pt  = sqrt( e*e - pow(0.511, 2.0)) / cosh(eta); // M_ele = 0.511, in MeV

      el->p.SetPtEtaPhiE(pt/1000.,eta,phi,e/1000.); //using et as proxy for pt... can fix later


      el->Set("cl_eta",  cl_eta);
      el->Set("cl_phi",  cl_phi);
      el->Set("cl_e",    cl_e);
      el->Set("trk_eta", trk_eta);
      el->Set("trk_theta", trk_theta);
      el->Set("trk_phi", trk_phi);
      el->Set("trk_pt",   trk_pt);
      el->Set("trk_qoverp",  Get<vector<float> >("el_trackqoverp").at(iEl));
      el->Set("d0_alt",          Get<vector<float> >("el_trackd0").at(iEl));
      el->Set("z0",          Get<vector<float> >("el_trackz0").at(iEl));
      
      if(!Rename){
        el->Set("Ptcone20",    Get<vector<float> >("el_ptcone20").at(iEl));
        el->Set("Ptcone30",    Get<vector<float> >("el_ptcone30").at(iEl));
      }
      else{
        el->Set("ptcone20",    Get<vector<float> >("el_ptcone20").at(iEl));
        el->Set("ptcone30",    Get<vector<float> >("el_ptcone30").at(iEl));
      }
            
      el->Set("Etcone20",    Get<vector<float> >("el_Etcone20").at(iEl));
      el->Set("Etcone30",    Get<vector<float> >("el_Etcone30").at(iEl));

      if(!Rename){
        el->Set("trackd0pvunbiased", Get<vector<float> >("el_trackd0pvunbiased").at(iEl));
        el->Set("trackz0pvunbiased", Get<vector<float> >("el_trackz0pvunbiased").at(iEl));
      }
      else{
        el->Set("d0pvunbiased", Get<vector<float> >("el_trackd0pvunbiased").at(iEl));
        el->Set("z0pvunbiased", Get<vector<float> >("el_trackz0pvunbiased").at(iEl));
      }

      //following variables needed to recompute electron PID quality, if needed
      el->Set("nPixHits", nPixHits);
      el->Set("nSCTHits", nSCTHits);
      el->Set("nSiHits", Get<vector<int> >("el_nSiHits").at(iEl));
      el->Set("nSCTOutliers", Get<vector<int> >("el_nSCTOutliers").at(iEl));
      el->Set("nPixelOutliers", Get<vector<int> >("el_nPixelOutliers").at(iEl));
      el->Set("nBLHits", Get<vector<int> >("el_nBLHits").at(iEl));
      el->Set("nBLayerOutliers", Get<vector<int> >("el_nBLayerOutliers").at(iEl));
      el->Set("expectHitInBLayer", Get<vector<float> >("el_expectHitInBLayer").at(iEl));
      el->Set("nTRTHits", Get<vector<int> >("el_nTRTHits").at(iEl));
      el->Set("TRTHighTOutliersRatio", Get<vector<float> >("el_TRTHighTOutliersRatio").at(iEl));
      el->Set("nTRTOutliers", Get<vector<int> >("el_nTRTOutliers").at(iEl));

      el->Set("etas2",  Get<vector<float> >("el_etas2").at(iEl));
      el->Set("Ethad", Get<vector<float> >("el_Ethad").at(iEl));
      el->Set("Ethad1", Get<vector<float> >("el_Ethad1").at(iEl));
      el->Set("reta", Get<vector<float> >("el_reta").at(iEl));
      el->Set("weta2", Get<vector<float> >("el_weta2").at(iEl));
      el->Set("f1", Get<vector<float> >("el_f1").at(iEl));
      el->Set("f3", Get<vector<float> >("el_f3").at(iEl)); //if not filled, use el_Es3/(el_Es0+el_Es1+el_Es2+el_Es3) 
      el->Set("f1core", Get<vector<float> >("el_f1core").at(iEl));
      el->Set("f3core", Get<vector<float> >("el_f3core").at(iEl)); 
      el->Set("wstot", Get<vector<float> >("el_wstot").at(iEl));
      el->Set("emaxs1", Get<vector<float> >("el_emaxs1").at(iEl));
      el->Set("Emax2", Get<vector<float> >("el_Emax2").at(iEl));
      el->Set("deltaeta1", Get<vector<float> >("el_deltaeta1").at(iEl));
      el->Set("deltaphi2", Get<vector<float> >("el_deltaphi2").at(iEl));
      el->Set("d0", Get<vector<float> >("el_trackd0_physics").at(iEl));

      //ph->Set("Etcone40_corrected",    Get<vector<float> >("ph_Etcone40_corrected")   .at(iEl));
      //ph->Set("Etcone40_ED_corrected", Get<vector<float> >("ph_Etcone40_ED_corrected").at(iEl));
      el->Set("author",                Get<vector<int> >("el_author")                .at(iEl));
      el->Set("isEM",                  Get<vector<unsigned int> >("el_isEM")         .at(iEl));
      el->Set("isEMLoose",             Get<vector<unsigned int> >("el_isEMLoose")    .at(iEl));
      el->Set("isEMMedium",            Get<vector<unsigned int> >("el_isEMMedium")   .at(iEl));
      el->Set("isEMTight",             Get<vector<unsigned int> >("el_isEMTight")    .at(iEl));
      el->Set("charge",                Get<vector<float> >       ("el_charge")       .at(iEl));

      
      el->Set("OQ",                    Get<vector<unsigned int> >("el_OQ")           .at(iEl));
      el->Set("convFlag",              Get<vector<int> >("el_convFlag")              .at(iEl));
      el->Set("isConv",                Get<vector<int> >("el_isConv")                .at(iEl));

      el->Set("mediumWithoutTrack",    Get<vector<int> >("el_mediumWithoutTrack")    .at(iEl));
      el->Set("mediumIsoWithoutTrack", Get<vector<int> >("el_mediumIsoWithoutTrack") .at(iEl));
      el->Set("tightWithoutTrack",     Get<vector<int> >("el_tightWithoutTrack")     .at(iEl));
      el->Set("tightIsoWithoutTrack",  Get<vector<int> >("el_tightIsoWithoutTrack")  .at(iEl));
      el->Set("loose",                 Get<vector<int> >("el_loose")                 .at(iEl));
      el->Set("looseIso",              Get<vector<int> >("el_looseIso")              .at(iEl));
      el->Set("medium",                Get<vector<int> >("el_medium")                .at(iEl));
      el->Set("mediumIso",             Get<vector<int> >("el_mediumIso")             .at(iEl));
      el->Set("tight",                 Get<vector<int> >("el_tight")                 .at(iEl));
      el->Set("tightIso",              Get<vector<int> >("el_tightIso")              .at(iEl));
      el->Set("loosePP",               Get<vector<int> >("el_loosePP")               .at(iEl));
      el->Set("loosePPIso",            Get<vector<int> >("el_loosePPIso")            .at(iEl));
      el->Set("mediumPP",              Get<vector<int> >("el_mediumPP")              .at(iEl));
      el->Set("mediumPPIso",           Get<vector<int> >("el_mediumPPIso")           .at(iEl));
      el->Set("tightPP",               Get<vector<int> >("el_tightPP")               .at(iEl));
      el->Set("tightPPIso",            Get<vector<int> >("el_tightPPIso")            .at(iEl));

      el->Set("goodOQ",                Get<vector<int> >("el_goodOQ")                .at(iEl));
      el->Set("refittedTrack_n",       Get<vector<int> >("el_refittedTrack_n")       .at(iEl));
      el->Set("vertweight",            Get<vector<float> >("el_vertweight")            .at(iEl));
      el->Set("hastrack",              Get<vector<int> >("el_hastrack")              .at(iEl));


      // int NewLoosePP= (int) isLoosePlusPlus(el->Float("etas2"), 
      //         el->Float("cl_e")/cosh(el->Float("etas2")),
      //         el->Float("Ethad")/(el->Float("cl_e")/cosh(el->Float("etas2"))), 
      //         el->Float("Ethad1")/(el->Float("cl_e")/cosh(el->Float("etas2"))), 
      //         el->Float("reta"), 
      //         el->Float("weta2"),
      //         el->Float("f1"),
      //         el->Float("wstot"), 
      //         ( el->Float("emaxs1") - el->Float("Emax2") )/ ( el->Float("emaxs1") + el->Float("Emax2") ), 
      //         el->Float("deltaeta1"), 
      //         el->Int("nSiHits"),
      //         el->Int("nSCTOutliers") + el->Int("nPixelOutliers"),
      //         el->Int("nPixHits"),
      //         el->Int("nPixelOutliers"));

      // int NewMediumPP = (int) isMediumPlusPlus(el->Float("etas2"), 
      //            el->Float("cl_e")/cosh(el->Float("etas2")),
      //            el->Float("f3"), //if not filled, use el_Es3/(el_Es0+el_Es1+el_Es2+el_Es3)
      //            el->Float("Ethad")/(el->Float("cl_e")/cosh(el->Float("etas2"))), 
      //            el->Float("Ethad1")/(el->Float("cl_e")/cosh(el->Float("etas2"))), 
      //            el->Float("reta"), 
      //            el->Float("weta2"),
      //            el->Float("f1"),
      //            el->Float("wstot"), 
      //            ( el->Float("emaxs1") - el->Float("Emax2") )/ ( el->Float("emaxs1") + el->Float("Emax2") ), 
      //            el->Float("deltaeta1"), 
      //            el->Float("d0"), 
      //            el->Float("TRTHighTOutliersRatio"),
      //            el->Int("nTRTHits"),
      //            el->Int("nTRTOutliers"),
      //            el->Int("nSiHits"),
      //            el->Int("nSCTOutliers") + el->Int("nPixelOutliers"),
      //            el->Int("nPixHits"),
      //            el->Int("nPixelOutliers"),
      //            el->Int("nBLHits"), 
      //            el->Int("nBLayerOutliers"), 
      //            el->Int("expectHitInBLayer"));

      // int NewTightPP = (int) isTightPlusPlus(el->Float("etas2"), 
      //          el->Float("cl_e")/cosh(el->Float("etas2")),
      //          el->Float("f3"), //if not filled, use el_Es3/(el_Es0+el_Es1+el_Es2+el_Es3)
      //          el->Float("Ethad")/(el->Float("cl_e")/cosh(el->Float("etas2"))), 
      //          el->Float("Ethad1")/(el->Float("cl_e")/cosh(el->Float("etas2"))), 
      //          el->Float("reta"), 
      //          el->Float("weta2"),
      //          el->Float("f1"),
      //          el->Float("wstot"), 
      //          ( el->Float("emaxs1") - el->Float("Emax2") )/ ( el->Float("emaxs1") + el->Float("Emax2") ), 
      //          el->Float("deltaeta1"), 
      //          el->Float("d0"), 
      //          el->Float("TRTHighTOutliersRatio"),
      //          el->Int("nTRTHits"),
      //          el->Int("nTRTOutliers"),
      //          el->Int("nSiHits"),
      //          el->Int("nSCTOutliers") + el->Int("nPixelOutliers"),
      //          el->Int("nPixHits"),
      //          el->Int("nPixelOutliers"),
      //          el->Int("nBLHits"), 
      //          el->Int("nBLayerOutliers"), 
      //          el->Int("expectHitInBLayer"),
      //          el->Float("cl_e") * fabs(el->Float("trk_qoverp")),
      //          el->Float("deltaphi2"),
      //          el->UInt("isEM") & (1 << egammaPID::ConversionMatch_Electron));
      

      // el->Set("loosePP_FromMacro", NewLoosePP);
      // el->Set("mediumPP_FromMacro", NewMediumPP);
      // el->Set("tightPP_FromMacro", NewTightPP);


      /*  
    eT: el_cl_E/cosh(el_etas2)
    eta: el_etas2
    rHad: el_Ethad/(el_cl_E/cosh(el_etas2))
    rHad1: el_Ethad1/(el_cl_E/cosh(el_etas2))
    Reta: el_reta
    w2: el_weta2
    f1: el_f1
    f3: el_f3 OR (if el_f3 is not filled) el_Es3/(el_Es0+el_Es1+el_Es2+el_Es3)
    wstot: el_wstot
    DEmaxs1: ( el_emaxs1 - el_Emax2 )/ ( el_emaxs1 + el_Emax2 )
    deltaEta: el_deltaeta1
    deltaPhi: el_deltaphi2
    eOverp: el_cl_E * fabs(el_trackqoverp)
    d0: el_trackd0_physics
    HTRatio, nTRT, nTRTOutliers: el_TRTHighTOutliersRatio , el_nTRTHits, el_nTRTOutliers
    nSi, nSiOutliers: el_nSiHits, el_nSCTOutliers + el_nPixelOutliers
    nPix, nPixOutliers: el_nPixHits, el_nPixelOutliers
    nBlayer, nBlayerOutliers: el_nBLHits, el_nBLayerOutliers
    expectBlayer: el_expectHitInBLayer
    ConvBit: el_isEM & (1 << egammaPID::ConversionMatch_Electron) NOTE: D3PD conversion quantities should not be used, you must use the isEM word) 
       */

      // add tracks associated to muon (TLorentzVector obj)
      if(fEvt->Debug()) cout << "In CopyElectrons: Adding tracks associated to electron" << endl;
      
      el->AddVec("AssocTrack");
      
      float track_phi    = Get<vector<float> >("el_trackphi").at(iEl);
      float track_eta    = Get<vector<float> >("el_tracketa").at(iEl);
      float track_pt     = Get<vector<float> >("el_trackpt"). at(iEl);
      TLorentzVector* track4V = new TLorentzVector();
      track4V->SetPtEtaPhiM(track_pt/1000., track_eta, track_phi, 0.000511);  // in GeV
      
      el->Add("AssocTrack", track4V);
     
      if(doLeptonTruth){
        // truth information
        if(fEvt->Debug()) cout << "In CopyElectrons: Adding electron truth information" << endl;
        
        el->Set("truth_matched",  Get<vector<int> >("el_truth_matched").at(iEl));
        el->Set("truth_index",    Get<vector<int> >("el_truth_index").at(iEl));
        el->Set("truth_barcode", Get<vector<int> >("el_truth_barcode").at(iEl));

        // it seems truth_matched is not reliable here ----> They do match. But according to the barcode, these are leptons out of Geant4 rather than MC generator. 
        bool el_truth_matched = doBarcode ? barcode_link->Exists(TString::Itoa(el->Int("truth_barcode"), 10)) : (el->Int("truth_index")!=-1);

        if(el_truth_matched){
          el->AddVec("TruthParticle", true);
          if(doBarcode){
            el->Add("TruthParticle", barcode_link->Obj<Particle>(TString::Itoa(el->Int("truth_barcode"), 10)));
          }
          else{
            el->Add("TruthParticle", &fEvt->truth(el->Int("truth_index")));
          }
        }

        if(fEvt->Debug()) cout << "In CopyElectrons: End of electron truth information copy" << endl;
      }

      fEvt->Add("electrons",el);
    }
  }
  
  if (fEvt->Debug()) cout << "EventBuilder_GRJETSv7_Hbb::CopyPhotons(): scheduled photons" << endl;
  return kTRUE;
}


//copied from BButler's hfsusy EventBuilder_susy1328.cxx, but should be generic
Bool_t EventBuilder_jetmet2012::DoHFOverlap(){

 if(fEvt->Debug()) cout << " doing hfor " << endl;
 if(doHFOR){
int datasetNumber = fEvt->ChannelNumber();
  
  static const MomKey mOutput("Output"), mEventWeight("EventWeight");
  static const BranchKey bmc_pt("mc_pt"), bmc_eta("mc_eta"), bmc_phi("mc_phi"), bmc_m("mc_m"), bmc_status("mc_status"), bmc_pdgId("mc_pdgId"), 
    bmc_vx_barcode("mc_vx_barcode"), bmc_parent_index("mc_parent_index"), bmc_child_index("mc_child_index"), bmc_n("mc_n");

  // Overlap removal for W/Z (alpgen) samples
  int hfor_type = -1;

  if ( (datasetNumber >= 110801 && datasetNumber <= 110804) ||   //wbb AlpgenPythia
    (datasetNumber >= 126601 && datasetNumber <= 126609) ||   //wc and wcc
    (datasetNumber >= 117680 && datasetNumber <= 117685) ||   //wenu
    (datasetNumber >= 117690 && datasetNumber <= 117695) ||   //wmunu
    (datasetNumber >= 117700 && datasetNumber <= 117705) ||   //wtaunu  
    (datasetNumber >= 110805 && datasetNumber <= 110828) ||   //zcc/zbb 
    (datasetNumber >= 117650 && datasetNumber <= 117655) ||   //zee
    (datasetNumber >= 117660 && datasetNumber <= 117665) ||   //zmumu
    (datasetNumber >= 117670 && datasetNumber <= 117675) ) {   //ztautau    

    myHforTool->setVerbosity(HforToolD3PD::ERROR);
        hfor_type = myHforTool->getDecision(datasetNumber, 
                          Get<Int_t>(bmc_n), 
                          GetP<vector<float> >(bmc_pt), 
                          GetP<vector<float> >(bmc_eta), 
                          GetP<vector<float> >(bmc_phi), 
                          GetP<vector<float> >(bmc_m), 
                                            GetP<vector<int> >(bmc_pdgId), 
                                            GetP<vector<int> >(bmc_status), 
                                            GetP<vector<int> >(bmc_vx_barcode), 
                                            GetP<vector<vector<int> > >(bmc_parent_index), 
                                            GetP<vector<vector<int> > >(bmc_child_index),
                                            HforToolD3PD::ALL);

         // Debugging purposes
      /*
      if (hfor_type ==  0) fEvt->Set("isBB",hfor_type);
      if (hfor_type ==  1) fEvt->Set("isCC",hfor_type); 
      if (hfor_type ==  2) fEvt->Set("isC",hfor_type); 
      if (hfor_type ==  3) fEvt->Set("isLightFlavor",hfor_type); 
      if (hfor_type ==  4) fEvt->Set("kill",hfor_type); 
      */
      
      if(hfor_type==4) {
        return kFALSE;
    }
  
    //Additional W weights only
      if ( (datasetNumber >= 110801 && datasetNumber <= 110804) ||   //wbb AlpgenPythia
    (datasetNumber >= 126601 && datasetNumber <= 126609) ||   //wc and wcc
    (datasetNumber >= 117680 && datasetNumber <= 117685) ||   //wenu
    (datasetNumber >= 117690 && datasetNumber <= 117695) ||   //wmunu
    (datasetNumber >= 117700 && datasetNumber <= 117705) ) {   //wtaunu   

        if (hfor_type == 0) {
          //fEvt->Set(mEventWeight,fEvt->Float(mEventWeight)*1.63*0.7); //additional k-factor for Wbb and Wcc, 0.7 added from stop analysis
          fEvt->Set(mOutput,"Wbb"); //redirect output, if file set up
        }
        else if (hfor_type == 1) {
          //fEvt->Set(mEventWeight,fEvt->Float(mEventWeight)*1.63*0.7); //additional k-factor for Wbb and Wcc, 0.7 added from stop analysis
          fEvt->Set(mOutput,"Wcc"); //redirect output, if file set up
        }
        else if (hfor_type ==  2) {
          //fEvt->Set(mEventWeight,fEvt->Float(mEventWeight)*1.11); //additional k-factor for Wc
          fEvt->Set(mOutput,"Wc"); //redirect output, if file set up
        }
    }
    else if ( (datasetNumber >= 110805 && datasetNumber <= 110828) ||   //zcc/zbb 
          (datasetNumber >= 117650 && datasetNumber <= 117655) ||   //zee
          (datasetNumber >= 117660 && datasetNumber <= 117665) ||   //zmumu
          (datasetNumber >= 117670 && datasetNumber <= 117675) ) {  //ztautau
      
        if (hfor_type == 0) {
          //fEvt->Set(mEventWeight,fEvt->Float(mEventWeight)*1.44); //additional k-factor Z+h.f.
          fEvt->Set(mOutput,"Zbb"); //redirect output, if file set up
        }
        else if (hfor_type == 1) {
          //fEvt->Set(mEventWeight,fEvt->Float(mEventWeight)*1.44); //additional k-factor Z+h.f.
          fEvt->Set(mOutput,"Zcc"); //redirect output, if file set up
        }
    }
  }
  
  //Sherpa
  if (datasetNumber >= 167740 && datasetNumber <= 167844) {
    bool isW = false;
    if ( (datasetNumber >= 167740 && datasetNumber <= 167748) ||
      (datasetNumber >= 167761 && datasetNumber <= 167796) ) isW = true;

    if (datasetNumber >= 167740 && datasetNumber <= 167760) { //inclusive
      //veto inclusive events with vector boson pT > 70
      
      int l1 = -1;
      int l2 = -1;  
      for(unsigned int i = 0; i<Get<vector<float> >(bmc_pt).size(); ++i) {
        if(Get<vector<int> >(bmc_status).at(i)!=3) continue;
        int pdgid = (int)fabs(Get<vector<int> >(bmc_pdgId).at(i));
        if(pdgid < 11 || pdgid > 16) continue;
        if(l1<0) {
          l1 = i;
          continue;
        }
        else {
          l2 = i;
          int pdgid1 = Get<vector<int> >(bmc_pdgId).at(l1);
          int pdgid2 = Get<vector<int> >(bmc_pdgId).at(l2);
          if((pdgid1*pdgid2) > 0) {
            Abort("EventBuilder_susy1328: ERROR lepton pair in V+jets cannot have same sign lepton number!");
          }
          if(isW) {
            if(fabs(fabs(pdgid1)-fabs(pdgid2))!=1) {
              Abort("EventBuilder_susy1328: ERROR lepton pair in W+jets must be the same flavor!");
            }
          }
          else {
            if(pdgid1 != -pdgid2) {
              Abort("EventBuilder_susy1328: ERROR lepton pair in Z+jets must be the same flavor!");
            }         
          }
          break;
        }
      }
      
      if(l1<0 || l2<0) Abort("EventBuilder_susy1328: ERROR two leptons not found!");
        
      TLorentzVector p1, p2;
      p1.SetPtEtaPhiM(Get<vector<float> >(bmc_pt).at(l1)/1000.,
        Get<vector<float> >(bmc_eta).at(l1),
        Get<vector<float> >(bmc_phi).at(l1),
        Get<vector<float> >(bmc_m).at(l1)/1000.);
      p2.SetPtEtaPhiM(Get<vector<float> >(bmc_pt).at(l2)/1000.,
        Get<vector<float> >(bmc_eta).at(l2),
        Get<vector<float> >(bmc_phi).at(l2),
        Get<vector<float> >(bmc_m).at(l2)/1000.);
        
      if((p1+p2).Pt()>70.) {
        return kFALSE;
      }
    }

    //Sort by output type, if desired
    int type = (datasetNumber - 167740) % 3;
    
    if(isW) {
      if(type == 0) {
          fEvt->Set(mOutput,"Wb"); //redirect output, if file set up    
      }
      else if(type == 1) {
          fEvt->Set(mOutput,"Wc"); //redirect output, if file set up    
      }
      //Type 2 is c and b veto, do not set  
    }
    else { //Z
    
      if(type == 0) {
          fEvt->Set(mOutput,"Zb"); //redirect output, if file set up    
      }
      else if(type == 1) {
          fEvt->Set(mOutput,"Zc"); //redirect output, if file set up    
      }
      //Type 2 is c and b veto, do not set
    }
  }

  }
  return kTRUE;
}


Bool_t EventBuilder_jetmet2012::MuonHasTriggerMatch(double mu_eta,
             double mu_phi,
             const vector<int> *trig_ef_trigmuonef_signature,
             int& efindex,
             int& eftrackindex,
             unsigned int trig_ef_trigmuonef_n,
             const vector<vector<float> > *trig_ef_trigmuonef_track_cb_eta,
             const vector<vector<float> > *trig_ef_trigmuonef_track_cb_phi,
             const vector<vector<int> > *trig_ef_trigmuonef_track_cb_hascb)
{
  /// code largely based on egammaanalysisutils/egammatriggermatching.cxx

  efindex = -1;
  eftrackindex = -1;

  if (trig_ef_trigmuonef_n==0) {
    //std::cout << "no ef object to be tested" << std::endl;
    return false;
  }
  
  /// checking consistency
  if (!trig_ef_trigmuonef_signature) {
    std::cout << "error: empty pointer provided for trig_ef_trigmuonef_signature .... returning false" << std::endl;
    return false;
  }
  if (trig_ef_trigmuonef_signature->size() != trig_ef_trigmuonef_n) {
    std::cout << "error: mismatch in the size of the trig_ef_trigmuonef_signature vector: " << trig_ef_trigmuonef_signature->size() << "  --> expected: " <<  trig_ef_trigmuonef_n << std::endl;
    return false;
  }

  if (!trig_ef_trigmuonef_track_cb_eta) {
    std::cout << "error: empty pointer provided for trig_ef_trigmuonef_track_cb_eta .... returning false" << std::endl;
    return false;
  }
  if (trig_ef_trigmuonef_track_cb_eta->size() != trig_ef_trigmuonef_n) {
    std::cout << "error: mismatch in the size of the trig_ef_trigmuonef_track_cb_eta vector: " << trig_ef_trigmuonef_track_cb_eta->size() << "  --> expected: " <<  trig_ef_trigmuonef_n << std::endl;
    return false;
  }

  if (trig_ef_trigmuonef_track_cb_phi==0) {
    std::cout << "error: empty pointer provided for trig_ef_trigmuonef_track_cb_phi .... returning false" << std::endl;
    return false;
  }
  if (trig_ef_trigmuonef_track_cb_phi->size() != trig_ef_trigmuonef_n) {
    std::cout << "error: mismatch in the size of the trig_ef_trigmuonef_track_cb_phi vector: " << trig_ef_trigmuonef_track_cb_phi->size() << "  --> expected: " <<  trig_ef_trigmuonef_n << std::endl;
    return false;
  }
  
  double drmax = 100;
 
  for (unsigned int imu=0; imu<trig_ef_trigmuonef_n; ++imu) {



    if ( !trig_ef_trigmuonef_signature->at(imu) ) continue;
    unsigned int ntracks = trig_ef_trigmuonef_track_cb_eta->at(imu).size();

    for (unsigned int itrack=0; itrack<ntracks; ++itrack) {
      if (!trig_ef_trigmuonef_track_cb_hascb->at(imu).at(itrack)) continue;
      double etaef = trig_ef_trigmuonef_track_cb_eta->at(imu).at(itrack);
      double phief = trig_ef_trigmuonef_track_cb_phi->at(imu).at(itrack);

      double dphi = TVector2::Phi_mpi_pi(phief - mu_phi);
      double deta = mu_eta - etaef;
      double deltar=sqrt(dphi*dphi + deta*deta);
      if ( deltar < drmax ) {
        drmax = deltar;
        efindex = imu;
        eftrackindex = itrack;
      }
    }
  }
 
  if ( drmax <= 0.15 ) {
    return true;
  } 

  efindex = -1;
  eftrackindex = -1;
  return false;
}



Bool_t EventBuilder_jetmet2012::CopyMuonTriggers()
{


  if(doMuonTriggers){

    fEvt->Set("EF_mu24i_tight", Get<bool>("EF_mu24i_tight") );
    fEvt->Set("EF_mu36_tight", Get<bool>("EF_mu36_tight") );
    
    static const MomKey m1muontrigger("1muontrigger");
    
    static const BranchKey bEF_mu24i_tight("EF_mu24i_tight");
    static const BranchKey bEF_mu36_tight("EF_mu36_tight");

    if(Get<Bool_t>(bEF_mu24i_tight)) fEvt->Set(m1muontrigger,true);
    else fEvt->Set(m1muontrigger,Get<Bool_t>(bEF_mu36_tight));
  }
  return kTRUE;
}

Bool_t EventBuilder_jetmet2012::CopyElectronTriggers()
{
  if(doElectronTriggers){
    //***   = =|||  ***//
  
  
    fEvt->Set("EF_e24vhi_medium1", Get<bool>("EF_e24vhi_medium1") );
    fEvt->Set("EF_e60_medium1", Get<bool>("EF_e60_medium1") );
  
  }

  return kTRUE;
}

Bool_t EventBuilder_jetmet2012::KillPileupFluctuations(){
  if(fEvt->isMC()){
    float maxpt6 = 0.;
    static const BranchKey bJet6PT("jet_AntiKt6TopoEM_pt");
    if(Get<vector<float> >(bJet6PT).size() > 0) maxpt6 = Get<vector<float> >(bJet6PT).at(0)/1000;   

    static const MomKey mChannelNumber("ChannelNumber");
    int chan = fEvt->Int(mChannelNumber);
    if(((chan == 105009 || chan == 113204 || chan == 126127) && maxpt6 > 0) || 
       ((chan == 105010 || chan == 113205 || chan == 126128) && maxpt6 > 70) || 
       ((chan == 105011 || chan == 113206 || chan == 126129) && maxpt6 > 140) ||
       ((chan == 105012 || chan == 113207 || chan == 126130) && maxpt6 > 240) || // max : 210 instead of 240
       ((chan == 105013 || chan == 113208 || chan == 126131) && maxpt6 > 440) ||
       ((chan == 105014 || chan == 113209 || chan == 126132) && maxpt6 > 800) ||
       ((chan == 105015 || chan == 113210 || chan == 126133) && maxpt6 > 1300) ||
       ((chan == 105016 || chan == 126134) && maxpt6 > 2400) ){
        if(fEvt->Debug()) 
          cout << "MC event with high-pt pileup jet!" << endl;
        return kFALSE;
    }

    if(chan >= 108081 && chan <= 108087){
      static const BranchKey bmcfEvt_pdf_scale("mcfEvt_pdf_scale");
      double pthat = Get<vector<double> >(bmcfEvt_pdf_scale)[0];

      if((chan == 108087 && (pthat > 45 || pthat < 0)) ||
         (chan == 108081 && (pthat > 85 || pthat < 45)) ||
         (chan == 108082 && (pthat > 150 || pthat < 85)) ||
         (chan == 108083 && (pthat > 310 || pthat < 150)) ||
         (chan == 108084 && (pthat > 1000000 || pthat < 310))){ // NEEDS WORK

        if(fEvt->Debug()) 
          cout << "killed event with pthat " << pthat << " for channel " << chan << endl;
        return kFALSE;
      } 
    }
    
    if(chan >= 147910 && chan <= 147917){
      bool status = true;
      static const BranchKey bJet4PT("jet_AntiKt4LCTopo_pt");
      BranchKey lead;
      if(doSMWZ)
        lead = "jet_";
      else
        lead = "";
      BranchKey bJet4PTT(lead+"AntiKt4Truth_pt");

      if(Get<vector<float> >(bJet4PT).size() == 0 || Get<vector<float> >(bJet4PT).size() == 1 || Get<vector<float> >(bJet4PTT).size() == 0){
        status = true;
      }
      else{
        float ptave = (Get<vector<float> >(bJet4PT).at(0)/1000. + Get<vector<float> >(bJet4PT).at(1)/1000.)/2.;
        float pttrue = (Get<vector<float> >(bJet4PTT).at(0)/1000.);

        if(ptave/pttrue < 1.4)
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
  return kTRUE;
}

Bool_t EventBuilder_jetmet2012::GetJetCalibrationTools(){
  cout << "In EventBuilder_jetmet2012::GetJetCalibrationTools(), doJetCalibrations: " << doJetCalibrations << endl;


    data4EMCalibrator = (JetCalibrationTool*) fInput->FindObject("data4EMCalibrator");
    data6EMCalibrator = (JetCalibrationTool*) fInput->FindObject("data6EMCalibrator");

    data4LCCalibrator = (JetCalibrationTool*) fInput->FindObject("data4LCCalibrator");
    data6LCCalibrator = (JetCalibrationTool*) fInput->FindObject("data6LCCalibrator");

    mc4EMCalibrator = (JetCalibrationTool*) fInput->FindObject("mc4EMCalibrator");
    mc6EMCalibrator = (JetCalibrationTool*) fInput->FindObject("mc6EMCalibrator");
    
    mc4LCCalibrator = (JetCalibrationTool*) fInput->FindObject("mc4LCCalibrator");
    mc6LCCalibrator = (JetCalibrationTool*) fInput->FindObject("mc6LCCalibrator");

    afiimc4EMCalibrator = (JetCalibrationTool*) fInput->FindObject("afiimc4EMCalibrator");
    afiimc4EMCalibrator = (JetCalibrationTool*) fInput->FindObject("afiimc4EMCalibrator");
    afiimc6EMCalibrator = (JetCalibrationTool*) fInput->FindObject("afiimc6EMCalibrator");
    afiimc6EMCalibrator = (JetCalibrationTool*) fInput->FindObject("afiimc6EMCalibrator");
  
    m_myJesNorm = (MultijetJESUncertaintyProvider*) fInput->FindObject("JES");
    m_myJesAFII = (MultijetJESUncertaintyProvider*) fInput->FindObject("JESA");

    m_myJER = (JetSmearingTool*) fInput->FindObject("JER");

  return kTRUE;

}

// function to copy only standard jet types
// use other copy functions to copy types not included here (using vectors to specify)
Bool_t EventBuilder_jetmet2012::CopyStandardJet(int JetFlag){
  MomKey mjets;
  MomKey mjetsNoJ;
  BranchKey bjets;
  JetCalibrationTool* myCal = 0;

  MomKey lead;

  // set up the type of jet we're dealing with based on the enum
  if(fEvt->Debug()) cout << "JetFlag " << JetFlag << endl;
  if(fEvt->Debug()) cout << "isMC " << fEvt->isMC() << endl;
  if(fEvt->Debug()) cout << "doAFII " << doAFII << endl;
  
  switch(JetFlag){
    case JET4LC:
      lead = MomKey("jet_");
      mjets = MomKey("jetsAntiKt4LCTopo");
      mjetsNoJ = MomKey("AntiKt4LCTopo");
      bjets = BranchKey("AntiKt4LCTopo");
      if(!fEvt->isMC()){
        myCal = data4LCCalibrator;
      } else if(doAFII){
        myCal = afiimc4LCCalibrator;
      } else {
        myCal = mc4LCCalibrator;
      }
      break;
    case JET6LC:
      lead = MomKey("jet_");
      mjets = MomKey("jetsAntiKt6LCTopo");
      mjetsNoJ = MomKey("AntiKt6LCTopo");
      bjets = BranchKey("AntiKt6LCTopo");
      if(!fEvt->isMC()){
        myCal = data6LCCalibrator;
      } else if(doAFII){
        myCal = afiimc6LCCalibrator;
      } else {
        myCal = mc6LCCalibrator;
      }
      break;
    case JET4EM:
      lead = MomKey("jet_");
      mjets = MomKey("jetsAntiKt4TopoEM");
      mjetsNoJ = MomKey("AntiKt4TopoEM");
      bjets = BranchKey("AntiKt4TopoEM");
      if(!fEvt->isMC()){
        myCal = data4EMCalibrator;
      } else if(doAFII){
        myCal = afiimc4EMCalibrator;
      } else {
        myCal = mc4EMCalibrator;
      }
      break;
    case JET6EM:
      lead = MomKey("jet_");
      mjets = MomKey("jetsAntiKt6TopoEM");
      mjetsNoJ = MomKey("AntiKt6TopoEM");
      bjets = BranchKey("AntiKt6TopoEM");
      if(!fEvt->isMC()){
        myCal = data6EMCalibrator;
      } else if(doAFII){
        myCal = afiimc6EMCalibrator;
      } else {
        myCal = mc6EMCalibrator;
      }
      break;
    case JET4TRUTH:
      if(doSMWZ)
        lead = "jet_";
      else
        lead = MomKey("");
      mjets = MomKey("jetsAntiKt4Truth");
      mjetsNoJ = MomKey("AntiKt4Truth");
      bjets = BranchKey("AntiKt4Truth");
      myCal = 0;
      break;
    case JET6TRUTH:
      if(doSMWZ)
        lead = "jet_";
      else
        lead = MomKey("");
      mjets = MomKey("jetsAntiKt6Truth");
      mjetsNoJ = MomKey("AntiKt6Truth");
      bjets = BranchKey("AntiKt6Truth");
      myCal = 0;
      break;
    case JET4OOTTRUTH:
      if(doSMWZ)
        lead = "jet_";
      else
        lead = MomKey("");
      mjets = MomKey("jetsOutOfTimeAntiKt4Truth");
      mjetsNoJ = MomKey("OutOfTimeAntiKt4Truth");
      bjets = BranchKey("OutOfTimeAntiKt4Truth");
      myCal = 0;
      break;
    case JET4INTIMETRUTH:
      if(doSMWZ)
        lead = "jet_";
      else
        lead = MomKey("");
      mjets = MomKey("jetsInTimeAntiKt4Truth");
      mjetsNoJ = MomKey("InTimeAntiKt4Truth");
      bjets = BranchKey("InTimeAntiKt4Truth");
      myCal = 0;
      break;
  }
 
  //cout << myCal << " for type " << JetFlag << endl;
  BranchKey fix;
  if(doCOMMON){
    fix = BranchKey("_allpv_1GeV");
  } else {
    fix = BranchKey("");
  }

  // ugly task of declaring all the branchkeys we want
  const BranchKey bjet_JetType_pt(lead+bjets+"_pt"), bjet_JetType_eta(lead+bjets+"_eta"), 
    bjet_JetType_phi(lead+bjets+"_phi"), bjet_JetType_E(lead+bjets+"_E"),
    bjet_JetType_emscale_pt(lead+bjets+"_emscale_pt"),
    bjet_JetType_emscale_eta(lead+bjets+"_emscale_eta"),
    bjet_JetType_emscale_phi(lead+bjets+"_emscale_phi"),
    bjet_JetType_emscale_E(lead+bjets+"_emscale_E"),
    bjet_JetType_constscale_eta(lead+bjets+"_constscale_eta"), 
    bjet_JetType_constscale_phi(lead+bjets+"_constscale_phi"),
    bjet_JetType_constscale_pt(lead+bjets+"_constscale_pt"),
    bjet_JetType_constscale_E(lead+bjets+"_constscale_E"), 
    bjet_JetType_constscale_m(lead+bjets+"_constscale_m"), 
    bjet_JetType_ActiveAreaPx(lead+bjets+"_ActiveAreaPx"), 
    bjet_JetType_ActiveAreaPy(lead+bjets+"_ActiveAreaPy"), 
    bjet_JetType_ActiveAreaPz(lead+bjets+"_ActiveAreaPz"), 
    bjet_JetType_ActiveAreaE(lead+bjets+"_ActiveAreaE"),
    bjet_JetType_ActiveArea(lead+bjets+"_ActiveArea"),
    bjet_JetType_emfrac(lead+bjets+"_emfrac"),
    bjet_JetType_hecf(lead+bjets+"_hecf"), bjet_JetType_LArQuality(lead+bjets+"_LArQuality"),
    bjet_JetType_HECQuality(lead+bjets+"_HECQuality"), bjet_JetType_AverageLArQF(lead+bjets+"_AverageLArQF"),
    bjet_JetType_Timing(lead+bjets+"_Timing"), bjet_JetType_sumPtTrk(lead+bjets+"_sumPtTrk"+fix),
    bjet_JetType_fracSamplingMax(lead+bjets+"_fracSamplingMax"), bjet_JetType_SamplingMax(lead+bjets+"_SamplingMax"),
    bjet_JetType_NegativeE(lead+bjets+"_NegativeE"),
    bjet_JetType_flavor_truth_label(lead+bjets+"_flavor_truth_label"), bjet_JetType_m(lead+bjets+"_m"), 
    bjet_JetType_YFlip12(lead+bjets+"_YFlip12"), bjet_JetType_YFlip23(lead+bjets+"_YFlip23"),
    bjet_JetType_flavor_weight_SV0(lead+bjets+"_flavor_weight_SV0"), 
    bjet_JetType_flavor_weight_JetFitterCOMBNN(lead+bjets+"_flavor_weight_JetFitterCOMBNN"), 
    bjet_JetType_flavor_weight_IP3D(lead+bjets+"_flavor_weight_IP3D"),
    bjet_JetType_flavor_weight_SV1(lead+bjets+"_flavor_weight_SV1"),
    bjet_JetType_flavor_weight_MV1(lead+bjets+"_flavor_weight_MV1"),
    bjet_JetType_jvtxf(lead+bjets+"_jvtxf"), 
    bjet_JetType_jvtxfFull(lead+bjets+"_jvtxfFull"), 
    bjet_JetType_isBadLoose(lead+bjets+"_isBadLoose"),
    bjet_JetType_isBadLooseMinus(lead+bjets+"_isBadLooseMinus"),
    bjet_JetType_EtaOrigin(lead+bjets+"_EtaOrigin"),
    bjet_JetType_PhiOrigin(lead+bjets+"_PhiOrigin"),
    bjet_JetType_MOrigin(lead+bjets+"_MOrigin"),
    bjet_JetType_BCH_CORR_JET(lead+bjets+"_BCH_CORR_JET"),
    bjet_JetType_nTrackAssoc(lead+bjets+"_TrackAssoc_n"),
    bjet_JetType_TrackAssoc_index(lead+bjets+"_TrackAssoc_index"),
    bjet_JetType_WIDTH(lead+bjets+"_WIDTH"), bjet_JetType_e_TileGap1(lead+bjets+"_e_TileGap1"),
    bjet_JetType_e_TileGap2(lead+bjets+"_e_TileGap2"),  bjet_JetType_e_TileGap3(lead+bjets+"_e_TileGap3"),
    bjet_JetType_nConstitAssoc(lead+bjets+"_constit_n"),
    bjet_JetType_ConstitAssoc_index(lead+bjets+"_constit_index"),
    bjet_JetType_TruthMFindex(lead+bjets+"_TruthMFindex"), bjet_JetType_TruthMF(lead+bjets+"_TruthMF"),
    bjet_JetType_NumTowers(lead+bjets+"_NumTowers"),
    bjet_JetType_sumPtTrk_pv0_500MeV(lead+bjets+"_sumPtTrk_pv0_500MeV");
   

    // ugly task of declaring all the moment keys

  static const MomKey emscale_pt("emscale_pt"),
    emscale_eta("emscale_eta"),
    emscale_phi("emscale_phi"),
    emscale_E("emscale_E"),
    constscale_eta("constscale_eta"), 
    constscale_phi("constscale_phi"),
    constscale_pt("constscale_pt"),
    constscale_E("constscale_E"), 
    constscale_m("constscale_m"), 
    ActiveAreaPx("ActiveAreaPx"), 
    ActiveAreaPy("ActiveAreaPy"), 
    ActiveAreaPz("ActiveAreaPz"), 
    ActiveAreaE("ActiveAreaE"),
    ActiveArea("ActiveArea"),
    emfrac("emfrac"),
    hecf("hecf"), LArQuality("LArQuality"),
    HECQuality("HECQuality"), AverageLArQF("AverageLArQF"),
    Timing("Timing"), sumPtTrk("sumPtTrk"),
    fracSamplingMax("fracSamplingMax"), SamplingMax("SamplingMax"),
    NegativeE("NegativeE"),
    flavor_truth_label("flavor_truth_label"), m("m"), 
    YFlip12("YFlip12"), YFlip23("YFlip23"),
    flavor_weight_SV0("flavor_weight_SV0"), 
    flavor_weight_JetFitterCOMBNN("flavor_weight_JetFitterCOMBNN"), 
    flavor_weight_IP3D("flavor_weight_IP3D"),
    flavor_weight_SV1("flavor_weight_SV1"),
    flavor_weight_MV1("flavor_weight_MV1"),
    JVF("JVF"), JVFFull("JVFfull"),
    isBadLoose("isBadLoose"),
    isBadLooseMinus("isBadLooseMinus"),
    EtaOrigin("EtaOrigin"),
    PhiOrigin("PhiOrigin"),
    MOrigin("MOrigin"),
    averageIntPerXing("averageIntPerXing"), Eventshape_rhoKt4LC("Eventshape_rhoKt4LC"),
    nTrackAssoc("nTrackAssoc"), GhostAssocTrack("GhostAssocTrack"), GhostAssocTrackSumPt("GhostAssocTrack_SumPt"),
    BCH_CORR_JET("BCH_CORR_JET"), 
    WIDTH("WIDTH"), TileGap1("TileGap1"), TileGap2("TileGap2"), TileGap3("TileGap3"),
    Pt_woJES("Pt_woJES"), Eta_woJES("Eta_woJES"), Phi_woJES("Phi_woJES"), E_woJES("E_woJES"), M_woJES("M_woJES"),
    Jet_woJES("Jet_woJES"),
    nConstitAssoc("nconstituents"), AssocConstit("constituents"), AssocConstitSumPt("constituents_SumPt"),
    TruthMFLinks("TruthMFLinks"), NumTowers("NumTowers"), sumPtTrk_pv0_500MeV("sumPtTrk_pv0_500MeV");

  fEvt->AddVec(mjets);
  if(fEvt->Debug()) cout << "In CopyStandardJet:: adding jet type " << mjets << endl ;

  // loop and copy
  for(unsigned int iJet = 0; iJet < Get<vector< float> >(bjet_JetType_m).size(); iJet++){
    if(fEvt->Debug()) cout << iJet << " ";

    Particle* jet = new Particle();


    if(fEvt->Debug()) cout << "doJetCalibrations " << doJetCalibrations << endl;
    if(fEvt->Debug()) cout << "myCal " << myCal << endl;
    // based on whether we're doing the calibration... do the calibration
    if(doJetCalibrations && myCal!=0){
      Particle* jet_woJES = new Particle();

      myCal->UseGeV(true);

      if(fEvt->Debug()) cout << "out myCal is type " << myCal->GetName() << endl;
    
      float ax  =     Get<vector<float> >(bjet_JetType_ActiveAreaPx)  .at(iJet);
      float ay  =     Get<vector<float> >(bjet_JetType_ActiveAreaPy)  .at(iJet);
      float az  =     Get<vector<float> >(bjet_JetType_ActiveAreaPz)  .at(iJet);
      float ae  =     Get<vector<float> >(bjet_JetType_ActiveAreaE)   .at(iJet);
      float pt  =     Get<vector<float> >(bjet_JetType_constscale_pt) .at(iJet);
      float eta =     Get<vector<float> >(bjet_JetType_constscale_eta).at(iJet);
      float phi =     Get<vector<float> >(bjet_JetType_constscale_phi).at(iJet);
      float m   =     Get<vector<float> >(bjet_JetType_constscale_m)  .at(iJet);
      float e   =     Get<vector<float> >(bjet_JetType_constscale_E)  .at(iJet);

      jet->Set(constscale_pt,        pt/1000.);
      jet->Set(constscale_eta,       eta);
      jet->Set(constscale_phi,       phi);
      jet->Set(constscale_m,         m/1000.);

      float rho = fEvt->Float(Eventshape_rhoKt4LC)/1000.;
      float mu  = (float) fEvt->Int(averageIntPerXing);
      float npv = (float) fEvt->vtxs();

      
      jet->p = myCal->ApplyJetAreaOffsetEtaJES(e/1000., eta, phi, m/1000., ax, ay, az, ae,
                                               rho, mu, npv); 


      if(whichsyste==JetMETSyst::JER)
      {
        if (fabs(jet->p.Eta()) < 4.5 && jet->p.Perp() > 20.){
          int seed = int(fabs(jet->p.Phi()*1.e+5));
          if(!seed) ++seed;
          m_myJER->SetSeed(seed);
          /// smear jet
          TLorentzVector tmp_tlv = jet->p;
          if(doAFII) m_myJER->SmearJet_Syst_AFII(tmp_tlv);
          else m_myJER->SmearJet_Syst(tmp_tlv);
          jet->p = tmp_tlv;
        }
      }



      jet_AntiKt4LCTopo_recalibrated_pt.push_back(jet->p.Pt());
      calJets.push_back(jet->p);




      // Store Pt without JES (but after area subtraction)
      TLorentzVector p_woJES = myCal->ApplyJetAreaOffset(e/1000., eta, phi, m/1000., ax, ay, az, ae,
                                                         fEvt->Float(Eventshape_rhoKt4LC)/1000.,
					                                               fEvt->Int(averageIntPerXing),
				                                                 fEvt->vtxs());

      TLorentzVector areaVector(0.,0.,0.,0.);
      areaVector.SetPxPyPzE(ax, ay, az, ae);
  
      if(fEvt->Debug()) cout << " before calib was " << pt/1000. << " " << eta  << " " << m/1000. << " after calib is " << jet->p.Pt() << " " << jet->p.Eta() <<  " " << jet->p.M() <<   " without JES is " << p_woJES.Pt() << " " << p_woJES.Eta() <<  " " <<  p_woJES.M() << " mu = " << fEvt->Float(averageIntPerXing) << " vtxs = " << fEvt->vtxs() <<  "rho is " << fEvt->Float("Eventshape_rhoKt4LC")/1000. << endl;

      if(fEvt->Debug()) cout << "ax = " << ax << " ay  = " << ay << " az = " << az << " ae = " << ae << "aPT = " << areaVector.Pt() << endl;
  
      jet_woJES->p = p_woJES;
      jet->AddVec(Jet_woJES);
      jet->Add(Jet_woJES, jet_woJES); 
      

    } else if(JetFlag!= JET6TRUTH && JetFlag!=JET4TRUTH && JetFlag!=JET4OOTTRUTH && JetFlag!=JET4INTIMETRUTH ){ // truth jets, and when we turned off cal, gets thrown here
      float ptJ  =     Get<vector<float> >(bjet_JetType_constscale_pt).at(iJet)/1000.;
      float etaJ =     Get<vector<float> >(bjet_JetType_constscale_eta).at(iJet);
      float phiJ =     Get<vector<float> >(bjet_JetType_constscale_phi).at(iJet);
      float mJ   =     Get<vector<float> >(bjet_JetType_constscale_m).at(iJet)/1000.;
      jet->p.SetPtEtaPhiM(ptJ, etaJ, phiJ, mJ);

      jet->Set(constscale_pt,  ptJ);
      jet->Set(constscale_eta, etaJ);
      jet->Set(constscale_phi, phiJ);
      jet->Set(constscale_m,   mJ);
    } else {
      float ptJ  =     Get<vector<float> >(bjet_JetType_pt).at(iJet)/1000.;
      float etaJ =     Get<vector<float> >(bjet_JetType_eta).at(iJet);
      float phiJ =     Get<vector<float> >(bjet_JetType_phi).at(iJet);
      float mJ   =     Get<vector<float> >(bjet_JetType_m).at(iJet)/1000.;
      jet->p.SetPtEtaPhiM(ptJ, etaJ, phiJ, mJ);

      jet->Set(constscale_pt,  ptJ);
      jet->Set(constscale_eta, etaJ);
      jet->Set(constscale_phi, phiJ);
      jet->Set(constscale_m,   mJ);


    }// end doJetCalibration blocks


    if(JetFlag != JET6TRUTH && JetFlag != JET4TRUTH && JetFlag!=JET4OOTTRUTH && JetFlag!=JET4INTIMETRUTH){
      jet->Set(emfrac,         Get<vector<float> >(bjet_JetType_emfrac)         .at(iJet));
      jet->Set(hecf,           Get<vector<float> >(bjet_JetType_hecf)           .at(iJet));
      jet->Set(LArQuality,     Get<vector<float> >(bjet_JetType_LArQuality)     .at(iJet));
      jet->Set(HECQuality,     Get<vector<float> >(bjet_JetType_HECQuality)     .at(iJet));
      jet->Set(AverageLArQF,   Get<vector<float> >(bjet_JetType_AverageLArQF)   .at(iJet));
      jet->Set(Timing,         Get<vector<float> >(bjet_JetType_Timing)         .at(iJet));
      jet->Set(sumPtTrk,       Get<vector<float> >(bjet_JetType_sumPtTrk)       .at(iJet));
      jet->Set(fracSamplingMax,Get<vector<float> >(bjet_JetType_fracSamplingMax).at(iJet));
      jet->Set(SamplingMax,    Get<vector<int> >  (bjet_JetType_SamplingMax)    .at(iJet));
      jet->Set(NegativeE,      Get<vector<float> >(bjet_JetType_NegativeE)      .at(iJet));
      jet->Set(JVF,            Get<vector<float> >(bjet_JetType_jvtxf)          .at(iJet));  
      jet->Set(TileGap1,       Get<vector<float> >(bjet_JetType_e_TileGap1)     .at(iJet) / 1000.);  
      jet->Set(TileGap2,       Get<vector<float> >(bjet_JetType_e_TileGap2)     .at(iJet) / 1000.);  
      jet->Set(TileGap3,       Get<vector<float> >(bjet_JetType_e_TileGap3)     .at(iJet) / 1000.);
      jet->Set(NumTowers,      Get<vector<float> >(bjet_JetType_NumTowers)      .at(iJet));

      jet->Set(ActiveAreaPx, Get<vector<float> >(bjet_JetType_ActiveAreaPx).at(iJet));
      jet->Set(ActiveAreaPy, Get<vector<float> >(bjet_JetType_ActiveAreaPy).at(iJet));
      jet->Set(ActiveAreaPz, Get<vector<float> >(bjet_JetType_ActiveAreaPz).at(iJet));
      jet->Set(ActiveAreaE, Get<vector<float> >(bjet_JetType_ActiveAreaE).at(iJet));
      jet->Set(ActiveArea,   Get<vector<float> >(bjet_JetType_ActiveArea).at(iJet));
      
      jet->Set(isBadLoose,     Get<vector<int> >  (bjet_JetType_isBadLoose)     .at(iJet));
      jet->Set(isBadLooseMinus,Get<vector<int> >  (bjet_JetType_isBadLooseMinus).at(iJet));
      
      // vector JVF
      if(doVectorJVFLinks){
          vector<float> jvtxfFull = Get<vector<vector<float> > >(bjet_JetType_jvtxfFull).at(iJet);
          jet->AddVec("JVFLinks");
          for (unsigned int iJVF=0; iJVF<jvtxfFull.size()-1; ++iJVF){ // last vectorJVF entry is for dummy vertex, which is not stored in fEvt->vtxs
                MomentObj* myLink = new MomentObj();
                myLink->Set("JVF", jvtxfFull[iJVF]);
                myLink->AddVec("vertex");
                if( iJVF < fEvt->vtxs()){  // fixing a very rare problem... 
                myLink->Add("vertex", &(fEvt->vtx(iJVF)));
                }
                jet->Add("JVFLinks", myLink);
          }
//          jet->Set(sumPtTrk_pv0_500MeV, Get<vector<float> > (bjet_JetType_sumPtTrk_pv0_500MeV).at(iJet));
      }

      // TruthMFLinks
      if(doTruthMFLinks){
        jet->AddVec(TruthMFLinks);

        MomentObj* myLink = new MomentObj();
        myLink->Set("TruthMFindex", Get<vector<float> >(bjet_JetType_TruthMFindex).at(iJet));
        myLink->Set("TruthMF",      Get<vector<float> >(bjet_JetType_TruthMF).at(iJet));    

        jet->Add(TruthMFLinks, myLink);
      }

      jet->Set(flavor_weight_MV1, Get<vector<float> >(bjet_JetType_flavor_weight_MV1).at(iJet));
      jet->Set(flavor_weight_SV1, Get<vector<float> >(bjet_JetType_flavor_weight_SV1).at(iJet));
      jet->Set(WIDTH, Get<vector<float> >(bjet_JetType_WIDTH).at(iJet));

      if(fEvt->isMC()){
        jet->Set(flavor_truth_label, Get<vector<int> >(bjet_JetType_flavor_truth_label).at(iJet));
      } else {
        jet->Set(flavor_truth_label, 0);
      }

      if(!doSMWZ){
        // Track association
        if(!Rename)
          jet->Set(nTrackAssoc, Get<vector<int> >(bjet_JetType_nTrackAssoc).at(iJet));
        else
          jet->Set("nconstits", Get<vector<int> >(bjet_JetType_nTrackAssoc).at(iJet));
        
        if(doTrack){
          vector<int> tracksIndices = Get<vector<vector<int> > >(bjet_JetType_TrackAssoc_index).at(iJet);
          
          if(!Rename) jet->AddVec(GhostAssocTrack);
          else        jet->AddVec("constit");
          
          double sumpt = 0.;

          for (unsigned int iTrack=0; iTrack<tracksIndices.size(); ++iTrack){
          
            if(!Rename) jet->Add(GhostAssocTrack, &(fEvt->track(tracksIndices[iTrack])));
            else        jet->Add("constit",       &(fEvt->track(tracksIndices[iTrack])));
            
            sumpt += fEvt->track(tracksIndices[iTrack]).p.Pt();
            
      	    // put a tag on tracks indicating which jet they belong to
	          // In case that a track is not associated with any jet, one should use 'Exists' to check.
       	    MomKey AssocJet_JetKey("AssocJet_");
	          AssocJet_JetKey += mjets; 
	          //fEvt->track(tracksIndices[iTrack]).AddVec(AssocJet_JetKey);
	          //fEvt->track(tracksIndices[iTrack]).Add(AssocJet_JetKey, jet); 
	          fEvt->track(tracksIndices[iTrack]).Set(AssocJet_JetKey, iJet); 
          }
        jet->Set(GhostAssocTrackSumPt, sumpt);
        }
      } // end check doSMWZ

    }

    // Add associated constit information: for truth jet, we add truth particle; for recon jet, we add cluster information.
    if(!doSMWZ) {

      jet->Set(nConstitAssoc, Get<vector<int> >(bjet_JetType_nConstitAssoc).at(iJet));
      if(doConstitLinks) {
        vector<int> constitIndices = Get<vector<vector<int> > >(bjet_JetType_ConstitAssoc_index).at(iJet);
        jet->AddVec(AssocConstit);
        double sumpt = 0.;
        for (unsigned int iConstit = 0; iConstit < constitIndices.size(); iConstit++) {
         if( (JetFlag == JET4TRUTH) || (JetFlag == JET6TRUTH) ) { 
           jet->Add(AssocConstit, &(fEvt->truth(constitIndices[iConstit], MomKey(""))));
           sumpt += fEvt->truth(constitIndices[iConstit], MomKey("")).p.Pt();
         }
         else { 
           jet->Add(AssocConstit, &(fEvt->cluster(constitIndices[iConstit], MomKey("LCTopo")))); 
           sumpt += fEvt->cluster(constitIndices[iConstit], MomKey("LCTopo")).p.Pt();
         }
       }
       jet->Set(AssocConstitSumPt, sumpt);
     }
   }    

   fEvt->Add(mjets, jet);
  } // end loop over jets

  // JES outside of the loop because of the need for fCloseby calculation, stupid
  if((whichsyste==JetMETSyst::JESUP || whichsyste ==JetMETSyst::JESDOWN) && fEvt->isMC() && !(JetFlag==JET4TRUTH || JetFlag==JET6TRUTH)){
    for(unsigned int iJet = 0; iJet<calJets.size(); ++iJet) {

      double fCloseby=0;
      for (unsigned int jJet=0; jJet<calJets.size(); ++jJet) {
        if (jJet==iJet) continue; // ignore probe jet
        if (calJets.at(jJet).Pt()<=12.) continue; // ignore jets with pT<=12 GeV
        if (calJets.at(jJet).DeltaR(calJets.at(iJet))>=1.1) continue; // only consider jets with DR<1.1
        // "fractional projection" of momentum along jet axis: p_jet*p_closeby / |p_jet|^2
        fCloseby += calJets.at(iJet).Vect().Dot(calJets.at(jJet).Vect()) / pow(calJets.at(iJet).P(),2) ; // close-by fraction
      }
      closeBy.push_back(fCloseby);
    }

    MultijetJESUncertaintyProvider* m_myJes = doAFII ? m_myJesAFII : m_myJesNorm;

    for(int iJet = 0; iJet < fEvt->jets(mjetsNoJ); iJet++){

      Particle* jet = &(fEvt->jet(iJet, mjetsNoJ));

      float fCloseby = closeBy[iJet];
      double uncert = 0.;
      double factor = 1.;

      float jetpt = 1000.* jet->p.Perp();
      float jetE  = 1000.* jet->p.E();
      float jeteta = jet->p.Eta();

      bool isBJet = Get<vector<int> >(bjet_JetType_flavor_truth_label).at(iJet) == 5;


      if(whichsyste==JetMETSyst::JESUP) {
        if(jetpt>20000 && fabs(jeteta)<4.5)
          uncert=m_myJes->getRelUncert(jetpt, jeteta,fCloseby,true,fEvt->vtxs(),fEvt->Float("averageIntPerXing"),isBJet);
      } else if(whichsyste==JetMETSyst::JESDOWN) {
        if(jetpt>20000 && fabs(jeteta)<4.5)
          uncert=m_myJes->getRelUncert(jetpt, jeteta,fCloseby,false,fEvt->vtxs(),fEvt->Float("averageIntPerXing"),isBJet);
        factor=-1;
      }

      jetpt=jetpt*(1+ (factor*uncert)) / 1000.;
      jetE=jetE*(1+ (factor*uncert)) / 1000.;
      
      jet_AntiKt4LCTopo_recalibrated_pt[iJet] = jetpt;
      jet->p.SetPtEtaPhiE(jetpt, jeteta, jet->p.Phi(), jetE);
    }
  }

  closeBy.clear();
  calJets.clear();
  if(fEvt->Debug()) cout << endl;

  return kTRUE;
} // end copy standard jets

// parse the non-standard jet strings
Bool_t EventBuilder_jetmet2012::ParseJetStrings(){
  if(!fEvt->Cfg()->Exists(mJETTYPEVEC)) {
    fEvt->Cfg()->AddVec(mJETTYPEVEC);
    //fEvt->AddVec(mJETTYPEVEC);
    TObjArray* arr = fEvt->Cfg()->String("JETTYPES").Tokenize(",");
    for(Int_t i = 0; i<arr->GetEntries(); ++i) {
      if (fEvt->Debug()) {
        TObjString* jetTypeObj = (TObjString*)arr->At(i);
    TString jetType        = jetTypeObj->GetString();
        cout << "EventBuilder_jetmet2012::CopyEvent(): DEBUG Adding jet type: " << jetType << endl; 
      }
      
      fEvt->Cfg()->Add(mJETTYPEVEC, (TObjString*)arr->At(i));
    }
    
    arr->SetOwner(kFALSE);  delete arr;
  }

  if(!fEvt->Cfg()->Exists(mBTAGVEC)) {
    fEvt->Cfg()->AddVec(mBTAGVEC);
    //fEvt->AddVec(mJETTYPEVEC);
    TObjArray* arr = fEvt->Cfg()->String("BTAGS").Tokenize(",");
    for(Int_t i = 0; i<arr->GetEntries(); ++i) {
      if (fEvt->Debug()) {
        TObjString* jetTypeObj = (TObjString*)arr->At(i);
    TString jetType        = jetTypeObj->GetString();
        cout << "EventBuilder_jetmet2012::CopyEvent(): DEBUG Adding btag type: " << jetType << endl; 
      }
      
      fEvt->Cfg()->Add(mBTAGVEC, (TObjString*)arr->At(i));
    }
    
    arr->SetOwner(kFALSE);  delete arr;
  }
  return kTRUE;
}

Bool_t EventBuilder_jetmet2012::AddPRWWeights(){
  if(doPRW){
    if(fEvt->Debug()) cout << " adding prw weights! " << endl;
    map<MomKey, Root::TPileupReweighting*>::iterator it;
    static const MomKey RunNumber("RunNumber"), ChannelNumber("ChannelNumber"), averageIntPerXing("averageIntPerXing");
    for(it = prwMap.begin(); it!=prwMap.end(); ++it){
      
      if(fEvt->isMC()){
        fEvt->Set(it->first, it->second->GetCombinedWeight(fEvt->Int(RunNumber),fEvt->Int(ChannelNumber),fEvt->Float(averageIntPerXing)));
//        cout << it->first << " " <<  it->second->GetCombinedWeight(fEvt->Int(RunNumber),fEvt->Int(ChannelNumber),fEvt->Float(averageIntPerXing)) << endl;
      }
      else{
        fEvt->Set(it->first, 1.);
      }
    }
    if(fEvt->Debug()) cout << " added prw weights! " << endl;

  }
  return kTRUE;
}

Bool_t EventBuilder_jetmet2012::CopyJetTriggers(){
    /// load branches that are necessary for trigger (if not already)
  if(doJetTriggers){
    const static MomKey Passed("Passed");

    Load("trig_DB_SMK");
    Load("trig_DB_L1PSK");
    Load("trig_DB_HLTPSK");
    Load("trig_L1_TAV");
    Load("trig_EF_passedPhysics");
    Load("trig_L2_passedPhysics");
    Load("trig_L1_TBP");
    Load("trig_L1_TAP");
    Load("trig_L2_passedRaw");
    Load("trig_EF_passedRaw");
    Load("trig_L2_resurrected");
    Load("trig_EF_resurrected");
    Load("trig_L2_passedThrough");
    Load("trig_EF_passedThrough");

    m_trigdefs->GetEntry(Entry());
    map<MomKey, Root::TPileupReweighting*>::iterator it;
    for(it = prwMap.begin(); it!=prwMap.end(); ++it){
      if(fEvt->Debug()) cout << "adding trigger " << it->first << " with Passed infomration!" << endl; 
      fEvt->Set(it->first+Passed, m_trigdefs->IsPassed(it->first.Data(), TrigDefs::Physics));
    }
  }
  return kTRUE;
}

