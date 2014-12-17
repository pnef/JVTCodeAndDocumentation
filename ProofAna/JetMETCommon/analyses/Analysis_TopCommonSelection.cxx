/**************************************************************************
 **
 **   File:         Analysis_TopCommonSelection.h
 **
 **   Description:  See header
 **                 
 **   Authors:      M. Swiatlowski, B. Nachman
 **
 **************************************************************************/

#define Analysis_TopCommonSelection_cxx

#include "Analysis_TopCommonSelection.h"
#include "AnaConfig.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "TKey.h"
#include "TObjString.h"
//#include "TObjArray.h"
#include "MuonEfficiencyCorrections/AnalysisMuonConfigurableScaleFactors.h"
#include "TrigMuonEfficiency/LeptonTriggerSF.h"
#include "TileTripReader/TTileTripReader.h" 



///=========================================
/// WorkerBegin: setup binnaing, etc
///=========================================
void Analysis_TopCommonSelection::WorkerBegin()
{
  if (Debug()) cout << "Analysis_TopCommonSelection: DEBUG In WorkerBegin()" << endl;

  Analysis_JetMET_Base::WorkerBegin();

  ChainCfg()->Get("DOSMWZ"        , doSMWZ       );
  ChainCfg()->Get("DOCOMMON"      , doCOMMON     );
  ChainCfg()->Get("SYSTTYPE"      , SystType     );
  doNOJVFcut = 0;
  Cfg()->Get("DONOJVFCUT"         , doNOJVFcut   );    


  BKey = GetBTagKey(SystType);

  MuonSFKey = GetMuonSFKey(SystType);

  MuonTrigSFKey = GetMuonTrigSFKey(SystType);

  //SetMaps(&m_h1d, &m_h2d);

  if (Debug()) cout << "Analysis_TopCommonSelection: DEBUG Finish WorkerBegin()" << endl;

  //MuonTriggerSF = new LeptonTriggerSF(2012,"../utils/SUSYTools/data/", "muon_trigger_sf_2012_AtoL.p1328.root", "../utils/ElectronEfficiencyCorrection/data/","rel17p2.v02");


  // std::string unit("GeV");
  // std::string muon_file_name("Muid_CB_plus_ST_2012_SF.txt");
  // std::string muon_sf_dir("../utils/MuonEfficiencyCorrections/share/"); /// Default path
  // Analysis::AnalysisMuonConfigurableScaleFactors::Configuration muon_configuration=Analysis::AnalysisMuonConfigurableScaleFactors::AverageOverRuns;
  // MuonSF = new Analysis::AnalysisMuonConfigurableScaleFactors(muon_sf_dir, muon_file_name, unit, muon_configuration);
  // MuonSF->Initialise();
    

  //MuonSF = (Analysis::AnalysisMuonConfigurableScaleFactors*) InputList()->FindObject("MuonSF");


    std::string unit("GeV");
    std::string muon_file_name("Muid_CB_plus_ST_2012_SF.txt");
    std::string maindir(gSystem->Getenv("PROOFANADIR"));
    if(maindir==".") maindir+="/libProofAna";
    std::string muon_sf_dir("/utils/MuonEfficiencyCorrections/share/"); /// Default path
    Analysis::AnalysisMuonConfigurableScaleFactors::Configuration muon_configuration=Analysis::AnalysisMuonConfigurableScaleFactors::AverageOverRuns;
    MuonSF = new Analysis::AnalysisMuonConfigurableScaleFactors(maindir+muon_sf_dir, muon_file_name, unit, muon_configuration);
    MuonSF->Initialise();

// * Pt                 | double       | 26.695848                   *
// * Eta                | double       | 2.460024                    *
// * Phi                | double       | 2.124920                    *
// * M                  | double       | 158.228891                  *
// *******************************************************************

   //TLorentzVector test(0.,0.,0.,0.);
   //float sf = 1.;
   //test.SetPtEtaPhiM(26.695848,2.460024,2.124920,158.228891);
   //float charge = -1.;
  //sf *= MuonSF->scaleFactor(charge, test);

  
    //cout << "MuonSF Init! " << MuonSF << " test sf = " << sf <<endl;


  //MuonTriggerSF = (LeptonTriggerSF*) InputList()->FindObject("MuonTriggerSF");
  MuonTriggerSF = new LeptonTriggerSF(2012,maindir+"/utils/SUSYTools/data/", "muon_trigger_sf_2012_AtoL.p1328.root", maindir+"/utils/ElectronEfficiencyCorrection/data/","rel17p2.v02");

  vector<TLorentzVector> dummy;
  TLorentzVector vec;
  vec.SetPtEtaPhiE(30.,0.5,0.5,50.);
  dummy.push_back(vec);
  MuonTriggerSF->GetTriggerSF(200804,true,dummy,combined,dummy,tightpp,0);

  
  myDataPeriod = (DataPeriod*) InputList()->FindObject("myDataPeriod");


  myBTagCalib = 0;
   myBTagCalib70 = 0;
   myBTagCalib80 = 0;
   myBTagCalib90 = 0;
   myBTagCalibLead = 0;
   myBTagCalibLead3 = 0;
   myJVFBTagCalib = 0;
   myJVFBTagCalib70 = 0;
   myJVFBTagCalib80 = 0;
   myJVFBTagCalib90 = 0;
   myJVFBTagCalibLead = 0;
   myJVFBTagCalibLead3 = 0;


  m_treader=new Root::TTileTripReader("myTripReader");

} 

///=========================================
/// WorkerTerminate: Usually nothing
///=========================================
void Analysis_TopCommonSelection::WorkerTerminate(){

}

///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_TopCommonSelection::ProcessEvent()
{
  if (Debug()) cout << "Analysis_TopCommonSelection: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;
  static const MomKey TopCommonJets("AntiKt4LCTopo");

  OutputDir()->cd();



  //MakeTruthPlots(); 
//-----------------------------------//
//              MUONS                //
//-----------------------------------//


  //if(Debug()){
  //  Show();
  //}

  AddVec("muonsgood");

  for (int i=0; i<muons(); i++){
   //Quality Criteria from https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TopCommonObjects#Muons
   //muon(i).Show();
   if (muon(i).Int("author")!=12) continue;
   if (muon(i).Int("tight")!=1) continue;
   if (muon(i).p.Pt() < 25) continue;
   if (fabs(muon(i).p.Eta())>2.5) continue;
   if (muon(i).Int("nPixHits")+muon(i).Int("nPixelDeadSensors") <= 0) continue;
   if (muon(i).Int("nSCTHits")+muon(i).Int("nSCTDeadSensors")<4) continue;
   if (muon(i).Int("nPixHoles")+muon(i).Int("nSCTHoles")>2) continue;
   if (!(muon(i).Int("expectBLayerHit") || muon(i).Int("nBLHits")) > 0) continue;
   if (fabs(muon(i).Float("z0wrtPV"))>2) continue;
   int n = muon(i).Int("nTRTOutliers")+muon(i).Int("nTRTHits");
   if (fabs(muon(i).p.Eta())<1.9 && fabs(muon(i).p.Eta())>0.1){
     if (n<=5) continue;
     if (muon(i).Int("nTRTOutliers")/n >=0.9) continue;
   }
   else if (n > 5){
     if (muon(i).Int("nTRTOutliers")/n >=0.9) continue;
   }
   //isolation cuts.
   if(muon(i).Float("Etcone20") > 4000) continue;
   if(muon(i).Float("Ptcone30") > 2500) continue; 
   bool isiso = true;
   for(int j=0; j<jets(TopCommonJets); j++){
     if (jet(j,TopCommonJets).p.Pt() < 25) continue;
     if (jet(j,TopCommonJets).Float("JVF") <0.25 && jet(j, TopCommonJets).p.Perp() < 50.) continue;
     if (jet(j,TopCommonJets).p.DeltaR(muon(i).p)<0.4){isiso=false;}
   }
   if (!isiso) continue;
   //if (!muon(i).Bool("TrigMatch") && i==0) continue; // not necessary, done in cut 5
   Add("muonsgood",&muon(i));
   //std::cout << muon(i).p.Pt() << std::endl;
 }

 if(Debug()) cout << "TopSelection: passed muons!" << endl;

  // Set up muon scale factors. Can only be done after good muons are defined.
  if(isMC()) 
    DoMuonSF();
  if(Debug()) cout << "TopSelection: setup muon sf!" << endl;


 //-----------------------------------//
 //                 ELECTRONS         //
 //-----------------------------------//

  for(int j=0; j<jets(TopCommonJets); j++){
    jet(j,TopCommonJets).Set("electron",0);
  }

 const static MomKey EGood("electronsgood");
 AddVec(EGood);

 for(int i=0; i < electrons(); i++){
   // From https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TopCommonObjects#Electrons
   if(!electron(i).Bool("tightPP")) continue;
   if(electron(i).Int("author")!=1 && electron(i).Int("author")!=3) continue;
   if(electron(i).Float("trackz0pvunbiased") >= 2.) continue;
   if(electron(i).p.Et() <= 25.) continue;
   float eta_c = fabs(electron(i).Float("cl_eta"));
   if( eta_c >= 2.47 || (eta_c < 1.52 && eta_c > 1.37) ) continue;
   if(electron(i).Float("Etcone20") <= 6000) continue;
   if(electron(i).Float("Ptcone30") <= 6000) continue; 


   bool isiso = true;
   for(int j=0; j<jets(TopCommonJets); j++){
     if (jet(j,TopCommonJets).p.Pt() < 25) continue;
     if (jet(j,TopCommonJets).Float("JVF") <0.25 && jet(j,TopCommonJets).p.Perp() < 50.) continue;
     if (jet(j,TopCommonJets).p.DeltaR(electron(i).p)<0.4){isiso=false; jet(j,TopCommonJets).Set("electron",1);}
   }
   if(!electron(i).Int("OQ")&1446) continue;
   Add(EGood, &electron(i));
 }

//-----------------------------------//
//                 JETS              //
//-----------------------------------//

 const static MomKey TopBJets = TopCommonJets + "B";
 const static MomKey JTopBJets = "jets"+TopBJets;
 const static MomKey TopGoodJets = TopCommonJets + "Good";
 const static MomKey JTopGoodJets ="jets"+ TopCommonJets + "Good";

 AddVec(JTopBJets);

 for(int i=0; i<jets(TopCommonJets); i++){
   jet(i,TopCommonJets).Set("btag",0);
   if(jet(i,TopCommonJets).Float("flavor_weight_MV1") > 0.7892){
     if (jet(i,TopCommonJets).p.Pt() < 25.) continue;
     if (jet(i,TopCommonJets).Bool("isBadLooseMinus")) continue;
     if (jet(i,TopCommonJets).p.Eta() >= 2.5) continue;
     if (!doNOJVFcut && jet(i,TopCommonJets).Float("JVF") < 0.25 && jet(i,TopCommonJets).p.Perp() < 50.) continue; 
     jet(i,TopCommonJets).Set("btag",1);
     Add(JTopBJets,&jet(i,TopCommonJets));
   }
 }

 if(Debug()) cout << "TopSelection: passed btags!" << endl;

 const static MomKey TopWJets = TopCommonJets + "W";
 const static MomKey JTopWJets = "jets"+TopWJets;

 AddVec(JTopGoodJets);
 for(int iJ = 0; iJ < jets(TopCommonJets); iJ++){
   if(jet(iJ,TopCommonJets).Bool("electron")){continue;}
   if(jet(iJ,TopCommonJets).p.Perp() < 25.) continue;
   if(!doNOJVFcut && jet(iJ,TopCommonJets).Float("JVF") < 0.25 && jet(iJ,TopCommonJets).p.Perp() < 50.) continue;
   if(jet(iJ,TopCommonJets).Bool("isBadLooseMinus")) continue;
   if(jet(iJ,TopCommonJets).p.Eta() >= 2.5) continue;
   Add(JTopGoodJets, &jet(iJ, TopCommonJets));
 }

 AddVec(JTopWJets);

 double chisqmin = 9999999999.;
 int jet1 = 0;
 int jet2 = 0;

 for (int j=0; j<jets(TopGoodJets); j++){
   jet(j,TopGoodJets).Set("W",0);
   if (jet(j,TopGoodJets).Bool("electron")){continue;}
   if( (jet(j,TopGoodJets).p.Eta() > 2.1) || (jet(j,TopGoodJets).p.Eta() < -2.1)){
     continue;
   }
   if (jet(j,TopGoodJets).Bool("btag")) continue;
   for (int k=0; k<jets(TopGoodJets); k++){
     if (j==k) continue;
     if (jet(k,TopGoodJets).Bool("btag")) continue;
     TLorentzVector W = jet(j,TopGoodJets).p+jet(k,TopGoodJets).p;
     double chisq = pow(W.M()-80.,2);
     if (chisq < chisqmin){
       chisqmin = chisq;
       jet1 = j;
       jet2 = k;
     }
   }
 }

 if (chisqmin < 900.){
   Add(JTopWJets,&jet(jet1,TopGoodJets));
   jet(jet1,TopGoodJets).Set("W",1);
   Add(JTopWJets,&jet(jet2,TopGoodJets));
   jet(jet2,TopGoodJets).Set("W",1);
 }


 if(Debug()) cout << "TopSelection: passed w!" << endl;

 if(isMC())
   DoBtaggingSF(true, TopGoodJets);

 if(Debug()){ cout << " Showing before selection " << endl; Show(); }

//-----------------------------------//
//             Selection             //
//-----------------------------------//

 BookCutflowNoPrefix("TopSelection_PreScale");
 //Weight(1.);
  if(isMC()){
    if(!Exists(BKey))
      Set(BKey,1.);
    if(!Exists(MuonSFKey))
      Set(MuonSFKey,1.);
    if(!Exists(MuonTrigSFKey))
      Set(MuonTrigSFKey,1.);
    if(!Exists("Full_lumi"))
      Set("Full_lumi",1.);
    Weight(DefaultWeight() * Float(MuonSFKey) * Float(MuonTrigSFKey) * Float(BKey) * Float("Full_lumi"));
    //cout << Float("Full_lumi") << endl;
  }

  BookCutflowNoPrefix("TopSelection_PostScale");

 bool all = true;

 Fill("TopSelection_MuPreCuts", Float("averageIntPerXing"), Weight(), 50., -0.5, 49.5);

  // -1: GRL
  if(!isMC() && !Bool("grl")){
    Set("TopSelection_GRL", false);
    all = false;
  } else {
    Set("TopSelection_GRL", true);
    if(all) BookCutflowNoPrefix("TopSelection_GRL");
  }

 // 0. lar error
 if(Int("larerror") > 1 || Int("tileerror") == 2 || (Int("coreflags")&0x40000) !=0 || m_treader->checkEvent(RunNumber(),Float("LBN"),EventNumber()) == 0) {
   Set("TopSelection_LarError", false);
   all = false;
 } else {
   Set("TopSelection_LarError", true);
   if(all) BookCutflowNoPrefix("TopSelection_LarError");
 }

 // 1. Good vertex
 if(vtxs() < 1){
   Set("TopSelection_Vtxs", false);
   all = false;
 } else {
   Set("TopSelection_Vtxs", true);
   if(all) BookCutflowNoPrefix("TopSelection_Vtxs");
 }

 if(Debug()) cout << "Passed cut 1" << endl;

  // 1.5 trigger
 if(doSMWZ || doCOMMON){
   if (!Bool("1muontrigger")){
     Set("TopSelection_MuonTrigger", false);
     all = false;
   } else {
     Set("TopSelection_MuonTrigger", true);
     if(all) BookCutflowNoPrefix("TopSelection_MuonTrigger");
   }
 }
 else { // else doSMWZ
   Set("TopSelection_MuTrigMatch",true);
 }

 // 2. at least 1 muon
 if (muons("good")<1){
   Set("TopSelection_AtLeastOneMuon", false);
   all = false;
 } else {
   Set("TopSelection_AtLeastOneMuon", true);
   if(all) BookCutflowNoPrefix("TopSelection_AtLeastOneMuon");
 }

 if(Debug()) cout << "Passed cut 2" << endl;

 // 3. exactly 1 muon
 if (muons("good")>1){
   Set("TopSelection_ExactlyOneMuon", false);
   all = false;
 } else {
   Set("TopSelection_ExactlyOneMuon", true);
   if(all) BookCutflowNoPrefix("TopSelection_ExactlyOneMuon");
 }

 if(Debug()) cout << "Passed cut 3" << endl;

 // 4. no good electrons
 if (electrons("good") > 0){
   Set("TopSelection_NoElectrons", false);
   all = false;
 } else {
  Set("TopSelection_NoElectrons", true);
  if(all) BookCutflowNoPrefix("TopSelection_NoElectrons");
 }

 if(Debug()) cout << "Passed cut 4" << endl;

 // 5. triggermatch
 if(doSMWZ || doCOMMON){
   if(muons("good") > 0 && Bool("1muontrigger")){
     if (muon(0,"good").Exists("TrigMatch")){
       if (muon(0,"good").Bool("TrigMatch")) {
         Set("TopSelection_MuTrigMatch", true);
         if(all) BookCutflowNoPrefix("TopSelection_MuTrigMatch");
       } else {
         all = false;
         Set("TopSelection_MuTrigMatch", false);
       }
     } else {
       all = false;
       Set("TopSelection_MuTrigMatch", false);
     }
   }  
   else { // else good muons and trigger  passed
     all = false;
     Set("TopSelection_MuTrigMatch",false);
   }
 } else { // else doSMWZ
   Set("TopSelection_MuTrigMatch",true);
 }

 if(Debug()) cout << "Passed cut 5" << endl;

 // 6. mu/e overlap removal
 bool phiOverlap = false, etaOverlap = false;
 for(int iE = 0; iE < electrons("good"); iE++){
   if(muons("good") > 0){
     if(fabs(electron(iE,"good").Float("trk_phi") - muon(0,"good").Float("id_phi")) < 0.005){
       phiOverlap = true;
     }
     if(fabs(electron(iE,"good").Float("trk_theta") - muon(0,"good").Float("id_theta")) < 0.005){
       etaOverlap = true;
     }
   }
 }
 if(phiOverlap && etaOverlap){
   all = false;
   Set("TopSelection_MuEOverlap", false);
 } else {
   Set("TopSelection_MuEOverlap", true);
   if(all) BookCutflowNoPrefix("TopSelection_MuEOverlap");
 }

 if(Debug()) cout << "Passed cut 6" << endl;

 // 7. reject if loose bad jet
 bool hasBadLoose = false;
 for(int iJ = 0; iJ < jets(TopCommonJets); iJ++){
   if(jet(iJ, TopCommonJets).Bool("isBadLoose")){
    hasBadLoose = true;
    break;
   }
 }
 if(hasBadLoose){
   all = false;
   Set("TopSelection_NoBadLooseJets", false);
 } else {
   Set("TopSelection_NoBadLooseJets", true);
   if(all) BookCutflowNoPrefix("TopSelection_NoBadLooseJets");
 }


 if(Debug()) cout << "Passed cut 7" << endl;

 // 8. at least 2 good jets
 if(jets(TopGoodJets) < 2){
   all = false;
   Set("TopSelection_2GoodJets", false);
 } else {
   Set("TopSelection_2GoodJets", true);
   if(all) BookCutflowNoPrefix("TopSelection_2GoodJets");
 }

 if(Debug()) cout << "Passed cut 8" << endl;

 // 9. at least 3 good jets
 if(jets(TopGoodJets) < 3){
   all = false;
   Set("TopSelection_3GoodJets", false);
 } else {
   Set("TopSelection_3GoodJets", true);
   if(all) BookCutflowNoPrefix("TopSelection_3GoodJets");
 }

 if(Debug()) cout << "Passed cut 9" << endl;

 // 10. at least 4 good jets
 if(jets(TopGoodJets) < 4){
   all = false;
   Set("TopSelection_4GoodJets", false);
 } else {
   Set("TopSelection_4GoodJets", true);
   if(all) BookCutflowNoPrefix("TopSelection_4GoodJets");
 }

 if(Debug()) cout << "Passed cut 10" << endl;

 ////////////////////////////////////
 ///  MET SETUP
 ////////////////////////////////////

 float MET_et  = Float("MET_RefFinal_et");
 float MET_phi = Float("MET_RefFinal_phi");
 float METx    = Float("MET_RefFinal_x");
 float METy    = Float("MET_RefFinal_y");

 // float MET_et = 0.; Float("NewMET");
 // float MET_phi = 0.; Float("NewMETphi");
 // float METx = 0.; Float("NewMETx");
 // float METy = 0.; Float("NewMETy");

 if(doSMWZ || doCOMMON){
   MET_et = Float("NewMET");
   MET_phi = Float("NewMETphi");
   float METx = Float("NewMETx");
   float METy = Float("NewMETy");
 }

 double MT = 0.;
 if(muons("good") > 0){
  double MTsquared = 2.*(MET_et*muon(0,"good").p.Pt()-METx*muon(0,"good").p.Px()-METy*muon(0,"good").p.Py());
  MT = sqrt(std::max(0.,MTsquared));
 }
 Set("MT",MT);


 // 11. met is > 20
 if(MET_et <= 20.){
   all = false;
   Set("TopSelection_MET20", false);
 } else {
   Set("TopSelection_MET20", true);
   if(all) BookCutflowNoPrefix("TopSelection_MET20");
 }


 if(Debug()) cout << "Passed cut 11" << endl;

 // 12. MET + MT > 60
 if(MET_et + MT <= 60.){
   all = false;
   Set("TopSelection_METMT60", false);
 } else {
   Set("TopSelection_METMT60", true);
   if(all) BookCutflowNoPrefix("TopSelection_METMT60");
 }
 
 if(Debug()) cout << "Passed cut 12" << endl;

 // 13. at least one b jet
 if(jets(TopBJets) < 1){
   all = false;
   Set("TopSelection_1BJet", false);
 } else {
   Set("TopSelection_1BJet", true);
   if(all) BookCutflowNoPrefix("TopSelection_1BJet");
 }

 Set("TopSelection_AllStandard", all);

 if(Debug()) cout << "Passed cut 13" << endl;

 Set("TopSelection_AllStandard", all);

 if(all){
  Fill("TopSelection_MuPostCuts", Float("averageIntPerXing"), Weight(), 50., -0.5, 49.5);
 }

 // 14. OUR OWN at least two b jets
 if(jets(TopBJets) < 2){
   all = false;
   Set("TopSelection_2BJet", false);
 } else {
   Set("TopSelection_2BJet", true);
   if(all) BookCutflowNoPrefix("TopSelection_2BJet");
 }

 Set("TopSelection_All2B",all);

 if(Debug()) cout << "Passed cut 14" << endl;
 // 15. OUR OWN 2 other jets in W mass window
 if(jets(TopWJets) < 2){
   all = false;
   Set("TopSelection_2WJets", false);
 } else {
   Set("TopSelection_2WJets", true);
   if(all) BookCutflowNoPrefix("TopSelection_2WJets");
 }
 if(Debug()) cout << "Passed cut 15" << endl;

 Set("TopSelection_All", all);

 if(Debug()){
  cout << "showing the full event now " << endl;
  Show();
 }
  

 // THIS IS ALL NOW DEFUNCT


 if (jets(TopBJets) <  2){
   Set("TopSelection_BJets", false);
   Set("TopSelection_BJetsJVF", false);
   all = false;
 } else {
   Set("TopSelection_BJets", true);
   if (jet(0,TopBJets).Float("JVF") < 0.5 || jet(1,TopBJets).Float("JVF") < 0.5){
     Set("TopSelection_BJetsJVF", false);
     all = false;
   } else {
     Set("TopSelection_BJetsJVF", true);
   }
 }

 if(Debug()) cout << "TopSelection: b-jet selection finished" << endl;

 if (jets(TopWJets) != 2){
   Set("TopSelection_WJets", false);
   Set("TopSelection_WJetsJVF", false);
   Set("TopSelection_WJetsPt", false);
   all = false;
 } else {
   Set("TopSelection_WJets", true);
   if (jet(0,TopWJets).Float("JVF") < 0.25 || jet(1,TopWJets).Float("JVF") < 0.25){
     Set("TopSelection_WJetsJVF", false);
     all = false;
   } else {
     Set("TopSelection_WJetsJVF", true);
   }

   if (jet(0,TopWJets).p.Pt() < 25. || jet(1,TopWJets).p.Pt() < 25.){
     Set("TopSelection_WJetsPt", false);
     all = false;
   } else {
     Set("TopSelection_WJetsPt", true);
   }
 }
 if(Debug()) cout << "TopSelection: w-jet selection finished" << endl;


 // END DEFUNCT SECTION


 if(Debug()) cout << "TopSelection: passed all selection!" << endl;

 return true;
}

