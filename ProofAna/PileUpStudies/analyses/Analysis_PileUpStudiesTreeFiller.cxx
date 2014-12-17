/**************************************************************************
 **
 **   File:         Analysis_PileUpStudiesTreeFiller.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      P. Nef
 **
 **************************************************************************/

#define Analysis_PileUpStudiesTreeFiller_cxx

#include "Analysis_PileUpStudiesTreeFiller.h"
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
#include "TMVA/Reader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"


///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_PileUpStudiesTreeFiller::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_PileUpStudiesTreeFiller: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();

  fTreeName = "TestTree";
  Cfg()->Get("TreeName", fTreeName);
  fTrkSel   = "tracksGoodSel0";
  Cfg()->Get("TrackSel", fTrkSel);

  // config
  fAddTacksToTree    = false;
  ChainCfg()->Get("ADDTRACKSTOTREE"   , fAddTacksToTree);
  fAddLeptonZInfoToTree = true;
  ChainCfg()->Get("ADDLEPTONZINFOTOTREE", fAddLeptonZInfoToTree);
  fisMuScan          = false;
  ChainCfg()->Get("isMuScan",           fisMuScan);

  // TMVA reader ----------------------------------------------------------------------------
  ChainCfg()->Get("MVA1name"        , fMVA1name );
  ChainCfg()->Get("MVA2name"        , fMVA2name );
  ChainCfg()->Get("MVA3name"        , fMVA3name );
  ChainCfg()->Get("MVA4name"        , fMVA4name );
  ChainCfg()->Get("MVA1file"        , fMVA1file );
  ChainCfg()->Get("MVA2file"        , fMVA2file );
  ChainCfg()->Get("MVA3file"        , fMVA3file );
  ChainCfg()->Get("MVA4file"        , fMVA4file );
  ChainCfg()->Get("doMVAEval"      , doMVAEval );
  ChainCfg()->Get("runCluster"     ,  runCluster);

  if (doMVAEval){
      cout << "-----------> TMVA ------------------" << endl;
      cout << "Analysis_PileUpStudies: Reading MVA " << endl
           <<  fMVA1name << " with input file " << fMVA1file << endl
           <<  fMVA2name << " with input file " << fMVA2file << endl
           <<  fMVA3name << " with input file " << fMVA3file << endl
           <<  fMVA4name << " with input file " << fMVA4file << endl;
      
      fTMVAreader = new TMVA::Reader( "!Color:!Silent");
      fTMVAreader->AddVariable( "JnPUTrkCorrJVF:=((JHStrkPtSum+JPUtrkPtSum>0)?JHStrkPtSum/(JHStrkPtSum+ JPUtrkPtSum/(TMath::Max(nPUTracks,1)*0.01)):-1)",  &fTJnPUTrkCorrJVF   );
      fTMVAreader->AddVariable( "JHSPVtrkSumOverPt",                                                                                          &fTJHSPVtrkSumOverPt);
      fTMVAreader->AddSpectator( "NPV",                                                                                                       &fTNPV);
      fTMVAreader->AddSpectator( "Mu",                                                                                                        &fTMu);
      fTMVAreader->AddSpectator( "Jpt",                                                                                                       &fTJpt);
      fTMVAreader->AddSpectator( "Jeta",                                                                                                      &fTJeta);
      fTMVAreader->AddSpectator( "JisPU",                                                                                                     &fTJisPU);
      fTMVAreader->AddSpectator( "JisHS",                                                                                                     &fTJisHS);
      fTMVAreader->AddSpectator( "Weight",                                                                                                    &fTWeight);
      fTMVAreader->AddSpectator( "VtxDzTruth",                                                                                                &fTVtxDzTruth);
      fTMVAreader->AddSpectator( "Jtruthpt",                                                                                                  &fTJtruthpt);
      fTMVAreader->AddSpectator( "JJVF",                                                                                                      &fTJJVF);
      fTMVAreader->AddSpectator( "JCorrPUtrkPtSumOverPt",                                                                                     &fTJCorr2PUtrkPtSumOverPt);
      fTMVAreader->AddSpectator( "JMuCorrJVF",                                                                                                &fTJMuCorrJVF);
      fTMVAreader->AddSpectator( "JNPVCorrJVF",                                                                                               &fTJNPVCorrJVF);
      fTMVAreader->AddSpectator( "nPUTracks",                                                                                                 &fTnPUTracks);
      fTMVAreader->AddSpectator( "JHStrkPtSum",                                                                                               &fTJHStrkPtSum);

      // Book MVA
      TString path;
      if(runCluster) path = "/tmp/pnef/";
      else           path = "/u/at/pnef/nfs/xml-files/PUID/"; 
      cout << TString::Format("%s%s",path.Data(),fMVA1file.Data()).Data() << endl;
      cout << TString::Format("%s%s",path.Data(),fMVA2file.Data()).Data() << endl;
      cout << TString::Format("%s%s",path.Data(),fMVA3file.Data()).Data() << endl;
      cout << TString::Format("%s%s",path.Data(),fMVA4file.Data()).Data() << endl;
      if(fMVA1file.Length()>0) fTMVAreader->BookMVA( fMVA1name, TString::Format("%s%s",path.Data(),fMVA1file.Data()).Data());
      if(fMVA2file.Length()>0) fTMVAreader->BookMVA( fMVA2name, TString::Format("%s%s",path.Data(),fMVA2file.Data()).Data());
      if(fMVA3file.Length()>0) fTMVAreader->BookMVA( fMVA3name, TString::Format("%s%s",path.Data(),fMVA3file.Data()).Data());
      if(fMVA4file.Length()>0) fTMVAreader->BookMVA( fMVA4name, TString::Format("%s%s",path.Data(),fMVA4file.Data()).Data());
  }
  
  // trees ------------------
  OutputDir()->cd();
  fJetTree = new TTree(fTreeName, "Tree with PU and HS jets");
  AddBranches(fJetTree);
  
} 


///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_PileUpStudiesTreeFiller::ProcessEvent()
{

  if (Debug()) cout << "Analysis_PileUpStudiesTreeFiller: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();
  
  if(! Exists("PileUpStudiesSelection")) return true;
  if(! Bool("PileUpStudiesSelection")  ) return true;

  // Trk definition ------------------------------------------------------------------------
  MomKey myTracksGhost(fTrkSel+"Ghost");
  MomKey myJVFLinks("JVFLinks_"+fTrkSel+"Ghost");

  // Fill Tree-----------------------------------------------------------------------------
  FillTree("AntiKt4LCTopo_ptgt20", myJVFLinks, myTracksGhost);


  // end----------------------------------------------------------------------
  if (Debug()) cout << "Analysis_PileUpStudiesTreeFiller: DEBUG End ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;
  return true;
}



///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_PileUpStudiesTreeFiller::WorkerTerminate()
{
  gDirectory->cd("/");
  fJetTree->Write();


  // Nothing more

}

///===========================================
/// Fill Tree
///========================================
void Analysis_PileUpStudiesTreeFiller::FillTree(const MomKey JetKey, const MomKey JVFKey, const MomKey trackKey){
  if(Debug()) cout <<"Analysis_PileUpStudiesTreeFiller::FillTree Start" << endl;
  ResetBranches(fJetTree);
  FillEventVars(fJetTree, JetKey, trackKey);
  fJetTree->Fill();
  for(int iJet=0; iJet<jets(JetKey); ++iJet){
    if(fisMuScan && !(jet(iJet,JetKey).Float("constscale_pt")>20 || jet(iJet,JetKey).Float("areacorr_pt")>20 || jet(iJet,JetKey).p.Pt()>20 )) continue; // cut for muscan samples
    ResetBranches(fJetTree);
    FillEventVars(fJetTree, JetKey, trackKey);
    FillJetVars(fJetTree, iJet, JetKey, JVFKey, trackKey);

    fJetTree->Fill();
  }
    
  if(Debug()) cout <<"Analysis_PileUpStudiesTreeFiller::FillTree End" << endl;
 return;

}

///=============================================
/// Add Branches To Tree
///=============================================
void Analysis_PileUpStudiesTreeFiller::AddBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_PileUpStudiesTreeFiller::AddBranches Start" << endl;
  
    // Event Info
    tree->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
    tree->Branch("RunNumber",                 &fTRunNumber,              "RunNumber/I");
    tree->Branch("Weight" ,                   &fTWeight,                 "Weight/F");
    tree->Branch("PUWeight" ,                 &fTPUWeight,               "PUWeight/F");
    tree->Branch("LumiXSecWeight" ,           &fTLumiXSecWeight,         "LumiXSecWeight/F");
    tree->Branch("NPV" ,                      &fTNPV,                    "NPV/F");
    tree->Branch("NPVtruth" ,                 &fTNPVtruth,               "NPVtruth/F");
    tree->Branch("VtxDzTruth",                &fTVtxDzTruth ,            "VtxDzTruth/F");
    tree->Branch("HSPVmindZPUVtx",            &fTHSPV_mindZPUVtx,        "HSPVmindZPUVtx/F");
    tree->Branch("NMCVertWithin1mmFromHS",    &fTNMCVertWithin1mmFromHS, "NMCVertWithin1mmFromHS/I");
    tree->Branch("HSVtxTruthminDz",           &fTHSVtxTruthminDz,        "HSVtxTruthminDz/F");
    tree->Branch("Mu" ,                       &fTMu,                     "Mu/F");
    tree->Branch("GRL",                       &fTGRL,                    "GRL/I");
    tree->Branch("NBJets" ,                   &fTNBJets,                 "NBJets/I");
    tree->Branch("RhoGoodPUTracks",           &fTRhoGoodPUTracks,        "RhoGoodPUTracks/F");
    tree->Branch("RhoPUTracks",               &fTRhoPUTracks,            "RhoPUTracks/F");
    tree->Branch("RhoKt4LC",                  &fTRhoKt4LC,               "RhoKt4LC/F");
    tree->Branch("nPUTracks",                 &fTnPUTracks,              "nPUTracks/I");
    tree->Branch("nTracksNotFromPV",          &fTnTracksNotFromPV,       "nTracksNotFromPV/I");
    tree->Branch("m1MuTrigScaleFactor" ,      &fTm1MuTrigScaleFactor,    "m1MuTrigScaleFactor/F");
    tree->Branch("m1MuScaleFactor" ,          &fTm1MuScaleFactor,        "m1MuScaleFactor/F");

    // Zboson
    if(! fisMuScan && fAddLeptonZInfoToTree){
    tree->Branch("ZPt",                       &fTZpt,                         "ZPt/F");
    tree->Branch("ZEta",                      &fTZeta,                        "ZEta/F");
    tree->Branch("ZPhi",                      &fTZphi,                        "ZPhi/F");
    tree->Branch("ZMass",                     &fTZmass,                       "ZMass/F");
    }
    
    // Muons
    if(! fisMuScan && fAddLeptonZInfoToTree){
    tree->Branch("NMuons",                   &fTNMuons,                      "NMuons/I");   
    tree->Branch("MuonPt",                   &fTMuonPt,                      "MuonPt[NMuons]/F");   
    tree->Branch("MuonPhi",                  &fTMuonPhi,                     "MuonPhi[NMuons]/F");   
    tree->Branch("MuonEta",                  &fTMuonEta,                     "MuonEta[NMuons]/F");   
    tree->Branch("MuonPtUnsmeared",          &fTMuonPtUnsmeared,             "MuonPtUnsmeared[NMuons]/F");   
    tree->Branch("MuonPhiUnsmeared",         &fTMuonPhiUnsmeared,            "MuonPhiUnsmeared[NMuons]/F");   
    tree->Branch("MuonEtaUnsmeared",         &fTMuonEtaUnsmeared,            "MuonEtaUnsmeared[NMuons]/F");   
    }

    // Jet multiplicity
    if(! fisMuScan){
    tree->Branch("NJetspt20eta2p1",          &fTNJetspt20eta2p1,             "NJetspt20eta2p1/I");   
    tree->Branch("NJetspt30eta2p1",          &fTNJetspt30eta2p1,             "NJetspt30eta2p1/I");   
    tree->Branch("NJetspt40eta2p1",          &fTNJetspt40eta2p1,             "NJetspt40eta2p1/I");   
    tree->Branch("NJetspt50eta2p1",          &fTNJetspt50eta2p1,             "NJetspt50eta2p1/I");   
    tree->Branch("NJetspt20eta2p1HS",        &fTNJetspt20eta2p1HS,           "NJetspt20eta2p1HS/I");   
    tree->Branch("NJetspt30eta2p1HS",        &fTNJetspt30eta2p1HS,           "NJetspt30eta2p1HS/I");   
    tree->Branch("NJetspt40eta2p1HS",        &fTNJetspt40eta2p1HS,           "NJetspt40eta2p1HS/I");   
    tree->Branch("NJetspt50eta2p1HS",        &fTNJetspt50eta2p1HS,           "NJetspt50eta2p1HS/I");   
    }


    //leading three jets
    tree->Branch("J0Pt",            &fTJ0Pt,             "J0Pt/F");
    tree->Branch("J1Pt",            &fTJ1Pt,             "J1Pt/F");
    tree->Branch("J2Pt",            &fTJ2Pt,             "J2Pt/F");
    tree->Branch("J0Eta",           &fTJ0Eta,            "J0Eta/F");
    tree->Branch("J1Eta",           &fTJ1Eta,            "J1Eta/F");
    tree->Branch("J2Eta",           &fTJ2Eta,            "J2Eta/F");
    tree->Branch("J0Phi",           &fTJ0Phi,            "J0Phi/F");
    tree->Branch("J1Phi",           &fTJ1Phi,            "J1Phi/F");
    tree->Branch("J2Phi",           &fTJ2Phi,            "J2Phi/F");
    tree->Branch("JLeadingTruthPt", &fTJLeadingTruthPt,  "JLeadingTruthPt/F");
    if(! fisMuScan && fAddLeptonZInfoToTree){
    tree->Branch("J0ZDPhi",         &fTJ0ZDPhi,          "J0ZDPhi/F");
    tree->Branch("J1ZDPhi",         &fTJ1ZDPhi,          "J1ZDPhi/F");
    tree->Branch("J2ZDPhi",         &fTJ2ZDPhi,          "J2ZDPhi/F");
    tree->Branch("J0MuonLeadDR",    &fTJ0MuonLeadDR,     "J0MuonLeadDR/F");
    tree->Branch("J1MuonLeadDR",    &fTJ1MuonLeadDR,     "J1MuonLeadDR/F");
    tree->Branch("J2MuonLeadDR",    &fTJ2MuonLeadDR,     "J2MuonLeadDR/F");
    tree->Branch("J0MuonTrailDR",   &fTJ0MuonTrailDR,    "J0MuonTrailDR/F");
    tree->Branch("J1MuonTrailDR",   &fTJ1MuonTrailDR,    "J1MuonTrailDR/F");
    tree->Branch("J2MuonTrailDR",   &fTJ2MuonTrailDR,    "J2MuonTrailDR/F");
    }

  
    // Jet vars --------------------------------------------------------------------------------------
    tree->Branch("Jindex",                    &fTJindex,                 "Jindex/I");
    tree->Branch("Jpt"  ,                     &fTJpt,                    "Jpt/F");
    tree->Branch("Jeta" ,                     &fTJeta,                   "Jeta/F");
    tree->Branch("Jphi" ,                     &fTJphi,                   "Jphi/F");
    tree->Branch("Jconstscalept"  ,           &fTJconstitpt,             "Jconstscalept/F");
    tree->Branch("Jareacorrpt"  ,             &fTJareacorrpt,            "Jareacorrpt/F");
    tree->Branch("Jtruthpt"  ,                &fTJtruthpt,               "Jtruthpt/F");
    tree->Branch("Jtrutheta" ,                &fTJtrutheta,              "Jtrutheta/F");
    tree->Branch("Jtruthphi" ,                &fTJtruthphi,              "Jtruthphi/F");
    tree->Branch("JtruthMatchDR" ,            &fTJtruthMatchDR,          "JtruthMatchDR/F");
    tree->Branch("JtruthInTimept"  ,          &fTJtruthInTimept,          "JtruthInTimept/F");
    tree->Branch("JtruthInTimeMatchDR" ,      &fTJtruthInTimeMatchDR,    "JtruthInTimeMatchDR/F");
    tree->Branch("JdRTruthPt4" ,              &fTJdRTruthPt4,            "dRTruthPt4/F");
    tree->Branch("JdRTruthPt6" ,              &fTJdRTruthPt6,            "dRTruthPt6/F");
    tree->Branch("JdRTruthPt8" ,              &fTJdRTruthPt8,            "dRTruthPt8/F");
    tree->Branch("JdRTruthPt10" ,             &fTJdRTruthPt10,           "dRTruthPt10/F");
    tree->Branch("JArea" ,                    &fTJArea,                  "JArea/F");
    tree->Branch("JisPU" ,                    &fTJisPU,                  "JisPU/I");
    tree->Branch("JisHS" ,                    &fTJisHS,                  "JisHS/I");
    tree->Branch("JTruthMF",                  &fTJTruthMF,               "JTruthMF/F");
    tree->Branch("JisBadLoose",               &fTJisBadLoose,            "JisBadLoose/I");
    tree->Branch("JisParton_uds",             &fTJisParton_uds,          "JisParton_uds/I");
    tree->Branch("JisParton_c",               &fTJisParton_c,            "JisParton_c/I");
    tree->Branch("JisParton_b",               &fTJisParton_b,            "JisParton_b/I");
    tree->Branch("JisParton_g",               &fTJisParton_g,            "JisParton_g/I");
    tree->Branch("JisBtagged",                &fTJisBtagged ,            "JisBtagged/I");
    tree->Branch("JNContribVtxs",             &fTJNContribVtxs,          "JNContribVtxs/I");
    tree->Branch("JClosebydR",                &fTClosebydR,              "JClosebydR/F");
    tree->Branch("JClosebykt",                &fTClosebykt,              "JClosebykt/F");
    tree->Branch("JClosebydRpT",              &fTClosebydRpT,            "JClosebydRpT/F");
    if(! fisMuScan){
    tree->Branch("JZDPhi",                    &fTJZDPhi,                 "JZDPhi/F");
    tree->Branch("JMuonLeadDR",               &fTJMuonLeadDR,            "JMuonLeadDR/F");
    tree->Branch("JMuonTrailDR",              &fTJMuonTrailDR,           "JMuonTrailDR/F");
    }

    if(! fisMuScan){
    tree->Branch("JJVF" ,                     &fTJJVF,                   "JJVF/F");
    tree->Branch("JstdJVF" ,                  &fTJstdJVF,                "JstdJVF/F");
    tree->Branch("JMaxJVF" ,                  &fTJMaxJVF,                "JMaxJVF/F");
    tree->Branch("JBetaStar",                 &fTJBetaStar               ,"JBetaStar/F");
    tree->Branch("JNPVCorrJVF",               &fTJNPVCorrJVF             , "JNPVCorrJVF/F");
    tree->Branch("JNContribVtxsCorrJVF",      &fTJNContribVtxsCorrJVF,     "JNContribVtxsCorrJVF/F");
    tree->Branch("JMuCorrJVF",                &fTJMuCorrJVF              , "JMuCorrJVF/F");
    tree->Branch("JMuCorrJVF0p1",             &fTJMuCorrJVF0p1           , "JMuCorrJVF0p1/F");
    tree->Branch("JMuCorrJVF0p14",            &fTJMuCorrJVF0p14          , "JMuCorrJVF0p14/F");
    tree->Branch("JnPUTrkCorrJVF",            &fTJnPUTrkCorrJVF,           "JnPUTrkCorrJVF/F");

    tree->Branch("JHStrkPtSum"    ,           &fTJHStrkPtSum,            "JHStrkPtSum/F");
    tree->Branch("JHStrkPtSumSquared",        &fTJHStrkPtSumSquared ,    "JHStrkPtSumSquared/F");
    tree->Branch("JHStrkPtSumdRWeight",       &fTJHStrkPtSumdRWeight ,   "JHStrkPtSumdRWeight/F");
    tree->Branch("JPUtrkPtSum",               &fTJPUtrkPtSum ,           "JPUtrkPtSum/F");
    tree->Branch("JPUtrkPtSumSquared",        &fTJPUtrkPtSumSquared ,    "JPUtrkPtSumSquared/F");
    tree->Branch("JPUtrkPtSumdRWeight",       &fTJPUtrkPtSumdRWeight ,   "JPUtrkPtSumdRWeight/F");
    tree->Branch("JPUtrkPtSumOverPt",         &fTJPUtrkPtSumOverPt,      "JPUtrkPtSumOverPt/F");
    tree->Branch("JCorrPUtrkPtSumOverPt",     &fTJCorrPUtrkPtSumOverPt,  "JCorrPUtrkPtSumOverPt/F");
    tree->Branch("JCorr2PUtrkPtSumOverPt",    &fTJCorr2PUtrkPtSumOverPt, "JCorr2PUtrkPtSumOverPt/F");
    tree->Branch("JHSPVtrkSumOverPt" ,        &fTJHSPVtrkSumOverPt,      "JHSPVtrkSumOverPt/F");
    tree->Branch("JHSPVtrkSumOverPtdRW" ,     &fTJHSPVtrkSumOverPtdRW,   "JHSPVtrkSumOverPtdRW/F");

    tree->Branch("JHSPVTrkWidth"     ,        &fTJHSPVTrkWidth,          "JHSPVTrkWidth/F");
    tree->Branch("JHSPVTrkMeanDR"    ,        &fTJHSPVTrkMeanDR,         "JHSPVTrkMeanDR/F");

    tree->Branch("JHSPVtrkPtlead" ,           &fTJHSPVtrkPtlead,         "JHSPVtrkPtlead/F");
    tree->Branch("JHSPVtrkPttrail" ,          &fTJHSPVtrkPttrail,        "JHSPVtrkPttrail/F");
    tree->Branch("JLeadTrkPt" ,               &fTJLeadTrkPt,             "JLeadTrkPt/F");
    tree->Branch("JHSPVTrk12dR" ,             &fTJHSPVTrk12dR,           "JHSPVTrk12dR/F");
    tree->Branch("JLeadTrkAllHSPVdZ",         &fTJLeadTrkAllHSPVdZ,      "JLeadTrkAllHSPVdZ/F");
    }
    
    tree->Branch("JNClusters",                &fTJNClusters,             "JNClusters/I");
    tree->Branch("JNumTowers",                &fTJNumTowers ,            "JNumTowers/F");
    tree->Branch("JCaloWidth",                &fTJCaloWidth,             "JCaloWidth/F");
    tree->Branch("JCaloWidth2",               &fTJCaloWidth2,            "JCaloWidth2/F");
    tree->Branch("JAnnulusPtRatio0to1",       &fTJAnnulusPtRatio0to1,    "JAnnulusPtRatio0to1/F");
    tree->Branch("JAnnulusPtRatio0to2",       &fTJAnnulusPtRatio0to2,    "JAnnulusPtRatio0to2/F");
    tree->Branch("JAnnulusPtRatio1to2",       &fTJAnnulusPtRatio1to2,    "JAnnulusPtRatio1to2/F");
    tree->Branch("JAnnulusPtRatio2to4",       &fTJAnnulusPtRatio2to4,    "JAnnulusPtRatio2to4/F");
    tree->Branch("JAverageDeltaRSquared",     &fTJAverageDeltaRSquared,  "JAverageDeltaRSquared/F");
    tree->Branch("JEMFrac",                   &fTJEMFrac,                "JEMFrac/F");
    tree->Branch("JTiming",                   &fTJTiming,                "JTiming/F");
    

    if(! fisMuScan){
    tree->Branch(TString::Format("J%s",fMVA1name.Data()).Data(), &fTJMVA1result,  TString::Format("J%s/F",fMVA1name.Data()).Data());
    tree->Branch(TString::Format("J%s",fMVA2name.Data()).Data(), &fTJMVA2result,  TString::Format("J%s/F",fMVA2name.Data()).Data());
    tree->Branch(TString::Format("J%s",fMVA3name.Data()).Data(), &fTJMVA3result,  TString::Format("J%s/F",fMVA3name.Data()).Data());
    tree->Branch(TString::Format("J%s",fMVA4name.Data()).Data(), &fTJMVA4result,  TString::Format("J%s/F",fMVA4name.Data()).Data());
    }

    // tracking
    tree->Branch("JNTracks"            , &fTNTracks,                "JNTracks/I");
    tree->Branch("JNTracksHS"          , &fTNTracksHS,              "JNTracksHS/I");
    tree->Branch("JNTracksPU"          , &fTNTracksPU,              "JNTracksPU/I");
    tree->Branch("JNTracksFilled"      , &fTNTracksFilled,          "JNTracksFilled/I");
    if(fAddTacksToTree){
    tree->Branch("TrkPt"               , &fTTrkPt,                  "TrkPt[JNTracksFilled]/F");
    tree->Branch("TrkEta"              , &fTTrkEta,                 "TrkEta[JNTracksFilled]/F");
    tree->Branch("TrkPhi"              , &fTTrkPhi,                 "TrkPhi[JNTracksFilled]/F");
    tree->Branch("Trkd0"               , &fTTrkd0,                  "Trkd0[JNTracksFilled]/F");
    tree->Branch("Trkd0errWrtPV"       , &fTTrkd0errWrtPV,          "Trkd0errWrtPV[JNTracksFilled]/F");
    tree->Branch("TrkZ0"               , &fTTrkZ0,                  "TrkZ0[JNTracksFilled]/F");
    tree->Branch("TrkZ0errWrtPV"       , &fTTrkZ0errWrtPV,          "TrkZ0errWrtPV[JNTracksFilled]/F");
    tree->Branch("TrkNPixHits"         , &fTTrkNPixHits,            "TrkNPixHits[JNTracksFilled]/I");
    tree->Branch("TrkNSCTHits"         , &fTTrkNSCTHits,            "TrkNSCTHits[JNTracksFilled]/I");
    tree->Branch("TrkNTRTHits"         , &fTTrkNTRTHits,            "TrkNTRTHits[JNTracksFilled]/I");
    tree->Branch("TrkChiSq"            , &fTTrkChiSq,               "TrkChiSq[JNTracksFilled]/F");
    tree->Branch("TrkNDof"             , &fTTrkNDof,                "TrkNDof[JNTracksFilled]/F");
    tree->Branch("TrkMCBarCode"        , &fTTrkMCBarCode,           "TrkMCBarCode[JNTracksFilled]/I");
    tree->Branch("TrkOrigin"           , &fTTrkOrigin,              "TrkOrigin[JNTracksFilled]/I");
    }

    // ROOTCore JVT
    tree->Branch("corrJVF_RootCore",    &fTJcorrJVF_RootCore,       "corrJVF_RootCore/F");
    tree->Branch("RpT_RootCore",        &fTJRpT_RootCore,           "RpT_RootCore/F");
    tree->Branch("JVT_RootCore",        &fTJJVT_RootCore,           "JVT_RootCore/F");

  if(Debug()) cout <<"Analysis_PileUpStudiesTreeFiller::AddBranches End" << endl;
    return;
}

void Analysis_PileUpStudiesTreeFiller::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_PileUpStudiesTreeFiller::ResetBranches Start" << endl;
  
    // Event Info
    fTEventNumber           = -9999;
    fTRunNumber             = -9999;
    fTWeight                = -999.99;
    fTPUWeight              = -999.99;
    fTLumiXSecWeight        = -999.99;
    fTm1MuScaleFactor       = -999.99;
    fTm1MuTrigScaleFactor   = -999.99;
    fTNPV                   = -99;
    fTNPVtruth              = -99;
    fTMu                    = -99;
    fTNBJets                = -99;
    fTGRL                   = -1;
    fTRhoKt4LC              = -999.99;
    fTnPUTracks             = -99;
    fTnTracksNotFromPV      = -99;

    // Zboson
    fTZpt         = -999.99;
    fTZeta        = -999.99;
    fTZphi        = -999.99;
    fTZmass       = -999.99;   

    // muons
    fTNMuons      = -99;
    for(int i=0;i<MaxNMuons; ++i){
        fTMuonPt[i]           = -999.99;
        fTMuonEta[i]          = -999.99;
        fTMuonPhi[i]          = -999.99;
        fTMuonPtUnsmeared[i]  = -999.99;
        fTMuonEtaUnsmeared[i] = -999.99;
        fTMuonPhiUnsmeared[i] = -999.99;
    }

    // njets
    fTNJetspt20eta2p1   = 0;
    fTNJetspt30eta2p1   = 0;
    fTNJetspt40eta2p1   = 0;
    fTNJetspt50eta2p1   = 0;
    fTNJetspt20eta2p1HS = 0;
    fTNJetspt30eta2p1HS = 0;
    fTNJetspt40eta2p1HS = 0;
    fTNJetspt50eta2p1HS = 0;
    

    // leading three jets
    fTJLeadingTruthPt= -999.99;
    fTJ0Pt          = -999.99;
    fTJ1Pt          = -999.99;
    fTJ2Pt          = -999.99;
    fTJ0Eta         = -999.99;
    fTJ1Eta         = -999.99;
    fTJ2Eta         = -999.99;
    fTJ0Phi         = -999.99;
    fTJ1Phi         = -999.99;
    fTJ2Phi         = -999.99;
    fTJ0ZDPhi       = -999.99;
    fTJ1ZDPhi       = -999.99;
    fTJ2ZDPhi       = -999.99;
    fTJ0MuonLeadDR  = -999.99;
    fTJ1MuonLeadDR  = -999.99;
    fTJ2MuonLeadDR  = -999.99;
    fTJ0MuonTrailDR = -999.99;
    fTJ1MuonTrailDR = -999.99;
    fTJ2MuonTrailDR = -999.99;

    
    // jet vars
    fTJindex                 = -999;
    fTJpt                    = -999.99;         
    fTJeta                   = -999.99;
    fTJphi                   = -999.99;
    fTJareacorrpt            = -999.99;
    fTJconstitpt             = -999.99;
    fTJtruthpt               = -999.99;         
    fTJtrutheta              = -999.99;
    fTJtruthphi              = -999.99;
    fTJtruthMatchDR          = -999.99;
    fTJtruthInTimept         = -999.99;      
    fTJtruthInTimeMatchDR    = -999.99;
    fTJdRTruthPt4            = -999.99;
    fTJdRTruthPt6            = -999.99;
    fTJdRTruthPt8            = -999.99;
    fTJdRTruthPt10           = -999.99;
    fTJisPU                  = -999;
    fTJisHS                  = -999;
    fTJisParton_uds          = -999;
    fTJisParton_c            = -999;
    fTJisParton_b            = -999;
    fTJisParton_g            = -999;
    fTJisBtagged             = -999;
    fTJArea                  = -999.99;
    fTJJVF                   = -999.99;
    fTJstdJVF                = -999.99;
    fTJMaxJVF                = -999.99;
    fTJHSPVtrkPtlead         = -999.99;
    fTJHSPVtrkPttrail        = -999.99; 
    fTJLeadTrkPt             = -999.99;
    fTJHSPVTrk12dR           = -999.99;
    fTJHSPVtrkSumOverPt      = -999.99;
    fTJHSPVtrkSumOverPtdRW   = -999.99;
    fTJHSPVTrkWidth          = -999.99;
    fTJHSPVTrkMeanDR         = -999.99;
    fTJMVA1result             = -999.99;
    fTJMVA2result             = -999.99;
    fTJMVA3result             = -999.99;
    fTJMVA4result             = -999.99;
    fTJTiming                = -999.99;
    fTJEMFrac                = -999.99;
    fTJTruthMF               = -999.99;
    fTJPUtrkPtSumOverPt      = -999.99;
    fTJCorrPUtrkPtSumOverPt  = -999.99;
    fTJCorr2PUtrkPtSumOverPt = -999.99;
    fTJPUtrkPtSum            = -999.99;
    fTJPUtrkPtSumdRWeight    = -999.99;
    fTJPUtrkPtSumSquared     = -999.99;
    fTJHStrkPtSum            = -999.99;
    fTJHStrkPtSumSquared     = -999.99;
    fTJHStrkPtSumdRWeight    = -999.99;
    fTJisBadLoose            = -999;
    fTRhoGoodPUTracks        = -999.99;
    fTRhoPUTracks            = -999.99;
    fTJZDPhi                 = -999.99;
    fTJMuonLeadDR            = -999.99;
    fTJMuonTrailDR           = -999.99;
    fTJNClusters             = -999;
    fTJNumTowers             = -999.99;
    fTJCaloWidth             = -999.99;
    fTJCaloWidth2            = -999.99;
    fTJAnnulusPtRatio0to1    = -999.99;
    fTJAnnulusPtRatio0to2    = -999.99;
    fTJAnnulusPtRatio1to2    = -999.99;
    fTJAnnulusPtRatio2to4    = -999.99;
    fTJAverageDeltaRSquared  = -999.99;
    fTJLeadTrkAllHSPVdZ      = -999.99;
    fTJBetaStar              = -999.99;
    fTJNContribVtxs          = -999;
    fTJNPVCorrJVF            = -999.99;
    fTJNContribVtxsCorrJVF   = -999.99;
    fTJMuCorrJVF             = -999.99;
    fTJMuCorrJVF0p1          = -999.99;
    fTJMuCorrJVF0p14         = -999.99;
    fTJnPUTrkCorrJVF         = -999.99;
    fTClosebykt              = -999.99;
    fTClosebydR              = -999.99;
    fTClosebydRpT            = -999.99;
    fTJJVT_RootCore          = -999.99;
    fTJcorrJVF_RootCore      = -999.99;
    fTJRpT_RootCore          = -999.99;

         

    fTNTracks               = 0;
    fTNTracksPU             = 0;
    fTNTracksHS             = 0;
    fTNTracksFilled         = 0;
    for(int i=0;i<MaxNTrks; ++i){
        fTTrkPt             [i] = -999.99;
        fTTrkEta            [i] = -999.99;
        fTTrkPhi            [i] = -999.99;
        fTTrkZ0             [i] = -999.99;
        fTTrkd0             [i] = -999.99;
        fTTrkChiSq          [i] = -999.99;
        fTTrkNDof           [i] = -999.99;
        fTTrkOrigin         [i] = -99;
        fTTrkNPixHits       [i] = -99;
        fTTrkNSCTHits       [i] = -99;
        fTTrkNTRTHits       [i] = -99;
        fTTrkMCBarCode      [i] = -99;
        fTTrkZ0errWrtPV     [i] = -999.99;
        fTTrkd0errWrtPV     [i] = -999.99;
    }

    // vertex
    fTHSVtxTruthminDz      = -999.99;
    fTVtxDzTruth           = -999.99;
    fTHSPV_mindZPUVtx      = -999.99;
    fTNMCVertWithin1mmFromHS = -999;

  if(Debug()) cout <<"Analysis_PileUpStudiesTreeFiller::ResetBranches End" << endl;
    return;
}


void Analysis_PileUpStudiesTreeFiller::FillEventVars(TTree *tree, const MomKey JetKey, const MomKey TrackKey){
  if(Debug()) cout <<"Analysis_PileUpStudiesTreeFiller::FillEventVars Begin" << endl;
  
    TString   simpleTrackKey = TrackKey.Data(); simpleTrackKey.ReplaceAll("Ghost", "");

    // Event Info
    fTEventNumber                 = Int("EventNumber");
    fTRunNumber                   = Int("RunNumber");
    fTNPV                         = Exists("NPV")? Int("NPV"):-1;
    fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
    fTMu                          = Float("averageIntPerXing");                  
    fTRhoKt4LC                    = Float("Eventshape_rhoKt4LC")/1000.;
    fTRhoGoodPUTracks             = Exists("RhoGoodPUTracks_"+simpleTrackKey)? Float("RhoGoodPUTracks_"+simpleTrackKey):-1;
    fTRhoPUTracks                 = Exists("RhoGoodPUTracks")? Float("RhoGoodPUTracks"):-1;
    fTnPUTracks                   = Exists("nPUTracks_"+simpleTrackKey)           ? Int("nPUTracks_"+simpleTrackKey):-1;
    fTnTracksNotFromPV            = Exists("nTrksNotFromPV_"+simpleTrackKey)      ? Int("nTrksNotFromPV_"+simpleTrackKey):-1;
    if(Exists("jetsAntiKt4LCTopoB")){
    fTNBJets                      = jets("AntiKt4LCTopoB");
    }
    fTWeight                      = Float("EventWeight");
    if(Exists("PUWeight")){
    fTPUWeight                    = Float("PUWeight");
    }
    fTLumiXSecWeight              = DefaultWeight();
    fTGRL                         = Bool("grl");
    if(isMC()){
    if(Exists("m1MuTrigScaleFactor")){
    fTm1MuTrigScaleFactor         = Float("m1MuTrigScaleFactor");
    }
    if(Exists("m1MuScaleFactor")){
    fTm1MuScaleFactor             = Float("m1MuScaleFactor");
    }
    }

    // Zboson
    if(Exists("recosZCandidate")){
        if(recos("ZCandidate")==1){
            fTZpt                      = reco(0,"ZCandidate").p.Pt();
            fTZeta                     = reco(0,"ZCandidate").p.Eta();
            fTZphi                     = reco(0,"ZCandidate").p.Phi();
            fTZmass                    = reco(0,"ZCandidate").p.M();
        }
    }

    if(Exists("muonsgood")){
        fTNMuons           = muons("good");
        for(int iMu=0; iMu<muons("good"); ++iMu){
            if (iMu>1) break;
            fTMuonPt [iMu]          = muon(iMu).p.Pt();
            fTMuonEta[iMu]          = muon(iMu).p.Eta();
            fTMuonPhi[iMu]          = muon(iMu).p.Phi();
            fTMuonPtUnsmeared [iMu] = muon(iMu).Float("ptUnsmeared");
            fTMuonEtaUnsmeared[iMu] = muon(iMu).Float("etaUnsmeared");
            fTMuonPhiUnsmeared[iMu] = muon(iMu).Float("phiUnsmeared");
        }
    }

    for(int iJet=0; iJet<jets(JetKey); ++iJet){
        if(jet(iJet,JetKey).p.Pt() > 20)                                      fTNJetspt20eta2p1++;
        if(jet(iJet,JetKey).p.Pt() > 30)                                      fTNJetspt30eta2p1++;
        if(jet(iJet,JetKey).p.Pt() > 40)                                      fTNJetspt40eta2p1++;
        if(jet(iJet,JetKey).p.Pt() > 50)                                      fTNJetspt50eta2p1++;
        if(jet(iJet,JetKey).p.Pt() > 20 && jet(iJet, JetKey).Bool("isHSJet")) fTNJetspt20eta2p1HS++;
        if(jet(iJet,JetKey).p.Pt() > 30 && jet(iJet, JetKey).Bool("isHSJet")) fTNJetspt30eta2p1HS++;
        if(jet(iJet,JetKey).p.Pt() > 40 && jet(iJet, JetKey).Bool("isHSJet")) fTNJetspt40eta2p1HS++;
        if(jet(iJet,JetKey).p.Pt() > 50 && jet(iJet, JetKey).Bool("isHSJet")) fTNJetspt50eta2p1HS++;
    }

    if(jets(JetKey)>0){
        fTJ0Pt  = jet(0, JetKey).p.Pt();
        fTJ0Eta = jet(0, JetKey).p.Eta();
        fTJ0Phi = jet(0, JetKey).p.Phi();
        if(Exists("recosZCandidate")){
        if(recos("ZCandidate")==1) fTJ0ZDPhi       = jet(0, JetKey).p.DeltaPhi(reco(0,"ZCandidate").p);
        }
        if(Exists("muonsgood")){
        if(muons("good")>0)        fTJ0MuonLeadDR  = jet(0, JetKey).p.DeltaR(muon(0, "good").p);
        if(muons("good")>1)        fTJ0MuonTrailDR = jet(0, JetKey).p.DeltaR(muon(1, "good").p);
        }
    }
    if(jets(JetKey)>1){
        fTJ1Pt  = jet(1, JetKey).p.Pt();
        fTJ1Eta = jet(1, JetKey).p.Eta();
        fTJ1Phi = jet(1, JetKey).p.Phi();
        if(Exists("recosZCandidate")){
        if(recos("ZCandidate")==1) fTJ1ZDPhi       = jet(1, JetKey).p.DeltaPhi(reco(0,"ZCandidate").p);
        }
        if(Exists("muonsgood")){
        if(muons("good")>0)        fTJ1MuonLeadDR  = jet(1, JetKey).p.DeltaR(muon(0, "good").p);
        if(muons("good")>1)        fTJ1MuonTrailDR = jet(1, JetKey).p.DeltaR(muon(1, "good").p);
        }
    }
    if(jets(JetKey)>2){
        fTJ2Pt  = jet(2, JetKey).p.Pt();
        fTJ2Eta = jet(2, JetKey).p.Eta();
        fTJ2Phi = jet(2, JetKey).p.Phi();
        if(Exists("recosZCandidate")){
        if(recos("ZCandidate")==1) fTJ2ZDPhi       = jet(2, JetKey).p.DeltaPhi(reco(0,"ZCandidate").p);
        }
        if(Exists("muonsgood")){
        if(muons("good")>0)        fTJ2MuonLeadDR  = jet(2, JetKey).p.DeltaR(muon(0, "good").p);
        if(muons("good")>1)        fTJ2MuonTrailDR = jet(2, JetKey).p.DeltaR(muon(1, "good").p);
        }
    }

    if(jets("AntiKt4Truth")>0){
        fTJLeadingTruthPt = jet(0,"AntiKt4Truth").p.Pt(); 
    }

  // Vertices --------------------------------
  // truth  ---
  if(Debug()) cout << "Analysis_PileUpStudiesTreeFiller :: Filling vertex info " << endl; 
  if(isMC() && vtxs("Truth")>0){

    // distance between truth HS and closest truth PU vertex
    float mindZTruthTruth=9999;
    for (int iV=1; iV<vtxs("Truth"); ++iV){ // index >=1 excludes HS 
        float dzTruthTruth = fabs(vtx(0,"Truth").x.z() - vtx(iV, "Truth").x.z());
        if(dzTruthTruth < mindZTruthTruth) {mindZTruthTruth = dzTruthTruth;}
    }
    fTHSVtxTruthminDz = (mindZTruthTruth < 9999 ? mindZTruthTruth : -1); // distance between truth HS and closest truth PU vertex
    

    // distance between reco  HS and truth HS PV
    float mindZTruthReco= vtxs()>0 ? fabs(vtx(0).x.z() - vtx(0, "Truth").x.z()): -1;
    fTVtxDzTruth  = mindZTruthReco; // distance between reco  HS and truth HS PV
  }

  // reco ---
  if(vtxs()>0){
    float mindZ=99999;
    for(int iV=1; iV<vtxs(); ++iV){
        if(fabs(vtx(iV).x.z() - vtx(0).x.z()) < mindZ){ 
            mindZ = fabs(vtx(iV).x.z() - vtx(0).x.z());
//            if (mindZ < 0.5) {vtx(0).Show(); vtx(iV).Show();}
        }
    }
    fTHSPV_mindZPUVtx = ((mindZ < 99999)? mindZ : -1);                   // distance between reco HS and closest PU vertex 
  }

  // Number of Truth Vertices within 1 mm of the HS PV
  //
  if(isMC() && vtxs("Truth")>0 && vtxs()>0 ){
      int NMerged=0;
      for(int iV=1; iV<vtxs("Truth"); ++iV){
            float dzTruthTruth = fabs(vtx(0).x.z() - vtx(iV, "Truth").x.z());  
            if(dzTruthTruth<1) NMerged++;
      }
      fTNMCVertWithin1mmFromHS = NMerged;
  }
  // -------------------------


  if(Debug()) cout <<"Analysis_PileUpStudiesTreeFiller::FillEventVars End" << endl;
  return;
}

///============================================================
/// Fill jet vars 
///============================================================
void Analysis_PileUpStudiesTreeFiller::FillJetVars(TTree* tree, int jindex, const MomKey JetKey, const MomKey JVFKey, const MomKey trackKey){
  if(Debug()) cout <<"Analysis_PileUpStudiesTreeFiller::FillJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex, JetKey));
  MomentObj *myJVFzero     = myjet->Objs(JVFKey)    >0   ? (MomentObj*)(myjet->Obj(JVFKey, 0)):NULL;
  MomentObj *stdJVFzero    = myjet->Exists("JVFLinks")   ? (MomentObj*)(myjet->Obj("JVFLinks", 0)):NULL;
  Particle* trklead        = NULL;
  Particle* trktrail       = NULL;
  if(myJVFzero!=NULL){
  trklead        = myJVFzero->Objs("tracks")>0    ?  (Particle*) myJVFzero->Obj("tracks",0) : NULL;
  trktrail       = myJVFzero->Objs("tracks")>1    ?  (Particle*) myJVFzero->Obj("tracks",1) : NULL;
  }
  TString   simpleTrackKey = trackKey.Data(); simpleTrackKey.ReplaceAll("Ghost", "");
  
  // set tree variables ---------------------------------
  fTJindex               = jindex;
  fTJpt                  = myjet->p.Pt();       
  fTJeta                 = myjet->p.Eta();
  fTJphi                 = myjet->p.Phi();
  fTJconstitpt           = myjet->Float("constscale_pt");
  fTJareacorrpt          = myjet->Float("areacorr_pt");
  fTJisPU                = myjet->Bool("isPUJet");
  fTJisHS                = myjet->Bool("isHSJet");
  fTJisBtagged           = myjet->Exists("btag")? myjet->Float("btag"): -1;
  fTJdRTruthPt4          = myjet->Float("dRTruthPt4");
  fTJdRTruthPt6          = myjet->Float("dRTruthPt6");
  fTJdRTruthPt8          = myjet->Float("dRTruthPt8");
  fTJdRTruthPt10         = myjet->Float("dRTruthPt10");
  fTJisParton_uds        = myjet->Bool("Parton_uds");
  fTJisParton_c          = myjet->Bool("Parton_c");
  fTJisParton_b          = myjet->Bool("Parton_b");
  fTJisParton_g          = myjet->Bool("Parton_g");
  fTJisBadLoose          = myjet->Int("isBadLoose");
  fTJArea                = myjet->Float("Area");
  fTJTiming              = myjet->Float("Timing");
  fTJEMFrac              = myjet->Float("emfrac");
  if(myjet->Bool("isHSJet") && myjet->Exists("AntiKt4Truth_match")){
      Particle *truth = (Particle*) myjet->Obj("AntiKt4Truth_match");
      fTJtruthpt             = truth->p.Pt();                                                                                         
      fTJtrutheta            = truth->p.Eta();
      fTJtruthphi            = truth->p.Phi();
      fTJtruthMatchDR        = truth->p.DeltaR(myjet->p);
  }
  if(myjet->Exists("InTimeAntiKt4Truth_match")){
      Particle *truth = (Particle*) myjet->Obj("InTimeAntiKt4Truth_match");
      fTJtruthInTimept      = truth->p.Pt();                                                                                         
      fTJtruthInTimeMatchDR = truth->p.DeltaR(myjet->p);
  }

  fTJJVF                    = myjet->Float(JVFKey+"_JVF");
  fTJMaxJVF                 = myjet->Float(JVFKey+"_MaxJVF");
  fTJstdJVF                 = stdJVFzero!=NULL? stdJVFzero->Float("JVF"): -1;
  fTJNPVCorrJVF             = myjet->Float(JVFKey+"_NPVCorrJVF");
  fTJNContribVtxsCorrJVF    = myjet->Float(JVFKey+"_NContribVtxsCorrJVF");
  fTJMuCorrJVF              = myjet->Float(JVFKey+"_MuCorrJVF");
  fTJMuCorrJVF0p1           = myjet->Float(JVFKey+"_MuCorrJVF0p1");
  fTJMuCorrJVF0p14          = myjet->Float(JVFKey+"_MuCorrJVF0p14");
  fTJnPUTrkCorrJVF          = myjet->Float(JVFKey+"_nPUTrkCorrJVF");
  fTJHSPVtrkPtlead          = trklead ==0? -1: trklead ->p.Pt();            
  fTJHSPVtrkPttrail         = trktrail==0? -1: trktrail->p.Pt();           
  fTJHSPVTrkWidth           = myjet->Float(JVFKey+"_TrkWidth");
  fTJHSPVTrkMeanDR          = myjet->Float(JVFKey+"_TrkMeanDR");
  fTJHSPVTrk12dR            = (trktrail > 0 && trklead >0)? trklead ->p.DeltaR(trktrail->p): -0.5;   
  fTJHSPVtrkSumOverPt       = myjet->Float(JVFKey+"_HSPVtrkSumOverPt")>0 ? myjet->Float(JVFKey+"_HSPVtrkSumOverPt"):0;
  fTJHSPVtrkSumOverPtdRW    = myjet->Float(JVFKey+"_HSPVtrkSumOverPtdRWeight");
  fTJPUtrkPtSumOverPt       = myjet->Exists(trackKey+"_PUTrkSumPt")?myjet->Float(trackKey+"_PUTrkSumPt")/myjet->p.Pt():-1;
  fTJCorrPUtrkPtSumOverPt   = myjet->Float(JVFKey+"_CorrPUtrkPtSumOverPt");
  fTJCorr2PUtrkPtSumOverPt  = myjet->Float(JVFKey+"_Corr2PUtrkPtSumOverPt");
  fTJBetaStar               = myjet->Float(JVFKey+"_BetaStar");
  fTJPUtrkPtSum             = myjet->Exists(trackKey+"_PUTrkSumPt")? myjet->Float(trackKey+"_PUTrkSumPt"):-1;
  fTJHStrkPtSum             = myjet->Exists(trackKey+"_HSTrkSumPt")? myjet->Float(trackKey+"_HSTrkSumPt"):-1;
  fTJPUtrkPtSumSquared      = myjet->Float(JVFKey+"_PUTrkSumPtSquared");
  fTJHStrkPtSumSquared      = myjet->Float(JVFKey+"_HSTrkSumPtSquared");         
  fTJPUtrkPtSumdRWeight     = myjet->Float(JVFKey+"_PUTrkSumPtdRWeight");
  fTJHStrkPtSumdRWeight     = myjet->Float(JVFKey+"_HSTrkSumPtdRWeight");
  fTJLeadTrkAllHSPVdZ       = myjet->Float(JVFKey+"_leadTrkPVdZ");
  fTJNContribVtxs           = myjet->Int(JVFKey+"_NContribVtxs");
  fTClosebydR               = myjet->Float(JVFKey+"_Closeby_dR");
  fTClosebykt               = myjet->Float(JVFKey+"_Closeby_kt");
  fTClosebydRpT             = (myjet->Float(JVFKey+"_Closeby_dRindex")>=0 && myjet->Float(JVFKey+"_Closeby_dRindex")<jets(JetKey))? jet(myjet->Float(JVFKey+"_Closeby_dRindex"), JetKey).p.Pt() : -1; 
  fTJJVT_RootCore           = myjet->Float("JVT_RootCoreJVT");
  fTJcorrJVF_RootCore       = myjet->Float("corrJVF_RootCoreJVT");
  fTJRpT_RootCore           = myjet->Float("RpT_RootCoreJVT");

  if(myjet->Exists(trackKey)){
    if(myjet->Objs(trackKey)>0) fTJLeadTrkPt  = ((Particle*) myjet->Obj(trackKey,0))->p.Pt();
    else                        fTJLeadTrkPt  = -0.5;
  }else                         fTJLeadTrkPt  = -0.5;


  // Leptons and Z-boson ---------------------------
  if(Exists("recosZCandidate")){
  if(recos("ZCandidate")==1){
  fTJZDPhi               = myjet->p.DeltaPhi(reco(0,"ZCandidate").p);
  }
  }
  if(Exists("muonsgood")){
  if(muons("good")>0){
  fTJMuonLeadDR          = myjet->p.DeltaR(muon(0,"good").p);
  }
  if(muons("good")>1){
  fTJMuonTrailDR         = myjet->p.DeltaR(muon(1,"good").p);
  }
  }

  // calo vars -------------------------------------------
  fTJNClusters             = myjet->Int("NClusters");
  fTJNumTowers             = myjet->Float("NumTowers");
  fTJCaloWidth             = myjet->Float("CaloWidth");
  fTJCaloWidth2            = myjet->Float("CaloWidth2");
  fTJAnnulusPtRatio0to1    = myjet->Float("AnnulusPtRatio0to1");
  fTJAnnulusPtRatio0to2    = myjet->Float("AnnulusPtRatio0to2");
  fTJAnnulusPtRatio1to2    = myjet->Float("AnnulusPtRatio1to2");
  fTJAnnulusPtRatio2to4    = myjet->Float("AnnulusPtRatio2to4");
  fTJAverageDeltaRSquared  = myjet->Float("AverageDeltaRSquared");


  // Tracks ----------------------
  if(myjet->Exists(trackKey)){
      for(int iTrk=0; iTrk<myjet->Objs(trackKey);++iTrk){
          Particle *trk         = (Particle*) myjet->Obj(trackKey,iTrk);
          fTNTracks++;
          if(trk->Int("origin") ==0) {fTNTracksHS++;} // associated tracks from the HS PV
          if(trk->Int("origin") > 0) {fTNTracksPU++;} // associated tracks from PU PVs

          if(iTrk>=MaxNTrks)       continue;
          fTNTracksFilled++;
          fTTrkPt[iTrk]            = trk->p.Pt(); 
          fTTrkEta[iTrk]           = trk->p.Eta();
          fTTrkPhi[iTrk]           = trk->p.Phi();
          fTTrkZ0[iTrk]            = trk->Float("z0_wrtPV");
          fTTrkd0[iTrk]            = trk->Float("d0_wrtPV");
          fTTrkZ0errWrtPV[iTrk]    = trk->Exists("err_z0_wrtPV")?trk->Float("err_z0_wrtPV"):-999;
          fTTrkd0errWrtPV[iTrk]    = trk->Exists("err_d0_wrtPV")?trk->Float("err_d0_wrtPV"):-999;
          fTTrkOrigin[iTrk]        = trk->Int("origin");
          fTTrkNPixHits[iTrk]      = trk->Int("nPixHits")+trk->Int("nPixelDeadSensors");
          fTTrkNSCTHits[iTrk]      = trk->Int("nSCTHits")+trk->Int("nSCTDeadSensors");
          fTTrkNTRTHits[iTrk]      = trk->Int("nTRTHits");
          fTTrkChiSq[iTrk]         = trk->Float("chi2");
          fTTrkNDof[iTrk]          = trk->Float("ndof");
//          fTTrkMCBarCode[iTrk]     = trk->Int("mc_barcode");
      }
  }


    // TMVA -----------------------------------------
    if(doMVAEval && myjet->p.Pt()>20 && myjet->p.Pt()<50 && fabs(myjet->p.Eta())<2.4 ){ // MVA only trained in 20 to 50 GeV jets with |eta|<2.4
        if(fMVA1name.Length()>0 && fMVA1file.Length()>0) {  if      (fTJPUtrkPtSum ==0 && fTJHStrkPtSum==0  ) fTJMVA1result = 1.1;
                                                            else if (fTJHStrkPtSum ==0                      ) fTJMVA1result = 0;
                                                            else                                              fTJMVA1result = fTMVAreader->EvaluateMVA(fMVA1name);}
        if(fMVA2name.Length()>0 && fMVA2file.Length()>0) {  if      (fTJPUtrkPtSum ==0 && fTJHStrkPtSum==0  ) fTJMVA2result = 1.1;
                                                            else if (fTJHStrkPtSum ==0                      ) fTJMVA2result = 0;
                                                            else                                              fTJMVA2result = fTMVAreader->EvaluateMVA(fMVA2name);}
        if(fMVA3name.Length()>0 && fMVA3file.Length()>0) {  if      (fTJPUtrkPtSum ==0 && fTJHStrkPtSum==0  ) fTJMVA3result = 1.1;
                                                            else if (fTJHStrkPtSum ==0                      ) fTJMVA3result = 0;
                                                            else                                              fTJMVA3result = fTMVAreader->EvaluateMVA(fMVA3name);}
        if(fMVA4name.Length()>0 && fMVA4file.Length()>0) {  if      (fTJPUtrkPtSum ==0 && fTJHStrkPtSum==0  ) fTJMVA4result = 1.1;
                                                            else if (fTJHStrkPtSum ==0                      ) fTJMVA4result = 0;
                                                            else                                              fTJMVA4result = fTMVAreader->EvaluateMVA(fMVA4name);}
    }
    

  if(Debug()) cout <<"Analysis_PileUpStudiesTreeFiller::FillJetVars End" << endl;
  return;
}
