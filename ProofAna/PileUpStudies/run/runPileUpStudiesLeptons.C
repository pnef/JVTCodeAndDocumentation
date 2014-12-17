#include "../scripts/runLocal.C"
#include "../scripts/runProof.C"
#include "../scripts/helperFunc.C"
#include "../scripts/helperJetMETCommon.C"
#include "../scripts/loadLibraries.C"
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

void runPileUpStudiesLeptons(TString mode       = "local",         // local, lite, or cluster
TString identifier = "PileUpStudies_DiLeptonSelection_rev368524",    
// --------------- Data -------------
//TString dataset    = "MuonsSMWZPeriod_All.jetmet2012pileupcustom",
//TString dataset    = "MuonsSMWZPeriodAB_short.jetmet2012pileupcustom",
//TString dataset    = "MuonsSMWZPeriodAB.jetmet2012pileupcustom",
//TString dataset    = "MuonsSMWZPeriodCD.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodE1.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodE2.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodE3.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodE4.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodE5.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodE6.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodG1.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodG2.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodG3.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodH.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodI.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodJ.jetmet2012pileupcustom",
//TString dataset      = "MuonsSMWZPeriodL.jetmet2012pileupcustom",
/// ---------- Zmumu ------------
//TString dataset    = "Zmumu_PowhegPythia8_MC12_COMMON_short.jetmet2012pileupcustom",
//TString dataset    = "Zmumu_PowhegPythia8_MC12_COMMON.jetmet2012pileupcustom",
//TString dataset    = "Zmumu_PowhegPythia8_MC12_JETMET.jetmet2012pileupcustom",
//TString dataset    = "Zmumu_PowhegPythia8_MC12_JETMET_short.jetmet2012pileupcustom",
//TString dataset    = "Zmumu_MC12_Sherpa_COMMON.jetmet2012pileupcustom",
TString dataset    = "Zmumu_MC12_Sherpa_COMMON_short.jetmet2012pileupcustom",
//TString dataset   = "Zmumu_MC12_AlpgenPythia_COMMON_short.jetmet2012pileupcustom",
//TString dataset   = "Zmumu_MC12_AlpgenPythia_COMMON.jetmet2012pileupcustom",
// >>> Overlay
//TString dataset    = "Zmumu_PowhegPythia8_Overlay_COMMON.jetmet2012pileupcustom",
//TString dataset    = "Zmumu_PowhegPythia8_Overlay_COMMON_short.jetmet2012pileupcustom",
//TString dataset    = "Zmumu_PowhegPythia8_Overlay_new_COMMON.jetmet2012pileupcustom",
// ------------ VBF Higgs -------------
//TString dataset   = "VBFH125_ZZ4lep_PowhegPythia8MC12_COMMON.jetmet2012pileupcustom",
// --------------- TTbar -----------
//TString dataset    = "PowhegPythia_ttbar_LeptonFilterSMWZ_short.jetmet2012pileupcustom",
// ------------------------------------------------------------------------------------------------
TString username   = "pnef",                               // username (e.g. swiatlow, fizisist)
bool mcweights     = true,                                 // use mc weights?
bool debug         = false,                                // debug mode
Long64_t nentries  = 500,                              // nevents
bool doPRW         = true,
int      nworkers  = -1,
bool fullLumi      = true
    ) 
{ 
    
    TString date(currentDateTime().c_str());
    identifier = date+"_"+identifier;

    ///----------------------------------------------------------------
    /// Load libraries , set the config file, treenam, and cluster info
    ///----------------------------------------------------------------

    cout << "trying to load libraries" << endl; 
    loadLibraries();

    cout << " Libraries loaded " << endl;

    // SetConfig 
    TString configfile("../config/pileupstudies.config");

    
    // Best to leave alone  
    TString pathLite("");
    TString pathCluster("root://atlprf01.slac.stanford.edu:2094//atlas/output/");
    pathCluster.Append(username);
    pathCluster.Append("/");

    // Determine eventbuilder from dataset name
    TString eventbuilder(dataset);
    eventbuilder.Remove(0,eventbuilder.Last('.')+1); 

    //Change if defaulting to wrong TTree, otherwise leave
    TString treename;
    if(dataset.Contains("SMWZ") || dataset.Contains("COMMON")){
        treename = "physics";
    } else{
        treename = "qcd";
    }
   
 
    ///----------------------------------------------------------------
    /// Filename paths, URLs for PROOF running
    ///----------------------------------------------------------------
    bool runCluster(false);
    TString url(mode);
    TString path("");
    if(mode.CompareTo("lite")==0) {
        url = "lite://";
        path = pathLite;
    }
    else if(mode.CompareTo("cluster")==0) {
        url = TString(username+"@atlprf01.slac.stanford.edu");
        path = pathCluster;
        runCluster = true;
    }
    
    // Make an options file, edit as needed
    TFile* options = new TFile("options.root","RECREATE");

    ///----------------------------------------------------------------
    /// Overall Configuration
    ///----------------------------------------------------------------
    bool doBasic             = true;
    bool doJetCalibrations   = true;
    bool doJetTriggers       = false;
    bool doMuonTriggers      = true;
    bool doMuonTriggersMatch = false;
    bool doParentChild       = false;
    bool doTrack             = true;
    bool doVertex            = true;
    bool doLCCluster         = true; 
    bool doLCJets            = true;
    bool doJet4              = true;
    bool doVectorJVFLinks    = true;
    bool doTruth             = true;
    bool doTruthJets         = true;
    bool doTruthLinks        = false;
    bool doInTimeTruthJet4   = false;
    bool doPhotons           = false;
    bool doElectrons         = true;
    bool doMuons             = true;
    bool doConstitLinks      = true;
    bool doMETRefFinal       = true;
    bool evtweights          = true;
    bool doMuSmear           = false;
    bool doSMWZ              = false;
    bool doCOMMON            = false;
    bool doSMWZfixes         = true;
    bool doCOMMONfixes       = false;
    bool doMCGRL             = false;
    if(dataset.Contains("SMWZ")){
        doSMWZ = true;
    }
    if(dataset.Contains("COMMON")){
        doCOMMON= true;
    }
    if(dataset.Contains("Overlay")){
        doMCGRL = true;
    }

    TString prwTypes      = "PeriodAB_lumi";
    if(fullLumi) prwTypes = "Full_lumi"; 


    float      LUMI                =1;
    int         counterMax         = -1;
    bool        doMVAEval          = false;
    TString     MVA1name           = "KNN50";
    TString     MVA1file           = "skimmed_Jpt20to50.20131108.10.44_PileUpStudies_firstJPt20Eta2p4_newStuff.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN20_50_100.nPUtrkCorrJVF_RpT.JHStrkPtSumGE0_KNN50.weights.xml";
    TString     MVA2name           = "KNN50cat";
    TString     MVA2file           = "";
    TString     MVA3name           = "KNN100";
    TString     MVA3file           = "";
    TString     MVA4name           = "KNN100cat";
    TString     MVA4file           = "";

    float maxJetTruthMatchdR  = 0.3;
    float minJetPUdR          = 0.6;
    bool  requireTrueHSvertex = false;
    bool  doLeptonSelection   = true;
    bool  doHZZSelection      = false;
    bool  doQCDSelection      = false;
    bool  RecoverBtracks      = true;
    bool  addTracksToTree     = false;
    bool  doTrkFromVxt        = true;



    ///----------------------------------------------------------------
    /// Nominal Configuration
    /// set up the actual analyses
    ///----------------------------------------------------------------

    
    Config* TopSelection = new Config("TopSelection",configfile);
    TopSelection->Set("ANALYSIS","TopCommonSelection");
    TopSelection->Set("DEBUG",debug);
    
    Config* PileUpStudiesBase = new Config("PileUpStudiesBase",configfile);
    PileUpStudiesBase->Set("ANALYSIS","PileUpStudiesBase");
    PileUpStudiesBase->Set("DEBUG",debug);
    
    Config* PileUpStudies = new Config("PileUpStudies",configfile);
    PileUpStudies->Set("ANALYSIS","PileUpStudies");
    PileUpStudies->Set("DEBUG",debug);
    PileUpStudies->Set("TrackSel", "tracksGoodSel0");
    
    Config* PileUpStudiesTreeFiller = new Config("PileUpStudiesTreeFiller",configfile);
    PileUpStudiesTreeFiller->Set("ANALYSIS","PileUpStudiesTreeFiller");
    PileUpStudiesTreeFiller->Set("TreeName","JetTree0");
    PileUpStudiesTreeFiller->Set("DEBUG",debug);

    Config* PileUpStudiesEventTreeFiller= new Config("PileUpStudiesEventTreeFiller",configfile);
    PileUpStudiesEventTreeFiller->Set("ANALYSIS","PileUpStudiesEventTreeFiller");
    PileUpStudiesEventTreeFiller->Set("TrackSel", "tracksGoodSel0");
    PileUpStudiesEventTreeFiller->Set("DEBUG",debug);

    Config* chain = new Config("chain",configfile);
    chain->AddVec("ANALYSIS");
    chain->Add("ANALYSIS",TopSelection);
    chain->Add("ANALYSIS",PileUpStudiesBase);
    chain->Add("ANALYSIS",PileUpStudies);
    chain->Add("ANALYSIS",PileUpStudiesTreeFiller);
    chain->Add("ANALYSIS",PileUpStudiesEventTreeFiller);


    // set up configurations, this overwrites configs from configfile
    chain->Set("DOJETCALIBRATIONS",doJetCalibrations);
    chain->Set("DOJETTRIGGERS"   , doJetTriggers   );
    chain->Set("DOMUONTRIGGERS"  , doMuonTriggers  );
    chain->Set("DOMUONTRIGGERSMATCH", doMuonTriggersMatch  );
    chain->Set("DOJET4"          , doJet4          );
    chain->Set("DOVECTORJVFLINKS", doVectorJVFLinks);
    chain->Set("DOLCJETS"        , doLCJets        );
    chain->Set("DOCONSTITLINKS"  , doConstitLinks  );
    chain->Set("COUNTERMAX"      , counterMax      );
    chain->Set("DEBUG"           , debug           );
    chain->Set("MCWEIGHTS"       , mcweights       );
    chain->Set("EVTWEIGHTS"      , evtweights      );
    chain->Set("PILE"            , doPRW           );
    chain->Set("DOBASIC"         , doBasic         );
    chain->Set("DOTRUTHLINKS"    , doTruthLinks    );
    chain->Set("DOPARENTCHILD"   , doParentChild   );
    chain->Set("DOTRACK"         , doTrack         );
    chain->Set("DOLCCLUSTER"     , doLCCluster     );
    chain->Set("DOTRUTH"         , doTruth         );
    chain->Set("DOTRUTHJETS"     , doTruthJets     );
    chain->Set("DOINTIMETRUTHJET4", doInTimeTruthJet4);
    chain->Set("DOVTX"           , doVertex        );
    chain->Set("DOPHOTON"        , doPhotons       );
    chain->Set("DOELECTRONS"     , doElectrons     );
    chain->Set("DOMUONS"         , doMuons         );
    chain->Set("JETTYPES"        , ""              );
    chain->Set("BTAGS"           , ""              ); 
    chain->Set("DOMETREFFINAL"   , doMETRefFinal   );
    chain->Set("PRWTYPES"        , prwTypes        );
    chain->Set("ADDTRACKSTOTREE" , addTracksToTree );
    chain->Set("RecoverBtracks",   RecoverBtracks   );
    chain->Set("LUMI",             LUMI);
    chain->Set("requireTrueHSvertex", requireTrueHSvertex);
    chain->Set("DOSMWZ"            , doSMWZ          );
    chain->Set("DOCOMMON"          , doCOMMON          );
    chain->Set("DOSMWZfixes"       , doSMWZfixes          );
    chain->Set("DOCOMMONfixes"     , doCOMMONfixes          );
    chain->Set("doLeptonSelection" , doLeptonSelection);
    chain->Set("doHZZSelection"    , doHZZSelection);
    chain->Set("doQCDSelection"    , doQCDSelection);
    chain->Set("DOMUSMEAR",          doMuSmear);
    chain->Set("DOMCGRL",            doMCGRL);
    chain->Set("DOTRKFROMVTX",       doTrkFromVxt);
    chain->Set("doMVAEval",        doMVAEval);
    chain->Set("runCluster",       runCluster);
    chain->Set("MVA1name"        , MVA1name         );
    chain->Set("MVA2name"        , MVA2name         );
    chain->Set("MVA3name"        , MVA3name         );
    chain->Set("MVA4name"        , MVA4name         );
    chain->Set("MVA1file"        , MVA1file         );
    chain->Set("MVA2file"        , MVA2file         );
    chain->Set("MVA3file"        , MVA3file         );
    chain->Set("MVA4file"        , MVA4file         );
  

    // set cuts and selections
    chain->Set("MAXJETTRUTHMATCHDR"       , maxJetTruthMatchdR);
    chain->Set("MINJETPUDR"               , minJetPUdR);
    chain->Write();


    
    if(fullLumi){
    WriteGRLObject("data12_8TeV.periodAllYear_HEAD_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml");
    WriteZllPRWO(options, "Full_lumi"); 
    }else{
    WriteGRLObject("data12_8TeV.periodAB_HEAD_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml");
    WriteZllPRWO(options, "PeriodAB_lumi"); 
    }
    WriteMuonUtilities(options);
    WriteBTagCalibObject(options,"MV1","0_7892");
    WriteJetCalibrationObjects(options);



    //if (chain->Exists("GRL")) WriteGRLObject(chain->String("GRL"));  

    ///----------------------------------------------------------------
    /// ProofAna global Config object
    ///----------------------------------------------------------------
    Config* confProofAna = new Config("ProofAna");
    
    confProofAna->Set("DEBUG"          , false        );  // "false", 0, "0" etc. also works
    confProofAna->Set("SAVETIMERS"     , false        );  // ProofAna timer histos in output file   
    confProofAna->Set("IDENTIFIER"     , identifier   );
    confProofAna->Set("DATASET"        , dataset      );
    confProofAna->Set("OUTPUTPATH"     , path         );
    confProofAna->Set("EVENTBUILDER"   , eventbuilder );
    confProofAna->Set("MERGE",true);     // enable dataset mode
   
    cout << "set eventbuilder to " << eventbuilder << endl;
 
    ///----------------------------------------------------------------
    /// Read information used in MC weighting, multi-dataset jobs
    ///----------------------------------------------------------------
    ReadDatasetInfo(dataset  , confProofAna  ); 
    confProofAna->Write();
    options->Close();
    delete options;
 
 
    cout << "All setup, ready to go " << endl; 
    int runNevents=1; 
 

    // Decide to run local or on the cluster
    if(mode.CompareTo("local")==0)     runLocal(dataset,treename,nentries);
    else{
        runProof(url,dataset,nworkers,treename);
    }
    gSystem->Unlink("options.root");

}


const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y%m%d.%H.%M", &tstruct);
    return buf;
}
