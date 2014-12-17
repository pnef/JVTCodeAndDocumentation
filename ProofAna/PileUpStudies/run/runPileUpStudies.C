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

void runPileUpStudies(TString mode       = "cluster",         // local, lite, or cluster
TString identifier = "PileUpStudies_rev368524_tracksSelITK5GeV",                      // tag 
//TString dataset  = "PythJ1and2mc12aJETMET.jetmet2012pileupcustom",
//TString dataset   = "PythJXmc12aJETMETshort.jetmet2012pileupcustom",  // dataset name
//TString dataset   = "mc12_14TeV_Pythia8_J2_ITK_80_80_COMMON.jetmet2012pileupcustom", // ITK Upgrade, 14 TeV, sigma(mu=80), mu=80, 25ns 
//TString dataset   = "mc12_14TeV_Pythia8_J2_ITK_80_80_COMMON_short.jetmet2012pileupcustom", // ITK Upgrade, 14 TeV, sigma(mu=80), mu=80, 25ns 
//TString dataset   = "mc12_14TeV_Pythia8_J2_ITK_140_140_COMMON.jetmet2012pileupcustom", // ITK Upgrade, 14 TeV, sigma(mu=140), mu=140, 25ns 
//TString dataset   = "mc12_14TeV_Pythia8_J2_ITK_140_140_COMMON_short.jetmet2012pileupcustom", // ITK Upgrade, 14 TeV, sigma(mu=140), mu=140, 25ns 
//TString dataset = "mc12_14TeV_Pythia8_J2_ITK_250_250_COMMON_short.jetmet2012pileupcustom",
TString dataset = "mc12_14TeV_Pythia8_J2_ITK_250_250_COMMON.jetmet2012pileupcustom",
//TString dataset = "mc12_14TeV_Pythia8_J2_ITK_300_300_COMMON_short.jetmet2012pileupcustom",
//TString dataset = "mc12_14TeV_Pythia8_J2_ITK_300_300_COMMON.jetmet2012pileupcustom",
TString username   = "pnef",                               // username (e.g. swiatlow, fizisist)
bool mcweights     = true,                                 // use mc weights?
bool debug         = false,                                // debug mode
bool doPRWGen      = false,                                 // PRWgen
Long64_t nentries  = 1000                              // nevents
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
    bool doBasic          = true;
    bool doJetCalibrations= true;
    bool doJetTriggers    = false;
    bool doPRW            = false;
    bool doParentChild    = false;
    bool doTrack          = true;
    bool doVertex         = true;
    bool doLCCluster      = true; 
    bool doEMCluster      = false;
    bool doEMJets         = false;
    bool doLCJets         = true;
    bool doJet4           = true;
    bool doJet6           = false; 
    bool doVectorJVFLinks = true;
    bool doTruth          = true;
    bool doTruthJets      = true;
    bool doOOTtruthJet4   = false;
    bool doTruthLinks     = true;
    bool doPhotons        = false;
    bool doElectrons      = false;
    bool doConstitLinks   = true;
    bool doMuons          = false;
    int  counterMax       = -1;
    TString prwTypes      = "EF_j15_a4tchad,EF_j25_a4tchad,EF_j35_a4tchad,EF_j45_a4tchad,EF_j45_a4tchad_L2FS_L1J15,EF_j55_a4tchad,EF_j80_a4tchad,EF_j110_a4tchad,EF_j145_a4tchad,EF_j180_a4tchad,EF_j220_a4tchad,EF_j360_a4tchad,EF_j460_a4tchad";

    bool        doMVAEval  = false;
    TString     MVA1name   = "KNN50";
    TString     MVA1file   = "skimmed_Jpt20to50.20131108.10.44_PileUpStudies_firstJPt20Eta2p4_newStuff.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN20_50_100.nPUtrkCorrJVF_RpT.JHStrkPtSumGE0_KNN50.weights.xml";
    TString     MVA2name   = "KNN50cat";
    TString     MVA2file   = "";
    TString     MVA3name   = "KNN100";
    TString     MVA3file   = "";
    TString     MVA4name   = "KNN100cat";
    TString     MVA4file   = "";

    bool doCOMMON = false;
    bool doITKFixes=false;
    if(dataset.Contains("COMMON")){
        doCOMMON= true;
    }
    if(dataset.Contains("ITK")){
        doITKFixes= true;
    }

    float maxJetTruthMatchdR = 0.3;
    float minJetPUdR         = 0.6;
    bool  requireTrueHSvertex= false;
    float LUMI               =1;
    bool  doSMWZfixes         = false;
    bool  doCOMMONfixes       = false;
    bool  addLeptonZInfo      = false;
    bool  doQCDSelection      = true;
    bool  doInTimeTruthJet4   = false;
    bool  addTracksToTree     = false;
    bool  doTrkFromVxt        = true;
    bool  RecoverBtracks      = true;



    ///----------------------------------------------------------------
    /// Nominal Configuration
    /// set up the actual analyses
    ///----------------------------------------------------------------
    Config* QCDSelection = new Config("QCDSelection",configfile);
    QCDSelection->Set("ANALYSIS","QCDCommonSelection_Forward");
    QCDSelection->Set("DEBUG",debug);
    
    Config* PileUpStudiesBase = new Config("PileUpStudiesBase",configfile);
    PileUpStudiesBase->Set("ANALYSIS","PileUpStudiesBase");
    PileUpStudiesBase->Set("DEBUG",debug);

    Config* PileUpStudies0 = new Config("PileUpStudies0",configfile);
    PileUpStudies0->Set("ANALYSIS","PileUpStudies");
    PileUpStudies0->Set("TrackSel", "tracksGoodSelITK5GeV");
    PileUpStudies0->Set("DEBUG",debug);
    
    Config* PileUpStudiesTreeFiller0 = new Config("PileUpStudiesTreeFiller0",configfile);
    PileUpStudiesTreeFiller0->Set("ANALYSIS","PileUpStudiesTreeFiller");
    PileUpStudiesTreeFiller0->Set("DEBUG",debug);
    PileUpStudiesTreeFiller0->Set("TreeName","JetTree0");
    PileUpStudiesTreeFiller0->Set("TrackSel", "tracksGoodSelITK5GeV");
    
    Config* PileUpStudiesTruthTreeFiller = new Config("PileUpStudiesTruthTreeFiller",configfile);
    PileUpStudiesTruthTreeFiller->Set("ANALYSIS","PileUpStudiesTruthTreeFiller");
    PileUpStudiesTruthTreeFiller->Set("DEBUG",debug);

    Config* PileUpStudiesEventTreeFiller= new Config("PileUpStudiesEventTreeFiller",configfile);
    PileUpStudiesEventTreeFiller->Set("ANALYSIS","PileUpStudiesEventTreeFiller");
    PileUpStudiesEventTreeFiller->Set("TrackSel", "tracksGoodSelITK5GeV");
    PileUpStudiesEventTreeFiller->Set("DEBUG",debug);

    Config* chain = new Config("chain",configfile);
    chain->AddVec("ANALYSIS");
    chain->Add("ANALYSIS",QCDSelection);
    chain->Add("ANALYSIS",PileUpStudiesBase);
    chain->Add("ANALYSIS",PileUpStudies0);
//  chain->Add("ANALYSIS",PileUpStudies1);
    chain->Add("ANALYSIS",PileUpStudiesTreeFiller0);
//  chain->Add("ANALYSIS",PileUpStudiesTreeFiller1);
    chain->Add("ANALYSIS",PileUpStudiesTruthTreeFiller);
    chain->Add("ANALYSIS",PileUpStudiesEventTreeFiller);


    // set up configurations, this overwrites configs from configfile
    chain->Set("DOJETCALIBRATIONS",doJetCalibrations);
    chain->Set("DOJETTRIGGERS"   , doJetTriggers);
    chain->Set("DOJET4"          , doJet4          );
    chain->Set("DOJET6"          , doJet6          );
    chain->Set("DOVECTORJVFLINKS", doVectorJVFLinks);
    chain->Set("DOLCJETS"        , doLCJets        );
    chain->Set("DOEMJETS"        , doEMJets        );
    chain->Set("COUNTERMAX"      , counterMax      );
    chain->Set("DEBUG"           , debug           );
    chain->Set("MCWEIGHTS"       , mcweights       );
    chain->Set("PILE"            , doPRW           );
    chain->Set("DOBASIC"         , doBasic         );
    chain->Set("DOTRUTHLINKS"    , doTruthLinks    );
    chain->Set("DOCONSTITLINKS"  , doConstitLinks  );
    chain->Set("DOPARENTCHILD"   , doParentChild   );
    chain->Set("DOTRACK"         , doTrack         );
    chain->Set("DOLCCLUSTER"     , doLCCluster     );
    chain->Set("DOEMCLUSTER"     , doEMCluster     );
    chain->Set("DOTRUTH"         , doTruth         );
    chain->Set("DOTRUTHJETS"     , doTruthJets     );
    chain->Set("DOOOTTRUTHJET4"  , doOOTtruthJet4  );
    chain->Set("DOINTIMETRUTHJET4", doInTimeTruthJet4);
    chain->Set("DOVTX"           , doVertex        );
    chain->Set("DOPHOTON"        , doPhotons       );
    chain->Set("DOELECTRONS"     , doElectrons     );
    chain->Set("DOMUONS"         , doMuons         );
    chain->Set("JETTYPES"        , ""              );
    chain->Set("BTAGS"           , ""              ); 
    chain->Set("PRWTYPES"        , prwTypes        );
    chain->Set("ADDTRACKSTOTREE" , addTracksToTree );
    chain->Set("ADDLEPTONZINFOTOTREE", addLeptonZInfo);
    chain->Set("LUMI",             LUMI);
    chain->Set("requireTrueHSvertex", requireTrueHSvertex);
    chain->Set("DOSMWZfixes"       , doSMWZfixes          );
    chain->Set("DOCOMMON"          , doCOMMON          );
    chain->Set("DOCOMMONfixes"     , doCOMMONfixes          );
    chain->Set("doQCDSelection"    , doQCDSelection); 
    chain->Set("DOTRKFROMVTX",       doTrkFromVxt);
    chain->Set("runCluster",       runCluster);
    chain->Set("doMVAEval",        doMVAEval);
    chain->Set("MVA1name"        , MVA1name         );
    chain->Set("MVA2name"        , MVA2name         );
    chain->Set("MVA3name"        , MVA3name         );
    chain->Set("MVA4name"        , MVA4name         );
    chain->Set("MVA1file"        , MVA1file         );
    chain->Set("MVA2file"        , MVA2file         );
    chain->Set("MVA3file"        , MVA3file         );
    chain->Set("MVA4file"        , MVA4file         );
    chain->Set("RecoverBtracks",   RecoverBtracks   );
    chain->Set("DOITKFIXES",   doITKFixes);
  

    // set cuts and selections
    chain->Set("MAXJETTRUTHMATCHDR"       , maxJetTruthMatchdR);
    chain->Set("MINJETPUDR"               , minJetPUdR);
    chain->Write();


    
    if(doPRWGen){
    WritePRWConfigure(options, "MC12a");
    }
    WriteJetCalibrationObjects(options);
    WriteGRLObject("data12_8TeV.periodAB_HEAD_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml");
    WriteDijetPRWO(options, "EF_j15_a4tchad");
    WriteDijetPRWO(options, "EF_j25_a4tchad");
    WriteDijetPRWO(options, "EF_j35_a4tchad");
    WriteDijetPRWO(options, "EF_j45_a4tchad");
    WriteDijetPRWO(options, "EF_j45_a4tchad_L2FS_L1J15");
    WriteDijetPRWO(options, "EF_j55_a4tchad");
    WriteDijetPRWO(options, "EF_j80_a4tchad");
    WriteDijetPRWO(options, "EF_j110_a4tchad");
    WriteDijetPRWO(options, "EF_j145_a4tchad");
    WriteDijetPRWO(options, "EF_j180_a4tchad");
    WriteDijetPRWO(options, "EF_j220_a4tchad");
    WriteDijetPRWO(options, "EF_j360_a4tchad");
    WriteDijetPRWO(options, "EF_j460_a4tchad");



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
    //WriteGroomedPRWO(options , "EF_j240_a4tc_EFFS" );
    confProofAna->Write();
    options->Close();
    delete options;
 
 
    cout << "All setup, ready to go " << endl; 
    int runNevents=1; 
 

    // Decide to run local or on the cluster
    if(mode.CompareTo("local")==0) runLocal(dataset,treename,nentries);
    else{
        runProof(url,dataset,-1,treename);
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
