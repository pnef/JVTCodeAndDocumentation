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

void runLargeR(TString mode       = "local",         // local, lite, or cluster
TString identifier = "LargeR",                      // tag 
// TString dataset    = "PythJXmc12aJETMET.jetmet2012pileupcustom",  // dataset name
// TString dataset    = "Pythia8_Wprime_WZqqqq_m2000_COMMON.jetmet2012pileupcustom",
// TString dataset    = "Pythia8_Wprime_WZqqqq_m1000_COMMON.jetmet2012pileupcustom",
 TString dataset    = "Pythia8_Wprime_WZqqqq_m1000_COMMON_1.jetmet2012pileupcustom",
// TString dataset   = "Zprime_14TeV_COMMON_r4770.jetmet2012pileupcustom",
TString username   = "pnef",                               // username (e.g. swiatlow, fizisist)
bool mcweights     = true,                                 // use mc weights?
bool debug         = false,                                // debug mode
Long64_t nentries  = 300                              // nevents
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
    TString configfile("../config/LargeR.config");

    
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
    TString url(mode);
    TString path("");
    if(mode.CompareTo("lite")==0) {
        url = "lite://";
        path = pathLite;
    }
    else if(mode.CompareTo("cluster")==0) {
        url = TString(username+"@atlprf01.slac.stanford.edu");
        path = pathCluster;
    }
    
    // Make an options file, edit as needed
    TFile* options = new TFile("options.root","RECREATE");

    ///----------------------------------------------------------------
    /// Overall Configuration
    ///----------------------------------------------------------------
    bool doBasic          = true;
    bool doJetCalibrations= true;
    bool doJetTriggers    = false;
    bool doConstitLinks   = true;
    bool doPRW            = false;
    bool doParentChild    = true;
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
    bool doMuons          = false;
    int  counterMax       = -1;
    bool doSMWZfixes      = false;
    bool doCOMMONfixes    = true;
    TString prwTypes      = "";


    float maxJetTruthMatchdR = 0.3;
    float minJetPUdR         = 0.6;
    bool  applylumiWeights   = false;
    bool  requireTrueHSvertex= false;
    float LUMI               = 1;
    bool  doTrkFromVxt       = true;
    bool  RecoverBtracks     = true;
    bool  doQCDSelection      = true;

    bool doCOMMON            = false;
    if(dataset.Contains("COMMON")){
        doCOMMON= true;
    }
   string TrimmedJetNames = 
"jet_0p10Trim,\
jet_0p09Trim,\
jet_0p08Trim,\
jet_0p07Trim,\
jet_0p06Trim,\
jet_0p05Trim,\
jet_0p04Trim,\
jet_0p03Trim,\
jet_0p02Trim,\
jet_0p01Trim,\
jet_0p10TrimCorrJVF,\
jet_0p09TrimCorrJVF,\
jet_0p08TrimCorrJVF,\
jet_0p07TrimCorrJVF,\
jet_0p06TrimCorrJVF,\
jet_0p05TrimCorrJVF,\
jet_0p04TrimCorrJVF,\
jet_0p03TrimCorrJVF,\
jet_0p02TrimCorrJVF,\
jet_0p01TrimCorrJVF,\
jet_0p00TrimCorrJVF,\
jet_JVFCleansedTrim0p00,\
jet_JVFCleansedTrim0p01,\
jet_JVFCleansedTrim0p02,\
jet_JVFCleansedTrim0p03,\
jet_JVFCleansedTrim0p04,\
jet_JVFCleansedTrim0p05,\
jet_JVFCleansedTrim0p06,\
jet_JVFCleansedTrim0p07,\
jet_JVFCleansedTrim0p08,\
jet_JVFCleansedTrim0p09,\
jet_JVFCleansedTrim0p010,\
jet_LinCleansed0p55Trim0p00,\
jet_LinCleansed0p55Trim0p01,\
jet_LinCleansed0p55Trim0p02,\
jet_LinCleansed0p55Trim0p03,\
jet_LinCleansed0p55Trim0p04,\
jet_LinCleansed0p55Trim0p05,\
jet_LinCleansed0p55Trim0p06,\
jet_LinCleansed0p55Trim0p07,\
jet_LinCleansed0p55Trim0p08,\
jet_LinCleansed0p55Trim0p09,\
jet_LinCleansed0p55Trim0p010,\
jet_MyClensing,\
jet_MyClensing0p01Trim";


    ///----------------------------------------------------------------
    /// Nominal Configuration
    /// set up the actual analyses
    ///----------------------------------------------------------------
    Config* QCDSelection = new Config("QCDSelection",configfile);
    QCDSelection->Set("ANALYSIS","QCDCommonSelection");
    QCDSelection->Set("DEBUG",debug);
    
    Config* PileUpStudiesBase = new Config("PileUpStudiesBase",configfile);
    PileUpStudiesBase->Set("ANALYSIS","PileUpStudiesBase");
    PileUpStudiesBase->Set("DEBUG",debug);

    Config* PileUpStudies0 = new Config("PileUpStudies0",configfile);
    PileUpStudies0->Set("ANALYSIS","PileUpStudies");
    PileUpStudies0->Set("TrackSel", "tracksGoodSel0");
    PileUpStudies0->Set("DEBUG",debug);

    Config* LargeR = new Config("LargeR",configfile);
    LargeR->Set("ANALYSIS","LargeR");
    LargeR->Set("TrackSel", "tracksGoodSel0");
    LargeR->Set("DEBUG",debug);
    
    Config* TreeFiller = new Config("LargeRTreeFiller",configfile);
    TreeFiller->Set("ANALYSIS","LargeRTreeFiller");
    TreeFiller->Set("TrackSel", "tracksGoodSel0");
    TreeFiller->Set("DEBUG",debug);
    
    Config* chain = new Config("chain",configfile);
    chain->AddVec("ANALYSIS");
    chain->Add("ANALYSIS",QCDSelection);
    chain->Add("ANALYSIS",PileUpStudiesBase);
    chain->Add("ANALYSIS",PileUpStudies0);
    chain->Add("ANALYSIS",LargeR);
    chain->Add("ANALYSIS",TreeFiller);

    // set up configurations, this overwrites configs from configfile
    chain->Set("DOJETCALIBRATIONS",doJetCalibrations);
    chain->Set("DOJETTRIGGERS"   , doJetTriggers);
    chain->Set("DOJET4"          , doJet4          );
    chain->Set("DOJET6"          , doJet6          );
    chain->Set("DOVECTORJVFLINKS", doVectorJVFLinks);
    chain->Set("DOCONSTITLINKS"  , doConstitLinks);
    chain->Set("DOLCJETS"        , doLCJets        );
    chain->Set("DOEMJETS"        , doEMJets        );
    chain->Set("COUNTERMAX"      , counterMax      );
    chain->Set("DEBUG"           , debug           );
    chain->Set("MCWEIGHTS"       , mcweights       );
    chain->Set("PILE"            , doPRW           );
    chain->Set("DOBASIC"         , doBasic         );
    chain->Set("DOTRUTHLINKS"    , doTruthLinks    );
    chain->Set("DOPARENTCHILD"   , doParentChild   );
    chain->Set("DOTRACK"         , doTrack         );
    chain->Set("DOLCCLUSTER"     , doLCCluster     );
    chain->Set("DOEMCLUSTER"     , doEMCluster     );
    chain->Set("DOTRUTH"         , doTruth         );
    chain->Set("DOTRUTHJETS"     , doTruthJets     );
    chain->Set("DOOOTTRUTHJET4"  , doOOTtruthJet4  );
    chain->Set("DOVTX"           , doVertex        );
    chain->Set("DOPHOTON"        , doPhotons       );
    chain->Set("DOELECTRONS"     , doElectrons     );
    chain->Set("DOMUONS"         , doMuons         );
    chain->Set("JETTYPES"        , ""              );
    chain->Set("BTAGS"           , ""              ); 
    chain->Set("PRWTYPES"        , prwTypes        );
    chain->Set("applylumiWeights", applylumiWeights);
    chain->Set("LUMI"            , LUMI);
    chain->Set("requireTrueHSvertex", requireTrueHSvertex);
    chain->Set("DOSMWZfixes",         doSMWZfixes);
    chain->Set("DOCOMMON"          , doCOMMON          );
    chain->Set("DOCOMMONfixes",       doCOMMONfixes);
    chain->Set("DOTRKFROMVTX",       doTrkFromVxt);
    chain->Set("TRIMMEDJETNAMES",    TrimmedJetNames);
    chain->Set("RecoverBtracks",   RecoverBtracks   );
    chain->Set("doQCDSelection"    , doQCDSelection); 
  

    // set cuts and selections
    chain->Set("MAXJETTRUTHMATCHDR"       , maxJetTruthMatchdR);
    chain->Set("MINJETPUDR"               , minJetPUdR);
    chain->Write();


    WriteJetCalibrationObjects(options);
    WriteGRLObject("data12_8TeV.periodAB_HEAD_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml");



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
