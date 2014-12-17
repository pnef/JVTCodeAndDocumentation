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

void runTreeFiller(TString mode       = "local",         // local, lite, or cluster
TString identifier = "ClustersAndTruth",                      // tag 
//TString dataset  = "PythJ1and2mc12aJETMET.jetmet2012pileupcustom",
//TString dataset   = "PythJ1to3mc12aJETMET.jetmet2012",
//TString dataset   = "PythJXmc12aJETMETshort.jetmet2012",  // dataset name
TString dataset   = "Zmumu_MC12_Sherpa_COMMON_short.jetmet2012",
TString username   = "pnef",                               // username (e.g. swiatlow, fizisist)
bool mcweights     = true,                                 // use mc weights?
bool debug         = false,                                // debug mode
bool doPRWGen      = false,                                 // PRWgen
Long64_t nentries  = 5000                              // nevents
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
    TString treename("physics");
   
 
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
    bool doJetCalibrations= false;
    bool doJetTriggers    = false;
    bool doPRW            = false;
    bool doParentChild    = false;
    bool doTrack          = true;
    bool doVertex         = true;
    bool doLCCluster      = true; 
    bool doEMCluster      = true;
    bool doTruth          = true;
    bool doEMJets         = false;
    bool doLCJets         = true;
    bool doJet4           = true;
    bool doJet6           = false; 
    bool doVectorJVFLinks = false;
    bool doTruthJets      = true;
    bool doOOTtruthJet4   = false;
    bool doTruthLinks     = false;
    bool doPhotons        = false;
    bool doElectrons      = false;
    bool doConstitLinks   = false;
    bool doMuons          = false;
    bool doCOMMON            = true;
    int  counterMax       = -1;
    TString prwTypes      = "EF_j15_a4tchad,EF_j25_a4tchad,EF_j35_a4tchad,EF_j45_a4tchad,EF_j45_a4tchad_L2FS_L1J15,EF_j55_a4tchad,EF_j80_a4tchad,EF_j110_a4tchad,EF_j145_a4tchad,EF_j180_a4tchad,EF_j220_a4tchad,EF_j360_a4tchad,EF_j460_a4tchad";

    float LUMI               =1;



    ///----------------------------------------------------------------
    /// Nominal Configuration
    /// set up the actual analyses
    ///----------------------------------------------------------------
    
    Config* EventTreeFiller_ClustersAndTruth = new Config("EventTreeFiller_ClustersAndTruth",configfile);
    EventTreeFiller_ClustersAndTruth->Set("ANALYSIS","EventTreeFiller_ClustersAndTruth");
    EventTreeFiller_ClustersAndTruth->Set("DEBUG",debug);
    EventTreeFiller_ClustersAndTruth->Set("TreeName","JetTree0");

    Config* chain = new Config("chain",configfile);
    chain->AddVec("ANALYSIS");
    chain->Add("ANALYSIS",EventTreeFiller_ClustersAndTruth);


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
    chain->Set("DOVTX"           , doVertex        );
    chain->Set("DOPHOTON"        , doPhotons       );
    chain->Set("DOELECTRONS"     , doElectrons     );
    chain->Set("DOMUONS"         , doMuons         );
    chain->Set("JETTYPES"        , ""              );
    chain->Set("BTAGS"           , ""              ); 
    chain->Set("PRWTYPES"        , prwTypes        );
    chain->Set("LUMI",             LUMI);
    chain->Set("DOCOMMON"          , doCOMMON          );
    chain->Set("runCluster",       runCluster);
  
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
