//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Nov 12 10:42:46 2011 by ROOT version 5.28/00f
// from TTree EventBuilder_jetmet2012/EventBuilder_jetmet2012
// found on file: NTUP_JETMET.543738._000084.root.2
//////////////////////////////////////////////////////////

#ifndef EventBuilder_jetmet2012_h
#define EventBuilder_jetmet2012_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
//include <TSelector.h>
#include "EventBuilderBase.h"
#include <TLorentzVector.h>
//#include "AnalysisBase.h"
#include "TriggerMenuNtuple/TMNUtil.h" //Called in Notify(), so need include, not just forward declaration
#include "TrigRootAnalysis/TrigDecisionToolD3PD.h"
//#include "applyjetcalibration/applyjetcalibration/applyjetcalibration.h"
//#include "utils/TriggerConfViewer/TrigDefsGJ.h"
#include "MomentObj.h"

namespace Root {class TGoodRunsList; class TPileupReweighting;}
class JetCalibrationTool;
class HforToolD3PD;
class MultijetJESUncertaintyProvider;
class JetSmearingTool;
namespace MuonSmear {class SmearingClass;}
using namespace D3PD;

class METUtility;
#include "MissingETUtility/METUtility.h"

namespace JetMETSyst{
  enum{
    NONE = 0,
    JESUP = 1,
    JESDOWN = 2,
    JER = 3,
    BD = 4,
    CD = 5,
    LD = 6,
    BU = 7,
    CU = 8,
    LU = 9,
    METSCALESTUP = 10,
    METSCALESTDOWN = 11,
    METRESOST = 12,
    MMSUP = 13,
    MMSDOWN = 14,
    MIDUP = 15,
    MIDDOWN = 16,
    MSCALEUP = 17,
    MSCALEDOWN = 18,
    MEFFUP = 19,
    MEFFDOWN = 20,
    MTRIGUP = 21,
    MTRIGDOWN = 22,
    MASSUP = 23,
    MASSDOWN = 24
  };
}

class EventBuilder_jetmet2012 : public EventBuilderBase {
public :

  enum{
    JET4LC,
    JET6LC,
    JET4EM,
    JET6EM,
    JET4TRUTH,
    JET6TRUTH,
    JET4OOTTRUTH,
    JET4INTIMETRUTH
  };

   Bool_t   CopyEvent(AnalysisBase* evt);
   void  WriteAuxTrees(TDirectory* outfile);

   virtual bool DoCustomSetup(){return true;}

   void Initialize();
   Bool_t   Notify();
 
   Bool_t   CopyTruth();
   Bool_t   CopyTruthFull();

   Bool_t   CopyPhotons();
   Bool_t   CopyVertices();	
   Bool_t   CopyJet(TString, TString);
   Bool_t		CopyJets();
   Bool_t   CopyStandardJet(int JetFlag);
   Bool_t		CopyTracks();
   Bool_t   CopyClusters();
   Bool_t   CopyMuons();
   Bool_t   CopyElectrons();
   Bool_t   CopyMET();

   Bool_t   DoBtaggingSF(bool sys);

   Bool_t   GetJetCalibrationTools();

   Bool_t   KillPileupFluctuations();
   Bool_t   DoHFOverlap();
   Bool_t   MuonHasTriggerMatch(double mu_eta,
             double mu_phi,
             const vector<int> *trig_ef_trigmuonef_signature,
             int& efindex,
             int& eftrackindex,
             unsigned int trig_ef_trigmuonef_n,
             const vector<vector<float> > *trig_ef_trigmuonef_track_cb_eta,
             const vector<vector<float> > *trig_ef_trigmuonef_track_cb_phi,
             const vector<vector<int> > *trig_ef_trigmuonef_track_cb_hascb);
   Bool_t   CopyMuonTriggers();
   Bool_t   CopyElectronTriggers();
   Bool_t   CopyJetTriggers();

   Bool_t   AddTruthParentChild();

   Bool_t   AddBtags();
   Bool_t   AddBtag(TString, TString);

   Bool_t   ParseJetStrings();
   Bool_t   AddPRWWeights();
   Bool_t   SetupPRW();


   bool doAFII;
   bool doSMWZ;
   bool doCOMMON;
   bool doOffset;
   bool doJet4;
   bool doJet6;
   bool doTruthJets;
   bool doOOTtruthJet4;
   bool doInTimeTruthJet4;
   bool doVectorJVFLinks;
   bool doTruthMFLinks;
   bool doConstitLinks;
   bool doTruthParentChild;
   bool doJetCalibrations;
   bool doLCJets;
   bool doEMJets;
   bool doTrack;
   bool doPhotons;
   bool doLCCluster;
   bool doEMCluster;
   bool doTruth;
   bool doVertex;
   bool doMuons;
   bool doElectrons;
   bool doMuonTriggers;
   bool doElectronTriggers;
   bool doMuonTriggersMatch;
   bool doJetTriggers;
   bool doPRW;
   bool doHFOR;
   bool doMETRefFinal;
   bool doBarcode;
   bool Rename; // for TrackJetBTag naming mechanism
   bool doTrkJetBTag;
   bool doLeptonTruth;
   bool doMuSmear;
   bool doMCGRL;
   bool doITKFixes;

   float muRef;
   float npvRef;

   TString mJETTYPEVEC;
   TString mBTAGVEC;
   TString mPRWVEC;

	METUtility *m_util;
	vector<float> jet_AntiKt4LCTopo_recalibrated_pt;
  vector<TLorentzVector> calJets;
  vector<float> closeBy;

   // jet calibration tools
   JetCalibrationTool*   data4LCCalibrator;
   JetCalibrationTool* data4EMCalibrator;

   JetCalibrationTool* data6LCCalibrator;
   JetCalibrationTool* data6EMCalibrator;

   JetCalibrationTool*  mc4LCCalibrator;
   JetCalibrationTool* mc4EMCalibrator;

   JetCalibrationTool* mc6LCCalibrator;
   JetCalibrationTool* mc6EMCalibrator;


   JetCalibrationTool* afiimc4LCCalibrator;
   JetCalibrationTool* afiimc4EMCalibrator;

   JetCalibrationTool* afiimc6LCCalibrator;
   JetCalibrationTool* afiimc6EMCalibrator;

   MultijetJESUncertaintyProvider* m_myJesNorm;
   MultijetJESUncertaintyProvider* m_myJesAFII;

   JetSmearingTool* m_myJER;

   // muon smear tools
   MuonSmear::SmearingClass* MCPSmear;


   //Syst tracking
   int whichsyste;

    // trigger things
   Long64_t m_TrigEntry;
   Long64_t m_TMNEntry;
   TTree* m_trigConfTree;
   TrigDecisionToolD3PD* m_trigdefs;
   Root::TGoodRunsList* myGRL;
   Root::TPileupReweighting* myPRW;

   map<MomKey, Root::TPileupReweighting*> prwMap;

   HforToolD3PD* myHforTool;

   EventBuilder_jetmet2012(TTree * /*tree*/ =0);
   virtual ~EventBuilder_jetmet2012();

   AnalysisBase* fEvt;

   MomentObj* barcode_link;

   ClassDef(EventBuilder_jetmet2012,0);

protected:
   	bool		m_isNewInit;
};

#endif

#ifdef EventBuilder_jetmet2012_cxx

Bool_t EventBuilder_jetmet2012::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

	//flag for skimming, do not want to copy auxiliary tree if we already have
  //cout << "starting notify" << endl;
 
  //if(fEvt->Debug()) cout << "setting up trees for trig decisions, etc." << endl;
 // cout << " setting up test cout " << endl;
	m_isNewInit = true;
	
	//flags for trigger navigation and trigger decision readers, whether we have to load branches
	m_TrigEntry = -1;
	m_TMNEntry = -1;

	if(!Tree()->GetTree()) return kTRUE;

  //cout << "going to do the casting bits" << endl;

	TFile* file = 0;
	TChain* chain = dynamic_cast< TChain* >( Tree() );
  //cout << "casted" << endl;
	if( chain ) {
		// We are running locally...
//    cout <<"localmode?" << endl;
		file = chain->GetFile();
	} else {
		// We are running on PROOF:
    //cout <<"proofmode?" << endl;
		file = Tree()->GetCurrentFile();
	}

  //cout << "grabbing the trigconftree" << endl;
	m_trigConfTree = (TTree*)( file->Get( "qcdMeta/TrigConfTree" ) );

  //cout << "grabbed the tree" << endl;

  if(m_trigdefs) delete m_trigdefs;

  //cout << "init the tool!" << endl;
	m_trigdefs = new TrigDecisionToolD3PD(Tree()->GetTree(), m_trigConfTree,"trig_");
	//m_trigdefs->initTrigMetaData(m_trigConfTree,Tree());
 
  //cout << "All good on NotifyOnDemand() " << endl; 
  //if(fEvt->Debug()) cout << " All good on NotifyOnDemand()" << endl;
	return kTRUE;
}


#endif // #ifdef EventBuilder_jetmet2012_cxx
