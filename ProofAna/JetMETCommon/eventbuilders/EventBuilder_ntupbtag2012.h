//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Nov 12 10:42:46 2011 by ROOT version 5.28/00f
// from TTree EventBuilder_ntupbtag2012/EventBuilder_ntupbtag2012
// found on file: NTUP_JETMET.543738._000084.root.2
//////////////////////////////////////////////////////////

#ifndef EventBuilder_ntupbtag2012_h
#define EventBuilder_ntupbtag2012_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
//include <TSelector.h>
#include "EventBuilderBase.h"
#include <TLorentzVector.h>
//#include "AnalysisBase.h"
#include "TriggerMenuNtuple/TMNUtil.h" //Called in Notify(), so need include, not just forward declaration
#include "TrigRootAnalysis/TrigDecisionToolD3PD.h"
#include "ApplyJetCalibration/ApplyJetCalibration/ApplyJetCalibration.h"
//#include "utils/TriggerConfViewer/TrigDefsGJ.h"

namespace root {class TPileupReweighting;}
namespace Root { class TGoodRunsList; }
using namespace D3PD;

class EventBuilder_ntupbtag2012 : public EventBuilderBase {
public :
   Bool_t   CopyEvent(AnalysisBase* evt);
   void  WriteAuxTrees(TDirectory* outfile);

   virtual bool DoCustomSetup(){return true;}

   void Initialize();
   Bool_t  Notify();
 
   Bool_t		CopyTruth();
   Bool_t    CopyTruthFull();

   Bool_t		CopyPhotons();
   Bool_t		CopyVertices();	
   Bool_t   CopyJet(TString, TString);
   Bool_t		CopyJets();
   Bool_t		CopyTracks();
   Bool_t   CopyClusters();
   Bool_t   CopyMuons();

   Bool_t   AddTruthParentChild();

   Bool_t   AddBtags();
   Bool_t   AddBtag(TString, TString);

   bool doOffset;
   bool doJet4;
   bool doJet6;
   bool doTruthJets;
   bool doLCJets;
   bool doEMJets;
   bool doTrack;
   bool doPhotons;
   bool doLCCluster;
   bool doEMCluster;
   bool doTruth;
   bool doTruthParentChild;
   bool doTruthHeavyQ;
   bool doTruthLightQ;
   bool doTruthStable;
   bool doVertex;
   bool doMuons;

   float muRef;
   float npvRef;

   TString mJETTYPEVEC;
   TString mBTAGVEC;

   JetCalibrationTool* calibAntiKt4TopoEM;
   JetCalibrationTool* calibAntiKt6TopoEM;
   JetCalibrationTool* calibAntiKt4LCTopo;
   JetCalibrationTool* calibAntiKt6LCTopo;
    // trigger things
   Long64_t m_TrigEntry;
   Long64_t m_TMNEntry;
   TTree* m_trigConfTree;
   TrigDecisionToolD3PD* m_trigdefs;
   Root::TGoodRunsList* myGRL;

   EventBuilder_ntupbtag2012(TTree * /*tree*/ =0);
   virtual ~EventBuilder_ntupbtag2012();

   AnalysisBase* fEvt;

   ClassDef(EventBuilder_ntupbtag2012,0);

protected:

   	bool		m_isNewInit;
};

#endif

#ifdef EventBuilder_ntupbtag2012_cxx

Bool_t EventBuilder_ntupbtag2012::Notify()
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


#endif // #ifdef EventBuilder_ntupbtag2012_cxx
