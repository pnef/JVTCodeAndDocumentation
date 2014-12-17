/**************************************************************************
 **
 **   File:         Analysis_LargeR.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      P. Nef
 **
 **************************************************************************/

#define Analysis_LargeR_cxx

#include "Analysis_LargeR.h"
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
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Subtractor.hh"

#if fjContribFLAG==1 // preprocessor flag defined in makefile
#include "JetCleanser/JetCleanser.hh"

using namespace fastjet::contrib;
#endif

///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_LargeR::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_LargeR: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();

  /// Set Pt and eta bounds for plots
  ptBounds.push_back(0.); 
  ptBounds.push_back(10); ptBounds.push_back(20.);
  ptBounds.push_back(30); ptBounds.push_back(40.);
  ptBounds.push_back(50); 
  ptBounds.push_back(60.); ptBounds.push_back(80.);
  ptBounds.push_back(100.); ptBounds.push_back(125.); 
  ptBounds.push_back(150.); ptBounds.push_back(200.); 
  ptBounds.push_back(350.); ptBounds.push_back(500.); 
  ptBounds.push_back(1000.); 
  ptBounds.push_back(4000.);

  etaBounds.push_back(0.);  etaBounds.push_back(1.8); 
  etaBounds.push_back(20.);

  // Specific Config
  fTrkSel = "tracksGoodSel0";
  Cfg()->Get("TrackSel", fTrkSel);

  if (Debug()) cout << "Analysis_LargeR: DEBUG Finish WorkerBegin()" << endl;  
} 

///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_LargeR::ProcessEvent()
{

  if (Debug()) cout << "Analysis_LargeR: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();
  
  // Event Selection and so on
  if(! Exists("PileUpStudiesSelection")) return true;
  if(! Bool("PileUpStudiesSelection")  ) return true;


  // START ---------------------------
  //  make truth jets
  AddStableParticles();
  MakeJets(fastjet::antikt_algorithm, 1,"truthsStable");
//  Show();
      

  // add truth W and Z
  AddVec("recosW");
  AddVec("recosZ");
  AddVec("recosTop");
  AddVec("recosAntiTop");
  for(int iT=0; iT<truths(); ++iT){
    if(    truth(iT).Int("pdgId") ==23 && truth(iT).Objs("children")==2 ) { Add("recosZ", &truth(iT));  }
    if(abs(truth(iT).Int("pdgId"))==24 && truth(iT).Objs("children")==2 ) { Add("recosW", &truth(iT));  }
    if(abs(truth(iT).Int("pdgId"))==6  && truth(iT).Objs("children")==2 ) {
        Particle* child1 = (Particle*) truth(iT).Obj("children",0);
        Particle* child2 = (Particle*) truth(iT).Obj("children",1);
        if(! (abs(child1->Int("pdgId"))==5 || abs(child1->Int("pdgId"))==24)) continue;
        if(! (abs(child2->Int("pdgId"))==5 || abs(child2->Int("pdgId"))==24)) continue;

        if(truth(iT).Int("pdgId")==6) Add("recosTop", &truth(iT));
        else                          Add("recosAntiTop", &truth(iT));
        }
  }

    AddVec("jetsKt3LCTopo");
    AddVec("jetsKt3LCTopoSubjets");

	

    // ---------------------------------------------------------------------------
	//Now we setup the fat jets, first trim and then find the subjets.
	AddVec("jetsFat10AntiKt");
	double jetRad = 1.0;

    // get vector of PS with clusters and ghost tracks
	vector<fastjet::PseudoJet> particles = ObjsToPJ("clustersLCTopo");
    AddAsGhosts(fTrkSel,   &particles, true);

	
    // cluster and calculate hard scattering jet areas
	fastjet::JetDefinition *fJetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, jetRad, fastjet::Best); 
	fastjet::GhostedAreaSpec areaSpec(5.0);
	fastjet::ClusterSequenceArea hardClustSeq(particles,*fJetDef, areaSpec);
	vector<fastjet::PseudoJet> UnGroomedJets = fastjet::sorted_by_pt(hardClustSeq.inclusive_jets(10.)); // ungroomed fat jets
	 
    // JetCleanser ---------
    #if fjContribFLAG==1 // preprocessor flag defined in makefile --------------------
    fastjet::JetDefinition subjet_def_A(fastjet::kt_algorithm, 0.3);                                                                                                                              
    // JVF cleanser                                                                                                                                                                               
    JetCleanser jvf_cleanser_A(subjet_def_A, JetCleanser::jvf_cleansing, JetCleanser::input_nc_together);                                                                                         
    // linear cleanser                                                                                                                                                                            
    JetCleanser linear_cleanser_A(subjet_def_A, JetCleanser::linear_cleansing, JetCleanser::input_nc_together);                                                                                   
    #endif

	for(unsigned int iJet = 0; iJet < UnGroomedJets.size(); iJet++){
	  Particle* UnGroomed = new Particle();
	  if (UnGroomedJets[iJet].pt() < 50) continue;
	  UnGroomed->p.SetPtEtaPhiE(UnGroomedJets[iJet].pt(), UnGroomedJets[iJet].eta(),
				    UnGroomedJets[iJet].phi(), UnGroomedJets[iJet].e());
	  UnGroomed->Set("ptConstScale", UnGroomedJets[iJet].pt());
	  UnGroomed->Set("mConstScale",  UnGroomedJets[iJet].m());
	  UnGroomed->Set("etaConstScale",UnGroomedJets[iJet].eta()); 
	  UnGroomed->Set("phiConstScale",UnGroomedJets[iJet].phi()); 
	  UnGroomed->Set("ActiveArea",   UnGroomedJets[iJet].area());
	  fastjet::PseudoJet area = UnGroomedJets[iJet].area_4vector();
	  UnGroomed->Set("ActiveArea_px", area.px());
	  UnGroomed->Set("ActiveArea_py", area.py());
	  UnGroomed->Set("ActiveArea_pz", area.pz());
	  UnGroomed->Set("ActiveArea_e" , area.e());
	  UnGroomed->Set("constscale_pt", UnGroomedJets[iJet].pt());

      // get ghostassociated tracks
	  vector<fastjet::PseudoJet> consts = UnGroomedJets[iJet].constituents();
      vector<fastjet::PseudoJet> tracksgoodGhostHS; // vector <PJ> HS tracks for jetCleansing                                                                                                     
      vector<fastjet::PseudoJet> tracksgoodGhostPU; // vector <PJ> PU tracks for jetCleansing             
	  UnGroomed->AddVec(fTrkSel+"Ghost");
      UnGroomed->AddVec(fTrkSel+"GhostPU");
      UnGroomed->AddVec(fTrkSel+"GhostHS");
	  UnGroomed->AddVec("constituents");
	  for (int k=0; k<consts.size(); k++){
	    if (consts[k].user_info<PJ_Info>().GhostTrack){
	      UnGroomed->Add(fTrkSel+"Ghost",consts[k].user_info<PJ_Info>().Pointer);
          if(consts[k].user_info<PJ_Info>().Pointer->Int("origin")==0) {
              UnGroomed->Add(fTrkSel+"GhostHS",consts[k].user_info<PJ_Info>().Pointer);
              tracksgoodGhostHS.push_back(consts[k]); // HS track                   
          }else{
              UnGroomed->Add(fTrkSel+"GhostPU",consts[k].user_info<PJ_Info>().Pointer);
              tracksgoodGhostPU.push_back(consts[k]); // PU track      
          }
	    }else{
	      UnGroomed->Add("constituents",   consts[k].user_info<PJ_Info>().Pointer);
        }
	  }
	  Add("jetsFat10AntiKt", UnGroomed);



    #if fjContribFLAG==1 // preprocessor flag defined in makefile --------------------
    // GetCleanedJet-----------------------------------------------------------
    // JVF cleansing
    for(int fCut =0; fCut<11; ++fCut){
        jvf_cleanser_A.SetTrimming(((float)fCut)/100.);
        fastjet::PseudoJet jvf_cleansed_jet = jvf_cleanser_A( UnGroomedJets[iJet], tracksgoodGhostHS, tracksgoodGhostPU );
        Particle* JVFCleansedJet = new Particle(); 
        JVFCleansedJet->p.SetPtEtaPhiE(jvf_cleansed_jet.pt(), jvf_cleansed_jet.eta(), jvf_cleansed_jet.phi(), jvf_cleansed_jet.e());  
        UnGroomed->AddVec(TString::Format("jet_JVFCleansedTrim0p0%i",fCut));
        UnGroomed->Add(TString::Format("jet_JVFCleansedTrim0p0%i",fCut), JVFCleansedJet);                       
    }
    // linear cleansing with \gamma_0=0.55 and various trimming fCuts
    for(int fCut =0; fCut<11; ++fCut){
        linear_cleanser_A.SetLinearParameters(0.55);   
        linear_cleanser_A.SetTrimming(((float)fCut)/100.);
        fastjet::PseudoJet linear_cleansed_jet = linear_cleanser_A( UnGroomedJets[iJet], tracksgoodGhostHS, tracksgoodGhostPU );
        Particle* LinCleansedJet = new Particle(); 
        LinCleansedJet->p.SetPtEtaPhiE(linear_cleansed_jet.pt(), linear_cleansed_jet.eta(), linear_cleansed_jet.phi(), linear_cleansed_jet.e());  
        UnGroomed->AddVec(TString::Format("jet_LinCleansed0p55Trim0p0%i",fCut));
        UnGroomed->Add(TString::Format("jet_LinCleansed0p55Trim0p0%i",fCut), LinCleansedJet);                       
    }
    //      cout << ">>> ------------- CLeansed --------------- <<< " << endl;           
    //      JVFCleansedJet->Show();                                                     
    //      cout << ">>> ------------- Orig     --------------- <<< " << endl;         
    //      UnGroomed->Show();                                                        
    #endif //--------------------------------------------------------------------------


      // setup new cluster sequence to get kt subjets
	  vector<fastjet::PseudoJet> newparticles = ObjsToPJ("constituents",UnGroomed);
      AddAsGhosts(fTrkSel, &newparticles, true);  // add tracks to constituents of ungroomed jet in ghost mode, so that everything can be clustered
	  fastjet::JetDefinition *JetDefSubjets = new fastjet::JetDefinition(fastjet::kt_algorithm, 0.3, fastjet::Best); 
	  fastjet::GhostedAreaSpec areaSpec(5.0);
	  fastjet::ClusterSequenceArea newClustSeq(newparticles, *JetDefSubjets, areaSpec);
	  vector<fastjet::PseudoJet> subj = fastjet::sorted_by_pt(newClustSeq.inclusive_jets(5.)); // subjets of fastjet
      double rho = Float("Eventshape_rhoKt4LC")>0? Float("Eventshape_rhoKt4LC")/1000.:0.0001;
      fastjet::Subtractor subtractor(rho); // rho for jet area corr (mass affected)

      UnGroomed->AddVec("Subjets");
      for(int iSub=0; iSub<subj.size(); ++iSub){
          Particle* sub = new Particle();
	      sub->Set("ptConstScale",  subj[iSub].pt());
	      sub->Set("mConstScale",   subj[iSub].m());
	      sub->Set("etaConstScale", subj[iSub].eta());
	      fastjet::PseudoJet area = subj[iSub].area_4vector();
	      sub->Set("ActiveArea_px", area.px());
	      sub->Set("ActiveArea_py", area.py());
	      sub->Set("ActiveArea_pz", area.pz());
	      sub->Set("ActiveArea_e" , area.e());
//        fastjet::PseudoJet sub_corr = subtractor(subj[iSub]);
//        sub->p.SetPtEtaPhiE(sub_corr.pt(), sub_corr.eta(),sub_corr.phi(),sub_corr.e()); // PU corrected momentum (also mass)
          sub->p.SetPtEtaPhiE(subj[iSub].pt(),subj[iSub].eta(),subj[iSub].phi(),subj[iSub].e()); // *NOT* area corrected 
          sub->AddVec(fTrkSel+"Ghost");
          sub->AddVec("constituents");
	      vector<fastjet::PseudoJet> constit = subj[iSub].constituents();
	      for (int k=0; k<constit.size(); k++){
	        if (constit[k].user_info<PJ_Info>().GhostTrack){
		      sub->Add(fTrkSel+"Ghost",constit[k].user_info<PJ_Info>().Pointer);
//            constit[k].user_info<PJ_Info>().Pointer->Show();
	        }else{
		      sub->Add("constituents",constit[k].user_info<PJ_Info>().Pointer);
            }
	      }
          UnGroomed->Add("Subjets", sub);
          Add("jetsKt3LCTopoSubjets",sub);
      }

      delete JetDefSubjets;

	}
	delete fJetDef;  
    

    CalculateDR("Fat10AntiKt", "recosW");
    CalculateDR("Fat10AntiKt", "recosZ");
    CalculateDR("Fat10AntiKt", "recosTop");
    CalculateDR("Fat10AntiKt", "recosAntiTop");
    CalculateJVF("Kt3LCTopoSubjets",fTrkSel+"Ghost");
    AddTrimmedJetCollections("Fat10AntiKt");
    MakeJetPlots("Fat10AntiKt");
//    MakeJetPlots("Kt3LCTopoSubjets");
      

    return true;
}

//---------------------------
// AddTrimmedJetCollections
//---------------------------
void Analysis_LargeR::AddTrimmedJetCollections(const MomKey JetKey){


    for(int j=0; j<jets(JetKey); ++j){
        Particle *myjet = &jet(j, JetKey);
        Particle *jet_0p10Trim                  = new Particle();  jet_0p10Trim                 ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p09Trim                  = new Particle();  jet_0p09Trim                 ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p08Trim                  = new Particle();  jet_0p08Trim                 ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p07Trim                  = new Particle();  jet_0p07Trim                 ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p06Trim                  = new Particle();  jet_0p06Trim                 ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p05Trim                  = new Particle();  jet_0p05Trim                 ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p04Trim                  = new Particle();  jet_0p04Trim                 ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p03Trim                  = new Particle();  jet_0p03Trim                 ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p02Trim                  = new Particle();  jet_0p02Trim                 ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p01Trim                  = new Particle();  jet_0p01Trim                 ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p10TrimCorrJVF           = new Particle();  jet_0p10TrimCorrJVF          ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p09TrimCorrJVF           = new Particle();  jet_0p09TrimCorrJVF          ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p08TrimCorrJVF           = new Particle();  jet_0p08TrimCorrJVF          ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p07TrimCorrJVF           = new Particle();  jet_0p07TrimCorrJVF          ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p06TrimCorrJVF           = new Particle();  jet_0p06TrimCorrJVF          ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p05TrimCorrJVF           = new Particle();  jet_0p05TrimCorrJVF          ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p04TrimCorrJVF           = new Particle();  jet_0p04TrimCorrJVF          ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p03TrimCorrJVF           = new Particle();  jet_0p03TrimCorrJVF          ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p02TrimCorrJVF           = new Particle();  jet_0p02TrimCorrJVF          ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p01TrimCorrJVF           = new Particle();  jet_0p01TrimCorrJVF          ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_0p00TrimCorrJVF           = new Particle();  jet_0p00TrimCorrJVF          ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_MyClensing0p01Trim        = new Particle();  jet_MyClensing0p01Trim       ->p.SetPxPyPzE(0,0,0,0);
        Particle *jet_MyClensing                = new Particle();  jet_MyClensing               ->p.SetPxPyPzE(0,0,0,0);


        for(int iSubj=0; iSubj<myjet->Objs("Subjets"); ++iSubj){

            Particle * subj = (Particle*) myjet->Obj("Subjets", iSubj);
            subj->Set("PtFrac", subj->p.Pt() / myjet->p.Pt());

            // Trimming           
            if(subj->Float("PtFrac") > 0.10)                                                      (jet_0p10Trim->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.09)                                                      (jet_0p09Trim->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.08)                                                      (jet_0p08Trim->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.07)                                                      (jet_0p07Trim->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.06)                                                      (jet_0p06Trim->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.05)                                                      (jet_0p05Trim->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.04)                                                      (jet_0p04Trim->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.03)                                                      (jet_0p03Trim->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.02)                                                      (jet_0p02Trim->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.01)                                                      (jet_0p01Trim->p)+= (subj->p);  
            
            // trimming & CorrJVF>0.6 (90% eff)
            if(subj->Float("PtFrac") > 0.10 && subj->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF")>0.6)      (jet_0p10TrimCorrJVF->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.09 && subj->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF")>0.6)      (jet_0p09TrimCorrJVF->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.08 && subj->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF")>0.6)      (jet_0p08TrimCorrJVF->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.07 && subj->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF")>0.6)      (jet_0p07TrimCorrJVF->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.06 && subj->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF")>0.6)      (jet_0p06TrimCorrJVF->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.05 && subj->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF")>0.6)      (jet_0p05TrimCorrJVF->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.04 && subj->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF")>0.6)      (jet_0p04TrimCorrJVF->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.03 && subj->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF")>0.6)      (jet_0p03TrimCorrJVF->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.02 && subj->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF")>0.6)      (jet_0p02TrimCorrJVF->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.01 && subj->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF")>0.6)      (jet_0p01TrimCorrJVF->p)+= (subj->p);  
            if(subj->Float("PtFrac") > 0.00 && subj->Float("JVFLinks_"+fTrkSel+"Ghost_nPUTrkCorrJVF")>0.6)      (jet_0p00TrimCorrJVF->p)+= (subj->p);  



            //My clensing
            float scalefact = (1-subj->Float(fTrkSel+"Ghost_PUTrkSumPt")/(0.55*subj->p.Pt()));
            if(subj->p.Pt() * scalefact / myjet->p.Pt() > 0.00)                                   jet_MyClensing        ->p += (subj->p * scalefact); 
            if(subj->p.Pt() * scalefact / myjet->p.Pt() > 0.01)                                   jet_MyClensing0p01Trim->p += (subj->p * scalefact); 

        }

        myjet->Add("jet_0p10Trim"             , jet_0p10Trim);
        myjet->Add("jet_0p09Trim"             , jet_0p09Trim);
        myjet->Add("jet_0p08Trim"             , jet_0p08Trim);
        myjet->Add("jet_0p07Trim"             , jet_0p07Trim);
        myjet->Add("jet_0p06Trim"             , jet_0p06Trim);
        myjet->Add("jet_0p05Trim"             , jet_0p05Trim);
        myjet->Add("jet_0p04Trim"             , jet_0p04Trim);
        myjet->Add("jet_0p03Trim"             , jet_0p03Trim);
        myjet->Add("jet_0p02Trim"             , jet_0p02Trim);
        myjet->Add("jet_0p01Trim"             , jet_0p01Trim);

        myjet->Add("jet_0p10TrimCorrJVF"      , jet_0p10TrimCorrJVF);
        myjet->Add("jet_0p09TrimCorrJVF"      , jet_0p09TrimCorrJVF);
        myjet->Add("jet_0p08TrimCorrJVF"      , jet_0p08TrimCorrJVF);
        myjet->Add("jet_0p07TrimCorrJVF"      , jet_0p07TrimCorrJVF);
        myjet->Add("jet_0p06TrimCorrJVF"      , jet_0p06TrimCorrJVF);
        myjet->Add("jet_0p05TrimCorrJVF"      , jet_0p05TrimCorrJVF);
        myjet->Add("jet_0p04TrimCorrJVF"      , jet_0p04TrimCorrJVF);
        myjet->Add("jet_0p03TrimCorrJVF"      , jet_0p03TrimCorrJVF);
        myjet->Add("jet_0p02TrimCorrJVF"      , jet_0p02TrimCorrJVF);
        myjet->Add("jet_0p01TrimCorrJVF"      , jet_0p01TrimCorrJVF);
        myjet->Add("jet_0p00TrimCorrJVF"      , jet_0p00TrimCorrJVF);

        myjet->Add("jet_MyClensing"           , jet_MyClensing);
        myjet->Add("jet_MyClensing0p01Trim"   , jet_MyClensing0p01Trim);

//        myjet->Show();
//        ((Particle*)myjet->Obj("jet_JVFCleansedTrim0p01"))->Show();
//        ((Particle*)myjet->Obj("jet_MyClensing"))         ->Show();
//        ((Particle*)myjet->Obj("jet_MyClensing0p01Trim")) ->Show();

    }
}


//--------------------------------
// Make Jet plots
// -----------------------------
void Analysis_LargeR::MakeJetPlots(const MomKey JetKey){

    TDirectory* start_dir = OutputDir()->Dir();
    if( OutputDir()->GetList()->Contains(TString(JetKey) ))
        SetOutputDir( OutputDir()->Dir()->GetDirectory(TString(JetKey))  );
    else
        SetOutputDir( OutputDir()->Dir()->mkdir(TString(JetKey)) );


    TString JetKeyS = JetKey.Data();
    for(int j=0; j<jets(JetKey); ++j){
        if(jet(j, JetKey).p.Pt() < 300       ) continue;
        if(fabs(jet(j, JetKey).p.Eta())>1.5  ) continue;
        Particle *jet_0p05Trim= (Particle*) jet(j, JetKey).Obj("jet_0p05Trim");


        if(j==0){ // Leading Jet
            Fill("J0Pt"   , jet(j, JetKey).p.Pt(),   Weight(), 100, 0, 2000);
            Fill("J0Mass",  jet(j, JetKey).p.M(),    Weight(), 100, 0,  500);
            Fill("J0PtVsUngroomedMassVsNPV",    jet(j, JetKey).p.Pt(), jet(j, JetKey).p.M(),      Int("NPV"), Weight(), 100, 200, 2000, 200, -50, 500, 51, -0.5, 50.5); 
            Fill("J0PtVs_MpT05_Trimmed_vsNPV"             ,jet(j, JetKey).p.Pt(),  jet_0p05Trim->p.M(), Int("NPV"), Weight(), 100, 200, 2000, 200, -50, 500, 51, -0.5, 50.5);
        }
        Fill("JPt"   , jet(j, JetKey).p.Pt(),   Weight(), 100, 0, 2000);
        Fill("JMass",  jet(j, JetKey).p.M(),    Weight(), 100, 0,  500);
    }

    SetOutputDir( start_dir );
}


/// --------------------------------------------------
// calculatedR
/// --------------------------------------------------
void Analysis_LargeR::CalculateDR(const MomKey JetKey, const MomKey reco){
    Particle* Reco   = (Objs(reco)==1)? (Particle*)Obj(reco,0):0;
    if(Reco==0)      {cout << reco << " " << Objs(reco) << endl; return;}

    Particle* child1 = (Reco->Objs("children")>0)? ((Particle*)Reco->Obj("children", 0)):0;
    Particle* child2 = (Reco->Objs("children")>1)? ((Particle*)Reco->Obj("children", 1)):0;
    for(int iJet = 0; iJet < jets(JetKey); iJet++){
        float dR = jet(iJet, JetKey).p.DeltaR(Reco->p);
        jet(iJet, JetKey).Set(reco+"_dR", dR);
        if(child1!=0){
            float child1dR = jet(iJet, JetKey).p.DeltaR(child1->p);
            jet(iJet, JetKey).Set(reco+"_child1_dR", child1dR);
        }else {jet(iJet, JetKey).Set(reco+"_child1_dR", -99); }
        if(child2!=0){
            float child2dR = jet(iJet, JetKey).p.DeltaR(child2->p);
            jet(iJet, JetKey).Set(reco+"_child2_dR", child2dR);
        }else {jet(iJet, JetKey).Set(reco+"_child2_dR", -99); }
    }
    return;
}

///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_LargeR::WorkerTerminate()
{


  // Nothing more

}



///===============================================
/// Calculate JVF 
///===============================================
void Analysis_LargeR::CalculateJVF(MomKey JetType, MomKey TrackType, bool requireCuts){
    for(int iJet = 0; iJet < jets(JetType); iJet++) {
        Particle* j = &(jet(iJet, JetType));
        

        // add scalar track pt sum of all tracks associated to this jet
        MomKey TrackType_SumTrackPt_Key(TrackType+"_SumTrackPt");
        j->Set(TrackType_SumTrackPt_Key, 0.);   

        // add vectorJVF momentObj to jet, this is where JVF will be stored
        MomKey vectorJVFKey("JVFLinks_"+TrackType);
        j->AddVec(vectorJVFKey);
        for(int iVtx = 0; iVtx < vtxs(); iVtx++) {
            MomentObj* myJVF = new MomentObj();
            myJVF->AddVec("tracks"); // added later
            myJVF->AddVec("vertex");
            myJVF->Add("vertex", &(vtx(iVtx)));
            myJVF->Set("JVF",  0.); 
            myJVF->Set("JVPt", 0.); 
            j->Add(vectorJVFKey, myJVF); // add vector JVF for each vertex
        }

        // loop over associated tracks
        if(j->Exists(TrackType)){
            for (int t=0; t < j->Objs(TrackType); ++t){
                Particle* track = (Particle*) (j->Obj(TrackType,t));
                if(requireCuts && (track->Bool("passedCuts")==false)) continue; 
                if( track->Exists("matchedVertex") ){
                    j->Set(TrackType_SumTrackPt_Key, j->Float(TrackType_SumTrackPt_Key)+track->p.Pt()); // add track to jet->trackpt

                    MomentObj* vertex = (MomentObj*)(track->Obj("matchedVertex"));
                    // find track vertex in vectorJVFKey vertices
                    bool vertexMatchOK(false);
                    for (int iJVF=0; iJVF < j->Objs(vectorJVFKey); ++iJVF){
                        MomentObj* JVF        = (MomentObj*)(j->Obj(vectorJVFKey,iJVF)); 
                        MomentObj* JVFvertex  = (MomentObj*)(JVF->Obj("vertex"));
                        if (JVFvertex==vertex) {
                            if(Debug()){ 
                                cout << " match found with vertex " << iJVF  << " track origin " << track->Int("origin") << endl;
                            }
                            JVF->Set("JVF",  JVF->Float("JVF")+track->p.Pt());
                            JVF->Set("JVPt", JVF->Float("JVPt")+track->p.Pt());
                            JVF->Add("tracks", track); // add tracklinks
                            vertexMatchOK=true;
                            break;
                        }
                    } 
                    if(! vertexMatchOK) {
                        cout << "EventBuilder_jetmet2012pileupcustom::CalculateJVF: ERROR: vertexMatchOK=false, there is something wrong. exiting. " << endl;
                        exit(1);
                    }
                } else {if(Debug()) cout << "EventBuilder_jetmet2012pileupcustom::CalculateJVF: WARNING: track " << t << " is not matched to vertex" << endl;}
            }
        }
        // loop over vectorJVF and normalize each JVF
        int nJVFTracks =0; //number of tracks actually used in the JVF calculation
        for(int iJVF=0; iJVF<j->Objs(vectorJVFKey); ++iJVF){
            if( j->Float(TrackType_SumTrackPt_Key) > 0){
                ((MomentObj*)j->Obj(vectorJVFKey,iJVF))->SortPt("tracks"); // sort tracks according to Pt
                float tmpJVF = ((MomentObj*)(j->Obj(vectorJVFKey, iJVF)))->Float("JVF");
                ((MomentObj*)(j->Obj(vectorJVFKey, iJVF)))->Set("JVF",     tmpJVF/(j->Float(TrackType_SumTrackPt_Key)));
            } else {
                ((MomentObj*)(j->Obj(vectorJVFKey, iJVF)))->Set("JVF"     , -1);
            }
            nJVFTracks+=((MomentObj*)(j->Obj(vectorJVFKey, iJVF)))->Objs("tracks");
        }
        j->Set("nJVFTracks"+TrackType, nJVFTracks);
        if(j->Objs(vectorJVFKey)>0){
            float JVFZero = ((MomentObj*)(j->Obj(vectorJVFKey, 0)))->Float("JVF");
            j->Set(TrackType+"_PUTrkSumPt", j->Float(TrackType_SumTrackPt_Key) - JVFZero*(j->Float(TrackType_SumTrackPt_Key)));
            j->Set(TrackType+"_HSTrkSumPt", JVFZero*(j->Float(TrackType_SumTrackPt_Key)));
            j->Set("JVFLinks_"+TrackType+"_JVF", JVFZero);
        }
        j->Set("JVFLinks_"+TrackType+"_nPUTrkCorrJVF",            GetNPUTrkCorrJVF(j, "JVFLinks_"+TrackType, TrackType, 0.01));
    }
}

float Analysis_LargeR::GetNPUTrkCorrJVF(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, float scale){
    // return HS / (HS + PU/(nPUTrk*scale))

    TString   simpleTrackKey = TrackKey.Data(); simpleTrackKey.ReplaceAll("Ghost", "");

    if(! jet->Exists(TrackKey+"_PUTrkSumPt")) return -1;
    if(! jet->Exists(TrackKey+"_HSTrkSumPt")) return -1;
    if(! Exists("nPUTracks_"+simpleTrackKey)) return -1;

    float PUTrkPtSum = jet->Float(TrackKey+"_PUTrkSumPt");
    float HSTrkPtSum = jet->Float(TrackKey+"_HSTrkSumPt");
    int   nPUTrk     = Int("nPUTracks_"+simpleTrackKey);

    if     (HSTrkPtSum <0 || PUTrkPtSum <0  )  return -1;    // this should be the case if JVF == -1
    else if(HSTrkPtSum    +  PUTrkPtSum ==0 )  return -1;
    else   {return HSTrkPtSum / (HSTrkPtSum + PUTrkPtSum/( TMath::Max(1, nPUTrk) * scale));}
}

