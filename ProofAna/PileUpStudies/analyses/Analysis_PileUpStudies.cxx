/**************************************************************************
 **
 **   File:         Analysis_PileUpStudies.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      P. Nef
 **
 **************************************************************************/

#define Analysis_PileUpStudies_cxx

#include "Analysis_PileUpStudies.h"
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
#include "JetVertexTagger/JetVertexTagger.h"


///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_PileUpStudies::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_PileUpStudies: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();

  // Specific Config
  fTrkSel = "tracksGoodSel0";
  Cfg()->Get("TrackSel", fTrkSel);

  // config -------------
  fisMuScan   = false;
  ChainCfg()->Get("isMuScan",           fisMuScan);
  

  if (Debug()) cout << "Analysis_PileUpStudies: DEBUG Finish WorkerBegin()" << endl;  
  pucounter=0;
  hscounter=0;
  
  std::string maindir(gSystem->Getenv("PROOFANADIR"));
  if(maindir==".") maindir+="/libProofAna";
  cout << maindir << " is the maindir we work from! " << endl;
  cout << gSystem->Exec(TString("ls "+maindir).Data()) << endl;
  
  jvt = new JetVertexTagger(0.2, maindir+"/utils/JetVertexTagger/data/JVTlikelihood_20140805.root");
} 


///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_PileUpStudies::ProcessEvent()
{

  if (Debug()) cout << "Analysis_PileUpStudies: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  
  OutputDir()->cd();
  
  if(! Exists("PileUpStudiesSelection")) return true;
  if(! Bool("PileUpStudiesSelection")  ) return true;

  // Trk definition ------------------------------------------------------------------------
  MomKey myTracks(fTrkSel);
  MomKey myTracksGhost(fTrkSel+"Ghost");
  MomKey myJVFLinks("JVFLinks_"+fTrkSel+"Ghost");
  TrackSelector(myTracks);

  // Add GhostAssociated Tracks --------------------------------------------------------------
  if(!fisMuScan){
   AddGhostMatch("AntiKt4LCTopo", myTracks);
  }

  // Calculate JVF with GhostAssociated Tracks ------------------------------------------------
  CalculateJVF        ("AntiKt4LCTopo_ptgt20", myTracksGhost);
  CalculateJetVars    ("AntiKt4LCTopo_ptgt20", myJVFLinks, myTracksGhost);
  CalculateJetVarsCalo("AntiKt4LCTopo_ptgt20");



  return true;
}

///==========================================
/// CalculateJetVars -- Tracking
///==========================================
void Analysis_PileUpStudies::CalculateJetVars(const MomKey JetKey, const MomKey JVFKey, const MomKey TrackKey){

    for(int iJet = 0; iJet < jets(JetKey); iJet++){
        Particle  *myjet         = &(jet(iJet, JetKey));
        MomentObj *myJVFzero     = myjet->Objs(JVFKey)>0          ?  (MomentObj*)(myjet->Obj(JVFKey, 0)): NULL;
        Particle* trklead =NULL;
        Particle* trktrail=NULL;
        if(myJVFzero!=NULL){
            trklead        = myJVFzero->Objs("tracks")>0    ?  (Particle*) myJVFzero->Obj("tracks",0) : NULL;
            trktrail       = myJVFzero->Objs("tracks")>1    ?  (Particle*) myJVFzero->Obj("tracks",1) : NULL;
        }
        TString   simpleTrackKey = TrackKey.Data(); simpleTrackKey.ReplaceAll("Ghost", "");

        // trk vars
        myjet->Set(JVFKey+"_JVF",                      myJVFzero!=NULL? myJVFzero->Float("JVF"):-1);
        myjet->Set(JVFKey+"_MaxJVF",                   GetValueOrderedJVF(myjet, JVFKey, 0));
        myjet->Set(JVFKey+"_TrkWidth",                 GetJetPVWidth(myjet, JVFKey, 0));
        myjet->Set(JVFKey+"_TrkMeanDR",                GetJetPVTrkMeanDR(myjet, JVFKey, 0));
        myjet->Set(JVFKey+"_HSPVTrk12dR",              (trktrail > 0 && trklead >0)? trklead ->p.DeltaR(trktrail->p): -0.5);
        myjet->Set(JVFKey+"_HSPVtrkSumOverPt",         (GetJetPVTrkSum(myjet, JVFKey, 0, false) >0)? GetJetPVTrkSum(myjet, JVFKey, 0, false)/myjet->p.Pt(): -0.5);
        myjet->Set(JVFKey+"_HSPVtrkSumOverPtdRWeight", (GetJetPVTrkSum(myjet, JVFKey, 0, true ) >0)? GetJetPVTrkSum(myjet, JVFKey, 0, true )/myjet->p.Pt(): -0.5);
        myjet->Set(JVFKey+"_AverageTrkDeltaZ",         GetJetAverageTrkDeltaZ(myjet, TrackKey, false,false));
        myjet->Set(JVFKey+"_AverageTrkDeltaZlog",      GetJetAverageTrkDeltaZ(myjet, TrackKey, true, false));
        myjet->Set(JVFKey+"_AverageTrkDeltaZsqrt",     GetJetAverageTrkDeltaZ(myjet, TrackKey, false,true));
        myjet->Set(JVFKey+"_leadTrkPVdZ",              GetJetLeadingTrkDeltaZ(myjet, TrackKey, false,false));
        myjet->Set(JVFKey+"_leadTrkPVdZlog",           GetJetLeadingTrkDeltaZ(myjet, TrackKey, true, false));
        myjet->Set(JVFKey+"_leadTrkPVdZsqrt",          GetJetLeadingTrkDeltaZ(myjet, TrackKey, false,true));
        myjet->Set(JVFKey+"_BetaStar",                 GetJetBetaStar(myjet, TrackKey));
        myjet->Set(JVFKey+"_NPVCorrJVF",               GetNPVCorrJVF(myjet, JVFKey, TrackKey, 1));
        myjet->Set(JVFKey+"_MuCorrJVF",                GetMuCorrJVF(myjet, JVFKey, TrackKey, 1));
        myjet->Set(JVFKey+"_MuCorrJVF0p1",             GetMuCorrJVF(myjet, JVFKey, TrackKey, 0.1));
        myjet->Set(JVFKey+"_MuCorrJVF0p14",            GetMuCorrJVF(myjet, JVFKey, TrackKey, 0.14));
        myjet->Set(JVFKey+"_NContribVtxsCorrJVF",      GetNContribVtxCorrJVF(myjet, JVFKey, TrackKey, 1));
        myjet->Set(JVFKey+"_nPUTrkCorrJVF",            GetNPUTrkCorrJVF(myjet, JVFKey, TrackKey, 0.01));
        myjet->Set(JVFKey+"_NContribVtxs",             GetNContribVtxs(myjet, JVFKey));
        myjet->Set(JVFKey+"_PUTrkSumPtSquared",        GetWeightedPUTrkSumPt(myjet, JVFKey, TrackKey, false, true));
        myjet->Set(JVFKey+"_HSTrkSumPtSquared",        GetWeightedHSTrkSumPt(myjet, JVFKey, TrackKey, false, true));
        myjet->Set(JVFKey+"_PUTrkSumPtdRWeight",       GetWeightedPUTrkSumPt(myjet, JVFKey, TrackKey, true, false));
        myjet->Set(JVFKey+"_HSTrkSumPtdRWeight",       GetWeightedHSTrkSumPt(myjet, JVFKey, TrackKey, true, false));
        myjet->Set(JVFKey+"_CorrPUtrkPtSumOverPt",     (myjet->Float(TrackKey+"_PUTrkSumPt")-Float("RhoGoodPUTracks")*myjet->Float("Area"))/myjet->p.Pt());
        myjet->Set(JVFKey+"_Corr2PUtrkPtSumOverPt",    (myjet->Float(TrackKey+"_PUTrkSumPt")-Float("RhoGoodPUTracks_"+simpleTrackKey)*myjet->Float("Area"))/myjet->p.Pt());
        // closeby effects
        int Jindex_dRmin = -1;
        myjet->Set(JVFKey+"_Closeby_kt",               GetClosebyDistance(myjet, JetKey, true,  Jindex_dRmin));
        myjet->Set(JVFKey+"_Closeby_dR",               GetClosebyDistance(myjet, JetKey, false, Jindex_dRmin));
        myjet->Set(JVFKey+"_Closeby_dRindex",          Jindex_dRmin);

/*        if(myjet->Objs("constituents") == 1){
        cout << "Width" << myjet->Float("WIDTH") << " my width" <<  GetConstitPtWeightedMeanOverDr(myjet) << endl;
        cout << "std dev " << GetConstitPtWeightedStdDevOverDr(myjet) << " nclus" << myjet->Objs("constituents") << endl;
        cout << "skewness " << GetConstitPtWeightedSkewnessOverDr(myjet) << " kurtosis " << GetConstitPtWeightedKurtosisOverDr(myjet) << endl;
        }
*/

        // TestRootCore Package ------------------------------------------------------------
        // trk pt and z0
        vector<float> trk_pt, trk_z0_wrtPV;
        for(int it=0; it< tracks(); ++it){
          trk_pt                    .push_back(track(it).p.Pt());
          trk_z0_wrtPV              .push_back(track(it).Float("z0_wrtPV"));
          track(it).Set("JVTindex", it);
        }

        // trk to vtx assoc
        vector<vector<int> > vxp_trk_index;
        for(int iv=0; iv<vtxs(); ++iv){
          vector<int> assoc_track_indices;
          for(int it=0; it<vtx(iv).Objs("vtxTracks"); ++it){
              Particle* trk = (Particle*) vtx(iv).Obj("vtxTracks",it);
              assoc_track_indices.push_back(trk->Int("JVTindex"));
          }
          vxp_trk_index.push_back(assoc_track_indices);
        }

        // JVT
        jvt->init(trk_pt, trk_z0_wrtPV, vxp_trk_index);
//        for(int it=0; it<trk_pt.size(); ++it){
//            cout << "JVT Tracks " << it << " pt " << trk_pt[it] << " " << trk_z0_wrtPV[it] << endl;
//        }
//        for(int iv=0; iv<vxp_trk_index.size(); ++iv){
//            cout << "vxp_trk_index " << iv << endl;
//            vector<int> test = vxp_trk_index[iv];
//            for(int it=0; it<test.size(); ++it){
//                cout << " >>> trk " << test[it] << endl;
//            }
//        }


        vector<int> assoc_trk_indices;
        float sum1=0;
        for(int it=0; it<myjet->Int("nTrackAssoc"); ++it){
            Particle* trk = (Particle*) myjet->Obj("GhostAssocTrack", it);
            assoc_trk_indices.push_back(trk->Int("JVTindex"));
            if(trk->Int("origin") == 0) sum1+=trk->p.Pt();
        }

        float sum2=0;
        for(int it=0; it<myjet->Objs(TrackKey); ++it){
            Particle* trk = (Particle*) myjet->Obj(TrackKey, it);
            if(trk->Int("origin") == 0) sum2+=trk->p.Pt();

        }
        bool pass = (*jvt)(myjet->p.Pt(), assoc_trk_indices); 
//        cout << " ---- " << endl;
//        cout << "New corrJVF " << jvt->corrJVF() << " rpt " << jvt->RpT()  << " hstrksumpt " << jvt->RpT() * myjet->p.Pt()  << endl;
//        cout << "Old corrJVF " << myjet->Float(JVFKey+"_nPUTrkCorrJVF") << " rpt " << myjet->Float(JVFKey+"_HSPVtrkSumOverPt")  << " hstrksumpt " << << endl; 
        myjet->Set("corrJVF_RootCoreJVT", jvt->corrJVF());
        myjet->Set("RpT_RootCoreJVT",     jvt->RpT());
        myjet->Set("JVT_RootCoreJVT",     jvt->JVT());

//        cout << " >>>>> JVT Rpt     " << jvt->RpT()               << " orig " << myjet->Float(JVFKey+"_HSPVtrkSumOverPt") << endl;
//        cout << " >>>>> JVT cJVF    " << jvt->corrJVF()           << " orig " << myjet->Float(JVFKey+"_nPUTrkCorrJVF")    << endl;
//        cout << " >>>>> nPUtrks     " << Int("nPUTracks_tracksGoodSel0") << endl;
//

    }
    return;

}
///==========================================
/// CalculateJetVars -- Tracking
///==========================================
void Analysis_PileUpStudies::CalculateJetVarsCalo(const MomKey JetKey){

    for(int iJet = 0; iJet < jets(JetKey); iJet++){
        Particle  *myjet         = &(jet(iJet, JetKey));

        // calo vars 
        myjet->Set("CaloWidth",                        GetJetCaloWidth(myjet, false));
        myjet->Set("CaloWidth2",                       GetJetCaloWidth(myjet, true));
        myjet->Set("AnnulusPtRatio0to1",               GetJetAnnulusPtRatio(myjet, 0,   0.1));
        myjet->Set("AnnulusPtRatio0to2",               GetJetAnnulusPtRatio(myjet, 0,   0.2));
        myjet->Set("AnnulusPtRatio1to2",               GetJetAnnulusPtRatio(myjet, 0.1, 0.2));
        myjet->Set("AnnulusPtRatio2to4",               GetJetAnnulusPtRatio(myjet, 0.2, 0.4));
        myjet->Set("AverageDeltaRSquared",             GetAverageDeltaRSquared(myjet));
        myjet->Set("NClusters",                        myjet->Exists("constituents")? myjet->Objs("constituents"): -1);
        myjet->Set("NumTowers",                        myjet->Exists("NumTowers")?    myjet->Float("NumTowers"): -1);

    }
    return;

}

//-----------------------------------------------------
// Print Event 
//-----------------------------------------------------
void Analysis_PileUpStudies::PrintEvent(){


    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << "->>> Printint Event " << EventNumber()                             << endl;
    
    cout << "Reco Jets " << endl;
    for(int iJ=0; iJ<jets("AntiKt4LCTopo_ptgt20_etalt2p4"); ++iJ){
        Particle *j = &jet(iJ,"AntiKt4LCTopo_ptgt20_etalt2p4");
        cout << "jet " << iJ << " pt " << j->p.Pt() << " eta " << j->p.Eta() << " phi " << j->p.Phi() 
             <<        " isHS " << j->Bool("isHSJet")  << " isPU " << j->Bool("isPUJet") << endl;
    }
    cout << "Truth jets " << endl;
    for(int iJ=0; iJ<jets("AntiKt4Truth"); ++iJ){
        Particle *j = &jet(iJ,"AntiKt4Truth");
        cout << "jet " << iJ << " pt " << j->p.Pt() << " eta " << j->p.Eta() << " phi " << j->p.Phi()  << endl;
    }
    cout << "********************** event ***************** " << endl;
    Show();

}

///---------------------------------------
/// JVF Helper
///----------------------------------------
float Analysis_PileUpStudies::GetValueOrderedJVF(Particle *jet, const MomKey JVFKey, unsigned int position){
    vector<float> vectorJVF;

    for(int i=0; i<jet->Objs(JVFKey); ++i){
        vectorJVF.push_back(((MomentObj*)(jet->Obj(JVFKey,i)))->Float("JVF"));
    }
    std::sort(vectorJVF.begin(), vectorJVF.end(), std::greater<float>());
    if(position >= vectorJVF.size()){
        return -99.99;
    }else{
        return vectorJVF[position];
    }
}

vector<float> Analysis_PileUpStudies::GetValueOrderedJVF(Particle *jet, const MomKey JVFKey){
    vector<float> vectorJVF;

    for(int i=0; i<jet->Objs(JVFKey); ++i){
        vectorJVF.push_back(((MomentObj*)(jet->Obj(JVFKey,i)))->Float("JVF"));
    }
    std::sort(vectorJVF.begin(), vectorJVF.end(), std::greater<float>());
    return vectorJVF;
}

vector<float> Analysis_PileUpStudies::GetValueOrderedJVPt(Particle *jet, const MomKey JVFKey){
    vector<float> vectorJVPt;

    for(int i=0; i<jet->Objs(JVFKey); ++i){
        vectorJVPt.push_back(((MomentObj*)(jet->Obj(JVFKey,i)))->Float("JVPt"));
    }
    std::sort(vectorJVPt.begin(), vectorJVPt.end(), std::greater<float>());
    return vectorJVPt;
}

float Analysis_PileUpStudies::GetJVPt(Particle *jet, const MomKey JVFKey, const unsigned int position){
    return ((MomentObj*)(jet->Obj(JVFKey,position)))->Float("JVPt");
}

vector<float> Analysis_PileUpStudies::GetJVF(Particle *jet, const MomKey JVFKey){
    vector<float> vectorJVF;
    for(int i=0; i<jet->Objs(JVFKey); ++i){
        vectorJVF.push_back(((MomentObj*)(jet->Obj(JVFKey,i)))->Float("JVF"));
    }
    return vectorJVF;
}

float Analysis_PileUpStudies::GetMeanJVF(Particle* jet, const MomKey JVFKey){ 
    float tmp=0; 
    int counter=0;
    for(int i=0; i<jet->Objs(JVFKey); ++i){
        float tmpJVF = ((MomentObj*)(jet->Obj(JVFKey,i)))->Float("JVF");
        if(tmpJVF>0) {
            tmp    += tmpJVF;
            counter++;
        }
    }
    return tmp/counter;
}

float Analysis_PileUpStudies::GetRMSJVF(Particle* jet, const MomKey JVFKey){ 
    float tmp=0; 
    int counter=0;
    for(int i=0; i<jet->Objs(JVFKey); ++i){
        float tmpJVF = ((MomentObj*)(jet->Obj(JVFKey,i)))->Float("JVF");
        if(tmpJVF>0) {
            tmp    += pow(tmpJVF,2);
            counter++;
        }
    }
    return sqrt(tmp/counter);
}

float Analysis_PileUpStudies::GetNxJVF(Particle* jet, const MomKey JVFKey, float Xvalue){
    // get Nx value of value ordered JVF vector
    // Nx=80: returns number of JVFs needed to sum up in order to get an integral of xvalue=0.8
    vector<float> JVFvec = GetValueOrderedJVF(jet, JVFKey);
    float sum=0;
    if (Xvalue > 1) { cout << "Analysis_PileUpStudies::GetNxJVF: Error Xvalue >1 " << endl; return -1; }
    for (int i=0; i<JVFvec.size(); ++i){
        sum += JVFvec[i];
        if (sum >= Xvalue) {return i+1;}
    }
    return -2;

}


vector<Particle*> Analysis_PileUpStudies::GetJetPVTracks(Particle* jet, const MomKey JVFKey, const int PV){
    vector<Particle*> tracks;

    if (jet->Objs(JVFKey) <= PV){
        if(Debug()) cout << "Analysis_PileUpStudies::GetJetPVTracks: Error: jet->Objs(JVFKey) < PV!" << endl; return tracks;
    }

    MomentObj* JVFObj =(MomentObj*)(jet->Obj(JVFKey,PV));
    for (int itrk=0; itrk <JVFObj->Objs("tracks"); ++itrk){
        tracks.push_back((Particle*) JVFObj->Obj("tracks",itrk));
    }
    return tracks;
}

float Analysis_PileUpStudies::GetJetPVTrkSum(Particle* jet, const MomKey JVFKey, const int PV, const bool dRWeight, const bool useTrackAxis){
    // calculate scalar pt sum of tracks from PV i 
    // optionally, trk pt is multiplied by 
    if (jet->Objs(JVFKey) <= PV){
        if(Debug())  cout << "Analysis_PileUpStudies::GetJetPVWidth: Error: jet->Objs(JVFKey) < PV!" << endl; return -1;
    }
    
    MomentObj* JVFObj =(MomentObj*)(jet->Obj(JVFKey,PV));

    if (! JVFObj->Exists("tracks")) { return -1;}
    if (JVFObj->Objs("tracks") == 0){ return -1; }

    // get axis  
    TLorentzVector axis(0,0,0,0);
    if(useTrackAxis){
        for (int itrk=0; itrk <JVFObj->Objs("tracks"); ++itrk){
            Particle* trk = (Particle*) JVFObj->Obj("tracks",itrk);
            axis += trk->p;
        }
    } else{
        axis = jet->p;
    }

    float ptsum = 0;
    for (int itrk=0; itrk <JVFObj->Objs("tracks"); ++itrk){
        Particle* trk = (Particle*) JVFObj->Obj("tracks",itrk);
        if (! dRWeight){
            ptsum += trk->p.Pt();
        } else{
            ptsum += trk->p.Pt() * fabs(TMath::Cos(trk->p.DeltaR(axis) *TMath::Pi()/(2.*0.4) ));
        }
    }

    return ptsum;
}

float Analysis_PileUpStudies::GetJetPVWidth(Particle* jet, const MomKey JVFKey, const int PV){
    // calculate jet width using the tracks associated to the jet from a given PV
    if (jet->Objs(JVFKey) <= PV){
        if(Debug()) cout << "Analysis_PileUpStudies::GetJetPVWidth: Error: jet->Objs(JVFKey) < PV!" << endl; return -1;
    }
    
    MomentObj* JVFObj =(MomentObj*)(jet->Obj(JVFKey,PV));

//    JVFObj->Show();
    if (JVFObj->Objs("tracks") == 0){ return -1; }
    
    float width = 0;
    float ptsum = 0;
    
    for (int itrk=0; itrk <JVFObj->Objs("tracks"); ++itrk){
        Particle* trk = (Particle*) JVFObj->Obj("tracks",itrk);
        width += trk->p.Pt() * jet->p.DeltaR(trk->p);
        ptsum += trk->p.Pt();
    }
    if (fabs(JVFObj->Float("JVPt") - ptsum )>0.01){
        cout << "Analysis_PileUpStudies::GetJetPVWidth : there's something wrong here..." << endl;
    }

    return width / (ptsum * jet->Float("radius")); 
} 

float Analysis_PileUpStudies::GetJetPVTrkMeanDR(Particle* jet, const MomKey JVFKey, const int PV){
    // calculate mean dR of tracks from a PV wrt the jet axis
    if (jet->Objs(JVFKey) <= PV){
        if(Debug())  cout << "Analysis_PileUpStudies::GetJetPVTrkMeanDR: Error: jet->Objs(JVFKey) < PV!" << endl; return -1;
    }
    
    MomentObj* JVFObj =(MomentObj*)(jet->Obj(JVFKey,PV));

//    JVFObj->Show();
    if (JVFObj->Objs("tracks") == 0){ return -1; }
    
    float meanDR = 0;
    
    for (int itrk=0; itrk <JVFObj->Objs("tracks"); ++itrk){
        Particle* trk = (Particle*) JVFObj->Obj("tracks",itrk);
        meanDR += jet->p.DeltaR(trk->p);
    }
    meanDR = meanDR / JVFObj->Objs("tracks");

    return meanDR;
} 

float Analysis_PileUpStudies::GetJetCaloWidth(Particle *jet, bool squared){
    // calculate width = sum(pTi * DeltaR_i)/(pT*R) using jet constituents (clusters)
    if(! jet->Exists("constituents")) return -0.5;
    float numerator    = 0;
    float numerator2   = 0;
    float denominator  = 0;
    float denominator2 = 0;
    for(int iC=0; iC<jet->Objs("constituents"); ++iC){
        Particle* constit = (Particle*) jet->Obj("constituents", iC);
        numerator    += (constit->p.Pt() * constit->p.DeltaR(jet->p));
        numerator2   += (constit->p.Pt() * constit->p.Pt() * constit->p.DeltaR(jet->p));
        denominator  += (constit->p.Pt());
        denominator2 += (constit->p.Pt() * constit->p.Pt());
    }

    if     (denominator >0 && ! squared) return numerator  / (denominator  * jet->Float("radius"));
    else if(denominator >0 &&   squared) return numerator2 / (denominator2 * jet->Float("radius"));
    else                                 return -0.5;
}

float Analysis_PileUpStudies::GetJetAnnulusPtRatio(Particle *jet, float dRmin, float dRmax){
    // calculate annulus pT ratio = 1/pT Sum(pT_i), for dRmin <= DeltaR(jet, pT_i) < dRmax
    if(! jet->Exists("constituents")) return -0.5;
    float denominator = 0;
    float numerator   = 0;
    for(int iC=0; iC<jet->Objs("constituents"); ++iC){
        Particle* constit = (Particle*) jet->Obj("constituents", iC);
        denominator  += (constit->p.Pt());
        float dR = constit->p.DeltaR(jet->p);
        if(dR >= dRmin && dR < dRmax) numerator += constit->p.Pt();
    }

    if     (denominator >0 ) return numerator  / denominator;
    else                     return -0.5;
}

float Analysis_PileUpStudies::GetAverageDeltaRSquared(Particle *jet){
    // compute <Delta R^2> = Sum(DeltaR_i^2 * pT_i^2)/Sum(pT_i^2)
    if(! jet->Exists("constituents")) return -0.5;
    float denominator = 0;
    float numerator   = 0;
    for(int iC=0; iC<jet->Objs("constituents"); ++iC){
        Particle* constit = (Particle*) jet->Obj("constituents", iC);
        denominator  += (pow(constit->p.Pt(),2));
        numerator    += (pow(constit->p.DeltaR(jet->p),2) * pow(constit->p.Pt(),2));
    }

    if     (denominator >0 ) return numerator  / denominator;
    else                     return -0.5;
}

float Analysis_PileUpStudies::GetJetLeadingTrkDeltaZ(Particle *jet, const MomKey TrackKey, bool log, bool sqrt){
    // return deltaZ(HS PV, leading trk)
    if      (log && sqrt) {cout << "Analysis_PileUpStudies::GetJetLeadingTrkDeltaZ: Warning: log && sqrt" << endl; return -1;}
    if      (! (jet->Exists(TrackKey)) )   return -0.5;
    else if (jet->Objs(TrackKey)==0    )   return -0.5;
    else {
        Particle *trkleadAllPV = (Particle*) jet->Obj(TrackKey,0);
        if     ( log) return TMath::Log10(max(fabs(trkleadAllPV->Float("z0_wrtPV")),(float)1.));
        else if(sqrt) return TMath::Sqrt(fabs(trkleadAllPV->Float("z0_wrtPV")));
        else          return fabs(trkleadAllPV->Float("z0_wrtPV"));
    }
}

float Analysis_PileUpStudies::GetJetAverageTrkDeltaZ(Particle *jet, const MomKey TrackKey, bool log, bool sqrt){
    // return Sum{deltaZ(HS PV, trk_i) * pT_i} / Sum{pT_i}
    if     (log && sqrt) {cout << "Analysis_PileUpStudies::GetJetAverageTrkDeltaZ: Warning: log && sqrt" << endl; return -1;}
    if     (! (jet->Exists(TrackKey)) ) return -0.5;
    float numerator  =0;
    float denominator=0;
    for(int iT=0; iT<jet->Objs(TrackKey); ++iT){
        Particle* trk = (Particle*) jet->Obj(TrackKey, iT);
        if(log){
        numerator    += TMath::Log10(max(fabs(trk->Float("z0_wrtPV")),(float) 1.)) * trk->p.Pt();
        }else if (sqrt){
        numerator    += TMath::Sqrt(max(fabs(trk->Float("z0_wrtPV")), (float) 1.)) * trk->p.Pt();
        }else{
        numerator    += fabs(trk->Float("z0_wrtPV")) * trk->p.Pt();
        }
        denominator  += trk->p.Pt(); 
    }
    if(denominator>0) return numerator / denominator;
    else              return -0.5;
}

float Analysis_PileUpStudies::GetJetBetaStar(Particle *jet, const MomKey TrackKey){

    if     (! (jet->Exists(TrackKey)) ) return -0.5;

    float numerator  =0;
    float denominator=0;
    for(int iT=0; iT<jet->Objs(TrackKey); ++iT){
        Particle* trk = (Particle*) jet->Obj(TrackKey, iT);
        if      (! trk->Exists("origin")) numerator += trk->p.Pt(); // tracks not associated with any VTX
        else if (  trk->Int("origin")!=0) numerator += trk->p.Pt(); // tracks from Vertex!=HS 
        denominator += trk->p.Pt(); 
    }
    
    if (denominator>0) return numerator/denominator;
    else               return -0.5;
}

float Analysis_PileUpStudies::GetMuCorrJVF(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, float scale){
    // return HS / (HS + PU/(mu*scale))

    if(! jet->Exists(TrackKey+"_PUTrkSumPt")) return -0.5;
    if(! jet->Exists(TrackKey+"_HSTrkSumPt")) return -0.5;

    float PUTrkPtSum = jet->Float(TrackKey+"_PUTrkSumPt");
    float HSTrkPtSum = jet->Float(TrackKey+"_HSTrkSumPt");

    if     (HSTrkPtSum <0 || PUTrkPtSum <0  )  return -0.5;    // this should be the case if JVF == -1
    else if(HSTrkPtSum    +  PUTrkPtSum ==0 )  return -0.5;
    else   {return HSTrkPtSum / (HSTrkPtSum + PUTrkPtSum/(Float("averageIntPerXing") *scale));}

}

float Analysis_PileUpStudies::GetNPVCorrJVF(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, float scale){
    // return HS / (HS + PU/(mu*scale))

    if(! jet->Exists(TrackKey+"_PUTrkSumPt")) return -0.5;
    if(! jet->Exists(TrackKey+"_HSTrkSumPt")) return -0.5;

    float PUTrkPtSum = jet->Float(TrackKey+"_PUTrkSumPt");
    float HSTrkPtSum = jet->Float(TrackKey+"_HSTrkSumPt");

    if     (HSTrkPtSum <0 || PUTrkPtSum <0  )  return -0.5;    // this should be the case if JVF == -1
    else if(HSTrkPtSum    +  PUTrkPtSum ==0 )  return -0.5;
    else   {return HSTrkPtSum / (HSTrkPtSum + PUTrkPtSum/(Int("NPV") *scale));}

}

float Analysis_PileUpStudies::GetNPUTrkCorrJVF(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, float scale){
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

float Analysis_PileUpStudies::GetNContribVtxCorrJVF(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, float scale){
    // return HS / (HS + PU/(mu*scale))

    if(! jet->Exists(TrackKey+"_PUTrkSumPt")) return -0.5;
    if(! jet->Exists(TrackKey+"_HSTrkSumPt")) return -0.5;

    float PUTrkPtSum = jet->Float(TrackKey+"_PUTrkSumPt");
    float HSTrkPtSum = jet->Float(TrackKey+"_HSTrkSumPt");
    int   NContribVtx= GetNContribVtxs(jet, JVFKey);


    if     (NContribVtx<0)                     return -0.5;
    else if(NContribVtx ==0)                   return -0.5;
    else if(HSTrkPtSum <0 || PUTrkPtSum <0  )  return -0.5;    // this should be the case if JVF == -1
    else if(HSTrkPtSum    +  PUTrkPtSum ==0 )  return -0.5;
    else   {return HSTrkPtSum / (HSTrkPtSum + PUTrkPtSum/(NContribVtx *scale));}

}


int Analysis_PileUpStudies::GetNContribVtxs(Particle *jet, const MomKey JVFKey){
    // count the number of vertices that contribute to the JVF calculation, i.e. with at least 1 trk associated to the jet

    if(! jet->Exists(JVFKey)) return -1;

    int nContrib =0;
    for (int iJVF=0; iJVF<jet->Objs(JVFKey); ++iJVF){
        MomentObj *myJVF     = (MomentObj*)(jet->Obj(JVFKey, iJVF));
        if (myJVF->Float("JVF") >0) nContrib++;
    }

    return nContrib;
}


float Analysis_PileUpStudies::GetWeightedPUTrkSumPt(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, bool dodRWeight, bool dopTWeight){
    // calculate Weighted HS trum Sum pt
    if(! jet->Exists(JVFKey)) return -1;

    float result = 0;
    for(int iJVF=1; iJVF<jet->Objs(JVFKey); ++iJVF){ // JVF[0] is HS
        MomentObj *myJVF     = (MomentObj*)(jet->Obj(JVFKey, iJVF));  
        for(int iTrk=0; iTrk<myJVF->Objs("tracks"); ++iTrk){
            Particle* trk = (Particle*) myJVF->Obj("tracks", iTrk);
            float dRWeight =  fabs(TMath::Cos(trk->p.DeltaR(jet->p) *TMath::Pi()/(2.*0.4) ));
            float pTWeight =  trk->p.Pt();
            float increment=  trk->p.Pt(); 
            if(dodRWeight)    increment*= dRWeight; 
            if(dopTWeight)    increment*= pTWeight; 
            result += increment;
        }
    }
    return result;
}

float Analysis_PileUpStudies::GetWeightedHSTrkSumPt(Particle *jet, const MomKey JVFKey, const MomKey TrackKey, bool dodRWeight, bool dopTWeight){
    // calculate Weighted HS trum Sum pt
    if(jet->Objs(JVFKey)==0) return -1;

    float result = 0;
    MomentObj *myJVF     = (MomentObj*)(jet->Obj(JVFKey, 0));  
    for(int iTrk=0; iTrk<myJVF->Objs("tracks"); ++iTrk){
        Particle* trk = (Particle*) myJVF->Obj("tracks", iTrk);
        float dRWeight =  fabs(TMath::Cos(trk->p.DeltaR(jet->p) *TMath::Pi()/(2.*0.4) ));
        float pTWeight =  trk->p.Pt();
        float increment=  trk->p.Pt(); 
        if(dodRWeight)    increment*= dRWeight; 
        if(dopTWeight)    increment*= pTWeight; 
        result += increment;
    }
    return result;
}

float Analysis_PileUpStudies::GetClosebyDistance(Particle *myjet, const MomKey JetKey, bool kt_metric, int &index){
    // returns mindR(j, jets) or min(dR(j,j1)*min(pTj,pTj1))  
    float kTmetric =    99999.99;
    float dRmetric =    99999.99;
    for(int iJet = 0; iJet < jets(JetKey); iJet++){
        if (myjet == &jet(iJet, JetKey)) continue;
        float dR = myjet->p.DeltaR(jet(iJet, JetKey).p);
        float kt = myjet->p.DeltaR(jet(iJet, JetKey).p) * TMath::Min(myjet->p.Pt(), jet(iJet, JetKey).p.Pt());

        if (dR < dRmetric) {dRmetric = dR; index = iJet;}
        if (kt < kTmetric) {kTmetric = kt;}
    }
    if (kTmetric == 99999.99) return -1;
    else return kt_metric? kTmetric : dRmetric;
}

float Analysis_PileUpStudies::GetConstitPtWeightedMeanOverDr(Particle* myjet){
    // calculates pt weighted average dr: jet width
    // see http://www.nematrian.com/R.aspx?p=WeightedMomentsAndCumulants

    if (! myjet->Exists("constituents")         ) return -999;
    if (! myjet->Exists("constituents_SumPt")   ) return -999;

    float drsum = 0; float sumweight = 0;
    for(int ic=0; ic<myjet->Objs("constituents"); ++ic){
        Particle* constit = (Particle*) myjet->Obj("constituents", ic);
        float weight = constit->p.Pt()/myjet->Float("constituents_SumPt");

        drsum     += weight * constit->p.DeltaR(myjet->p);
        sumweight += weight;
    }
    return drsum / sumweight;
}

float Analysis_PileUpStudies::GetConstitPtWeightedStdDevOverDr(Particle* myjet){
    // calculates pt weighted standard deviation of dr
    // sigma^2 = sum (w_i (dr_i - <dr>)^2) / sum(w_i), w_i = pt_i / sum(pT), <dr> = pt weighted mean (jet width) 
    // see http://www.nematrian.com/R.aspx?p=WeightedMomentsAndCumulants

    if (! myjet->Exists("constituents")         ) return -999;
    if (! myjet->Exists("constituents_SumPt")   ) return -999;
    if(   myjet->Objs("constituents") ==1       ) return 0; // code below typically return 10^{-10} or something like that rather the 0 for only one cluster
                                                            // so, handle separately. 

    float mean = GetConstitPtWeightedMeanOverDr(myjet);

    float drsum = 0; float sumweight = 0;
    for(int ic=0; ic<myjet->Objs("constituents"); ++ic){
        Particle* constit = (Particle*) myjet->Obj("constituents", ic);

        float weight = constit->p.Pt()/myjet->Float("constituents_SumPt");

        drsum     += weight * pow(constit->p.DeltaR(myjet->p) - mean,2);
        sumweight += weight;
    }
    return sqrt(drsum / sumweight);
}

float Analysis_PileUpStudies::GetConstitPtWeightedSkewnessOverDr(Particle* myjet){
    // calculates pt weighted skewness of dr
    // gamma  = sum (w_i (dr_i - <dr>)^3 / sigma^3) / sum(w_i), w_i = pt_i / sum(pT), <dr> = pt weighted mean (jet width) , sigma = std dev 
    // see http://www.nematrian.com/R.aspx?p=WeightedMomentsAndCumulants

    if (! myjet->Exists("constituents")         ) return -999;
    if (! myjet->Exists("constituents_SumPt")   ) return -999;
    if(   myjet->Objs("constituents") ==1       ) return -10; // only one constit would give Inf as result. 

    float mean   = GetConstitPtWeightedMeanOverDr  (myjet);
    float sigma  = GetConstitPtWeightedStdDevOverDr(myjet);

    float drsum = 0; float sumweight = 0;
    for(int ic=0; ic<myjet->Objs("constituents"); ++ic){
        Particle* constit = (Particle*) myjet->Obj("constituents", ic);

        float weight = constit->p.Pt()/myjet->Float("constituents_SumPt");

        drsum     += weight * pow(constit->p.DeltaR(myjet->p) - mean,3)/pow(sigma,3);
        sumweight += weight;
    }
    return drsum / sumweight;
}

float Analysis_PileUpStudies::GetConstitPtWeightedKurtosisOverDr(Particle* myjet){
    // calculates pt weighted kurtosis of dr
    // kappa  = sum (w_i (dr_i - <dr>)^4 / sigma^4) / sum(w_i) - 3,
    // w_i = pt_i / sum(pT), <dr> = pt weighted mean (jet width) , sigma = std dev 
    // "-3" is so that normal distribution has kurtosis = 0
    // see http://www.nematrian.com/R.aspx?p=WeightedMomentsAndCumulants

    if (! myjet->Exists("constituents")         ) return -999;
    if (! myjet->Exists("constituents_SumPt")   ) return -999;
    if(   myjet->Objs("constituents") ==1       ) return -10; // only one constit would give Inf as result. 

    float mean   = GetConstitPtWeightedMeanOverDr  (myjet);
    float sigma  = GetConstitPtWeightedStdDevOverDr(myjet);

    float drsum = 0; float sumweight = 0;
    for(int ic=0; ic<myjet->Objs("constituents"); ++ic){
        Particle* constit = (Particle*) myjet->Obj("constituents", ic);

        float weight = constit->p.Pt()/myjet->Float("constituents_SumPt");

        drsum     += weight * pow(constit->p.DeltaR(myjet->p) - mean,4)/pow(sigma,4);
        sumweight += weight;
    }
    return drsum / sumweight -3; 
}

///===============================================
/// Calculate JVF 
///===============================================
void Analysis_PileUpStudies::CalculateJVF(MomKey JetType, MomKey TrackType, bool requireCuts){
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
        for(int iJVF=0; iJVF<j->Objs(vectorJVFKey); ++iJVF){
            if( j->Float(TrackType_SumTrackPt_Key) > 0){
                ((MomentObj*)j->Obj(vectorJVFKey,iJVF))->SortPt("tracks"); // sort tracks according to Pt
                float tmpJVF = ((MomentObj*)(j->Obj(vectorJVFKey, iJVF)))->Float("JVF");
                ((MomentObj*)(j->Obj(vectorJVFKey, iJVF)))->Set("JVF",     tmpJVF/(j->Float(TrackType_SumTrackPt_Key)));
            } else {
                ((MomentObj*)(j->Obj(vectorJVFKey, iJVF)))->Set("JVF"     , -1);
            }
        }
        if(j->Objs(vectorJVFKey)>0){
            float JVFZero = ((MomentObj*)(j->Obj(vectorJVFKey, 0)))->Float("JVF");
            j->Set(TrackType+"_PUTrkSumPt", j->Float(TrackType_SumTrackPt_Key) - JVFZero*(j->Float(TrackType_SumTrackPt_Key)));
            j->Set(TrackType+"_HSTrkSumPt", JVFZero*(j->Float(TrackType_SumTrackPt_Key)));
        }
    }
}



///-------------------------------------------------------------------------
// TrackToVertex Associator 
///------------------------------------------------------------------------
void Analysis_PileUpStudies::TrackSelector(MomKey TrkType){
    if(Debug()) cout << "Analysis_PileUpStudies::TrackSelector(): DEBUG Starting track selection..." << endl;

    if(Exists(TrkType)){cout << "Analysis_PileUpStudies::TrackSelector:: TrkType " << TrkType << " already exists. Aborting" << endl; exit(1);}

    float RhoGoodPUTracks    =0;
    int   nPUTracks          =0;
    int   nTrksNotFromPV     =0;
    float RhoGoodPUTracksSel =0;
    int   nPUTracksSel       =0;
    int   nTrksNotFromPVSel  =0;

    AddVec(TrkType);
    TString TrkTypeString = TrkType.Data();
    for(int i = 0; i < tracks(); i++){
        bool good(true);
        track(i).Set("GoodSel0", false); // tmp hack 
        if(TrkTypeString == "tracksGoodSel0"){
            if(fabs(track(i).Float("d0_wrtPV")) > 2.5 )                          {good = false;}
            if(track(i).p.Pt() < 0.5)                                            {good = false;}
            if(fabs(track(i).p.Eta()) > 2.5)                                     {good = false;}
            if(track(i).Int("nPixHits")+track(i).Int("nPixelDeadSensors") < 1)   {good = false;}
            if(track(i).Int("nSCTHits")+track(i).Int("nSCTDeadSensors") < 6)     {good = false;}
            if(track(i).Float("chi2")/((float)track(i).Int("ndof")) > 5.)        {good = false;}
            if(good) track(i).Set("GoodSel0", true); // tmp hack
        }else if (TrkTypeString == "tracksGoodSel1"){
            if(fabs(track(i).Float("d0_wrtPV")) > 2.5 )                          {good = false;}
            if(track(i).p.Pt() < 0.4)                                            {good = false;}
            if(fabs(track(i).p.Eta()) > 2.5)                                     {good = false;}
            if(track(i).Int("nPixHits")+track(i).Int("nPixelDeadSensors") < 1)   {good = false;}
            if(track(i).Int("nSCTHits")+track(i).Int("nSCTDeadSensors") < 6)     {good = false;}
            if(track(i).Float("chi2")/((float)track(i).Int("ndof")) > 5.)        {good = false;}
        }else if (TrkTypeString == "tracksGoodSel2"){
            if(track(i).p.Pt() < 0.4)                                            {good = false;}
            if(fabs(track(i).p.Eta()) > 2.5)                                     {good = false;}
        }else if (TrkTypeString == "tracksGoodSel3"){
            if(fabs(track(i).Float("d0_wrtPV")) > 3 )                            {good = false;}
            if(track(i).p.Pt() < 0.4)                                            {good = false;}
            if(fabs(track(i).p.Eta()) > 2.5)                                     {good = false;}
            if(track(i).Int("nPixHits")+track(i).Int("nPixelDeadSensors") < 1)   {good = false;}
            if(track(i).Int("nSCTHits")+track(i).Int("nSCTDeadSensors") < 4)     {good = false;}
            if(track(i).Float("chi2")/((float)track(i).Int("ndof")) > 5.)        {good = false;}
        }else if (TrkTypeString == "tracksGoodSel4"){
            if(fabs(track(i).Float("d0_wrtPV")) > 2.5 )                          {good = false;}
            if(track(i).p.Pt() < 1)                                              {good = false;}
            if(fabs(track(i).p.Eta()) > 2.5)                                     {good = false;}
            if(track(i).Int("nPixHits")+track(i).Int("nPixelDeadSensors") < 1)   {good = false;}
            if(track(i).Int("nSCTHits")+track(i).Int("nSCTDeadSensors") < 6)     {good = false;}
            if(track(i).Float("chi2")/((float)track(i).Int("ndof")) > 5.)        {good = false;}
        }else if (TrkTypeString == "tracksGoodSelITK"){
            if(fabs(track(i).Float("d0_wrtPV")) > 2.5 )                          {good = false;}
            if(track(i).p.Pt() < 1)                                              {good = false;}
            if(fabs(track(i).p.Eta()) > 2.7)                                     {good = false;}
            if(track(i).Int("nPixHits")+track(i).Int("nPixelDeadSensors") 
               +track(i).Int("nSCTHits")+track(i).Int("nSCTDeadSensors") < 11)   {good = false;}
            if(track(i).Int("nPixHoles") >1                                  )   {good = false;}
        }else if (TrkTypeString == "tracksGoodSelITK5GeV"){
            if(fabs(track(i).Float("d0_wrtPV")) > 2.5 )                          {good = false;}
            if(track(i).p.Pt() < 5)                                              {good = false;}
            if(fabs(track(i).p.Eta()) > 2.7)                                     {good = false;}
            if(track(i).Int("nPixHits")+track(i).Int("nPixelDeadSensors") 
               +track(i).Int("nSCTHits")+track(i).Int("nSCTDeadSensors") < 11)   {good = false;}
            if(track(i).Int("nPixHoles") >1                                  )   {good = false;}
        }else if (TrkTypeString == "tracksGoodSelITK3GeV"){
            if(fabs(track(i).Float("d0_wrtPV")) > 2.5 )                          {good = false;}
            if(track(i).p.Pt() < 3)                                              {good = false;}
            if(fabs(track(i).p.Eta()) > 2.7)                                     {good = false;}
            if(track(i).Int("nPixHits")+track(i).Int("nPixelDeadSensors") 
               +track(i).Int("nSCTHits")+track(i).Int("nSCTDeadSensors") < 11)   {good = false;}
            if(track(i).Int("nPixHoles") >1                                  )   {good = false;}
        }else if (TrkTypeString == "tracksGoodSelITK2GeV"){
            if(fabs(track(i).Float("d0_wrtPV")) > 2.5 )                          {good = false;}
            if(track(i).p.Pt() < 2)                                              {good = false;}
            if(fabs(track(i).p.Eta()) > 2.7)                                     {good = false;}
            if(track(i).Int("nPixHits")+track(i).Int("nPixelDeadSensors") 
               +track(i).Int("nSCTHits")+track(i).Int("nSCTDeadSensors") < 11)   {good = false;}
            if(track(i).Int("nPixHoles") >1                                  )   {good = false;}
        }else{
            cout << "Analysis_PileUpStudies::TrackSelector: unrecognized Track Selection " << TrkTypeString << ".   Aborting" << endl; exit(1);  
        }

        if(good)  Add(TrkType, &(track(i)));

        if ( (track(i).Int("origin")>0) && track(i).p.Perp()< 30)  {
            RhoGoodPUTracks +=track(i).p.Pt(); // tracks attached to PU vertices only
            nPUTracks       +=1;
            if(good){
                RhoGoodPUTracksSel +=track(i).p.Pt(); // tracks attached to PU vertices only
                nPUTracksSel       +=1;
            }  
        }
        if ( (track(i).Int("origin")!=0)   && track(i).p.Perp()< 30)  {
            nTrksNotFromPV       +=1;
            if(good){
                nTrksNotFromPVSel      +=1;
            }  
        }
    }

    Set("RhoGoodPUTracks",           RhoGoodPUTracks     /(2*TMath::Pi()*2.5*2));
    Set("RhoGoodPUTracks_"+TrkType,  RhoGoodPUTracksSel  /(2*TMath::Pi()*2.5*2));
    Set("nPUTracks",                 nPUTracks);
    Set("nPUTracks_"+TrkType,        nPUTracksSel);
    Set("nTrksNotFromPV",            nTrksNotFromPV);
    Set("nTrksNotFromPV_"+TrkType,   nTrksNotFromPVSel);

    if(Debug()) cout << "Analysis_PileUpStudies: TrackSelector End " << endl;
    return;
}

///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_PileUpStudies::WorkerTerminate()
{
    cout << "pu counter " << pucounter << " hs counter " << hscounter << endl;

  // Nothing more

}
