/**************************************************************************
 **
 **   File:         EventBuilder_jetmet2012pileupcustom
 **
 **   Description:  EventBuilder custom extension 
 **                 
 **
 **   Author:       P. Nef
 **
 **   Created:      4/12/2013
 **   Modified:
 **
 **************************************************************************/
#define EventBuilder_jetmet2012pileupcustom_cxx

#include "EventBuilder_jetmet2012pileupcustom.h"
#include "AnaConfig.h"
#include "AnalysisBase.h"
#include "TrigRootAnalysis/Conditions.h"
#include "TRandom3.h"

///----------------------------------------------
/// DoCustomSetup is called from within EventBuilder_jetmet2012::CopyEvent()
///  -> calles all methods to do calculations on top of EventBuilder_jetmet2012
///----------------------------------------------
bool EventBuilder_jetmet2012pileupcustom::DoCustomSetup(){
    // copy selection cuts
    fEvt->Set("MAXJETTRUTHMATCHDR"  , (float) fEvt->Cfg()->Float("MAXJETTRUTHMATCHDR"));
    fEvt->Set("MINJETPUDR"          , (float) fEvt->Cfg()->Float("MINJETPUDR"));
    fEvt->Set("doSMWZfixes"         , (float) fEvt->Cfg()->Float("DOSMWZfixes"));
    fEvt->Set("doCOMMONfixes"       , (float) fEvt->Cfg()->Float("DOCOMMONfixes"));
    fEvt->Set("doTrkFromVxt",         (bool)  fEvt->Cfg()->Bool("DOTRKFROMVTX"));
  
    fRecoverBtracks  = false;
    fEvt->ChainCfg()->Get("RecoverBtracks",    fRecoverBtracks);


    // select tracks and add vertexlinks
	if (doTrack)  {
        TrackVertexLinks();
        if( doTruth) DoTrackTruthMatching();
    }


    // setup jets
    if(doJet4 && doLCCluster){
        const static MomKey JetKey4LC("AntiKt4LCTopo"); 
        const static MomKey TruthJetKey4("AntiKt4Truth");
        const static MomKey OOTTruthJet4Key("OutOfTimeAntiKt4Truth");
        
        // select jets
        JetSelector(JetKey4LC);
        
        // truth match
        if(doTruthJets){
            AddTruthMatch(JetKey4LC,TruthJetKey4, true, true);
            if (fEvt->Exists("jetsOutOfTimeAntiKt4Truth")) AddTruthMatch(JetKey4LC,OOTTruthJet4Key, false, true);
            if (fEvt->Exists("jetsInTimeAntiKt4Truth"))    AddTruthMatch(JetKey4LC,"InTimeAntiKt4Truth", false, true);
        }
        PartonFlavorMatch(JetKey4LC);
    } 
    return true;
}


///----------------------------------------------
/// AddTruthMatch
///  -> matches jets to truth jets
///  -> if closest truth jet is withing dR<MAXJETTRUTHMATCHDR, truth jet link is added
//   -> if dR > MINJETPUDR, jet is flagged as PileUp Jet
///----------------------------------------------
void EventBuilder_jetmet2012pileupcustom::AddTruthMatch(const MomKey JetType, const MomKey TruthJetType, const bool doLabel, bool highestPt){
    if(! doTruthJets) {
        cout << "EventBuilder_jetmet2012pileupcustom::AddTruthMatch: WARNING: cannot do truth matching without doTruthJets=true" << endl;
        return;
    }
    for(int iJet = 0; iJet < fEvt->jets(JetType); iJet++) {
        Particle* jet = &(fEvt->jet(iJet, JetType));

        // defaults
        if(doLabel) jet->Set("isPUJet", false);
        if(doLabel) jet->Set("isHSJet", false);  
        if(doLabel) jet->Set("dRTruthPt4",  999.);  
        if(doLabel) jet->Set("dRTruthPt6",  999.);  
        if(doLabel) jet->Set("dRTruthPt8",  999.);  
        if(doLabel) jet->Set("dRTruthPt10", 999.);  
        
        if(! fEvt->isMC()){ // data ! -----------------------------
            continue;
        }

        float mindR= 999.99; 
        float maxPt=-999.99; 
        int minDRindex =-1;
        int maxPtIndex =-1;
        for(int iTrueJ=0; iTrueJ < fEvt->jets(TruthJetType); ++iTrueJ){
            Particle* trueJ = &(fEvt->jet(iTrueJ, TruthJetType));

            float dR = (jet->p).DeltaR(trueJ->p);
            if(dR<mindR)                                                            { mindR = dR;            minDRindex = iTrueJ;}
            if(dR < (fEvt->Float("MAXJETTRUTHMATCHDR")) && maxPt < trueJ->p.Pt())   { maxPt = trueJ->p.Pt(); maxPtIndex = iTrueJ;} 
            if(doLabel){
                // get dR to truth jet with certain pT threshold
                if(trueJ->p.Pt()>4  && dR < jet->Float("dRTruthPt4" ))  jet->Set("dRTruthPt4" ,dR);
                if(trueJ->p.Pt()>6  && dR < jet->Float("dRTruthPt6" ))  jet->Set("dRTruthPt6" ,dR);
                if(trueJ->p.Pt()>8  && dR < jet->Float("dRTruthPt8" ))  jet->Set("dRTruthPt8" ,dR);
                if(trueJ->p.Pt()>10 && dR < jet->Float("dRTruthPt10"))  jet->Set("dRTruthPt10",dR);
            }
        }
        if(highestPt){
            if(maxPtIndex != -1){
                jet->Add(TruthJetType+"_match", &(fEvt->jet(maxPtIndex, TruthJetType)));
                if(doLabel) jet->Set("isHSJet", true);  
                if(doLabel) jet->Set("isPUJet", false); 
                if (!(fEvt->jet(maxPtIndex, TruthJetType).Exists(JetType+"_match"))) fEvt->jet(maxPtIndex, TruthJetType).Add(JetType+"_match", jet, true);
            } else if (mindR> (fEvt->Float("MINJETPUDR"))){
                if(doLabel) jet->Set("isPUJet", true);
                if(doLabel) jet->Set("isHSJet", false);  
            } else {
                if(doLabel) jet->Set("isPUJet", false);
                if(doLabel) jet->Set("isHSJet", false);  
            }
        } else {
            if(mindR< (fEvt->Float("MAXJETTRUTHMATCHDR")) && minDRindex!=-1) {
                jet->Add(TruthJetType+"_match", &(fEvt->jet(minDRindex, TruthJetType)));
                if(doLabel) jet->Set("isHSJet", true);  
                if(doLabel) jet->Set("isPUJet", false); 
                if(!(fEvt->jet(minDRindex, TruthJetType).Exists(JetType+"_match"))) fEvt->jet(minDRindex, TruthJetType).Add(JetType+"_match", jet, true);
            } else if (mindR> (fEvt->Float("MINJETPUDR")) && minDRindex!=-1){
                if(doLabel) jet->Set("isPUJet", true);
                if(doLabel) jet->Set("isHSJet", false);  
            } else {
                if(doLabel) jet->Set("isPUJet", false);
                if(doLabel) jet->Set("isHSJet", false);  
            }
        }
    }
}

void EventBuilder_jetmet2012pileupcustom::PartonFlavorMatch(const MomKey JetType){
    // match truth particles to jets within dR < 0.3. 
    // determine flavor based on pdg_id of highest pT matched parton

    for(int iJet = 0; iJet < fEvt->jets(JetType); iJet++) {
        Particle* jet = &(fEvt->jet(iJet, JetType));
        // default values
        jet->Set("Parton_uds"  , false);
        jet->Set("Parton_c"    , false);
        jet->Set("Parton_b"    , false);
        jet->Set("Parton_g"    , false);

        if(! fEvt->isMC())    continue; // do nothing for data
        double maxPt      = -999;
        int    maxPtPDGid = -999;
        for(int iT=0; iT<fEvt->truths(); ++iT){
            if (fEvt->truth(iT).p.Pt() < 1)                                              continue;
            if (fEvt->truth(iT).p.DeltaR(jet->p)>0.3 )                                   continue;
            if (abs(fEvt->truth(iT).Int("pdgId"))>5 && fEvt->truth(iT).Int("pdgId")!=21) continue;
            
            if ( fEvt->truth(iT).p.Pt() > maxPt) {
                maxPt      = fEvt->truth(iT).p.Pt(); 
                maxPtPDGid = abs(fEvt->truth(iT).Int("pdgId"));
            }
        }
        if(maxPt > -999){
            if     (maxPtPDGid == 1 || maxPtPDGid == 2 || maxPtPDGid == 3) {jet->Set("Parton_uds"  , true);}
            else if(maxPtPDGid == 4                                      ) {jet->Set("Parton_c"    , true);}
            else if(maxPtPDGid == 5                                      ) {jet->Set("Parton_b"    , true);}
            else if(maxPtPDGid ==21                                      ) {jet->Set("Parton_g"    , true);}
        }
    }
    return;
}

void EventBuilder_jetmet2012pileupcustom::DoTrackTruthMatching(){
    if(fEvt->Debug()) cout << "Begin EventBuilder_jetmet2012pileupcustom::DoTrackTruthMatching " << endl;
    if(! fEvt->isMC()){return;}
    for(int iTrk=0; iTrk<fEvt->tracks(); ++iTrk){
        Particle *truthmatch =0;
        float minDR=9999.99;
        for(int iTr=0; iTr<fEvt->truths(); ++iTr){
            if(fEvt->truth(iTr).Int("charge")==0 || fEvt->truth(iTr).p.Pt()<0.001) continue;
            float dR=fEvt->truth(iTr).p.DeltaR(fEvt->track(iTrk).p);
            if (dR < minDR) {minDR = dR; truthmatch = &(fEvt->truth(iTr));}
        }
        (fEvt->track(iTrk)).AddVec("truthmatch");
        (fEvt->track(iTrk)).Add("truthmatch", truthmatch);
    }
    if(fEvt->Debug()) cout << "End EventBuilder_jetmet2012pileupcustom::DoTrackTruthMatching " << endl;
}

///----------------------------------------------
/// JetSelector
///  -> labels jet as isGood if basic kinematic selection is satisfied
///----------------------------------------------
void EventBuilder_jetmet2012pileupcustom::JetSelector(const MomKey JetType){
    for(int iJet = 0; iJet < fEvt->jets(JetType); iJet++) {
        Particle* jet = &(fEvt->jet(iJet, JetType));
        
        // add moments
        if( JetIsIsolated(JetType, jet, 1.) ){
            jet->Set("isIsolated_dR1", true);
        } else {
            jet->Set("isIsolated_dR1", false);
        }

        TString type(JetType);
        if      ( type.First("4") !=-1){jet->Set("radius", 0.4);}
        else if ( type.First("6") !=-1){jet->Set("radius", 0.6);}
        else    { cout << "EventBuilder_jetmet2012pileupcustom::JetSelector: could not determine jet radius" << endl;}

    }
}


///----------------------------------------------
/// TrackSelector
///  -> labells good tracks 
//   -> adds vertexLinks to tracks
///----------------------------------------------
Bool_t EventBuilder_jetmet2012pileupcustom::TrackVertexLinks(){
    if(fEvt->Debug()) cout << "EventBuilder_jetmet2012pileupcustom::TrackVertexLinks(): DEBUG Starting track selection..." << endl;


    for(int i = 0; i < fEvt->tracks(); i++){

        // add vertex links to tracks
        fEvt->track(i).Set("origin", -1);
        if (fEvt->Bool("doTrkFromVxt")){
            for(int iVtx = 0; iVtx < fEvt->vtxs(); iVtx++) { 
                // method 2 ---------------
                if(fEvt->track(i).Int("origin")!=-1) break;
                for(int iTrk=0; iTrk<fEvt->vtx(iVtx).Objs("vtxTracks"); ++iTrk){
                    if(&(fEvt->track(i)) == (fEvt->vtx(iVtx).Obj("vtxTracks",iTrk))){
                        fEvt->track(i).Add("matchedVertex", &(fEvt->vtx(iVtx)));
                        fEvt->track(i).Set("origin",        iVtx);
                        break;
                    }
                }
            }
        }else{
            for(int iVtx = 0; iVtx < fEvt->vtxs(); iVtx++) { 
                if(fEvt->track(i).Int("origin")!=-1) break;
                // method 1 -----
                float dzSinTheta = 999.99;
                if (fEvt->Bool("doSMWZfixes") || fEvt->Bool("doCOMMONfixes")){ // nasty fix
                    dzSinTheta = fEvt->track(i).Float("z0_wrtPV")+fEvt->vtx(0).x.z() - fEvt->vtx(iVtx).x.z();
                }else{
                    dzSinTheta = fEvt->track(i).Float("z0") - fEvt->vtx(iVtx).x.z();
                }
                dzSinTheta *= TMath::Sin(fEvt->track(i).Float("theta"));
                if(fabs(dzSinTheta) < 1.0 && fEvt->track(i).Int("origin") == -1) { // match found, if fabs(dzSinTheta) < 1.0 satisfied for multiple vertices, first match is taken
                    fEvt->track(i).Add("matchedVertex", &(fEvt->vtx(iVtx)));
                    fEvt->track(i).Set("origin", iVtx);
                }
                // end method 1 ------------------------------
            }
       }
       // add tracks to HS PV if origin == -1 and |dz sin(theta)|<3 mm  (Qi selection)
       if(fRecoverBtracks && fEvt->track(i).Int("origin") == -1 && fEvt->vtxs()>0){
           float dzSinTheta = fEvt->track(i).Float("z0_wrtPV");
           if(fabs(dzSinTheta) < 3) {
                fEvt->track(i).Add("matchedVertex", &(fEvt->vtx(0)));
                fEvt->track(i).Set("origin", 0);
           }
       }

    }

    if(fEvt->Debug()) cout << "EventBuilder_jetmet2012pileupcustom ::TrackVertexLinks End " << endl;
    return kTRUE;
}

///----------------------------------------------
/// JetIsIsolated
///  -> returns true if no jet of type JetType is within DeltaR < minDR
///----------------------------------------------
bool EventBuilder_jetmet2012pileupcustom::JetIsIsolated( const MomKey JetType, Particle* jet,  float minDR){
    for(int iJet = 0; iJet < fEvt->jets(JetType); iJet++) {
        if(&(fEvt->jet(iJet, JetType)) == jet  ) continue;   
        if(fEvt->jet(iJet, JetType).p.Pt() < 20) continue;
        float dR = jet->p.DeltaR((fEvt->jet(iJet, JetType)).p);
        if(dR < minDR){
            return false;
        }
    }
    return true;
}

