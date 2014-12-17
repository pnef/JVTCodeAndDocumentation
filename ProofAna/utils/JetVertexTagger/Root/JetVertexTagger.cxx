#include "JetVertexTagger/JetVertexTagger.h"
#include <algorithm>    // std::find


// default Constructor
JetVertexTagger::JetVertexTagger()
    : debug_(false), goodfile_(true), goodhisto_(true),
      likelihood_name_("JVTRootCore_kNN100trim_pt20to50_Likelihood"), jvt_cut_(-1),
      tot_n_pu_tracks_(-1), corrJVF_(-99.), RpT_(-99.), JVT_(-99.),
      k_z0_cut_(3.), k_corrJVF_scale_fact_(0.01){
    
      // nothing to be done here. This is a dummy constructor only. 
}

// Constructor 
JetVertexTagger::JetVertexTagger(const float &jvt_cut, const TString filein)
    : debug_(false), goodfile_(true), goodhisto_(true),
      likelihood_name_("JVTRootCore_kNN100trim_pt20to50_Likelihood"), jvt_cut_(jvt_cut),
      tot_n_pu_tracks_(-1), corrJVF_(-99.), RpT_(-99.), JVT_(-99.),
      k_z0_cut_(3.), k_corrJVF_scale_fact_(0.01){

      // load JVT likelihood
      file_ = new TFile(filein, "open");
      if(! file_->IsOpen()){
          std::cerr << "JetVertexTagger::JetVertexTagger: cannot open file " << filein << std::endl;
          goodfile_ = false;
      }
      likelihood_ = (TH2F*) file_->Get(likelihood_name_);
      if(likelihood_ == 0){
          std::cerr << "JetVertexTagger::JetVertexTagger: can't find TH2F " << likelihood_name_ << std::endl;
          goodhisto_ = false;
      }

}


// init method: to be called event-by-events
// WARNING: trk_pt is in GeV!! 
void JetVertexTagger::init(const std::vector<float> &trk_pt, const std::vector<float> &trk_z0_wrtPV, const std::vector<std::vector<int> > &vxp_trk_index)
{

    // clear tot_n_pu_tracks_
    tot_n_pu_tracks_ = 0;

    // clear tracks for new events 
    for(int iT=0; iT<trks_.size(); ++iT){ delete trks_[iT];}
    trks_.clear();

    // do track to vertex association
    // step 1: use association from vtx reconstruction
    // step 2: use |z0|<3 mm requirement for tracks not assiciated to any vertex after step 1
    for(int iT=0; iT<trk_pt.size(); ++iT){
        JVTTrack *trk = new JVTTrack(trk_pt[iT], trk_z0_wrtPV[iT]);
        trks_.push_back(trk);
        
        // step 1
        // loop over vertices and see if trak is associated to vtx 0 (HS) or any other vertex (PU)
        for(int iVtx=0; iVtx<vxp_trk_index.size(); ++iVtx){
            if(std::find(vxp_trk_index[iVtx].begin(), vxp_trk_index[iVtx].end(), iT) != vxp_trk_index[iVtx].end()) {
                if(iVtx == 0) trk->origin_ = 0;
                else          trk->origin_ = 1;
                break; 
            }
        }

        // step 2
        if(trk->origin_ >=0   ) continue; // look at track not associated with any vertex
        if(fabs(trk->z0_wrtPV_) < k_z0_cut_ && vxp_trk_index.size()>0) trk->origin_ = 0;

    }


    // count total number of pileup tracks in the events. 
    // this is needed for corrJVF
    for(int iT=0; iT<trks_.size(); ++iT){
        if( trks_[iT]->origin_ == 1 && trks_[iT]->pt_ < 30 ) tot_n_pu_tracks_++;
    }

    


}

// Destructor 
JetVertexTagger::~JetVertexTagger(){
    // delete tfile
    delete file_;
}


// This is where the magic happens
// WARNING: jet_pt is in GeV! 
bool JetVertexTagger::operator()(const float jet_pt, const std::vector<int> assoc_trk_indices){
    
    // sanity checks
    if(! (goodfile_ && goodhisto_)) {
        std::cerr << "JetVertexTagger::operator(): setup invalid. check input file and histogram. returning false" << std::endl;
        return false;
    }

    // calculate scalar pT sum of associated HS and PU tracks
    float hs_trk_pt_sum = 0;
    float pu_trk_pt_sum = 0;
    for(int i=0; i<assoc_trk_indices.size(); ++i){
        int trk_index = assoc_trk_indices[i];
        if(trk_index >= trks_.size() || trk_index < 0) {
            std::cerr << "JetVertexTagger::operator(): received invalid track index. ignoring track." << std::endl; 
            continue;
        }

        if     (trks_[trk_index]->origin_ == 0) hs_trk_pt_sum += trks_[trk_index]->pt_; // HS
        else if(trks_[trk_index]->origin_ == 1) pu_trk_pt_sum += trks_[trk_index]->pt_; // PU 

    }

    // calculate corrJVF 
    corrJVF_ =  -1;
    if(hs_trk_pt_sum + pu_trk_pt_sum >0){
        corrJVF_ = hs_trk_pt_sum/ (hs_trk_pt_sum+ (pu_trk_pt_sum/ (std::max(1, tot_n_pu_tracks_) * k_corrJVF_scale_fact_) ) );
    }

    // calculate RpT
    RpT_ = hs_trk_pt_sum / jet_pt;

    // calculate JVT
    if      (corrJVF_ == -1.) JVT_ = -0.1;        // no assocaited tracks
    else if (corrJVF_ ==  0.) JVT_ =   0.;         // no hard-scatter tracks
    else {
        float RpT_forJVT = std::min(RpT_, (float) 1.);         // max RpT
        int bin = likelihood_->FindBin(corrJVF_, RpT_forJVT);
        JVT_    = likelihood_->GetBinContent(bin);
//      JVT_    = likelihood_->Interpolate(corrJVF_, RpT_forJVT); // interpolation
    }

    // return result
    if(JVT_ > jvt_cut_) return true;
    else                return false;
}


