// JetVertexTagger package
// Questions / Comments: pascal.nef@cern.ch, sch@slac.stanford.edu, Andreas.Wildauer@cern.ch



#ifndef _JETVERTEXTAGGER_
#define _JETVERTEXTAGGER_


#include "TNamed.h"
#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include <iostream>
#include <vector>
#include <math.h> 

class JVTTrack;

//----------------------------------------------------
/// \class JetVertexTagger
/// Class for calculation of the JVT variable to tag 
/// and suppress pileup jets.
///
/// See ATLAS-CONF-2014-018 for details. 
//----------------------------------------------------
class JetVertexTagger : public TNamed{
    public:

        // constructor
        JetVertexTagger();
        JetVertexTagger(const float &jvt_cut, const TString filein);
        // destructor
        ~JetVertexTagger();

        // init method: call this for each new event
        // WARNING: trk_pt is in GeV, NOT MeV
        void init(const std::vector<float> &trk_pt, const std::vector<float> &trk_z0_wrtPV, const std::vector<std::vector<int> > &vxp_trk_index);

        // actual worker
        // WARNING: jet_pt is fully calibrated and in GeV, NOT MeV
        bool operator()(const float jet_pt, const std::vector<int> assoc_trk_indices);

        // setters 
        inline void setDebug(bool debug){
            debug_ = debug;
        }

        // helpers and info
        float corrJVF(){return corrJVF_;}
        float JVT()    {return JVT_;}
        float RpT()    {return RpT_;}


    private: 
        bool                    debug_;
        bool                    goodfile_;
        bool                    goodhisto_;
        TFile*                  file_;
        TH2F*                   likelihood_;
        TString                 likelihood_name_;
        std::vector<JVTTrack*>  trks_;
        float                   jvt_cut_;
        int                     tot_n_pu_tracks_;
        float                   corrJVF_;
        float                   RpT_;
        float                   JVT_;

        // consts
        const float k_z0_cut_;
        const float k_corrJVF_scale_fact_;
        
        #ifdef JETVERTEXTAGGER_STANDALONE
        ClassDef(JetVertexTagger,1);
        #endif
};

#endif

//----------------------------------------------------
/// \class JVTTrack
/// Helper class for tracks used in JetVertexTagger
//----------------------------------------------------
#ifndef _JVTTrack_
#define _JVTTrack_

#include "TNamed.h"
class JVTTrack : public TNamed  {
    public:
        JVTTrack()                        : pt_(0) ,   z0_wrtPV_(-9999),    origin_(-1){};
        JVTTrack(float pt, float z0_wrtPV): pt_(pt) ,  z0_wrtPV_(z0_wrtPV), origin_(-1){};

        float pt_;
        float z0_wrtPV_;
        float origin_;
        
        #ifdef JETVERTEXTAGGER_STANDALONE
        ClassDef(JVTTrack,1);
        #endif
};
#endif

