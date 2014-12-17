ProofAna
=======================
This repository contains the most important parts of ProofAna code I've been using. 

This ProofAna code is used to:
 * Read NTUP_COMMON D3PDs 
 * Perform a matching of reco jets to truth jets to define pileup and hard-scatter jets 
 * Compute JetVertex Tagging information and related variables: [JetVertexTagger RootCore package](utils/JetVertexTagger)
 * Write out small ROOT ntuples that are later used to produce plots, using e.g. the plotting macros from [JVTConfNotePlots/macros](../JVTConfNotePlots/macros)


## Truth Jet Matching & the definition of hard-scatter and pileup jets


## Jet Vertex Tagger

The JetVertexTagger RootCore package is [here](utils/JetVertexTagger). This also containes Standalone installation instructions. The ATLAS JVT twiki for run-I is [here](https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetVertexTagger).

### Usage within ProofAna

In a ProofAna analysis, the JVT RootCore package can be used after a few simple steps:

#### Check out the package:
In your main ProofAna directory:

```
cd utils
svn co $SVNOFF/Reconstruction/Jet/JetAnalysisTools/JetVertexTagger/tags/JetVertexTagger-00-00-01 JetVertexTagger
```

#### in the Analysis code

Don't forget to add
```
#include "JetVertexTagger/JetVertexTagger.h"
```

and initialize the tool in the Begin() method:

```
void Analysis_PileUpStudies::WorkerBegin(){
...
jvt = new JetVertexTagger(0.2, maindir+"/utils/JetVertexTagger/data/JVTlikelihood_20140805.root");
...
}
```


In the ProcessEvent() method (Event Loop), do something like:
```
        // trk pt and z0 ------------------------------------------------------------
        vector<float> trk_pt, trk_z0_wrtPV;
        for(int it=0; it< tracks(); ++it){
          trk_pt                    .push_back(track(it).p.Pt());
          trk_z0_wrtPV              .push_back(track(it).Float("z0_wrtPV"));
          track(it).Set("JVTindex", it);
        }

        // trk to vtx assoc ---------------------------------------------------------
        vector<vector<int> > vxp_trk_index;
        for(int iv=0; iv<vtxs(); ++iv){
          vector<int> assoc_track_indices;
          for(int it=0; it<vtx(iv).Objs("vtxTracks"); ++it){
              Particle* trk = (Particle*) vtx(iv).Obj("vtxTracks",it);
              assoc_track_indices.push_back(trk->Int("JVTindex"));
          }
          vxp_trk_index.push_back(assoc_track_indices);
        }

        // init JVT ----------------------------------------------------------------
        jvt->init(trk_pt, trk_z0_wrtPV, vxp_trk_index);

        // JVT Info for a given jet: myjet --------------------------------------
        vector<int> assoc_trk_indices;
        for(int it=0; it<myjet->Int("nTrackAssoc"); ++it){
            Particle* trk = (Particle*) myjet->Obj("GhostAssocTrack", it);
            assoc_trk_indices.push_back(trk->Int("JVTindex"));
        }

        bool pass = (*jvt)(myjet->p.Pt(), assoc_trk_indices); 
        myjet->Set("corrJVF_RootCoreJVT", jvt->corrJVF());
        myjet->Set("RpT_RootCoreJVT",     jvt->RpT());
        myjet->Set("JVT_RootCoreJVT",     jvt->JVT());
```
