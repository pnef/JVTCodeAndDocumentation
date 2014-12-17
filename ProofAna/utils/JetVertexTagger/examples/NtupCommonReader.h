//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 27 16:14:36 2014 by ROOT version 5.34/05
// from TTree physics/physics
// found on file: /atlas/dq2//mc12_8TeV/NTUP_COMMON/e1434_s1499_s1504_r3658_r3549_p1562/mc12_8TeV.147771.Sherpa_CT10_Zmumu.merge.NTUP_COMMON.e1434_s1499_s1504_r3658_r3549_p1562_tid01313919_00/NTUP_COMMON.01313919._000001.root.1
//////////////////////////////////////////////////////////

#ifndef NtupCommonReader_h
#define NtupCommonReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class NtupCommonReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EventNumber;
   Int_t           jet_AntiKt4LCTopo_n;
   vector<float>   *jet_AntiKt4LCTopo_E;
   vector<float>   *jet_AntiKt4LCTopo_pt;
   vector<float>   *jet_AntiKt4LCTopo_m;
   vector<float>   *jet_AntiKt4LCTopo_eta;
   vector<float>   *jet_AntiKt4LCTopo_phi;
   vector<float>   *jet_AntiKt4LCTopo_EtaOrigin;
   vector<float>   *jet_AntiKt4LCTopo_PhiOrigin;
   vector<float>   *jet_AntiKt4LCTopo_MOrigin;
   vector<float>   *jet_AntiKt4LCTopo_WIDTH;
   vector<float>   *jet_AntiKt4LCTopo_Timing;
   vector<float>   *jet_AntiKt4LCTopo_LArQuality;
   vector<float>   *jet_AntiKt4LCTopo_OriginIndex;
   vector<float>   *jet_AntiKt4LCTopo_HECQuality;
   vector<float>   *jet_AntiKt4LCTopo_NegativeE;
   vector<float>   *jet_AntiKt4LCTopo_AverageLArQF;
   vector<float>   *jet_AntiKt4LCTopo_BCH_CORR_CELL;
   vector<float>   *jet_AntiKt4LCTopo_BCH_CORR_DOTX;
   vector<float>   *jet_AntiKt4LCTopo_BCH_CORR_JET;
   vector<float>   *jet_AntiKt4LCTopo_BCH_CORR_JET_FORCELL;
   vector<float>   *jet_AntiKt4LCTopo_ENG_BAD_CELLS;
   vector<float>   *jet_AntiKt4LCTopo_N_BAD_CELLS;
   vector<float>   *jet_AntiKt4LCTopo_N_BAD_CELLS_CORR;
   vector<float>   *jet_AntiKt4LCTopo_BAD_CELLS_CORR_E;
   vector<float>   *jet_AntiKt4LCTopo_NumTowers;
   vector<float>   *jet_AntiKt4LCTopo_ootFracCells5;
   vector<float>   *jet_AntiKt4LCTopo_ootFracCells10;
   vector<float>   *jet_AntiKt4LCTopo_ootFracClusters5;
   vector<float>   *jet_AntiKt4LCTopo_ootFracClusters10;
   vector<int>     *jet_AntiKt4LCTopo_SamplingMax;
   vector<float>   *jet_AntiKt4LCTopo_fracSamplingMax;
   vector<float>   *jet_AntiKt4LCTopo_hecf;
   vector<float>   *jet_AntiKt4LCTopo_tgap3f;
   vector<int>     *jet_AntiKt4LCTopo_isUgly;
   vector<int>     *jet_AntiKt4LCTopo_isBadLooseMinus;
   vector<int>     *jet_AntiKt4LCTopo_isBadLoose;
   vector<int>     *jet_AntiKt4LCTopo_isBadMedium;
   vector<int>     *jet_AntiKt4LCTopo_isBadTight;
   vector<float>   *jet_AntiKt4LCTopo_emfrac;
   vector<float>   *jet_AntiKt4LCTopo_emscale_E;
   vector<float>   *jet_AntiKt4LCTopo_emscale_pt;
   vector<float>   *jet_AntiKt4LCTopo_emscale_m;
   vector<float>   *jet_AntiKt4LCTopo_emscale_eta;
   vector<float>   *jet_AntiKt4LCTopo_emscale_phi;
   vector<float>   *jet_AntiKt4LCTopo_ActiveArea;
   vector<float>   *jet_AntiKt4LCTopo_ActiveAreaPx;
   vector<float>   *jet_AntiKt4LCTopo_ActiveAreaPy;
   vector<float>   *jet_AntiKt4LCTopo_ActiveAreaPz;
   vector<float>   *jet_AntiKt4LCTopo_ActiveAreaE;
   vector<float>   *jet_AntiKt4LCTopo_jvtxf;
   vector<vector<float> > *jet_AntiKt4LCTopo_jvtxfFull;
   vector<float>   *jet_AntiKt4LCTopo_jvtx_x;
   vector<float>   *jet_AntiKt4LCTopo_jvtx_y;
   vector<float>   *jet_AntiKt4LCTopo_jvtx_z;
   vector<float>   *jet_AntiKt4LCTopo_TruthMFindex;
   vector<float>   *jet_AntiKt4LCTopo_TruthMF;
   vector<float>   *jet_AntiKt4LCTopo_e_PreSamplerB;
   vector<float>   *jet_AntiKt4LCTopo_e_EMB1;
   vector<float>   *jet_AntiKt4LCTopo_e_EMB2;
   vector<float>   *jet_AntiKt4LCTopo_e_EMB3;
   vector<float>   *jet_AntiKt4LCTopo_e_PreSamplerE;
   vector<float>   *jet_AntiKt4LCTopo_e_EME1;
   vector<float>   *jet_AntiKt4LCTopo_e_EME2;
   vector<float>   *jet_AntiKt4LCTopo_e_EME3;
   vector<float>   *jet_AntiKt4LCTopo_e_HEC0;
   vector<float>   *jet_AntiKt4LCTopo_e_HEC1;
   vector<float>   *jet_AntiKt4LCTopo_e_HEC2;
   vector<float>   *jet_AntiKt4LCTopo_e_HEC3;
   vector<float>   *jet_AntiKt4LCTopo_e_TileBar0;
   vector<float>   *jet_AntiKt4LCTopo_e_TileBar1;
   vector<float>   *jet_AntiKt4LCTopo_e_TileBar2;
   vector<float>   *jet_AntiKt4LCTopo_e_TileGap1;
   vector<float>   *jet_AntiKt4LCTopo_e_TileGap2;
   vector<float>   *jet_AntiKt4LCTopo_e_TileGap3;
   vector<float>   *jet_AntiKt4LCTopo_e_TileExt0;
   vector<float>   *jet_AntiKt4LCTopo_e_TileExt1;
   vector<float>   *jet_AntiKt4LCTopo_e_TileExt2;
   vector<float>   *jet_AntiKt4LCTopo_e_FCAL0;
   vector<float>   *jet_AntiKt4LCTopo_e_FCAL1;
   vector<float>   *jet_AntiKt4LCTopo_e_FCAL2;
   vector<float>   *jet_AntiKt4LCTopo_constscale_E;
   vector<float>   *jet_AntiKt4LCTopo_constscale_pt;
   vector<float>   *jet_AntiKt4LCTopo_constscale_m;
   vector<float>   *jet_AntiKt4LCTopo_constscale_eta;
   vector<float>   *jet_AntiKt4LCTopo_constscale_phi;
   vector<float>   *jet_AntiKt4LCTopo_CENTER_LAMBDA;
   vector<float>   *jet_AntiKt4LCTopo_SECOND_R;
   vector<float>   *jet_AntiKt4LCTopo_SECOND_LAMBDA;
   vector<float>   *jet_AntiKt4LCTopo_CLUSTER_PT;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_Comb;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_GbbNN;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_IP2D;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_IP3D;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_JetFitterCharm;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_MV1;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_SV0;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_SV1;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_SV2;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_SecondSoftMuonTagChi2;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_SoftMuonTagChi2;
   vector<int>     *jet_AntiKt4LCTopo_flavor_truth_label;
   vector<float>   *jet_AntiKt4LCTopo_flavor_truth_dRminToB;
   vector<float>   *jet_AntiKt4LCTopo_flavor_truth_dRminToC;
   vector<float>   *jet_AntiKt4LCTopo_flavor_truth_dRminToT;
   vector<int>     *jet_AntiKt4LCTopo_flavor_truth_BHadronpdg;
   vector<float>   *jet_AntiKt4LCTopo_flavor_truth_vx_x;
   vector<float>   *jet_AntiKt4LCTopo_flavor_truth_vx_y;
   vector<float>   *jet_AntiKt4LCTopo_flavor_truth_vx_z;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_ip2d_pu;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_ip2d_pb;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_ip2d_isValid;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_ip2d_ntrk;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_ip3d_pu;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_ip3d_pb;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_ip3d_pc;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_ip3d_isValid;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv1_pu;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv1_pb;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv1_pc;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv2_pu;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv2_pb;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_pu;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_pb;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_pc;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitcomb_pu;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitcomb_pb;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitcomb_pc;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_jfit_nvtx;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_jfit_ntrkAtVx;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_efrc;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_mass;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_sig3d;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_deltaPhi;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_deltaEta;
   vector<vector<string> > *jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropName;
   vector<vector<double> > *jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropValue;
   vector<vector<string> > *jet_AntiKt4LCTopo_flavor_component_jfitc_intPropName;
   vector<vector<int> > *jet_AntiKt4LCTopo_flavor_component_jfitc_intPropValue;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitc_pu;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitc_pb;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitc_pc;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_ipplus_trk_n;
   vector<vector<int> > *jet_AntiKt4LCTopo_flavor_component_ipplus_trk_index;
   vector<vector<float> > *jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0val;
   vector<vector<float> > *jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0sig;
   vector<vector<float> > *jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0val;
   vector<vector<float> > *jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0sig;
   vector<vector<float> > *jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w2D;
   vector<vector<float> > *jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w3D;
   vector<vector<int> > *jet_AntiKt4LCTopo_flavor_component_ipplus_trk_grade;
   vector<vector<int> > *jet_AntiKt4LCTopo_flavor_component_ipplus_trk_isFromV0;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_svp_isValid;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_svp_ntrkv;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_svp_n2t;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_mass;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_efrc;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_x;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_y;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_z;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_err_x;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_err_y;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_err_z;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_cov_xy;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_cov_xz;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_cov_yz;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_chi2;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_svp_ndof;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_svp_ntrk;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_svp_trk_n;
   vector<vector<int> > *jet_AntiKt4LCTopo_flavor_component_svp_trk_index;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_sv0p_isValid;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_sv0p_n2t;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_mass;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_efrc;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_x;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_y;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_z;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_err_x;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_err_y;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_err_z;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xy;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xz;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_cov_yz;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_chi2;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_sv0p_ndof;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_sv0p_trk_n;
   vector<vector<int> > *jet_AntiKt4LCTopo_flavor_component_sv0p_trk_index;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_n;
   vector<vector<int> > *jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_index;
   vector<vector<float> > *jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_w;
   vector<vector<float> > *jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_pTRel;
   vector<vector<float> > *jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_dRJet;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_msvp_isValid;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_msvp_ntrkv;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_msvp_n2t;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_msvp_nVtx;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_msvp_normWeightedDist;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_n;
   vector<vector<int> > *jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_index;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitvx_theta;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitvx_phi;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitvx_errtheta;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitvx_errphi;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitvx_chi2;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfitvx_ndof;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_n;
   vector<vector<int> > *jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_index;
   vector<int>     *jet_AntiKt4LCTopo_bbDecision;
   vector<float>   *jet_AntiKt4LCTopo_TrackMFindex;
   vector<float>   *jet_AntiKt4LCTopo_TrackMF;
   vector<float>   *jet_AntiKt4LCTopo_nTrk_pv0_1GeV;
   vector<float>   *jet_AntiKt4LCTopo_sumPtTrk_pv0_1GeV;
   vector<float>   *jet_AntiKt4LCTopo_nTrk_allpv_1GeV;
   vector<float>   *jet_AntiKt4LCTopo_sumPtTrk_allpv_1GeV;
   vector<float>   *jet_AntiKt4LCTopo_nTrk_pv0_500MeV;
   vector<float>   *jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV;
   vector<float>   *jet_AntiKt4LCTopo_trackWIDTH_pv0_1GeV;
   vector<float>   *jet_AntiKt4LCTopo_trackWIDTH_allpv_1GeV;
   vector<float>   *jet_AntiKt4LCTopo_LikeLihood_0;
   vector<float>   *jet_AntiKt4LCTopo_LowEtConstituentsFrac;
   vector<float>   *jet_AntiKt4LCTopo_KtDr;
   vector<float>   *jet_AntiKt4LCTopo_SPLIT12;
   vector<float>   *jet_AntiKt4LCTopo_SPLIT23;
   vector<float>   *jet_AntiKt4LCTopo_SPLIT34;
   vector<vector<int> > *jet_AntiKt4LCTopo_truthAssoc_index;
   vector<int>     *jet_AntiKt4LCTopo_constit_n;
   vector<vector<int> > *jet_AntiKt4LCTopo_constit_index;
   vector<int>     *jet_AntiKt4LCTopo_TrackAssoc_n;
   vector<vector<int> > *jet_AntiKt4LCTopo_TrackAssoc_index;
   Int_t           trk_n;
   vector<float>   *trk_pt;
   vector<float>   *trk_eta;
   vector<float>   *trk_d0_wrtPV;
   vector<float>   *trk_z0_wrtPV;
   vector<float>   *trk_phi_wrtPV;
   vector<float>   *trk_theta_wrtPV;
   vector<float>   *trk_qoverp_wrtPV;
   vector<float>   *trk_cov_d0_wrtPV;
   vector<float>   *trk_cov_z0_wrtPV;
   vector<float>   *trk_cov_phi_wrtPV;
   vector<float>   *trk_cov_theta_wrtPV;
   vector<float>   *trk_cov_qoverp_wrtPV;
   vector<float>   *trk_d0_wrtBL;
   vector<float>   *trk_d0_err_wrtBL;
   vector<float>   *trk_eta_atCalo;
   vector<float>   *trk_phi_atCalo;
   vector<float>   *trk_eta_atCalo2ndLayer;
   vector<float>   *trk_phi_atCalo2ndLayer;
   vector<float>   *trk_chi2;
   vector<int>     *trk_ndof;
   vector<int>     *trk_nBLHits;
   vector<int>     *trk_nPixHits;
   vector<int>     *trk_nSCTHits;
   vector<int>     *trk_nTRTHits;
   vector<int>     *trk_nTRTHighTHits;
   vector<int>     *trk_nPixHoles;
   vector<int>     *trk_nSCTHoles;
   vector<int>     *trk_nTRTHoles;
   vector<int>     *trk_nPixelDeadSensors;
   vector<int>     *trk_nSCTDeadSensors;
   vector<int>     *trk_expectBLayerHit;
   vector<float>   *trk_mc_probability;
   vector<int>     *trk_mc_barcode;
   vector<int>     *trk_mc_index;
   Int_t           vxp_n;
   vector<float>   *vxp_x;
   vector<float>   *vxp_y;
   vector<float>   *vxp_z;
   vector<float>   *vxp_cov_x;
   vector<float>   *vxp_cov_y;
   vector<float>   *vxp_cov_z;
   vector<float>   *vxp_cov_xy;
   vector<float>   *vxp_cov_xz;
   vector<float>   *vxp_cov_yz;
   vector<int>     *vxp_type;
   vector<float>   *vxp_chi2;
   vector<int>     *vxp_ndof;
   vector<float>   *vxp_px;
   vector<float>   *vxp_py;
   vector<float>   *vxp_pz;
   vector<float>   *vxp_E;
   vector<float>   *vxp_m;
   vector<int>     *vxp_nTracks;
   vector<float>   *vxp_sumPt;
   vector<vector<float> > *vxp_trk_weight;
   vector<int>     *vxp_trk_n;
   vector<vector<int> > *vxp_trk_index;

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_jet_AntiKt4LCTopo_n;   //!
   TBranch        *b_jet_AntiKt4LCTopo_E;   //!
   TBranch        *b_jet_AntiKt4LCTopo_pt;   //!
   TBranch        *b_jet_AntiKt4LCTopo_m;   //!
   TBranch        *b_jet_AntiKt4LCTopo_eta;   //!
   TBranch        *b_jet_AntiKt4LCTopo_phi;   //!
   TBranch        *b_jet_AntiKt4LCTopo_EtaOrigin;   //!
   TBranch        *b_jet_AntiKt4LCTopo_PhiOrigin;   //!
   TBranch        *b_jet_AntiKt4LCTopo_MOrigin;   //!
   TBranch        *b_jet_AntiKt4LCTopo_WIDTH;   //!
   TBranch        *b_jet_AntiKt4LCTopo_Timing;   //!
   TBranch        *b_jet_AntiKt4LCTopo_LArQuality;   //!
   TBranch        *b_jet_AntiKt4LCTopo_OriginIndex;   //!
   TBranch        *b_jet_AntiKt4LCTopo_HECQuality;   //!
   TBranch        *b_jet_AntiKt4LCTopo_NegativeE;   //!
   TBranch        *b_jet_AntiKt4LCTopo_AverageLArQF;   //!
   TBranch        *b_jet_AntiKt4LCTopo_BCH_CORR_CELL;   //!
   TBranch        *b_jet_AntiKt4LCTopo_BCH_CORR_DOTX;   //!
   TBranch        *b_jet_AntiKt4LCTopo_BCH_CORR_JET;   //!
   TBranch        *b_jet_AntiKt4LCTopo_BCH_CORR_JET_FORCELL;   //!
   TBranch        *b_jet_AntiKt4LCTopo_ENG_BAD_CELLS;   //!
   TBranch        *b_jet_AntiKt4LCTopo_N_BAD_CELLS;   //!
   TBranch        *b_jet_AntiKt4LCTopo_N_BAD_CELLS_CORR;   //!
   TBranch        *b_jet_AntiKt4LCTopo_BAD_CELLS_CORR_E;   //!
   TBranch        *b_jet_AntiKt4LCTopo_NumTowers;   //!
   TBranch        *b_jet_AntiKt4LCTopo_ootFracCells5;   //!
   TBranch        *b_jet_AntiKt4LCTopo_ootFracCells10;   //!
   TBranch        *b_jet_AntiKt4LCTopo_ootFracClusters5;   //!
   TBranch        *b_jet_AntiKt4LCTopo_ootFracClusters10;   //!
   TBranch        *b_jet_AntiKt4LCTopo_SamplingMax;   //!
   TBranch        *b_jet_AntiKt4LCTopo_fracSamplingMax;   //!
   TBranch        *b_jet_AntiKt4LCTopo_hecf;   //!
   TBranch        *b_jet_AntiKt4LCTopo_tgap3f;   //!
   TBranch        *b_jet_AntiKt4LCTopo_isUgly;   //!
   TBranch        *b_jet_AntiKt4LCTopo_isBadLooseMinus;   //!
   TBranch        *b_jet_AntiKt4LCTopo_isBadLoose;   //!
   TBranch        *b_jet_AntiKt4LCTopo_isBadMedium;   //!
   TBranch        *b_jet_AntiKt4LCTopo_isBadTight;   //!
   TBranch        *b_jet_AntiKt4LCTopo_emfrac;   //!
   TBranch        *b_jet_AntiKt4LCTopo_emscale_E;   //!
   TBranch        *b_jet_AntiKt4LCTopo_emscale_pt;   //!
   TBranch        *b_jet_AntiKt4LCTopo_emscale_m;   //!
   TBranch        *b_jet_AntiKt4LCTopo_emscale_eta;   //!
   TBranch        *b_jet_AntiKt4LCTopo_emscale_phi;   //!
   TBranch        *b_jet_AntiKt4LCTopo_ActiveArea;   //!
   TBranch        *b_jet_AntiKt4LCTopo_ActiveAreaPx;   //!
   TBranch        *b_jet_AntiKt4LCTopo_ActiveAreaPy;   //!
   TBranch        *b_jet_AntiKt4LCTopo_ActiveAreaPz;   //!
   TBranch        *b_jet_AntiKt4LCTopo_ActiveAreaE;   //!
   TBranch        *b_jet_AntiKt4LCTopo_jvtxf;   //!
   TBranch        *b_jet_AntiKt4LCTopo_jvtxfFull;   //!
   TBranch        *b_jet_AntiKt4LCTopo_jvtx_x;   //!
   TBranch        *b_jet_AntiKt4LCTopo_jvtx_y;   //!
   TBranch        *b_jet_AntiKt4LCTopo_jvtx_z;   //!
   TBranch        *b_jet_AntiKt4LCTopo_TruthMFindex;   //!
   TBranch        *b_jet_AntiKt4LCTopo_TruthMF;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_PreSamplerB;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_EMB1;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_EMB2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_EMB3;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_PreSamplerE;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_EME1;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_EME2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_EME3;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_HEC0;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_HEC1;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_HEC2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_HEC3;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_TileBar0;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_TileBar1;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_TileBar2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_TileGap1;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_TileGap2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_TileGap3;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_TileExt0;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_TileExt1;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_TileExt2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_FCAL0;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_FCAL1;   //!
   TBranch        *b_jet_AntiKt4LCTopo_e_FCAL2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_constscale_E;   //!
   TBranch        *b_jet_AntiKt4LCTopo_constscale_pt;   //!
   TBranch        *b_jet_AntiKt4LCTopo_constscale_m;   //!
   TBranch        *b_jet_AntiKt4LCTopo_constscale_eta;   //!
   TBranch        *b_jet_AntiKt4LCTopo_constscale_phi;   //!
   TBranch        *b_jet_AntiKt4LCTopo_CENTER_LAMBDA;   //!
   TBranch        *b_jet_AntiKt4LCTopo_SECOND_R;   //!
   TBranch        *b_jet_AntiKt4LCTopo_SECOND_LAMBDA;   //!
   TBranch        *b_jet_AntiKt4LCTopo_CLUSTER_PT;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_Comb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_GbbNN;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_IP2D;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_IP3D;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_JetFitterCharm;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_MV1;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_SV0;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_SV1;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_SV2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_SecondSoftMuonTagChi2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_SoftMuonTagChi2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_truth_label;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_truth_dRminToB;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_truth_dRminToC;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_truth_dRminToT;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_truth_BHadronpdg;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_truth_vx_x;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_truth_vx_y;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_truth_vx_z;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip2d_pu;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip2d_pb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip2d_isValid;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip2d_ntrk;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip3d_pu;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip3d_pb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip3d_pc;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip3d_isValid;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv1_pu;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv1_pb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv1_pc;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv2_pu;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv2_pb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_pu;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_pb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_pc;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitcomb_pu;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitcomb_pb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitcomb_pc;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_nvtx;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_ntrkAtVx;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_efrc;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_mass;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_sig3d;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_deltaPhi;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_deltaEta;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropName;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropValue;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitc_intPropName;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitc_intPropValue;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitc_pu;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitc_pb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitc_pc;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_n;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_index;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0val;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0sig;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0val;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0sig;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w2D;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w3D;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_grade;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_isFromV0;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_isValid;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_ntrkv;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_n2t;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_mass;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_efrc;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_x;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_y;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_z;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_err_x;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_err_y;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_err_z;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_cov_xy;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_cov_xz;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_cov_yz;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_chi2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_ndof;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_ntrk;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_trk_n;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_trk_index;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_isValid;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_n2t;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_mass;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_efrc;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_x;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_y;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_z;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_err_x;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_err_y;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_err_z;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xy;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xz;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_cov_yz;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_chi2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_ndof;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_trk_n;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_trk_index;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_n;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_index;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_w;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_pTRel;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_dRJet;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_msvp_isValid;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_msvp_ntrkv;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_msvp_n2t;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_msvp_nVtx;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_msvp_normWeightedDist;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_n;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_index;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitvx_theta;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitvx_phi;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitvx_errtheta;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitvx_errphi;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitvx_chi2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfitvx_ndof;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_n;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_index;   //!
   TBranch        *b_jet_AntiKt4LCTopo_bbDecision;   //!
   TBranch        *b_jet_AntiKt4LCTopo_TrackMFindex;   //!
   TBranch        *b_jet_AntiKt4LCTopo_TrackMF;   //!
   TBranch        *b_jet_AntiKt4LCTopo_nTrk_pv0_1GeV;   //!
   TBranch        *b_jet_AntiKt4LCTopo_sumPtTrk_pv0_1GeV;   //!
   TBranch        *b_jet_AntiKt4LCTopo_nTrk_allpv_1GeV;   //!
   TBranch        *b_jet_AntiKt4LCTopo_sumPtTrk_allpv_1GeV;   //!
   TBranch        *b_jet_AntiKt4LCTopo_nTrk_pv0_500MeV;   //!
   TBranch        *b_jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV;   //!
   TBranch        *b_jet_AntiKt4LCTopo_trackWIDTH_pv0_1GeV;   //!
   TBranch        *b_jet_AntiKt4LCTopo_trackWIDTH_allpv_1GeV;   //!
   TBranch        *b_jet_AntiKt4LCTopo_LikeLihood_0;   //!
   TBranch        *b_jet_AntiKt4LCTopo_LowEtConstituentsFrac;   //!
   TBranch        *b_jet_AntiKt4LCTopo_KtDr;   //!
   TBranch        *b_jet_AntiKt4LCTopo_SPLIT12;   //!
   TBranch        *b_jet_AntiKt4LCTopo_SPLIT23;   //!
   TBranch        *b_jet_AntiKt4LCTopo_SPLIT34;   //!
   TBranch        *b_jet_AntiKt4LCTopo_truthAssoc_index;   //!
   TBranch        *b_jet_AntiKt4LCTopo_constit_n;   //!
   TBranch        *b_jet_AntiKt4LCTopo_constit_index;   //!
   TBranch        *b_jet_AntiKt4LCTopo_TrackAssoc_n;   //!
   TBranch        *b_jet_AntiKt4LCTopo_TrackAssoc_index;   //!
   TBranch        *b_trk_n;   //!
   TBranch        *b_trk_pt;   //!
   TBranch        *b_trk_eta;   //!
   TBranch        *b_trk_d0_wrtPV;   //!
   TBranch        *b_trk_z0_wrtPV;   //!
   TBranch        *b_trk_phi_wrtPV;   //!
   TBranch        *b_trk_theta_wrtPV;   //!
   TBranch        *b_trk_qoverp_wrtPV;   //!
   TBranch        *b_trk_cov_d0_wrtPV;   //!
   TBranch        *b_trk_cov_z0_wrtPV;   //!
   TBranch        *b_trk_cov_phi_wrtPV;   //!
   TBranch        *b_trk_cov_theta_wrtPV;   //!
   TBranch        *b_trk_cov_qoverp_wrtPV;   //!
   TBranch        *b_trk_d0_wrtBL;   //!
   TBranch        *b_trk_d0_err_wrtBL;   //!
   TBranch        *b_trk_eta_atCalo;   //!
   TBranch        *b_trk_phi_atCalo;   //!
   TBranch        *b_trk_eta_atCalo2ndLayer;   //!
   TBranch        *b_trk_phi_atCalo2ndLayer;   //!
   TBranch        *b_trk_chi2;   //!
   TBranch        *b_trk_ndof;   //!
   TBranch        *b_trk_nBLHits;   //!
   TBranch        *b_trk_nPixHits;   //!
   TBranch        *b_trk_nSCTHits;   //!
   TBranch        *b_trk_nTRTHits;   //!
   TBranch        *b_trk_nTRTHighTHits;   //!
   TBranch        *b_trk_nPixHoles;   //!
   TBranch        *b_trk_nSCTHoles;   //!
   TBranch        *b_trk_nTRTHoles;   //!
   TBranch        *b_trk_nPixelDeadSensors;   //!
   TBranch        *b_trk_nSCTDeadSensors;   //!
   TBranch        *b_trk_expectBLayerHit;   //!
   TBranch        *b_trk_mc_probability;   //!
   TBranch        *b_trk_mc_barcode;   //!
   TBranch        *b_trk_mc_index;   //!
   TBranch        *b_vxp_n;   //!
   TBranch        *b_vxp_x;   //!
   TBranch        *b_vxp_y;   //!
   TBranch        *b_vxp_z;   //!
   TBranch        *b_vxp_cov_x;   //!
   TBranch        *b_vxp_cov_y;   //!
   TBranch        *b_vxp_cov_z;   //!
   TBranch        *b_vxp_cov_xy;   //!
   TBranch        *b_vxp_cov_xz;   //!
   TBranch        *b_vxp_cov_yz;   //!
   TBranch        *b_vxp_type;   //!
   TBranch        *b_vxp_chi2;   //!
   TBranch        *b_vxp_ndof;   //!
   TBranch        *b_vxp_px;   //!
   TBranch        *b_vxp_py;   //!
   TBranch        *b_vxp_pz;   //!
   TBranch        *b_vxp_E;   //!
   TBranch        *b_vxp_m;   //!
   TBranch        *b_vxp_nTracks;   //!
   TBranch        *b_vxp_sumPt;   //!
   TBranch        *b_vxp_trk_weight;   //!
   TBranch        *b_vxp_trk_n;   //!
   TBranch        *b_vxp_trk_index;   //!

   NtupCommonReader(TTree *tree=0);
   virtual ~NtupCommonReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NtupCommonReader_cxx
NtupCommonReader::NtupCommonReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/atlas/dq2//mc12_8TeV/NTUP_COMMON/e1434_s1499_s1504_r3658_r3549_p1562/mc12_8TeV.147771.Sherpa_CT10_Zmumu.merge.NTUP_COMMON.e1434_s1499_s1504_r3658_r3549_p1562_tid01313919_00/NTUP_COMMON.01313919._000001.root.1");
      if (!f || !f->IsOpen()) {
         f = new TFile("/atlas/dq2//mc12_8TeV/NTUP_COMMON/e1434_s1499_s1504_r3658_r3549_p1562/mc12_8TeV.147771.Sherpa_CT10_Zmumu.merge.NTUP_COMMON.e1434_s1499_s1504_r3658_r3549_p1562_tid01313919_00/NTUP_COMMON.01313919._000001.root.1");
      }
      f->GetObject("physics",tree);

   }
   Init(tree);
}

NtupCommonReader::~NtupCommonReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NtupCommonReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NtupCommonReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NtupCommonReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   EventNumber = 0;
   jet_AntiKt4LCTopo_E = 0;
   jet_AntiKt4LCTopo_pt = 0;
   jet_AntiKt4LCTopo_m = 0;
   jet_AntiKt4LCTopo_eta = 0;
   jet_AntiKt4LCTopo_phi = 0;
   jet_AntiKt4LCTopo_EtaOrigin = 0;
   jet_AntiKt4LCTopo_PhiOrigin = 0;
   jet_AntiKt4LCTopo_MOrigin = 0;
   jet_AntiKt4LCTopo_WIDTH = 0;
   jet_AntiKt4LCTopo_Timing = 0;
   jet_AntiKt4LCTopo_LArQuality = 0;
   jet_AntiKt4LCTopo_OriginIndex = 0;
   jet_AntiKt4LCTopo_HECQuality = 0;
   jet_AntiKt4LCTopo_NegativeE = 0;
   jet_AntiKt4LCTopo_AverageLArQF = 0;
   jet_AntiKt4LCTopo_BCH_CORR_CELL = 0;
   jet_AntiKt4LCTopo_BCH_CORR_DOTX = 0;
   jet_AntiKt4LCTopo_BCH_CORR_JET = 0;
   jet_AntiKt4LCTopo_BCH_CORR_JET_FORCELL = 0;
   jet_AntiKt4LCTopo_ENG_BAD_CELLS = 0;
   jet_AntiKt4LCTopo_N_BAD_CELLS = 0;
   jet_AntiKt4LCTopo_N_BAD_CELLS_CORR = 0;
   jet_AntiKt4LCTopo_BAD_CELLS_CORR_E = 0;
   jet_AntiKt4LCTopo_NumTowers = 0;
   jet_AntiKt4LCTopo_ootFracCells5 = 0;
   jet_AntiKt4LCTopo_ootFracCells10 = 0;
   jet_AntiKt4LCTopo_ootFracClusters5 = 0;
   jet_AntiKt4LCTopo_ootFracClusters10 = 0;
   jet_AntiKt4LCTopo_SamplingMax = 0;
   jet_AntiKt4LCTopo_fracSamplingMax = 0;
   jet_AntiKt4LCTopo_hecf = 0;
   jet_AntiKt4LCTopo_tgap3f = 0;
   jet_AntiKt4LCTopo_isUgly = 0;
   jet_AntiKt4LCTopo_isBadLooseMinus = 0;
   jet_AntiKt4LCTopo_isBadLoose = 0;
   jet_AntiKt4LCTopo_isBadMedium = 0;
   jet_AntiKt4LCTopo_isBadTight = 0;
   jet_AntiKt4LCTopo_emfrac = 0;
   jet_AntiKt4LCTopo_emscale_E = 0;
   jet_AntiKt4LCTopo_emscale_pt = 0;
   jet_AntiKt4LCTopo_emscale_m = 0;
   jet_AntiKt4LCTopo_emscale_eta = 0;
   jet_AntiKt4LCTopo_emscale_phi = 0;
   jet_AntiKt4LCTopo_ActiveArea = 0;
   jet_AntiKt4LCTopo_ActiveAreaPx = 0;
   jet_AntiKt4LCTopo_ActiveAreaPy = 0;
   jet_AntiKt4LCTopo_ActiveAreaPz = 0;
   jet_AntiKt4LCTopo_ActiveAreaE = 0;
   jet_AntiKt4LCTopo_jvtxf = 0;
   jet_AntiKt4LCTopo_jvtxfFull = 0;
   jet_AntiKt4LCTopo_jvtx_x = 0;
   jet_AntiKt4LCTopo_jvtx_y = 0;
   jet_AntiKt4LCTopo_jvtx_z = 0;
   jet_AntiKt4LCTopo_TruthMFindex = 0;
   jet_AntiKt4LCTopo_TruthMF = 0;
   jet_AntiKt4LCTopo_e_PreSamplerB = 0;
   jet_AntiKt4LCTopo_e_EMB1 = 0;
   jet_AntiKt4LCTopo_e_EMB2 = 0;
   jet_AntiKt4LCTopo_e_EMB3 = 0;
   jet_AntiKt4LCTopo_e_PreSamplerE = 0;
   jet_AntiKt4LCTopo_e_EME1 = 0;
   jet_AntiKt4LCTopo_e_EME2 = 0;
   jet_AntiKt4LCTopo_e_EME3 = 0;
   jet_AntiKt4LCTopo_e_HEC0 = 0;
   jet_AntiKt4LCTopo_e_HEC1 = 0;
   jet_AntiKt4LCTopo_e_HEC2 = 0;
   jet_AntiKt4LCTopo_e_HEC3 = 0;
   jet_AntiKt4LCTopo_e_TileBar0 = 0;
   jet_AntiKt4LCTopo_e_TileBar1 = 0;
   jet_AntiKt4LCTopo_e_TileBar2 = 0;
   jet_AntiKt4LCTopo_e_TileGap1 = 0;
   jet_AntiKt4LCTopo_e_TileGap2 = 0;
   jet_AntiKt4LCTopo_e_TileGap3 = 0;
   jet_AntiKt4LCTopo_e_TileExt0 = 0;
   jet_AntiKt4LCTopo_e_TileExt1 = 0;
   jet_AntiKt4LCTopo_e_TileExt2 = 0;
   jet_AntiKt4LCTopo_e_FCAL0 = 0;
   jet_AntiKt4LCTopo_e_FCAL1 = 0;
   jet_AntiKt4LCTopo_e_FCAL2 = 0;
   jet_AntiKt4LCTopo_constscale_E = 0;
   jet_AntiKt4LCTopo_constscale_pt = 0;
   jet_AntiKt4LCTopo_constscale_m = 0;
   jet_AntiKt4LCTopo_constscale_eta = 0;
   jet_AntiKt4LCTopo_constscale_phi = 0;
   jet_AntiKt4LCTopo_CENTER_LAMBDA = 0;
   jet_AntiKt4LCTopo_SECOND_R = 0;
   jet_AntiKt4LCTopo_SECOND_LAMBDA = 0;
   jet_AntiKt4LCTopo_CLUSTER_PT = 0;
   jet_AntiKt4LCTopo_flavor_weight_Comb = 0;
   jet_AntiKt4LCTopo_flavor_weight_GbbNN = 0;
   jet_AntiKt4LCTopo_flavor_weight_IP2D = 0;
   jet_AntiKt4LCTopo_flavor_weight_IP3D = 0;
   jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN = 0;
   jet_AntiKt4LCTopo_flavor_weight_JetFitterCharm = 0;
   jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN = 0;
   jet_AntiKt4LCTopo_flavor_weight_MV1 = 0;
   jet_AntiKt4LCTopo_flavor_weight_SV0 = 0;
   jet_AntiKt4LCTopo_flavor_weight_SV1 = 0;
   jet_AntiKt4LCTopo_flavor_weight_SV2 = 0;
   jet_AntiKt4LCTopo_flavor_weight_SecondSoftMuonTagChi2 = 0;
   jet_AntiKt4LCTopo_flavor_weight_SoftMuonTagChi2 = 0;
   jet_AntiKt4LCTopo_flavor_truth_label = 0;
   jet_AntiKt4LCTopo_flavor_truth_dRminToB = 0;
   jet_AntiKt4LCTopo_flavor_truth_dRminToC = 0;
   jet_AntiKt4LCTopo_flavor_truth_dRminToT = 0;
   jet_AntiKt4LCTopo_flavor_truth_BHadronpdg = 0;
   jet_AntiKt4LCTopo_flavor_truth_vx_x = 0;
   jet_AntiKt4LCTopo_flavor_truth_vx_y = 0;
   jet_AntiKt4LCTopo_flavor_truth_vx_z = 0;
   jet_AntiKt4LCTopo_flavor_component_ip2d_pu = 0;
   jet_AntiKt4LCTopo_flavor_component_ip2d_pb = 0;
   jet_AntiKt4LCTopo_flavor_component_ip2d_isValid = 0;
   jet_AntiKt4LCTopo_flavor_component_ip2d_ntrk = 0;
   jet_AntiKt4LCTopo_flavor_component_ip3d_pu = 0;
   jet_AntiKt4LCTopo_flavor_component_ip3d_pb = 0;
   jet_AntiKt4LCTopo_flavor_component_ip3d_pc = 0;
   jet_AntiKt4LCTopo_flavor_component_ip3d_isValid = 0;
   jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk = 0;
   jet_AntiKt4LCTopo_flavor_component_sv1_pu = 0;
   jet_AntiKt4LCTopo_flavor_component_sv1_pb = 0;
   jet_AntiKt4LCTopo_flavor_component_sv1_pc = 0;
   jet_AntiKt4LCTopo_flavor_component_sv2_pu = 0;
   jet_AntiKt4LCTopo_flavor_component_sv2_pb = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_pu = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_pb = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_pc = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitcomb_pu = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitcomb_pb = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitcomb_pc = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_nvtx = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_ntrkAtVx = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_efrc = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_mass = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_sig3d = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_deltaPhi = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_deltaEta = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropName = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropValue = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitc_intPropName = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitc_intPropValue = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitc_pu = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitc_pb = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitc_pc = 0;
   jet_AntiKt4LCTopo_flavor_component_ipplus_trk_n = 0;
   jet_AntiKt4LCTopo_flavor_component_ipplus_trk_index = 0;
   jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0val = 0;
   jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0sig = 0;
   jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0val = 0;
   jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0sig = 0;
   jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w2D = 0;
   jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w3D = 0;
   jet_AntiKt4LCTopo_flavor_component_ipplus_trk_grade = 0;
   jet_AntiKt4LCTopo_flavor_component_ipplus_trk_isFromV0 = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_isValid = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_ntrkv = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_n2t = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_mass = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_efrc = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_x = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_y = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_z = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_err_x = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_err_y = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_err_z = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_cov_xy = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_cov_xz = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_cov_yz = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_chi2 = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_ndof = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_ntrk = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_trk_n = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_trk_index = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_isValid = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_n2t = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_mass = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_efrc = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_x = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_y = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_z = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_err_x = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_err_y = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_err_z = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xy = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xz = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_cov_yz = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_chi2 = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_ndof = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_trk_n = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_trk_index = 0;
   jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_n = 0;
   jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_index = 0;
   jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_w = 0;
   jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_pTRel = 0;
   jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_dRJet = 0;
   jet_AntiKt4LCTopo_flavor_component_msvp_isValid = 0;
   jet_AntiKt4LCTopo_flavor_component_msvp_ntrkv = 0;
   jet_AntiKt4LCTopo_flavor_component_msvp_n2t = 0;
   jet_AntiKt4LCTopo_flavor_component_msvp_nVtx = 0;
   jet_AntiKt4LCTopo_flavor_component_msvp_normWeightedDist = 0;
   jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_n = 0;
   jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_index = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitvx_theta = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitvx_phi = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitvx_errtheta = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitvx_errphi = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitvx_chi2 = 0;
   jet_AntiKt4LCTopo_flavor_component_jfitvx_ndof = 0;
   jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_n = 0;
   jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_index = 0;
   jet_AntiKt4LCTopo_bbDecision = 0;
   jet_AntiKt4LCTopo_TrackMFindex = 0;
   jet_AntiKt4LCTopo_TrackMF = 0;
   jet_AntiKt4LCTopo_nTrk_pv0_1GeV = 0;
   jet_AntiKt4LCTopo_sumPtTrk_pv0_1GeV = 0;
   jet_AntiKt4LCTopo_nTrk_allpv_1GeV = 0;
   jet_AntiKt4LCTopo_sumPtTrk_allpv_1GeV = 0;
   jet_AntiKt4LCTopo_nTrk_pv0_500MeV = 0;
   jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV = 0;
   jet_AntiKt4LCTopo_trackWIDTH_pv0_1GeV = 0;
   jet_AntiKt4LCTopo_trackWIDTH_allpv_1GeV = 0;
   jet_AntiKt4LCTopo_LikeLihood_0 = 0;
   jet_AntiKt4LCTopo_LowEtConstituentsFrac = 0;
   jet_AntiKt4LCTopo_KtDr = 0;
   jet_AntiKt4LCTopo_SPLIT12 = 0;
   jet_AntiKt4LCTopo_SPLIT23 = 0;
   jet_AntiKt4LCTopo_SPLIT34 = 0;
   jet_AntiKt4LCTopo_truthAssoc_index = 0;
   jet_AntiKt4LCTopo_constit_n = 0;
   jet_AntiKt4LCTopo_constit_index = 0;
   jet_AntiKt4LCTopo_TrackAssoc_n = 0;
   jet_AntiKt4LCTopo_TrackAssoc_index = 0;
   trk_pt = 0;
   trk_eta = 0;
   trk_d0_wrtPV = 0;
   trk_z0_wrtPV = 0;
   trk_phi_wrtPV = 0;
   trk_theta_wrtPV = 0;
   trk_qoverp_wrtPV = 0;
   trk_cov_d0_wrtPV = 0;
   trk_cov_z0_wrtPV = 0;
   trk_cov_phi_wrtPV = 0;
   trk_cov_theta_wrtPV = 0;
   trk_cov_qoverp_wrtPV = 0;
   trk_d0_wrtBL = 0;
   trk_d0_err_wrtBL = 0;
   trk_eta_atCalo = 0;
   trk_phi_atCalo = 0;
   trk_eta_atCalo2ndLayer = 0;
   trk_phi_atCalo2ndLayer = 0;
   trk_chi2 = 0;
   trk_ndof = 0;
   trk_nBLHits = 0;
   trk_nPixHits = 0;
   trk_nSCTHits = 0;
   trk_nTRTHits = 0;
   trk_nTRTHighTHits = 0;
   trk_nPixHoles = 0;
   trk_nSCTHoles = 0;
   trk_nTRTHoles = 0;
   trk_nPixelDeadSensors = 0;
   trk_nSCTDeadSensors = 0;
   trk_expectBLayerHit = 0;
   trk_mc_probability = 0;
   trk_mc_barcode = 0;
   trk_mc_index = 0;
   vxp_x = 0;
   vxp_y = 0;
   vxp_z = 0;
   vxp_cov_x = 0;
   vxp_cov_y = 0;
   vxp_cov_z = 0;
   vxp_cov_xy = 0;
   vxp_cov_xz = 0;
   vxp_cov_yz = 0;
   vxp_type = 0;
   vxp_chi2 = 0;
   vxp_ndof = 0;
   vxp_px = 0;
   vxp_py = 0;
   vxp_pz = 0;
   vxp_E = 0;
   vxp_m = 0;
   vxp_nTracks = 0;
   vxp_sumPt = 0;
   vxp_trk_weight = 0;
   vxp_trk_n = 0;
   vxp_trk_index = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

//   fChain->SetBranchStatus("*", 0);
//   fChain->SetBranchStatus("jet_AntiKt4LCTopo*", 1);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_n", &jet_AntiKt4LCTopo_n, &b_jet_AntiKt4LCTopo_n);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_E", &jet_AntiKt4LCTopo_E, &b_jet_AntiKt4LCTopo_E);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_pt", &jet_AntiKt4LCTopo_pt, &b_jet_AntiKt4LCTopo_pt);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_m", &jet_AntiKt4LCTopo_m, &b_jet_AntiKt4LCTopo_m);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_eta", &jet_AntiKt4LCTopo_eta, &b_jet_AntiKt4LCTopo_eta);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_phi", &jet_AntiKt4LCTopo_phi, &b_jet_AntiKt4LCTopo_phi);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_EtaOrigin", &jet_AntiKt4LCTopo_EtaOrigin, &b_jet_AntiKt4LCTopo_EtaOrigin);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_PhiOrigin", &jet_AntiKt4LCTopo_PhiOrigin, &b_jet_AntiKt4LCTopo_PhiOrigin);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_MOrigin", &jet_AntiKt4LCTopo_MOrigin, &b_jet_AntiKt4LCTopo_MOrigin);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_WIDTH", &jet_AntiKt4LCTopo_WIDTH, &b_jet_AntiKt4LCTopo_WIDTH);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_Timing", &jet_AntiKt4LCTopo_Timing, &b_jet_AntiKt4LCTopo_Timing);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_LArQuality", &jet_AntiKt4LCTopo_LArQuality, &b_jet_AntiKt4LCTopo_LArQuality);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_OriginIndex", &jet_AntiKt4LCTopo_OriginIndex, &b_jet_AntiKt4LCTopo_OriginIndex);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_HECQuality", &jet_AntiKt4LCTopo_HECQuality, &b_jet_AntiKt4LCTopo_HECQuality);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_NegativeE", &jet_AntiKt4LCTopo_NegativeE, &b_jet_AntiKt4LCTopo_NegativeE);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_AverageLArQF", &jet_AntiKt4LCTopo_AverageLArQF, &b_jet_AntiKt4LCTopo_AverageLArQF);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_BCH_CORR_CELL", &jet_AntiKt4LCTopo_BCH_CORR_CELL, &b_jet_AntiKt4LCTopo_BCH_CORR_CELL);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_BCH_CORR_DOTX", &jet_AntiKt4LCTopo_BCH_CORR_DOTX, &b_jet_AntiKt4LCTopo_BCH_CORR_DOTX);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_BCH_CORR_JET", &jet_AntiKt4LCTopo_BCH_CORR_JET, &b_jet_AntiKt4LCTopo_BCH_CORR_JET);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_BCH_CORR_JET_FORCELL", &jet_AntiKt4LCTopo_BCH_CORR_JET_FORCELL, &b_jet_AntiKt4LCTopo_BCH_CORR_JET_FORCELL);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ENG_BAD_CELLS", &jet_AntiKt4LCTopo_ENG_BAD_CELLS, &b_jet_AntiKt4LCTopo_ENG_BAD_CELLS);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_N_BAD_CELLS", &jet_AntiKt4LCTopo_N_BAD_CELLS, &b_jet_AntiKt4LCTopo_N_BAD_CELLS);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_N_BAD_CELLS_CORR", &jet_AntiKt4LCTopo_N_BAD_CELLS_CORR, &b_jet_AntiKt4LCTopo_N_BAD_CELLS_CORR);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_BAD_CELLS_CORR_E", &jet_AntiKt4LCTopo_BAD_CELLS_CORR_E, &b_jet_AntiKt4LCTopo_BAD_CELLS_CORR_E);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_NumTowers", &jet_AntiKt4LCTopo_NumTowers, &b_jet_AntiKt4LCTopo_NumTowers);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ootFracCells5", &jet_AntiKt4LCTopo_ootFracCells5, &b_jet_AntiKt4LCTopo_ootFracCells5);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ootFracCells10", &jet_AntiKt4LCTopo_ootFracCells10, &b_jet_AntiKt4LCTopo_ootFracCells10);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ootFracClusters5", &jet_AntiKt4LCTopo_ootFracClusters5, &b_jet_AntiKt4LCTopo_ootFracClusters5);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ootFracClusters10", &jet_AntiKt4LCTopo_ootFracClusters10, &b_jet_AntiKt4LCTopo_ootFracClusters10);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_SamplingMax", &jet_AntiKt4LCTopo_SamplingMax, &b_jet_AntiKt4LCTopo_SamplingMax);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_fracSamplingMax", &jet_AntiKt4LCTopo_fracSamplingMax, &b_jet_AntiKt4LCTopo_fracSamplingMax);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_hecf", &jet_AntiKt4LCTopo_hecf, &b_jet_AntiKt4LCTopo_hecf);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_tgap3f", &jet_AntiKt4LCTopo_tgap3f, &b_jet_AntiKt4LCTopo_tgap3f);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_isUgly", &jet_AntiKt4LCTopo_isUgly, &b_jet_AntiKt4LCTopo_isUgly);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_isBadLooseMinus", &jet_AntiKt4LCTopo_isBadLooseMinus, &b_jet_AntiKt4LCTopo_isBadLooseMinus);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_isBadLoose", &jet_AntiKt4LCTopo_isBadLoose, &b_jet_AntiKt4LCTopo_isBadLoose);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_isBadMedium", &jet_AntiKt4LCTopo_isBadMedium, &b_jet_AntiKt4LCTopo_isBadMedium);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_isBadTight", &jet_AntiKt4LCTopo_isBadTight, &b_jet_AntiKt4LCTopo_isBadTight);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emfrac", &jet_AntiKt4LCTopo_emfrac, &b_jet_AntiKt4LCTopo_emfrac);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_E", &jet_AntiKt4LCTopo_emscale_E, &b_jet_AntiKt4LCTopo_emscale_E);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_pt", &jet_AntiKt4LCTopo_emscale_pt, &b_jet_AntiKt4LCTopo_emscale_pt);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_m", &jet_AntiKt4LCTopo_emscale_m, &b_jet_AntiKt4LCTopo_emscale_m);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_eta", &jet_AntiKt4LCTopo_emscale_eta, &b_jet_AntiKt4LCTopo_emscale_eta);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_phi", &jet_AntiKt4LCTopo_emscale_phi, &b_jet_AntiKt4LCTopo_emscale_phi);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ActiveArea", &jet_AntiKt4LCTopo_ActiveArea, &b_jet_AntiKt4LCTopo_ActiveArea);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ActiveAreaPx", &jet_AntiKt4LCTopo_ActiveAreaPx, &b_jet_AntiKt4LCTopo_ActiveAreaPx);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ActiveAreaPy", &jet_AntiKt4LCTopo_ActiveAreaPy, &b_jet_AntiKt4LCTopo_ActiveAreaPy);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ActiveAreaPz", &jet_AntiKt4LCTopo_ActiveAreaPz, &b_jet_AntiKt4LCTopo_ActiveAreaPz);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ActiveAreaE", &jet_AntiKt4LCTopo_ActiveAreaE, &b_jet_AntiKt4LCTopo_ActiveAreaE);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_jvtxf", &jet_AntiKt4LCTopo_jvtxf, &b_jet_AntiKt4LCTopo_jvtxf);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_jvtxfFull", &jet_AntiKt4LCTopo_jvtxfFull, &b_jet_AntiKt4LCTopo_jvtxfFull);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_jvtx_x", &jet_AntiKt4LCTopo_jvtx_x, &b_jet_AntiKt4LCTopo_jvtx_x);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_jvtx_y", &jet_AntiKt4LCTopo_jvtx_y, &b_jet_AntiKt4LCTopo_jvtx_y);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_jvtx_z", &jet_AntiKt4LCTopo_jvtx_z, &b_jet_AntiKt4LCTopo_jvtx_z);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_TruthMFindex", &jet_AntiKt4LCTopo_TruthMFindex, &b_jet_AntiKt4LCTopo_TruthMFindex);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_TruthMF", &jet_AntiKt4LCTopo_TruthMF, &b_jet_AntiKt4LCTopo_TruthMF);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_PreSamplerB", &jet_AntiKt4LCTopo_e_PreSamplerB, &b_jet_AntiKt4LCTopo_e_PreSamplerB);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_EMB1", &jet_AntiKt4LCTopo_e_EMB1, &b_jet_AntiKt4LCTopo_e_EMB1);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_EMB2", &jet_AntiKt4LCTopo_e_EMB2, &b_jet_AntiKt4LCTopo_e_EMB2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_EMB3", &jet_AntiKt4LCTopo_e_EMB3, &b_jet_AntiKt4LCTopo_e_EMB3);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_PreSamplerE", &jet_AntiKt4LCTopo_e_PreSamplerE, &b_jet_AntiKt4LCTopo_e_PreSamplerE);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_EME1", &jet_AntiKt4LCTopo_e_EME1, &b_jet_AntiKt4LCTopo_e_EME1);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_EME2", &jet_AntiKt4LCTopo_e_EME2, &b_jet_AntiKt4LCTopo_e_EME2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_EME3", &jet_AntiKt4LCTopo_e_EME3, &b_jet_AntiKt4LCTopo_e_EME3);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_HEC0", &jet_AntiKt4LCTopo_e_HEC0, &b_jet_AntiKt4LCTopo_e_HEC0);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_HEC1", &jet_AntiKt4LCTopo_e_HEC1, &b_jet_AntiKt4LCTopo_e_HEC1);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_HEC2", &jet_AntiKt4LCTopo_e_HEC2, &b_jet_AntiKt4LCTopo_e_HEC2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_HEC3", &jet_AntiKt4LCTopo_e_HEC3, &b_jet_AntiKt4LCTopo_e_HEC3);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_TileBar0", &jet_AntiKt4LCTopo_e_TileBar0, &b_jet_AntiKt4LCTopo_e_TileBar0);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_TileBar1", &jet_AntiKt4LCTopo_e_TileBar1, &b_jet_AntiKt4LCTopo_e_TileBar1);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_TileBar2", &jet_AntiKt4LCTopo_e_TileBar2, &b_jet_AntiKt4LCTopo_e_TileBar2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_TileGap1", &jet_AntiKt4LCTopo_e_TileGap1, &b_jet_AntiKt4LCTopo_e_TileGap1);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_TileGap2", &jet_AntiKt4LCTopo_e_TileGap2, &b_jet_AntiKt4LCTopo_e_TileGap2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_TileGap3", &jet_AntiKt4LCTopo_e_TileGap3, &b_jet_AntiKt4LCTopo_e_TileGap3);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_TileExt0", &jet_AntiKt4LCTopo_e_TileExt0, &b_jet_AntiKt4LCTopo_e_TileExt0);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_TileExt1", &jet_AntiKt4LCTopo_e_TileExt1, &b_jet_AntiKt4LCTopo_e_TileExt1);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_TileExt2", &jet_AntiKt4LCTopo_e_TileExt2, &b_jet_AntiKt4LCTopo_e_TileExt2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_FCAL0", &jet_AntiKt4LCTopo_e_FCAL0, &b_jet_AntiKt4LCTopo_e_FCAL0);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_FCAL1", &jet_AntiKt4LCTopo_e_FCAL1, &b_jet_AntiKt4LCTopo_e_FCAL1);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_e_FCAL2", &jet_AntiKt4LCTopo_e_FCAL2, &b_jet_AntiKt4LCTopo_e_FCAL2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_constscale_E", &jet_AntiKt4LCTopo_constscale_E, &b_jet_AntiKt4LCTopo_constscale_E);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_constscale_pt", &jet_AntiKt4LCTopo_constscale_pt, &b_jet_AntiKt4LCTopo_constscale_pt);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_constscale_m", &jet_AntiKt4LCTopo_constscale_m, &b_jet_AntiKt4LCTopo_constscale_m);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_constscale_eta", &jet_AntiKt4LCTopo_constscale_eta, &b_jet_AntiKt4LCTopo_constscale_eta);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_constscale_phi", &jet_AntiKt4LCTopo_constscale_phi, &b_jet_AntiKt4LCTopo_constscale_phi);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_CENTER_LAMBDA", &jet_AntiKt4LCTopo_CENTER_LAMBDA, &b_jet_AntiKt4LCTopo_CENTER_LAMBDA);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_SECOND_R", &jet_AntiKt4LCTopo_SECOND_R, &b_jet_AntiKt4LCTopo_SECOND_R);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_SECOND_LAMBDA", &jet_AntiKt4LCTopo_SECOND_LAMBDA, &b_jet_AntiKt4LCTopo_SECOND_LAMBDA);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_CLUSTER_PT", &jet_AntiKt4LCTopo_CLUSTER_PT, &b_jet_AntiKt4LCTopo_CLUSTER_PT);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_Comb", &jet_AntiKt4LCTopo_flavor_weight_Comb, &b_jet_AntiKt4LCTopo_flavor_weight_Comb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_GbbNN", &jet_AntiKt4LCTopo_flavor_weight_GbbNN, &b_jet_AntiKt4LCTopo_flavor_weight_GbbNN);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_IP2D", &jet_AntiKt4LCTopo_flavor_weight_IP2D, &b_jet_AntiKt4LCTopo_flavor_weight_IP2D);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_IP3D", &jet_AntiKt4LCTopo_flavor_weight_IP3D, &b_jet_AntiKt4LCTopo_flavor_weight_IP3D);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN", &jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN, &b_jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_JetFitterCharm", &jet_AntiKt4LCTopo_flavor_weight_JetFitterCharm, &b_jet_AntiKt4LCTopo_flavor_weight_JetFitterCharm);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN", &jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN, &b_jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_MV1", &jet_AntiKt4LCTopo_flavor_weight_MV1, &b_jet_AntiKt4LCTopo_flavor_weight_MV1);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_SV0", &jet_AntiKt4LCTopo_flavor_weight_SV0, &b_jet_AntiKt4LCTopo_flavor_weight_SV0);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_SV1", &jet_AntiKt4LCTopo_flavor_weight_SV1, &b_jet_AntiKt4LCTopo_flavor_weight_SV1);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_SV2", &jet_AntiKt4LCTopo_flavor_weight_SV2, &b_jet_AntiKt4LCTopo_flavor_weight_SV2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_SecondSoftMuonTagChi2", &jet_AntiKt4LCTopo_flavor_weight_SecondSoftMuonTagChi2, &b_jet_AntiKt4LCTopo_flavor_weight_SecondSoftMuonTagChi2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_SoftMuonTagChi2", &jet_AntiKt4LCTopo_flavor_weight_SoftMuonTagChi2, &b_jet_AntiKt4LCTopo_flavor_weight_SoftMuonTagChi2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_truth_label", &jet_AntiKt4LCTopo_flavor_truth_label, &b_jet_AntiKt4LCTopo_flavor_truth_label);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_truth_dRminToB", &jet_AntiKt4LCTopo_flavor_truth_dRminToB, &b_jet_AntiKt4LCTopo_flavor_truth_dRminToB);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_truth_dRminToC", &jet_AntiKt4LCTopo_flavor_truth_dRminToC, &b_jet_AntiKt4LCTopo_flavor_truth_dRminToC);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_truth_dRminToT", &jet_AntiKt4LCTopo_flavor_truth_dRminToT, &b_jet_AntiKt4LCTopo_flavor_truth_dRminToT);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_truth_BHadronpdg", &jet_AntiKt4LCTopo_flavor_truth_BHadronpdg, &b_jet_AntiKt4LCTopo_flavor_truth_BHadronpdg);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_truth_vx_x", &jet_AntiKt4LCTopo_flavor_truth_vx_x, &b_jet_AntiKt4LCTopo_flavor_truth_vx_x);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_truth_vx_y", &jet_AntiKt4LCTopo_flavor_truth_vx_y, &b_jet_AntiKt4LCTopo_flavor_truth_vx_y);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_truth_vx_z", &jet_AntiKt4LCTopo_flavor_truth_vx_z, &b_jet_AntiKt4LCTopo_flavor_truth_vx_z);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip2d_pu", &jet_AntiKt4LCTopo_flavor_component_ip2d_pu, &b_jet_AntiKt4LCTopo_flavor_component_ip2d_pu);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip2d_pb", &jet_AntiKt4LCTopo_flavor_component_ip2d_pb, &b_jet_AntiKt4LCTopo_flavor_component_ip2d_pb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip2d_isValid", &jet_AntiKt4LCTopo_flavor_component_ip2d_isValid, &b_jet_AntiKt4LCTopo_flavor_component_ip2d_isValid);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip2d_ntrk", &jet_AntiKt4LCTopo_flavor_component_ip2d_ntrk, &b_jet_AntiKt4LCTopo_flavor_component_ip2d_ntrk);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip3d_pu", &jet_AntiKt4LCTopo_flavor_component_ip3d_pu, &b_jet_AntiKt4LCTopo_flavor_component_ip3d_pu);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip3d_pb", &jet_AntiKt4LCTopo_flavor_component_ip3d_pb, &b_jet_AntiKt4LCTopo_flavor_component_ip3d_pb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip3d_pc", &jet_AntiKt4LCTopo_flavor_component_ip3d_pc, &b_jet_AntiKt4LCTopo_flavor_component_ip3d_pc);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip3d_isValid", &jet_AntiKt4LCTopo_flavor_component_ip3d_isValid, &b_jet_AntiKt4LCTopo_flavor_component_ip3d_isValid);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk", &jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk, &b_jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv1_pu", &jet_AntiKt4LCTopo_flavor_component_sv1_pu, &b_jet_AntiKt4LCTopo_flavor_component_sv1_pu);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv1_pb", &jet_AntiKt4LCTopo_flavor_component_sv1_pb, &b_jet_AntiKt4LCTopo_flavor_component_sv1_pb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv1_pc", &jet_AntiKt4LCTopo_flavor_component_sv1_pc, &b_jet_AntiKt4LCTopo_flavor_component_sv1_pc);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv2_pu", &jet_AntiKt4LCTopo_flavor_component_sv2_pu, &b_jet_AntiKt4LCTopo_flavor_component_sv2_pu);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv2_pb", &jet_AntiKt4LCTopo_flavor_component_sv2_pb, &b_jet_AntiKt4LCTopo_flavor_component_sv2_pb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_pu", &jet_AntiKt4LCTopo_flavor_component_jfit_pu, &b_jet_AntiKt4LCTopo_flavor_component_jfit_pu);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_pb", &jet_AntiKt4LCTopo_flavor_component_jfit_pb, &b_jet_AntiKt4LCTopo_flavor_component_jfit_pb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_pc", &jet_AntiKt4LCTopo_flavor_component_jfit_pc, &b_jet_AntiKt4LCTopo_flavor_component_jfit_pc);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitcomb_pu", &jet_AntiKt4LCTopo_flavor_component_jfitcomb_pu, &b_jet_AntiKt4LCTopo_flavor_component_jfitcomb_pu);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitcomb_pb", &jet_AntiKt4LCTopo_flavor_component_jfitcomb_pb, &b_jet_AntiKt4LCTopo_flavor_component_jfitcomb_pb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitcomb_pc", &jet_AntiKt4LCTopo_flavor_component_jfitcomb_pc, &b_jet_AntiKt4LCTopo_flavor_component_jfitcomb_pc);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_nvtx", &jet_AntiKt4LCTopo_flavor_component_jfit_nvtx, &b_jet_AntiKt4LCTopo_flavor_component_jfit_nvtx);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t", &jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t, &b_jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_ntrkAtVx", &jet_AntiKt4LCTopo_flavor_component_jfit_ntrkAtVx, &b_jet_AntiKt4LCTopo_flavor_component_jfit_ntrkAtVx);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_efrc", &jet_AntiKt4LCTopo_flavor_component_jfit_efrc, &b_jet_AntiKt4LCTopo_flavor_component_jfit_efrc);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_mass", &jet_AntiKt4LCTopo_flavor_component_jfit_mass, &b_jet_AntiKt4LCTopo_flavor_component_jfit_mass);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_sig3d", &jet_AntiKt4LCTopo_flavor_component_jfit_sig3d, &b_jet_AntiKt4LCTopo_flavor_component_jfit_sig3d);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_deltaPhi", &jet_AntiKt4LCTopo_flavor_component_jfit_deltaPhi, &b_jet_AntiKt4LCTopo_flavor_component_jfit_deltaPhi);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_deltaEta", &jet_AntiKt4LCTopo_flavor_component_jfit_deltaEta, &b_jet_AntiKt4LCTopo_flavor_component_jfit_deltaEta);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropName", &jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropName, &b_jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropName);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropValue", &jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropValue, &b_jet_AntiKt4LCTopo_flavor_component_jfitc_doublePropValue);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitc_intPropName", &jet_AntiKt4LCTopo_flavor_component_jfitc_intPropName, &b_jet_AntiKt4LCTopo_flavor_component_jfitc_intPropName);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitc_intPropValue", &jet_AntiKt4LCTopo_flavor_component_jfitc_intPropValue, &b_jet_AntiKt4LCTopo_flavor_component_jfitc_intPropValue);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitc_pu", &jet_AntiKt4LCTopo_flavor_component_jfitc_pu, &b_jet_AntiKt4LCTopo_flavor_component_jfitc_pu);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitc_pb", &jet_AntiKt4LCTopo_flavor_component_jfitc_pb, &b_jet_AntiKt4LCTopo_flavor_component_jfitc_pb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitc_pc", &jet_AntiKt4LCTopo_flavor_component_jfitc_pc, &b_jet_AntiKt4LCTopo_flavor_component_jfitc_pc);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ipplus_trk_n", &jet_AntiKt4LCTopo_flavor_component_ipplus_trk_n, &b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_n);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ipplus_trk_index", &jet_AntiKt4LCTopo_flavor_component_ipplus_trk_index, &b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_index);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0val", &jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0val, &b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0val);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0sig", &jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0sig, &b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_d0sig);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0val", &jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0val, &b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0val);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0sig", &jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0sig, &b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_z0sig);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w2D", &jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w2D, &b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w2D);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w3D", &jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w3D, &b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_w3D);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ipplus_trk_grade", &jet_AntiKt4LCTopo_flavor_component_ipplus_trk_grade, &b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_grade);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ipplus_trk_isFromV0", &jet_AntiKt4LCTopo_flavor_component_ipplus_trk_isFromV0, &b_jet_AntiKt4LCTopo_flavor_component_ipplus_trk_isFromV0);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_isValid", &jet_AntiKt4LCTopo_flavor_component_svp_isValid, &b_jet_AntiKt4LCTopo_flavor_component_svp_isValid);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_ntrkv", &jet_AntiKt4LCTopo_flavor_component_svp_ntrkv, &b_jet_AntiKt4LCTopo_flavor_component_svp_ntrkv);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_n2t", &jet_AntiKt4LCTopo_flavor_component_svp_n2t, &b_jet_AntiKt4LCTopo_flavor_component_svp_n2t);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_mass", &jet_AntiKt4LCTopo_flavor_component_svp_mass, &b_jet_AntiKt4LCTopo_flavor_component_svp_mass);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_efrc", &jet_AntiKt4LCTopo_flavor_component_svp_efrc, &b_jet_AntiKt4LCTopo_flavor_component_svp_efrc);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_x", &jet_AntiKt4LCTopo_flavor_component_svp_x, &b_jet_AntiKt4LCTopo_flavor_component_svp_x);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_y", &jet_AntiKt4LCTopo_flavor_component_svp_y, &b_jet_AntiKt4LCTopo_flavor_component_svp_y);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_z", &jet_AntiKt4LCTopo_flavor_component_svp_z, &b_jet_AntiKt4LCTopo_flavor_component_svp_z);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_err_x", &jet_AntiKt4LCTopo_flavor_component_svp_err_x, &b_jet_AntiKt4LCTopo_flavor_component_svp_err_x);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_err_y", &jet_AntiKt4LCTopo_flavor_component_svp_err_y, &b_jet_AntiKt4LCTopo_flavor_component_svp_err_y);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_err_z", &jet_AntiKt4LCTopo_flavor_component_svp_err_z, &b_jet_AntiKt4LCTopo_flavor_component_svp_err_z);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_cov_xy", &jet_AntiKt4LCTopo_flavor_component_svp_cov_xy, &b_jet_AntiKt4LCTopo_flavor_component_svp_cov_xy);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_cov_xz", &jet_AntiKt4LCTopo_flavor_component_svp_cov_xz, &b_jet_AntiKt4LCTopo_flavor_component_svp_cov_xz);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_cov_yz", &jet_AntiKt4LCTopo_flavor_component_svp_cov_yz, &b_jet_AntiKt4LCTopo_flavor_component_svp_cov_yz);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_chi2", &jet_AntiKt4LCTopo_flavor_component_svp_chi2, &b_jet_AntiKt4LCTopo_flavor_component_svp_chi2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_ndof", &jet_AntiKt4LCTopo_flavor_component_svp_ndof, &b_jet_AntiKt4LCTopo_flavor_component_svp_ndof);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_ntrk", &jet_AntiKt4LCTopo_flavor_component_svp_ntrk, &b_jet_AntiKt4LCTopo_flavor_component_svp_ntrk);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_trk_n", &jet_AntiKt4LCTopo_flavor_component_svp_trk_n, &b_jet_AntiKt4LCTopo_flavor_component_svp_trk_n);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_trk_index", &jet_AntiKt4LCTopo_flavor_component_svp_trk_index, &b_jet_AntiKt4LCTopo_flavor_component_svp_trk_index);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_isValid", &jet_AntiKt4LCTopo_flavor_component_sv0p_isValid, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_isValid);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv", &jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_n2t", &jet_AntiKt4LCTopo_flavor_component_sv0p_n2t, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_n2t);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_mass", &jet_AntiKt4LCTopo_flavor_component_sv0p_mass, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_mass);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_efrc", &jet_AntiKt4LCTopo_flavor_component_sv0p_efrc, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_efrc);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_x", &jet_AntiKt4LCTopo_flavor_component_sv0p_x, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_x);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_y", &jet_AntiKt4LCTopo_flavor_component_sv0p_y, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_y);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_z", &jet_AntiKt4LCTopo_flavor_component_sv0p_z, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_z);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_err_x", &jet_AntiKt4LCTopo_flavor_component_sv0p_err_x, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_err_x);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_err_y", &jet_AntiKt4LCTopo_flavor_component_sv0p_err_y, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_err_y);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_err_z", &jet_AntiKt4LCTopo_flavor_component_sv0p_err_z, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_err_z);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xy", &jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xy, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xy);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xz", &jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xz, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_cov_xz);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_cov_yz", &jet_AntiKt4LCTopo_flavor_component_sv0p_cov_yz, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_cov_yz);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_chi2", &jet_AntiKt4LCTopo_flavor_component_sv0p_chi2, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_chi2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_ndof", &jet_AntiKt4LCTopo_flavor_component_sv0p_ndof, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_ndof);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk", &jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_trk_n", &jet_AntiKt4LCTopo_flavor_component_sv0p_trk_n, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_trk_n);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_trk_index", &jet_AntiKt4LCTopo_flavor_component_sv0p_trk_index, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_trk_index);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_n", &jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_n, &b_jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_n);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_index", &jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_index, &b_jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_index);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_w", &jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_w, &b_jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_w);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_pTRel", &jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_pTRel, &b_jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_pTRel);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_dRJet", &jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_dRJet, &b_jet_AntiKt4LCTopo_flavor_component_softmuonchi2info_muon_dRJet);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_msvp_isValid", &jet_AntiKt4LCTopo_flavor_component_msvp_isValid, &b_jet_AntiKt4LCTopo_flavor_component_msvp_isValid);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_msvp_ntrkv", &jet_AntiKt4LCTopo_flavor_component_msvp_ntrkv, &b_jet_AntiKt4LCTopo_flavor_component_msvp_ntrkv);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_msvp_n2t", &jet_AntiKt4LCTopo_flavor_component_msvp_n2t, &b_jet_AntiKt4LCTopo_flavor_component_msvp_n2t);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_msvp_nVtx", &jet_AntiKt4LCTopo_flavor_component_msvp_nVtx, &b_jet_AntiKt4LCTopo_flavor_component_msvp_nVtx);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_msvp_normWeightedDist", &jet_AntiKt4LCTopo_flavor_component_msvp_normWeightedDist, &b_jet_AntiKt4LCTopo_flavor_component_msvp_normWeightedDist);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_n", &jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_n, &b_jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_n);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_index", &jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_index, &b_jet_AntiKt4LCTopo_flavor_component_msvp_msvinjet_index);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitvx_theta", &jet_AntiKt4LCTopo_flavor_component_jfitvx_theta, &b_jet_AntiKt4LCTopo_flavor_component_jfitvx_theta);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitvx_phi", &jet_AntiKt4LCTopo_flavor_component_jfitvx_phi, &b_jet_AntiKt4LCTopo_flavor_component_jfitvx_phi);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitvx_errtheta", &jet_AntiKt4LCTopo_flavor_component_jfitvx_errtheta, &b_jet_AntiKt4LCTopo_flavor_component_jfitvx_errtheta);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitvx_errphi", &jet_AntiKt4LCTopo_flavor_component_jfitvx_errphi, &b_jet_AntiKt4LCTopo_flavor_component_jfitvx_errphi);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitvx_chi2", &jet_AntiKt4LCTopo_flavor_component_jfitvx_chi2, &b_jet_AntiKt4LCTopo_flavor_component_jfitvx_chi2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfitvx_ndof", &jet_AntiKt4LCTopo_flavor_component_jfitvx_ndof, &b_jet_AntiKt4LCTopo_flavor_component_jfitvx_ndof);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_n", &jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_n, &b_jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_n);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_index", &jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_index, &b_jet_AntiKt4LCTopo_flavor_component_jftwotrackvertex_index);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_bbDecision", &jet_AntiKt4LCTopo_bbDecision, &b_jet_AntiKt4LCTopo_bbDecision);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_TrackMFindex", &jet_AntiKt4LCTopo_TrackMFindex, &b_jet_AntiKt4LCTopo_TrackMFindex);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_TrackMF", &jet_AntiKt4LCTopo_TrackMF, &b_jet_AntiKt4LCTopo_TrackMF);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_nTrk_pv0_1GeV", &jet_AntiKt4LCTopo_nTrk_pv0_1GeV, &b_jet_AntiKt4LCTopo_nTrk_pv0_1GeV);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_sumPtTrk_pv0_1GeV", &jet_AntiKt4LCTopo_sumPtTrk_pv0_1GeV, &b_jet_AntiKt4LCTopo_sumPtTrk_pv0_1GeV);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_nTrk_allpv_1GeV", &jet_AntiKt4LCTopo_nTrk_allpv_1GeV, &b_jet_AntiKt4LCTopo_nTrk_allpv_1GeV);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_sumPtTrk_allpv_1GeV", &jet_AntiKt4LCTopo_sumPtTrk_allpv_1GeV, &b_jet_AntiKt4LCTopo_sumPtTrk_allpv_1GeV);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_nTrk_pv0_500MeV", &jet_AntiKt4LCTopo_nTrk_pv0_500MeV, &b_jet_AntiKt4LCTopo_nTrk_pv0_500MeV);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV", &jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV, &b_jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_trackWIDTH_pv0_1GeV", &jet_AntiKt4LCTopo_trackWIDTH_pv0_1GeV, &b_jet_AntiKt4LCTopo_trackWIDTH_pv0_1GeV);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_trackWIDTH_allpv_1GeV", &jet_AntiKt4LCTopo_trackWIDTH_allpv_1GeV, &b_jet_AntiKt4LCTopo_trackWIDTH_allpv_1GeV);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_LikeLihood_0", &jet_AntiKt4LCTopo_LikeLihood_0, &b_jet_AntiKt4LCTopo_LikeLihood_0);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_LowEtConstituentsFrac", &jet_AntiKt4LCTopo_LowEtConstituentsFrac, &b_jet_AntiKt4LCTopo_LowEtConstituentsFrac);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_KtDr", &jet_AntiKt4LCTopo_KtDr, &b_jet_AntiKt4LCTopo_KtDr);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_SPLIT12", &jet_AntiKt4LCTopo_SPLIT12, &b_jet_AntiKt4LCTopo_SPLIT12);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_SPLIT23", &jet_AntiKt4LCTopo_SPLIT23, &b_jet_AntiKt4LCTopo_SPLIT23);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_SPLIT34", &jet_AntiKt4LCTopo_SPLIT34, &b_jet_AntiKt4LCTopo_SPLIT34);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_truthAssoc_index", &jet_AntiKt4LCTopo_truthAssoc_index, &b_jet_AntiKt4LCTopo_truthAssoc_index);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_constit_n", &jet_AntiKt4LCTopo_constit_n, &b_jet_AntiKt4LCTopo_constit_n);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_constit_index", &jet_AntiKt4LCTopo_constit_index, &b_jet_AntiKt4LCTopo_constit_index);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_TrackAssoc_n", &jet_AntiKt4LCTopo_TrackAssoc_n, &b_jet_AntiKt4LCTopo_TrackAssoc_n);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_TrackAssoc_index", &jet_AntiKt4LCTopo_TrackAssoc_index, &b_jet_AntiKt4LCTopo_TrackAssoc_index);
   fChain->SetBranchAddress("trk_n", &trk_n, &b_trk_n);
   fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
   fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
   fChain->SetBranchAddress("trk_d0_wrtPV", &trk_d0_wrtPV, &b_trk_d0_wrtPV);
   fChain->SetBranchAddress("trk_z0_wrtPV", &trk_z0_wrtPV, &b_trk_z0_wrtPV);
   fChain->SetBranchAddress("trk_phi_wrtPV", &trk_phi_wrtPV, &b_trk_phi_wrtPV);
   fChain->SetBranchAddress("trk_theta_wrtPV", &trk_theta_wrtPV, &b_trk_theta_wrtPV);
   fChain->SetBranchAddress("trk_qoverp_wrtPV", &trk_qoverp_wrtPV, &b_trk_qoverp_wrtPV);
   fChain->SetBranchAddress("trk_cov_d0_wrtPV", &trk_cov_d0_wrtPV, &b_trk_cov_d0_wrtPV);
   fChain->SetBranchAddress("trk_cov_z0_wrtPV", &trk_cov_z0_wrtPV, &b_trk_cov_z0_wrtPV);
   fChain->SetBranchAddress("trk_cov_phi_wrtPV", &trk_cov_phi_wrtPV, &b_trk_cov_phi_wrtPV);
   fChain->SetBranchAddress("trk_cov_theta_wrtPV", &trk_cov_theta_wrtPV, &b_trk_cov_theta_wrtPV);
   fChain->SetBranchAddress("trk_cov_qoverp_wrtPV", &trk_cov_qoverp_wrtPV, &b_trk_cov_qoverp_wrtPV);
   fChain->SetBranchAddress("trk_d0_wrtBL", &trk_d0_wrtBL, &b_trk_d0_wrtBL);
   fChain->SetBranchAddress("trk_d0_err_wrtBL", &trk_d0_err_wrtBL, &b_trk_d0_err_wrtBL);
   fChain->SetBranchAddress("trk_eta_atCalo", &trk_eta_atCalo, &b_trk_eta_atCalo);
   fChain->SetBranchAddress("trk_phi_atCalo", &trk_phi_atCalo, &b_trk_phi_atCalo);
   fChain->SetBranchAddress("trk_eta_atCalo2ndLayer", &trk_eta_atCalo2ndLayer, &b_trk_eta_atCalo2ndLayer);
   fChain->SetBranchAddress("trk_phi_atCalo2ndLayer", &trk_phi_atCalo2ndLayer, &b_trk_phi_atCalo2ndLayer);
   fChain->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
   fChain->SetBranchAddress("trk_ndof", &trk_ndof, &b_trk_ndof);
   fChain->SetBranchAddress("trk_nBLHits", &trk_nBLHits, &b_trk_nBLHits);
   fChain->SetBranchAddress("trk_nPixHits", &trk_nPixHits, &b_trk_nPixHits);
   fChain->SetBranchAddress("trk_nSCTHits", &trk_nSCTHits, &b_trk_nSCTHits);
   fChain->SetBranchAddress("trk_nTRTHits", &trk_nTRTHits, &b_trk_nTRTHits);
   fChain->SetBranchAddress("trk_nTRTHighTHits", &trk_nTRTHighTHits, &b_trk_nTRTHighTHits);
   fChain->SetBranchAddress("trk_nPixHoles", &trk_nPixHoles, &b_trk_nPixHoles);
   fChain->SetBranchAddress("trk_nSCTHoles", &trk_nSCTHoles, &b_trk_nSCTHoles);
   fChain->SetBranchAddress("trk_nTRTHoles", &trk_nTRTHoles, &b_trk_nTRTHoles);
   fChain->SetBranchAddress("trk_nPixelDeadSensors", &trk_nPixelDeadSensors, &b_trk_nPixelDeadSensors);
   fChain->SetBranchAddress("trk_nSCTDeadSensors", &trk_nSCTDeadSensors, &b_trk_nSCTDeadSensors);
   fChain->SetBranchAddress("trk_expectBLayerHit", &trk_expectBLayerHit, &b_trk_expectBLayerHit);
   fChain->SetBranchAddress("trk_mc_probability", &trk_mc_probability, &b_trk_mc_probability);
   fChain->SetBranchAddress("trk_mc_barcode", &trk_mc_barcode, &b_trk_mc_barcode);
   fChain->SetBranchAddress("trk_mc_index", &trk_mc_index, &b_trk_mc_index);
   fChain->SetBranchAddress("vxp_n", &vxp_n, &b_vxp_n);
   fChain->SetBranchAddress("vxp_x", &vxp_x, &b_vxp_x);
   fChain->SetBranchAddress("vxp_y", &vxp_y, &b_vxp_y);
   fChain->SetBranchAddress("vxp_z", &vxp_z, &b_vxp_z);
   fChain->SetBranchAddress("vxp_cov_x", &vxp_cov_x, &b_vxp_cov_x);
   fChain->SetBranchAddress("vxp_cov_y", &vxp_cov_y, &b_vxp_cov_y);
   fChain->SetBranchAddress("vxp_cov_z", &vxp_cov_z, &b_vxp_cov_z);
   fChain->SetBranchAddress("vxp_cov_xy", &vxp_cov_xy, &b_vxp_cov_xy);
   fChain->SetBranchAddress("vxp_cov_xz", &vxp_cov_xz, &b_vxp_cov_xz);
   fChain->SetBranchAddress("vxp_cov_yz", &vxp_cov_yz, &b_vxp_cov_yz);
   fChain->SetBranchAddress("vxp_type", &vxp_type, &b_vxp_type);
   fChain->SetBranchAddress("vxp_chi2", &vxp_chi2, &b_vxp_chi2);
   fChain->SetBranchAddress("vxp_ndof", &vxp_ndof, &b_vxp_ndof);
   fChain->SetBranchAddress("vxp_px", &vxp_px, &b_vxp_px);
   fChain->SetBranchAddress("vxp_py", &vxp_py, &b_vxp_py);
   fChain->SetBranchAddress("vxp_pz", &vxp_pz, &b_vxp_pz);
   fChain->SetBranchAddress("vxp_E", &vxp_E, &b_vxp_E);
   fChain->SetBranchAddress("vxp_m", &vxp_m, &b_vxp_m);
   fChain->SetBranchAddress("vxp_nTracks", &vxp_nTracks, &b_vxp_nTracks);
   fChain->SetBranchAddress("vxp_sumPt", &vxp_sumPt, &b_vxp_sumPt);
   fChain->SetBranchAddress("vxp_trk_weight", &vxp_trk_weight, &b_vxp_trk_weight);
   fChain->SetBranchAddress("vxp_trk_n", &vxp_trk_n, &b_vxp_trk_n);
   fChain->SetBranchAddress("vxp_trk_index", &vxp_trk_index, &b_vxp_trk_index);
   Notify();
}

Bool_t NtupCommonReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NtupCommonReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NtupCommonReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NtupCommonReader_cxx
