void JVTTest(){

    gSystem->CompileMacro("Loader.C");
    gSystem->CompileMacro("NtupCommonReader.C");
    gSystem->Load("../StandAlone/libJetVertexTagger.so");
    JetVertexTagger *jvt = new JetVertexTagger(-1, "../data/JVTlikelihood_20140805.root");

    TFile *fT = new TFile("mc12_8TeV.147771.Sherpa_CT10_Zmumu.merge.NTUP_COMMON.e1434_s1499_s1504_r3658_r3549_p1562.NTUP_COMMON.01313919._000880.skimmed100Events.root");
    TTree *tree = (TTree*) fT->Get("physics");
    NtupCommonReader reader(tree);

    Long64_t nentries = reader.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    // loop over 200 events 
    for (Long64_t jentry=0; jentry<100;jentry++) {
      Long64_t ientry = reader.LoadTree(jentry);
      if (ientry < 0) break;
      nb = reader.fChain->GetEntry(jentry);   nbytes += nb;


      vector<float> trk_pt                           = *(reader.trk_pt);
      vector<float> trk_z0_wrtPV                     = *(reader.trk_z0_wrtPV);
      vector<vector<int> > vxp_trk_index             = *(reader.vxp_trk_index);
      vector<vector<int> > assoc_trk_indices_alljets = *(reader.jet_AntiKt4LCTopo_TrackAssoc_index);

      // trk pt needs to be in GeV!!!
      for(int it=0; it<trk_pt.size(); ++it){trk_pt[it] = trk_pt[it]/1000.;}

      // initialize jvt-tool
      jvt->init(trk_pt, trk_z0_wrtPV, vxp_trk_index); 

      // loop over jets
      for(int ij=0; ij<reader.jet_AntiKt4LCTopo_n; ++ij){

         // associated tracks
         vector<int> assoc_trk_indices = assoc_trk_indices_alljets[ij];

         // pt of current jet
         // !!!! WARNING: for an analysis usage, use the fully corrected jet pT and not the constscale pT, needs to be in GeV
         float jetpt  = (reader.jet_AntiKt4LCTopo_constscale_pt)->at(ij) / 1000.;  // WARNING (see above)
         float jeteta = (reader.jet_AntiKt4LCTopo_eta)->at(ij);

         if( jetpt< 20 || jetpt > 50 || fabs(jeteta)>2.4 )    continue;

         // calculate JVT
         bool pass = (*jvt)(jetpt, assoc_trk_indices);
         cout << setprecision(3) << "JVTTest: >>> EventNumber " << (reader.EventNumber)  
                                 <<  " jet pt " << std::setw(10) <<  jetpt 
                                 <<  " eta "    << std::setw(10) <<  jeteta 
                                 << " jvt "     << std::setw(10) <<  jvt->JVT() 
                                 << " corrJVF " << std::setw(10) <<  jvt->corrJVF() 
                                 << " RpT "     << std::setw(10) <<  jvt->RpT()  << " <<< "  << endl;


      }
   }


   fT->Close();
   delete jvt;
   delete fT;
}
