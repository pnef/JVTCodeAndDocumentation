/**************************************************************************
 **
 **   File:         Analysis_JetMET_Base.h
 **
 **   Description:  See header
 **                 
 **   Authors:      M. Swiatlowski
 **
 **************************************************************************/

#define Analysis_JetMET_Base_cxx

#include "Analysis_JetMET_Base.h"
#include "AnaConfig.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <iostream>
#include <cstdlib>
#include "TKey.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "fastjet/tools/Pruner.hh"
#include "Nsubjettiness.h"
#include "MuonEfficiencyCorrections/AnalysisMuonConfigurableScaleFactors.h"
#include "CalibrationDataInterface/CalibrationDataInterfaceROOT.h"
#include "TrigMuonEfficiency/LeptonTriggerSF.h"
#include "SUSYTools/DataPeriod.h"
#include "EventBuilder_jetmet2012.h"

struct BTagIndex {
  unsigned int b[2];
  unsigned int c[2];
  unsigned int t[2];
  unsigned int l[2];
};

bool Analysis_JetMET_Base::ManuallyAddConstituents(const MomKey JetKey){
  const static MomKey jetl("jets");

  Sort(jetl+JetKey, "constscale_pt");

  TString JetKeyString(JetKey.Data());
  MomKey constType;

  if(JetKeyString.Contains("LCTopo")){
    constType = "clustersLCTopo";
  } else if (JetKeyString.Contains("Track")){
    constType = "tracksgood";
  } else if (JetKeyString.Contains("Truth")){
    constType = "truthsStable";
  } else {
    cout << " THIS JET TYPE IS NOT SUPPORTED FOR MANUAL CONSTITUENTS " << JetKey << endl;
    return false;
  }

  vector<fastjet::PseudoJet> inputConst = ObjsToPJ(constType);

  float jetR = 0.4;

  if(JetKeyString.Contains("4")){
    jetR = 0.4;
  } else if (JetKeyString.Contains("6")){
    jetR = 0.6;
  } else {
    cout << "THIS JET RAD CURRENTLY NOT SUPPORTED " << JetKey << endl; 
    return false;
  }

  fastjet::JetAlgorithm algo;

  if(JetKeyString.Contains("AntiKt")){
    algo = fastjet::antikt_algorithm;
  } else if (JetKeyString.Contains("CamKt")){
    algo = fastjet::cambridge_algorithm;
  } else if (JetKeyString.Contains("Kt")){
    algo = fastjet::kt_algorithm;
  } else {
    cout << "THIS JET ALGO CURRENTLY NOT SUPPORTED " << JetKey << endl;
    return false;
  }

  fastjet::JetDefinition jetDef(algo, jetR,fastjet::E_scheme, fastjet::Best);
  fastjet::ClusterSequence clustSeq(inputConst, jetDef);  

  vector<fastjet::PseudoJet> inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets());

  const static MomKey mConstituents("constituents");

  for(int iJet = 0; iJet < inclusiveJets.size() && iJet < jets(JetKey); iJet++){
    if(Debug()){
      cout << "testing compatibility of jets" << endl;
      cout << inclusiveJets[iJet].pt() << " " << inclusiveJets[iJet].eta() << endl;
      cout << jet(iJet, JetKey).Float("constscale_pt") << " " << jet(iJet, JetKey).Float("constscale_eta") << endl;
    }
    jet(iJet, JetKey).AddVec(mConstituents);
    vector<fastjet::PseudoJet> tempConst = inclusiveJets[iJet].constituents();
    for(int iC = 0; iC < tempConst.size(); iC++){
      jet(iJet, JetKey).Add(mConstituents, tempConst[iC].user_info<PJ_Info>().Pointer);
    }
    if(Debug()) cout << " added " << tempConst.size() << " constituents to the " << iJet << "th jet" << endl;
  }

  SortPt(jetl+JetKey);
  return true;
}

bool Analysis_JetMET_Base::SumMatchedCollectionJets(const MomKey JetKey, const MomKey collection){
  for(int iJet = 0; iJet < jets(JetKey); iJet++){
    SumMatchedCollection(&jet(iJet, JetKey), collection);
  }
  return true;
}

bool Analysis_JetMET_Base::SumMatchedCollection(Particle* theJet, const MomKey collection){
  const MomKey sumKey = collection + MomKey("Summed");
  theJet->AddVec(sumKey);
  TLorentzVector pSum(0,0,0,0); 
  for(int iC = 0; iC < theJet->Objs(collection); iC++){
    pSum += ((Particle*) theJet->Obj(collection, iC))->p;
  }
  Particle* theSum = new Particle(pSum);
  theJet->Add(sumKey, theSum);
  return true;
}

bool Analysis_JetMET_Base::ReduceTruthBHadronCollection(const MomKey TruthKey){
  vector<vector<double> > links;

  // Prepare Link Vector -- links[i][j] represents the link between particle i and particle (i+1+j)
  // Compute the average on the flight
  double average = 0.;

  for(int i = 0; i < truths(TruthKey); i++){
    vector<double> row_i;
    Particle* p_i = &truth(i, TruthKey);

    for(int j = i+1; j < truths(TruthKey); j++){
      Particle* p_j = &truth(j, TruthKey);
      double dR = p_i->p.DeltaR(p_j->p);

      row_i.push_back(dR);
      average += dR;
    }

    links.push_back(row_i);
  }

  average = average/(truths(TruthKey)*(truths(TruthKey)-1)/2);

  // Initialize status of each input particle
  vector<int> inputs_status;
  for(unsigned int i = 0; i < truths(TruthKey); i++)
    inputs_status.push_back(1);

  // Go through the links again to decide boundary
  for(unsigned int i = 0; i < links.size(); i++){
    vector<double> links_row_i = links.at(i);
    Particle* p_i = &truth(i, TruthKey);

    for(unsigned int j = 0; j < links_row_i.size(); j++){
      Particle* p_j = &truth(i+1+j, TruthKey);
    
      if(links_row_i.at(j) >= average) continue;
      if(p_i->Int("status") == p_j->Int("status")) continue;
      //if(p_i->Int("pdgId") == (-1*p_j->Int("pdgId"))) continue; // abandon this because of B mixing(oscillation)

      // merge those below boundary -- remain the one with highest mc_status
      if(p_j->Int("status") > p_i->Int("status")){ 
        inputs_status.at(i) = 0;

        //cout << "Merging two particles:" << endl;
        //PrintParticleInfo(p_i, false);
        //PrintParticleInfo(p_j, false);

      }
      else{
        inputs_status.at(i+1+j) = 0;

        //cout << "Merging two particles: " << endl;
        //PrintParticleInfo(p_i, false);
        //PrintParticleInfo(p_j, false);
     }
    }
  } 

  static const MomKey mTruths("truths");
  // re-organize output
  for(int i = 0; i < truths(TruthKey); i++){
    if(inputs_status.at(i) == 0){
      if(Debug()) cout << " we are removing a b hadron! " << endl;
      Remove(mTruths+TruthKey, i);
      inputs_status.erase(inputs_status.begin()+ i);
      i--;
    }
  }

  /*/
  if(outputs.size() == 3){
    cout << "average dR: " << average << endl;
    for(unsigned int i = 0; i < outputs.size(); i++){
      PrintParticleInfo(outputs.at(i), false);
    }
  } */
  return true;
}

// label jets as b-jets.
// if there is b-hadron > 5 GeV inside it (via DR)
// boolean "IsBJet" will be set true; otherwise false.
bool Analysis_JetMET_Base::LabelBJetsDR(const MomKey JetKey){
  if(Debug()) cout << " In LabelBJetsDR " << endl;
  if(!isMC())
    return false;

  const static MomKey mBHadrons("BHadrons");
  const static MomKey mtruthsBHadrons("truthsBHadrons");
  const static MomKey mtruthsBHadronsDR("truthsBHadronsDR");
  const static MomKey IsBJetTruth("IsBJetTruth");

  if(Debug()) cout << " Initializing truth B hadron collection if necessary " << endl;
  if(!Exists(mtruthsBHadrons)){
    MakeTruthBHadronCollection();
    if(Debug()) cout << truths(mBHadrons) << " is the number of B hadrons before cleanup ";
    //ReduceTruthBHadronCollection();
    if(Debug()) cout << truths(mBHadrons) << " is the number after " << endl;
  }

  //if(truths(mBHadrons) == 0)
  //  return false;

  //typically want a 5 GeV cut, remove things with less than this
  for(int iTruth = 0; iTruth < truths(mBHadrons); iTruth++){
    if(truth(iTruth,mBHadrons).p.Perp() < 5.){
      Remove(mtruthsBHadrons, iTruth);
      iTruth--;
    }
  }

  AddDRMatch(JetKey, mtruthsBHadrons);

  for(int iJet = 0; iJet < jets(JetKey); iJet++){
    if(jet(iJet, JetKey).Objs(mtruthsBHadronsDR)){
      jet(iJet, JetKey).Set(IsBJetTruth, true);
    } else {
      jet(iJet, JetKey).Set(IsBJetTruth, false);
    }
  }

  return true;
}

///=========================================
/// Print truth info for 1 particle
///=========================================
void Analysis_JetMET_Base::PrintParticle(Particle* part){
  if(!isMC()) 
    return;

  if(!part)
    return;
    
  std::cout<<"Particle index ("<<part->Int("index")<<") {pT,Eta,Phi,M}={"<<part->p.Pt()<<","<<part->p.Eta()<<","<<part->p.Phi()<<","<<part->p.M()<<"} "<<"[barcode="<<part->Int("barcode")<<"] pdgId="<<part->Int("pdgId")<<" status="<<part->Int("status");
  
  std::cout<<" parents=";
  for(unsigned int iPar=0; iPar<part->ObjVec("parents").size(); iPar++)
    std::cout<<dynamic_cast<Particle*>(part->Obj("parents",iPar))->Int("index")<<"("<<dynamic_cast<Particle*>(part->Obj("parents",iPar))->Int("pdgId")<<"), ";
  
  std::cout<<" child=";
  for(unsigned int iChi=0; iChi<part->ObjVec("children").size(); iChi++)
    std::cout<<dynamic_cast<Particle*>(part->Obj("children",iChi))->Int("index")<<"("<< dynamic_cast<Particle*>(part->Obj("children",iChi))->Int("pdgId")<<","<<dynamic_cast<Particle*>(part->Obj("children",iChi))->Int("status")<<"), ";
  
  std::cout<<std::endl;
  
  return;

}


bool Analysis_JetMET_Base::GetBHadronDecayProducts(Particle* part, bool printFinalProducts, bool printTree)
{
  if(!isMC()) 
    return false;

  if(!part)
    return false;

  if(part->Exists("FinalDecayProducts"))
    return true;
  
  
  Particle* initial =  part;
  initial->AddVec("FinalDecayProducts", true);
  initial->AddVec("ChargedDecayProducts", true);


  std::vector< std::set<Particle*> > decay_tree;
  std::set< Particle* > start_node; 
  start_node.insert( initial );
  decay_tree.push_back( start_node );

  std::set< Particle* > final_particles;

  //int tree_level=0;
  std::set< Particle* > next_node = start_node;
  std::set< Particle* >::iterator it;
  bool allFinished=false;
  while(!allFinished)
    {
      
      std::set< Particle* > curr_node = next_node;
      next_node.clear();
      
      for(it=curr_node.begin(); it!=curr_node.end(); it++)
  {
    Particle* curr_part = (*it);

    if(curr_part->ObjVec("children").size()==0)
      {
        final_particles.insert(curr_part);
        continue;
      }
    
    for(unsigned int iChi=0; iChi<curr_part->ObjVec("children").size(); iChi++)
      {
        Particle* child = dynamic_cast<Particle*>(curr_part->Obj("children",iChi));
        next_node.insert(child);
      }//end loop over current tree level particle children 
  }//end loop over current tree level particles
      
      if(next_node.size() > 0)
  decay_tree.push_back(next_node);
      else
  allFinished=true;
      
    }//check if all done

  for(it=final_particles.begin(); it!=final_particles.end(); it++){
      initial->Add("FinalDecayProducts", (*it));
      if((*it)->Float("charge") != 0)
  initial->Add("ChargedDecayProducts", (*it));
  }

  float dr_ave = 0;
  float dr2_ave = 0;
  float npart_ave=0;

  float dr_ptave = 0;
  float dr2_ptave = 0;
  float sum_pt = 0;

  for(unsigned int iChd=0; iChd < initial->ObjVec("ChargedDecayProducts").size(); iChd++){
    Particle* chdProd = dynamic_cast<Particle*>(initial->Obj("ChargedDecayProducts",iChd));

    dr_ave += initial->p.DeltaR(chdProd->p);
    dr2_ave += initial->p.DeltaR(chdProd->p) * initial->p.DeltaR(chdProd->p);
    npart_ave += 1.0;

    dr_ptave += chdProd->p.Pt() * initial->p.DeltaR(chdProd->p);
    dr2_ptave += chdProd->p.Pt() * initial->p.DeltaR(chdProd->p) * initial->p.DeltaR(chdProd->p);
    sum_pt += chdProd->p.Pt();
  }
  if(npart_ave > 0){
    dr_ave /= npart_ave;
    dr2_ave /= npart_ave;

    initial->Set("ChargedDecay_DR_ave", dr_ave );
    initial->Set("ChargedDecay_DR_rms", sqrt(dr2_ave - dr_ave*dr_ave) );

    Fill("ChargedDecay_DR_ave", dr_ave, 1.0, 100, 0, 2);
    Fill("ChargedDecay_DR_rms", sqrt(dr2_ave - dr_ave*dr_ave), 1.0, 100, 0, 2);

  }
  if(sum_pt > 0){
    dr_ptave /= sum_pt;
    dr2_ptave /= sum_pt;
    
    initial->Set("ChargedDecay_DR_ptave", dr_ptave );
    initial->Set("ChargedDecay_DR_ptrms", sqrt(dr2_ptave - dr_ptave*dr_ptave) );

    Fill("ChargedDecay_DR_ptave", dr_ptave, 1.0, 100, 0, 2);
    Fill("ChargedDecay_DR_ptrms", sqrt(dr2_ptave - dr_ptave*dr_ptave), 1.0, 100, 0, 2);
  }

  
  if(printFinalProducts)
    {
      std::cout<<"-------------------------------------------------------"<<std::endl;
      std::cout<<"FINAL DECAY PRODUCTS (particle pdgId="<<part->Int("pdgId")<<" index = "<<part->Int("index")<<")"<<std::endl;
      std::cout<<"-------------------------------------------------------"<<std::endl;
      for(it=final_particles.begin(); it!=final_particles.end(); it++)
  PrintParticle( (*it) );

      std::cout<<"-------------------------------------------------------"<<std::endl;
      std::cout<<std::endl;
    }

  if(printTree)
    {
      std::cout<<"-------------------------------------------------------"<<std::endl;
      std::cout<<"DECAY TREE (particle pdgId="<<part->Int("pdgId")<<" index="<<part->Int("index")<<")"<<std::endl;
      std::cout<<"-------------------------------------------------------"<<std::endl;

      for(unsigned int i=0; i<decay_tree.size(); i++)
  {
    std::cout<<"LEVEL = "<<i<<std::endl;
    std::cout<<"-----------"<<std::endl;
    
    for(it=decay_tree.at(i).begin(); it!=decay_tree.at(i).end(); it++)
      PrintParticle( (*it) );
    
    std::cout<<"-------------------------------------------------------"<<std::endl;
    std::cout<<std::endl;
  }
    }
  

  return true;
}



// Make a truth bhadron collection directly //
bool Analysis_JetMET_Base::MakeTruthBHadronCollection()
{

  if(Debug()) cout << "Entering Analysis_JetMET_Base::MakeTruthBHadronCollection" << endl;
  
  if(!isMC()) 
    return false;
  
  if(truths()<=0){
    std::cout<<"Analysis_JetMET_Base::MakeTruthBHadronCollection - No Truth particles found???"<<std::endl;
    return false;
  }    

  if(Exists("truthsBHadrons") || Exists("truthsWeaklyDecayingBHadrons"))
    return true;

  
  AddVec("truthsBHadrons");
  AddVec("truthsWeaklyDecayingBHadrons");
  AddVec("truthsWeaklyDecayW");
  
  for(int iTruth=0; iTruth<truths(); iTruth++){
    if(Debug()) cout << "scanning " << iTruth << " of " << truths() << endl;
   
    int pdgId = truth(iTruth).Int("pdgId");
    //int status = truth(iTruth).Int("status");

    bool isBHadron= IsBHadron(pdgId);
    if(Debug()) cout << "Checking if this truth's pdgId is bhadron: " << isBHadron << endl;
    bool isWeaklyDecayingBHadron=IsWeaklyDecayingBHadron(pdgId);
    if(Debug()) cout << "Checking if this turth's pdgId is weakly decaying bhadron: " << isWeaklyDecayingBHadron << endl;

    if(!isBHadron && !isWeaklyDecayingBHadron)
      continue;

    bool isEndOfChain=true;
    bool childIsBHadron=false;
    for(unsigned int iChi=0; iChi<truth(iTruth).ObjVec("children").size(); iChi++){
      if(Debug()) cout << "Scanning child " << iChi << " of " << truth(iTruth).ObjVec("children").size() << endl;
    
      Particle* child = dynamic_cast<Particle*>(truth(iTruth).Obj("children",iChi));
      if(child->Int("pdgId") == pdgId){
        if(Debug()) cout << "Still not the end" << endl;
      
        isEndOfChain=false;
         break;
      }
      if(IsBHadron(child->Int("pdgId"))){
        if(Debug()) cout << "Before the end: a bhadron child" << endl;
      
        childIsBHadron=true;
      }
    }
    
    if(Debug()){
      cout << "End of child loop" << endl;
      cout << "isEndOfChain ? " << isEndOfChain << endl;
      cout << "childIsBHadron ? " << childIsBHadron << endl;
    }

    if(!isEndOfChain)
      continue;

    if(childIsBHadron)
      continue;

    bool isFromTop = IsBHadronFromTop(iTruth);

    if(Debug()) cout << "Does the BHadron come from top ? " << isFromTop << endl;

    truth(iTruth).Set("IsFromTop",isFromTop);

    if(isBHadron){
      if(Debug()) cout << "This is truth BHadron!" << endl;
      Add("truthsBHadrons",&truth(iTruth));
    }
    
    if(isWeaklyDecayingBHadron){
      if(Debug()) cout << "This is weakly decaying BHadron!" << endl;
      GetBHadronDecayProducts(&truth(iTruth), false, false);
      Add("truthsWeaklyDecayingBHadrons",&truth(iTruth));
    }
    
  }//end loop through truths
  return true;
  
}

///=========================================
/// truth B-hadron pdg id's
///=========================================
bool Analysis_JetMET_Base::IsBHadron(int pdgId)
{
  int abs_pdgId = abs(pdgId);
  int abs_pdgId_mod10k = (abs(pdgId)%10000);
  if( (abs_pdgId_mod10k>=500 && abs_pdgId_mod10k<600) /*mesons*/  ||
      (abs_pdgId>=5000      && abs_pdgId<6000)      /*baryons*/   )
    return true;

  return false;
}

///=========================================
/// truth weakly decaying B-hadron pdg id's
///=========================================
bool Analysis_JetMET_Base::IsWeaklyDecayingBHadron(int pdgId)
{
  /* From Gio Zevi Della Porta, from Gavin
     stableB.push_back(511);
     stableB.push_back(521);
     stableB.push_back(531);
     stableB.push_back(541);
     //stableB.push_back(551); //<- eta 1s, excluded!
     //stableB.push_back(553); //<-upsilon, excluded!
     stableB.push_back(5112);
     stableB.push_back(5122);
     stableB.push_back(5132);
     stableB.push_back(5142);
     stableB.push_back(5212);
     stableB.push_back(5222);
     stableB.push_back(5232);
     stableB.push_back(5242);
     stableB.push_back(5332);
  */
  int abs_pdgId = abs(pdgId);
  if( abs_pdgId==511  || abs_pdgId==521  || abs_pdgId==531  || abs_pdgId==541  || 
      abs_pdgId==5112 || abs_pdgId==5122 || abs_pdgId==5132 || abs_pdgId==5142 ||
      abs_pdgId==5212 || abs_pdgId==5222 || abs_pdgId==5232 || abs_pdgId==5242 ||
      abs_pdgId==5332)
    return true;

  return false;
}

bool Analysis_JetMET_Base::IsBHadronFromTop(int iTruth)
{
  if(Debug()) cout << "Entering Analysis_JetMET_Base::IsBHadronFromTop" << endl;

  if(!isMC())
    return false;
  
  if(truths()<=0){
    std::cout<<"Analysis_JetMET_Base::MakeTruthBHadronCollection - No Truth particles found???"<<std::endl;
    return false;
  }    

  int pdgId = truth(iTruth).Int("pdgId");
  if(!IsBHadron(pdgId)){
    return false;
  }

  //loop up through parents.  Parent must be b-hadron, b-quark, or transition particle (pdgId==91)
  Particle* curr_particle = &truth(iTruth);
  bool foundTop=false;
  truth(iTruth).AddVec("TopParent",true);

  //Particle* init = curr_particle;

  while(!foundTop)
    {
      bool goodParent=false;
      for(unsigned int iPar=0; iPar<curr_particle->ObjVec("parents").size(); iPar++)
  {
    Particle* parent = dynamic_cast<Particle*>(curr_particle->Obj("parents",iPar));
    int parent_pdgId = parent->Int("pdgId");
    
    //if(parent == init){ // MARK a potential bug that is very very rare abd bizzare (but happens!)
    //  cout << "Warning in Analysis_JetMET_Base::IsBHadronFromTop: The BHadron goes back to itself through parent chain!" << endl;
    //  goodParent = false;
    //break;
    //}
    
    if(IsBHadron(parent_pdgId) || abs(parent_pdgId)==5 || abs(parent_pdgId)==91){
      if(Debug()) cout << "Good parent found. Change current obj " << curr_particle << " to parent obj " << parent << endl;
    
      goodParent=true;
      curr_particle = parent;
      break;
    }
    if(abs(parent_pdgId)==6){
      foundTop  = true;
        truth(iTruth).Add("TopParent", parent);
      break;
    }
  }
      
      //If none of the parents are top and none of the parents are b-hadron related
      //then this particle can not come from the top.
      if(!foundTop && !goodParent)
        break;
    }
  
  if(Debug()) cout << "Leaving Analysis_JetMET_Base::IsBHadronFromTop" << endl;
  
  return foundTop;
}








MomKey Analysis_JetMET_Base::GetBTagKey(const int SystType){
  MomKey mSys;
  MomKey mBtagScaleFactor70 = MomKey("BtagScaleFactor70");

  switch(SystType) {
      case JetMETSyst::BD:
         mSys = MomKey("Bdown");

        mBtagScaleFactor70+=mSys;
        break;
      case JetMETSyst::CD:
         mSys = MomKey("Cdown");

        mBtagScaleFactor70+=mSys;
        break;
      case JetMETSyst::LD:
        mSys = MomKey("Ldown");

        mBtagScaleFactor70+=mSys;
        break;
      case JetMETSyst::BU:
        mSys = MomKey("Bup");

        mBtagScaleFactor70+=mSys;
        break;
      case JetMETSyst::CU:
        mSys = MomKey("Cup");

        mBtagScaleFactor70+=mSys;
        break;
      case JetMETSyst::LU:
        mSys = MomKey("Lup");

        mBtagScaleFactor70+=mSys;
        break;
      default:
        mSys = MomKey("Nom");

        mBtagScaleFactor70+=mSys;
        break;
    }

    return mBtagScaleFactor70;
}

MomKey Analysis_JetMET_Base::GetMuonSFKey(const int SystType){
  MomKey mSys; 
  MomKey m1MuScaleFactor("m1MuScaleFactor");
  switch(SystType){
    case JetMETSyst::MEFFUP:
      mSys = MomKey("Up");
      m1MuScaleFactor+=mSys;
      break;
    case JetMETSyst::MEFFDOWN:
      mSys = MomKey("Down");
      m1MuScaleFactor+=mSys;
      break;
    default:

      break;
  }
  return m1MuScaleFactor;
}

MomKey Analysis_JetMET_Base::GetMuonTrigSFKey(const int SystType){
  MomKey mSys; 
  MomKey m1MuTrigScaleFactor("m1MuTrigScaleFactor");
  switch(SystType){
    case JetMETSyst::MTRIGUP:
      mSys = MomKey("Up");
      m1MuTrigScaleFactor+=mSys;
      break;
    case JetMETSyst::MTRIGDOWN:
      mSys = MomKey("Down");
      m1MuTrigScaleFactor+=mSys;
      break;
    default:

      break;
  }
  return m1MuTrigScaleFactor;
}

void Analysis_JetMET_Base::DoBtaggingSF(bool sys, const MomKey JetType)
{

  if(Debug()) cout << "DoBtaggingSF: starting! " << sys << " " << JetType << endl;

  static const AnaKey a_Process_CopyEvent_DoBtaggingSF(":Process:CopyEvent:DoBtaggingSF");
  const AnaKey aName_Process_CopyEvent_DoBtaggingSF(Name()+a_Process_CopyEvent_DoBtaggingSF);
  Time()->Start(aName_Process_CopyEvent_DoBtaggingSF);

  static const MomKey mBtagScaleFactor("BtagScaleFactor"), mBtagScaleFactor70("BtagScaleFactor70"),
    mBtagScaleFactor80("BtagScaleFactor80"), mBtagScaleFactor90("BtagScaleFactor90"),
    mBtagScaleFactorLead("BtagScaleFactorLead"), mBtagScaleFactorLead3("BtagScaleFactorLead3"), mMV1("flavor_weight_MV1"), mflavor("flavor_truth_label");

  //This function stolen and modified from BTagCalib.cxx in SUSYTools-00-00-57

  // static vector<MomKey> keys;
  // static vector<float> opvals;
  // static vector<std::string> OPs;
  // static vector<Analysis::CalibrationDataInterfaceROOT*> BTagCalibs;
  // static vector<Analysis::CalibrationDataInterfaceROOT*> JVFBTagCalibs;

  if(Debug()) cout << "DoBtaggingSF: have vectors, should be blank! " << keys.size() << " " << opvals.size() <<  " " << OPs.size() << " " << BTagCalibs.size() << " " << JVFBTagCalibs.size() << endl;
  /*
  if(!myBTagCalib) {
    myBTagCalib = (Analysis::CalibrationDataInterfaceROOT*)fInput->FindObject(TString("myBTagCalib").Append("MV1").Append("0_9827"));
    if(!myBTagCalib) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");
    keys.push_back(mBtagScaleFactor);
    opvals.push_back(0.9827);
    OPs.push_back("0_9827");
    BTagCalibs.push_back(myBTagCalib);

    myJVFBTagCalib = (Analysis::CalibrationDataInterfaceROOT*)fInput->FindObject(TString("myJVFBTagCalib").Append("MV1").Append("0_9827"));
    if(!myJVFBTagCalib) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");

    JVFBTagCalibs.push_back(myJVFBTagCalib);
  }
  */
  if(!myBTagCalib70) {
    myBTagCalib70 = (Analysis::CalibrationDataInterfaceROOT*)InputList()->FindObject(TString("myBTagCalib").Append("MV1").Append("0_7892"));
    if(!myBTagCalib70) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");
    keys.push_back(mBtagScaleFactor70);
    opvals.push_back(0.7892);
    OPs.push_back("0_7892");
    BTagCalibs.push_back(myBTagCalib70);

    myJVFBTagCalib70 = (Analysis::CalibrationDataInterfaceROOT*)InputList()->FindObject(TString("myJVFBTagCalib").Append("MV1").Append("0_7892"));
    if(!myJVFBTagCalib70) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");
    JVFBTagCalibs.push_back(myJVFBTagCalib70);
  }
    if(Debug()) cout << "DoBtaggingSF: have vectors, should be 1! " << keys.size() << " " << opvals.size() <<  " " << OPs.size() << " " << BTagCalibs.size() << " " << JVFBTagCalibs.size() << endl;


  //cout << "myBTagCalib70 " << myBTagCalib70 << endl;
  /*
  if(!myBTagCalib80) {
    myBTagCalib80 = (Analysis::CalibrationDataInterfaceROOT*)fInput->FindObject(TString("myBTagCalib").Append("MV1").Append("0_3511"));
    if(!myBTagCalib80) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");
    keys.push_back(mBtagScaleFactor80);
    opvals.push_back(0.3511);
    OPs.push_back("0_3511");
    BTagCalibs.push_back(myBTagCalib80);

    myJVFBTagCalib80 = (Analysis::CalibrationDataInterfaceROOT*)fInput->FindObject(TString("myJVFBTagCalib").Append("MV1").Append("0_3511"));
    if(!myJVFBTagCalib80) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");
    JVFBTagCalibs.push_back(myJVFBTagCalib80);
  }

  if(!myBTagCalib90) {
    myBTagCalib90 = (Analysis::CalibrationDataInterfaceROOT*)fInput->FindObject(TString("myBTagCalib").Append("MV1").Append("0_0617"));
    if(!myBTagCalib90) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");
    keys.push_back(mBtagScaleFactor90);
    opvals.push_back(0.0617);
    OPs.push_back("0_0617");
    BTagCalibs.push_back(myBTagCalib90);

    myJVFBTagCalib90 = (Analysis::CalibrationDataInterfaceROOT*)fInput->FindObject(TString("myJVFBTagCalib").Append("MV1").Append("0_0617"));
    if(!myJVFBTagCalib90) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");
    JVFBTagCalibs.push_back(myJVFBTagCalib90);
  }

  if(!myBTagCalibLead) { //tight OP, but for 2 leading jets only
    myBTagCalibLead = (Analysis::CalibrationDataInterfaceROOT*)fInput->FindObject(TString("myBTagCalib").Append("MV1").Append("0_9827"));
    if(!myBTagCalibLead) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");
    keys.push_back(mBtagScaleFactorLead);
    opvals.push_back(0.9827);
    OPs.push_back("0_9827");
    BTagCalibs.push_back(myBTagCalibLead);

    myJVFBTagCalibLead = (Analysis::CalibrationDataInterfaceROOT*)fInput->FindObject(TString("myJVFBTagCalib").Append("MV1").Append("0_9827"));
    if(!myJVFBTagCalibLead) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");
    JVFBTagCalibs.push_back(myJVFBTagCalibLead);
  }

  if(!myBTagCalibLead3) { //tight OP, but for 3 leading jets only
    myBTagCalibLead3 = (Analysis::CalibrationDataInterfaceROOT*)fInput->FindObject(TString("myBTagCalib").Append("MV1").Append("0_9827"));
    if(!myBTagCalibLead3) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");
    keys.push_back(mBtagScaleFactorLead3);
    opvals.push_back(0.9827);
    OPs.push_back("0_9827");
    BTagCalibs.push_back(myBTagCalibLead3);

    myJVFBTagCalibLead3 = (Analysis::CalibrationDataInterfaceROOT*)fInput->FindObject(TString("myJVFBTagCalib").Append("MV1").Append("0_9827"));
    if(!myJVFBTagCalibLead3) Abort("EventBuilder_susy1328: ERROR cannot retrieve CalibrationDataInterfaceROOT object from input list!");
    JVFBTagCalibs.push_back(myJVFBTagCalibLead3);
  }
  */
  static Analysis::CalibrationDataVariables* BTagVars = 0;
  if(!BTagVars) {
    BTagVars = new Analysis::CalibrationDataVariables();
        BTagVars->jetAuthor = "AntiKt4TopoLCnoJVF";
  }
  static Analysis::CalibrationDataVariables* BTagVarsJVF = 0;
  if(!BTagVarsJVF) {
    BTagVarsJVF = new Analysis::CalibrationDataVariables();
        BTagVarsJVF->jetAuthor = "AntiKt4TopoLCJVF";
  }

  if(Debug()) cout << "DoBtaggingSF: have vectors! " << keys.size() << " " << opvals.size() <<  " " << OPs.size() << " " << BTagCalibs.size() << " " << JVFBTagCalibs.size() << endl;

  //do subkeys, indices for speed
  static vector<vector<MomKey> > subkeys;
  static vector<BTagIndex> index;
  static vector<BTagIndex> indexJVF;

  if(!subkeys.size()) {
    for(unsigned int j = 0; j<BTagCalibs.size(); ++j) {

      //subkeys
      vector<MomKey> types;
      types.push_back(keys.at(j)+"Nom");
      types.push_back(keys.at(j)+"Bdown");
      types.push_back(keys.at(j)+"Cdown");
      types.push_back(keys.at(j)+"Ldown");
      types.push_back(keys.at(j)+"Bup");
      types.push_back(keys.at(j)+"Cup");
      types.push_back(keys.at(j)+"Lup");
      subkeys.push_back(types);

      if(Debug()) cout << "DoBtaggingSF: have types/subkeys!" << endl;


      //non-JVF indices
      BTagIndex ind;
      bool result = true;
      for (unsigned int i = 0; i!=2; ++i) {
        //0 is efficiency, 1 is SF
        result = result && BTagCalibs.at(j)->retrieveCalibrationIndex("B",OPs.at(j),BTagVars->jetAuthor,i,ind.b[i]);
        result = result && BTagCalibs.at(j)->retrieveCalibrationIndex("C",OPs.at(j),BTagVars->jetAuthor,i,ind.c[i]);
        result = result && BTagCalibs.at(j)->retrieveCalibrationIndex("T",OPs.at(j),BTagVars->jetAuthor,i,ind.t[i]);
        result = result && BTagCalibs.at(j)->retrieveCalibrationIndex("Light",OPs.at(j),BTagVars->jetAuthor,i,ind.l[i]);
      }
      index.push_back(ind);

      if(Debug()) cout << "DoBtaggingSF: have indices!" << endl;

      //JVF indices
      for (unsigned int i = 0; i!=2; ++i) {
        //0 is efficiency, 1 is SF
        result = result && JVFBTagCalibs.at(j)->retrieveCalibrationIndex("B",OPs.at(j),BTagVarsJVF->jetAuthor,i,ind.b[i]);
        result = result && JVFBTagCalibs.at(j)->retrieveCalibrationIndex("C",OPs.at(j),BTagVarsJVF->jetAuthor,i,ind.c[i]);
        result = result && JVFBTagCalibs.at(j)->retrieveCalibrationIndex("T",OPs.at(j),BTagVarsJVF->jetAuthor,i,ind.t[i]);
        result = result && JVFBTagCalibs.at(j)->retrieveCalibrationIndex("Light",OPs.at(j),BTagVarsJVF->jetAuthor,i,ind.l[i]);
      }
      indexJVF.push_back(ind);
      if(!result) Abort("EventBuilder_susy1328: ERROR retrieving b-tagging calibration indices");
      if(Debug()) cout << "DoBtaggingSF: have JVF indices!" << endl;

    }
  }

  bool noSFs = false;
  
  float btagPt = 25.;
  
  //cout << subkeys.size() << " is subkey size" << endl;

    // ***************************************************************************************************************
    // ****                 Function for applying Scale Factors and associated Errors to Analyses                 ****
    // ***************************************************************************************************************
    //
    // - Input Jet information: Pt, Eta, bweight, PDG ID. First three are vectors of floats and last one is vector of ints
    // - Input option for applying shift due to errors. Takes integer value: 0 = Default, 1 = Shift Up, 2 = Shift Down
    // - Function then uses input values to feed into root file, which provides Scale Factors and associated errors
    //   as a function of Pt and Eta. These Scale Factors are then used to calculate the weights applied to individual
    //   jets and the event as a whole.
    // - Function then returns a pair of vector :
    //   first : vector of the weights in the following order :
    // 0 : standard
    // 1 : B-jet efficiency down
    // 2 : C-jet efficiency down
    // 3 : mistag rate efficiency down
    // 4 : B-jet efficiency up
    // 5 : C-jet efficiency up
    // 6 : mistag rate efficiency up
    // The B, C and mistag rate efficiencies have to be added in quadrature at the end of the analysis
    // second : vector of floats containing the weights for the individual jets in the event.
    //
    // - IMPORTANT: SCALE FACTORS NOT DEFINED BELOW JET PT OF 20 GeV OR |ETA| ABOVE 2.5!!!
    //
    // - The efficiency and inefficiency SFs for b and c-jets are valid from 20 to 200 GeV.
    // - The efficiency SFs for light-jets are valid from 20 to 750 GeV.
    // - The inefficiency SFs for light-jets are valid from 20 to 200 GeV
    //   For jets with pT above these limits, values for either 200 GeV or 750 GeV will be used (handled internally now).
    //
    // ***************************************************************************************************************

    float evtWeight = 1.0;
  float evtWeight_B_down = 1.0;
  float evtWeight_B_up = 1.0;
  float evtWeight_C_down = 1.0;
  float evtWeight_C_up = 1.0;
  float evtWeight_L_down = 1.0;
  float evtWeight_L_up = 1.0;

    Analysis::Uncertainty BTagSys = Analysis::Total; //! None, Total, Statistical, Systematic

  //iterate over operating points
  for(unsigned int opiter = 0; opiter<opvals.size(); ++opiter) {

    // BTagCalib objects
      Analysis::CalibrationDataInterfaceROOT* CalibObj = 0;
    Analysis::CalibrationDataVariables* CalibVars = 0;
    BTagIndex* CalibIndex = 0;

      float opval = opvals.at(opiter);
      //cout << opval << " is the opval " << endl;

      evtWeight = 1.0;

      if(sys) {
      evtWeight_B_down = 1.0;
      evtWeight_B_up = 1.0;
      evtWeight_C_down = 1.0;
      evtWeight_C_up = 1.0;
      evtWeight_L_down = 1.0;
      evtWeight_L_up = 1.0;
    }


      //! Perform actions for each individual jet in an event
      int max = jets(JetType);
      if((max>2)&&(keys.at(opiter)==mBtagScaleFactorLead)) max = 2;
      else if((max>3)&&(keys.at(opiter)==mBtagScaleFactorLead3)) max = 3;

      for (int jitr = 0; jitr < max; jitr++) {
      
        if(jet(jitr,JetType).p.Pt()<btagPt) continue;

          if(true) {//if((fabs(jet(jitr,JetType).p.Eta())<2.4) && (jet(jitr,JetType).p.Pt()<50.)) { needs JVF criteria!
            CalibObj = JVFBTagCalibs.at(opiter);
            CalibIndex = &indexJVF.at(opiter);
            CalibVars = BTagVarsJVF;
          }
          else {
            CalibObj = BTagCalibs.at(opiter);
            CalibIndex = &index.at(opiter);
            CalibVars = BTagVars;
          }

      unsigned int* ind = 0;

          //! Set quark flavour to be used
          if (jet(jitr,JetType).Int(mflavor) == 4) ind = CalibIndex->c;
          else if (jet(jitr,JetType).Int(mflavor) == 5) ind = CalibIndex->b;
          else if (jet(jitr,JetType).Int(mflavor) == 15) ind = CalibIndex->t;
          else ind = CalibIndex->l;

          //! Set necessary b-tag variables for retrieving SF + errors
          CalibVars->jetPt  = jet(jitr,JetType).p.Pt()*1000.; //MeV
          CalibVars->jetEta = jet(jitr,JetType).p.Eta();

      Analysis::CalibResult BTagCalibResult;

      if(jet(jitr,JetType).Float(mMV1)>opval) { //b-tagged jet
        BTagCalibResult = CalibObj->getScaleFactor(*CalibVars, ind[1], BTagSys);
      }
        else { //not b-tagged jet
          BTagCalibResult = CalibObj->getInefficiencyScaleFactor(*CalibVars, ind[1], ind[0], BTagSys);
        }

        float jetWeight = noSFs ? 1.0 : (float)BTagCalibResult.first;
        float scale = noSFs ? (float)BTagCalibResult.first : 1.0;

      float jetWeight_down = jetWeight;
      float jetWeight_up = jetWeight;

        evtWeight *= jetWeight;

      jet(jitr,JetType).Set(subkeys.at(opiter).at(0),jetWeight);

      //Treat efficiency and the inefficiency scale factor uncertainties as anti-correlated
      if(sys) {
        if(jet(jitr,JetType).Int(mflavor)==5) {
          if(jet(jitr,JetType).Float(mMV1)>opval) {
            jetWeight_down -= (float)BTagCalibResult.second/scale;
            jetWeight_up += (float)BTagCalibResult.second/scale;
          }
          else {
            jetWeight_down += (float)BTagCalibResult.second/scale;
            jetWeight_up -= (float)BTagCalibResult.second/scale;
          }

          evtWeight_B_down *= jetWeight_down;
          evtWeight_B_up *= jetWeight_up;
          evtWeight_C_down *= jetWeight;
          evtWeight_C_up *= jetWeight;
          evtWeight_L_down *= jetWeight;
          evtWeight_L_up *= jetWeight;

          jet(jitr,JetType).Set(subkeys.at(opiter).at(1),jetWeight_down);
          jet(jitr,JetType).Set(subkeys.at(opiter).at(4),jetWeight_up);
        }
        else if(jet(jitr,JetType).Int(mflavor)==4 || jet(jitr,JetType).Int(mflavor)==15) {
          if(jet(jitr,JetType).Float(mMV1)>opval) {
            jetWeight_down -= (float)BTagCalibResult.second/scale;
            jetWeight_up += (float)BTagCalibResult.second/scale;
          }
          else {
            jetWeight_down += (float)BTagCalibResult.second/scale;
            jetWeight_up -= (float)BTagCalibResult.second/scale;
          }

          evtWeight_B_down *= jetWeight;
          evtWeight_B_up *= jetWeight;
          evtWeight_C_down *= jetWeight_down;
          evtWeight_C_up *= jetWeight_up;
          evtWeight_L_down *= jetWeight;
          evtWeight_L_up *= jetWeight;

          jet(jitr,JetType).Set(subkeys.at(opiter).at(2),jetWeight_down);
          jet(jitr,JetType).Set(subkeys.at(opiter).at(5),jetWeight_up);
        }
        else {
          if(jet(jitr,JetType).Float(mMV1)>opval) {
            jetWeight_down -= (float)BTagCalibResult.second/scale;
            jetWeight_up += (float)BTagCalibResult.second/scale;
          }
          else {
            jetWeight_down += (float)BTagCalibResult.second/scale;
            jetWeight_up -= (float)BTagCalibResult.second/scale;
          }

          evtWeight_B_down *= jetWeight;
          evtWeight_B_up *= jetWeight;
          evtWeight_C_down *= jetWeight;
          evtWeight_C_up *= jetWeight;
          evtWeight_L_down *= jetWeight_down;
          evtWeight_L_up *= jetWeight_up;

          jet(jitr,JetType).Set(subkeys.at(opiter).at(3),jetWeight_down);
          jet(jitr,JetType).Set(subkeys.at(opiter).at(6),jetWeight_up);
        }
      }

      } //end of jet loop

    Set(subkeys.at(opiter).at(0),evtWeight);

    //cout << "btag sf = " <<  subkeys.at(opiter).at(0) << " " << evtWeight << endl;
    if(!sys) continue;

    Set(subkeys.at(opiter).at(1),evtWeight_B_down);
    Set(subkeys.at(opiter).at(2),evtWeight_C_down);
    Set(subkeys.at(opiter).at(3),evtWeight_L_down);
    Set(subkeys.at(opiter).at(4),evtWeight_B_up);
    Set(subkeys.at(opiter).at(5),evtWeight_C_up);
    Set(subkeys.at(opiter).at(6),evtWeight_L_up);

    } //end of OP loop

  Time()->Stop(aName_Process_CopyEvent_DoBtaggingSF);

  if(Debug()){
    cout << "After BTag SF, event is : "<< endl;
    Show();
    cout << opvals.size() << " is the number of operating points" << endl;
  }


  return;
}



void Analysis_JetMET_Base::DoMuonSF(){
  vector<TLorentzVector> sigMuons;
  vector<TLorentzVector> sigElectrons;
  vector<muon_quality> sigMuonsq;
  vector<electron_quality> sigElectronsq; 

  float trigscalefactor = 1.;
  float trigscalefactorup = 1.;
  float trigscalefactordown = 1.;

  int runnumber = -1;
  //set random seed
  if(muons("good")) {
    int seed = (int)(1.e+5*fabs(muon(0,"good").Float("phiUnsmeared")));
      if(!seed) ++seed;
    runnumber = myDataPeriod->getRunNumber_MC12(RunNumber(),seed);
  }

  if(Debug()) cout << " set up init for DoMuonSF " << endl;

  for( int iMu = 0; iMu < muons("good"); iMu++){
    TLorentzVector unsmeared(0.,0.,0.,0.);
    unsmeared.SetPtEtaPhiM( muon(iMu,"good").Float("ptUnsmeared"),
                            muon(iMu,"good").Float("etaUnsmeared"),
                            muon(iMu,"good").Float("phiUnsmeared"),
                            105.6/1000. );
    sigMuons.push_back(muon(iMu,"good").p);

    if(muon(iMu,"good").Bool("isCombinedMuon")){
      sigMuonsq.push_back(combined);
    } else {
      sigMuonsq.push_back(loose);
    }

  }
  if(Debug()) cout << " set up muons for DoMuonSF " << endl;


  if(muons("good")){
    trigscalefactor = MuonTriggerSF->GetTriggerSF(runnumber, true, sigMuons, sigMuonsq, sigElectrons, sigElectronsq, 0).first;
    trigscalefactorup = MuonTriggerSF->GetTriggerSF(runnumber, true, sigMuons, sigMuonsq, sigElectrons, sigElectronsq, 5).first;
    trigscalefactordown = MuonTriggerSF->GetTriggerSF(runnumber, true, sigMuons, sigMuonsq, sigElectrons, sigElectronsq, 6).first;

  }
  const static MomKey m1MuTrigScaleFactor("m1MuTrigScaleFactor");
  const static MomKey m1MuTrigScaleFactorUp("m1MuTrigScaleFactorUp");
  const static MomKey m1MuTrigScaleFactorDown("m1MuTrigScaleFactorDown");

  if(Debug()) cout << " did muon trig scale factor for DoMuonSF , val " << trigscalefactor << endl;


  Set(m1MuTrigScaleFactor, trigscalefactor);
  Set(m1MuTrigScaleFactorUp, trigscalefactorup);
  Set(m1MuTrigScaleFactorDown, trigscalefactordown);

  float sf = 1.;
  float sf_up = 1.;
  float sf_down = 1.;

  if(Debug()) cout << " going to loop over signal muons for for DoMuonSF " << endl;

  for(int iMu = 0; iMu < muons("good"); iMu++){
    float charge = muon(iMu,"good").Float("charge");
    TLorentzVector muVec = muon(iMu,"good").p;
    if(Debug()) muon(iMu,"good").Show();

    float val = MuonSF->scaleFactor(charge,muVec);

    sf *= val;

    double stat = MuonSF->scaleFactorUncertainty(charge,muVec);
    double syst = MuonSF->scaleFactorSystematicUncertainty(charge,muVec);
      /// The systematic error has to be added linearly to the statistical error because it describes how the efficiency scale factors is shifted after a modified tag-and-probe selection.
    float sf_unc = stat + syst;

    sf_up   *= (val + sf_unc);
    sf_down *= (val - sf_unc);

    if(Debug()) cout << "In MuonSF, Actually did a MuonSF!" << endl;
  }

  if(Debug()) cout << " did muon scale factor for DoMuonSF " << endl;


  const static MomKey m1MuScaleFactor("m1MuScaleFactor");
  const static MomKey m1MuScaleFactorUp("m1MuScaleFactorUp");
  const static MomKey m1MuScaleFactorDown("m1MuScaleFactorDown");
  Set(m1MuScaleFactor,sf);
  Set(m1MuScaleFactorUp,sf_up);
  Set(m1MuScaleFactorDown,sf_down);

  if(Debug()) cout << " did everything for DoMuonSF " << endl;
  
}

void Analysis_JetMET_Base::MakeDRMatch(const MomKey JetKey1, const MomKey JetKey2){
  if(Debug()) cout << "making DR Match between " << JetKey1 << " and " << JetKey2 << endl;
  static const MomKey DRMatch("DRMatch");
  for(int iJet1 = 0; iJet1 < jets(JetKey1); iJet1++){

    jet(iJet1, JetKey1).AddVec(DRMatch+JetKey2);
    for(int iJet2 = 0; iJet2 < jets(JetKey2); iJet2++){

      if(jet(iJet1, JetKey1).p.DeltaR(jet(iJet2, JetKey2).p) < 0.3){
        if(!jet(iJet2,JetKey2).Exists(DRMatch+JetKey1)){
          jet(iJet1, JetKey1).Add(DRMatch+JetKey2, &jet(iJet2, JetKey2));
          jet(iJet2, JetKey2).AddVec(DRMatch+JetKey1, true);
          jet(iJet2, JetKey2).Add(DRMatch+JetKey1, &jet(iJet1, JetKey1));
        }
      }

    }

  }

  for(int iJet2 = 0; iJet2 < jets(JetKey2); iJet2++){
    if(!jet(iJet2,JetKey2).Exists(DRMatch+JetKey1)){
      jet(iJet2,JetKey2).AddVec(DRMatch+JetKey1, true);
    }
  }
  if(Debug()) cout << "finished making DR Match between " << JetKey1 << " and " << JetKey2 << endl;

}

void Analysis_JetMET_Base::AddGhostMatch(const MomKey JetKey, const MomKey MatchObject){

  const TString JetTypeS(JetKey.Data());
  MomKey baseCons = MomKey("clustersLCTopo");
  fastjet::JetAlgorithm baseAlgo = fastjet::antikt_algorithm;
  float baseR = 0.4;

  if(JetTypeS.Contains("7"))
    baseR = 0.7;
  else if(JetTypeS.Contains("6"))
    baseR = 0.6;
  else if(!JetTypeS.Contains("4")){
    cout << "YOUR TYPE OF JET IS NOT BEING PARSED PROPERLY, DEFAULTING TO 0.4" << endl;
  }


  if(!JetTypeS.Contains("AntiKt") && JetTypeS.Contains("CamKt"))
    baseAlgo = fastjet::cambridge_algorithm;
  else if(!JetTypeS.Contains("AntiKt") && JetTypeS.Contains("Kt"))
    baseAlgo = fastjet::kt_algorithm;
  else if(!JetTypeS.Contains("AntiKt"))
    cout << "YOUR TYPE OF JET IS NOT BEING PARSED PROPERLY, DEFAULTING TO ANTIKT" << endl;
  
  if(JetTypeS.Contains("Truth"))
    baseCons = MomKey("truthsStable");
  else if(JetTypeS.Contains("Track"))
    baseCons = MomKey("tracksgood");
  else if(!JetTypeS.Contains("LCTopo"))
    cout << "YOUR TYPE OF JET IS NOT BEING PARSED PROPERLY, DEFAULTING TO LCTOPO" << endl;


  AddGhostMatch(JetKey, MatchObject, baseCons, baseAlgo, baseR);

}


bool Analysis_JetMET_Base::AddDRMatch(const MomKey JetKey, const MomKey MatchObject){
  float baseR = 0.4;

  TString JetTypeS = TString(JetKey.Data());

  if(JetTypeS.Contains("7"))
    baseR = 0.7;
  else if(JetTypeS.Contains("6"))
    baseR = 0.6;
  else if(!JetTypeS.Contains("4")){
    cout << "YOUR TYPE OF JET IS NOT BEING PARSED PROPERLY, DEFAULTING TO 0.4" << endl;
  }

  const static MomKey DR("DR");

  for(int iJet = 0; iJet < jets(JetKey); iJet++){
    jet(iJet,JetKey).AddVec(MatchObject+DR);
    for(int iTr = 0; iTr < Objs(MatchObject); iTr++){
      if(((Particle*) Obj(MatchObject,iTr))->p.DeltaR(jet(iJet,JetKey).p) < baseR){
        jet(iJet,JetKey).Add(MatchObject+DR,Obj(MatchObject,iTr));
      }
    }
  }
  return true;
}

void Analysis_JetMET_Base::AddGhostMatch(const MomKey JetKey, const MomKey MatchObject, const MomKey ConstType, const fastjet::JetAlgorithm AlgoIn, const float Radius){

  if(Debug()) cout << "gonna match you some ghosts" << endl;
  // make a vector to be clustered
  vector<fastjet::PseudoJet> AllTheThings = ObjsToPJ(ConstType);


  bool doTrack = false;

  if(TString(MatchObject.Data()).Contains("track"))
    doTrack = true;
  AddAsGhosts(MatchObject, &AllTheThings, doTrack);

  fastjet::JetDefinition jetDef(AlgoIn, Radius, fastjet::E_scheme, fastjet::Best);
  fastjet::ClusterSequence clustSeq(AllTheThings, jetDef);  

  vector<fastjet::PseudoJet> inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets());

  const static MomKey constscale_pt("constscale_pt");
  const static MomKey Ghost("Ghost");
  const static MomKey jetl("jets");
  Sort(jetl+JetKey, constscale_pt);

  if(Debug()) cout << " in ghost matching, sorted " << endl;

  for(int iJet = 0 ; iJet < jets(JetKey) ; iJet++){
    jet(iJet,JetKey).AddVec(MatchObject+Ghost);

    if(iJet > inclusiveJets.size() - 1)
      break;
    if(inclusiveJets.size() ==0)  // extra protection for case jJet == 0 and inclusiveJets.size() ==0 
      break;                      // where above statement fails due to unsigned int comparison with negative number....
    fastjet::PseudoJet jetFJ = inclusiveJets[iJet];

   

    vector<fastjet::PseudoJet> constituents = jetFJ.constituents();

     //cout << " sanity check time! pt of fj jet = " << jetFJ.perp() << " pt of pa jet " << jet(iJet,JetKey).Float(constscale_pt) << endl;

    for(unsigned int iConst  = 0; iConst < constituents.size(); iConst++){
      if(doTrack){
        if(constituents[iConst].user_info<PJ_Info>().GhostTrack){
          jet(iJet,JetKey).Add(MatchObject+Ghost, constituents[iConst].user_info<PJ_Info>().Pointer);
        }
      } else {
        if(constituents[iConst].user_info<PJ_Info>().GhostParticle){
          jet(iJet,JetKey).Add(MatchObject+Ghost, constituents[iConst].user_info<PJ_Info>().Pointer);
        }
      }
    }
    

  }// end loop over jets

  SortPt(jetl+JetKey);



}

void Analysis_JetMET_Base::AddAsGhosts(const MomKey GhostObject, vector<fastjet::PseudoJet>* TheThings, bool IsTrack){

  for(int iGhost = 0; iGhost < Objs(GhostObject); iGhost++){
    Particle* part = (Particle*) Obj(GhostObject,iGhost);
    TLorentzVector* p = &(part->p);
    if(p->E() > 0){
      //fastjet::PseudoJet clus(p->Px()/10000000., p->Py()/1000000., p->Pz()/1000000., p->E()/1000000.);
      fastjet::PseudoJet clus(0., 0., 0., 0.);
      float wrappedPhi = p->Phi();
      const float pi = 3.14159265358979323846;
      if(wrappedPhi > pi){
        wrappedPhi = -1.* (wrappedPhi - pi);
      }
      clus.reset_PtYPhiM(p->Perp()/1000000., p->Rapidity(), wrappedPhi, p->M());

      //cout << " new eta = " << clus.eta() << " old eta is " << part->p.Eta() << endl;
      //cout << " new phi = " << clus.phi_std() << " old phi is " << part->p.Phi() << endl;

      if(IsTrack)
        clus.set_user_info(new PJ_Info(true, false, part));
      else
        clus.set_user_info(new PJ_Info(false, true, part));

      TheThings->push_back(clus);
    } // energy check
  }// ghost loop

}


// specific function for truth jets
// add a new list, constituentsCharged
// this willl contain only charged constitutents of truth jets
void Analysis_JetMET_Base::AddTruthStableCharged(const MomKey JetKey){

  if(TString(JetKey.Data()).Contains("Truth")==0){
    cout << " ONLY USE THIS FUNCTION ON TRUTH JETS" << endl;
    return;
  }

  static const MomKey mCC("constituentsCharged");
  static const MomKey mC("constituents");
  static const MomKey mCharge("charge");

  for(int iJet = 0; iJet < jets(JetKey); iJet++){
    if(!jet(iJet,JetKey).Exists(mC)){
      cout << "NEED CONSTITUENTS OF TRUTH JETS SET UP" << endl;
      return;
    }
    jet(iJet, JetKey).AddVec(mCC);
    for(int iC = 0; iC < jet(iJet,JetKey).Objs(mC); iC++){
      Particle* theConst = (Particle*) jet(iJet, JetKey).Obj(mC, iC);
      if(fabs(theConst->Float(mCharge)) < 1E-13){
        continue;
      }
      jet(iJet,JetKey).Add(mCC, theConst);
    }
  }

}


void Analysis_JetMET_Base::CalculateResponseOffset(const MomKey JetKey, const MomKey MatchKey){
  static const MomKey Offset("Offset");
  static const MomKey Response("Response");
  static const MomKey TruthPT("TruthPT");
  if(Debug()) cout << " making Response/Offset for " << JetKey << " and " << MatchKey << endl;

  for(int iJet = 0; iJet < jets(JetKey); iJet++){

    if(jet(iJet, JetKey).Objs(MatchKey) > 0){
      Particle* truth = (Particle*) jet(iJet, JetKey).Obj(MatchKey, 0);
      Particle* reco  = &jet(iJet, JetKey);
      jet(iJet, JetKey).Set(Offset + MatchKey, reco->p.Perp() - truth->p.Perp());
      jet(iJet, JetKey).Set(Response + MatchKey, reco->p.Perp() / truth->p.Perp());
    } else {
      jet(iJet,JetKey).Set(Offset + MatchKey, -999999);
      jet(iJet,JetKey).Set(Response + MatchKey, -999999);
    }

  }
  if(Debug()) cout << "finished making Response/Offset for " << JetKey << " and " << MatchKey << endl;

}

void Analysis_JetMET_Base::MakeDRResponseOffset(const MomKey JetKey, const MomKey TruthJetKey){
  MakeDRMatch(JetKey, TruthJetKey);
  const static MomKey DRMatch("DRMatch");
  CalculateResponseOffset(JetKey, DRMatch+TruthJetKey);
}

void Analysis_JetMET_Base::WorkerBegin(){

}

void Analysis_JetMET_Base::RemoveJetsPtCut(const MomKey JetType, float PtCut){
  const MomKey FullJetType = MomKey("jets")+JetType;
  for(int iJet = 0; iJet < jets(JetType); iJet++){
    if(jet(iJet,JetType).p.Perp() < PtCut){
      Remove(FullJetType, iJet);
      iJet--;
    } // end check on jet pt
  } // end loop over jets
  return;
}

MomKey Analysis_JetMET_Base::AddSubjets(const MomKey JetType, const fastjet::JetAlgorithm algo, const double fcut, const double ptcut){
  const static MomKey ConsKey("constituents");
  const static MomKey kt("KT");
  const static MomKey ca("CA");
  MomKey NSubKey("NSub");

  switch(algo){
    case fastjet::kt_algorithm:
      NSubKey += kt;
      break;
    case fastjet::cambridge_algorithm:
      NSubKey += ca;
      break;
    default:
      cout << " this jet algorithm is not supported! quitting " << endl;
      exit(-1);
      break;
  }

  NSubKey += TString::Format("FCut%.0fPtCut%.0f",100.*fcut,ptcut);

  fastjet::JetDefinition jetDef(algo,1.0);
  for(int iJet = 0; iJet < jets(JetType); iJet++){

    vector<fastjet::PseudoJet> jetConst = ObjsToPJ(ConsKey, &jet(iJet,JetType)); 
    fastjet::ClusterSequence subCluster(jetConst,jetDef);

    vector<fastjet::PseudoJet> subs = subCluster.exclusive_jets(TMath::Power(fcut * jet(iJet,JetType).p.Perp(), 2));
    int nSub = 0;
    for(unsigned int iSub = 0; iSub < subs.size(); iSub++){
      if(subs[iSub].pt() > ptcut){
        nSub++;
      }
    }

    if(Debug()) cout << " jet " << iJet << " type " << JetType << " started with " << subs.size() << " final number " << nSub <<  " dcut is " << TMath::Power(fcut * jet(iJet,JetType).p.Perp(),2) << endl;
    jet(iJet, JetType).Set(NSubKey, nSub);
  }

  if(Debug()) cout << "subjets " << NSubKey << " are set up now " << endl;

  return NSubKey;

}

MomKey Analysis_JetMET_Base::SetEventSubjets(const MomKey JetType, const MomKey SubjetType, int NJets){
  int nsub = 0;

  for(int iJet = 0; iJet < jets(JetType) && iJet < NJets; iJet++){
    Particle* j = &jet(iJet,JetType);
    if(!j->Exists(SubjetType)){
      cout << " Set up nsubjets for this jet type, " << JetType << ", first! " << endl;
      return MomKey("");
    }
    nsub += j->Int(SubjetType);
  }

  MomKey NSubKey = MomKey("NSubjet") + JetType + SubjetType + TString::Format("NJets%i",NJets);

  Set(NSubKey, nsub);

  if(Debug()) cout << " finished setting up " << NSubKey << endl;

  return NSubKey;
}


MomKey Analysis_JetMET_Base::SetEventSubjettiness(const MomKey JetType, const MomKey nsubType, const int NJets){
  float nsubprod = 1.;

  for(int iJet = 0; iJet < jets(JetType) && iJet < NJets; iJet++){
    Particle* j = &jet(iJet,JetType);
    if(!j->Exists(nsubType)){
      cout << " Set up nsubjettiness for this jet type, " << JetType << ", first! " << endl;
      return MomKey("");
    }
    nsubprod *= j->Float(nsubType);
  }

  double ensub = TMath::Power(nsubprod, 1./((float) NJets));

  MomKey ESubKey = MomKey("EventSubjettiness") + JetType + nsubType + TString::Format("NJets%i",NJets);

  Set(ESubKey, nsubprod);

  if(Debug()) cout << " finished setting up " << ESubKey << endl;

  return ESubKey;
}

MomKey Analysis_JetMET_Base::SetEventMass(const MomKey JetType, int NJets){
  float massTotal = 0.;

  for(int iJet = 0; iJet < jets(JetType) && iJet < NJets; iJet++){
    massTotal += jet(iJet,JetType).p.M();
  }

  MomKey MassKey = MomKey("TotalJetMass")+JetType+TString::Format("NJets%i",NJets);

  Set(MassKey, massTotal);

  if(Debug()) cout << " finished setting up " << MassKey << endl;

  return MassKey;
}

MomKey Analysis_JetMET_Base::SetHt(const MomKey JetType, int NJets){

  float ht = 0.;

  int lim = ((NJets > 0) && (jets(JetType) > NJets)) ? NJets : jets(JetType);

  MomKey HtKey = MomKey("Ht") + JetType;

  if(NJets > 0){
    HtKey +=  TString::Format("NJets%i",NJets);
  }

  for(int iJet = 0; iJet < lim; iJet++){
    ht += jet(iJet,JetType).p.Perp(); 
  }

  Set(HtKey, ht);

  if(Debug()) cout << " finished setting up " << HtKey << endl;

  return HtKey;
}

void Analysis_JetMET_Base::SetFatJetOverlap(const MomKey FatKey, const MomKey SmallKey, float dR){
	const MomKey OverlapKey = "Overlap_" + SmallKey;

	for(int iFat = 0; iFat < jets(FatKey); iFat++){
		int NOverlap = 0;
		for(int iSmall = 0; iSmall < jets(SmallKey); iSmall++){
			if(jet(iFat, FatKey).p.DeltaR(jet(iSmall, SmallKey).p) < dR){
				NOverlap++;
			}// end if DeltaR
		} // end for iSmall
		jet(iFat, FatKey).Set(OverlapKey, NOverlap);
	} // end for iFat
	return;
}


///=========================================
/// Generic function to return a vector of pseudojets which are the objects to calculate
/// Some fast-jet related properties. If the last argument is not supplied,
/// Will give full event level clusters/tracks/etc. Otherwise for a specified jet
///=========================================
vector<fastjet::PseudoJet> Analysis_JetMET_Base::ObjsToPJ(const MomKey type, MomentObj* jet){
  jet = (jet!=0)?jet:this;
  vector<fastjet::PseudoJet> out;
  for(int i = 0; i < jet->Objs(type); i++){
    Particle* part = (Particle*) jet->Obj(type,i);
    TLorentzVector* p = &(part->p);
    if(p->E() > 0){
      fastjet::PseudoJet clus(p->Px(), p->Py(), p->Pz(), p->E());
      // label this as not a ghost track, not a ghost particle, and add the pointer to the particle
      // note that if you use this method to get tracks (or truth), they'll be labeled as 'clusters', but that's fine
      clus.set_user_info(new PJ_Info(false, false, part));
      out.push_back(clus);
    }
  }
  return out;
}

MomKey Analysis_JetMET_Base::MakePrunedJetsName(const MomKey JetType, const fastjet::JetAlgorithm algo, const double zcut, const double dcut){

 const static MomKey SJetKey("jets");
  TString AlgoKeyS;
  switch(algo){
    case fastjet::kt_algorithm:
      AlgoKeyS = "Kt";
      break;

    case fastjet::cambridge_algorithm:
      AlgoKeyS = "CamKt";
      break;

    default:
      cout << " this jet algorithm is not supported! quitting " << endl;
      exit(-1);
      break;
  }
  MomKey ParamKey;
  if(dcut < 0.){
    ParamKey = AlgoKeyS + TString::Format("ZCut%.0fDCutMPT", zcut*100.);
  } else {
    ParamKey = AlgoKeyS + TString::Format("ZCut%.0fDCut%.0f", zcut*100., dcut*100.);
  } 
  const static MomKey PrunedKeyLiteral("Pruned");
  MomKey PrunedKey = JetType+PrunedKeyLiteral+ParamKey;
  return PrunedKey;
}


//add pruned jets to an existing collection. need to have constituents of that jet.
MomKey Analysis_JetMET_Base::MakePrunedJets(const MomKey JetType, const fastjet::JetAlgorithm algo, const double zcut, const double dcut){
  
  if(Debug()) cout << "going to make some pruned jets" << endl;
  const static MomKey SJetKey("jets");
  TString AlgoKeyS;
  switch(algo){
  	case fastjet::kt_algorithm:
  	  AlgoKeyS = "Kt";
  	  break;

  	case fastjet::cambridge_algorithm:
  	  AlgoKeyS = "CamKt";
  	  break;

    default:
      cout << " this jet algorithm is not supported! quitting " << endl;
      exit(-1);
      break;
  }
  MomKey ParamKey;
  if(dcut < 0.){
    ParamKey = AlgoKeyS + TString::Format("ZCut%.0fDCutMPT", zcut*100.);
  } else {
    ParamKey = AlgoKeyS + TString::Format("ZCut%.0fDCut%.0f", zcut*100., dcut*100.);
  } 
  const static MomKey PrunedKeyLiteral("Pruned");
  MomKey PrunedKey = JetType+PrunedKeyLiteral+ParamKey;
  if(Debug()) cout << "PrunedKey = " << PrunedKey << endl;
  MomKey FPrunedKey = SJetKey + PrunedKey;

  if(Debug()) cout << "constructing pruned jets " << FPrunedKey << endl;

  AddVec(FPrunedKey);


  if(Debug()) cout << "pruner ready" << endl;

  const TString JetTypeS = JetType.Data();
  double baseR;
  fastjet::JetAlgorithm baseAlgo;
  MomKey baseCons;

  if(JetTypeS.Contains("7"))
  	baseR = 0.7;
  else if(JetTypeS.Contains("4"))
    baseR = 0.4;
  if(JetTypeS.Contains("AntiKt"))
  	baseAlgo = fastjet::antikt_algorithm;
  if(JetTypeS.Contains("LCTopo"))
  	baseCons = MomKey("clustersLCTopo");
  if(JetTypeS.Contains("Truth"))
    baseCons = MomKey("truthsStable");
  if(JetTypeS.Contains("Track"))
    baseCons = MomKey("tracksgood");

  vector<fastjet::PseudoJet> inputConst = ObjsToPJ(baseCons);
  fastjet::JetDefinition jetDef(baseAlgo, baseR,fastjet::E_scheme, fastjet::Best);
  fastjet::ClusterSequence clustSeq(inputConst, jetDef);  

  vector<fastjet::PseudoJet> normalJets = sorted_by_pt(clustSeq.inclusive_jets(0.));


  for(unsigned int iJet = 0; iJet < normalJets.size() && iJet < jets(JetType); iJet++){

    // if dcut is negative, use the m/pt for dcut
    // otherwise choose the explicit dcut provided
    fastjet::Pruner* pruner;
    if(dcut < 0.){
      pruner = new fastjet::Pruner(algo, zcut, normalJets[iJet].m() / normalJets[iJet].pt());
    } else{
      pruner = new fastjet::Pruner(algo, zcut, dcut);
    }

  	fastjet::PseudoJet theJet = normalJets[iJet];
  	fastjet::PseudoJet prunedJetP = (*pruner)(theJet);
  	Particle* prunedJet = new Particle();
  	prunedJet->p.SetPtEtaPhiE(prunedJetP.pt(), prunedJetP.eta(), prunedJetP.phi(), prunedJetP.e());
  	const static MomKey ParentKey("parent");
  	prunedJet->AddVec(ParentKey);
  	prunedJet->Add(ParentKey, &jet(iJet, JetType));

  	Add(FPrunedKey, prunedJet);

  	jet(iJet, JetType).AddVec(FPrunedKey, true); // true for circular link prevention! gc
  	jet(iJet, JetType).Add(FPrunedKey, prunedJet);

  	vector<fastjet::PseudoJet> constituentsP = prunedJetP.constituents();
  	static const MomKey ConsKey("constituents");
  	prunedJet->AddVec(ConsKey);
  	for(unsigned int iCons = 0; iCons < constituentsP.size(); iCons++){
  		const PJ_Info* info = &(constituentsP[iCons].user_info<PJ_Info>());
  		prunedJet->Add(ConsKey, info->Pointer);
  	} // end loop over cons

    delete pruner;
  } // end loop over normal jets

  if(Debug()) cout << "pruned jets " << FPrunedKey << " finished!" << endl;
  if(Debug()) cout << "there are " << jets(PrunedKey) << " of them " << endl;
  return PrunedKey;
}

//make normal jets, no trimming, etc.
MomKey Analysis_JetMET_Base::MakeJets(const fastjet::JetAlgorithm algo, const double jetR, const MomKey constType, const MomKey extra, const double minpt){

  const static MomKey SJetKey("jets");
  vector<fastjet::PseudoJet> inputConst = ObjsToPJ(constType);

  fastjet::JetDefinition jetDef(algo, jetR,fastjet::E_scheme, fastjet::Best);
  fastjet::ClusterSequence clustSeq(inputConst, jetDef);  

  vector<fastjet::PseudoJet> inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(minpt));

  TString key;

  switch(algo){
    case fastjet::antikt_algorithm:
      key = "AntiKt";
      break;
    
    case fastjet::kt_algorithm:
      key = "Kt";
      break;
    
    case fastjet::cambridge_algorithm:
      key = "CamKt";
      break;
    
    default:
      cout << " this jet algorithm is not supported! quitting " << endl;
      exit(-1);
      break;
  }
  
  key+= TString::Format("%.0f",10.*jetR);

  const static MomKey LCKey("clustersLCTopo");
  const static MomKey TrackKey("tracksgood");
  const static MomKey TruthKey("truthsStable");

  if(constType==LCKey){
  	key+="LCTopo";
  } else if (constType==TrackKey){
  	key+="TrackZ";
  } else if(constType==TruthKey){
    key+="Truth";
  }
  key+=extra;

  if(Debug()) cout << "MakeJets with key " << key << endl;

  MomKey FinalKey(key);
  MomKey FFinalKey = SJetKey + FinalKey;

  AddVec(SJetKey+FinalKey);

  for(unsigned int iJet = 0 ; iJet < inclusiveJets.size() ; iJet++){
  	fastjet::PseudoJet jet = inclusiveJets[iJet];
  	vector<fastjet::PseudoJet> constituents = jet.constituents();
  	Particle* jetP = new Particle();
  	jetP->p.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.e());
  	static const MomKey ConsKey("constituents");
  	jetP->AddVec(ConsKey);
  	for(unsigned int iCons = 0; iCons < constituents.size(); iCons++){
  		const PJ_Info* info = &(constituents[iCons].user_info<PJ_Info>());
  		jetP->Add(ConsKey, info->Pointer);
  	} // end loop over cons
  	Add(FFinalKey, jetP);
  }// end loop over jets

  return FinalKey;
}

///=========================================
/// Adds a calculation for the minimized n-subjettiness moments (up to \tau_6)
/// Uses the plugin version 0.5.1 from jthaler.net
/// Plugin slightly modified; located in utils of JetMETCommon
///=========================================
void Analysis_JetMET_Base::AddNsub(const MomKey constType, const MomKey jetCollection, const bool min, const float beta){
  if(Debug()) cout << "In AddMinNsub with constType " << constType << " and jetCollection " << jetCollection << endl;

  const TString jetCollectionS = jetCollection.Data();

  double jetRad = 1.0;

  if(jetCollectionS.Contains("10"))
    jetRad=1.0;
  else if(jetCollectionS.Contains("12"))
    jetRad=1.2;
  else if(jetCollectionS.Contains("4"))
    jetRad=0.4; 
  else if(jetCollectionS.Contains("6"))
    jetRad=0.6; 
  else if(jetCollectionS.Contains("15"))
    jetRad=1.5;
  else if(jetCollectionS.Contains("7"))
    jetRad=0.7;

  for(int iJet = 0; iJet < jets(jetCollection); iJet++){
    Particle* myjet = &jet(iJet, jetCollection);

    if(Debug()){
      cout << "adding nsubjettiness for the jet " << iJet << " of jet collection " << jetCollectionS << endl;
      myjet->Show();
    }
    
    CalcNsub(myjet, constType, min, beta, jetRad);

  }
  return;
}

void Analysis_JetMET_Base::CalcNsub(Particle* myjet, const MomKey constType, const bool min, const float beta, const float jetRad){
	//if(Debug()) myjet->Show();
    vector<fastjet::PseudoJet> particles = ObjsToPJ(constType, myjet);

    double tau_1, tau_2, tau_3, tau_4, tau_5, tau_6;
    NsubParameters myNsubParam(beta, jetRad);

    Njettiness* NJetCalc;
    if(min){
      NJetCalc = new Njettiness(Njettiness::onepass_kt_axes, myNsubParam);
    }
    else{
      NJetCalc = new Njettiness(Njettiness::kt_axes, myNsubParam);
    }

    tau_1 = NJetCalc->getTau(1, particles);
    tau_2 = NJetCalc->getTau(2, particles);
    tau_3 = NJetCalc->getTau(3, particles);
    tau_4 = NJetCalc->getTau(4, particles);
    tau_5 = NJetCalc->getTau(5, particles);
    tau_6 = NJetCalc->getTau(6, particles);

    if(Debug()) cout << "In AddMinNsub, did calculations" << endl;

    delete NJetCalc;

    MomKey extra = "";
    if(min) extra += "Min";
    if(beta > 1.00001 || beta < .999999) extra += TString::Format("Beta%.1f",beta);

    myjet->Set(constType+"tau1"+extra, tau_1);
    myjet->Set(constType+"tau2"+extra, tau_2);
    myjet->Set(constType+"tau3"+extra, tau_3);
    myjet->Set(constType+"tau4"+extra, tau_4);
    myjet->Set(constType+"tau5"+extra, tau_5);
    myjet->Set(constType+"tau6"+extra, tau_6);

    myjet->Set(constType+"tau32"+extra, tau_3/tau_2);
    myjet->Set(constType+"tau31"+extra, tau_3/tau_1);
    myjet->Set(constType+"tau21"+extra, tau_2/tau_1);

    myjet->Set(constType+"tau63"+extra, tau_6/tau_3);
    myjet->Set(constType+"tau43"+extra, tau_4/tau_3);
    myjet->Set(constType+"tau42"+extra, tau_4/tau_2);
    myjet->Set(constType+"tau41"+extra, tau_4/tau_1);
    myjet->Set(constType+"tau61"+extra, tau_6/tau_1);

    if(Debug()) cout << "In CalcNsub, set all values!" << endl;
    return;
}

////////////////////////////////////////////////////////
///AddStableParticles: simple method to copy stable particles somewhere known
////////////////////////////////////////////////////////
MomKey Analysis_JetMET_Base::AddStableParticles(){
  const static MomKey tstable("truthsStable");
  const static MomKey statkey("status");
  const static MomKey pdgId("pdgId"), barcode("barcode"), vx_barcode("vx_barcode");
  AddVec(tstable);
  for(int iTruth = 0; iTruth < truths(); iTruth++){
    int status = truth(iTruth).Int(statkey);
    int bar    = truth(iTruth).Int(barcode);
    int vx_bar = truth(iTruth).Int(vx_barcode);
    int pdg = fabs(truth(iTruth).Int(pdgId));
    if((((status % 1000 == 1) ||
        (status % 1000 == 2 && status > 1000) ||
        (status==2 && vx_bar<-200000)) && (bar<200000)) &&
       !(pdg==12 || pdg==14 || pdg==16 || pdg==13 ||
        (pdg==1000022 &&  status%1000==1 ) ||
        (pdg==5100022 &&  status%1000==1 ) ||
        (pdg==1000024 &&  status%1000==1 ) ||
        (pdg==39      &&  status%1000==1 ) ||
        (pdg==1000039 &&  status%1000==1 ) ||
        (pdg==5000039 &&  status%1000==1 ))){ // then we have a stable particle
          // sort of a misnomer: we actually want stable without nu/mu... might clean this later
          Add(tstable, &truth(iTruth));
          //cout << "adding particle with pdgId " << pdg << " and pt " << truth(iTruth).p.Pt() << " and eta " << truth(iTruth).p.Eta() <<  endl;
    } // end if
  } // end for
  return tstable;
}// end function

////////////////////////////////////////////////////////
///AddStableParticles: simple method to add the good tracks
////////////////////////////////////////////////////////
MomKey Analysis_JetMET_Base::AddGoodTracks(){
  const static MomKey tracksgood("tracksgood");
  AddVec(tracksgood);
  for(int iTr = 0; iTr < tracks(); iTr++){
    if(fabs(track(iTr).Float("z0_wrtPV") * sin(track(iTr).Float("theta")) ) > 1.5) continue;
    if(fabs(track(iTr).Float("d0_wrtPV")) > 1.) continue;
    if(track(iTr).p.Pt() < 0.5) continue;
    if(fabs(track(iTr).p.Eta()) > 2.5) continue;
    if(track(iTr).Int("nPixHits")+track(iTr).Int("nPixelDeadSensors") < 1) continue;
    if(track(iTr).Int("nSCTHits")+track(iTr).Int("nSCTDeadSensors") < 6) continue;
    if(track(iTr).Float("chi2")/((float)track(iTr).Int("ndof")) > 3.) continue;
    Add(tracksgood,&track(iTr));
  }

  return tracksgood;
}


TH1* Analysis_JetMET_Base::Fill(AnaKey name, double XVal, double weight, int NumBinsX, double XMin, double XMax){
  
  	static TClass* th1cl = TH1::Class();
	
	AnaKey histoname(name);

	TObject* obj = OutputDir()->Get(histoname);
	if(!obj) {
		TDirectory* bkp = gDirectory;
		OutputDir()->cd();
    	TH1D* th1 = new TH1D(histoname.Data(), histoname.Data(), NumBinsX, XMin, XMax);
    	gDirectory = bkp;
      cout << "Booked new histo " << histoname.Data() << endl;	
    	
		th1->StatOverflows(1);
    	th1->Sumw2();
		th1->Fill(XVal,weight);
		return th1;	
	}
	else if(obj->InheritsFrom(th1cl)) {
		TH1* th1 = static_cast<TH1*>(obj);
		th1->Fill(XVal,weight);
		return th1;
	}

	TString msg("Analysis_JetMET_Base: ERROR No valid object ");
	msg.Append(name);
	msg.Append(" found");
	Abort(msg);		

	return 0;
}


TH1* Analysis_JetMET_Base::Fill(AnaKey name, double XVal, double weight, int NumBinsX, const float* Xbins){

  	static TClass* th1cl = TH1::Class();
	
	AnaKey histoname(name);

	TObject* obj = OutputDir()->Get(histoname);
	if(!obj) {
		TDirectory* bkp = gDirectory;
		OutputDir()->cd();
    	TH1D* th1 = new TH1D(histoname.Data(), histoname.Data(), NumBinsX, Xbins);
    	gDirectory = bkp;
      cout << "Booked new histo " << histoname.Data() << endl;	
    	
		th1->StatOverflows(1);
    	th1->Sumw2();
		th1->Fill(XVal,weight);
		return th1;	
	}
	else if(obj->InheritsFrom(th1cl)) {
		TH1* th1 = static_cast<TH1*>(obj);
		th1->Fill(XVal,weight);
		return th1;
	}

	TString msg("Analysis_JetMET_Base: ERROR No valid object ");
	msg.Append(name);
	msg.Append(" found");
	Abort(msg);		

	return 0;
}



TH1* Analysis_JetMET_Base::Fill(AnaKey name, std::string label, double weight, std::vector<std::string> x_labels){

  	static TClass* th1cl = TH1::Class();
	
	AnaKey histoname(name);

	TObject* obj = OutputDir()->Get(histoname);
	if(!obj) {
		TDirectory* bkp = gDirectory;
		OutputDir()->cd();
    	TH1D* th1 = new TH1D(histoname.Data(), histoname.Data(), x_labels.size(), 0.5, x_labels.size()+0.5);
    	gDirectory = bkp;
    	
      cout << "Booked new histo " << histoname.Data() << endl;	
		for(unsigned int i=0; i<x_labels.size(); i++)
		th1->GetXaxis()->SetBinLabel(i+1, x_labels.at(i).c_str());
    	
		th1->StatOverflows(1);
    	th1->Sumw2();
    	
      int bin=0;
      for(unsigned int i=0; i<x_labels.size(); i++){
	if(label==x_labels.at(i)){
	  bin=i+1;
	  th1->Fill(bin, weight);
	}
      }
      if(bin==0) th1->Fill(bin, weight);

		return th1;	
	}
	else if(obj->InheritsFrom(th1cl)) {
		TH1* th1 = static_cast<TH1*>(obj);
      int bin=0;
      for(unsigned int i=0; i<x_labels.size(); i++){
	if(label==x_labels.at(i)){
	  bin=i+1;
	  th1->Fill(bin, weight);
	}
      }
      if(bin==0) th1->Fill(bin, weight);
		return th1;
	}

	TString msg("Analysis_JetMET_Base: ERROR No valid object ");
	msg.Append(name);
	msg.Append(" found");
	Abort(msg);		

	return 0;
}


TH2* Analysis_JetMET_Base::Fill(AnaKey name, double XVal, double YVal, double weight, 
	    int NumBinsX, double XMin, double XMax, int NumBinsY, double YMin, double YMax){

  	static TClass* th2cl = TH2::Class();
	
	AnaKey histoname(name);

	TObject* obj = OutputDir()->Get(histoname);
	if(!obj) {
		TDirectory* bkp = gDirectory;
		OutputDir()->cd();
    	TH2D* th2 = new TH2D(histoname.Data(), histoname.Data(), NumBinsX, XMin, XMax, NumBinsY, YMin, YMax);
    	gDirectory = bkp;
      cout << "Booked new histo " << histoname.Data() << endl;	
    	
		th2->StatOverflows(1);
    	th2->Sumw2();
		th2->Fill(XVal, YVal, weight);
		return th2;	
	}
	else if(obj->InheritsFrom(th2cl)) {
		TH2* th2 = static_cast<TH2*>(obj);
		th2->Fill(XVal, YVal, weight);
		return th2;
	}

	TString msg("Analysis_JetMET_Base: ERROR No valid object ");
	msg.Append(name);
	msg.Append(" found");
	Abort(msg);		

	return 0;
}



TH2* Analysis_JetMET_Base::Fill(AnaKey name, double XVal, double YVal, double weight, 
	    int NumBinsX, const double* Xbins, int NumBinsY, double YMin, double YMax){

  	static TClass* th2cl = TH2::Class();
	
	AnaKey histoname(name);

	TObject* obj = OutputDir()->Get(histoname);
	if(!obj) {
		TDirectory* bkp = gDirectory;
		OutputDir()->cd();
    	TH2D* th2 = new TH2D(histoname.Data(), histoname.Data(), NumBinsX, Xbins, NumBinsY, YMin, YMax);
    	gDirectory = bkp;
      cout << "Booked new histo " << histoname.Data() << endl;	
    	
		th2->StatOverflows(1);
    	th2->Sumw2();
		th2->Fill(XVal, YVal, weight);
		return th2;	
	}
	else if(obj->InheritsFrom(th2cl)) {
		TH2* th2 = static_cast<TH2*>(obj);
		th2->Fill(XVal, YVal, weight);
		return th2;
	}

	TString msg("Analysis_JetMET_Base: ERROR No valid object ");
	msg.Append(name);
	msg.Append(" found");
	Abort(msg);		

	return 0;
}


TH2* Analysis_JetMET_Base::Fill(AnaKey name, double XVal, double YVal, double weight, 
		      int NumBinsX, double XMin, double XMax, int NumBinsY, const double* Ybins){

  	static TClass* th2cl = TH2::Class();
	
	AnaKey histoname(name);

	TObject* obj = OutputDir()->Get(histoname);
	if(!obj) {
		TDirectory* bkp = gDirectory;
		OutputDir()->cd();
    	TH2D* th2 = new TH2D(histoname.Data(), histoname.Data(), NumBinsX, XMin, XMax, NumBinsY, Ybins);
    	gDirectory = bkp;
      cout << "Booked new histo " << histoname.Data() << endl;	
		
		th2->StatOverflows(1);
    	th2->Sumw2();
		th2->Fill(XVal, YVal, weight);
		return th2;	
	}
	else if(obj->InheritsFrom(th2cl)) {
		TH2* th2 = static_cast<TH2*>(obj);
		th2->Fill(XVal, YVal, weight);
		return th2;
	}

	TString msg("Analysis_JetMET_Base: ERROR No valid object ");
	msg.Append(name);
	msg.Append(" found");
	Abort(msg);		

	return 0;
}


TH2* Analysis_JetMET_Base::Fill(AnaKey name, double XVal, double YVal, double weight, 
			int NumBinsX, const float* Xbins, int NumBinsY, const float* Ybins){

  	static TClass* th2cl = TH2::Class();
	
	AnaKey histoname(name);

	TObject* obj = OutputDir()->Get(histoname);
	if(!obj) {
		TDirectory* bkp = gDirectory;
		OutputDir()->cd();
    	TH2D* th2 = new TH2D(histoname.Data(), histoname.Data(), NumBinsX, Xbins, NumBinsY, Ybins);
    	gDirectory = bkp;
    	
		th2->StatOverflows(1);
    	th2->Sumw2();
		th2->Fill(XVal, YVal, weight);
		return th2;	
	}
	else if(obj->InheritsFrom(th2cl)) {
		TH2* th2 = static_cast<TH2*>(obj);
		th2->Fill(XVal, YVal, weight);
		return th2;
	}

	TString msg("Analysis_JetMET_Base: ERROR No valid object ");
	msg.Append(name);
	msg.Append(" found");
	Abort(msg);		

	return 0;
}


TH2* Analysis_JetMET_Base::Fill(AnaKey name, std::string xlabel, std::string ylabel, double weight,
			std::vector<std::string> x_labels, std::vector<std::string> y_labels){


  	static TClass* th2cl = TH2::Class();
	
	AnaKey histoname(name);

	TObject* obj = OutputDir()->Get(histoname);
	if(!obj) {
		TDirectory* bkp = gDirectory;
		OutputDir()->cd();
    	TH2D* th2 = new TH2D(histoname.Data(), histoname.Data(), x_labels.size(), 0.5, x_labels.size()+0.5, y_labels.size(), 0.5, y_labels.size()+0.5);
		for(unsigned int i=0; i<x_labels.size(); i++)
		th2->GetXaxis()->SetBinLabel(i+1, x_labels.at(i).c_str());

		for(unsigned int i=0; i<y_labels.size(); i++)
		th2->GetYaxis()->SetBinLabel(i+1, y_labels.at(i).c_str());
		
    	gDirectory = bkp;
    	
		th2->StatOverflows(1);
    	th2->Sumw2();

      int xbin=0;
      int ybin=0;
      for(unsigned int i=0; i<x_labels.size(); i++){
	if(xlabel==x_labels.at(i))
	  xbin=i+1;
      }
      for(unsigned int i=0; i<y_labels.size(); i++){
	if(ylabel==y_labels.at(i))
	  ybin=i+1;
      }
      th2->Fill(xbin, ybin, weight);
		return th2;	
	}
	else if(obj->InheritsFrom(th2cl)) {
		TH2* th2 = static_cast<TH2*>(obj);
      int xbin=0;
      int ybin=0;
      for(unsigned int i=0; i<x_labels.size(); i++){
	if(xlabel==x_labels.at(i))
	  xbin=i+1;
      }
      for(unsigned int i=0; i<y_labels.size(); i++){
	if(ylabel==y_labels.at(i))
	  ybin=i+1;
      }
      th2->Fill(xbin, ybin, weight);
		return th2;
	}

	TString msg("Analysis_JetMET_Base: ERROR No valid object ");
	msg.Append(name);
	msg.Append(" found");
	Abort(msg);		

	return 0;
}


TH3* Analysis_JetMET_Base::Fill(AnaKey name, double XVal, double YVal, double ZVal, double weight, 
	    int NumBinsX, double XMin, double XMax, int NumBinsY, double YMin, double YMax, int NumBinsZ, double ZMin, double ZMax){

  	static TClass* th3cl = TH3::Class();
	
	AnaKey histoname(name);

	TObject* obj = OutputDir()->Get(histoname);
	if(!obj) {
		TDirectory* bkp = gDirectory;
		OutputDir()->cd();
    	TH3F* th3 = new TH3F(histoname.Data(), histoname.Data(), NumBinsX, XMin, XMax, NumBinsY, YMin, YMax, NumBinsZ, ZMin, ZMax );
    	gDirectory = bkp;
    	
    	th3->Sumw2();
		th3->Fill(XVal, YVal, ZVal, weight);
		return th3;	
	}
	else if(obj->InheritsFrom(th3cl)) {
		TH3* th3 = static_cast<TH3*>(obj);
		th3->Fill(XVal, YVal, ZVal,  weight);
		return th3;
	}

	TString msg("Analysis_JetMET_Base: ERROR No valid object ");
	msg.Append(name);
	msg.Append(" found");
	Abort(msg);		

	return 0;
}

