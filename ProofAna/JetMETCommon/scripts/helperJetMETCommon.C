#include <fstream>
#include <vector>
//#include "../utils/TrigMuonEfficiency/TrigMuonEfficiency/LeptonTriggerSF.h"
//#include "ApplyJetCalibration/ApplyJetCalibration.h"


void WriteBTagCalibObject(TFile* options, TString tagger = "SV0", TString OP = "5_85")
{
	map<std::string,std::string> names;
	names["B"]="default";
	names["C"]="default";
	names["T"]="default";
	names["Light"]="default";

	Analysis::CalibrationDataInterfaceROOT* myObj = new Analysis::CalibrationDataInterfaceROOT(tagger.Data(),"../utils/SUSYTools/data/BTagCalibration.env","../utils/SUSYTools/data/"); // Class used for b-tagging scale factors
	myObj->setEffCalibrationNames(names);
	myObj->setSFCalibrationNames(names);
	myObj->initialize("AntiKt4TopoLCnoJVF",OP.Data(),Analysis::Total);
	//myObj->initialize("AntiKt4Topo",OP.Data(),Analysis::Systematic);

	TString name("myBTagCalib");
	name.Append(tagger);
	name.Append(OP);

	myObj->SetName(name.Data());
	options->cd();
	myObj->Write();

	myObj = new Analysis::CalibrationDataInterfaceROOT(tagger.Data(),"../utils/SUSYTools/data/BTagCalibration.env","../utils/SUSYTools/data/"); // Class used for b-tagging scale factors
	myObj->setEffCalibrationNames(names);
	myObj->setSFCalibrationNames(names);
	myObj->initialize("AntiKt4TopoLCJVF",OP.Data(),Analysis::Total);
	//myObj->initialize("AntiKt4Topo",OP.Data(),Analysis::Systematic);

	name = TString("myJVFBTagCalib");
	name.Append(tagger);
	name.Append(OP);

	myObj->SetName(name.Data());
	options->cd();
	myObj->Write();
}

void WriteGRLObject(TString xmlfile)
{
	///++++++++++++++++++++++++++++++++++++++++++++++++
	/// Good Run List reader
	///++++++++++++++++++++++++++++++++++++++++++++++++
	Root::TGoodRunsListReader reader;   
	reader.SetXMLFile(TString("../config/grl/").Append(xmlfile));
	reader.Interpret();
	Root::TGoodRunsList* m_GRL = new Root::TGoodRunsList(reader.GetMergedGoodRunsList());  
	m_GRL->SetName("myGRL");
	m_GRL->Write();
}

void WriteMuonUtilities(TFile* options){
	MuonSmear::SmearingClass* m_mcp_smear = new MuonSmear::SmearingClass("Data12","muid","q_pT","Rel17.2Repro", "../utils/MuonMomentumCorrections/share/");
	//m_mcp_smear->Initialize("Data12","muid","q_pT","Rel17.2Repro", "../utils/MuonMomentumCorrections/share/");
	//options->WriteObject(m_mcp_smear,"MCPSmear");
	






    std::string unit("GeV");
    std::string muon_file_name("Muid_CB_plus_ST_2012_SF.txt");
    std::string muon_sf_dir("../utils/MuonEfficiencyCorrections/share/"); /// Default path
    Analysis::AnalysisMuonConfigurableScaleFactors::Configuration muon_configuration=Analysis::AnalysisMuonConfigurableScaleFactors::AverageOverRuns;
    cout << "about to create" << endl;
    Analysis::AnalysisMuonConfigurableScaleFactors* MuonSF = new Analysis::AnalysisMuonConfigurableScaleFactors(muon_sf_dir, muon_file_name, unit, muon_configuration);
    cout << "about to init" << endl;
    MuonSF->Initialise();
    //options->WriteObject(MuonSF,"MuonSF");
    
    cout << "done with muon sf's" << endl;

	//MuonTriggerSFTool* myMuonSF = new MuonTriggerSFTool();
	//myMuonSF->initialize_moriond2013("../utils/SUSYTools/data/muon_trigger_sf_2012_AtoL.p1328.root");


    //options->WriteObject(myMuonSF,"LeptonTriggerSF");

    cout << "done with muon trigger sf's " << endl;

   	DataPeriod* myDataPeriod = new DataPeriod();
	myDataPeriod->initialize();
	options->cd();
	myDataPeriod->SetName("myDataPeriod");
	myDataPeriod->Write();


}

void WriteJetCalibrationObjectsMuScan(TFile* options){
	TString mainConfig0 ("../utils/ApplyJetCalibration/data/CalibrationConfigs/muScan2013/Rel17.2_JES_muScan_25ns_0_30.config");
	TString mainConfig40("../utils/ApplyJetCalibration/data/CalibrationConfigs/muScan2013/Rel17.2_JES_muScan_25ns_40_40.config");
	TString mainConfig60("../utils/ApplyJetCalibration/data/CalibrationConfigs/muScan2013/Rel17.2_JES_muScan_25ns_60_60.config");
	TString mainConfig80("../utils/ApplyJetCalibration/data/CalibrationConfigs/muScan2013/Rel17.2_JES_muScan_25ns_80_80.config");
	TString mainConfig140("../utils/ApplyJetCalibration/data/CalibrationConfigs/muScan2013/Rel17.2_JES_muScan_25ns_140_140.config");
	TString mainConfig200("../utils/ApplyJetCalibration/data/CalibrationConfigs/muScan2013/Rel17.2_JES_muScan_25ns_200_200.config");
	options->cd();
	
	JetCalibrationTool* mc4LCCalibrator_25ns_0_30 = new JetCalibrationTool();
	mc4LCCalibrator_25ns_0_30->init("AntiKt4LCTopo", mainConfig0, false);
	mc4LCCalibrator_25ns_0_30->SetName("mc4LCCalibrator_25ns_0_30");
	options->cd();
	mc4LCCalibrator_25ns_0_30->Write();

	JetCalibrationTool* mc4LCCalibrator_25ns_40_40 = new JetCalibrationTool();
	mc4LCCalibrator_25ns_40_40->init("AntiKt4LCTopo", mainConfig40, false);
	mc4LCCalibrator_25ns_40_40->SetName("mc4LCCalibrator_25ns_40_40");
	options->cd();
	mc4LCCalibrator_25ns_40_40->Write();
	
    JetCalibrationTool* mc4LCCalibrator_25ns_60_60 = new JetCalibrationTool();
	mc4LCCalibrator_25ns_60_60->init("AntiKt4LCTopo", mainConfig60, false);
	mc4LCCalibrator_25ns_60_60->SetName("mc4LCCalibrator_25ns_60_60");
	options->cd();
	mc4LCCalibrator_25ns_60_60->Write();

	JetCalibrationTool* mc4LCCalibrator_25ns_80_80 = new JetCalibrationTool();
	mc4LCCalibrator_25ns_80_80->init("AntiKt4LCTopo", mainConfig80, false);
	mc4LCCalibrator_25ns_80_80->SetName("mc4LCCalibrator_25ns_80_80");
	options->cd();
	mc4LCCalibrator_25ns_80_80->Write();
    
    JetCalibrationTool* mc4LCCalibrator_25ns_140_140 = new JetCalibrationTool();
	mc4LCCalibrator_25ns_140_140->init("AntiKt4LCTopo", mainConfig140, false);
	mc4LCCalibrator_25ns_140_140->SetName("mc4LCCalibrator_25ns_140_140");
	options->cd();
	mc4LCCalibrator_25ns_140_140->Write();

	JetCalibrationTool* mc4LCCalibrator_25ns_200_200 = new JetCalibrationTool();
	mc4LCCalibrator_25ns_200_200->init("AntiKt4LCTopo", mainConfig200, false);
	mc4LCCalibrator_25ns_200_200->SetName("mc4LCCalibrator_25ns_200_200");
	options->cd();
	mc4LCCalibrator_25ns_200_200->Write();

}

void WriteJetCalibrationObjects(TFile* options){
	TString mainConfig("../utils/ApplyJetCalibration/data/CalibrationConfigs/JES_Full2012dataset_Preliminary_Jan13.config");
	TString afiiConfig("../utils/ApplyJetCalibration/data/CalibrationConfigs/JES_Full2012dataset_Preliminary_AFII_Jan13.config");
	options->cd();

	JetCalibrationTool* data4LCCalibrator = new JetCalibrationTool();
	data4LCCalibrator->init("AntiKt4LCTopo", mainConfig, true);
	data4LCCalibrator->SetName("data4LCCalibrator");
	options->cd();
	data4LCCalibrator->Write();

	JetCalibrationTool* data4EMCalibrator = new JetCalibrationTool();
	data4EMCalibrator->init("AntiKt4TopoEM", mainConfig, true);
	data4EMCalibrator->SetName("data4EMCalibrator");
	options->cd();
	data4EMCalibrator->Write();

	JetCalibrationTool* data6LCCalibrator = new JetCalibrationTool();
	data6LCCalibrator->init("AntiKt6LCTopo", mainConfig, true);
	data6LCCalibrator->SetName("data6LCCalibrator");
	options->cd();
	data6LCCalibrator->Write();

	JetCalibrationTool* data6EMCalibrator = new JetCalibrationTool();
	data6EMCalibrator->init("AntiKt6TopoEM", mainConfig, true);
	data6EMCalibrator->SetName("data6EMCalibrator");
	options->cd();
	data6EMCalibrator->Write();




	JetCalibrationTool* mc4LCCalibrator = new JetCalibrationTool();
	mc4LCCalibrator->init("AntiKt4LCTopo", mainConfig, false);
	mc4LCCalibrator->SetName("mc4LCCalibrator");
	options->cd();
	mc4LCCalibrator->Write();


	JetCalibrationTool* mc4EMCalibrator = new JetCalibrationTool();
	mc4EMCalibrator->init("AntiKt4TopoEM", mainConfig, false);
	mc4EMCalibrator->SetName("mc4EMCalibrator");
	options->cd();
	mc4EMCalibrator->Write();

	JetCalibrationTool* mc6LCCalibrator = new JetCalibrationTool();
	mc6LCCalibrator->init("AntiKt6LCTopo", mainConfig, false);
	mc6LCCalibrator->SetName("mc6LCCalibrator");
	options->cd();
	mc6LCCalibrator->Write();

	JetCalibrationTool* mc6EMCalibrator = new JetCalibrationTool();
	mc6EMCalibrator->init("AntiKt6TopoEM", mainConfig, false);
	mc6EMCalibrator->SetName("mc6EMCalibrator");
	options->cd();
	mc6EMCalibrator->Write();


 // cout << "testing with Pascal mc4LC" << endl;

  /*TLorentzVector vec = mc4LCCalibrator->ApplyJetAreaOffsetEtaJES(127.39, -0.857549, 1.0, 18.4123., -0.405253, 0.219108, -0.45514, 0.660416,
                                                                 10.4646,
                                                                 17,
                                                                 14);

*/

 // mc4LCCalibrator->UseGeV(true);
  //TLorentzVector vec = mc4LCCalibrator->ApplyJetAreaOffsetEtaJES(50., 0., 0., 10., -0.405253, 0.219108, -0.45514, 0.660416,
                                                                 //10.,
                                                                 //20.,
                                                                 //15.);


  //cout << "test is " << vec.Pt() << endl;
  

    JetCalibrationTool* mcAFII4LCCalibrator = new JetCalibrationTool();
    mcAFII4LCCalibrator->init("AntiKt4LCTopo", afiiConfig, false);
	mcAFII4LCCalibrator->SetName("mcAFII4LCCalibrator");
	options->cd();
	mcAFII4LCCalibrator->Write();

	JetCalibrationTool* mcAFII4EMCalibrator = new JetCalibrationTool();
	mcAFII4EMCalibrator->init("AntiKt4TopoEM", afiiConfig, false);
	mcAFII4EMCalibrator->SetName("mcAFII4EMCalibrator");
	options->cd();
	mcAFII4EMCalibrator->Write();

	// JetCalibrationTool* mcAFII6LCCalibrator = new JetCalibrationTool();
	// mcAFII6LCCalibrator->init("AntiKt6LCTopo", afiiConfig, false);
	// mcAFII6LCCalibrator->SetName("mcAFII6LCCalibrator");
	// mcAFII6LCCalibrator->Write();

	// JetCalibrationTool* mcAFII6EMCalibrator = new JetCalibrationTool();
	// mcAFII6EMCalibrator->init("AntiKt6TopoEM", afiiConfig, false);
	// mcAFII6EMCalibrator->SetName("mcAFII6EMCalibrator");
	// mcAFII6EMCalibrator->Write();
 
 	TString multijes_file = "JES_2012/Moriond2013/MultijetJES_2012.config";
  	TString jes_file = "JES_2012/Moriond2013/InsituJES2012_14NP.config";
  	TString jes_path = "../utils/JetUncertainties/share/";

  	TString mcType = "MC12a";
  	TString mcTypeA = "AFII";

  	TString m_jetAlgo = "AntiKt4LCTopo";
    MultijetJESUncertaintyProvider* m_myJes=new MultijetJESUncertaintyProvider(multijes_file,jes_file,m_jetAlgo,mcType,jes_path);
    m_myJes->SetName("JES");
    options->cd();
    m_myJes->Write();


    MultijetJESUncertaintyProvider* m_myJesA=new MultijetJESUncertaintyProvider(multijes_file,jes_file,m_jetAlgo,mcTypeA,jes_path);
    m_myJesA->SetName("JESA");
    options->cd();
    m_myJesA->Write();

  TString jetAlgo   = "AntiKt4LCTopo";
  TString jer_file = "../utils/JetResolution/share/JERProviderPlots_2012.root";
JetSmearingTool* m_myJER = new JetSmearingTool(jetAlgo,jer_file);
  m_myJER->init();
  m_myJER->SetName("JER");
  options->cd();
  m_myJER->Write();




}

void WritePRWConfigure(TFile* options, TString period){
	Root::TPileupReweighting* m_pileupTool = new Root::TPileupReweighting();
    m_pileupTool->UsePeriodConfig(period); //for MC11b change the string to "MC11b". For MC12a, change to "MC12a"
    m_pileupTool->Initialize();
    m_pileupTool->SetName("PRWConfigure");
    options->cd();
    m_pileupTool->Write();
}

void WriteTopPRWO(TFile* options, TString file){
  Root::TPileupReweighting* myPile = new Root::TPileupReweighting();
//  myPile->AddConfigFile("../config/lumicalcs_R17/gj.prw.root");
  myPile->AddConfigFile("../utils/PileupReweighting/share/mc12a_defaults.prw.root");
  myPile->AddLumiCalcFile("../config/"+file+".root");
  myPile->SetUnrepresentedDataAction(2); //Action needs investigation
  myPile->DisableWarnings(true);
  myPile->Initialize();
  myPile->SetName(file);
  options->cd();
  myPile->Write();
}

void WriteZllPRWO(TFile* options, TString file){
  Root::TPileupReweighting* myPile = new Root::TPileupReweighting();
//  myPile->AddConfigFile("../config/lumicalcs_R17/gj.prw.root");
  myPile->SetDataScaleFactors(1./1.11) ;
  myPile->AddConfigFile("../utils/PileupReweighting/share/mc12a_defaults.prw.root");
  myPile->AddLumiCalcFile("../config/"+file+".root");
  myPile->SetUnrepresentedDataAction(2); //Action needs investigation
  myPile->DisableWarnings(true);
  myPile->Initialize();
  myPile->SetName(file);
  options->cd();
  myPile->Write();
}

void WriteDijetPRWO(TFile* options, TString file, TString prepend = "PeriodAB_"){
  Root::TPileupReweighting* myPile = new Root::TPileupReweighting();
//  myPile->AddConfigFile("../config/lumicalcs_R17/gj.prw.root");
  myPile->AddConfigFile("../utils/PileupReweighting/share/mc12a_defaults.prw.root");
  myPile->AddLumiCalcFile("../config/"+prepend+file+"_lumi.root");
  myPile->SetUnrepresentedDataAction(2); //Action needs investigation
  myPile->DisableWarnings(true);
  myPile->Initialize();
  myPile->SetName(file);
  options->cd();
  myPile->Write();
}

void WriteGroomedPRWO(TFile* options, TString file){
  Root::TPileupReweighting* myPile = new Root::TPileupReweighting();
//  myPile->AddConfigFile("../config/lumicalcs_R17/gj.prw.root");
  myPile->AddConfigFile("../utils/PileupReweighting/share/mc12a_defaults.prw.root");
  myPile->AddLumiCalcFile("../config/groomedPRW/"+file+".root");
  myPile->SetUnrepresentedDataAction(2); //Action needs investigation
  myPile->DisableWarnings(true);
  myPile->Initialize();
  myPile->SetName(file);
  options->cd();
  myPile->Write();
}

