#include "BoostedTop.h"
#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

int main(int argc, char** argv)
{
	UInt_t jobnum = atoi(argv[1])+1;
	UInt_t jobcount = atoi(argv[2]);
	string lumifile = string(argv[3]);
	//string outfilename = string(argv[4]);
	stringstream ssoutfilename;
	ssoutfilename << "out" << "_" << jobnum-1 << ".root";

	BoostedTop ana(ssoutfilename.str());
	
	if(lumifile.find("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8") != string::npos) {ana.TTCENTRAL = true;}
	if(lumifile.find("TTToSemiLeptonic") != string::npos) {ana.TTMC_ = true;}
	if(lumifile.find("TTToHadronic") != string::npos) {ana.TTMC_ = true;}
	if(lumifile.find("TTTo2L2Nu") != string::npos) {ana.TTMC_ = true;}
	if(lumifile.find("TT_") != string::npos) {ana.TTMC_ = true;}
	if(lumifile.find("TTJets_") != string::npos) {ana.TTMC_ = true;}
	if(lumifile.find("TT_AK8HT650") != string::npos) {ana.TTPT650MC_ = true;}
	if(lumifile.find("SingleMuon") != string::npos) {ana.DAMU_ = true;}
	if(lumifile.find("SingleElectron") != string::npos) {ana.DAEL_ = true;}
	if(lumifile.find("EGamma") != string::npos) {ana.DAEL_ = true;}
	if(lumifile.find("JetHT") != string::npos) {ana.DAEL_ = true; ana.DAMU_ = true;}
	if(lumifile.find("QCD_HT") != string::npos) {ana.QCDMC_ = true;}
	if(lumifile.find("QCD_Pt") != string::npos && lumifile.find("Enriched") == string::npos) {ana.QCDMC_ = true;}

	if(lumifile.find("hdampDOWN") != string::npos) {ana.mcname_ = "hddown";}
	if(lumifile.find("hdampUP") != string::npos) {ana.mcname_ = "hdup";}
	if(lumifile.find("TuneCP5down") != string::npos) {ana.mcname_ = "tunedown";}
	if(lumifile.find("TuneCP5up") != string::npos) {ana.mcname_ = "tuneup";}
	if(lumifile.find("isrdown") != string::npos) {ana.mcname_ = "isrdown";}
	if(lumifile.find("isrup") != string::npos) {ana.mcname_ = "isrup";}
	if(lumifile.find("fsrdown") != string::npos) {ana.mcname_ = "fsrdown";}
	if(lumifile.find("fsrup") != string::npos) {ana.mcname_ = "fsrup";}
	if(lumifile.find("TuneCUETP8M2T4down") != string::npos) {ana.mcname_ = "tunedown";}
	if(lumifile.find("TuneCUETP8M2T4up") != string::npos) {ana.mcname_ = "tuneup";}

	ana.SetPrintInfo(1000);
    //ana.EnableDuplicateCheck();
	ana.AddLumiFile(lumifile);
	for(int n = 4 ; n < argc ; ++n)
	{
		ana.AddLumiFile(argv[n]);
	}
	ana.Batch_Prepare(jobnum, jobcount);

////	R2 Trigger
//	ana.AddTriggerSelection("IsoMu24", {"HLT_IsoMu24_v.*"});
//	ana.GetTriggerSelection("IsoMu24")->PrintInfo();
//	ana.AddTriggerSelection("IsoTkMu24", {"HLT_IsoTkMu24_v.*"});
//	ana.GetTriggerSelection("IsoTkMu24")->PrintInfo();
//	ana.AddTriggerSelection("IsoMu27", {"HLT_IsoMu27_v.*"});
//	ana.GetTriggerSelection("IsoMu27")->PrintInfo();
//	//ana.AddTriggerSelection("IsoTkMu27", {"HLT_IsoTkMu27_v.*"});
//	//ana.GetTriggerSelection("IsoTkMu27")->PrintInfo();
//	ana.AddTriggerSelection("Mu50", {"HLT_Mu50_v.*"});
//	ana.GetTriggerSelection("Mu50")->PrintInfo();
//	ana.AddTriggerSelection("Mu27", {"HLT_Mu27_v.*"}, true);
//	ana.GetTriggerSelection("Mu27")->PrintInfo();
//
//	ana.AddTriggerSelection("Ele27", {"HLT_Ele27_WPTight_Gsf_v.*"});
//	ana.GetTriggerSelection("Ele27")->PrintInfo();
//	ana.AddTriggerSelection("Ele32", {"HLT_Ele32_WPTight_Gsf_v.*"});
//	ana.GetTriggerSelection("Ele32")->PrintInfo();
//	//ana.AddTriggerSelection("Ele30jet", {"HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v.*"});
//	//ana.GetTriggerSelection("Ele30jet")->PrintInfo();
//	ana.AddTriggerSelection("ELJET", {"HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v.*"});
//	ana.GetTriggerSelection("ELJET")->PrintInfo();
//
//	ana.AddTriggerSelection("JET500", {"HLT_AK8PFJet500_v.*"},true);
//	ana.GetTriggerSelection("JET500")->PrintInfo();
//	ana.AddTriggerSelection("JET450", {"HLT_AK8PFJet450_v.*"},true);
//	ana.GetTriggerSelection("JET450")->PrintInfo();
//	ana.AddTriggerSelection("JET400", {"HLT_AK8PFJet400_v.*"},true);
//	ana.GetTriggerSelection("JET400")->PrintInfo();
//	ana.AddTriggerSelection("JET320", {"HLT_AK8PFJet320_v.*"},true);
//	ana.GetTriggerSelection("JET320")->PrintInfo();
//	ana.AddTriggerSelection("JET260", {"HLT_AK8PFJet260_v.*"},true);
//	ana.GetTriggerSelection("JET260")->PrintInfo();
//
////	ana.AddTriggerSelection("PH50", {"HLT_Photon50_v.*"},true);
////	ana.GetTriggerSelection("PH50")->PrintInfo();
////	ana.AddTriggerSelection("PH120", {"HLT_Photon120_v.*"},true);
////	ana.GetTriggerSelection("PH120")->PrintInfo();
////	ana.AddTriggerSelection("PH175", {"HLT_Photon175_v.*"},true);
////	ana.GetTriggerSelection("PH175")->PrintInfo();



	//vector<string> dmtrigger;
    ////dmtrigger.push_back("HLT_Mu17_Mu8_DZ_v.*");
    //dmtrigger.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v.*");
	//ana.AddTriggerSelection("DMT", dmtrigger);
	//ana.GetTriggerSelection("DMT")->PrintInfo();

//	vector<string> detrigger;
//	detrigger.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v.*");
//	detrigger.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v.*");
//	ana.AddTriggerSelection("DET", detrigger);
//	ana.GetTriggerSelection("DET")->PrintInfo();
	
	ana.Loop();
}

