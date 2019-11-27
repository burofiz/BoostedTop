#include "BoostedTop.h"
#include <iostream>
#include <sstream>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

#include <LV.h>
#include <ConfigParser.h>
#include <NeutrinoSolver.h>
#include "NNJet.h"


#include "BAnalysis.h"

using namespace std;
using namespace TMath;

bool GenSelect(const TTEvent& ev)
{
	//if(ev.LLep()->Pt() < 30 || abs(ev.LLep()->Eta()) > 2.4) {return false;}
	//if(ev.BLep()->Pt() < 30 || abs(ev.BLep()->Eta()) > 2.4) {return false;}

	//if(ev.THad()->Pt() > 200. && abs(ev.THad()->Eta()) < 2.4) {return true;}

	if(ev.JaHad()->Pt() < 0 || abs(ev.JaHad()->Eta()) > 2.4) {return false;}
	if(ev.JbHad()->Pt() < 0 || abs(ev.JbHad()->Eta()) > 2.4) {return false;}
	if(ev.BHad()->Pt() < 0 || abs(ev.BHad()->Eta()) > 2.4) {return false;}

	return true;

}

BoostedTop::BoostedTop(string outfilename) :
		unfolding_all_("unfolding"),
		unfolding_ps_("unfolding_ps"),
		hgenPS("genPS"),
		hgenAll("genAll"),
		hgenLHE("genLHE"),
		hgenHadAll("genHadAll"),
		hgenHadRes("genHadRes"),
		hgenHadBoosted("genHadBoosted"),
		hgenHadBoostedMerged("genHadBoostedMerged"),
		hgenHadBroken("genHadBroken"),
		hgenHadBrokenRes("genHadBrokenRes"),
		hgenHadNo("genHadNo"),
		hgenLepRes("genLepRes"),
		hgenLepBoosted("genLepBoosted"),
		hgenLepNo("genLepNo"),
		hgenRES("genRES"),
		hgenOther("genOther"),
		hgenBH("genBH"),
		hgenBL("genBL"),
		hgenBHBL("genBHBL"),
		hrecRESright("recRESright"),
		hrecRESwrong("recRESwrong"),
		hrecRESnonreco("recRESnonreco"),
		hrecRESother("recRESother"),
		hrecBHright("recBHright"),
		hrecBHwrong("recBHwrong"),
		hrecBHnonreco("recBHnonreco"),
		hrecBHother("recBHother"),
		hrecBLright("recBLright"),
		hrecBLwrong("recBLwrong"),
		hrecBLnonreco("recBLnonreco"),
		hrecBLother("recBLother"),
		hrecBHBLright("recBHBLright"),
		hrecBHBLwrong("recBHBLwrong"),
		hrecBHBLnonreco("recBHBLnonreco"),
		hrecBHBLother("recBHBLother"),
		lepjetmu_right("lepjetmu_right"),
		lepjetel_right("lepjetel_right"),
		lepjetmu_wrong("lepjetmu_wrong"),
		lepjetel_wrong("lepjetel_wrong"),
		hadjet_right("hadjet_right"),
		hadjet_wrong("hadjet_wrong"),
		hadjet_other("hadjet_other"),
		genhadjet_right("genhadjet_right"),
		genhadjet_wrong("genhadjet_wrong"),
		genhadjet_other("genhadjet_other"),
		hadjetisolep_right("hadjetisolep_right"),
		hadjetisolep_wrong("hadjetisolep_wrong"),
		hadjetisolep_other("hadjetisolep_other"),
		hadjetnonisolep_right("hadjetnonisolep_right"),
		hadjetnonisolep_wrong("hadjetnonisolep_wrong"),
		hadjetnonisolep_other("hadjetnonisolep_other"),
		hadjetbkgPH("hadjetbkgPH"),
		hadjetbkgB("hadjetbkgB"),
		_hists1d("extra"),
		_hists2d("extra")
{
	std::cout<<"BT 97"<<endl;
	CP = make_unique<ConfigParser>("config.cfg");

	BAn = make_unique<BAnalysis>();

	BTAG_EFF_MODE_ = CP->Get<bool>("BTag_Eff_Mode");
	ISRunc_ = CP->Get<int>("ISRunc");
	if(ISRunc_ == -1.){mcname_ = "isrdown";}
	else if(ISRunc_ == 1.){mcname_ = "isrup";}
	FSRunc_ = CP->Get<int>("FSRunc");
	if(FSRunc_ == -1.){mcname_ = "fsrdown";}
	else if(FSRunc_ == 1.){mcname_ = "fsrup";}
	FSunc_ = CP->Get<int>("FSunc");
	RSunc_ = CP->Get<int>("RSunc");
	PS_MODE_ = CP->Get<bool>("PseudoTop_Mode");
	ISOSIDEBAND_ = CP->Get<bool>("IsoSideband_Mode");
	SKIPPT650_ = CP->Get<bool>("SkipPT650");
	HadTopJet::PS_MODE = PS_MODE_;
	ctopptweight_ = CP->Get<double>("TopPtWeight");
	cttmweight_ = CP->Get<double>("TTMassWeight");
	TRAIN = CP->Get<bool>("FillTrainTrees");
	DELPHES = CP->Get<bool>("DELPHES");
	PDFUNCS_ = CP->Get<bool>("PDFUNCS");
	OJet::DELPHES = DELPHES;
	UserTime(CP->Get<int>("Time"));
	string likelihoodfile = "likelihoods_pa.root";
	if(PS_MODE_){ likelihoodfile = "likelihoods_ps.root";}


	if(CP->Get<double>("Unc_BDecay") == -1.) {mcname_ = "bdecaydown";} 
	else if(CP->Get<double>("Unc_BDecay") == 1.) {mcname_ = "bdecayup";} 
	if(CP->Get<double>("Unc_BFrag") == -1.) {mcname_ = "bfragdown";} 
	else if(CP->Get<double>("Unc_BFrag") == 1.) {mcname_ = "bfragup";} 

	cfillsigtree_ = CP->Get<int>("FillSigTree");
	cfillbkgtree_ = CP->Get<int>("FillBkgTree");
	btag_loose_ = CP->Get<double>("BTag_Loose");
	btag_medium_ = CP->Get<double>("BTag_Medium");
	btag_tight_ = CP->Get<double>("BTag_Tight");

	cbsidebandmin_ = CP->Get<double>("BTagSideband_Min");
	cbsidebandmax_ = CP->Get<double>("BTagSideband_Max");
	if(cbsidebandmin_ == cbsidebandmax_) //no sideband selection use default
	{
		cbsidebandmin_ = btag_medium_;
		cbsidebandmax_ = 100.;
	}
	std::cout<<"BT 144"<<endl;

	cuncjer_ = CP->Get<double>("JEC_UncJER");
	cuncjes_ = CP->Get<double>("JEC_UncJES");

	bdecayweights.Init(CP->Get<double>("Unc_BDecay"));
    bfragweights.Init("bfragweights.root", CP->Get<double>("Unc_BFrag"));


	string puhistname("pu_central");
    if(CP->Get<int>("PU_Unc") == -1) puhistname = "pu_minus";
    if(CP->Get<int>("PU_Unc") == 1) puhistname = "pu_plus";
    TFile* f = TFile::Open(CP->Get<string>("PU_WeightFile").c_str());
    puhist_ = dynamic_cast<TH1D*>(f->Get(puhistname.c_str()));

    if(UserTime() == 2017 || UserTime() == 2016){lepcorrector.init2(CP->Get<string>("LEP_SFMu"), CP->Get<string>("LEP_SFEl"));}
	else {lepcorrector.init(CP->Get<string>("LEP_SFMu"), CP->Get<string>("LEP_SFEl"));}
    //elcorrector.init(CP->Get<string>("LEP_SFEl"), 0.0);
	csigmamu_ = CP->Get<double>("LEP_UncMu");
	csigmael_ = CP->Get<double>("LEP_UncEl");

	trgres.emplace("IsoMu24" , "HLT_IsoMu24_v.*");
	trgres.emplace("IsoTkMu24" , "HLT_IsoTkMu24_v.*");
	trgres.emplace("IsoMu27" , "HLT_IsoMu27_v.*");
	trgres.emplace("Mu50" , "HLT_Mu50_v.*");
	trgres.emplace("Mu27" , "HLT_Mu27_v.*");
	trgres.emplace("Ele27" , "HLT_Ele27_WPTight_Gsf_v.*");
	trgres.emplace("Ele32" , "HLT_Ele32_WPTight_Gsf_v.*");
	trgres.emplace("ElJet" , "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v.*");
	trgres.emplace("Jet260" , "HLT_AK8PFJet260_v.*");
	trgres.emplace("Jet320" , "HLT_AK8PFJet320_v.*");
	trgres.emplace("Jet400" , "HLT_AK8PFJet400_v.*");
	trgres.emplace("Jet450" , "HLT_AK8PFJet450_v.*");
	trgres.emplace("Jet500" , "HLT_AK8PFJet500_v.*");

	netmu1 = new NeuralNet("NNmu1.root");
	netmu2 = new NeuralNet("NNmu2.root");

	nethad1 = new NeuralNet("NNhad1.root");
	nethad2 = new NeuralNet("NNhad2.root");
	nethad3 = new NeuralNet("NNhad3.root");
	nethad4 = new NeuralNet("NNhad4.root");
	//nethad5 = new NeuralNet("NNhad5.root");

	if(UserTime() == 2016 || UserTime() == 2017) {prefireprob.Init(CP->Get<string>("prefiremap"), CP->Get<double>("prefire_sigma"));}

	TFile* tf = new TFile(likelihoodfile.c_str());
	h_tmwm = dynamic_cast<TH2D*>(tf->Get("had_tmwm"));
	h_dnu = dynamic_cast<TH1D*>(tf->Get("had_dnu"));

	if(UserTime() == 2023)
	{
		//bins_toppt = {0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 430.0, 500., 800.,};
		//bins_topy = {0.0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.0, 2.5};
		//bins_ttm = {300.0, 360.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 1200.0, 1500.0, 2500.0};
		//bins_ttpt = {0.0, 40.0, 80.0, 150.0, 220.0, 300.0, 380.0, 500.0, 1000.0};	
		//bins_tty = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4};
		bins_toppt_large = {0, 160, 350, 450, 7000.};
		bins_toppt = {0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 430.0, 500., 800.,};
		bins_topy = {0.0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0, 3.4};
		bins_ttm = {300.0, 360.0, 430.0, 500.0, 580.0, 680.0, 800.0, 1000.0, 1200.0, 1500.0, 2500.0};
		bins_ttpt = {0.0, 40.0, 80.0, 150.0, 220.0, 300.0, 380.0, 500.0, 1000.0};	
		bins_tty = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.4, 2.6, 3.0};
		bins_cts = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
	}
	else
	{
		bins_toppt_large = {0, 160, 350, 450, 7000.};
		//bins_toppt = {0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500, 600, 700, 800, 900, 1000, 1200, 1500};
		bins_toppt = {0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500, 600, 700, 800, 900, 1000, 1500};
		bins_topptsoft = {0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500, 600, 700, 800, 1500};
		bins_topy = {0.0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.0, 2.5};
		//bins_ttm = {320.0, 400., 500.0, 600.0, 700.0, 800.0, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 3500};
		bins_ttm = {320.0, 400., 500.0, 600.0, 700.0, 800.0, 1000, 1200, 1400, 1700, 2000, 3500};
		//bins_ttm = {300.0, 380.0, 480.0, 580.0, 680.0, 800.0, 1000, 1200, 1400, 1700, 2000, 3500};
		bins_ttpt = {0.0, 40.0, 80.0, 150.0, 220.0, 300.0, 380.0, 500.0, 600., 1000.0};
		bins_tty = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4};
		bins_cts = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
	}
	std::cout<<"BT 223"<<endl;

	histfile = new TFile(outfilename.c_str(), "recreate");
	histfile->cd();
	btageff.Init(btag_medium_, CP->Get<double>("BTagSideband_Min"), CP->Get<double>("BTagSideband_Max"));

	topdatasize = 5;
	topdata = new float[topdatasize];
	topmutreesig = new TTree("topmutreesig", "topmutreesig");
	topmutreesig->Branch("data", topdata, ("data["+to_string(topdatasize)+"]/F").c_str());
	topmutreesig->Branch("pt", &phpt, "pt/F");
	topmutreesig->Branch("eta", &pheta, "eta/F");
	topmutreebkg = new TTree("topmutreebkg", "topmutreebkg");
	topmutreebkg->Branch("data", topdata, ("data["+to_string(topdatasize)+"]/F").c_str());
	topmutreebkg->Branch("pt", &phpt, "pt/F");
	topmutreebkg->Branch("eta", &pheta, "eta/F");

	topeltreesig = new TTree("topeltreesig", "topeltreesig");
	topeltreesig->Branch("data", topdata, ("data["+to_string(topdatasize)+"]/F").c_str());
	topeltreesig->Branch("pt", &phpt, "pt/F");
	topeltreesig->Branch("eta", &pheta, "eta/F");
	topeltreebkg = new TTree("topeltreebkg", "topeltreebkg");
	topeltreebkg->Branch("data", topdata, ("data["+to_string(topdatasize)+"]/F").c_str());
	topeltreebkg->Branch("pt", &phpt, "pt/F");
	topeltreebkg->Branch("eta", &pheta, "eta/F");

	tophaddatasize = 21;
	tophaddata = new float[tophaddatasize];
	tophadtreesig = new TTree("tophadtreesig", "tophadtreesig");
	tophadtreesig->Branch("data", tophaddata, ("data["+to_string(tophaddatasize)+"]/F").c_str());
	tophadtreesig->Branch("pt", &phpt, "pt/F");
	tophadtreesig->Branch("eta", &pheta, "eta/F");
	tophadtreesigbroken = new TTree("tophadtreesigbroken", "tophadtreesigbroken");
	tophadtreesigbroken->Branch("data", tophaddata, ("data["+to_string(tophaddatasize)+"]/F").c_str());
	tophadtreesigbroken->Branch("pt", &phpt, "pt/F");
	tophadtreesigbroken->Branch("eta", &pheta, "eta/F");
	tophadtreebkg = new TTree("tophadtreebkg", "tophadtreebkg");
	tophadtreebkg->Branch("data", tophaddata, ("data["+to_string(tophaddatasize)+"]/F").c_str());
	tophadtreebkg->Branch("pt", &phpt, "pt/F");
	tophadtreebkg->Branch("eta", &pheta, "eta/F");

	ttbartree = new TTree("ttbartree", "ttbartree");
	ttbartree->Branch("tleppt", &tleppt, "tleppt/F");
	ttbartree->Branch("thadpt", &thadpt, "thadpt/F");
	ttbartree->Branch("thardpt", &thardpt, "thardpt/F");
	ttbartree->Branch("tsoftpt", &tsoftpt, "tsoftpt/F");
	ttbartree->Branch("toppt", &toppt, "toppt/F");
	ttbartree->Branch("topbarpt", &topbarpt, "topbarpt/F");
	ttbartree->Branch("thady", &thady, "thady/F");
	ttbartree->Branch("tlepy", &tlepy, "tlepy/F");
	ttbartree->Branch("topy", &topy, "topy/F");
	ttbartree->Branch("topbary", &topbary, "topbary/F");
	ttbartree->Branch("cts", &cts, "cts/F");
	ttbartree->Branch("ttm", &ttm, "ttm/F");
	ttbartree->Branch("ttpt", &ttpt, "ttpt/F");
	ttbartree->Branch("tty", &tty, "tty/F");
	ttbartree->Branch("thadpt+ttm+tty", &thadpt_ttm_tty, "thadpt+ttm+tty/F");
	ttbartree->Branch("thadpt+ttm+cts", &thadpt_ttm_cts, "thadpt+ttm+cts/F");
	ttbartree->Branch("thadptref", &thadptref, "thadptref/F");
	ttbartree->Branch("thadetaref", &thadetaref, "thadetaref/F");
	ttbartree->Branch("btag", &btag, "btaf/F");
	ttbartree->Branch("nnres", &nnres, "nnres/F");
	ttbartree->Branch("type", &type, "type/I");
	ttbartree->Branch("nvtx", &nvtx, "nvtx/I");
	ttbartree->Branch("weight", &weight, "weight/D");

	qcdbkgtree = new TTree("qcdbkgtree", "qcdbkgtree");
	qcdbkgtree->Branch("thadptref", &thadptref, "thadptref/F");
	qcdbkgtree->Branch("thadetaref", &thadetaref, "thadetaref/F");
	qcdbkgtree->Branch("nnres", &nnres, "nnres/F");
	qcdbkgtree->Branch("weight", &weight, "weight/D");

	wbkgtree = new TTree("wbkgtree", "wbkgtree");
	wbkgtree->Branch("thadptref", &thadptref, "thadptref/F");
	wbkgtree->Branch("thadetaref", &thadetaref, "thadetaref/F");
	wbkgtree->Branch("nnres", &nnres, "nnres/F");
	wbkgtree->Branch("weight", &weight, "weight/D");

	
	_hists1d.AddHist("genttdecay", 10, 0, 10, "decay type", "");

	_hists1d.AddHist("counter", 100, 0, 100, "counter", "");
	_hists1d.AddHist("mcinfo", 5, 0, 5, "mcinfo", "");
	_hists1d.AddHist("nhadjets_iso", 10, 0, 10, "hadjets iso", "");
	_hists1d.AddHist("nhadjets_noniso", 10, 0, 10, "hadjets noniso", "");
	_hists1d.AddHist("mcweight", 2000, -10, 10, "MC weight", "");
	_hists1d.AddHist("puweight", 1000, 0, 10, "PU weight", "");
	_hists1d.AddHist("lepweight", 1000, 0, 2, "lepton weight", "");
	_hists1d.AddHist("prefireweight", 1000, 0, 2, "lepton weight", "");
	_hists1d.AddHist("btagweight", 1000, 0, 2, "b-tagging weight", "");
	_hists2d.AddHist("btaghigh_deepcsv", 100, 0, 1, {30, 50, 70, 100, 140, 200, 300, 600}, "deepcsv_{high}", "P_{T}(b jet) [GeV]");
	_hists2d.AddHist("btaglow_deepcsv", 100, 0, 1, {30, 50, 70, 100, 140, 200, 300, 600}, "deepcsv_{low}", "P_{T}(b jet) [GeV]");
	_hists1d.AddHist("mu", 200, 0, 200, "#mu", "");
	_hists1d.AddHist("mu_weighted", 200, 0, 200, "#mu", "");
	_hists1d.AddHist("muon_iso", 100, 0, 1, "muon iso", "");
	_hists1d.AddHist("el_iso", 100, 0, 1, "el iso", "");
	_hists1d.AddHist("jet_pt_btagged0", 100, 0, 500, "btagged jets", "");
	_hists1d.AddHist("jet_pt_btagged1", 100, 0, 500, "btagged jets", "");
	_hists1d.AddHist("jet_pt_btagged2", 100, 0, 500, "btagged jets", "");
	_hists1d.AddHist("jet_pt_btagged3", 100, 0, 500, "btagged jets", "");
	_hists1d.AddHist("jet_pt_btagged4", 100, 0, 500, "btagged jets", "");
	_hists1d.AddHist("jet_pt_allbtagged", 100, 0, 500, "nonbtagged jets", "");
	_hists1d.AddHist("jet_pt_missbtagged0", 100, 0, 500, "missbtagged jets", "");
	_hists1d.AddHist("jet_pt_missbtagged1", 100, 0, 500, "missbtagged jets", "");
	_hists1d.AddHist("jet_pt_missbtagged2", 100, 0, 500, "missbtagged jets", "");
	_hists1d.AddHist("jet_pt_missbtagged3", 100, 0, 500, "missbtagged jets", "");
	_hists1d.AddHist("jet_pt_missbtagged4", 100, 0, 500, "missbtagged jets", "");
	_hists1d.AddHist("jet_pt_allmissbtagged", 100, 0, 500, "nonmissbtagged jets", "");
	_hists2d.AddHist("jet_nbjets", 10, 0, 10, 10, 0, 10, "nbjets", "nljets");
	_hists1d.AddHist("selected_bjets_pt", 100, 0, 500, "p_{T}(b-jets) [GeV]", "");
	_hists1d.AddHist("selected_ljets_pt", 100, 0, 500, "p_{T}(l-jets) [GeV]", "");
	_hists2d.AddHist("selected_bljets", 10, 0, 10, 10, 0, 10, "num. b-jets", "num. l-jets");

	_hists2d.AddHist("jet_Eb_pt_res", 120, 30, 630, 100, -1,1 , "pt", "eta");
	_hists2d.AddHist("jet_Bb_pt_res", 120, 30, 630, 100, -1,1 , "pt", "eta");
	//_hists2d.AddHist("jet_Ec_pt_res", 120, 30, 630, 100, -1,1 , "pt", "eta");
	//_hists2d.AddHist("jet_Bc_pt_res", 120, 30, 630, 100, -1,1 , "pt", "eta");
	_hists2d.AddHist("jet_El_pt_res", 120, 30, 630, 100, -1,1 , "pt", "eta");
	_hists2d.AddHist("jet_Bl_pt_res", 120, 30, 630, 100, -1,1 , "pt", "eta");

	_hists2d.AddHist("trel_eta_pt_all", 50, -2.5, 2.5, 150, 30, 330, "eta", "pt");
	_hists2d.AddHist("trel_eta_pt_passnoniso", 50, -2.5, 2.5, 150, 30, 330, "eta", "pt");
	_hists2d.AddHist("trel_eta_pt_pass", 50, -2.5, 2.5, 150, 30, 330, "eta", "pt");
	_hists2d.AddHist("trel_eta_pt_passiso", 50, -2.5, 2.5, 150, 30, 330, "eta", "pt");
	_hists2d.AddHist("treliso_eta_pt_all", 50, -2.5, 2.5, 150, 30, 330, "eta", "pt");
	_hists2d.AddHist("treliso_eta_pt_passnoniso", 50, -2.5, 2.5, 150, 30, 330, "eta", "pt");
	_hists2d.AddHist("treliso_eta_pt_pass", 50, -2.5, 2.5, 150, 30, 330, "eta", "pt");
	_hists2d.AddHist("trjet_eta_pt_all", 50, -2.5, 2.5, 150, 100, 400, "eta", "pt");
	_hists2d.AddHist("trjet_eta_pt_pass", 50, -2.5, 2.5, 150, 100, 400, "eta", "pt");
	_hists2d.AddHist("trmu_eta_pt_all", 50, -2.5, 2.5, 150, 30, 330, "eta", "pt");
	_hists2d.AddHist("trmu_eta_pt_passnoniso", 50, -2.5, 2.5, 150, 30, 330, "eta", "pt");
	_hists2d.AddHist("trmu_eta_pt_pass", 50, -2.5, 2.5, 150, 30, 330, "eta", "pt");

	_hists2d.AddHist("ptgen_drmaxgen", 400, 0, 2000, 300, 0,3 , "pt", "drmax");
	_hists2d.AddHist("ptgen_drmingen", 400, 0, 2000, 300, 0,3 , "pt", "drmin");
	_hists2d.AddHist("BLmu_pt_Dpt", 150, 0, 1500, 200, -1, 1, "", "");
	_hists2d.AddHist("BLmu_pt_Ddr", 150, 0, 1500, 150, 0, 1.5, "", "");
	_hists2d.AddHist("BLmu_pt_Dpz", 150, 0, 1500, 200, -1, 1, "", "");
	_hists2d.AddHist("BLmu_pt_Dnu", 150, 0, 1500, 500, 0, 2500, "", "");
	_hists2d.AddHist("BLel_pt_Dpt", 150, 0, 1500, 200, -1, 1, "", "");
	_hists2d.AddHist("BLel_pt_Ddr", 150, 0, 1500, 150, 0, 1.5, "", "");
	_hists2d.AddHist("BLel_pt_Dpz", 150, 0, 1500, 200, -1, 1, "", "");
	_hists2d.AddHist("BLel_pt_Dnu", 150, 0, 1500, 500, 0, 2500, "", "");
	_hists2d.AddHist("RES_Dpt_0", 120, 0, 600, 200, 0., 2., "", "");
	_hists2d.AddHist("RES_Dpt_1", 120, 0, 600, 200, 0., 2., "", "");
	_hists2d.AddHist("RES_Dpt_2", 120, 0, 600, 200, 0., 2., "", "");
	_hists2d.AddHist("RES_Dpt_3", 120, 0, 600, 200, 0., 2., "", "");
	_hists2d.AddHist("RES_Dpt_4", 120, 0, 600, 200, 0., 2., "", "");
	_hists2d.AddHist("RES_Dpt_5", 120, 0, 600, 200, 0., 2., "", "");
	_hists2d.AddHist("RES_Dpt", 120, 0, 600, 200, 0., 2., "", "");
	_hists2d.AddHist("BH_Dpt_0", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BH_Dpt_1", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BH_Dpt_2", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BH_Dpt_3", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BH_Dpt_4", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BH_Dpt_5", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BH_Dpt", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BH_Ddr", 150, 0, 1500, 150, 0, 1.5, "", "");
	_hists2d.AddHist("BHbroken_Dpt_0", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BHbroken_Dpt_1", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BHbroken_Dpt_2", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BHbroken_Dpt_3", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BHbroken_Dpt_4", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BHbroken_Dpt_5", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BHbroken_Dpt", 150, 0, 1500, 200, 0., 2., "", "");
	_hists2d.AddHist("BHbroken_Ddr", 150, 0, 1500, 150, 0, 1.5, "", "");
	_hists2d.AddHist("genBH_Dpt", 150, 0, 1500, 200, -1, 1., "", "");
	_hists2d.AddHist("genBH_Ddr", 150, 0, 1500, 150, 0, 1.5, "", "");
	_hists2d.AddHist("genBHbroken_Dpt", 150, 0, 1500, 200, -1, 1., "", "");
	_hists2d.AddHist("genBHbroken_Ddr", 150, 0, 1500, 150, 0, 1.5, "", "");
	_hists2d.AddHist("papshad_600_dpt_m", 200, -1, 1, 200, 0, 400, "", "");
	_hists2d.AddHist("papshad_600_dr_m", 200, 0, 2, 200, 0, 400, "", "");
	//_hists2d.AddHist("papshad_605_dpt_m", 200, -1, 1, 200, 0, 400, "", "");
	//_hists2d.AddHist("papshad_605_dr_m", 200, 0, 2, 200, 0, 400, "", "");
	//_hists2d.AddHist("papshad_604_dpt_m", 200, -1, 1, 200, 0, 400, "", "");
	//_hists2d.AddHist("papshad_604_dr_m", 200, 0, 2, 200, 0, 400, "", "");
	//_hists2d.AddHist("papshad_601_dpt_m", 200, -1, 1, 200, 0, 400, "", "");
	//_hists2d.AddHist("papshad_601_dr_m", 200, 0, 2, 200, 0, 400, "", "");
	_hists2d.AddHist("papslep_dpt_m", 200, -1, 1, 200, 0, 400, "", "");
	_hists2d.AddHist("papslep_dr_m", 200, 0, 2, 200, 0, 400, "", "");
	//_hists2d.AddHist("BH_righthyp", 10, 0, 10, 100, 0, 2000, "", "");
	//_hists2d.AddHist("BH_wronghyp", 10, 0, 10, 100, 0, 2000, "", "");

	_hists2d.AddHist("BB_ptlep_genrec", 200, 0, 2000, 200, 0, 2000, "", "");
	_hists2d.AddHist("BB_pthad_genrec", 200, 0, 2000, 200, 0, 2000, "", "");
	//_hists2d.AddHist("BB_lnet_hnet_right", 100, 0, 1, 100, 0, 1, "", "");
	//_hists2d.AddHist("BB_lnet_hnet_wrong", 100, 0, 1, 100, 0, 1, "", "");
	//_hists2d.AddHist("BB_lnet_hnet_other", 100, 0, 1, 100, 0, 1, "", "");
	_hists2d.AddHist("Mt_MW_right", 500, 0, 500, 500, 0, 500, "M_{t} [GeV]", "M_{W} [GeV]");
	_hists2d.AddHist("nusolver_tleppt_respt", 60, 0, 600, 200, -2, 2, "tleppt", "respt");
	_hists2d.AddHist("nusolver_tleppt_resmet", 60, 0, 600, 200, -2, 2, "tleppt", "respt");
	_hists2d.AddHist("nusolver_tleppt_respz", 60, 0, 600, 200, -2, 2, "tleppt", "respz");
	_hists1d.AddHist("Dnu_right", 600, 0, 150, "D_{#nu}", "Events");
	_hists2d.AddHist("Mt_MW_wrong", 500, 0, 500, 500, 0, 500, "M_{t} [GeV]", "M_{W} [GeV]");
	_hists1d.AddHist("Dnu_wrong", 600, 0, 150, "D_{#nu}", "Events");

	_hists1d.AddHist("melel", 500, 0, 500, "melel", "Events");
	_hists2d.AddHist("melel_jet", 10, 0, 10, 250, 0, 500, "jets", "melel");
	_hists2d.AddHist("melel_bjet", 10, 0, 10, 250, 0, 500, "jets", "melel");
	_hists2d.AddHist("melel_genbjet", 10, 0, 10, 250, 0, 500, "jets", "melel");
	_hists1d.AddHist("mmumu", 100, 0, 200, "mmumu", "Events");
	_hists2d.AddHist("mmumu_jet", 10, 0, 10, 250, 0, 500, "jets", "mmumu");
	_hists2d.AddHist("mmumu_bjet", 10, 0, 10, 250, 0, 500, "jets", "mmumu");
	_hists2d.AddHist("mmumu_genbjet", 10, 0, 10, 250, 0, 500, "jets", "mmumu");
	_hists2d.AddHist("prefire_bxm1_ph", 160, -4, 4, 200, 0, 1000, "#eta(jet)", "P_{T}(jet)");
	_hists2d.AddHist("prefire_bxm1_jet", 160, -4, 4, 200, 0, 1000, "#eta(jet)", "P_{T}(jet)");
	_hists2d.AddHist("prefire_all_ph", 160, -4, 4, 200, 0, 1000, "#eta(jet)", "P_{T}(jet)");
	_hists2d.AddHist("prefire_all_jet", 160, -4, 4, 200, 0, 1000, "#eta(jet)", "P_{T}(jet)");

	_hists2d.AddHist("bkg_nnres_nnres", 200, 0, 1, 200, 0, 1, "", "");
	_hists2d.AddHist("pt_jetpt_B", 100, 300, 1300, 100, -1, 1, "P_{T,gen} [GeV]", "Barrel: #Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("pt_brokentopjetpt_B", 100, 300, 1300, 100, -1, 1, "P_{T,gen} [GeV]", "#Barrel: Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("pt_topjetpt_B", 100, 300, 1300, 100, -1, 1, "P_{T,gen} [GeV]", "Barrel: #Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("pt_jetpt_E", 100, 300, 1300, 100, -1, 1, "P_{T,gen} [GeV]", "Endcap: #Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("pt_brokentopjetpt_E", 100, 300, 1300, 100, -1, 1, "P_{T,gen} [GeV]", "Endcap: #Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("pt_topjetpt_E", 100, 300, 1300, 100, -1, 1, "P_{T,gen} [GeV]", "Endcap: #Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("nvtx_pta_jetpt", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("nvtx_pta_jetm", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_pta_jetmtop", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptb_jetpt", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("nvtx_ptb_jetm", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptb_jetmtop", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptc_jetpt", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("nvtx_ptc_jetm", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptc_jetmtop", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_pta_brokentopjetpt", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("nvtx_pta_brokentopjetm", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_pta_brokentopjetmtop", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptb_brokentopjetpt", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("nvtx_ptb_brokentopjetm", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptb_brokentopjetmtop", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptc_brokentopjetpt", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("nvtx_ptc_brokentopjetm", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptc_brokentopjetmtop", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_pta_topjetpt", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("nvtx_pta_topjetm", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_pta_topjetmtop", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptb_topjetpt", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("nvtx_ptb_topjetm", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptb_topjetmtop", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptc_topjetpt", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta P_{T}/P_{T}"); 
	_hists2d.AddHist("nvtx_ptc_topjetm", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 
	_hists2d.AddHist("nvtx_ptc_topjetmtop", 75, 0, 300, 100, -1, 1, "nvtx", "#Delta M/M"); 

	//_hists2d.AddHist("njet_nnres", 10, -0.5, 9.5, 100, 0, 1, "njet", "H_{NN}"); 
	//_hists2d.AddHist("njet_b_nnres", 10, -0.5, 9.5, 100, 0, 1, "njet", "H_{NN}"); 
	//_hists2d.AddHist("njet_ph_nnres", 10, -0.5, 9.5, 100, 0, 1, "njet", "H_{NN}"); 
	_hists2d.AddHist("compare_pt", 100, 300, 1300, 200, -2, 2, "pt", "dpt"); 
	_hists1d.AddHist("compare_dr", 90, 0, 3, "dr", ""); 


	_hists1d.AddHist("jet_xb", 100, 0, 1, "x_{b}", "");

	//_hists1d.AddHist("noniso_b_nnres", 100, 0, 1, "nnres b", "");
	//_hists1d.AddHist("noniso_c_nnres", 100, 0, 1, "nnres c", "");
	//_hists1d.AddHist("noniso_l_nnres", 100, 0, 1, "nnres l", "");
	//_hists1d.AddHist("noniso_o_nnres", 100, 0, 1, "nnres o", "");
	//_hists1d.AddHist("iso_b_nnres", 100, 0, 1, "nnres b", "");
	//_hists1d.AddHist("iso_c_nnres", 100, 0, 1, "nnres c", "");
	//_hists1d.AddHist("iso_l_nnres", 100, 0, 1, "nnres l", "");
	//_hists1d.AddHist("iso_o_nnres", 100, 0, 1, "nnres o", "");
	//_hists1d.AddHist("iso_drjet", 300, 0, 6, "nnres o", "");
	//_hists1d.AddHist("noniso_drjet", 300, 0, 6, "nnres o", "");

	//_hists2d.AddHist("ptbmus", 200, 0, 400, 200, 0, 400, "pt1", "pt2"); 
	histfile->cd();
	TDirectory* dir_BA = histfile->mkdir("BAnalysis");
	dir_BA->cd();
	std::cout<<"BT 490"<<endl;
	BAn->Init();

	histfile->cd();
	TDirectory* dir_gen = histfile->mkdir("BT");
	dir_gen->cd();

	hgenPS.Init();
	hgenAll.Init();
	hgenLHE.Init();
	hgenHadAll.Init();
	hgenHadRes.Init();
	hgenHadBoosted.Init();
	hgenHadBoostedMerged.Init();
	hgenHadBroken.Init();
	hgenHadBrokenRes.Init();
	hgenHadNo.Init();
	hgenLepRes.Init();
	hgenLepBoosted.Init();
	hgenLepNo.Init();
	hgenRES.Init();
	hgenOther.Init();
	hgenBH.Init();
	hgenBL.Init();
	hgenBHBL.Init();


	hrecRESright.Init();
	hrecRESwrong.Init();
	hrecRESnonreco.Init();
	hrecRESother.Init();
	hrecBHright.Init();
	hrecBHwrong.Init();
	hrecBHnonreco.Init();
	hrecBHother.Init();
	hrecBLright.Init();
	hrecBLwrong.Init();
	hrecBLnonreco.Init();
	hrecBLother.Init();
	hrecBHBLright.Init();
	hrecBHBLwrong.Init();
	hrecBHBLnonreco.Init();
	hrecBHBLother.Init();

	lepjetmu_right.Init();
	lepjetel_right.Init();
	lepjetmu_wrong.Init();
	lepjetel_wrong.Init();
	hadjet_right.Init();
	hadjet_wrong.Init();
	hadjet_other.Init();
	genhadjet_right.Init();
	genhadjet_wrong.Init();
	genhadjet_other.Init();
	hadjetisolep_right.Init();
	hadjetisolep_wrong.Init();
	hadjetisolep_other.Init();
	hadjetnonisolep_right.Init();
	hadjetnonisolep_wrong.Init();
	hadjetnonisolep_other.Init();
	hadjetbkgPH.Init();
	hadjetbkgB.Init();

	unfolding_all_.Init(PDFUNCS_);
	unfolding_ps_.Init();

	ElJet_El.Init(&RunInfo::HLTElectronAllNames, "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v.*", {"hltEle50CaloIdVTGsfTrkIdTGsfDphiFilter"});
	ElJet_Jet.Init(&RunInfo::HLTJetAllNames, "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v.*", {"hltMonoPFJet165"});

	if(UserTime() == 2016)
	{
		MuIso.Init(&RunInfo::HLTMuonAllNames, "HLT_IsoMu24_v.*", {"hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09"});
		MuTkIso.Init(&RunInfo::HLTMuonAllNames, "HLT_IsoTkMu24_v.*", {"hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09"});
		ElIso27.Init(&RunInfo::HLTElectronAllNames, "HLT_Ele27_WPTight_Gsf_v.*", {"hltEle27WPTightGsfTrackIsoFilter"});
	}
	else
	{
		MuIso.Init(&RunInfo::HLTMuonAllNames, "HLT_IsoMu27_v.*", {"hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"});
		ElIso32.Init(&RunInfo::HLTElectronAllNames, "HLT_Ele32_WPTight_Gsf_v.*", {"hltEle32WPTightGsfTrackIsoFilter"});
		ElIso27.Init(&RunInfo::HLTElectronAllNames, "HLT_Ele27_WPTight_Gsf_v.*", {"hltEle27WPTightGsfTrackIsoFilter"});
	}
	Mu50.Init(&RunInfo::HLTMuonAllNames, "HLT_Mu50_v.*", {"hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"});

}

BoostedTop::~BoostedTop()
{
	histfile->Write();
	histfile->Close();
}




Int_t BoostedTop::AnalyseEvent()
{
	ALLHAD_ = false;
	DILEP_ = false;
	SEMILEP_ = false;
	GenMus.clear();
	GenEls.clear();
	gps.clear();
	HadTops.clear();
	SBasicParticles.clear();
	RecMusLoose.clear();
	RecMusVeto.clear();
	RecMus.clear();
	//SMuons.clear();
	SMuons.reset(NumIOMuons());
	RecElsLoose.clear();
	RecElsVeto.clear();
	RecEls.clear();
	//SElectrons.clear();
	SElectrons.reset(NumIOElectrons());
	RecPhs.clear();
	RecPrefirePhs.clear();
	//SPhotons.clear();
	SPhotons.reset(NumIOPhotons());
	//SGenAK8s.clear();
	SGenAK8s.reset(NumIOGenAK8Jets());
	GenAK8s.clear();
	//SAK4s.clear();
	SAK4s.reset(NumIOPFAK4Jets());
	RecAK4s.clear();
	RecPrefireAK4s.clear();
	//SAK8s.clear();
	SAK8s.reset(NumIOPFAK8Jets());
	RecAK8s.clear();
	RecHadJets.clear();
	SelectedHadJets.clear();
	SHadJets.clear();
	SGenHadJets.clear();
	GenHadJets.clear();
	RecLepJets.clear();
	SelectedLepJets.clear();
	SLepJets.clear();
	pstphotons.clear();
	GenAllJets.clear();
	GenBJets.clear();
	GenCJets.clear();
	GenLJets.clear();
	gent = nullptr;
	gentbar = nullptr;

	ran_BoostedTH_ = false;
	ran_BoostedTL_ = false;
	ran_SelectRECOParticles_ = false;
	ran_SelectGENParticles_ = false;
	ran_GenBoostedTH_ = false;

	ttgen.Reset();
	ttgenps.Reset();
	ttgenpa.Reset();
	ttgenlhe.Reset();
	tttruthres.Reset();
	tttruthboosted.Reset();

	//if(NumGoodVertices() < 30) {return 1;}
	run = Run();
	//if(!IsMC() && run <= 304671) {return 1;}
	//if(!IsMC() && run > 304671) {return 1;}
	rho = Rho();
	weight = 1.;
	puweight_ = 1.;
	//NTCweight = 1.;
	lepweight_ = 1.;
	prefireweight_ = 1.;
	btagweight_ = 1.;
	_hists1d["counter"]->Fill(0.5, weight);
	if(IsData() && CurrentLBInfo()->LumiValue() == -1.){cout << "NOLUMI " << Run() << " " << LumiBlock() << endl; return(1);}
	_hists1d["counter"]->Fill(1.5, weight);
	//if(CheckDuplicate() != 0) {cout << "DUPLICATE" << endl; return(1);}


	if(IsMC())
	{
		if(weighttype_ == 1 && GetGenInfo(0).Num_PSWeights() > 1) {weight = GetGenInfo(0).PSWeights(ttmcweightnum_);}
		else if(weighttype_ == 0 && GetGenInfo(0).Num_ScaleWeights() > 0) {weight = GetGenInfo(0).ScaleWeights(ttmcweightnum_)/abs(GetGenInfo(0).ScaleWeights(0));}
		else {weight = GetGenInfo(0).PSWeights(0);}

		defaultpdfweight_ = GetGenInfo(0).PSWeights(0);

		_hists1d["mcweight"]->Fill(weight);


//cout << "WEI " <<  weight << " " << GetGenInfo(0).PSWeights(0) << " " << GetGenInfo(0).ScaleWeights(0) << " " << GetGenInfo(0).NNPDFWeights(0) << endl;

	}

	_hists1d["counter"]->Fill(2.5, weight);
	if(IsPT650() && TTMC_) 
	{
		//cout << "HT650" << endl;
		_hists1d["mcinfo"]->Fill(0.5);
		_hists1d["mcinfo"]->Fill(1.5, weight);
		if(SKIPPT650_ && !TTPT650MC_) {return 1;}
	}
	else 
	{
		//cout << "NoHT650" << endl;
		_hists1d["mcinfo"]->Fill(2.5);
		_hists1d["mcinfo"]->Fill(3.5, weight);
	}

	_hists1d["counter"]->Fill(3.5, weight);

if(abs(weight) > 100 || weight != weight) {cout << "Extreme MC weight " << weight << endl; weight = 1.;}

	if(IsMC() && NumIOMETs() != 0 && UserTime() != 2023)
	{
		double npu = GetGenInfo(0).NumTrueInteractions();
		_hists1d["mu"]->Fill(npu, weight);
		puweight_ = puhist_->GetBinContent(puhist_->FindFixBin(npu));
		weight *= puweight_;
		_hists1d["puweight"]->Fill(puweight_);
		_hists1d["mu_weighted"]->Fill(npu, weight);
	}

	_hists1d["counter"]->Fill(4.5, weight);

	SelectGENParticles();
	if(DELPHES && gent != nullptr && gentbar != nullptr)
	{
		bool skip = false;
		//if(gent->DeltaR(pret) < 1E-5 && gentbar->DeltaR(pretb) < 1E-5) {cout << "LIKEPRE " << Run()*100000+LumiBlock() << endl; skip = true;}
		if(gent->DeltaR(pret) < 1E-5 && gentbar->DeltaR(pretb) < 1E-5) {skip = true;}
		pret=*gent;
		pretb=*gentbar;

		if(gent->P() == 0. && gentbar->P() == 0.) {skip = true;}

		if(skip) {return 1;}
	}


	if(TTMC_ && PS_MODE_)
	{
		if(BTAG_EFF_MODE_) { SelectPseudoTop(); }
		else { SelectPseudoTopBoosted(); }
	}
	FillGenJets();
	//if(TTCENTRAL && (*gent + *gentbar).M() > 700){return 1;}

	ttgen = (PS_MODE_ ? ttgenps : ttgenpa);
	if(TTMC_) {ttgen.Calculate();}

	GenBoostedTH();

	if(ttgen.IsComplete()){unfolding_all_.FillTruth(ttgen, weight);}
	if(ttgenpa.IsComplete()){unfolding_ps_.FillTruth(ttgenpa, weight);}
	if(ttgenps.IsComplete()){unfolding_ps_.FillReco(ttgenps, "RES", weight);}
	if(ttgenpa.IsComplete() && ttgenps.IsComplete()){unfolding_ps_.FillTruthReco(ttgenpa, ttgenps, "RES", weight);}

	if(NumIOMETs() == 0) {return 1;}

	if(UserTime() == 2023) {SelectRECOParticles2023();} else{SelectRECOParticles();}
	BoostedTL();
	BoostedTH();
	GenRecoMatching();

	bool muisotrg = false;
	bool muhighpttrg = false;
	bool elisotrg = false;
	bool elhighpttrg = false;

	if(UserTime() == 2016)
	{
		muisotrg = trgres["IsoMu24"].GetResult() == 1 || trgres["IsoTkMu24"].GetResult() == 1;
		muhighpttrg = trgres["Mu50"].GetResult() == 1;
		elisotrg = trgres["Ele27"].GetResult() == 1;
		elhighpttrg = trgres["ElJet"].GetResult() == 1;
	}
	else if(UserTime() == 2017)
	{
		muisotrg = trgres["IsoMu27"].GetResult() == 1;
		muhighpttrg = trgres["Mu50"].GetResult() == 1;
		elisotrg = trgres["Ele27"].GetResult() == 1 || trgres["Ele32"].GetResult() == 1;
		elhighpttrg = trgres["ElJet"].GetResult() == 1;

		if(IsMC() && !elisotrg && !muisotrg && !muhighpttrg && elhighpttrg) {weight *= 0.8845;}
	}
	else if(UserTime() == 2018)
	{
		muisotrg = trgres["IsoMu27"].GetResult() == 1;
		muhighpttrg = trgres["Mu50"].GetResult() == 1;
		elisotrg = trgres["Ele32"].GetResult() == 1;
		elhighpttrg = trgres["ElJet"].GetResult() == 1;
	}
	else if(UserTime() == 2023)
	{
		muisotrg = true;
		muhighpttrg = true;
		elisotrg = true;
		elhighpttrg = true;
	}

	bool mutrg = muisotrg || muhighpttrg;
	bool eltrg = elisotrg || elhighpttrg;

	bool mettest = (GetIOEventInfo(0).METFilterPassed() & metfiltertest_) == metfiltertest_;
	//if(UserTime() == 2016 && GetIOEventInfo(0).METFilterPassed() == 159) {mettest = true;} //temporary hack

	if(mettest && ( (IsMC() && (mutrg || eltrg)) || (DAMU_ && mutrg) || (DAEL_ && eltrg && !mutrg)))
	{

		_hists1d["counter"]->Fill(5.5, weight);

		if(IsMC())
		{
			if(UserTime() == 2017 || UserTime() == 2016)
			{
				prefireweight_ =  prefireprob.getProb();
				weight *= prefireweight_;
				_hists1d["prefireweight"]->Fill(prefireweight_);
			}

			if(UserTime() == 2017 || UserTime() == 2016) { lepweight_ = lepcorrector.correctionEvent2(RecEls, RecElsLoose,  RecMus, RecMusLoose, RecAK4s, csigmael_, csigmamu_); }
			else { lepweight_ = lepcorrector.correctionEvent(RecEls, RecMus, RecAK4s, csigmael_, csigmamu_); }
			_hists1d["lepweight"]->Fill(lepweight_);
			weight *= lepweight_;

			btagweight_ = btagweighter.SF(RecAK4s, GenBJets, GenCJets, weight);
			_hists1d["btagweight"]->Fill(btagweight_);
			if(!BTAG_EFF_MODE_) weight *= btagweight_;
		}

		if(mutrg && RecMus.size() == 1) { ElJetTriggerAnalysis();}
		PrefireAnalysis();

	
		if(!ISOSIDEBAND_ && (RecAK4s.size() >= 4 || RecHadJets.size() >= 1))
		{

if(weight != weight) {cout << "WEIGHT ALARM " << weight << " " << puweight_ << " " << lepweight_ << " " << btagweight_ << " " << prefireweight_ << endl;}
			if(RecMusLoose.size() + RecElsLoose.size() > 0 && RecMusVeto.size() + RecElsVeto.size() == 0)
			{
				if(RecMusLoose.size() > 0) _hists1d["counter"]->Fill(13.5, weight);
				if(RecElsLoose.size() > 0) _hists1d["counter"]->Fill(14.5, weight);
				AnalysisNonIsoLep();
			}
			else if(RecMus.size() + RecEls.size() == 1 && RecMusVeto.size() + RecElsVeto.size() == 1)
			{
//cout << elisotrg << " " << elhighpttrg  << " - " <<  NumIOElectrons() << " " << RecEls.size() << " j " << NumIOPFAK4Jets()  << " " << RecAK4s.size() << endl;
				if(RecMus.size() > 0) _hists1d["counter"]->Fill(15.5, weight);
				if(RecEls.size() > 0) _hists1d["counter"]->Fill(16.5, weight);
				AnalysisIsoLep();
			}
		}
	}

	if(ISOSIDEBAND_ && (RecAK4s.size() >= 4 || RecHadJets.size() >= 1))
	{
		if(RecMusLoose.size() == 1 && RecMusVeto.size() + RecElsVeto.size() == 0)
		{
			int isotrg = trgres["IsoMu27"].GetResult();
			int highpttrg = trgres["Mu50"].GetResult();
			int lowpttrg = trgres["Mu27"].GetResult();
			
			isotrg = isotrg <= 0 ? 1000 : isotrg;
			highpttrg = highpttrg <= 0 ? 1000 : highpttrg;
			lowpttrg = lowpttrg <= 0 ? 1000 : lowpttrg;
			int prescale = min({isotrg, highpttrg, lowpttrg});

			if(prescale < 1000)
			{
				weight *= prescale;
				OMuon* mu = RecMusLoose[0];
				RecMus.push_back(mu);
				for(OJet* jet : RecAK4s)
				{
					for(size_t nmu = 0 ; nmu < jet->Num_MatchedMuons() ; ++nmu)
					{
						if(jet->MatchedMuons(nmu) == mu->Num()) 
						{
							jet->SetPxPyPzE(jet->Px()-mu->Px(), jet->Py()-mu->Py(), jet->Pz()-mu->Pz(), jet->E()-mu->E());
						}
					}
				}
				AnalysisIsoLep();
				weight /= prescale;
			}
		}
	}

	AnalysisBKG(1);


	return(1);
}



void BoostedTop::FileChanged()
{

	//cout << "File Changed" << endl;
	//LoadPrimaryVertex(false);
	//LoadIOEventInfo(false);
	//LoadIOPFJet(false);
	//LoadGenInfo(false);
	//LoadIOMET(false);
	//LoadIOGenAK4Jet(false);
	//LoadSelectedGenParticle(false);
	//LoadAllGenParticle(false);
	//LoadIOElectron(false);
	//LoadIOTrack(false);
	//LoadIOBeamSpot(false);
	//LoadIOPhoton(false);
}


void BoostedTop::BeginLoop()
{
	histfile->cd();

	if(TTMC_)
	{
		if(FSRunc_ == 1) {weighttype_ = 1; ttmcweightnum_ = 3;}
		else if(FSRunc_ == -1) {weighttype_ = 1; ttmcweightnum_ = 5;}
		else if(FSRunc_ == -2) {weighttype_ = 1; ttmcweightnum_ = 9;}
		else if(FSRunc_ == 2) {weighttype_ = 1; ttmcweightnum_ = 7;}
		else if(ISRunc_ == 1) {weighttype_ = 1; ttmcweightnum_ = 6;}
		else if(ISRunc_ == -1) {weighttype_ = 1; ttmcweightnum_ = 8;}
		else if(RSunc_ == -1 && FSunc_ == -1) {weighttype_ = 0; ttmcweightnum_ = 8;}
		else if(RSunc_ == 0 && FSunc_ == -1) {weighttype_ = 0; ttmcweightnum_ = 2;}
		else if(RSunc_ == 1 && FSunc_ == -1) {weighttype_ = 0; ttmcweightnum_ = 5;}
		else if(RSunc_ == -1 && FSunc_ == 0) {weighttype_ = 0; ttmcweightnum_ = 6;}
		else if(RSunc_ == 0 && FSunc_ == 0) {weighttype_ = 0; ttmcweightnum_ = 0;}
		else if(RSunc_ == 1 && FSunc_ == 0) {weighttype_ = 0; ttmcweightnum_ = 3;}
		else if(RSunc_ == -1 && FSunc_ == 1) {weighttype_ = 0; ttmcweightnum_ = 7;}
		else if(RSunc_ == 0 && FSunc_ == 1) {weighttype_ = 0; ttmcweightnum_ = 1;}
		else if(RSunc_ == 1 && FSunc_ == 1) {weighttype_ = 0; ttmcweightnum_ = 4;}
	}

	cout << "MCNAME: " << mcname_ << endl;
	if(CP->Get<bool>("MCSpecificBTAG")){
		btagweighter.Init(CP->Get<string>("BTag_FileName"), CP->Get<string>("BTag_Efficiencies"), btag_loose_, btag_medium_, mcname_, CP->Get<double>("BTag_BUnc"), CP->Get<double>("BTag_LUnc"));
	}
	else
	{
		btagweighter.Init(CP->Get<string>("BTag_FileName"), CP->Get<string>("BTag_Efficiencies"), btag_loose_, btag_medium_, "P8", CP->Get<double>("BTag_BUnc"), CP->Get<double>("BTag_LUnc"));
	}

	if(DELPHES)
	{
		delphesjec.Init("HL_YR_JEC.root");
	}
	else
	{

		jetscalerPUPPI.Init(CP->Get<string>("JEC_UncertaintyFilePUPPI"), CP->Get<string>("JEC_UncB"), true);
		jetscalerPUPPI.Init(CP->Get<string>("JEC_UncertaintyFilePUPPI"), CP->Get<string>("JEC_UncLIGHT"), false);
		jetscalerPUPPI.InitResolution(CP->Get<string>("JEC_ResFilePUPPI"), CP->GetVector<string>("JEC_ResSFFilePUPPI"));
		//jetscalerPUPPI.InitMCrescale(mcname_, CP->Get<string>("JEC_MCSpecificFile"));

		jetscaler.Init(CP->Get<string>("JEC_UncertaintyFile"), CP->Get<string>("JEC_UncB"), true);
		jetscaler.Init(CP->Get<string>("JEC_UncertaintyFile"), CP->Get<string>("JEC_UncLIGHT"), false);
		jetscaler.InitResolution(CP->Get<string>("JEC_ResFile"), CP->GetVector<string>("JEC_ResSFFile"));
		if(CP->Get<bool>("MCSpecificJEC")) {jetscaler.InitMCrescale(mcname_, CP->Get<string>("JEC_MCSpecificFile"));}
	}
}

void BoostedTop::RunChanged()
{
	vector<string> mettests;

	if(IsMC() && UserTime() == 2016)
	{
		mettests = {"Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter"};
	}
	else if(!IsMC() && UserTime() == 2016)
	{
		mettests = {"Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter", "Flag_eeBadScFilter"};
	}
	else if(IsMC() && UserTime() >= 2017)
	{
		mettests = {"Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter", "Flag_ecalBadCalibFilter"};
	}
	else if(!IsMC() && UserTime() >= 2017)
	{
		mettests = {"Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter", "Flag_ecalBadCalibFilter", "Flag_eeBadScFilter"};
	}
	else
	{
		cout << "WARNING: no MET filter specified." << endl;
	}


	metfiltertest_ = 0;
	for(const string& filter : mettests)
	{
		for(size_t b = 0 ; b < CurrentRunInfo()->METFilterNames().size() ; ++b)
		{
			if(filter == CurrentRunInfo()->METFilterNames()[b])
			{
				metfiltertest_ |= 1<<b;
				break;
			}
			if(b == CurrentRunInfo()->METFilterNames().size()-1) {cout << "WARNING: could not find MET filter " << filter << endl;}
		}
	}


	//for(int n = 0 ; n < GetNumHLTriggers() ; ++n)
	//{
	//	cout << n << " " << GetHLTNames(n) << endl;
	//}


}

