#ifndef BOOSTEDTOP
#define BOOSTEDTOP
#include <Analyse.h>
#include <RLObject.h>
#include <OMuon.h>
#include <OElectron.h>
#include <OPhoton.h>
#include <OJet.h>
#include <Track.h>
#include <helper.h>
#include <NN.h>
#include <BTagWeight.h> 
#include <JetScaler.h>
#include <GenBasicParticle.h>
#include <GenSelectedParticle.h>
#include <GenSimpleParticle.h>
#include <GenJet.h>
#include <BHadronDecayWeights.h>

#include <vector>
#include <list>
#include <set>
#include <unordered_map>

#include "TTEvent.h" 
#include "TTPlots.h" 
#include "TopJets.h" 
#include "LepEffCorrection.h" 
#include "Unfolding.h" 
#include "BtagEff.h" 
#include "NNJet.h" 
//#include "THResidualCorrector.h"
//#include "PrefireCheck.h"
#include "PrefireProb.h"
#include "DelphesJEC.h"

class TFile;
class ConfigParser;
class BAnalysis;

using namespace std;


class BoostedTop : public Analyse
{
	private:
		bool ALLHAD_ = false;
		bool DILEP_ = false;
		bool SEMILEP_ = false;

		Unfolding unfolding_all_;
		Unfolding unfolding_ps_;

		bool ran_BoostedTH_ = false;
		bool ran_BoostedTL_ = false;
		bool ran_SelectRECOParticles_ = false;
		bool ran_SelectGENParticles_ = false;
		bool ran_GenBoostedTH_ = false;

		TrgPathMatch ElJet_El;
		TrgPathMatch ElJet_Jet;
		TrgPathMatch ElIso32;
		TrgPathMatch ElIso27;
		TrgPathMatch MuIso;
		TrgPathMatch MuTkIso;
		TrgPathMatch Mu50;
	public:
		bool IsUnPrefireable() { return GetIOEventInfo(0).METFilterPassed() & (1<<CurrentRunInfo()->METFilterNames().size());}
		bool IsPT650() { return GetIOEventInfo(0).METFilterPassed() & (1<<(CurrentRunInfo()->METFilterNames().size()+1));}
		TFile* histfile;

		BtagEff btageff;

		TTEvent ttgen;
		TTEvent ttgenps;
		TTEvent ttgenpa;
		TTEvent ttgenlhe;
		TTEvent tttruthres;
		TTEvent tttruthboosted;
		TLorentzVector met;

		TLorentzVector nu;
		vector<GenSimpleParticle> gps; 
		list<GenSimpleParticle> SBasicParticles; 
		list<GenSimpleParticle*> HadTops; 
		vector<GenSimpleParticle*> GenMus; 
		vector<GenSimpleParticle*> GenEls; 
		vector<GenSimpleParticle*> GenAllJets; 
		vector<GenSimpleParticle*> GenBJets; 
		vector<GenSimpleParticle*> GenCJets; 
		vector<GenSimpleParticle*> GenLJets; 
		GenSimpleParticle* gent;
		GenSimpleParticle* gentbar;
		vector<TLorentzVector*> pstphotons;

		ocontainer<OMuon> SMuons;
		vector<OMuon*> RecMusLoose;
		vector<OMuon*> RecMusVeto;
		vector<OMuon*> RecMus;
		ocontainer<OElectron> SElectrons;
		vector<OElectron*> RecElsLoose;
		vector<OElectron*> RecElsVeto;
		vector<OElectron*> RecEls;
		ocontainer<OPhoton> SPhotons;
		vector<OPhoton*> RecPhs;
		vector<OPhoton*> RecPrefirePhs;
		ocontainer<OJet> SAK8s;
		vector<OJet*> RecAK8s;
		ocontainer<OJet> SAK4s;
		vector<OJet*> RecAK4s;
		vector<OJet*> RecPrefireAK4s;
		ocontainer<GenJet> SGenAK8s;
		vector<GenJet*> GenAK8s;

		list<HadTopJet> SHadJets;
		list<HadTopJet> SHadJetsRight;
		vector<HadTopJet*> RecHadJets;
		vector<HadTopJet*> SelectedHadJets;
		list<LepTopJet> SLepJets;
		vector<LepTopJet*> RecLepJets;
		vector<LepTopJet*> SelectedLepJets;
		list<HadTopJet> SGenHadJets;
		vector<HadTopJet*> GenHadJets;

		void SelectGENParticles();
		void FillGenJets();
		void SelectPseudoTop();
		void SelectPseudoTopBoosted();
		void SelectRECOParticles2023();
		void SelectRECOParticles();
		void GenRecoMatching();
		void AnalyseResolved(TTEvent& ttrec);
		bool BoostedTL();
		void BoostedTH();
		void GenBoostedTH();
		int isResolved(vector<GenSimpleParticle*>& decays);
		bool AnalysisIsoLep();
		bool AnalysisNonIsoLep();
		bool AnalysisBKG(int prescale);
		
		void ElJetTriggerAnalysis();
		void PrefireAnalysis();
		
		double rho;
		int run;
		string mcname_ = "P8";

		double drmaxtop;
		double drmaxantitop;

		vector<string> isomuvec;
		vector<string> dmfvec;
		vector<string> smfvec;

		TTPlots hgenPS;
		TTPlots hgenAll;
		TTPlots hgenLHE;

		TTPlots hgenHadAll;
		TTPlots hgenHadRes;
		TTPlots hgenHadBoosted;
		TTPlots hgenHadBoostedMerged;
		TTPlots hgenHadBroken;
		TTPlots hgenHadBrokenRes;
		TTPlots hgenHadNo;
		TTPlots hgenLepRes;
		TTPlots hgenLepBoosted;
		TTPlots hgenLepNo;

		TTPlots hgenRES;
		TTPlots hgenOther;
		TTPlots hgenBH;
		TTPlots hgenBL;
		TTPlots hgenBHBL;

		TTPlots hrecRESright;
		TTPlots hrecRESwrong;
		TTPlots hrecRESnonreco;
		TTPlots hrecRESother;
		TTPlots hrecBHright;
		TTPlots hrecBHwrong;
		TTPlots hrecBHnonreco;
		TTPlots hrecBHother;
		TTPlots hrecBLright;
		TTPlots hrecBLwrong;
		TTPlots hrecBLnonreco;
		TTPlots hrecBLother;
		TTPlots hrecBHBLright;
		TTPlots hrecBHBLwrong;
		TTPlots hrecBHBLnonreco;
		TTPlots hrecBHBLother;

		LepJetHists lepjetmu_right;
		LepJetHists lepjetel_right;
		LepJetHists lepjetmu_wrong;
		LepJetHists lepjetel_wrong;
		HadJetHists hadjet_right;
		HadJetHists hadjet_wrong;
		HadJetHists hadjet_other;
		HadJetHists genhadjet_right;
		HadJetHists genhadjet_wrong;
		HadJetHists genhadjet_other;
		HadJetHists hadjetisolep_right;
		HadJetHists hadjetisolep_wrong;
		HadJetHists hadjetisolep_other;
		HadJetHists hadjetnonisolep_right;
		HadJetHists hadjetnonisolep_wrong;
		HadJetHists hadjetnonisolep_other;
		HadJetHists hadjetbkgPH;
		HadJetHists hadjetbkgB;
		bool promptph;
		bool QCDMC_ = false;

		TH1DCollection _hists1d;
		TH2DCollection _hists2d;

		//JetPlots H_alljets;
		//JetPlots H_bjets;
		//JetPlots H_cjets;
		//JetPlots H_wjets;
		//JetPlots H_ljets;

		int topdatasize = 180;
		float phpt, pheta;
		float* topdata = nullptr;
		TTree* topmutreesig = nullptr;
		TTree* topmutreebkg = nullptr;
		TTree* topeltreesig = nullptr;
		TTree* topeltreebkg = nullptr;

		int tophaddatasize = 180;
		float* tophaddata = nullptr;
		TTree* tophadtreesig = nullptr;
		TTree* tophadtreesigbroken = nullptr;
		TTree* tophadtreebkg = nullptr;
		
		NeuralNet* netmu1 = nullptr;
		NeuralNet* netmu2 = nullptr;
		NeuralNet* nethad1 = nullptr;
		NeuralNet* nethad2 = nullptr;
		NeuralNet* nethad3 = nullptr;
		NeuralNet* nethad4 = nullptr;
		//NeuralNet* nethad5 = nullptr;

		bool DELPHES = false;
		bool PDFUNCS_ = false;
		bool TRAIN = false;
		double cuncjer_ = 0.;
		double cuncjes_ = 0.;

		unique_ptr<ConfigParser> CP;

		bool BTAG_EFF_MODE_;
		bool PS_MODE_;
		bool ISOSIDEBAND_;
		bool SKIPPT650_ = false;
		int FSRunc_ = 0;
		int ISRunc_ = 0;
		int FSunc_ = 0;
		int RSunc_ = 0;

		BTagWeight btagweighter; 
		double btag_loose_;
		double btag_medium_;
		double btag_tight_;

		JetScaler jetscaler;
		JetScaler jetscalerPUPPI;

		BHadronDecayWeights bdecayweights;
        BFragWeights bfragweights;

		TH1D* puhist_ = nullptr;

        LepEffCorrection lepcorrector;
		double csigmamu_ = 1.;
		double csigmael_ = 1.;
		BoostedTop(string outfile);
		virtual ~BoostedTop();
		virtual Int_t AnalyseEvent();
		virtual void FileChanged();
		virtual void RunChanged();
		virtual void BeginLoop();

		bool TTCENTRAL = false;
		bool TTMC_ = false;
		bool TTPT650MC_ = false;
		bool DAMU_ = false;
		bool DAEL_ = false;
		int ttmcweightnum_ = 0;
		int weighttype_ = 0;
		bool first = true;
		double PuWeight() const {return puweight_;}
		const string& MCName() const {return mcname_;}

		vector<double> bins_toppt_large;
		vector<double> bins_toppt;
		vector<double> bins_topptsoft;
		vector<double> bins_topy;
		vector<double> bins_ttm;
		vector<double> bins_ttpt;
		vector<double> bins_tty;
		vector<double> bins_cts;


		TH1D* h_hadnnres = nullptr;
		TH2D* h_tmwm = nullptr;
		TH1D* h_dnu = nullptr;

		double cbsidebandmin_;
		double cbsidebandmax_;

int cmu = 0;
int cel = 0;

	int cfillsigtree_;
	int cfillbkgtree_;
	TTree* ttbartree;
	TTree* wbkgtree;
	TTree* qcdbkgtree;
	float thadptref;
	float thadetaref;
	float tleppt;
	float thadpt;
	float toppt;
	float topbarpt;
	float topy;
	float topbary;
	float thardpt;
	float tsoftpt;
	float thady;
	float tlepy;
	float ttm;
	float ttpt;
	float tty;
	float cts;
	float thadpt_ttm_tty;
	float thadpt_ttm_cts;
	float nnres;
	float btag;
	int type;
	int nvtx;

	double ctopptweight_;
	double cttmweight_;

	TLorentzVector pret;
	TLorentzVector pretb;

	PrefireProb prefireprob;

	double weight = 1.;	
	double puweight_ = 1.;	
	double lepweight_ = 1.;
	double prefireweight_ = 1.;
	double btagweight_ = 1.;
	double defaultpdfweight_ = 1.;

	DelphesJEC delphesjec;

	uint32_t metfiltertest_;

	unordered_map<string, TrgResult> trgres;

	L1Prefire prefiretest_;
	unique_ptr<BAnalysis> BAn;
};

#endif
