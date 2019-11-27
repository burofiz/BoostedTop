#include <ConfigParser.h>
#include <Analyse.h>
#include <Run.h>
#include <OPhoton.h>
#include <OJet.h>
#include <LV.h>

#include <regex>

#include "PrefireProb.h"
#include "BoostedTop.h"

using namespace std;

void L1Prefire::UpdateRun()
{
	if(GLAN->IsMC()) {return;}
	RunInfo* run = GLAN->CurrentRunInfo();

	regex reg_isoEGer("L1_SingleIsoEG\\d*er(2p1)*");
	regex reg_isoEG("L1_SingleIsoEG\\d*(er2p5)*");
	regex reg_EGer("L1_SingleEG\\d*er(2p1)*");
	regex reg_EG("L1_SingleEG\\d*(er2p5)*");

	
	in_EG.clear();
	in_EG.emplace_back(0, -1);
	for(size_t l1 = 0 ; l1 < run->NumL1() ; ++l1)
	{
		string ls = run->L1Name(l1);
		if(regex_match(ls, reg_EG))
		{
			vector<string> sdata = string_split(ls, {"EG"});
			int pt = stof(sdata[1]);
			//cout << ls << " " << pt << endl;
			in_EG.emplace_back(pt, l1);
		}

	}
	if(in_EG.size() > 0) {in_EG.emplace_back(10000, in_EG.back().second);}
	sort(in_EG.begin(), in_EG.end(), [](auto& a, auto& b){return a.first < b.first;});

	in_EGer.clear();
	in_EGer.emplace_back(0, -1);
	for(size_t l1 = 0 ; l1 < run->NumL1() ; ++l1)
	{
		string ls = run->L1Name(l1);
		if(regex_match(ls, reg_EGer))
		{
			vector<string> sdata = string_split(ls, {"EG"});
			int pt = stof(sdata[1]);
			//cout << ls << " " << pt << endl;
			in_EGer.emplace_back(pt, l1);
		}
	}
	if(in_EGer.size() > 0) {in_EGer.emplace_back(10000, in_EGer.back().second);}
	sort(in_EGer.begin(), in_EGer.end(), [](auto& a, auto& b){return a.first < b.first;});

	in_isoEG.clear();
	in_isoEG.emplace_back(0, -1);
	for(size_t l1 = 0 ; l1 < run->NumL1() ; ++l1)
	{
		string ls = run->L1Name(l1);
		if(regex_match(ls, reg_isoEG))
		{
			vector<string> sdata = string_split(ls, {"EG"});
			int pt = stof(sdata[1]);
			//cout << ls << " " << pt << endl;
			in_isoEG.emplace_back(pt, l1);
		}

	}
	if(in_isoEG.size() > 0) {in_isoEG.emplace_back(10000, in_isoEG.back().second);}
	sort(in_isoEG.begin(), in_isoEG.end(), [](auto& a, auto& b){return a.first < b.first;});

	in_isoEGer.clear();
	in_isoEGer.emplace_back(0, -1);
	for(size_t l1 = 0 ; l1 < run->NumL1() ; ++l1)
	{
		string ls = run->L1Name(l1);
		if(regex_match(ls, reg_isoEGer))
		{
			vector<string> sdata = string_split(ls, {"EG"});
			int pt = stof(sdata[1]);
			//cout << ls << " " << pt << endl;
			in_isoEGer.emplace_back(pt, l1);
		}

	}
	if(in_isoEGer.size() > 0) {in_isoEGer.emplace_back(10000, in_isoEGer.back().second);}
	sort(in_isoEGer.begin(), in_isoEGer.end(), [](auto& a, auto& b){return a.first < b.first;});
}

void L1Prefire::UpdateLB()
{
	if(GLAN->IsMC()) {return;}
	Int_t pscolumn  = GLAN->CurrentLBInfo()->HLTTable();
	const double psmax = 1.E100;
	double ps = psmax;
	for(size_t n = 0 ; n < in_EG.size()-1 ; ++n)
	{
		Bin ptbin(in_EG[n].first, in_EG[n+1].first);
		if(in_EG[n].second == -1){ps = psmax;}		
		else { double psnew = GLAN->CurrentRunInfo()->L1Prescale(in_EG[n].second, pscolumn);
		if(psnew != 0) {ps = psnew;}}
		//cout << "aEG: " << in_EG[n].first << " " <<  in_EG[n+1].first << " " << ps << endl;
		ps_EG[ptbin] = ps;
	}

	ps = psmax;
	for(size_t n = 0 ; n < in_EGer.size()-1 ; ++n)
	{
		Bin ptbin(in_EGer[n].first, in_EGer[n+1].first);
		if(in_EGer[n].second == -1){ps = psmax;}		
		else { double psnew = GLAN->CurrentRunInfo()->L1Prescale(in_EGer[n].second, pscolumn);
		if(psnew != 0) {ps = psnew;}}
		//cout << "bEG: " << in_EGer[n].first << " " <<  in_EGer[n+1].first << " " << ps << endl;
		ps_EGer[ptbin] = ps;
	}

	ps = psmax;
	for(size_t n = 0 ; n < in_isoEG.size()-1 ; ++n)
	{
		Bin ptbin(in_isoEG[n].first, in_isoEG[n+1].first);
		if(in_isoEG[n].second == -1){ps = psmax;}		
		else { double psnew = GLAN->CurrentRunInfo()->L1Prescale(in_isoEG[n].second, pscolumn);
		if(psnew != 0) {ps = psnew;}}
		//cout << "cEG: " << in_isoEG[n].first << " " <<  in_isoEG[n+1].first << " " << ps << endl;
		ps_isoEG[ptbin] = ps;
	}

	ps = psmax;
	for(size_t n = 0 ; n < in_isoEGer.size()-1 ; ++n)
	{
		Bin ptbin(in_isoEGer[n].first, in_isoEGer[n+1].first);
		if(in_isoEGer[n].second == -1){ps = psmax;}		
		else { double psnew = GLAN->CurrentRunInfo()->L1Prescale(in_isoEGer[n].second, pscolumn);
		if(psnew != 0) {ps = psnew;}}
		//cout << "dEG: " << in_isoEGer[n].first << " " <<  in_isoEGer[n+1].first << " " << ps << endl;
		ps_isoEGer[ptbin] = ps;
	}
}

double L1Prefire::FireProb(const L1Object& l1eg)
{
	double eta = l1eg.Eta();
	double pt = l1eg.Pt();
	bool iso = (l1eg.Iso() & 1); //tight iso bit 
	double psmin = 1.E100;
	if(abs(eta) <= 2.1)
	{
		psmin = min(psmin, ps_EGer[pt]);
		if(iso) {psmin = min(psmin, ps_isoEGer[pt]);}	
	}
	else
	{
		psmin = min(psmin,ps_EG[pt]);
		if(iso) {psmin = min(psmin,ps_isoEG[pt]);}	
	}
//cout << "PSMIN: " << psmin << " " << ps_EG[pt] << " " << ps_EGer[pt]<< " " << ps_isoEG[pt]<< " " << ps_isoEGer[pt]<< endl;
	return 1./psmin;
}

void PrefireProb::Init(const string& prefiremapfile, double sigma)
{
	TDirectory* cdir = gDirectory;

	TFile* tf = TFile::Open(prefiremapfile.c_str());
	prefiremapph = dynamic_cast<TH2D*>(tf->Get("prefiremap_ph"));
	prefiremapjet = dynamic_cast<TH2D*>(tf->Get("prefiremap_jet"));

	for(int x = 1 ; x <= prefiremapph->GetNbinsX() ; ++x)
	{
		for(int y = 1 ; y <= prefiremapph->GetNbinsY() ; ++y)
		{
			prefiremapph->SetBinContent(x,y, prefiremapph->GetBinContent(x,y)+sigma*prefiremapph->GetBinError(x,y));
		}
	}

	for(int x = 1 ; x <= prefiremapjet->GetNbinsX() ; ++x)
	{
		for(int y = 1 ; y <= prefiremapjet->GetNbinsY() ; ++y)
		{
			prefiremapjet->SetBinContent(x,y, prefiremapjet->GetBinContent(x,y)+sigma*prefiremapjet->GetBinError(x,y));
		}
	}


	cdir->cd();

	ConfigParser CP("config.cfg");

	ps_EG = CP.GetVector<double>("ps_EG");
	ps_isoEG = CP.GetVector<double>("ps_isoEG");
	ps_EGer = CP.GetVector<double>("ps_EGer");
	ps_isoEGer = CP.GetVector<double>("ps_isoEGer");

}

PrefireProb::~PrefireProb()
{
	if(tf != nullptr) {tf->Close();}
}




double PrefireProb::getProb() const
{

	double preprob = 1.;

	for(OPhoton* ph : dynamic_cast<BoostedTop*>(GLAN)->RecPrefirePhs)
	{
		preprob *= 1.-prefiremapph->Interpolate(ph->Eta(), min(150., ph->Pt()));
	}
	for(OJet* jet : dynamic_cast<BoostedTop*>(GLAN)->RecPrefireAK4s)
	{
		preprob *= 1.-prefiremapjet->Interpolate(jet->Eta(), min(250., jet->Pt()));
	}

	return preprob;

}


double PrefireProb::FireProb(const L1Object& l1eg) const
{
	double eta = l1eg.Eta();
	double pt = l1eg.Pt();
	int bpt = min(49, (int)pt);
	bool iso = (l1eg.Iso() & 1); //tight iso bit 
	double fireprob = 0.;
	if(abs(eta) <= 2.1)
	{
		fireprob = max(fireprob, ps_EGer[bpt]);
		if(iso) {fireprob = max(fireprob, ps_isoEGer[bpt]);}	
	}
	else
	{
		fireprob = max(fireprob,ps_EG[bpt]);
		if(iso) {fireprob = max(fireprob,ps_isoEG[bpt]);}	
	}
	return fireprob;
}

