#include "LepEffCorrection.h"

using namespace std;

void LepEffScales::init(const string& filename, const string& type)
{
	type_ = type;
	TDirectory* dir = gDirectory;
	efffile = TFile::Open(filename.c_str());
	heta = dynamic_cast<TH1D*>(efffile->Get(("eta_"+type).c_str()));
	for(int i = 0 ; i < heta->GetNbinsX() ; ++i)
	{
		effhists_mc.push_back(dynamic_cast<TH1D*>(efffile->Get(("MC_"+type+"_"+to_string(i)).c_str())));
		effhists_da.push_back(dynamic_cast<TH1D*>(efffile->Get(("DA_"+type+"_"+to_string(i)).c_str())));
	}
	dir->cd();
}

const vector<double>& LepEffScales::getvals(TLorentzVector* vec)
{
	double eta = vec->Eta();
	OElectron* el = dynamic_cast<OElectron*>(vec);
	if(el != nullptr) {eta = SuperCluster(el->SC()).Eta();}

	int etabin = heta->FindFixBin(eta)-1;
	double pt = min(vec->Pt(), effhists_mc[etabin]->GetXaxis()->GetXmax()-0.01);
	pt = max(pt, effhists_mc[etabin]->GetXaxis()->GetXmin()+0.01);
	int ptbin = effhists_mc[etabin]->FindFixBin(pt);
	//results[0] = effhists_mc[etabin]->Interpolate(pt);
	//results[1] = effhists_da[etabin]->Interpolate(pt);
	results[0] = effhists_mc[etabin]->GetBinContent(ptbin);
	results[1] = effhists_da[etabin]->GetBinContent(ptbin);
	results[2] = effhists_da[etabin]->GetBinError(ptbin);

	return results;

}


////////////LepEffCorrection//////////////////////////
void LepEffCorrection::init(const string& mufile, const string& elfile)
{
	for(const string& type : {"ElReco_ElTOT", "ElISOTRG", "ElORTRG"})
	{
		corrections_[type] = LepEffScales();
		corrections_[type].init(elfile, type);
	}
	for(const string& type : {"MuTRK_MuTOT", "MuISOTRG", "MuORTRG"})
	{
		corrections_[type] = LepEffScales();
		corrections_[type].init(mufile, type);
	}
}

void LepEffCorrection::init2(const string& mufile, const string& elfile)
{
	for(const string& type : {"ElReco_ElTOT", "ElReco_ElID", "ElISOTRG", "treliso", "trel", "trisoel", "trjet", "trmu"})
	{
		corrections_[type] = LepEffScales();
		corrections_[type].init(elfile, type);
	}
	for(const string& type : {"MuTRK_MuTOT", "MuTRK_MuID", "MuORTRG"})
	{
		corrections_[type] = LepEffScales();
		corrections_[type].init(mufile, type);
	}
}


double LepEffCorrection::correctionEvent(vector<OElectron*> els, vector<OMuon*> mus, vector<OJet*> jets, double sigmael, double sigmamu)
{
	if(mus.size() == 0 && els.size() == 0) {return 1.;}
	double recoeffmc = 1.;
	double trgeffmc = 1.;
	double recoeffda = 1.;
	double trgeffda = 1.;
	if(els.size() > 0)
	{
		double maxjetpt = 0.;
		for(OJet* j : jets) {maxjetpt = max(maxjetpt, j->Pt());}

		string trgname = "ElISOTRG";
		if(maxjetpt > 165.) {trgname = "ElORTRG";}

		for(OElectron* el : els)
		{
			const vector<double>& totres = corrections_["ElReco_ElTOT"].getvals(el);
			recoeffmc *= totres[0];
			recoeffda *= totres[1]+sigmael*totres[2];

			const vector<double>& trgres = corrections_[trgname].getvals(el);
			trgeffmc *= 1.-trgres[0];
			trgeffda *= 1.-trgres[1]-sigmael*trgres[2];

		}
	}
	if(mus.size() > 0)
	{
		for(OMuon* mu : mus)
		{
			const vector<double>& totres = corrections_["MuTRK_MuTOT"].getvals(mu);
			recoeffmc *= totres[0];
			recoeffda *= totres[1]+sigmamu*totres[2];

			const vector<double>& trgres = corrections_["MuORTRG"].getvals(mu);
			trgeffmc *= 1.-trgres[0];
			trgeffda *= 1.-trgres[1]-sigmamu*trgres[2];

		}
	}
	return recoeffda/recoeffmc * (1.-trgeffda)/(1.-trgeffmc);
}

double LepEffCorrection::correctionEvent2(vector<OElectron*> els, vector<OElectron*> looseels, vector<OMuon*> mus, vector<OMuon*> loosemus, vector<OJet*> jets, double sigmael, double sigmamu)
{
	if(mus.size()+loosemus.size() == 0 && els.size()+looseels.size() == 0) {return 1.;}
	double recoeffmc = 1.;
	double trgeffmc = 1.;
	double recoeffda = 1.;
	double trgeffda = 1.;
	if(els.size()+looseels.size() > 0)
	{
		double jetprobda = 1.;
		double jetprobmc = 1.;
		for(OJet* jet : jets)
		{
			if(jet->Pt() < 155.) {continue;}
			const vector<double>& jetres = corrections_["trjet"].getvals(jet);
			jetprobmc *= 1.-jetres[0];
			jetprobda *= 1.-jetres[1]-sigmael*jetres[2];
		}
		jetprobda = 1.-jetprobda;
		jetprobmc = 1.-jetprobmc;

		for(OElectron* el : els)
		{
			const vector<double>& totres = corrections_["ElReco_ElTOT"].getvals(el);
			recoeffmc *= totres[0];
			recoeffda *= totres[1]+sigmael*totres[2];

			const vector<double>& trgor = corrections_["treliso"].getvals(el);
			const vector<double>& trgiso = corrections_["ElISOTRG"].getvals(el);
			trgeffmc *= 1.-(jetprobmc*trgor[0] + (1.-jetprobmc)*trgiso[0]);
			trgeffda *= 1.-(jetprobda*(trgor[1]+sigmael*trgor[2]) + (1.-jetprobda)*(trgiso[1]+sigmael*trgor[2]));

		}
		for(OElectron* el : looseels)
		{
			if(els.size() == 0)
			{
				const vector<double>& totres = corrections_["ElReco_ElID"].getvals(el);
				recoeffmc *= totres[0];
				recoeffda *= totres[1]+sigmael*totres[2];
			}
			const vector<double>& trgor = corrections_["trel"].getvals(el);
			const vector<double>& trgiso = corrections_["trisoel"].getvals(el);
			trgeffmc *= 1.-(jetprobmc*trgor[0] + (1.-jetprobmc)*trgiso[0]);
			trgeffda *= 1.-(jetprobda*(trgor[1]+sigmael*trgor[2]) + (1.-jetprobda)*(trgiso[1]+sigmael*trgor[2]));

		}
	}

	for(OMuon* mu : mus)
	{
		const vector<double>& totres = corrections_["MuTRK_MuTOT"].getvals(mu);
		recoeffmc *= totres[0];
		recoeffda *= totres[1]+sigmamu*totres[2];

		const vector<double>& trgres = corrections_["MuORTRG"].getvals(mu);
		trgeffmc *= 1.-trgres[0];
		trgeffda *= 1.-trgres[1]-sigmamu*trgres[2];

	}
	for(OMuon* mu : loosemus)
	{
		if(mus.size() == 0)
		{
			const vector<double>& totres = corrections_["MuTRK_MuID"].getvals(mu);
			recoeffmc *= totres[0];
			recoeffda *= totres[1]+sigmamu*totres[2];
		}
		const vector<double>& trgres = corrections_["trmu"].getvals(mu);
		trgeffmc *= 1.-trgres[0];
		trgeffda *= 1.-trgres[1]-sigmamu*trgres[2];

	}
	return recoeffda/recoeffmc * (1.-trgeffda)/(1.-trgeffmc);
}

