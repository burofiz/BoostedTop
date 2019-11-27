#include "Unfolding.h"
#include "TTEvent.h"
#include "BoostedTop.h"


vector<double> getnvec(int n)
{
	vector<double> vec(n+1);
	double val = 0.;
	for(double& v : vec) {v = val; val++;}
	return vec;
}

void Unfolding::Init(bool PDFUNCS)
{
	PDFUNCS_ = PDFUNCS;
	AN = dynamic_cast<BoostedTop*>(GLAN);
	const vector<double>& bins_toppt_large = AN->bins_toppt_large;
	const vector<double>& bins_toppt = AN->bins_toppt;
	const vector<double>& bins_topptsoft = AN->bins_topptsoft;
	const vector<double>& bins_topy = AN->bins_topy;
	const vector<double>& bins_ttm = AN->bins_ttm;
	const vector<double>& bins_ttpt = AN->bins_ttpt;
	const vector<double>& bins_tty = AN->bins_tty;
	const vector<double>& bins_cts = AN->bins_cts;
	vector<string> categories = {"RES", "BH", "BL", "BHBL"};

	if(GLAN->UserTime() == 2023)
	{
		vector<double> jetbins_large = {-0.5, 0.5, 1.5, 2.5, 3.5};
		vector<double> topptbins_large = {0.0, 90.0, 180.0, 270.0, 800.0};
		vector<double> topybins_large = {0.0, 0.5, 1.0, 1.5, 3.5};
		//vector<double> topybins_large = {0.0, 0.5, 1.0, 1.5, 2.5};
		vector<double> ttmbins_large = {300.0, 450.0, 625.0, 850.0, 2000.0};
		//vector<double> ttybins_large = {0.0, 0.4, 0.8, 1.2, 2.3};
		vector<double> ttybins_large = {0.0, 0.5, 1.0, 1.5, 3.0};
		vector<double> ttptbins_large = {0.0, 35.0, 80.0, 140.0, 500.0};

		unfolding1d_.AddMatrix("thadpt", bins_toppt, bins_toppt, categories);
		unfolding1d_.AddMatrix("thady", bins_topy, bins_topy, categories);
		unfolding1d_.AddMatrix("tleppt", bins_toppt, bins_toppt, categories);
		unfolding1d_.AddMatrix("tlepy", bins_topy, bins_topy, categories);
		unfolding1d_.AddMatrix("ttm", bins_ttm, bins_ttm, categories);
		unfolding1d_.AddMatrix("ttpt", bins_ttpt, bins_ttpt, categories);
		unfolding1d_.AddMatrix("tty", bins_tty, bins_tty, categories);

		vector< vector<double> > bins_thady_thadpt = {
			{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
			{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
			{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0},
			{0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 330.0, 380.0, 450.0, 800.0}
		};
		unfolding2d_.AddMatrix("thady+thadpt", topybins_large, bins_thady_thadpt, categories);

		vector< vector<double> > bins_ttm_tty = {
			{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0},
			{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.0},
			{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.0},
			{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 3.0}
		};
		unfolding2d_.AddMatrix("ttm+tty", ttmbins_large, bins_ttm_tty, categories);


	}
	else
	{
		vector<double> bins_toppt_unfolding = {0., 350., 7000.};
		unfolding2d_.AddMatrix("thadpt", {0, 7000.}, {bins_toppt}, categories);
		unfolding2d_.AddMatrix("tleppt", bins_toppt_unfolding, {bins_toppt, bins_toppt}, categories);
		unfolding2d_.AddMatrix("thardpt", bins_toppt_unfolding, {bins_toppt, bins_toppt}, categories);
		unfolding2d_.AddMatrix("tsoftpt", bins_toppt_unfolding, {bins_topptsoft, bins_topptsoft}, categories);
		unfolding2d_.AddMatrix("toppt", bins_toppt_unfolding, {bins_toppt, bins_toppt}, categories);
		unfolding2d_.AddMatrix("topy", bins_toppt_large, {bins_topy, bins_topy, bins_topy, bins_topy}, categories);
		unfolding2d_.AddMatrix("topbarpt", bins_toppt_unfolding, {bins_toppt, bins_toppt}, categories);
		unfolding2d_.AddMatrix("topbary", bins_toppt_large, {bins_topy, bins_topy, bins_topy, bins_topy}, categories);
		unfolding2d_.AddMatrix("thady", bins_toppt_large, {bins_topy, bins_topy, bins_topy, bins_topy}, categories);
		unfolding2d_.AddMatrix("tlepy", bins_toppt_large, {bins_topy, bins_topy, bins_topy, bins_topy}, categories);
		unfolding2d_.AddMatrix("ttm", bins_toppt_unfolding, {bins_ttm, bins_ttm}, categories);
		unfolding2d_.AddMatrix("ttpt", bins_toppt_unfolding, {bins_ttpt, bins_ttpt}, categories);
		unfolding2d_.AddMatrix("tty", bins_toppt_large, {bins_tty, bins_tty, bins_tty, bins_tty}, categories);
		unfolding2d_.AddMatrix("cts", bins_toppt_large, {bins_cts, bins_cts, bins_cts, bins_cts}, categories);

		vector<double> jetbins_large = {-0.5, 0.5, 1.5, 2.5, 3.5};
		vector<double> topptbins_large = {0.0, 90.0, 180.0, 270.0, 6500.0};
		vector<double> topybins_large = {0.0, 0.5, 1.0, 1.5, 2.5};
		//vector<double> ttmbins_large = {300.0, 450.0, 625.0, 850.0, 2000.0};
		vector<double> ttmbins_large = {300.0, 500.0, 700.0, 900.0, 13000.0};
		vector<double> ttybins_large = {0.3, 0.6, 0.9, 1.2, 1.5, 2.5};
		vector<double> ctsbins_large = {-1, -0.6, -0.2, 0.2, 0.6, 1.};
		vector<double> ttptbins_large = {0.0, 35.0, 80.0, 140.0, 500.0};
		unfolding2d_.AddMatrix("ttm+tty", ttmbins_large, {ttybins_large, ttybins_large, ttybins_large, ttybins_large}, categories);
		vector<double> nvecttm_tty = getnvec(unfolding2d_.GetNBins("ttm+tty"));
		unfolding2d_.AddMatrix("thadpt+ttm+tty", bins_toppt_unfolding, {nvecttm_tty, nvecttm_tty}, categories);
		unfolding2d_.AddMatrix("ttm+cts", ttmbins_large, {ctsbins_large, ctsbins_large, ctsbins_large, ctsbins_large}, categories);
		vector<double> nvecttm_cts = getnvec(unfolding2d_.GetNBins("ttm+cts"));
		unfolding2d_.AddMatrix("thadpt+ttm+cts", bins_toppt_unfolding, {nvecttm_cts, nvecttm_cts}, categories);
		if(PDFUNCS_)
		{
			pdf1d_.AddArray("nnpdf31", &GenInfo::Num_NNPDFWeights, &GenInfo::NNPDFWeights);
			pdf1d_.AddArray("nnpdf31as", &GenInfo::Num_NNPDFasWeights, &GenInfo::NNPDFasWeights);
			pdf1d_.AddArray("scales", &GenInfo::Num_ScaleWeights, &GenInfo::ScaleWeights);
			for(const string& name : {"thadpt", "thady", "tleppt", "tlepy", "thardpt", "tsoftpt", "toppt", "topy", "topbarpt", "topbary", "ttm", "ttpt", "tty", "cts", "thadpt+ttm+tty", "thadpt+ttm+cts"})
			{
				int nbins = unfolding2d_.GetNBins(name);
				pdf1d_.AddHIST(TH1D((name+"_gen").c_str(), (name+"_gen").c_str(), nbins, 0, nbins));
				pdf1d_.AddHIST(TH1D((name+"_rec").c_str(), (name+"_rec").c_str(), nbins, 0, nbins));
			}
		}
	}
}

void Unfolding::FillTruth(const TTEvent& evt, double weight)
{
	if(GLAN->UserTime() == 2023)
	{
		unfolding2d_.FillTruth("thady+thadpt", abs(evt.THad()->Rapidity()), evt.THad()->Pt(), weight);
		unfolding2d_.FillTruth("ttm+tty", evt.TT()->M(), abs(evt.TT()->Rapidity()), weight);
		unfolding1d_.FillTruth("thadpt", evt.THad()->Pt(), weight);
		unfolding1d_.FillTruth("thady", abs(evt.THad()->Rapidity()), weight);
		unfolding1d_.FillTruth("tleppt", evt.TLep()->Pt(), weight);
		unfolding1d_.FillTruth("tlepy", abs(evt.TLep()->Rapidity()), weight);
		unfolding1d_.FillTruth("ttm", evt.TT()->M(), weight);
		unfolding1d_.FillTruth("ttpt", evt.TT()->Pt(), weight);
		unfolding1d_.FillTruth("tty", abs(evt.TT()->Rapidity()), weight);
	}
	else
	{
		unfolding2d_.FillTruth("thadpt", evt.THad()->Pt(), evt.THad()->Pt(), weight);
		unfolding2d_.FillTruth("tleppt", evt.THad()->Pt(), evt.TLep()->Pt(), weight);
		unfolding2d_.FillTruth("thardpt", evt.THad()->Pt(), evt.THard()->Pt(), weight);
		unfolding2d_.FillTruth("tsoftpt", evt.THad()->Pt(), evt.TSoft()->Pt(), weight);
		unfolding2d_.FillTruth("toppt", evt.THad()->Pt(), evt.T()->Pt(), weight);
		unfolding2d_.FillTruth("topy", evt.THad()->Pt(), abs(evt.T()->Rapidity()), weight);
		unfolding2d_.FillTruth("topbarpt", evt.THad()->Pt(), evt.TBar()->Pt(), weight);
		unfolding2d_.FillTruth("topbary", evt.THad()->Pt(), abs(evt.TBar()->Rapidity()), weight);
		unfolding2d_.FillTruth("thady", evt.THad()->Pt(), abs(evt.THad()->Rapidity()), weight);
		unfolding2d_.FillTruth("tlepy", evt.THad()->Pt(), abs(evt.TLep()->Rapidity()), weight);
		unfolding2d_.FillTruth("ttm", evt.THad()->Pt(), evt.TT()->M(), weight);
		unfolding2d_.FillTruth("ttpt", evt.THad()->Pt(), evt.TT()->Pt(), weight);
		unfolding2d_.FillTruth("tty", evt.THad()->Pt(), abs(evt.TT()->Rapidity()), weight);
		unfolding2d_.FillTruth("cts", evt.THad()->Pt(), evt.CosThetaCMS(), weight);
		unfolding2d_.FillTruth("thadpt+ttm+tty", evt.THad()->Pt(), GetBin2D("ttm+tty", evt.TT()->M(), abs(evt.TT()->Rapidity())), weight);
		unfolding2d_.FillTruth("thadpt+ttm+cts", evt.THad()->Pt(), GetBin2D("ttm+cts", evt.TT()->M(), evt.CosThetaCMS()), weight);
		if(PDFUNCS_)
		{
			double defaultpdfweight = AN->defaultpdfweight_;
			pdf1d_.Fill("thadpt_gen", GetBin2D("thadpt", evt.THad()->Pt(), evt.THad()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("tleppt_gen", GetBin2D("tleppt", evt.THad()->Pt(), evt.TLep()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("thardpt_gen", GetBin2D("thardpt", evt.THad()->Pt(), evt.THard()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("tsoftpt_gen", GetBin2D("tsoftpt", evt.THad()->Pt(), evt.TSoft()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("toppt_gen", GetBin2D("toppt", evt.THad()->Pt(), evt.T()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("topy_gen", GetBin2D("topy", evt.THad()->Pt(), abs(evt.T()->Rapidity())), weight, defaultpdfweight);
			pdf1d_.Fill("topbarpt_gen", GetBin2D("topbarpt", evt.THad()->Pt(), evt.TBar()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("topbary_gen", GetBin2D("topbary", evt.THad()->Pt(), abs(evt.TBar()->Rapidity())), weight, defaultpdfweight);
			pdf1d_.Fill("thady_gen", GetBin2D("thady", evt.THad()->Pt(), abs(evt.THad()->Rapidity())), weight, defaultpdfweight);
			pdf1d_.Fill("tlepy_gen", GetBin2D("tlepy", evt.THad()->Pt(), abs(evt.TLep()->Rapidity())), weight, defaultpdfweight);
			pdf1d_.Fill("ttm_gen", GetBin2D("ttm", evt.THad()->Pt(), evt.TT()->M()), weight, defaultpdfweight);
			pdf1d_.Fill("ttpt_gen", GetBin2D("ttpt", evt.THad()->Pt(), evt.TT()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("tty_gen", GetBin2D("tty", evt.THad()->Pt(), abs(evt.TT()->Rapidity())), weight, defaultpdfweight);
			pdf1d_.Fill("cts_gen", GetBin2D("cts", evt.THad()->Pt(), evt.CosThetaCMS()), weight, defaultpdfweight);
			pdf1d_.Fill("thadpt+ttm+tty_gen", GetBin2D("thadpt+ttm+tty", evt.THad()->Pt(), GetBin2D("ttm+tty", evt.TT()->M(), abs(evt.TT()->Rapidity()))), weight, defaultpdfweight);
			pdf1d_.Fill("thadpt+ttm+cts_gen", GetBin2D("thadpt+ttm+cts", evt.THad()->Pt(), GetBin2D("ttm+cts", evt.TT()->M(), evt.CosThetaCMS())), weight, defaultpdfweight);
		}
	}
}

void Unfolding::FillTruthReco(const TTEvent& truth_evt, const TTEvent& reco_evt, const string& category, double weight)
{
	if(GLAN->UserTime() == 2023)
	{
		unfolding2d_.FillTruthReco("thady+thadpt_"+category, abs(truth_evt.THad()->Rapidity()), truth_evt.THad()->Pt(), abs(reco_evt.THad()->Rapidity()), reco_evt.THad()->Pt(), weight);
		unfolding2d_.FillTruthReco("ttm+tty_"+category, truth_evt.TT()->M(), abs(truth_evt.TT()->Rapidity()), reco_evt.TT()->M(), abs(reco_evt.TT()->Rapidity()), weight);
		unfolding1d_.FillTruthReco("thadpt_"+category, truth_evt.THad()->Pt(), reco_evt.THad()->Pt(), weight);
		unfolding1d_.FillTruthReco("thady_"+category, abs(truth_evt.THad()->Rapidity()), abs(reco_evt.THad()->Rapidity()), weight);
		unfolding1d_.FillTruthReco("tleppt_"+category, truth_evt.TLep()->Pt(), reco_evt.TLep()->Pt(), weight);
		unfolding1d_.FillTruthReco("tlepy_"+category, abs(truth_evt.TLep()->Rapidity()), abs(reco_evt.TLep()->Rapidity()), weight);
		unfolding1d_.FillTruthReco("ttm_"+category, truth_evt.TT()->M(), reco_evt.TT()->M(), weight);
		unfolding1d_.FillTruthReco("ttpt_"+category, truth_evt.TT()->Pt(), reco_evt.TT()->Pt(), weight);
		unfolding1d_.FillTruthReco("tty_"+category, abs(truth_evt.TT()->Rapidity()), abs(reco_evt.TT()->Rapidity()), weight);

	}
	else
	{
		unfolding2d_.FillTruthReco("thadpt_"+category, truth_evt.THad()->Pt(), truth_evt.THad()->Pt(), reco_evt.THad()->Pt(), reco_evt.THad()->Pt(), weight);
		unfolding2d_.FillTruthReco("tleppt_"+category, truth_evt.THad()->Pt(), truth_evt.TLep()->Pt(), reco_evt.THad()->Pt(), reco_evt.TLep()->Pt(), weight);
		unfolding2d_.FillTruthReco("thardpt_"+category, truth_evt.THad()->Pt(), truth_evt.THard()->Pt(), reco_evt.THad()->Pt(), reco_evt.THard()->Pt(), weight);
		unfolding2d_.FillTruthReco("tsoftpt_"+category, truth_evt.THad()->Pt(), truth_evt.TSoft()->Pt(), reco_evt.THad()->Pt(), reco_evt.TSoft()->Pt(), weight);
		unfolding2d_.FillTruthReco("toppt_"+category, truth_evt.THad()->Pt(), truth_evt.T()->Pt(), reco_evt.THad()->Pt(), reco_evt.T()->Pt(), weight);
		unfolding2d_.FillTruthReco("topy_"+category, truth_evt.THad()->Pt(), abs(truth_evt.T()->Rapidity()), reco_evt.THad()->Pt(), abs(reco_evt.T()->Rapidity()), weight);
		unfolding2d_.FillTruthReco("topbarpt_"+category, truth_evt.THad()->Pt(), truth_evt.TBar()->Pt(), reco_evt.THad()->Pt(), reco_evt.TBar()->Pt(), weight);
		unfolding2d_.FillTruthReco("topbary_"+category, truth_evt.THad()->Pt(), abs(truth_evt.TBar()->Rapidity()), reco_evt.THad()->Pt(), abs(reco_evt.TBar()->Rapidity()), weight);
		unfolding2d_.FillTruthReco("thady_"+category, truth_evt.THad()->Pt(), abs(truth_evt.THad()->Rapidity()), reco_evt.THad()->Pt(), abs(reco_evt.THad()->Rapidity()), weight);
		unfolding2d_.FillTruthReco("tlepy_"+category, truth_evt.THad()->Pt(), abs(truth_evt.TLep()->Rapidity()), reco_evt.THad()->Pt(), abs(reco_evt.TLep()->Rapidity()), weight);
		unfolding2d_.FillTruthReco("ttm_"+category, truth_evt.THad()->Pt(), truth_evt.TT()->M(), reco_evt.THad()->Pt(), reco_evt.TT()->M(), weight);
		unfolding2d_.FillTruthReco("ttpt_"+category, truth_evt.THad()->Pt(), truth_evt.TT()->Pt(), reco_evt.THad()->Pt(), reco_evt.TT()->Pt(), weight);
		unfolding2d_.FillTruthReco("tty_"+category, truth_evt.THad()->Pt(), abs(truth_evt.TT()->Rapidity()), reco_evt.THad()->Pt(), abs(reco_evt.TT()->Rapidity()), weight);
		unfolding2d_.FillTruthReco("cts_"+category, truth_evt.THad()->Pt(), truth_evt.CosThetaCMS(), reco_evt.THad()->Pt(), reco_evt.CosThetaCMS(), weight);
		unfolding2d_.FillTruthReco("thadpt+ttm+tty_"+category, truth_evt.THad()->Pt(), GetBin2D("ttm+tty", truth_evt.TT()->M(), abs(truth_evt.TT()->Rapidity())), reco_evt.THad()->Pt(), GetBin2D("ttm+tty", reco_evt.TT()->M(), abs(reco_evt.TT()->Rapidity())), weight);
		unfolding2d_.FillTruthReco("thadpt+ttm+cts_"+category, truth_evt.THad()->Pt(), GetBin2D("ttm+cts", truth_evt.TT()->M(), truth_evt.CosThetaCMS()), reco_evt.THad()->Pt(), GetBin2D("ttm+cts", reco_evt.TT()->M(), reco_evt.CosThetaCMS()), weight);
		if(PDFUNCS_)
		{
			double defaultpdfweight = AN->defaultpdfweight_;
			pdf1d_.Fill("thadpt_rec", GetBin2D("thadpt", truth_evt.THad()->Pt(), truth_evt.THad()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("tleppt_rec", GetBin2D("tleppt", truth_evt.THad()->Pt(), truth_evt.TLep()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("thardpt_rec", GetBin2D("thardpt", truth_evt.THad()->Pt(), truth_evt.THard()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("tsoftpt_rec", GetBin2D("tsoftpt", truth_evt.THad()->Pt(), truth_evt.TSoft()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("toppt_rec", GetBin2D("toppt", truth_evt.THad()->Pt(), truth_evt.T()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("topy_rec", GetBin2D("topy", truth_evt.THad()->Pt(), abs(truth_evt.T()->Rapidity())), weight, defaultpdfweight);
			pdf1d_.Fill("topbarpt_rec", GetBin2D("topbarpt", truth_evt.THad()->Pt(), truth_evt.TBar()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("topbary_rec", GetBin2D("topbary", truth_evt.THad()->Pt(), abs(truth_evt.TBar()->Rapidity())), weight, defaultpdfweight);
			pdf1d_.Fill("thady_rec", GetBin2D("thady", truth_evt.THad()->Pt(), abs(truth_evt.THad()->Rapidity())), weight, defaultpdfweight);
			pdf1d_.Fill("tlepy_rec", GetBin2D("tlepy", truth_evt.THad()->Pt(), abs(truth_evt.TLep()->Rapidity())), weight, defaultpdfweight);
			pdf1d_.Fill("ttm_rec", GetBin2D("ttm", truth_evt.THad()->Pt(), truth_evt.TT()->M()), weight, defaultpdfweight);
			pdf1d_.Fill("ttpt_rec", GetBin2D("ttpt", truth_evt.THad()->Pt(), truth_evt.TT()->Pt()), weight, defaultpdfweight);
			pdf1d_.Fill("tty_rec", GetBin2D("tty", truth_evt.THad()->Pt(), abs(truth_evt.TT()->Rapidity())), weight, defaultpdfweight);
			pdf1d_.Fill("cts_rec", GetBin2D("cts", truth_evt.THad()->Pt(), truth_evt.CosThetaCMS()), weight, defaultpdfweight);
			pdf1d_.Fill("thadpt+ttm+tty_rec", GetBin2D("thadpt+ttm+tty", truth_evt.THad()->Pt(), GetBin2D("ttm+tty", truth_evt.TT()->M(), abs(truth_evt.TT()->Rapidity()))), weight, defaultpdfweight);
			pdf1d_.Fill("thadpt+ttm+cts_rec", GetBin2D("thadpt+ttm+cts", truth_evt.THad()->Pt(), GetBin2D("ttm+cts", truth_evt.TT()->M(), truth_evt.CosThetaCMS())), weight, defaultpdfweight);
		}
	}
}

void Unfolding::FillReco(const TTEvent& evt, const string& category, double weight)
{
	if(GLAN->UserTime() == 2023)
	{
		unfolding2d_.FillReco("thady+thadpt_"+category, abs(evt.THad()->Rapidity()), evt.THad()->Pt(), weight);
		unfolding2d_.FillReco("ttm+tty_"+category, evt.TT()->M(), abs(evt.TT()->Rapidity()), weight);
		unfolding1d_.FillReco("thadpt_"+category, evt.THad()->Pt(), weight);
		unfolding1d_.FillReco("thady_"+category, abs(evt.THad()->Rapidity()), weight);
		unfolding1d_.FillReco("tleppt_"+category, evt.TLep()->Pt(), weight);
		unfolding1d_.FillReco("tlepy_"+category, abs(evt.TLep()->Rapidity()), weight);
		unfolding1d_.FillReco("ttm_"+category, evt.TT()->M(), weight);
		unfolding1d_.FillReco("ttpt_"+category, evt.TT()->Pt(), weight);
		unfolding1d_.FillReco("tty_"+category, abs(evt.TT()->Rapidity()), weight);

	}
	else
	{
		unfolding2d_.FillReco("thadpt_"+category, evt.THad()->Pt(), evt.THad()->Pt(), weight);
		unfolding2d_.FillReco("tleppt_"+category, evt.THad()->Pt(), evt.TLep()->Pt(), weight);
		unfolding2d_.FillReco("thardpt_"+category, evt.THad()->Pt(), evt.THard()->Pt(), weight);
		unfolding2d_.FillReco("tsoftpt_"+category, evt.THad()->Pt(), evt.TSoft()->Pt(), weight);
		unfolding2d_.FillReco("toppt_"+category, evt.THad()->Pt(), evt.T()->Pt(), weight);
		unfolding2d_.FillReco("topy_"+category, evt.THad()->Pt(), abs(evt.T()->Rapidity()), weight);
		unfolding2d_.FillReco("topbarpt_"+category, evt.THad()->Pt(), evt.TBar()->Pt(), weight);
		unfolding2d_.FillReco("topbary_"+category, evt.THad()->Pt(), abs(evt.TBar()->Rapidity()), weight);
		unfolding2d_.FillReco("thady_"+category, evt.THad()->Pt(), abs(evt.THad()->Rapidity()), weight);
		unfolding2d_.FillReco("tlepy_"+category, evt.THad()->Pt(), abs(evt.TLep()->Rapidity()), weight);
		unfolding2d_.FillReco("ttm_"+category, evt.THad()->Pt(), evt.TT()->M(), weight);
		unfolding2d_.FillReco("ttpt_"+category, evt.THad()->Pt(), evt.TT()->Pt(), weight);
		unfolding2d_.FillReco("tty_"+category, evt.THad()->Pt(), abs(evt.TT()->Rapidity()), weight);
		unfolding2d_.FillReco("cts_"+category, evt.THad()->Pt(), evt.CosThetaCMS(), weight);
		unfolding2d_.FillReco("thadpt+ttm+tty_"+category, evt.THad()->Pt(), GetBin2D("ttm+tty", evt.TT()->M(), abs(evt.TT()->Rapidity())), weight);
		unfolding2d_.FillReco("thadpt+ttm+cts_"+category, evt.THad()->Pt(), GetBin2D("ttm+cts", evt.TT()->M(), evt.CosThetaCMS()), weight);
	}
}

void Unfolding::FillRes(const TTEvent& truth_evt, const TTEvent& reco_evt, const string& category, double weight)
{
	if(GLAN->UserTime() == 2023)
	{
		unfolding2d_.FillRes("thady+thadpt_"+category, abs(truth_evt.THad()->Rapidity()), truth_evt.THad()->Pt(), abs(reco_evt.THad()->Rapidity()), reco_evt.THad()->Pt(), weight);
		unfolding2d_.FillRes("ttm+tty_"+category, truth_evt.TT()->M(), abs(truth_evt.TT()->Rapidity()), reco_evt.TT()->M(), abs(reco_evt.TT()->Rapidity()), weight);
		unfolding1d_.FillRes("thadpt_"+category, truth_evt.THad()->Pt(), reco_evt.THad()->Pt(), weight);
		unfolding1d_.FillRes("thady_"+category, abs(truth_evt.THad()->Rapidity()), abs(reco_evt.THad()->Rapidity()), weight);
		unfolding1d_.FillRes("tleppt_"+category, truth_evt.TLep()->Pt(), reco_evt.TLep()->Pt(), weight);
		unfolding1d_.FillRes("tlepy_"+category, abs(truth_evt.TLep()->Rapidity()), abs(reco_evt.TLep()->Rapidity()), weight);
		unfolding1d_.FillRes("ttm_"+category, truth_evt.TT()->M(), reco_evt.TT()->M(), weight);
		unfolding1d_.FillRes("ttpt_"+category, truth_evt.TT()->Pt(), reco_evt.TT()->Pt(), weight);
		unfolding1d_.FillRes("tty_"+category, abs(truth_evt.TT()->Rapidity()), abs(reco_evt.TT()->Rapidity()), weight);

	}
	else
	{
		unfolding2d_.FillRes("thadpt_"+category, truth_evt.THad()->Pt(), truth_evt.THad()->Pt(), reco_evt.THad()->Pt(), reco_evt.THad()->Pt(), weight);
		unfolding2d_.FillRes("tleppt_"+category, truth_evt.THad()->Pt(), truth_evt.TLep()->Pt(), reco_evt.THad()->Pt(), reco_evt.TLep()->Pt(), weight);
		unfolding2d_.FillRes("thardpt_"+category, truth_evt.THad()->Pt(), truth_evt.THard()->Pt(), reco_evt.THad()->Pt(), reco_evt.THard()->Pt(), weight);
		unfolding2d_.FillRes("tsoftpt_"+category, truth_evt.THad()->Pt(), truth_evt.TSoft()->Pt(), reco_evt.THad()->Pt(), reco_evt.TSoft()->Pt(), weight);
		unfolding2d_.FillRes("toppt_"+category, truth_evt.THad()->Pt(), truth_evt.T()->Pt(), reco_evt.THad()->Pt(), reco_evt.T()->Pt(), weight);
		unfolding2d_.FillRes("topy_"+category, truth_evt.THad()->Pt(), abs(truth_evt.T()->Rapidity()), reco_evt.THad()->Pt(), abs(reco_evt.T()->Rapidity()), weight);
		unfolding2d_.FillRes("topbarpt_"+category, truth_evt.THad()->Pt(), truth_evt.TBar()->Pt(), reco_evt.THad()->Pt(), reco_evt.TBar()->Pt(), weight);
		unfolding2d_.FillRes("topbary_"+category, truth_evt.THad()->Pt(), abs(truth_evt.TBar()->Rapidity()), reco_evt.THad()->Pt(), abs(reco_evt.TBar()->Rapidity()), weight);
		unfolding2d_.FillRes("thady_"+category, truth_evt.THad()->Pt(), abs(truth_evt.THad()->Rapidity()), reco_evt.THad()->Pt(), abs(reco_evt.THad()->Rapidity()), weight);
		unfolding2d_.FillRes("tlepy_"+category, truth_evt.THad()->Pt(), abs(truth_evt.TLep()->Rapidity()), reco_evt.THad()->Pt(), abs(reco_evt.TLep()->Rapidity()), weight);
		unfolding2d_.FillRes("ttm_"+category, truth_evt.THad()->Pt(), truth_evt.TT()->M(), reco_evt.THad()->Pt(), reco_evt.TT()->M(), weight);
		unfolding2d_.FillRes("ttpt_"+category, truth_evt.THad()->Pt(), truth_evt.TT()->Pt(), reco_evt.THad()->Pt(), reco_evt.TT()->Pt(), weight);
		unfolding2d_.FillRes("tty_"+category, truth_evt.THad()->Pt(), abs(truth_evt.TT()->Rapidity()), reco_evt.THad()->Pt(), abs(reco_evt.TT()->Rapidity()), weight);
		unfolding2d_.FillRes("cts_"+category, truth_evt.THad()->Pt(), truth_evt.CosThetaCMS(), reco_evt.THad()->Pt(), reco_evt.CosThetaCMS(), weight);
		unfolding2d_.FillRes("ttm+tty_"+category, truth_evt.TT()->M(), abs(truth_evt.TT()->Rapidity()), reco_evt.TT()->M(), abs(reco_evt.TT()->Rapidity()), weight);
		unfolding2d_.FillRes("ttm+cts_"+category, truth_evt.TT()->M(), truth_evt.CosThetaCMS(), reco_evt.TT()->M(), reco_evt.CosThetaCMS(), weight);
	}
}
