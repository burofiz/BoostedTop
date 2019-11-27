#include <TLorentzVector.h>
#include <OJet.h>

#include "TTPlots.h"
#include "TTEvent.h"
#include "TopJets.h"
#include "BoostedTop.h"

Plots::Plots(const string& prefix) : _prefix(prefix), _hists1d(""), _hists2d("")
{

}

TTPlots::TTPlots(const string& prefix) : Plots(prefix)
{

}

void TTPlots::Init()
{
	TDirectory* olddir = gDirectory;
	TDirectory* newdir = gDirectory->mkdir(_prefix.c_str());
	newdir->cd();
	const vector<double>& bins_toppt_large = dynamic_cast<BoostedTop*>(GLAN)->bins_toppt_large;
	const vector<double>& bins_toppt = dynamic_cast<BoostedTop*>(GLAN)->bins_toppt;
	for(const string& leptype : leptypes_)
	{
		_hists1d.AddHist(leptype+"_nvtx_unweighted", 150, 0, 300, "number of vertices (unweighted)", "");
		_hists1d.AddHist(leptype+"_nvtx", 150, 0, 300, "number of vertices (weighted)", "");
		_hists1d.AddHist(leptype+"_njets", 10, 0, 10, "additional jets", "");
		_hists1d.AddHist(leptype+"_nhtjets", 10, 0, 10, "additional boosted t_{h}", "");
		_hists1d.AddHist(leptype+"_nhtjets_nnres", 100, 0, 1, "additional t_{h} N_{h}", "");
		_hists2d.AddHist(leptype+"_thadpt_njets", bins_toppt_large, 10, 0, 10, "p_{T}(t_{h}) [GeV]", "additional jets");
		_hists1d.AddHist(leptype+"_wm", 100, 0, 200, "M_{W} [GeV]", "");
		_hists1d.AddHist(leptype+"_tm", 200, 0, 400, "M_{t} [GeV]", "");
		_hists2d.AddHist(leptype+"_tm_thadpt", 150, 0, 300, 150, 0, 1500, "M_{t} [GeV]", "p_{T}(t_{h}) [GeV]");
		_hists2d.AddHist(leptype+"_tm_thadpt_eta0", 150, 0, 300, bins_toppt_large, "M_{t} [GeV]", "p_{T}(t_{h}) [GeV]");
		_hists2d.AddHist(leptype+"_tm_thadpt_eta1", 150, 0, 300, bins_toppt_large, "M_{t} [GeV]", "p_{T}(t_{h}) [GeV]");
		_hists2d.AddHist(leptype+"_tm_thadpt_eta2", 150, 0, 300, bins_toppt_large, "M_{t} [GeV]", "p_{T}(t_{h}) [GeV]");
		_hists2d.AddHist(leptype+"_tm_thadpt_eta3", 150, 0, 300, bins_toppt_large, "M_{t} [GeV]", "p_{T}(t_{h}) [GeV]");
		_hists1d.AddHist(leptype+"_thadpt", 400, 0, 2000, "p_{T}(t_{h}) [GeV]", "");
		_hists1d.AddHist(leptype+"_thardpt", 400, 0, 2000, "p_{T}(t_{high}) [GeV]", "");
		_hists1d.AddHist(leptype+"_tsoftpt", 400, 0, 2000, "p_{T}(t_{low}) [GeV]", "");
		_hists1d.AddHist(leptype+"_toppt", 400, 0, 2000, "p_{T}(t) [GeV]", "");
		_hists1d.AddHist(leptype+"_topy", 100, 0, 5, "|y(t)|", "");
		_hists1d.AddHist(leptype+"_topeta", 100, 0, 5, "|#eta(t)|", "");
		_hists1d.AddHist(leptype+"_topbarpt", 400, 0, 2000, "p_{T}(#bar{t}) [GeV]", "");
		_hists1d.AddHist(leptype+"_topbary", 50, 0, 5, "|y(#bar{t})|", "");
		_hists1d.AddHist(leptype+"_topbareta", 50, 0, 5, "|#eta(#bar{t})|", "");
		_hists2d.AddHist(leptype+"_thadpt_thadeta", bins_toppt_large, 50, 0, 5., "p_{T}(t_{h}) [GeV]", "#eta(t_{h})");
		_hists1d.AddHist(leptype+"_thady", 100, 0, 5, "|y(t_{h})|", "");
		_hists1d.AddHist(leptype+"_tleppt", 400, 0, 2000, "p_{T}(t_{l}) [GeV]", "");
		_hists2d.AddHist(leptype+"_tleppt_thadpt", 75, 0, 1500, 75, 0, 1500,  "p_{T}(t_{l}) [GeV]", "p_{T}(t_{h}) [GeV]");
		_hists1d.AddHist(leptype+"_tlepy", 100, 0, 5, "|y(t_{l})|", "");
		_hists1d.AddHist(leptype+"_ttm", 500, 0, 5000, "M(t#bar{t}) [GeV]", "");
		_hists1d.AddHist(leptype+"_ttpt", 200, 0, 1000, "p_{T}(t#bar{t}) [GeV]", "");
		_hists1d.AddHist(leptype+"_tty", 250, 0, 5, "|y(t#bar{t})|", "");
		_hists2d.AddHist(leptype+"_ttdy_FBthadpt", 100, -5, 5, bins_toppt_large, "#Deltay(t#bar{t})", "p_{T}(t_{h}) [GeV]");
		_hists2d.AddHist(leptype+"_ttdphi_FBthadpt", 100, 0, Pi(), bins_toppt_large, "#Delta #phi(t#bar{t})", "p_{T}(t_{h}) [GeV]");
		_hists2d.AddHist(leptype+"_ht_FBthadpt", 100, 0, 1000, bins_toppt_large, "H_{T}", "p_{T}(t_{h}) [GeV]");
		_hists2d.AddHist(leptype+"_cts_FBthadpt", 100, -1, 1, bins_toppt_large, "cos(#theta)", "p_{T}(t_{h}) [GeV]");
		vector<double> ttmbins_large = {300.0, 500.0, 700.0, 900.0, 13000.0};
		_hists2d.AddHist(leptype+"_bstar_FBttm", 100, -1, 1, ttmbins_large, "cos(#theta_{b})", "M(t#bar{t}) [GeV]");
		_hists2d.AddHist(leptype+"_bbarstar_FBttm", 100, -1, 1, ttmbins_large, "cos(#theta_{#bar{b}})", "M(t#bar{t}) [GeV]");
		_hists2d.AddHist(leptype+"_deltabbstar_FBttm", 100, -1, 1, ttmbins_large, "cos(#theta_{b#bar{b}})", "M(t#bar{t}) [GeV]");
		_hists1d.AddHist(leptype+"_nupt", 250, 0, 1000, "p_{T}(#nu) [GeV]", "");
		_hists1d.AddHist(leptype+"_nueta", 100, -2.5, 2.5, "#eta(#nu)", "");
		_hists1d.AddHist(leptype+"_blpt", 200, 0, 1000, "p_{T}(bl) [GeV]", "");
		_hists1d.AddHist(leptype+"_bhpt", 200, 0, 1000, "p_{T}(bh) [GeV]", "");
		_hists1d.AddHist(leptype+"_jet1pt", 200, 0, 1000, "p_{T}(j_{1}) [GeV]", "");
		_hists1d.AddHist(leptype+"_jet2pt", 200, 0, 800, "p_{T}(j_{2}) [GeV]", "");
		_hists1d.AddHist(leptype+"_jet3pt", 200, 0, 400, "p_{T}(j_{3}) [GeV]", "");
		_hists1d.AddHist(leptype+"_wlowpt", 200, 0, 1000, "p_{T}(Wlow) [GeV]", "");
		_hists1d.AddHist(leptype+"_whighpt", 200, 0, 1000, "p_{T}(Whigh) [GeV]", "");
		_hists1d.AddHist(leptype+"_bleta", 200, -5, 5, "#eta(bl)", "");
		_hists1d.AddHist(leptype+"_bheta", 200, -5, 5, "#eta(bh)", "");
		_hists2d.AddHist(leptype+"_bhetaphi", 50, -2.5, 2.5, 90, -Pi(), Pi(), "#eta(bh)", "#phi(bh)");
		_hists1d.AddHist(leptype+"_wloweta", 200, -5, 5, "#eta(Wlow)", "");
		_hists1d.AddHist(leptype+"_whigheta", 200, -5, 5, "#eta(Whigh)", "");
		_hists1d.AddHist(leptype+"_leppt", 200, 0, 400, "p_{T}(l) [GeV]", "");
		_hists1d.AddHist(leptype+"_lepeta", 200, -5.0, 5.0, "#eta(l)", "");
		_hists1d.AddHist(leptype+"_lepptFW", 200, 0, 400, "p_{T}(l) [GeV] for |y(t#bar{t})| > 2.4", "");
		_hists1d.AddHist(leptype+"_lepetaFW", 200, -5.0, 5.0, "#eta(l) for |y(t#bar{t})| > 2.4", "");
		_hists2d.AddHist(leptype+"_lepetapt", 48, -2.4, 2.4, 15, 30, 180, "#eta(l)", "p_{T}(l) [GeV]");
		_hists2d.AddHist(leptype+"_lepetaphi", {-2.4, -1.5, 0., 1.5, 2.4}, 180, -Pi(), Pi(), "#eta(l)", "#phi(l)");
		_hists2d.AddHist(leptype+"_thadpt_WMT", bins_toppt, 100, 0, 200, "p_{T}(t_{h}) [GeV]", "M_{T}(W) [GeV]");
		_hists1d.AddHist(leptype+"_met", 250, 0, 1000, "p_{T}^{miss} [GeV]", "");
		_hists1d.AddHist(leptype+"_metphi", 180, -Pi(), Pi(), "#phi^{miss}", "");
		_hists1d.AddHist(leptype+"_thadpt_FB", bins_toppt, "p_{T}(t_{h}) [GeV]", "");
		_hists1d.AddHist(leptype+"_tleppt_FB", bins_toppt, "p_{T}(t_{l}) [GeV]", "");
		_hists1d.AddHist(leptype+"_ttm_FB", dynamic_cast<BoostedTop*>(GLAN)->bins_ttm, "M(t#bar{t}) [GeV]", "");
		_hists1d.AddHist(leptype+"_ttpt_FB", dynamic_cast<BoostedTop*>(GLAN)->bins_ttpt, "p_{T}(t#bar{t}) [GeV]", "");
		_hists2d.AddHist(leptype+"_thadpt_prob", bins_toppt, 60, 10, 25, "p_{T}(t_{h}) [GeV]", "-log(#lambda)");
		_hists2d.AddHist(leptype+"_thadpt_probhad", bins_toppt, 40, 5, 15, "p_{T}(t_{h}) [GeV]", "-log(#lambda_{had})");
		_hists2d.AddHist(leptype+"_tleppt_problep", bins_toppt, 40, 0, 10, "p_{T}(t_{l}) [GeV]", "-log(#lambda_{lep})");
		_hists1d.AddHist(leptype+"_dnu", 100, 0, 200, "D_{#nu}", "");
		_hists1d.AddHist(leptype+"_isonear", 100, 0, 200, "I_{near}", "");
	}
	olddir->cd();
}

void TTPlots::Fill(const string& leptype, TTEvent& ev, double weight)
{
	BoostedTop* BT = dynamic_cast<BoostedTop*>(GLAN);
	
	_hists1d[leptype+"_nvtx"]->Fill(BT->NumGoodVertices(), weight);
	_hists1d[leptype+"_nvtx_unweighted"]->Fill(BT->NumGoodVertices(), BT->defaultpdfweight_*BT->lepweight_*BT->btagweight_*BT->prefireweight_);

	_hists1d[leptype+"_njets"]->Fill(ev.NAddJets(), weight);
	_hists2d[leptype+"_thadpt_njets"]->Fill(ev.THad()->Pt(), ev.NAddJets(), weight);

	const vector<const HadTopJet* > addhtjets = ev.AddHTJets();
	_hists1d[leptype+"_nhtjets"]->Fill(addhtjets.size(), weight);
	if(addhtjets.size() >= 2) {_hists1d[leptype+"_nhtjets_nnres"]->Fill(addhtjets[0]->NNRes()*addhtjets[1]->NNRes(), weight);}
	if(ev.IsCompleteRH())
	{
		_hists1d[leptype+"_wm"]->Fill((*ev.WJPtmax() + *ev.WJPtmin()).M(), weight);
		_hists1d[leptype+"_bhpt"]->Fill(ev.BHad()->Pt(), weight);
		_hists1d[leptype+"_wlowpt"]->Fill(ev.WJPtmin()->Pt(), weight);
		_hists1d[leptype+"_whighpt"]->Fill(ev.WJPtmax()->Pt(), weight);
		_hists1d[leptype+"_bheta"]->Fill(ev.BHad()->Eta(), weight);
		_hists2d[leptype+"_bhetaphi"]->Fill(ev.BHad()->Eta(), ev.BHad()->Phi(), weight);
		_hists1d[leptype+"_wloweta"]->Fill(ev.WJPtmin()->Eta(), weight);
		_hists1d[leptype+"_whigheta"]->Fill(ev.WJPtmax()->Eta(), weight);
	}
	if(ev.IsCompleteRL())
	{
		_hists1d[leptype+"_blpt"]->Fill(ev.BLep()->Pt(), weight);
		_hists1d[leptype+"_bleta"]->Fill(ev.BLep()->Eta(), weight);
	}
	if(ev.Prob() > 0)
	{
		_hists2d[leptype+"_thadpt_prob"]->Fill(ev.THad()->Pt(),ev.Prob(), weight);
		_hists2d[leptype+"_thadpt_probhad"]->Fill(ev.THad()->Pt(), ev.ProbHad(), weight);
		_hists2d[leptype+"_tleppt_problep"]->Fill(ev.TLep()->Pt(), ev.ProbLep(), weight);
		_hists1d[leptype+"_dnu"]->Fill(ev.Dnu(), weight);
	}
	_hists1d[leptype+"_tm"]->Fill(ev.THad()->M(), weight);
	_hists2d[leptype+"_tm_thadpt"]->Fill(ev.THad()->M(), ev.THad()->Pt(), weight);
	if(ev.THad()->Eta() < -1.) {_hists2d[leptype+"_tm_thadpt_eta0"]->Fill(ev.THad()->M(), ev.THad()->Pt(), weight);}
	else if(ev.THad()->Eta() < 0.) {_hists2d[leptype+"_tm_thadpt_eta1"]->Fill(ev.THad()->M(), ev.THad()->Pt(), weight);}
	else if(ev.THad()->Eta() < 1.) {_hists2d[leptype+"_tm_thadpt_eta2"]->Fill(ev.THad()->M(), ev.THad()->Pt(), weight);}
	else {_hists2d[leptype+"_tm_thadpt_eta3"]->Fill(ev.THad()->M(), ev.THad()->Pt(), weight);}
	_hists1d[leptype+"_thadpt"]->Fill(ev.THad()->Pt(), weight);
	_hists1d[leptype+"_thardpt"]->Fill(ev.THard()->Pt(), weight);
	_hists1d[leptype+"_tsoftpt"]->Fill(ev.TSoft()->Pt(), weight);
	_hists1d[leptype+"_toppt"]->Fill(ev.T()->Pt(), weight);
	_hists1d[leptype+"_topy"]->Fill(abs(ev.T()->Rapidity()), weight);
	_hists1d[leptype+"_topeta"]->Fill(abs(ev.T()->Eta()), weight);
	_hists1d[leptype+"_topbarpt"]->Fill(ev.TBar()->Pt(), weight);
	_hists1d[leptype+"_topbary"]->Fill(abs(ev.TBar()->Rapidity()), weight);
	_hists1d[leptype+"_topbareta"]->Fill(abs(ev.TBar()->Eta()), weight);
	_hists2d[leptype+"_thadpt_thadeta"]->Fill(ev.THad()->Pt(), abs(ev.THad()->Eta()), weight);
	_hists1d[leptype+"_thadpt_FB"]->Fill(ev.THad()->Pt(), weight);
	_hists1d[leptype+"_thady"]->Fill(Abs(ev.THad()->Rapidity()), weight);
	_hists1d[leptype+"_tleppt"]->Fill(ev.TLep()->Pt(), weight);
	_hists1d[leptype+"_tleppt_FB"]->Fill(ev.TLep()->Pt(), weight);
	_hists2d[leptype+"_tleppt_thadpt"]->Fill(ev.TLep()->Pt(), ev.THad()->Pt(), weight);
	_hists1d[leptype+"_tlepy"]->Fill(Abs(ev.TLep()->Rapidity()), weight);
	_hists1d[leptype+"_ttm"]->Fill(ev.TT()->M(), weight);
	_hists1d[leptype+"_ttm_FB"]->Fill(ev.TT()->M(), weight);
	_hists1d[leptype+"_ttpt"]->Fill(ev.TT()->Pt(), weight);
	_hists1d[leptype+"_tty"]->Fill(Abs(ev.TT()->Rapidity()), weight);
	_hists2d[leptype+"_ttdy_FBthadpt"]->Fill(abs(ev.T()->Rapidity()) - abs(ev.TBar()->Rapidity()), ev.THad()->Pt(), weight);
	_hists2d[leptype+"_ttdphi_FBthadpt"]->Fill(abs(ev.T()->DeltaPhi(*ev.TBar())), ev.THad()->Pt(), weight);
	_hists2d[leptype+"_cts_FBthadpt"]->Fill(ev.CosThetaCMS(), ev.THad()->Pt(), weight);
	_hists2d[leptype+"_bstar_FBttm"]->Fill(ev.BPhiCMS(), ev.TT()->M(), weight);
	_hists2d[leptype+"_bbarstar_FBttm"]->Fill(ev.BBarPhiCMS(), ev.TT()->M(), weight);
	_hists2d[leptype+"_deltabbstar_FBttm"]->Fill(ev.DeltaBBPhiCMS(), ev.TT()->M(), weight);
	_hists2d[leptype+"_ht_FBthadpt"]->Fill(ev.Ht(), ev.THad()->Pt(), weight);
	if(ev.MET() != nullptr)
	{
		_hists1d[leptype+"_met"]->Fill(ev.MET()->Pt(), weight);
		_hists1d[leptype+"_metphi"]->Fill(ev.MET()->Phi(), ev.MET()->Pt()*weight);
	}	
	if(ev.NLep() != nullptr)
	{
		_hists1d[leptype+"_nupt"]->Fill(ev.NLep()->Pt(), weight);
		_hists1d[leptype+"_nueta"]->Fill(ev.NLep()->Eta(), weight);
	}
	if(ev.LLep() != nullptr)
	{
		_hists1d[leptype+"_leppt"]->Fill(ev.LLep()->Pt(), weight);
		_hists1d[leptype+"_lepeta"]->Fill(ev.LLep()->Eta(), weight);
		_hists2d[leptype+"_lepetapt"]->Fill(ev.LLep()->Eta(), ev.LLep()->Pt(),  weight);
		_hists2d[leptype+"_lepetaphi"]->Fill(ev.LLep()->Eta(), ev.LLep()->Phi(),  weight);

		if(abs(ev.TT()->Rapidity()) > 2.4)
		{
			_hists1d[leptype+"_lepetaFW"]->Fill(ev.LLep()->Eta(), weight);
			_hists1d[leptype+"_lepptFW"]->Fill(ev.LLep()->Pt(), weight);
		}
		if(ev.MET() != nullptr)
		{
			_hists2d[leptype+"_thadpt_WMT"]->Fill(ev.THad()->Pt(), sqrt(pow(ev.LLep()->Pt() + ev.MET()->Pt(), 2) - pow(ev.LLep()->Px() + ev.MET()->Px(),2) - pow(ev.LLep()->Py() + ev.MET()->Py(),2)) ,  weight);
		}	
	}


	if(ev.AddJets().size() > 0)
	{
		_hists1d[leptype+"_jet1pt"]->Fill(ev.AddJets()[0]->Pt(), weight);
	}
	if(ev.AddJets().size() > 1)
	{
		_hists1d[leptype+"_jet2pt"]->Fill(ev.AddJets()[1]->Pt(), weight);
	}
	if(ev.AddJets().size() > 2)
	{
		_hists1d[leptype+"_jet3pt"]->Fill(ev.AddJets()[2]->Pt(), weight);
	}
}

void TTPlots::Fill(TTEvent& ev, double weight)
{
	if(abs(ev.LPDGID()) == 13)
	{
		OMuon* mu = dynamic_cast<OMuon*>(ev.LLep());
		if(dynamic_cast<BoostedTop*>(GLAN)->DELPHES == false &&  mu != nullptr) {_hists1d["mu_isonear"]->Fill(mu->IsoNearAll()/mu->IsoCentralAll(), weight);}
		Fill("mu", ev, weight);
	}
	else if(abs(ev.LPDGID()) == 11)
	{
		OElectron* el = dynamic_cast<OElectron*>(ev.LLep());
		if(dynamic_cast<BoostedTop*>(GLAN)->DELPHES == false && el != nullptr) {_hists1d["el_isonear"]->Fill(el->IsoNearAll()/el->IsoCentralAll(), weight);}
		Fill("el", ev, weight);
	}
}

///////////////////////////////////////////////////////////////////7
//HadJetHists
HadJetHists::HadJetHists(const string& prefix) : Plots(prefix)
{

}

void HadJetHists::Init()
{
	TDirectory* olddir = gDirectory;
	TDirectory* newdir = gDirectory->mkdir(_prefix.c_str());
	newdir->cd();
	_hists2d.AddHist("tm_thadpt_eta0", 150, 0, 300, 150, 0, 1500, "M_{t} [GeV]", "p_{T}(t_{h}) [GeV]");
	_hists2d.AddHist("tm_thadpt_eta1", 150, 0, 300, 150, 0, 1500, "M_{t} [GeV]", "p_{T}(t_{h}) [GeV]");
	_hists2d.AddHist("tm_thadpt_eta2", 150, 0, 300, 150, 0, 1500, "M_{t} [GeV]", "p_{T}(t_{h}) [GeV]");
	_hists2d.AddHist("tm_thadpt_eta3", 150, 0, 300, 150, 0, 1500, "M_{t} [GeV]", "p_{T}(t_{h}) [GeV]");
	int ptbins = 150;
	double ptmin = 0;
	double ptmax = 1500;
	_hists2d.AddHist("ptall_mtop", ptbins, ptmin, ptmax, 100, 0., 400, "p_{T} [GeV]", "M_{top} [GeV]");
	_hists2d.AddHist("ptall_mw", ptbins, ptmin, ptmax, 100, 0., 400, "p_{T} [GeV]", "M_{W} [GeV]");
	_hists2d.AddHist("ptall_m3a", ptbins, ptmin, ptmax, 100, 0., 400, "p_{T} [GeV]", "M_{3,a} [GeV]");
	_hists2d.AddHist("ptall_m3aovermall", ptbins, ptmin, ptmax, 100, 0., 1., "p_{T} [GeV]", "M_{3,a}/M_{all}");
	_hists2d.AddHist("ptall_m3b", ptbins, ptmin, ptmax, 100, 0., 400, "p_{T} [GeV]", "M_{3,b} [GeV]");
	_hists2d.AddHist("ptall_m2a", ptbins, ptmin, ptmax, 100, 0., 400, "p_{T} [GeV]", "M_{2,a} [GeV]");
	_hists2d.AddHist("ptall_m2b", ptbins, ptmin, ptmax, 100, 0., 400, "p_{T} [GeV]", "M_{2,b} [GeV]");
	_hists2d.AddHist("ptall_m2c", ptbins, ptmin, ptmax, 100, 0., 400, "p_{T} [GeV]", "M_{2,c} [GeV]");
	_hists2d.AddHist("ptall_hardpairs", ptbins, ptmin, ptmax, 50, -0.5, 49.5, "p_{T} [GeV]", "pairs");
	_hists2d.AddHist("ptall_tau21", ptbins, ptmin, ptmax, 100, 0., 1, "p_{T} [GeV]", "#tau_{2}/#tau_{1}");
	_hists2d.AddHist("ptall_tau32", ptbins, ptmin, ptmax, 100, 0., 1, "p_{T} [GeV]", "#tau_{3}/#tau_{2}");
	_hists2d.AddHist("ptall_tau43", ptbins, ptmin, ptmax, 100, 0., 1, "p_{T} [GeV]", "#tau_{4}/#tau_{3}");
	_hists2d.AddHist("ptall_tau54", ptbins, ptmin, ptmax, 100, 0., 1, "p_{T} [GeV]", "#tau_{5}/#tau_{4}");
	_hists2d.AddHist("ptall_ja", ptbins, ptmin, ptmax, 100, 0., 400, "p_{T} [GeV]", "E(j_{a}) [GeV]");
	_hists2d.AddHist("ptall_jb", ptbins, ptmin, ptmax, 100, 0., 400, "p_{T} [GeV]", "E(j_{b}) [GeV]");
	_hists2d.AddHist("ptall_jc", ptbins, ptmin, ptmax, 100, 0., 400, "p_{T} [GeV]", "E(j_{c}) [GeV]");
	_hists2d.AddHist("ptall_jd", ptbins, ptmin, ptmax, 100, 0., 400, "p_{T} [GeV]", "E(j_{d}) [GeV]");
	_hists2d.AddHist("ptall_hardjets", ptbins, ptmin, ptmax, 20, 0., 20, "p_{T} [GeV]", "hard jets");
	_hists2d.AddHist("ptall_sphericity", ptbins, ptmin, ptmax, 100, 0., 1., "p_{T} [GeV]", "sphericity");
	_hists2d.AddHist("ptall_m2boverm2a", ptbins, ptmin, ptmax, 100, 0., 1., "p_{T} [GeV]", "M_{2,b}/M_{2,a}");
	_hists2d.AddHist("ptall_plan", ptbins, ptmin, ptmax, 100, 0., 1., "p_{T} [GeV]", "planarity");
	_hists2d.AddHist("ptall_fraca", ptbins, ptmin, ptmax, 100, 0., 1., "p_{T} [GeV]", "f_{a}");
	_hists2d.AddHist("ptall_fracb", ptbins, ptmin, ptmax, 100, 0., 1., "p_{T} [GeV]", "f_{b}");
	_hists2d.AddHist("ptall_fracc", ptbins, ptmin, ptmax, 100, 0., 1., "p_{T} [GeV]", "f_{c}");
	_hists2d.AddHist("ptall_nnres", ptbins, ptmin, ptmax, 100, 0., 1, "p_{T} [GeV]", "H_{NN}");
	_hists2d.AddHist("ptall_btag", ptbins, ptmin, ptmax, 100, 0., 1, "p_{T} [GeV]", "DeepCSV");

	_hists2d.AddHist("ptlow_eta_nnres", 48, -2.4, 2.4, 50, 0., 1, "#eta", "H_{NN}");
	_hists2d.AddHist("ptmedium_eta_nnres", 48, -2.4, 2.4, 50, 0., 1, "#eta", "H_{NN}");
	_hists2d.AddHist("pthigh_eta_nnres", 48, -2.4, 2.4, 50, 0., 1, "#eta", "H_{NN}");

	_hists2d.AddHist("bjets_pt_nnres", ptbins, ptmin, ptmax, 100, 0., 1, "p_{T}(b-jets) [GeV]", "H_{NN}");
	_hists2d.AddHist("ljets_pt_nnres", ptbins, ptmin, ptmax, 100, 0., 1, "p_{T}(l-jets) [GeV]", "H_{NN}");
	_hists2d.AddHist("bjets_pt_btag", ptbins, ptmin, ptmax, 100, 0., 1, "p_{T}(b-jets) [GeV]", "DeepCSV");
	_hists2d.AddHist("ljets_pt_btag", ptbins, ptmin, ptmax, 100, 0., 1, "p_{T}(l-jets) [GeV]", "DeepCSV");

	int npvbins = 40;
	double npvmin = 0;
	double npvmax = 80;
	_hists2d.AddHist("nvertices_ptall", npvbins, npvmin, npvmax, ptbins, ptmin, ptmax, "# vertices", "p_{T} [GeV]");
	_hists2d.AddHist("nvertices_mtop", npvbins, npvmin, npvmax, 100, 0., 400, "# vertices", "M_{top} [GeV]");
	_hists2d.AddHist("nvertices_mw", npvbins, npvmin, npvmax, 100, 0., 400, "# vertices", "M_{W} [GeV]");
	_hists2d.AddHist("nvertices_m3a", npvbins, npvmin, npvmax, 100, 0., 400, "# vertices", "M_{3,a} [GeV]");
	_hists2d.AddHist("nvertices_m3aovermall", npvbins, npvmin, npvmax, 100, 0., 1, "# vertices", "M_{3,a}/M_{all}");
	_hists2d.AddHist("nvertices_m3b", npvbins, npvmin, npvmax, 100, 0., 400, "# vertices", "M_{3,b} [GeV]");
	_hists2d.AddHist("nvertices_m2a", npvbins, npvmin, npvmax, 100, 0., 400, "# vertices", "M_{2,a} [GeV]");
	_hists2d.AddHist("nvertices_m2b", npvbins, npvmin, npvmax, 100, 0., 400, "# vertices", "M_{2,b} [GeV]");
	_hists2d.AddHist("nvertices_m2c", npvbins, npvmin, npvmax, 100, 0., 400, "# vertices", "M_{2,c} [GeV]");
	_hists2d.AddHist("nvertices_hardpairs", npvbins, npvmin, npvmax, 50, -0.5, 49.5, "# vertices", "pairs");
	_hists2d.AddHist("nvertices_tau21", npvbins, npvmin, npvmax, 100, 0., 1, "# vertices", "#tau_{2}/#tau_{1}");
	_hists2d.AddHist("nvertices_tau32", npvbins, npvmin, npvmax, 100, 0., 1, "# vertices", "#tau_{3}/#tau_{2}");
	_hists2d.AddHist("nvertices_tau43", npvbins, npvmin, npvmax, 100, 0., 1, "# vertices", "#tau_{4}/#tau_{3}");
	_hists2d.AddHist("nvertices_tau54", npvbins, npvmin, npvmax, 100, 0., 1, "# vertices", "#tau_{5}/#tau_{4}");
	_hists2d.AddHist("nvertices_ja", npvbins, npvmin, npvmax, 100, 0., 400, "# vertices", "E(j_{a}) [GeV]");
	_hists2d.AddHist("nvertices_jb", npvbins, npvmin, npvmax, 100, 0., 400, "# vertices", "E(j_{b}) [GeV]");
	_hists2d.AddHist("nvertices_jc", npvbins, npvmin, npvmax, 100, 0., 400, "# vertices", "E(j_{c}) [GeV]");
	_hists2d.AddHist("nvertices_jd", npvbins, npvmin, npvmax, 100, 0., 400, "# vertices", "E(j_{d}) [GeV]");
	_hists2d.AddHist("nvertices_hardjets", npvbins, npvmin, npvmax, 20, 0., 20, "p_{T} [GeV]", "hard jets");
	_hists2d.AddHist("nvertices_sphericity", npvbins, npvmin, npvmax, 100, 0., 1., "p_{T} [GeV]", "sphericity");
	_hists2d.AddHist("nvertices_m2boverm2a", npvbins, npvmin, npvmax, 100, 0., 1., "# vertices", "M_{2,b}/M_{2,a}");
	_hists2d.AddHist("nvertices_plan", npvbins, npvmin, npvmax, 100, 0., 1., "# vertices", "planarity");
	_hists2d.AddHist("nvertices_nnres", npvbins, npvmin, npvmax, 100, 0., 1., "# vertices", "H_{NN}");
	_hists2d.AddHist("nvertices_btag", npvbins, npvmin, npvmax, 100, 0., 1., "# vertices", "DeepCSV");

	int etabins = 25;
	double etamin = 0;
	double etamax = 2.5;
	_hists2d.AddHist("etaall_mtop", etabins, etamin, etamax, 100, 0., 400, "|#eta|", "M_{top} [GeV]");
	_hists2d.AddHist("etaall_mw", etabins, etamin, etamax, 100, 0., 400, "|#eta|", "M_{W} [GeV]");
	_hists2d.AddHist("etaall_m3a", etabins, etamin, etamax, 100, 0., 400, "|#eta|", "M_{3,a} [GeV]");
	_hists2d.AddHist("etaall_m3aovermall", etabins, etamin, etamax, 100, 0., 1., "|#eta|", "M_{3,a}/M_{all}");
	_hists2d.AddHist("etaall_m3b", etabins, etamin, etamax, 100, 0., 400, "|#eta|", "M_{3,b} [GeV]");
	_hists2d.AddHist("etaall_m2a", etabins, etamin, etamax, 100, 0., 400, "|#eta|", "M_{2,a} [GeV]");
	_hists2d.AddHist("etaall_m2b", etabins, etamin, etamax, 100, 0., 400, "|#eta|", "M_{2,b} [GeV]");
	_hists2d.AddHist("etaall_m2c", etabins, etamin, etamax, 100, 0., 400, "|#eta|", "M_{2,c} [GeV]");
	_hists2d.AddHist("etaall_hardpairs", etabins, etamin, etamax, 50, -0.5, 49.5, "|#eta|", "pairs");
	_hists2d.AddHist("etaall_tau21", etabins, etamin, etamax, 100, 0., 1, "|#eta|", "#tau_{2}/#tau_{1}");
	_hists2d.AddHist("etaall_tau32", etabins, etamin, etamax, 100, 0., 1, "|#eta|", "#tau_{3}/#tau_{2}");
	_hists2d.AddHist("etaall_tau43", etabins, etamin, etamax, 100, 0., 1, "|#eta|", "#tau_{4}/#tau_{3}");
	_hists2d.AddHist("etaall_tau54", etabins, etamin, etamax, 100, 0., 1, "|#eta|", "#tau_{5}/#tau_{4}");
	_hists2d.AddHist("etaall_ja", etabins, etamin, etamax, 100, 0., 400, "|#eta|", "E(j_{a}) [GeV]");
	_hists2d.AddHist("etaall_jb", etabins, etamin, etamax, 100, 0., 400, "|#eta|", "E(j_{b}) [GeV]");
	_hists2d.AddHist("etaall_jc", etabins, etamin, etamax, 100, 0., 400, "|#eta|", "E(j_{c}) [GeV]");
	_hists2d.AddHist("etaall_jd", etabins, etamin, etamax, 100, 0., 400, "|#eta|", "E(j_{d}) [GeV]");
	_hists2d.AddHist("etaall_hardjets", etabins, etamin, etamax, 20, 0., 20, "|#eta|", "hard jets");
	_hists2d.AddHist("etaall_sphericity", etabins, etamin, etamax, 100, 0., 1., "|#eta|", "sphericity");
	_hists2d.AddHist("etaall_m2boverm2a", etabins, etamin, etamax, 100, 0., 1., "|#eta|", "M_{2,b}/M_{2,a}");
	_hists2d.AddHist("etaall_plan", etabins, etamin, etamax, 100, 0., 1., "|#eta|", "planarity");
	_hists2d.AddHist("etaall_fraca", etabins, etamin, etamax, 100, 0., 1., "|#eta|", "f_{a}");
	_hists2d.AddHist("etaall_fracb", etabins, etamin, etamax, 100, 0., 1., "|#eta|", "f_{b}");
	_hists2d.AddHist("etaall_fracc", etabins, etamin, etamax, 100, 0., 1., "|#eta|", "f_{c}");
	_hists2d.AddHist("etaall_nnres", etabins, etamin, etamax, 100, 0., 1, "|#eta|", "H_{NN}");
	_hists2d.AddHist("etaall_btag", etabins, etamin, etamax, 100, 0., 1, "|#eta|", "DeepCSV");
	olddir->cd();
}

void HadJetHists::Fill(HadTopJet* tj, double weight)
{
	BoostedTop* BT = dynamic_cast<BoostedTop*>(GLAN);
	int npv = BT->NumGoodVertices();
	double ptall = tj->Pt();
	double etaall = abs(tj->Eta());

	//For JECS
	if(tj->Eta() < -1.) {_hists2d["tm_thadpt_eta0"]->Fill(tj->M(), tj->Pt(), weight);}
	else if(tj->Eta() < 0.) {_hists2d["tm_thadpt_eta1"]->Fill(tj->M(), tj->Pt(), weight);}
	else if(tj->Eta() < 1.) {_hists2d["tm_thadpt_eta2"]->Fill(tj->M(), tj->Pt(), weight);}
	else {_hists2d["tm_thadpt_eta3"]->Fill(tj->M(), tj->Pt(), weight);}

	TVector3 cross(SJet(tj->SubJets(0)).Vect().Cross(SJet(tj->SubJets(1)).Vect()));
	cross *= 1./cross.Mag();
	double plan = abs(cross * SJet(tj->SubJets(2)).Vect())/SJet(tj->SubJets(2)).P();
	_hists2d["nvertices_ptall"]->Fill(npv, ptall, weight);
	_hists2d["nvertices_mtop"]->Fill(npv, tj->T().M(), weight);
	_hists2d["nvertices_mw"]->Fill(npv, tj->W().M(), weight);
	_hists2d["nvertices_m3a"]->Fill(npv, tj->M3s()[0], weight);
	_hists2d["nvertices_m3aovermall"]->Fill(npv, tj->M3s()[0]/tj->M(), weight);
	if(tj->M3s().size() > 1) {_hists2d["nvertices_m3b"]->Fill(npv, tj->M3s().size() > 1 ? tj->M3s()[1] : 0., weight);}
	_hists2d["nvertices_m2a"]->Fill(npv, tj->M2s()[0], weight);
	_hists2d["nvertices_m2b"]->Fill(npv, tj->M2s()[1], weight);
	_hists2d["nvertices_m2c"]->Fill(npv, tj->M2s()[2], weight);
	_hists2d["nvertices_hardpairs"]->Fill(npv, tj->HardPairs(), weight);
	_hists2d["nvertices_tau21"]->Fill(npv, tj->Taus(1)/tj->Taus(0), weight);
	_hists2d["nvertices_tau32"]->Fill(npv, tj->Taus(2)/tj->Taus(1), weight);
	_hists2d["nvertices_tau43"]->Fill(npv, tj->Num_Taus() > 3 ? tj->Taus(3)/tj->Taus(2) : 0., weight);
	_hists2d["nvertices_tau54"]->Fill(npv, tj->Num_Taus() > 4 ? tj->Taus(4)/tj->Taus(3) : 0., weight);
	_hists2d["nvertices_ja"]->Fill(npv, SJet(tj->SubJets(0)).E(), weight);
	_hists2d["nvertices_jb"]->Fill(npv, SJet(tj->SubJets(1)).E(), weight);
	_hists2d["nvertices_jc"]->Fill(npv, SJet(tj->SubJets(2)).E(), weight);
	_hists2d["nvertices_jd"]->Fill(npv, tj->Num_SubJets() > 3 ? SJet(tj->SubJets(3)).E() : 0., weight);
	_hists2d["nvertices_hardjets"]->Fill(npv, tj->HardJets(), weight);
	_hists2d["nvertices_sphericity"]->Fill(npv, tj->Sphericity(), weight);
	_hists2d["nvertices_m2boverm2a"]->Fill(npv, tj->M2s()[1]/tj->M2s()[0], weight);
	_hists2d["nvertices_plan"]->Fill(npv, plan, weight);
	_hists2d["nvertices_nnres"]->Fill(npv, tj->NNRes(), weight);
	_hists2d["nvertices_btag"]->Fill(npv, tj->BTag(), weight);

	_hists2d["ptall_mtop"]->Fill(ptall, tj->T().M(), weight);
	_hists2d["ptall_mw"]->Fill(ptall, tj->W().M(), weight);
	_hists2d["ptall_m3a"]->Fill(ptall, tj->M3s()[0], weight);
	_hists2d["ptall_m3aovermall"]->Fill(ptall, tj->M3s()[0]/tj->M(), weight);
	if(tj->M3s().size() > 1) {_hists2d["ptall_m3b"]->Fill(ptall, tj->M3s().size() > 1 ? tj->M3s()[1] : 0., weight);}
	_hists2d["ptall_m2a"]->Fill(ptall, tj->M2s()[0], weight);
	_hists2d["ptall_m2b"]->Fill(ptall, tj->M2s()[1], weight);
	_hists2d["ptall_m2c"]->Fill(ptall, tj->M2s()[2], weight);
	_hists2d["ptall_hardpairs"]->Fill(ptall, tj->HardPairs(), weight);
	_hists2d["ptall_tau21"]->Fill(ptall, tj->Taus(1)/tj->Taus(0), weight);
	_hists2d["ptall_tau32"]->Fill(ptall, tj->Taus(2)/tj->Taus(1), weight);
	_hists2d["ptall_tau43"]->Fill(ptall, tj->Num_Taus() > 3 ? tj->Taus(3)/tj->Taus(2) : 0., weight);
	_hists2d["ptall_tau54"]->Fill(ptall, tj->Num_Taus() > 4 ? tj->Taus(4)/tj->Taus(3) : 0., weight);
	_hists2d["ptall_ja"]->Fill(ptall, SJet(tj->SubJets(0)).E(), weight);
	_hists2d["ptall_jb"]->Fill(ptall, SJet(tj->SubJets(1)).E(), weight);
	_hists2d["ptall_jc"]->Fill(ptall, SJet(tj->SubJets(2)).E(), weight);
	_hists2d["ptall_jd"]->Fill(ptall, tj->Num_SubJets() > 3 ? SJet(tj->SubJets(3)).E() : 0., weight);
	_hists2d["ptall_hardjets"]->Fill(ptall, tj->HardJets(), weight);
	_hists2d["ptall_sphericity"]->Fill(ptall, tj->Sphericity(), weight);
	_hists2d["ptall_m2boverm2a"]->Fill(ptall, tj->M2s()[1]/tj->M2s()[0], weight);
	_hists2d["ptall_plan"]->Fill(ptall, plan, weight);
	_hists2d["ptall_fraca"]->Fill(ptall, tj->Frac(0), weight);
	_hists2d["ptall_fracb"]->Fill(ptall, tj->Frac(1), weight);
	_hists2d["ptall_fracc"]->Fill(ptall, tj->Frac(2), weight);
	_hists2d["ptall_nnres"]->Fill(ptall, tj->NNRes(), weight);
	_hists2d["ptall_btag"]->Fill(ptall, tj->BTag(), weight);

	_hists2d["etaall_mtop"]->Fill(etaall, tj->T().M(), weight);
	_hists2d["etaall_mw"]->Fill(etaall, tj->W().M(), weight);
	_hists2d["etaall_m3a"]->Fill(etaall, tj->M3s()[0], weight);
	_hists2d["etaall_m3aovermall"]->Fill(etaall, tj->M3s()[0]/tj->M(), weight);
	if(tj->M3s().size() > 1) {_hists2d["etaall_m3b"]->Fill(etaall, tj->M3s().size() > 1 ? tj->M3s()[1] : 0., weight);}
	_hists2d["etaall_m2a"]->Fill(etaall, tj->M2s()[0], weight);
	_hists2d["etaall_m2b"]->Fill(etaall, tj->M2s()[1], weight);
	_hists2d["etaall_m2c"]->Fill(etaall, tj->M2s()[2], weight);
	_hists2d["etaall_hardpairs"]->Fill(etaall, tj->HardPairs(), weight);
	_hists2d["etaall_tau21"]->Fill(etaall, tj->Taus(1)/tj->Taus(0), weight);
	_hists2d["etaall_tau32"]->Fill(etaall, tj->Taus(2)/tj->Taus(1), weight);
	_hists2d["etaall_tau43"]->Fill(etaall, tj->Num_Taus() > 3 ? tj->Taus(3)/tj->Taus(2) : 0., weight);
	_hists2d["etaall_tau54"]->Fill(etaall, tj->Num_Taus() > 4 ? tj->Taus(4)/tj->Taus(3) : 0., weight);
	_hists2d["etaall_ja"]->Fill(etaall, SJet(tj->SubJets(0)).E(), weight);
	_hists2d["etaall_jb"]->Fill(etaall, SJet(tj->SubJets(1)).E(), weight);
	_hists2d["etaall_jc"]->Fill(etaall, SJet(tj->SubJets(2)).E(), weight);
	_hists2d["etaall_jd"]->Fill(etaall, tj->Num_SubJets() > 3 ? SJet(tj->SubJets(3)).E() : 0., weight);
	_hists2d["etaall_hardjets"]->Fill(etaall, tj->HardJets(), weight);
	_hists2d["etaall_sphericity"]->Fill(etaall, tj->Sphericity(), weight);
	_hists2d["etaall_m2boverm2a"]->Fill(etaall, tj->M2s()[1]/tj->M2s()[0], weight);
	_hists2d["etaall_plan"]->Fill(etaall, plan, weight);
	_hists2d["etaall_fraca"]->Fill(etaall, tj->Frac(0), weight);
	_hists2d["etaall_fracb"]->Fill(etaall, tj->Frac(1), weight);
	_hists2d["etaall_fracc"]->Fill(etaall, tj->Frac(2), weight);
	_hists2d["etaall_nnres"]->Fill(etaall, tj->NNRes(), weight);
	_hists2d["etaall_btag"]->Fill(etaall, tj->BTag(), weight);


	if(ptall < 500 ) {_hists2d["ptlow_eta_nnres"]->Fill(tj->Eta(), tj->NNRes(), weight);}
	else if(ptall < 700 ) {_hists2d["ptmedium_eta_nnres"]->Fill(tj->Eta(), tj->NNRes(), weight);}
	else {_hists2d["pthigh_eta_nnres"]->Fill(tj->Eta(), tj->NNRes(), weight);}

	auto gjb = find_if(BT->GenBJets.begin(), BT->GenBJets.end(), [&](TLorentzVector* bp){return tj->CloseToMember(bp);});
	if(gjb != BT->GenBJets.end())
	{
		_hists2d["bjets_pt_nnres"]->Fill(tj->Pt(), tj->NNRes(), weight);
		_hists2d["bjets_pt_btag"]->Fill(tj->Pt(), tj->BTag(), weight);
	}
	else
	{
		_hists2d["ljets_pt_nnres"]->Fill(tj->Pt(), tj->NNRes(), weight);
		_hists2d["ljets_pt_btag"]->Fill(tj->Pt(), tj->BTag(), weight);

	}
}

///////////////////////////////////////////////////////////////////7
//LepJetHists
LepJetHists::LepJetHists(const string& prefix) : Plots(prefix)
{

}

void LepJetHists::Init()
{
	TDirectory* olddir = gDirectory;
	TDirectory* newdir = gDirectory->mkdir(_prefix.c_str());
	newdir->cd();
	int ptbins = 200;
	double ptmin = 0;
	double ptmax = 2000;
	_hists2d.AddHist("Mt_Mw", 500, 0, 500, 500, 0, 500, "", "");
	_hists2d.AddHist("met_lep_pt", 500, 0, 500, 500, 0, 500, "", "");
	_hists2d.AddHist("lepfrac_metfrac", 100, 0,1, 100, 0, 1, "", "");
	_hists1d.AddHist("dphimetlep", 100, 0,3, "", "");
	_hists1d.AddHist("Mtmetlep", 100, 0,200, "", "");
	_hists2d.AddHist("jetpt_lepinjet", ptbins, ptmin, ptmax, 2, 0, 2,"", "");
	_hists2d.AddHist("jetpt_jetmass", ptbins, ptmin, ptmax, 150, 0, 300,"", "");
	_hists2d.AddHist("jetpt_jetredmass", ptbins, ptmin, ptmax, 100, 0, 1,"", "");
	_hists2d.AddHist("jetpt_lepfrac", ptbins, ptmin, ptmax, 100, 0, 1,"", "");
	_hists2d.AddHist("jetpt_metfrac", ptbins, ptmin, ptmax, 100, 0, 1,"", "");
	_hists2d.AddHist("jetpt_lepdrin", ptbins, ptmin, ptmax, 100, 0, 2,"", "");
	_hists2d.AddHist("jetpt_lepdrout", ptbins, ptmin, ptmax, 100, 0, 2,"", "");
	_hists2d.AddHist("jetpt_nnres", ptbins, ptmin, ptmax, 100, 0, 1, "L_{NN}", "");
	_hists2d.AddHist("jetpt_iso", ptbins, ptmin, ptmax, 100, 0, 1, "", "");
	_hists2d.AddHist("jetpt_isofar", ptbins, ptmin, ptmax, 100, 0, 0.2, "", "");
	_hists2d.AddHist("jetpt_isonear", ptbins, ptmin, ptmax, 100, 0, 200, "", "");
	_hists2d.AddHist("jetpt_isochargedfar", ptbins, ptmin, ptmax, 100, 0, 0.2, "", "");
	_hists2d.AddHist("jetpt_isochargednear", ptbins, ptmin, ptmax, 100, 0, 400, "", "");
	olddir->cd();
}

void LepJetHists::Fill(LepTopJet* tj, double weight)
{
	TLorentzVector* lep = tj->Lep();
	TLorentzVector* met = tj->MET();
	_hists2d["Mt_Mw"]->Fill(tj->M(), (*met+*lep).M(), weight);
	_hists2d["met_lep_pt"]->Fill(met->Pt(), lep->Pt(), weight);
	_hists2d["lepfrac_metfrac"]->Fill(lep->Pt()/tj->Pt(), met->Pt()/tj->Pt(), weight);
	_hists1d["dphimetlep"]->Fill(abs(lep->DeltaPhi(*met)), weight);
	_hists1d["Mtmetlep"]->Fill(sqrt(lep->Pt()*met->Pt()-lep->Px()*met->Px() - lep->Py()*met->Py()), weight);
	_hists2d["jetpt_jetmass"]->Fill(tj->Pt(), tj->LepJet().M(), weight);
	_hists2d["jetpt_jetredmass"]->Fill(tj->Pt(), (tj->LepJet() - *lep).M()/tj->LepJet().M(), weight);
	_hists2d["jetpt_lepfrac"]->Fill(tj->Pt(), lep->Pt()/tj->LepJet().Pt(), weight);
	_hists2d["jetpt_metfrac"]->Fill(tj->Pt(), met->Pt()/tj->Pt(), weight);
	if(tj->LepInJet()) {_hists2d["jetpt_lepdrin"]->Fill(tj->Pt(), tj->Jet()->DeltaR(*lep), weight);}
	else{_hists2d["jetpt_lepdrout"]->Fill(tj->Pt(), tj->Jet()->DeltaR(*lep), weight);}
	_hists2d["jetpt_nnres"]->Fill(tj->Pt(), tj->NNRes(), weight);
	if(abs(tj->LepPDGID()) == 13)
	{
		OMuon* mu = dynamic_cast<OMuon*>(lep);
		_hists2d["jetpt_lepinjet"]->Fill(tj->Pt(), tj->LepInJet() ? 1.5 : 0.5, weight);
		_hists2d["jetpt_iso"]->Fill(tj->Pt(), mu->PFIsolationDB(), weight);
		_hists2d["jetpt_isofar"]->Fill(tj->Pt(), mu->IsoFarAll()/mu->IsoCentralAll(), weight);
		_hists2d["jetpt_isonear"]->Fill(tj->Pt(), mu->IsoNearAll()/mu->IsoCentralAll(), weight);
		_hists2d["jetpt_isochargedfar"]->Fill(tj->Pt(), mu->IsoFar().Charged()/mu->IsoCentral().Charged(), weight);
		_hists2d["jetpt_isochargednear"]->Fill(tj->Pt(), mu->IsoNear().Charged()/mu->IsoCentral().Charged(), weight);
	}
	if(abs(tj->LepPDGID()) == 11)
	{
		OElectron* el = dynamic_cast<OElectron*>(lep);
		_hists2d["jetpt_lepinjet"]->Fill(tj->Pt(), tj->LepInJet() ? 1.5 : 0.5, weight);
		_hists2d["jetpt_iso"]->Fill(tj->Pt(), el->CorPFIsolation2016(), weight);
		_hists2d["jetpt_isofar"]->Fill(tj->Pt(), el->IsoFarAll()/el->IsoCentralAll(), weight);
		_hists2d["jetpt_isonear"]->Fill(tj->Pt(), el->IsoNearAll()/el->IsoCentralAll(), weight);
		_hists2d["jetpt_isochargedfar"]->Fill(tj->Pt(), el->IsoFar().Charged()/el->IsoCentral().Charged(), weight);
		_hists2d["jetpt_isochargednear"]->Fill(tj->Pt(), el->IsoNear().Charged()/el->IsoCentral().Charged(), weight);
	}
}

