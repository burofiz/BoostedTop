#include "BtagEff.h"
#include <TTree.h>
#include <TRandom3.h>
#include <iostream>

#include <algorithm>

#include "TTEvent.h"
#include "OJet.h"

#include "BoostedTop.h"

BtagEff::BtagEff(){}

void BtagEff::Init(double btagsel, double bkgcutmin, double bkgcutmax)
{
	bkgcutmin_ = bkgcutmin;
	bkgcutmax_ = bkgcutmax;
	btagselection = btagsel;
	btagtree = new TTree("btagtree", "btagtree", 1);
	btagtree->Branch("j", j, "j[17]/F");
	btagtree->Branch("prob", &prob, "prob/F");
	btagtree->Branch("problep", &problep, "problep/F");
	btagtree->Branch("probhad", &probhad, "probhad/F");
	btagtree->Branch("met", &met, "met/F");
	btagtree->Branch("weight", &weight, "weight/F");
	btagtree->Branch("nvtx", &nvtx, "nvtx/F");
	btagtree->Branch("run", &run, "run/I");
	btagtree->Branch("typ", &typ, "typ/I");
	btagtree->Branch("test", &test, "test/I");

	_hists1d.AddHist("prob_pass_sig", 50, 0, 50, "prob", "Events");
	_hists1d.AddHist("prob_fail_sig", 50, 0, 50, "prob", "Events");
	_hists1d.AddHist("prob_pass_bkg", 50, 0, 50, "prob", "Events");
	_hists1d.AddHist("prob_fail_bkg", 50, 0, 50, "prob", "Events");
}

void BtagEff::PrepareFill(OJet* jet)
{
	j[0] = jet->Px();
	j[1] = jet->Py();
	j[2] = jet->Pz();
	j[3] = jet->E();

	j[4] = jet->FlavorInfo(0).CSV();
	j[5] = jet->FlavorInfo(0).DeepCSVbb();
	j[6] = jet->FlavorInfo(0).DeepCSVb();
	j[7] = jet->FlavorInfo(0).DeepCSVc();
	j[8] = jet->FlavorInfo(0).DeepCSVudsg();
	j[9] = jet->FlavorInfo(0).DeepFlavourbb();
	j[10] = jet->FlavorInfo(0).DeepFlavourb();
	j[11] = jet->FlavorInfo(0).DeepFlavourlepb();
	j[12] = jet->FlavorInfo(0).DeepFlavourc();
	j[13] = jet->FlavorInfo(0).DeepFlavourg();
	j[14] = jet->FlavorInfo(0).DeepFlavouruds();
	j[15] = jet->FlavorInfo(0).ComCvsL();
	j[16] = jet->FlavorInfo(0).ComCvsB();


//	j[4] = jet->FlavorInfo(0).CSV();
//	j[5] = jet->FlavorInfo(0).DeepCSVb()+jet->FlavorInfo(0).DeepCSVbb(); //DeepCSV
//	j[6] = 0.;
//	if(jet->FlavorInfo(0).ComCvsL() > -0.53 &&  jet->FlavorInfo(0).ComCvsB() > -0.26) {j[6]=1;}
//	j[7] = 0.;
//	if(jet->FlavorInfo(0).ComCvsL() > 0.07 &&  jet->FlavorInfo(0).ComCvsB() > -0.10) {j[7]=1;}
//	j[8] = 0.;
//	if(jet->FlavorInfo(0).ComCvsL() > 0.87 &&  jet->FlavorInfo(0).ComCvsB() > -0.30) {j[8]=1;}
//
//	double csvb = jet->FlavorInfo(0).DeepCSVc() <= -1 ? -10. : jet->FlavorInfo(0).DeepCSVc() / (jet->FlavorInfo(0).DeepCSVc() + jet->FlavorInfo(0).DeepCSVb() + jet->FlavorInfo(0).DeepCSVbb());
//	double csvl = jet->FlavorInfo(0).DeepCSVc() <= -1 ? -10. : jet->FlavorInfo(0).DeepCSVc() / (jet->FlavorInfo(0).DeepCSVc() + jet->FlavorInfo(0).DeepCSVudsg());
//	j[9] = 0.;
//	if(csvl > 0.05 &&  csvb > 0.33) {j[9]=1;}
//	j[10] = 0.;
//	if(csvl > 0.15 &&  csvb > 0.28) {j[10]=1;}
//	j[11] = 0.;
//	if(csvl > 0.8 &&  csvb > 0.1) {j[11]=1;}
}

bool BtagEff::Fill(TTEvent& per, TTEvent& genper, double theweight)
{
	bool use = false;
	if(per.BHad()->DeltaR(*per.JaHad()) < 0.8 || per.BHad()->DeltaR(*per.JbHad()) < 0.8 || per.BHad()->DeltaR(*per.BLep()) < 0.8 || per.BHad()->DeltaR(*per.LLep()) < 0.5 || per.BLep()->DeltaR(*per.LLep()) < 0.5) {return use;}
	//if(per.NAddJets() > 0 && per.THad()->Pt() < 200 && per.TLep()->Pt() < 200) {return use;}
	if(per.NAddJets() > 0) {return use;}

	OJet* blep = dynamic_cast<OJet*>(per.BLep());
	OJet* bhad = dynamic_cast<OJet*>(per.BHad());
	OJet* wja = dynamic_cast<OJet*>(per.JaHad());
	OJet* wjb = dynamic_cast<OJet*>(per.JbHad());

	typ = 0;
	//if(genper.IsComplete())
	{
		BoostedTop* BT = dynamic_cast<BoostedTop*>(GLAN);

		auto genjetbhad = find_if(BT->GenBJets.begin(), BT->GenBJets.end(), [&](TLorentzVector* bp){return bhad->DeltaR(*bp) < 0.3;});
        if(genjetbhad != BT->GenBJets.end())
        { typ |= 1<<0; }
		auto genjetblep = find_if(BT->GenBJets.begin(), BT->GenBJets.end(), [&](TLorentzVector* bp){return blep->DeltaR(*bp) < 0.3;});
        if(genjetblep != BT->GenBJets.end())
        { typ |= 1<<1; }
		auto genjetchad = find_if(BT->GenCJets.begin(), BT->GenCJets.end(), [&](TLorentzVector* bp){return bhad->DeltaR(*bp) < 0.3;});
        if(genjetchad != BT->GenCJets.end())
        { typ |= 1<<2; }
		auto genjetclep = find_if(BT->GenCJets.begin(), BT->GenCJets.end(), [&](TLorentzVector* bp){return blep->DeltaR(*bp) < 0.3;});
        if(genjetclep != BT->GenCJets.end())
        { typ |= 1<<3; }

		//if(per.IsTHadCorrectKin(genper)) {typ |= 1<<0;}
		//if(per.IsTLepCorrectKin(genper)) {typ |= 1<<1;}
		//if(per.IsCorrectKin(genper)) typ = 3;
		//if(per.BHad()->DeltaR(*genper.BHad())<0.2) typ |= 1<<0;
		//if(per.BLep()->DeltaR(*genper.BLep())<0.2) typ |= 1<<1;
		//if(per.IsTLepCorrectKin(genper)) typ = 2;
		//if(per.IsCorrectKin(genper)) typ = 3;
	}
	weight = theweight;
	run = GLAN->Run();
	prob = per.Prob();
	problep = per.ProbLep();
	probhad = per.ProbHad();
	met = per.MET()->Pt();

	if(blep->BTag() > btagselection)
	{
		test = 1;
		PrepareFill(bhad);
		btagtree->Fill();
		use = true;
	}

	if(bhad->BTag() > btagselection)
	{
		test = 2;
		PrepareFill(blep);
		btagtree->Fill();
		use = true;
	}

//	if(bhad->BTag() > btagselection && blep->BTag() > btagselection)
//	{
//		test = 3;
//		PrepareFill(wja);
//		btagtree->Fill();
//		PrepareFill(wjb);
//		btagtree->Fill();
//	}

	double bmax = max({bhad->BTag(), blep->BTag(), wja->BTag(), wjb->BTag()});
	if(bmax > bkgcutmin_ && bmax < bkgcutmax_)
	{
		test = -1;
		PrepareFill(bhad);
		btagtree->Fill();
		test = -2;
		PrepareFill(blep);
		btagtree->Fill();
		//test = -3;
		//PrepareFill(wja);
		//btagtree->Fill();
		//PrepareFill(wjb);
		//btagtree->Fill();
	}
	return use;
}

