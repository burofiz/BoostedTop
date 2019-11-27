#include "BoostedTop.h"

void BoostedTop::SelectPseudoTopBoosted()
{
    int lc = 0;
    nu = TLorentzVector(0.,0.,0.,0.);
    vector<GenSimpleParticle*> pstbjets;
    vector<GenSimpleParticle*> pstljets;
    vector<GenSimpleParticle*> promptleps;
    vector<GenSimpleParticle*> nonisoleps;
    vector<GenSimpleParticle*> tlepboosted;
    vector<GenSimpleParticle*> thadboosted;
    vector<GenSimpleParticle*> promptnus;
	for(size_t n = 0 ; n < NumPLObjects() ; ++n)
    {
		GenSimpleParticle pl(GetPLObject(n));	
        if(abs(pl.PDGID()) < 6) continue;

        if(pl.PDGID() == 22)
        {
            if(pl.IsoR3() < 0.25 && abs(pl.Eta()) < 2.4 && pl.Pt() > 30.)
            {
                SBasicParticles.push_back(pl);
                pstphotons.push_back(&(SBasicParticles.back()));
            }
        }
        else if(abs(pl.PDGID()) == 13 || abs(pl.PDGID()) == 11)
        {
			if(((int)pl.Mother() & (1<<0)) == 0) //not from Hadron
            {
                if(abs(pl.Eta()) < 2.4 && pl.Pt() > 15)
                {
                    lc++;
                }
                if(abs(pl.Eta()) < 2.4 && pl.Pt() > 30.)
                {
                    SBasicParticles.push_back(pl);
					promptleps.push_back(&(SBasicParticles.back()));
                }
            }
            //if(pl.IsoR3() < 0.35)
            //{
            //    if(abs(pl.Eta()) < 2.4 && pl.Pt() > 15)
            //    {
            //        lc++;
            //    }
            //    if(abs(pl.Eta()) < 2.4 && pl.Pt() > 30.)
            //    {
            //        SBasicParticles.push_back(pl);
			//		isoleps.push_back(&(SBasicParticles.back()));
            //    }
            //}
			//else
			//{
			//	if(abs(pl.Eta()) < 2.4 && pl.Pt() > 50.)
			//	{
			//		if(((int)pl.Mother() & (1<<0)) == 0) //not from Hadron
			//		{
			//			SBasicParticles.push_back(pl);
			//			nonisoleps.push_back(&(SBasicParticles.back()));
			//		}
			//	}
			//}
		}
        else if(abs(pl.PDGID()) == 12 || abs(pl.PDGID()) == 14 || abs(pl.PDGID()) == 16)
        {
			if(pl.Mother() == 0)
			{
				SBasicParticles.push_back(pl);
				promptnus.push_back(&(SBasicParticles.back()));
			}
            nu += pl;
        }
	}
    if(lc >= 2) {return;}

	for(size_t n = 0 ; n < NumPLObjects() ; ++n)
    {
		GenSimpleParticle pl(GetPLObject(n));	
		if((abs(pl.PDGID()) == 611 || abs(pl.PDGID()) == 613) && abs(pl.Eta()) < 2.4 && pl.Pt() > 400.)
		{
			for(TLorentzVector* nil : promptleps)
			{
				if(nil->Pt() > 50. && nil->DeltaR(pl) < 0.6)
				{
					SBasicParticles.push_back(pl);
					tlepboosted.push_back(&(SBasicParticles.back()));
					break;
				}
			}
		}
		else if(abs(pl.PDGID()) == 605 && abs(pl.Eta()) < 2.4 && pl.Pt() > 350. && pl.M() > 120.)
		{
			if(find_if(promptleps.begin(), promptleps.end(), [&](const TLorentzVector* nil){return nil->DeltaR(pl) < 0.8;}) == promptleps.end())
			{
				SBasicParticles.push_back(pl);
				thadboosted.push_back(&(SBasicParticles.back()));
			}
		}
	}

	for(size_t n = 0 ; n < NumPLObjects() ; ++n)
    {
		GenSimpleParticle pl(GetPLObject(n));	
        if(abs(pl.PDGID()) > 6) {continue;}
        if(pl.Pt() < 25 || abs(pl.Eta()) > 2.4) {continue;}
        bool islep = false;
        for(TLorentzVector* lep : promptleps) {if(lep->DeltaR(pl) < 0.4) {islep = true;}}
        if(islep) {continue;}
        bool isph = false;
        for(TLorentzVector* ph : pstphotons) {if(ph->DeltaR(pl) < 0.4) {isph = true;}}
        if(isph) {continue;}
        SBasicParticles.push_back(pl);
        if(abs(pl.PDGID()) == 5)
        {
            pstbjets.push_back(&(SBasicParticles.back()));
        }
        else
        {
            pstljets.push_back(&(SBasicParticles.back()));
        }
    }

	double MW = 80.4;
	double Mt = 172.5;
	double bestchi2 = -1.;
	if(promptleps.size() == 1)
	{
		GenSimpleParticle* lepton = promptleps[0];
		for(TLorentzVector* th : thadboosted)
		{
			for(TLorentzVector* bl : pstbjets)
			{
				if(th->DeltaR(*bl) < 1.2 || th->DeltaR(*lepton) < 1.2) continue;
				TLorentzVector tl(*bl + *lepton + nu);
				double chi2 = pow(tl.M()-Mt, 2) + pow(th->M()-Mt, 2);
				if(bestchi2 < 0 || chi2 < bestchi2)
				{
					bestchi2 = chi2;
					ttgenps.Reset();
					ttgenps.THad(th);
					ttgenps.BLep(bl);
					ttgenps.LLep(lepton, lepton->PDGID());
					ttgenps.NLep(&nu);
				}
			}
		}

		if(pstbjets.size() >= 2 && pstljets.size() >= 2)
		{
			TLorentzVector wl(*lepton + nu);
			double bestchi2had = 1.E100;
			for(size_t wa = 0 ; wa < pstljets.size() ; wa++)
			{
				for(size_t wb = 0 ; wb < wa ; wb++)
				{
					TLorentzVector wh(*pstljets[wa] + *pstljets[wb]);
					for(size_t bl = 0 ; bl < pstbjets.size() ; bl++)
					{
						if(pstbjets[bl] == pstljets[wa] || pstbjets[bl] == pstljets[wb]) {continue;}
						for(size_t bh = 0 ; bh < pstbjets.size() ; bh++)
						{
							if(pstbjets[bh] == pstljets[wa] || pstbjets[bh] == pstljets[wb] || pstbjets[bh] == pstbjets[bl]){continue;}
							TLorentzVector th(wh + *pstbjets[bh]);
							TLorentzVector tl(wl + *pstbjets[bl]);
							double chi2had = pow(th.M() - Mt, 2) + pow(tl.M() - Mt, 2) + pow(wh.M() - MW, 2);
							double chi2 = pow(th.M() - Mt, 2) + pow(tl.M() - Mt, 2);
							if(bestchi2 < 0 || (chi2 < bestchi2 && chi2had < bestchi2had))
							{
								bestchi2 = chi2;
								bestchi2had = chi2had;
								ttgenps.Reset();
								ttgenps.JaHad(pstljets[wa]);
								ttgenps.JbHad(pstljets[wb]);
								ttgenps.BHad(pstbjets[bh]);
								ttgenps.BLep(pstbjets[bl]);
								ttgenps.LLep(lepton, lepton->PDGID());
								ttgenps.NLep(&nu);
							}
						}
					}
				}
			}
		}
	}

	for(GenSimpleParticle* tl : tlepboosted)
	{
		for(TLorentzVector* th : thadboosted)
		{
			if(th->DeltaR(*tl) < 1.2) {continue;}

			double chi2 = pow(th->M()-Mt, 2)+pow(tl->M() - Mt, 2);
			if(bestchi2 < 0 || chi2 < bestchi2)
			{
				bestchi2 = chi2;
				ttgenps.Reset();
				ttgenps.LLep(nullptr, (abs(tl->PDGID()) -600)*tl->PDGID()/abs(tl->PDGID()));
				ttgenps.THad(th);
				ttgenps.TLep(tl);
			}
		}
	}		

	double bestchi2had = 1.E100;
	for(GenSimpleParticle* tl : tlepboosted)
	{
		for(size_t wa = 0 ; wa < pstljets.size() ; wa++)
		{
			if(pstljets[wa]->DeltaR(*tl) < 1.2) {continue;} 
			for(size_t wb = 0 ; wb < wa ; wb++)
			{
				if(pstljets[wb]->DeltaR(*tl) < 1.2) {continue;} 
				TLorentzVector wh(*pstljets[wa] + *pstljets[wb]);
				for(size_t bh = 0 ; bh < pstbjets.size() ; bh++)
				{
					if(pstbjets[bh]->DeltaR(*tl) < 1.2){continue;}
					TLorentzVector th(wh + *pstbjets[bh]);
					double chi2had = pow(th.M() - Mt, 2) + pow(wh.M() - MW, 2);
					double chi2 = pow(th.M()-Mt, 2)+pow(tl->M() - Mt, 2);
					if(bestchi2 < 0 || (chi2 < bestchi2 && chi2had < bestchi2had))
					{
						bestchi2 = chi2;
						bestchi2had = chi2had;
						ttgenps.Reset();
						ttgenps.LLep(nullptr, (abs(tl->PDGID()) -600)*tl->PDGID()/abs(tl->PDGID()));
						ttgenps.JaHad(pstljets[wa]);
						ttgenps.JbHad(pstljets[wb]);
						ttgenps.BHad(pstbjets[bh]);
						ttgenps.TLep(tl);
					}
				}
			}
		}
	}

	if(ttgenps.IsComplete())
	{
		vector<TLorentzVector*> genaddjets;
		for(TLorentzVector* gj : pstbjets)
		{
			if(ttgenps.IsJetIn(gj) == -1)
			{
				genaddjets.push_back(gj);
			}
		}
		for(TLorentzVector* gj : pstljets)
		{
			if(ttgenps.IsJetIn(gj) == -1)
			{
				genaddjets.push_back(gj);
			}
		}
		ttgenps.SetAdditionalJets(genaddjets);
		ttgenps.Calculate();
		hgenPS.Fill(ttgenps, weight);
		if(ttgenps.IsCompleteRH() && ttgenps.IsCompleteRL())
		{
			hgenRES.Fill(ttgenps, weight);
		}
		else if(ttgenps.IsCompleteBH() && ttgenps.IsCompleteRL())
		{
			hgenBH.Fill(ttgenps, weight);
		}
		else if(ttgenps.IsCompleteRH() && ttgenps.IsCompleteBL())
		{
			hgenBL.Fill(ttgenps, weight);
		}
		else if(ttgenps.IsCompleteBH() && ttgenps.IsCompleteBL())
		{
			hgenBHBL.Fill(ttgenps, weight);
		}
		else
		{
			hgenOther.Fill(ttgenps, weight);
		}
	}
}


void BoostedTop::SelectPseudoTop()
{
    GenSimpleParticle* lepton = nullptr;
    int lc = 0;
    nu = TLorentzVector(0.,0.,0.,0.);
    vector<TLorentzVector*> pstbjets;
    vector<TLorentzVector*> pstljets;
	for(size_t n = 0 ; n < NumPLObjects() ; ++n)
    {
		GenSimpleParticle pl(GetPLObject(n));	
        if(abs(pl.PDGID()) < 6) continue;

        if(abs(pl.PDGID()) == 12 || abs(pl.PDGID()) == 14 || abs(pl.PDGID()) == 16)
        {
            nu += pl;
            continue;
        }

        if(pl.PDGID() == 22)
        {
            if(pl.IsoR3() < 0.25 && abs(pl.Eta()) < 2.4 && pl.Pt() > 15.)
            {
                SBasicParticles.push_back(pl);
                pstphotons.push_back(&(SBasicParticles.back()));
            }
        }
        if(abs(pl.PDGID()) == 13 || abs(pl.PDGID()) == 11)
        {
            if(pl.IsoR3() < 0.35)
            {
                if(abs(pl.Eta()) < 2.4 && pl.Pt() > 15)
                {
                    lc++;
                }
                if(abs(pl.Eta()) < 2.4 && pl.Pt() > 30.)
                {
                    SBasicParticles.push_back(pl);
                    lepton = &(SBasicParticles.back());
                }
            }
        }
    }
    if(lc >= 2 || lepton == nullptr) {return;}
	for(size_t n = 0 ; n < NumPLObjects() ; ++n)
    {
		GenSimpleParticle pl(GetPLObject(n));	
        if(abs(pl.PDGID()) > 6) {continue;}
        if(pl.Pt() < 25 || abs(pl.Eta()) > 2.4) {continue;}
        if(lepton->DeltaR(pl) < 0.4) {continue;}
        bool isph = false;
        for(TLorentzVector* ph : pstphotons) {if(ph->DeltaR(pl) < 0.2) {isph = true;}}
        if(isph) {continue;}
        SBasicParticles.push_back(pl);
        if(abs(pl.PDGID()) == 5)
        {
            pstbjets.push_back(&(SBasicParticles.back()));
        }
        else
        {
            pstljets.push_back(&(SBasicParticles.back()));
        }
    }
    if(pstbjets.size() < 2 || pstljets.size() < 2){return;}

    TLorentzVector wl(*lepton + nu);
    double chi2min = 1.E100;
    double MW = 80.4;
    double Mt = 172.5;
    for(size_t wa = 0 ; wa < pstljets.size() ; wa++)
    {
        for(size_t wb = 0 ; wb < wa ; wb++)
        {
            //if((pstljets[wa]->Pt() < 25. && pstljets[wa]->Pt() < 25.)){continue;}
            TLorentzVector wh(*pstljets[wa] + *pstljets[wb]);
            for(size_t bl = 0 ; bl < pstbjets.size() ; bl++)
            {
                if(pstbjets[bl] == pstljets[wa] || pstbjets[bl] == pstljets[wb]) {continue;}
                for(size_t bh = 0 ; bh < pstbjets.size() ; bh++)
                {
                    if(pstbjets[bh] == pstljets[wa] || pstbjets[bh] == pstljets[wb] || pstbjets[bh] == pstbjets[bl]){continue;}
                    //if(pstbjets[bh]->Pt() < 25. && pstbjets[bl]->Pt() < 25.) {continue;}
                    TLorentzVector th(wh + *pstbjets[bh]);
                    TLorentzVector tl(wl + *pstbjets[bl]);
                    double chi2 = Power(th.M() - Mt, 2) + Power(tl.M() - Mt, 2) + Power(wh.M() - MW, 2);
                    if(chi2 < chi2min)
                    {
						chi2min = chi2;
						ttgenps.JaHad(pstljets[wa]);
						ttgenps.JbHad(pstljets[wb]);
						ttgenps.BHad(pstbjets[bh]);
						ttgenps.BLep(pstbjets[bl]);
						ttgenps.LLep(lepton, lepton->PDGID());
						ttgenps.NLep(&nu);
					}
                }
            }
        }
    }
	if(ttgenps.IsComplete())
	{
		vector<TLorentzVector*> genaddjets;
		for(TLorentzVector* gj : pstbjets)
		{
			if(ttgenps.IsJetIn(gj) == -1)
			{
				genaddjets.push_back(gj);
			}
		}
		for(TLorentzVector* gj : pstljets)
		{
			if(ttgenps.IsJetIn(gj) == -1)
			{
				genaddjets.push_back(gj);
			}
		}
		ttgenps.SetAdditionalJets(genaddjets);
		ttgenps.Calculate();
		hgenPS.Fill(ttgenps, weight);
	}
}

void BoostedTop::SelectGENParticles()
{
	if(ran_SelectGENParticles_ == true){return;}
	ran_SelectGENParticles_ = true;
	if(NumPartonicTops() == 0) {return;}
	for(UInt_t i = 0 ; i < NumPartonicTops() ; i++)
	{
		gps.push_back(GenSimpleParticle(GetPartonicTop(i)));
		//cout << i << ": " << gps.back().PDGID() << " " << gps.back().Px() << " " <<gps.back().Py() << " " <<gps.back().Pz() << " " <<gps.back().E() << endl;
	}


	bool tlep = false;
	bool tbarlep = false;
	bool tem = false;
	bool tbarem = false;
	if(gps[2].PDGID() == -11 || gps[2].PDGID() == -13){tem = true;}
	if(gps[7].PDGID() == 11 || gps[7].PDGID() == 13){tbarem = true;}
	if(gps[2].PDGID() == -11 || gps[2].PDGID() == -13 || gps[2].PDGID() == -15){tlep = true;}
	if(gps[7].PDGID() == 11 || gps[7].PDGID() == 13 || gps[7].PDGID() == 15){tbarlep = true;}

	gent = &gps[0];
	gentbar = &gps[4];

	weight *=  1. + ctopptweight_*(gent->Pt()/750. - 1.);
	weight *=  1. + cttmweight_*(gent->Pt()/1500. - 1.);

	if(tem && !tbarlep)
	{
		ttgenpa.TLep(&gps[0]);
		ttgenpa.BLep(&gps[1]);
		ttgenpa.LLep(&gps[2], gps[2].PDGID());
		ttgenpa.NLep(&gps[3]);
		ttgenpa.THad(&gps[4]);
		ttgenpa.BHad(&gps[5]);
		ttgenpa.JaHad(&gps[6]);
		ttgenpa.JbHad(&gps[7]);
		_hists1d["genttdecay"]->Fill(3.5, weight);

		ttgenlhe.LLep(&gps[2], gps[2].PDGID());
		ttgenlhe.THad(&gps[8]);
		ttgenlhe.TLep(&gps[9]);
	}
	else if(!tlep && tbarem)
	{
		ttgenpa.TLep(&gps[4]);
		ttgenpa.BLep(&gps[5]);
		ttgenpa.LLep(&gps[7], gps[7].PDGID());
		ttgenpa.NLep(&gps[6]);
		ttgenpa.THad(&gps[0]);
		ttgenpa.BHad(&gps[1]);
		ttgenpa.JaHad(&gps[2]);
		ttgenpa.JbHad(&gps[3]);
		_hists1d["genttdecay"]->Fill(3.5, weight);

		ttgenlhe.LLep(&gps[7], gps[7].PDGID());
		ttgenlhe.THad(&gps[9]);
		ttgenlhe.TLep(&gps[8]);
	}

	if(tlep && tbarlep)
	{
		_hists1d["genttdecay"]->Fill(2.5, weight);
		DILEP_ = true;
	}
	else if(!tlep && !tbarlep)
	{
		_hists1d["genttdecay"]->Fill(0.5, weight);
		ALLHAD_ = true;
		HadTops.push_back(&gps[0]);
		HadTops.push_back(&gps[4]);
	}
	else if(!tlep)
	{
		_hists1d["genttdecay"]->Fill(1.5, weight);
		SEMILEP_ = true;
		HadTops.push_back(&gps[0]);

	}
	else if(!tbarlep)
	{
		_hists1d["genttdecay"]->Fill(1.5, weight);
		SEMILEP_ = true;
		HadTops.push_back(&gps[4]);
	}


	if(ttgenpa.IsComplete())
	{
		ttgenlhe.Calculate();
		hgenLHE.Fill(ttgenlhe, weight);
		ttgenpa.Calculate();
		hgenAll.Fill(ttgenpa, weight);
		if(abs(ttgenpa.THad()->Eta()) < 2.4)
		{
			_hists2d["ptgen_drmaxgen"]->Fill(ttgenpa.THad()->Pt(), ttgenpa.DRmaxHad());
			_hists2d["ptgen_drmingen"]->Fill(ttgenpa.THad()->Pt(), ttgenpa.DRminHad());
		}
	}
}

void BoostedTop::FillGenJets()
{
	//Gen Jets
	double bweights = 1.;
	promptph = false;
	for(size_t n = 0 ; n < NumPLObjects() ; ++n)
	{
		GenSimpleParticle pl(GetPLObject(n));	
        if(pl.PDGID() == 22 && ((int)pl.Mother() & (1<<0)) == 0)
        {
			promptph = true;
		}
		if(abs(pl.PDGID()) > 6) continue;
		if(pl.Pt() < 25 || Abs(pl.Eta()) > 5.0) continue;
		SBasicParticles.push_back(pl);
		GenAllJets.push_back(&(SBasicParticles.back()));
		if(abs(pl.PDGID()) == 5)
		{
			GenBJets.push_back(&(SBasicParticles.back()));
			double fragweight = bfragweights.Weight(pl);
			bweights *= fragweight;
            bweights *= bdecayweights.Weight(pl);
			_hists1d["jet_xb"]->Fill(pl.Xb(), fragweight);
		}
		else if(abs(pl.PDGID()) == 4)
		{
			GenCJets.push_back(&(SBasicParticles.back()));
		}
		else
		{
			GenLJets.push_back(&(SBasicParticles.back()));
		}
    }
	weight *= bweights;


}

