#include "BoostedTop.h"

#include <Eigen/Dense>

using Eigen::MatrixXd;

void BoostedTop::AnalyseResolved(TTEvent& ttrec)
{
	ttrec.Prob(1E-50, 1E-50);
	if(RecAK4s.size() < 4) {return;}
	sort(RecAK4s.begin(), RecAK4s.end(), [](OJet* a, OJet* b){return a->BTag() > b->BTag();});
	//if(!BTAG_EFF_MODE_ && RecAK4s[1]->BTag() < btag_medium_) {return;}	
	if(!BTAG_EFF_MODE_ && RecAK4s[1]->BTag() < cbsidebandmin_) {return;}	
	if(!BTAG_EFF_MODE_ && RecAK4s[0]->BTag() > cbsidebandmax_) {return;}	
	//if(BTAG_EFF_MODE_ && RecAK4s.size() != 4) {return;}


	TTEvent tttest;
	if(RecMus.size() == 1)
	{
		tttest.LLep(RecMus[0], -13*RecMus[0]->Charge());
	}		
	else if(RecEls.size() == 1)
	{
		tttest.LLep(RecEls[0], -11*RecEls[0]->Charge());
	}		
	if(tttest.LLep() == nullptr) {return;}

	tttest.MET(&met);
	for(size_t nbl = 0 ; nbl < (BTAG_EFF_MODE_ ? RecAK4s.size() : 2) ; ++nbl)
	{
		OJet* bl = RecAK4s[nbl];
		tttest.BLep(bl);
		tttest.CalculateTLep();
		if(tttest.Dnu() != tttest.Dnu() || tttest.Dnu() < 0 || tttest.Dnu() > 150.) continue;
		double dnu = tttest.Dnu();
		if(ttgen.IsComplete() && tttest.IsTLepCorrectKin(ttgen))
		{
			_hists1d["Dnu_right"]->Fill(dnu, weight);
		}
		else if(ttgen.IsComplete())
		{
			_hists1d["Dnu_wrong"]->Fill(dnu, weight);
		}
		double probnu = h_dnu->Interpolate(dnu);
		for(size_t nbh = 0 ; nbh < (BTAG_EFF_MODE_ ? RecAK4s.size() : 2) ; ++nbh)
		{
			OJet* bh = RecAK4s[nbh];
			if(bl == bh) continue;
			for(size_t nja = (BTAG_EFF_MODE_ ? 0 : 2) ; nja < RecAK4s.size() ; ++nja)
			{
				OJet* ja = RecAK4s[nja];
				if(ja == bl || ja == bh) continue;
				for(size_t njb = (BTAG_EFF_MODE_ ? 0 : 2) ; njb < nja ; ++njb)
				{
					OJet* jb = RecAK4s[njb];
					if(jb == bl || jb == bh || jb == ja) continue;
					tttest.BHad(bh);
					tttest.JaHad(ja);
					tttest.JbHad(jb);
					tttest.CalculateTHad();
					if(ttgen.IsCompleteRES() && tttest.IsTHadCorrectKin(ttgen))
					{
						_hists2d["Mt_MW_right"]->Fill(tttest.THad()->M(), (*tttest.JaHad() + *tttest.JbHad()).M(), weight);
					}
					else if(ttgen.IsCompleteRES())
					{
						_hists2d["Mt_MW_wrong"]->Fill(tttest.THad()->M(), (*tttest.JaHad() + *tttest.JbHad()).M(), weight);
					}
					double probmass = h_tmwm->Interpolate(min(499., tttest.THad()->M()), min(499., (*tttest.JaHad() + *tttest.JbHad()).M()));
					//if(probmass <= 0.) continue;
					if(probmass <= 1.67E-5) continue;
					//tttest.Prob(-1.*log(probmass*probnu));
					tttest.Prob(probmass, probnu);
					if(tttest.Prob() < ttrec.Prob())
					{
						ttrec = tttest;
					}
				}
			}
		}
	}

	if(ttrec.IsComplete())
	{
		ttrec.Calculate();
		ttrec.SetAdditionalJets(RecAK4s);
		ttrec.SetAdditionalHTJets(RecHadJets);
	}

}

void BoostedTop::GenBoostedTH()
{
	if(ran_GenBoostedTH_ == true){return;}
	ran_GenBoostedTH_ = true;
	//Fill ALL! AK8 Jets
	for(UInt_t i = 0 ; i < NumIOGenAK8Jets() ; i++)
	{   
		GenJet jet(GetIOGenAK8Jet(i));
		SGenAK8s.push_back(jet);
		GenAK8s.push_back(&SGenAK8s.back());
	}

	if(gps.size() > 0)
	{
		for(size_t j = 1 ; j < 4 ; ++j)
		{
			if(abs(gps[j].PDGID()) > 5) continue;
			double drmin = 0.8;
			int posmin = -1;
			for(size_t i = 0 ; i < GenAK8s.size() ; i++)
			{
				double dr = gps[j].DeltaR(*GenAK8s[i]);
				if(dr < drmin)
				{
					posmin = i;
					drmin = dr;
				}
			}
			if(posmin != -1)
			{
				GenAK8s[posmin]->matchedtpartons+=1;
				GenAK8s[posmin]->ptpartons += gps[j];
			}
		}
		for(size_t j = 5 ; j < 8 ; ++j)
		{
			if(abs(gps[j].PDGID()) > 5) continue;

			double drmin = 0.8;
			int posmin = -1;
			for(size_t i = 0 ; i < GenAK8s.size() ; i++)
			{
				double dr = gps[j].DeltaR(*GenAK8s[i]);
				if(dr < drmin)
				{
					posmin = i;
					drmin = dr;
				}
			}
			if(posmin != -1)
			{
				GenAK8s[posmin]->matchedtbarpartons+=1;
				GenAK8s[posmin]->ptbarpartons += gps[j];
			}
		}
	}

	for(size_t i = 0 ; i < NumIOGenAK8TopJets() ; i++)
	{
		HadTopJet htjet(GetIOGenAK8TopJet(i), GenAK8s);
		if(htjet.Pt() < 350. || abs(htjet.Eta()) >  2.4){continue;}
		if(htjet.Num_SubJets() < 3) { continue;}
		//if(htjet.M() < 60.) {continue;}
		//if(htjet.Num_MemberJets() != 1) { continue;}
		bool clean = true;
		for(size_t mj = 0 ; mj < htjet.Num_MemberJets() ; ++mj)
		{
			size_t mjet = htjet.MemberJets(mj);
			if((ttgen.LLep() != nullptr && GenAK8s[mjet]->DeltaR(*ttgen.LLep()) < 0.8) || abs(GenAK8s[mjet]->Eta()) > 2.4){clean = false; break;}
		}
		if(!clean) {continue;}

		phpt = htjet.Pt(); pheta = htjet.Eta(); 
		htjet.SetTree(tophaddata);
		if(htjet.IsCleanHadTopJet())
		{
			genhadjet_right.Fill(&htjet, weight);	
			htjet.NNRes(1);
			if(ttgen.IsCompleteHad())
			{
				_hists2d["genBH_Ddr"]->Fill(htjet.Pt(),  htjet.DeltaR(*(ttgen.THad())), weight);
				_hists2d["genBH_Dpt"]->Fill(htjet.Pt(), (htjet.Pt()-ttgen.THad()->Pt())/ttgen.THad()->Pt(), weight);
			}
		}
		else if(htjet.IsCleanBrokenHadTopJet())
		{
			genhadjet_wrong.Fill(&htjet, weight);	
			htjet.NNRes(0.5);
			if(ttgen.IsCompleteHad())
			{
				_hists2d["genBHbroken_Ddr"]->Fill(htjet.Pt(),  htjet.DeltaR(*(ttgen.THad())), weight);
				_hists2d["genBHbroken_Dpt"]->Fill(htjet.Pt(), (htjet.Pt()-ttgen.THad()->Pt())/ttgen.THad()->Pt(), weight);
			}
		}
		else
		{
			htjet.NNRes(0.);
			genhadjet_other.Fill(&htjet, weight);	
		}

		SGenHadJets.push_back(move(htjet));

	}

	SGenHadJets.sort([](const HadTopJet& A, const HadTopJet& B) {return A.NNRes() > B.NNRes();});
	for(auto a = SGenHadJets.begin() ; a != SGenHadJets.end(); ++a)
	{
		for(auto b = next(a) ; b != SGenHadJets.end(); )
		{
			if(a->Overlap(*b)) {b = SGenHadJets.erase(b);}
			else {++b;}
		}
	}

	for(HadTopJet& htjet : SGenHadJets)
	{
		GenHadJets.push_back(&htjet);
	}
}


void BoostedTop::BoostedTH()
{
	if(ran_BoostedTH_ == true){return;}
	ran_BoostedTH_ = true;
	//Fill ALL! AK8 Jets
	for(UInt_t i = 0 ; i < NumIOPFAK8Jets() ; i++)
	{   
		OJet jet(GetIOPFAK8Jet(i));
		SAK8s.push_back(jet);
		RecAK8s.push_back(&SAK8s.back());
	}

	if(gps.size() > 0)
	{
		for(size_t j = 1 ; j < 4 ; ++j)
		{
			if(abs(gps[j].PDGID()) > 5) continue;

			double drmin = 0.8;
			int posmin = -1;
			for(size_t i = 0 ; i < RecAK8s.size() ; ++i)
			{
				double dr = gps[j].DeltaR(*RecAK8s[i]);
				if(dr < drmin)
				{
					posmin = i;
					drmin = dr;
				}
			}
			if(posmin != -1)
			{
				RecAK8s[posmin]->matchedtpartons+=1;
				RecAK8s[posmin]->ptpartons += gps[j];
			}
		}
		for(size_t j = 5 ; j < 8 ; ++j)
		{
			if(abs(gps[j].PDGID()) > 5) continue;

			double drmin = 0.8;
			int posmin = -1;
			for(size_t i = 0 ; i < RecAK8s.size() ; i++)
			{
				double dr = gps[j].DeltaR(*RecAK8s[i]);
				if(dr < drmin)
				{
					posmin = i;
					drmin = dr;
				}
			}
			if(posmin != -1)
			{
				RecAK8s[posmin]->matchedtbarpartons+=1;
				RecAK8s[posmin]->ptbarpartons += gps[j];
			}
		}


	}

	for(size_t i = 0 ; i < NumIOPFAK8TopJets() ; i++)
	{
		HadTopJet htjet(GetIOPFAK8TopJet(i), RecAK8s);

		if(PS_MODE_ && ttgen.IsCompleteHad() && ttgen.THad()->DeltaR(htjet) < 0.8)
		{
			htjet.IsPS(true);
		}

		HadTopJet* genhtjet = nullptr;
		for(HadTopJet* genhtjettmp : GenHadJets)
		{
			if(genhtjettmp->DeltaR(htjet) < 0.4)
			{
				genhtjet = genhtjettmp;
			}
		}

		TLorentzVector jmom(0.,0.,0.,0.);
		for(size_t mj = 0 ; mj < htjet.Num_MemberJets() ; ++mj)
		{
			OJet* mjet = RecAK8s[htjet.MemberJets(mj)];
			if(IsMC())
			{
				double sf = jetscaler.GetRes(mjet, GenAK8s, Rho(), cuncjer_);
				sf *= jetscalerPUPPI.GetScale(mjet, GenBJets, cuncjes_);
				jmom += (*mjet)*sf;
			}
			else
			{
				jmom += *mjet;
			}
			//cout << mj << ": " << mjet->EnergyCorrection() << endl;
		}
		double sc = jmom.Pt()/htjet.Pt();
		htjet.ScaleRes(sc, 1., genhtjet);

		if(htjet.Pt() < 350. || abs(htjet.Eta()) >  2.4){continue;}
		if(htjet.Num_SubJets() < 3) {continue;}
		if(PS_MODE_ && htjet.M() < 120.) {continue;}
		if(htjet.M() < 60.) {continue;}
		//if(htjet.Num_MemberJets() != 1) { continue;}
		bool clean = true;
		//cout << "Jet: " << htjet.Num_MemberJets() << " " << htjet.Pt() << " " << htjet.Eta() << endl;
		for(size_t mj = 0 ; mj < htjet.Num_MemberJets() ; ++mj)
		{
			OJet* ak8jet = RecAK8s[htjet.MemberJets(mj)];
			//if(!ak8jet->Clean(RecMusVeto, RecElsVeto, RecPhs, 0.8) || abs(ak8jet->Eta()) > 2.4){clean = false; break;}
			if(!ak8jet->Clean(RecMusVeto, RecElsVeto, {}, 0.8) || abs(ak8jet->Eta()) > 2.4){clean = false; break;}
			htjet.BTag(max(htjet.BTag(), ak8jet->BTag()));
		}
		if(!clean) {continue;}

		phpt = htjet.Pt(); pheta = htjet.Eta(); 
		htjet.SetTree(tophaddata);
		Layer layin(tophaddatasize);
		for(int n = 0 ; n < tophaddatasize ; ++n)
		{
			layin.output()(n, 0) = tophaddata[n];
		}

		if(htjet.Pt() < 500)
		{
			double nnoutcut = nethad1->apply(layin)(0,0);
			htjet.NNRes((nnoutcut+1)*0.5);
		}
		else if(htjet.Pt() < 700)
		{
			double nnoutcut = nethad2->apply(layin)(0,0);
			htjet.NNRes((nnoutcut+1)*0.5);
		}
		else if(htjet.Pt() < 1000)
		{
			double nnoutcut = nethad3->apply(layin)(0,0);
			htjet.NNRes((nnoutcut+1)*0.5);
		}
		else
		{
			double nnoutcut = nethad4->apply(layin)(0,0);
			htjet.NNRes((nnoutcut+1)*0.5);
		}

		if(htjet.IsCleanHadTopJet())
		{
			if(TRAIN) {tophadtreesig->Fill();}
			hadjet_right.Fill(&htjet, weight);	
			if(ttgen.IsComplete()) 
			{
				_hists2d["BH_Ddr"]->Fill(htjet.Pt(),  htjet.DeltaR(*(ttgen.THad())), weight);
				_hists2d["BH_Dpt"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);
				if(htjet.Eta() < -1.5){_hists2d["BH_Dpt_0"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
				else if(htjet.Eta() < -0.7){_hists2d["BH_Dpt_1"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
				else if(htjet.Eta() < 0.0){_hists2d["BH_Dpt_2"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
				else if(htjet.Eta() < 0.7){_hists2d["BH_Dpt_3"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
				else if(htjet.Eta() < 1.5){_hists2d["BH_Dpt_4"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
				else {_hists2d["BH_Dpt_5"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
			}
		}
		else if(htjet.IsCleanBrokenHadTopJet())
		{
			if(TRAIN) {tophadtreesig->Fill();}
			hadjet_wrong.Fill(&htjet, weight);	
			if(ttgen.IsComplete())
			{
				_hists2d["BHbroken_Ddr"]->Fill(htjet.Pt(),  htjet.DeltaR(*(ttgen.THad())), weight);
				_hists2d["BHbroken_Dpt"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);
				if(htjet.Eta() < -1.5){_hists2d["BHbroken_Dpt_0"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
				else if(htjet.Eta() < -0.7){_hists2d["BHbroken_Dpt_1"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
				else if(htjet.Eta() < 0.0){_hists2d["BHbroken_Dpt_2"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
				else if(htjet.Eta() < 0.7){_hists2d["BHbroken_Dpt_3"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
				else if(htjet.Eta() < 1.5){_hists2d["BHbroken_Dpt_4"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
				else {_hists2d["BHbroken_Dpt_5"]->Fill(htjet.Pt(), ttgen.THad()->Pt()/htjet.Pt(), weight);}
			}
		}
		else
		{
			if(TRAIN) {tophadtreebkg->Fill();}
			hadjet_other.Fill(&htjet, weight);	
		}

		SHadJets.push_back(move(htjet));
	}

	SHadJets.sort([](const HadTopJet& A, const HadTopJet& B) {return A.NNRes() > B.NNRes();});
	for(auto a = SHadJets.begin() ; a != SHadJets.end(); ++a)
	{
		for(auto b = next(a) ; b != SHadJets.end(); )
		{
			if(a->Overlap(*b)) {b = SHadJets.erase(b);}
			else {++b;}
		}
	}

	for(HadTopJet& htjet : SHadJets)
	{
		RecHadJets.push_back(&htjet);

		if(htjet.MatchedGenJet() == nullptr) {continue;}

		HadTopJet* genhtjet = htjet.MatchedGenJet();
		double ptres = (htjet.Pt()-genhtjet->Pt())/genhtjet->Pt();
		double mres = (htjet.M()-genhtjet->M())/genhtjet->M();
		double mtopres = (htjet.T().M()-genhtjet->T().M())/genhtjet->T().M();
		if(htjet.IsCleanHadTopJet())
		{
			if(abs(htjet.Eta()) < 1.5) {_hists2d["pt_topjetpt_B"]->Fill(htjet.Pt(), ptres);}
			else {_hists2d["pt_topjetpt_E"]->Fill(htjet.Pt(), ptres);}
			if(genhtjet->Pt() < 500)
			{
				_hists2d["nvtx_pta_topjetpt"]->Fill(GetGenInfo(0).NumTrueInteractions(), ptres); 
				_hists2d["nvtx_pta_topjetm"]->Fill(GetGenInfo(0).NumTrueInteractions(), mres); 
				_hists2d["nvtx_pta_topjetmtop"]->Fill(GetGenInfo(0).NumTrueInteractions(), mtopres); 
			}
			else if(genhtjet->Pt() < 800)
			{
				_hists2d["nvtx_ptb_topjetpt"]->Fill(GetGenInfo(0).NumTrueInteractions(), ptres); 
				_hists2d["nvtx_ptb_topjetm"]->Fill(GetGenInfo(0).NumTrueInteractions(), mres); 
				_hists2d["nvtx_ptb_topjetmtop"]->Fill(GetGenInfo(0).NumTrueInteractions(), mtopres); 
			}
			else
			{
				_hists2d["nvtx_ptc_topjetpt"]->Fill(GetGenInfo(0).NumTrueInteractions(), ptres); 
				_hists2d["nvtx_ptc_topjetm"]->Fill(GetGenInfo(0).NumTrueInteractions(), mres); 
				_hists2d["nvtx_ptc_topjetmtop"]->Fill(GetGenInfo(0).NumTrueInteractions(), mtopres); 
			}
		}
		else if(htjet.IsCleanBrokenHadTopJet())
		{
			if(abs(htjet.Eta()) < 1.5) {_hists2d["pt_brokentopjetpt_B"]->Fill(htjet.Pt(), ptres);}
			else {_hists2d["pt_brokentopjetpt_E"]->Fill(htjet.Pt(), ptres);}
			if(genhtjet->Pt() < 500)
			{
				_hists2d["nvtx_pta_brokentopjetpt"]->Fill(GetGenInfo(0).NumTrueInteractions(), ptres); 
				_hists2d["nvtx_pta_brokentopjetm"]->Fill(GetGenInfo(0).NumTrueInteractions(), mres); 
				_hists2d["nvtx_pta_brokentopjetmtop"]->Fill(GetGenInfo(0).NumTrueInteractions(), mtopres); 
			}
			else if(genhtjet->Pt() < 800)
			{
				_hists2d["nvtx_ptb_brokentopjetpt"]->Fill(GetGenInfo(0).NumTrueInteractions(), ptres); 
				_hists2d["nvtx_ptb_brokentopjetm"]->Fill(GetGenInfo(0).NumTrueInteractions(), mres); 
				_hists2d["nvtx_ptb_brokentopjetmtop"]->Fill(GetGenInfo(0).NumTrueInteractions(), mtopres); 
			}
			else
			{
				_hists2d["nvtx_ptc_brokentopjetpt"]->Fill(GetGenInfo(0).NumTrueInteractions(), ptres); 
				_hists2d["nvtx_ptc_brokentopjetm"]->Fill(GetGenInfo(0).NumTrueInteractions(), mres); 
				_hists2d["nvtx_ptc_brokentopjetmtop"]->Fill(GetGenInfo(0).NumTrueInteractions(), mtopres); 
			}
		}
		else if(htjet.NumMatchPartons() == 0)
		{
			if(abs(htjet.Eta()) < 1.5) {_hists2d["pt_jetpt_B"]->Fill(htjet.Pt(), ptres);}
			else {_hists2d["pt_jetpt_E"]->Fill(htjet.Pt(), ptres);}
			if(genhtjet->Pt() < 500)
			{
				_hists2d["nvtx_pta_jetpt"]->Fill(GetGenInfo(0).NumTrueInteractions(), ptres); 
				_hists2d["nvtx_pta_jetm"]->Fill(GetGenInfo(0).NumTrueInteractions(), mres); 
				_hists2d["nvtx_pta_jetmtop"]->Fill(GetGenInfo(0).NumTrueInteractions(), mtopres); 
			}
			else if(genhtjet->Pt() < 800)
			{
				_hists2d["nvtx_ptb_jetpt"]->Fill(GetGenInfo(0).NumTrueInteractions(), ptres); 
				_hists2d["nvtx_ptb_jetm"]->Fill(GetGenInfo(0).NumTrueInteractions(), mres); 
				_hists2d["nvtx_ptb_jetmtop"]->Fill(GetGenInfo(0).NumTrueInteractions(), mtopres); 
			}
			else
			{
				_hists2d["nvtx_ptc_jetpt"]->Fill(GetGenInfo(0).NumTrueInteractions(), ptres); 
				_hists2d["nvtx_ptc_jetm"]->Fill(GetGenInfo(0).NumTrueInteractions(), mres); 
				_hists2d["nvtx_ptc_jetmtop"]->Fill(GetGenInfo(0).NumTrueInteractions(), mtopres); 
			}
		}
	}
	sort(RecHadJets.begin(), RecHadJets.end(), [](HadTopJet* A, HadTopJet* B){return A->NNRes() > B->NNRes();});
}

bool BoostedTop::BoostedTL()
{
	if(ran_BoostedTL_ == true){return false;}
	ran_BoostedTL_ = true;

	for(OJet* jet : RecAK4s)
	{
		if(jet->BTag() < btag_loose_) {continue;}

		TLorentzVector* lep = nullptr;
		int pdgid = 0;
		bool injet = false;
		for(OMuon* mu : RecMusLoose)
		{
			injet = false;
			if(jet->DeltaR(*mu) < 0.8 && (lep == nullptr || mu->Pt() > lep->Pt()))
			{
				lep = mu;
				pdgid = -13*mu->Charge();
				for(size_t nmu = 0 ; nmu < jet->Num_MatchedMuons() ; ++nmu)
				{
					if(jet->MatchedMuons(nmu) == mu->Num()) {injet = true;}
				}
			}
		}
		for(OElectron* el : RecElsLoose)
		{
			injet = false;
			if(jet->DeltaR(*el) < 0.8 && (lep == nullptr || el->Pt() > lep->Pt()))
			{
				lep =el;
				pdgid = -11*el->Charge();
				for(size_t nel = 0 ; nel < jet->Num_MatchedElectrons() ; ++nel)
				{
					if(jet->MatchedElectrons(nel) == el->Num()) {injet = true;}
				}
			}
		}

		if(lep == nullptr) {continue;}

		SLepJets.push_back(LepTopJet(jet, &met, lep, pdgid, injet));
		LepTopJet* bljet = &(SLepJets.back());
		if(bljet->Pt() < 400.) {continue;}

		if(PS_MODE_)
		{
			if(ttgen.IsCompleteLep() && ttgen.TLep()->DeltaR(*bljet) < 0.4) {bljet->IsTop(true);}	
		}
		else if(gps.size() > 0)
		{
			int matchlep = 0;
			int matchb = 0;
			for(size_t j = 1 ; j < 4 ; ++j)
			{
				if(abs(gps[j].PDGID()) == 5 && jet->DeltaR(gps[j]) < 0.4) {matchb++;}
				if(abs(gps[j].PDGID()) == 13 && lep->DeltaR(gps[j]) < 0.2) {matchlep++;}
				if(abs(gps[j].PDGID()) == 11 && lep->DeltaR(gps[j]) < 0.2) {matchlep++;}
			}
			if(matchlep == 1 && matchb == 1) {bljet->IsTop(true);}
			matchlep = 0;
			matchb = 0;
			for(size_t j = 5 ; j < 8 ; ++j)
			{
				if(abs(gps[j].PDGID()) == 5 && jet->DeltaR(gps[j]) < 0.4) {matchb++;}
				if(abs(gps[j].PDGID()) == 13 && lep->DeltaR(gps[j]) < 0.2) {matchlep++;}
				if(abs(gps[j].PDGID()) == 11 && lep->DeltaR(gps[j]) < 0.2) {matchlep++;}
			}
			if(matchlep == 1 && matchb == 1) {bljet->IsTop(true);}
		}

		//if(ttgen.IsComplete() && jet->DeltaR(*ttgen.BLep()) < 0.3 && lep->DeltaR(*ttgen.LLep()) < 0.2)
		if(bljet->IsTop())
		{
			if(abs(pdgid) == 13)
			{
				_hists2d["BLmu_pt_Dnu"]->Fill(bljet->Pt(), bljet->Dnu());
			}
			if(abs(pdgid) == 11)
			{
				_hists2d["BLel_pt_Dnu"]->Fill(bljet->Pt(), bljet->Dnu());
			}
		}

		if(bljet->Dnu() < 0 || bljet->Dnu() != bljet->Dnu()) {continue;}

		phpt = bljet->Pt(); pheta = bljet->Eta(); 
		if(abs(pdgid) == 13)
		{
			//OMuon* mu = dynamic_cast<OMuon*>(lep);
			bljet->SetTree<OMuon>(topdata);

			Layer layin(topdatasize);
			for(int n = 0 ; n < topdatasize ; ++n)
			{
				layin.output()(n, 0) = topdata[n];
			}
			if(bljet->Pt() < 650) {bljet->NNRes(0.5*(netmu1->apply(layin)(0,0)+1));}
			else {bljet->NNRes(0.5*(netmu2->apply(layin)(0,0)+1));}
			//if(mu->IsoNearAll()/mu->IsoCentralAll() > 200.) {continue;}

		}
		if(abs(pdgid) == 11)
		{
			//OElectron* el = dynamic_cast<OElectron*>(lep);
			bljet->SetTree<OElectron>(topdata);

			Layer layin(topdatasize);
			for(int n = 0 ; n < topdatasize ; ++n)
			{
				layin.output()(n, 0) = topdata[n];
			}
			if(bljet->Pt() < 650) {bljet->NNRes(0.5*(netmu1->apply(layin)(0,0)+1));}
			else {bljet->NNRes(0.5*(netmu2->apply(layin)(0,0)+1));}
			//if(el->IsoNearAll()/el->IsoCentralAll() > 200.) {continue;}

		}
		
		if(bljet->Dnu() < 250.)	{RecLepJets.push_back(&(SLepJets.back()));}

		if(bljet->IsTop() && ttgen.IsComplete())
		{
			if(abs(pdgid) == 13)
			{
				_hists2d["BLmu_pt_Dpt"]->Fill(bljet->Pt(), (bljet->Pt()-ttgen.TLep()->Pt())/ttgen.TLep()->Pt());
				_hists2d["BLmu_pt_Ddr"]->Fill(bljet->Pt(), bljet->DeltaR(*(ttgen.TLep())));
				_hists2d["BLmu_pt_Dpz"]->Fill(bljet->Pt(), (bljet->Pz()-ttgen.TLep()->Pz())/ttgen.TLep()->Pz());
				if(TRAIN) {topmutreesig->Fill();}
			}
			if(abs(pdgid) == 11)
			{
				_hists2d["BLel_pt_Dpt"]->Fill(bljet->Pt(), (bljet->Pt()-ttgen.TLep()->Pt())/ttgen.TLep()->Pt());
				_hists2d["BLel_pt_Ddr"]->Fill(bljet->Pt(), bljet->DeltaR(*(ttgen.TLep())));
				_hists2d["BLel_pt_Dpz"]->Fill(bljet->Pt(), (bljet->Pz()-ttgen.TLep()->Pz())/ttgen.TLep()->Pz());
				if(TRAIN) {topeltreesig->Fill();}
			}
		}
		else
		{
			if(abs(pdgid) == 13)
			{
				if(TRAIN) {topmutreebkg->Fill();}
			}
			if(abs(pdgid) == 11)
			{
				if(TRAIN) {topeltreebkg->Fill();}
			}
		}
	}

	sort(RecLepJets.begin(), RecLepJets.end(), [](LepTopJet* jeta, LepTopJet* jetb) {return jeta->NNRes() > jetb->NNRes();});
	return true;
}

