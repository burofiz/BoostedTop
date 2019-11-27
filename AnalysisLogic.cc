#include "BoostedTop.h"
#include "BAnalysis.h"

#include <TRandom3.h>

bool BoostedTop::AnalysisNonIsoLep()
{
	BoostedTL();

	for(LepTopJet* ltj : RecLepJets)
	{
		if(ltj->IsTop() &&  abs(ltj->LepPDGID()) == 13) {lepjetmu_right.Fill(ltj, weight);}
		else if(ltj->IsTop() &&  abs(ltj->LepPDGID()) == 11) {lepjetel_right.Fill(ltj, weight);}
		else if(!ltj->IsTop() &&  abs(ltj->LepPDGID()) == 13) {lepjetmu_wrong.Fill(ltj, weight);}
		else if(!ltj->IsTop() &&  abs(ltj->LepPDGID()) == 11) {lepjetel_wrong.Fill(ltj, weight);}
	}

	sort(RecLepJets.begin(), RecLepJets.end(), [](LepTopJet* A, LepTopJet* B) {return A->NNRes() > B->NNRes();});

	int goodtl = count_if(RecLepJets.begin(), RecLepJets.end(), [](LepTopJet* ltj) {return ltj->NNRes() > 0.8;});
	if(goodtl == 1)
	{
		LepTopJet* ltj = RecLepJets[0];
		TTEvent ttrecres;
		double kmin = -0.5;
		for(OJet* bh : RecAK4s)
		{
			if(bh->BTag() < btag_medium_) {continue;}
			if(bh->DeltaR(*ltj->Jet()) < 0.4) {continue;}
			for(OJet* ja : RecAK4s)
			{
				if(ja == bh) continue;
				if(ja->DeltaR(*ltj->Jet()) < 0.4) {continue;}
				for(OJet* jb : RecAK4s)
				{
					if(jb == bh || jb == ja) continue;
					if(jb->DeltaR(*ltj->Jet()) < 0.4) {continue;}
					TLorentzVector W(*ja + *jb);
					TLorentzVector T(W+*bh);
					double K = h_tmwm->Interpolate(min(499., T.M()), min(499., W.M()));
					if(K > kmin)
					{
						kmin = K;
						ttrecres.TLep(ltj);
						ttrecres.LLep(ltj->Lep(), ltj->LepPDGID());
						ttrecres.BHad(bh);
						ttrecres.JaHad(ja);
						ttrecres.JbHad(jb);
					}
				}
			}
		}


		_hists1d["counter"]->Fill(18.5, weight);
		_hists1d["nhadjets_noniso"]->Fill(RecHadJets.size()-1, weight);
		HadTopJet* htj = nullptr;
		sort(RecHadJets.begin(), RecHadJets.end(), [](HadTopJet* A, HadTopJet* B){return A->NNRes() > B->NNRes();});
		for(HadTopJet* testhtj : RecHadJets)
		{
			if(!testhtj->CloseToMember(&(ltj->PureJet()), 0.8) && !testhtj->CloseToMember(ltj->Lep(), 0.8))
			{
				htj = testhtj;
				break;
			}
		}

		if(htj != nullptr)
		{
			//if(TRAIN) {htj->SetTree(tophaddata);}
			if(htj->IsHadTopJet())
			{
				hadjetnonisolep_right.Fill(htj, weight);
				//if(TRAIN) {tophadtreesig->Fill();}
			}
			else if(htj->IsBrokenHadTopJet())
			{
				hadjetnonisolep_wrong.Fill(htj, weight);
				//if(TRAIN) {tophadtreesig->Fill();}
			}
			else
			{
				hadjetnonisolep_other.Fill(htj, weight);
				//if(TRAIN) {tophadtreebkg->Fill();}
			}
			TTEvent ttrec;
			ttrec.THad(htj);
			ttrec.TLep(ltj);
			ttrec.LLep(ltj->Lep(), ltj->LepPDGID());
			ttrec.Calculate();
			ttrec.SetAdditionalJets(RecAK4s);
			ttrec.SetAdditionalHTJets(RecHadJets);
			tleppt = ttrec.TLep()->Pt();
			thadpt  = ttrec.THad()->Pt();
			thardpt = ttrec.THard()->Pt();
			tsoftpt = ttrec.TSoft()->Pt();
			toppt = ttrec.T()->Pt();
			topy = abs(ttrec.T()->Rapidity());
			topbarpt = ttrec.TBar()->Pt();
			topbary = abs(ttrec.TBar()->Rapidity());
			thady  = abs(ttrec.THad()->Rapidity());
			tlepy  = abs(ttrec.TLep()->Rapidity());
			ttm = ttrec.TT()->M();
			ttpt = ttrec.TT()->Pt();
			tty = abs(ttrec.TT()->Rapidity());
			thadpt_ttm_tty = unfolding_all_.GetBin2D("ttm+tty", ttrec.TT()->M(), abs(ttrec.TT()->Rapidity()));
			thadpt_ttm_cts = unfolding_all_.GetBin2D("ttm+cts", ttrec.TT()->M(), ttrec.CosThetaCMS());
			thadptref  = ttrec.THad()->Pt();
			thadetaref  = ttrec.THad()->Eta();
			cts  = ttrec.CosThetaCMS();
			nnres = htj->NNRes();
			type = htj->IsSignalJet();
			nvtx = NumGoodVertices();
			btag = htj->BTag();
			if(cfillsigtree_) {ttbartree->Fill();}
			//std::cout<<"AL.cc, 11"<<endl;
			_hists1d["counter"]->Fill(19.5, weight);
			if(htj->IsSignalJet())
			{
				if(ttgen.IsComplete())
				{
					unfolding_all_.FillTruthReco(ttgen, ttrec, "BHBL", weight);
					unfolding_all_.FillRes(ttgen, ttrec, "BHBL", weight);
				}
				unfolding_all_.FillReco(ttrec, "BHBL", weight);
			}

			if(ttgen.IsComplete() && ttrec.IsCorrectKin(ttgen))
			{
				_hists2d["BB_pthad_genrec"]->Fill(ttgen.THad()->Pt(), ttrec.THad()->Pt());
				_hists2d["BB_ptlep_genrec"]->Fill(ttgen.TLep()->Pt(), ttrec.TLep()->Pt());
				hrecBHBLright.Fill(ttrec, weight);
				//cout << kmin << ": " <<  ttrec.IsTHadCorrectKin(ttgen) << " " << htj->DeltaR(*ttgen.THad()) << " " << ttgen.THad()->Pt() << " CORRECT " << ttrec.IsTLepCorrectKin(ttgen) << endl; 
			}
			else if(ttgen.IsComplete() && tttruthboosted.IsComplete())
			{
				//cout << kmin << ": " << ttrec.IsTHadCorrectKin(ttgen) << " " << htj->DeltaR(*ttgen.THad()) << " " << ttgen.THad()->Pt() << " WRONG " << ttrec.IsTLepCorrectKin(ttgen) << endl; 
				hrecBHBLwrong.Fill(ttrec, weight);
			}
			else if(ttgen.IsComplete())
			{
				hrecBHBLnonreco.Fill(ttrec, weight);
			}
			else
			{
				hrecBHBLother.Fill(ttrec, weight);
			}
		}
		else if(kmin > 0)
		{
			std::cout<<"AL.cc, 151"<<endl;
			_hists1d["counter"]->Fill(20.5, weight);

			ttrecres.Calculate();
			ttrecres.SetAdditionalJets(RecAK4s);
			ttrecres.SetAdditionalHTJets(RecHadJets);
			unfolding_all_.FillReco(ttrecres, "BL", weight);
			if(ttgen.IsComplete() && ttrecres.IsCorrectKin(ttgen))
			{
				hrecBLright.Fill(ttrecres, weight);
				unfolding_all_.FillTruthReco(ttgen, ttrecres, "BL", weight);
				unfolding_all_.FillRes(ttgen, ttrecres, "BL", weight);
			}
			else if(ttgen.IsComplete() && tttruthboosted.IsComplete())
			{
				hrecBLwrong.Fill(ttrecres, weight);
			}
			else if(ttgen.IsComplete())
			{
				hrecBLnonreco.Fill(ttrecres, weight);
				unfolding_all_.FillTruthReco(ttgen, ttrecres, "BL", weight);
			}
			else
			{
				hrecBLother.Fill(ttrecres, weight);
			}
		}
	}
	return true;
}



bool BoostedTop::AnalysisIsoLep()
{
	sort(RecAK4s.begin(), RecAK4s.end(), [](OJet* a, OJet* b){return a->BTag() > b->BTag();});
	if(RecAK4s.size() >= 4)
	{
		_hists1d["counter"]->Fill(25.5, weight);
		_hists2d["btaghigh_deepcsv"]->Fill(RecAK4s[0]->BTag(), RecAK4s[0]->Pt(), weight);
		_hists2d["btaglow_deepcsv"]->Fill(RecAK4s[1]->BTag(), RecAK4s[1]->Pt(), weight);

		if(RecAK4s[0]->BTag() > btag_medium_) _hists1d["counter"]->Fill(26.5, weight);
		if(RecAK4s[1]->BTag() > btag_medium_) _hists1d["counter"]->Fill(27.5, weight);
	}

	TTEvent ttrec_boost;
	if(RecMus.size() == 1)
	{
		ttrec_boost.LLep(RecMus[0], -13*RecMus[0]->Charge());
	}		
	if(RecEls.size() == 1)
	{
		ttrec_boost.LLep(RecEls[0], -11*RecEls[0]->Charge());
	}		
	ttrec_boost.MET(&met);
	sort(RecHadJets.begin(), RecHadJets.end(), [](HadTopJet* A, HadTopJet* B){return A->NNRes() > B->NNRes();});
	_hists1d["nhadjets_iso"]->Fill(RecHadJets.size(), weight);
	for(HadTopJet* htj : RecHadJets)
	{
		if(htj->CloseToMember(ttrec_boost.LLep(), 0.8)) continue;
		float btagmax = btag_medium_;
		OJet* bestbl = nullptr;
		for(OJet* bl : RecAK4s)
		{
			if(!htj->CloseToMember(bl, 0.8) && bl->BTag() > btagmax)
			{
				btagmax = bl->BTag();
				bestbl = bl;
			}
		}

		if(bestbl != nullptr)
		{
			ttrec_boost.THad(htj);
			ttrec_boost.BLep(bestbl);
			ttrec_boost.Calculate();
			if(ttrec_boost.Dnu() != ttrec_boost.Dnu() || ttrec_boost.Dnu() < 0)
			{
				ttrec_boost.THad(nullptr);
				ttrec_boost.BLep(nullptr);
				//continue;
			}
			//break;
		}
		break;
	}

	if(ttrec_boost.THad() != nullptr)
	{
		_hists1d["counter"]->Fill(30.5, weight);
		ttrec_boost.SetAdditionalJets(RecAK4s);
		ttrec_boost.SetAdditionalHTJets(RecHadJets);
	}

	TTEvent ttrec_res;
	AnalyseResolved(ttrec_res);

//	if(ttrec_res.IsComplete())
//	{
//		_hists1d["counter"]->Fill(31.5, weight);
//		ttrec_res.Calculate();
//		ttrec_res.SetAdditionalJets(RecAK4s);
//		ttrec_res.SetAdditionalHTJets(RecHadJets);
//	}

	if(ttrec_boost.IsComplete() && ttrec_res.IsComplete())
	{
		double mdr = ttrec_boost.THad()->DeltaR(*ttrec_res.THad());
		_hists1d["compare_dr"]->Fill(mdr);
		if(mdr < 1.)
		{
		_hists2d["compare_pt"]->Fill(ttrec_boost.THad()->Pt(), ttrec_boost.THad()->Pt()/ttrec_res.THad()->Pt()-1);
		}
	}


	if(!BTAG_EFF_MODE_ && ttrec_boost.IsComplete() && !ttrec_res.IsComplete())
	{
		_hists1d["counter"]->Fill(32.5, weight);
		HadTopJet* htj = dynamic_cast<HadTopJet*>(ttrec_boost.THad());

		//if(TRAIN) {htj->SetTree(tophaddata);}
		if(htj->IsHadTopJet())
		{
			hadjetisolep_right.Fill(htj, weight);
			//if(TRAIN) {tophadtreesig->Fill();}
		}
		else if(htj->IsBrokenHadTopJet())
		{
			hadjetisolep_wrong.Fill(htj, weight);
			//if(TRAIN) {tophadtreesig->Fill();}
		}
		else
		{
			hadjetisolep_other.Fill(htj, weight);
			//if(TRAIN) {tophadtreebkg->Fill();}
		}

		tleppt = ttrec_boost.TLep()->Pt();
		thadpt  = ttrec_boost.THad()->Pt();
		thardpt = ttrec_boost.THard()->Pt();
		tsoftpt = ttrec_boost.TSoft()->Pt();
		toppt = ttrec_boost.T()->Pt();
		topy = abs(ttrec_boost.T()->Rapidity());
		topbarpt = ttrec_boost.TBar()->Pt();
		topbary = abs(ttrec_boost.TBar()->Rapidity());
		thady  = abs(ttrec_boost.THad()->Rapidity());
		tlepy  = abs(ttrec_boost.TLep()->Rapidity());
		ttm = ttrec_boost.TT()->M();
		ttpt = ttrec_boost.TT()->Pt();
		tty = abs(ttrec_boost.TT()->Rapidity());
		thadptref  = ttrec_boost.THad()->Pt();
		thadetaref  = ttrec_boost.THad()->Eta();
		cts  = ttrec_boost.CosThetaCMS();
		thadpt_ttm_tty = unfolding_all_.GetBin2D("ttm+tty", ttrec_boost.TT()->M(), abs(ttrec_boost.TT()->Rapidity()));
		thadpt_ttm_cts = unfolding_all_.GetBin2D("ttm+cts", ttrec_boost.TT()->M(), ttrec_boost.CosThetaCMS());
		nnres = htj->NNRes();
		type = htj->IsSignalJet();
		nvtx = NumGoodVertices();
		btag = htj->BTag();

		if(cfillsigtree_) {ttbartree->Fill();}

		//if(htj->NNCut())
		{
			//std::cout<<"AL.cc, 317"<<endl;
			_hists1d["counter"]->Fill(33.5, weight);
			if(htj->IsSignalJet())
			{
				if(ttgen.IsComplete())
				{
					unfolding_all_.FillTruthReco(ttgen, ttrec_boost, "BH", weight);
					unfolding_all_.FillRes(ttgen, ttrec_boost, "BH", weight);
				}
				unfolding_all_.FillReco(ttrec_boost, "BH", weight);
			}

			if(ttgen.IsComplete() && ttrec_boost.IsCorrectKin(ttgen))
			{
				hrecBHright.Fill(ttrec_boost, weight);
			}
			else if(ttgen.IsComplete() && tttruthboosted.IsComplete())
			{
				hrecBHwrong.Fill(ttrec_boost, weight);
			}
			else if(ttgen.IsComplete())
			{
				hrecBHnonreco.Fill(ttrec_boost, weight);
			}
			else
			{
				hrecBHother.Fill(ttrec_boost, weight);
			}

		}
	}
	else if(ttrec_res.IsComplete())
	{	
		std::cout<<"350"<<endl;
		//std::cout<<"AL.cc, 350"<<endl;
		BAn->RunEvent(ttrec_res);//, ttgen);
		_hists1d["counter"]->Fill(34.5, weight);
		bool usebtag = false;
		if(BTAG_EFF_MODE_) { usebtag = btageff.Fill(ttrec_res, ttgen, weight);}
		std::cout<<"AL.cc, 356"<<endl;
		if(!BTAG_EFF_MODE_ || usebtag)
		{
			std::cout<<"AL.cc, 358"<<endl;
			std::cout<<"AL.cc, 359"<<endl;
			unfolding_all_.FillReco(ttrec_res, "RES", weight);

			double drmax = 0.2;
			int nbjets = 0;
			int nljets = 0;
			for(TLorentzVector* jet : {ttrec_res.BLep(), ttrec_res.BHad()})
			{
				auto gjb = find_if(GenBJets.begin(), GenBJets.end(), [&](TLorentzVector* bp){return jet->DeltaR(*bp) < drmax;});
				if(gjb != GenBJets.end())
				{
					nbjets++;
					_hists1d["selected_bjets_pt"]->Fill(jet->Pt(), weight);
				}
				else
				{
					nljets++;
					_hists1d["selected_ljets_pt"]->Fill(jet->Pt(), weight);

				}
			}	
			_hists2d["selected_bljets"]->Fill(nbjets, nljets, weight);

			std::cout<<"AL.cc, 382"<<endl;
			if(ttgen.IsComplete() && ttrec_res.IsCorrectKin(ttgen))
			{
				_hists2d["RES_Dpt"]->Fill(ttrec_res.THad()->Pt(), ttgen.THad()->Pt()/ttrec_res.THad()->Pt(), weight);
				if(ttrec_res.THad()->Eta() < -1.5){_hists2d["RES_Dpt_0"]->Fill(ttrec_res.THad()->Pt(), ttgen.THad()->Pt()/ttrec_res.THad()->Pt(), weight);}
				else if(ttrec_res.THad()->Eta() < -0.7){_hists2d["RES_Dpt_1"]->Fill(ttrec_res.THad()->Pt(), ttgen.THad()->Pt()/ttrec_res.THad()->Pt(), weight);}
				else if(ttrec_res.THad()->Eta() < 0.0){_hists2d["RES_Dpt_2"]->Fill(ttrec_res.THad()->Pt(), ttgen.THad()->Pt()/ttrec_res.THad()->Pt(), weight);}
				else if(ttrec_res.THad()->Eta() < 0.7){_hists2d["RES_Dpt_3"]->Fill(ttrec_res.THad()->Pt(), ttgen.THad()->Pt()/ttrec_res.THad()->Pt(), weight);}
				else if(ttrec_res.THad()->Eta() < 1.5){_hists2d["RES_Dpt_4"]->Fill(ttrec_res.THad()->Pt(), ttgen.THad()->Pt()/ttrec_res.THad()->Pt(), weight);}
				else {_hists2d["RES_Dpt_5"]->Fill(ttrec_res.THad()->Pt(), ttgen.THad()->Pt()/ttrec_res.THad()->Pt(), weight);}

				if(ttgen.NLep() != nullptr)
				{
					_hists2d["nusolver_tleppt_respt"]->Fill(ttrec_res.TLep()->Pt(), (ttrec_res.NLep()->Pt() - ttgen.NLep()->Pt())/ttgen.NLep()->Pt(), weight);
					_hists2d["nusolver_tleppt_resmet"]->Fill(ttrec_res.TLep()->Pt(), (met.Pt() - ttgen.NLep()->Pt())/ttgen.NLep()->Pt(), weight);
					_hists2d["nusolver_tleppt_respz"]->Fill(ttrec_res.TLep()->Pt(), (ttrec_res.NLep()->Pz() - ttgen.NLep()->Pz())/ttgen.NLep()->Pz(), weight);
				}
				hrecRESright.Fill(ttrec_res, weight);
				unfolding_all_.FillTruthReco(ttgen, ttrec_res, "RES", weight);
				unfolding_all_.FillRes(ttgen, ttrec_res, "RES", weight);
			}
			else if(ttgen.IsComplete() && tttruthres.IsComplete())
			{
				hrecRESwrong.Fill(ttrec_res, weight);
				unfolding_all_.FillTruthReco(ttgen, ttrec_res, "RES", weight);
				unfolding_all_.FillRes(ttgen, ttrec_res, "RES", weight);
			}
			else if(ttgen.IsComplete())
			{
				hrecRESnonreco.Fill(ttrec_res, weight);
				unfolding_all_.FillTruthReco(ttgen, ttrec_res, "RES", weight);
			}
			else
			{
				hrecRESother.Fill(ttrec_res, weight);
			}
		}
	}

	return true;
}


bool BoostedTop::AnalysisBKG(int prescale)
{

	if(!(RecAK4s.size() > 0 && RecHadJets.size() > 0 && RecMusVeto.size() + RecElsVeto.size() == 0 && (RecLepJets.size() == 0 || RecLepJets[0]->NNRes() < 0.8) && (!QCDMC_ || !promptph))){ return false;}

	double trgprob = 1.;
	if(trgres["Jet500"].GetResult() > 0) {int trgpre = trgres["Jet500"].GetResult(); trgprob *= 1./trgpre;}
	else if(trgres["Jet450"].GetResult() > 0) {int trgpre = trgres["Jet450"].GetResult(); trgprob *= 1./trgpre;}
	else if(trgres["Jet400"].GetResult() > 0) {int trgpre = trgres["Jet400"].GetResult(); trgprob *= 1./trgpre;}
	else if(trgres["Jet320"].GetResult() > 0) {int trgpre = trgres["Jet320"].GetResult(); trgprob *= 1./trgpre;}
	else if(trgres["Jet260"].GetResult() > 0) {int trgpre = trgres["Jet260"].GetResult(); trgprob *= 1./trgpre;}
	else {return false;}
	if(!IsMC()) {weight/=trgprob;}	

	vector<OPhoton*> cleanphLs;
	//vector<OPhoton*> cleanphMs;
	//vector<OPhoton*> cleanphTs;
	for(OPhoton* ph : RecPhs)
	{
		//if(ph->Pt() < 45.) {continue;}
		double mindrph =  (*min_element(RecAK4s.begin(), RecAK4s.end(), [&ph](OJet* A, OJet* B){return ph->DeltaR(*A) < ph->DeltaR(*B);}))->DeltaR(*ph);
		//bool closejet = find_if(RecAK4s.begin(), RecAK4s.end(), [&ph](OJet* jet){return ph->DeltaR(*jet) < 0.8;}) != RecAK4s.end();
		//if(closejet) {continue;}
		if(ph->Pt() > 35 && mindrph > 0.5) {cleanphLs.push_back(ph);}
		//if(ph->Pt() > 45 && mindrph > 0.7) {cleanphMs.push_back(ph);}
		//if(ph->Pt() > 45 && mindrph > 0.7) {cleanphTs.push_back(ph);}
	}

	bool L = true;
	bool M = true;
	//bool T = true;

	sort(RecHadJets.begin(), RecHadJets.end(), [](HadTopJet* A, HadTopJet* B){return A->NNRes() > B->NNRes();});
	for(HadTopJet* htj : RecHadJets)
	{
		bool othertop = find_if(RecHadJets.begin(), RecHadJets.end(), [&htj](HadTopJet* otherhtj){return htj != otherhtj && otherhtj->NNRes() > 0.4;}) != RecHadJets.end();
		if(othertop) {continue;}
		//bool hashighptjet = find_if(RecAK4s.begin(), RecAK4s.end(), [&](OJet* jet){return jet->Pt() > 400.;}) != RecAK4s.end();
		//if(!hashighptjet) {continue;}

		//int naddjets = count_if(RecAK4s.begin(), RecAK4s.end(), [&](OJet* jet){return !htj->CloseToMember(jet, 1.2);});

		bool hasb = find_if(RecAK4s.begin(), RecAK4s.end(), [&](OJet* jet){return jet->BTag() > btag_medium_ && !htj->CloseToMember(jet);}) != RecAK4s.end();
		bool hasphL = find_if(cleanphLs.begin(), cleanphLs.end(), [&](OPhoton* ph){return !htj->CloseToMember(ph);}) != cleanphLs.end();
		//bool hasphM = find_if(cleanphMs.begin(), cleanphMs.end(), [&](OPhoton* ph){return !htj->CloseToMember(ph);}) != cleanphMs.end();
		//bool hasphT = find_if(cleanphTs.begin(), cleanphTs.end(), [&](OPhoton* ph){return !htj->CloseToMember(ph);}) != cleanphTs.end();

		if(L && hasphL)
		{
			hadjetbkgPH.Fill(htj, weight);
			thadptref = htj->Pt();
			thadetaref = htj->Eta();
			nnres = htj->NNRes();
			nvtx = NumGoodVertices();
			if(cfillbkgtree_) {wbkgtree->Fill();}
			//if(UserTime() == 2016) {L = false;}
			//cout << "Q" << endl; 
		}
		if(M && hasb)
		{
			hadjetbkgB.Fill(htj, weight);
			thadptref = htj->Pt();
			thadetaref = htj->Eta();
			nnres = htj->NNRes();
			nvtx = NumGoodVertices();
			if(cfillbkgtree_) {qcdbkgtree->Fill();}
			//M = false;
			//cout << "Q" << endl; 
		}
		//if(T && hasphT)
		//{
		//	hadjetbkgPHT.Fill(htj, weight);
		//	T = false;
		//	//cout << "Q" << endl; 
		//}
	}

	return true;
}
