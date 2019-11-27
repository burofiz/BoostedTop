#include "BoostedTop.h"

void BoostedTop::GenRecoMatching()
{

	for(OJet* jet : RecAK4s)
	{
		auto genjetb = find_if(GenBJets.begin(), GenBJets.end(), [&](TLorentzVector* bp){return jet->DeltaR(*bp) < 0.2 && abs(jet->E()/bp->E()-1) < 1.0;});
		if(genjetb != GenBJets.end())
		{
			double res = (jet->Pt()-(*genjetb)->Pt())/(*genjetb)->Pt();
			if(abs(jet->Eta()) < 1.5)
			{
				_hists2d["jet_Bb_pt_res"]->Fill(jet->Pt(), res, weight);
			}
			else if(abs(jet->Eta()) < 2.4)
			{
				_hists2d["jet_Eb_pt_res"]->Fill(jet->Pt(), res, weight);
			}
		}
		else
		{
			auto genjetl = find_if(GenAllJets.begin(), GenAllJets.end(), [&](TLorentzVector* bp){return jet->DeltaR(*bp) < 0.2 && abs(jet->E()/bp->E()-1) < 1.0;});
			if(genjetl != GenAllJets.end())
			{
				double res = (jet->Pt()-(*genjetl)->Pt())/(*genjetl)->Pt();
				if(abs(jet->Eta()) < 1.5)
				{
					_hists2d["jet_Bl_pt_res"]->Fill(jet->Pt(), res, weight);
				}
				else if(abs(jet->Eta()) < 2.4)
				{
					_hists2d["jet_El_pt_res"]->Fill(jet->Pt(), res, weight);
				}
			}
		}	
	}

	double drmax = 0.2;
	int nbjets = 0;
	int nljets = 0;
	for(OJet* jet : RecAK4s)
	{
		auto gjb = find_if(GenBJets.begin(), GenBJets.end(), [&](TLorentzVector* bp){return jet->DeltaR(*bp) < drmax;});
		if(gjb != GenBJets.end())
		{
			nbjets++;
			if(DELPHES)
			{
				if(int(jet->FlavorInfo(0).DeepCSVb()) & 1<<0) {_hists1d["jet_pt_btagged0"]->Fill(jet->Pt());}
				if(int(jet->FlavorInfo(0).DeepCSVb()) & 1<<1) {_hists1d["jet_pt_btagged1"]->Fill(jet->Pt());}
				if(int(jet->FlavorInfo(0).DeepCSVb()) & 1<<2) {_hists1d["jet_pt_btagged2"]->Fill(jet->Pt());}
				if(int(jet->FlavorInfo(0).DeepCSVb()) & 1<<3) {_hists1d["jet_pt_btagged3"]->Fill(jet->Pt());}
				if(int(jet->FlavorInfo(0).DeepCSVb()) & 1<<4) {_hists1d["jet_pt_btagged4"]->Fill(jet->Pt());}
			}
			else
			{
				if(jet->BTag() > btag_loose_) {_hists1d["jet_pt_btagged0"]->Fill(jet->Pt());}
				if(jet->BTag() > btag_medium_) {_hists1d["jet_pt_btagged1"]->Fill(jet->Pt());}
				if(jet->BTag() > btag_tight_) {_hists1d["jet_pt_btagged2"]->Fill(jet->Pt());}
			}
			_hists1d["jet_pt_allbtagged"]->Fill(jet->Pt());
		}
		else
		{
			nljets++;
			if(DELPHES)
			{
				if(int(jet->FlavorInfo(0).DeepCSVb()) & 1<<0) {_hists1d["jet_pt_missbtagged0"]->Fill(jet->Pt());}
				if(int(jet->FlavorInfo(0).DeepCSVb()) & 1<<1) {_hists1d["jet_pt_missbtagged1"]->Fill(jet->Pt());}
				if(int(jet->FlavorInfo(0).DeepCSVb()) & 1<<2) {_hists1d["jet_pt_missbtagged2"]->Fill(jet->Pt());}
				if(int(jet->FlavorInfo(0).DeepCSVb()) & 1<<3) {_hists1d["jet_pt_missbtagged3"]->Fill(jet->Pt());}
				if(int(jet->FlavorInfo(0).DeepCSVb()) & 1<<4) {_hists1d["jet_pt_missbtagged4"]->Fill(jet->Pt());}
			}
			else
			{
				if(jet->BTag() > btag_loose_) {_hists1d["jet_pt_missbtagged0"]->Fill(jet->Pt());}
				if(jet->BTag() > btag_medium_) {_hists1d["jet_pt_missbtagged1"]->Fill(jet->Pt());}
				if(jet->BTag() > btag_tight_) {_hists1d["jet_pt_missbtagged2"]->Fill(jet->Pt());}
			}
			_hists1d["jet_pt_allmissbtagged"]->Fill(jet->Pt());

		}
	}	
	_hists2d["jet_nbjets"]->Fill(nbjets, nljets);

	if(!ttgen.IsCompleteRES()) return;
	tttruthres.MET(&met);


	if(Abs(ttgen.LPDGID()) == 13)
	{
		vector<OMuon*>::const_iterator mu = min_element(RecMus.begin(), RecMus.end(), [&](OMuon* a, OMuon* b){return ttgen.LLep()->DeltaR(*a) < ttgen.LLep()->DeltaR(*b);});
		if(mu != RecMus.end() && (*mu)->DeltaR(*ttgen.LLep()) < 0.2)
		{
			tttruthres.LLep(*mu, -13*(*mu)->Charge());
		}
	}
	if(Abs(ttgen.LPDGID()) == 11)
	{
		vector<OElectron*>::iterator el = min_element(RecEls.begin(), RecEls.end(), [&](OElectron* a, OElectron* b){return ttgen.LLep()->DeltaR(*a) < ttgen.LLep()->DeltaR(*b);});
		if(el != RecEls.end() && (*el)->DeltaR(*ttgen.LLep()) < 0.2)
		{
			tttruthres.LLep(*el, -11* (*el)->Charge());
		}
	}



	for(OJet* jet : RecAK4s)
	{
		if(jet->DeltaR(*ttgen.BLep()) < drmax && (tttruthres.BLep() == nullptr || jet->Pt() > tttruthres.BLep()->Pt()))
		{
			tttruthres.BLep(jet);
		}
		if(jet->DeltaR(*ttgen.BHad()) < drmax && (tttruthres.BHad() == nullptr || jet->Pt() > tttruthres.BHad()->Pt()))
		{
			tttruthres.BHad(jet);
		}
		if(jet->DeltaR(*ttgen.JaHad()) < drmax && (tttruthres.JaHad() == nullptr || jet->Pt() > tttruthres.JaHad()->Pt()))
		{
			tttruthres.JaHad(jet);
		}
		if(jet->DeltaR(*ttgen.JbHad()) < drmax && (tttruthres.JbHad() == nullptr || jet->Pt() > tttruthres.JbHad()->Pt()))
		{
			tttruthres.JbHad(jet);
		}
	}

	tttruthboosted = tttruthres;

	for(HadTopJet* htjet : RecHadJets)
	{
		if(htjet->IsSignalJet())
		{
			tttruthboosted.THad(htjet);
			break;
		}
	}

	drmax = 0.4;
	for(LepTopJet* ltjet : RecLepJets)
	{
		if(ltjet->DeltaR(*ttgen.TLep()) < drmax)
		{
			drmax = ltjet->DeltaR(*ttgen.TLep());
			tttruthboosted.TLep(ltjet);
		}
	}

	if(tttruthboosted.IsCompleteBL())
	{
		hgenLepBoosted.Fill(ttgen, weight);
	}
	else if(tttruthres.IsCompleteRL())
	{
		hgenLepRes.Fill(ttgen, weight);
	}
	else
	{
		hgenLepNo.Fill(ttgen, weight);
	}


	if(tttruthres.IsCompleteRH() && tttruthboosted.IsCompleteBH() && dynamic_cast<HadTopJet*>(tttruthboosted.THad())->NumMatchPartons() == 3)
	{
		hgenHadAll.Fill(ttgen, weight);
	}
	else if(tttruthboosted.IsCompleteBH() && dynamic_cast<HadTopJet*>(tttruthboosted.THad())->NumMatchPartons() == 3 && dynamic_cast<HadTopJet*>(tttruthboosted.THad())->Num_MemberJets() == 1)
	{
		hgenHadBoosted.Fill(ttgen, weight);
	}
	else if(tttruthboosted.IsCompleteBH() && dynamic_cast<HadTopJet*>(tttruthboosted.THad())->NumMatchPartons() == 3 && dynamic_cast<HadTopJet*>(tttruthboosted.THad())->Num_MemberJets() == 2)
	{
		hgenHadBoostedMerged.Fill(ttgen, weight);
	}
	else if(tttruthboosted.IsCompleteBH() && dynamic_cast<HadTopJet*>(tttruthboosted.THad())->NumMatchPartons() == 2 && tttruthres.IsCompleteRH())
	{
		hgenHadBrokenRes.Fill(ttgen, weight);
	}
	else if(tttruthres.IsCompleteRH())
	{
		hgenHadRes.Fill(ttgen, weight);
	}
	else if(tttruthboosted.IsCompleteBH() && dynamic_cast<HadTopJet*>(tttruthboosted.THad())->NumMatchPartons() == 2)
	{
		hgenHadBroken.Fill(ttgen, weight);
	}
	else
	{
		hgenHadNo.Fill(ttgen, weight);
	}

}
