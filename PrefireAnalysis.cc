#include "BoostedTop.h"

void BoostedTop::PrefireAnalysis()
{

	bool barrelel = false;
	for(OElectron* el : RecEls) {if(abs(el->Eta()) < 1.8 ) {barrelel = true; break;}}

	if(IsData() && (UserTime() == 2016 || UserTime() == 2017) && (RecMus.size() > 0 || barrelel))
	{
		if(IsUnPrefireable())
		{
			_hists1d["counter"]->Fill(98.5, weight);
			for(OPhoton* ph : RecPrefirePhs)
			{
				double preprob = 1.;
				for(size_t nl1eg = 0 ; nl1eg < NumL1EGammas() ; ++nl1eg)
				{
					L1Object l1eg = GetL1EGamma(nl1eg);
					if(l1eg.BX() == -1 && ph->DeltaR(l1eg) < 0.2)
					{
						//preprob *= 1.-prefireprob.FireProb(l1eg);
						preprob *= 1.-prefiretest_.FireProb(l1eg);
					}
				}
				_hists2d["prefire_all_ph"]->Fill(ph->Eta(), ph->Pt());
				_hists2d["prefire_bxm1_ph"]->Fill(ph->Eta(), ph->Pt(), 1.-preprob);
			}

			for(OJet* jet : RecPrefireAK4s)
			{   
				double preprob = 1.;
				for(size_t nl1eg = 0 ; nl1eg < NumL1EGammas() ; ++nl1eg)
				{
					L1Object l1eg = GetL1EGamma(nl1eg);
					if(l1eg.BX() == -1 && jet->DeltaR(l1eg) < 0.4)
					{
						//preprob *= 1.-prefireprob.FireProb(l1eg);
						preprob *= 1.-prefiretest_.FireProb(l1eg);
					}
				}
				_hists2d["prefire_all_jet"]->Fill(jet->Eta(), jet->Pt());
				_hists2d["prefire_bxm1_jet"]->Fill(jet->Eta(), jet->Pt(), 1.-preprob);
			}

		}
	}


	//Two lepton plots

	int nbjets = count_if(RecAK4s.begin(), RecAK4s.end(), [&](OJet* j) {return j->BTag() > btag_medium_;});
	if(RecEls.size() == 2 && RecEls[0]->Charge() != RecEls[1]->Charge() && (*RecEls[0] + *RecEls[1]).M() > 50.)
	{
		_hists1d["counter"]->Fill(6.5, weight);
		_hists1d["melel"]->Fill((*RecEls[0] + *RecEls[1]).M(), weight);
		_hists2d["melel_jet"]->Fill(RecAK4s.size()+0.5, (*RecEls[0] + *RecEls[1]).M(), weight);
		_hists2d["melel_bjet"]->Fill(nbjets+0.5, (*RecEls[0] + *RecEls[1]).M(), weight);
		_hists2d["melel_genbjet"]->Fill(GenBJets.size()+0.5, (*RecEls[0] + *RecEls[1]).M(), weight);
		if(nbjets == 1) {_hists1d["counter"]->Fill(7.5, weight);}
		if(nbjets == 2) {_hists1d["counter"]->Fill(8.5, weight);}
	}

	if(RecMus.size() == 2 && RecMus[0]->Charge() != RecMus[1]->Charge() && (*RecMus[0] + *RecMus[1]).M() > 50.)
	{
		_hists1d["counter"]->Fill(9.5, weight);
		TLorentzVector Z(*RecMus[0] + *RecMus[1]);
		_hists1d["mmumu"]->Fill(Z.M(), weight);
		_hists2d["mmumu_jet"]->Fill(RecAK4s.size()+0.5, Z.M(), weight);
		_hists2d["mmumu_bjet"]->Fill(nbjets+0.5, (*RecMus[0] + *RecMus[1]).M(), weight);
		_hists2d["mmumu_genbjet"]->Fill(GenBJets.size()+0.5, (*RecMus[0] + *RecMus[1]).M(), weight);
		if(nbjets == 1) {_hists1d["counter"]->Fill(10.5, weight);}
		if(nbjets == 2) {_hists1d["counter"]->Fill(11.5, weight);}
	}

}
