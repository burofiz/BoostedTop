#include "BoostedTop.h"

void BoostedTop::ElJetTriggerAnalysis()
{

	int nbjets = count_if(RecAK4s.begin(), RecAK4s.end(), [&](OJet* j) {return j->BTag() > btag_medium_;});
	double nlweight = weight/lepweight_;
	if(nbjets >= 2)
	{
		if(MuIso.GetResult(RecMus[0]) == 1)
		{
			for(OMuon* mu : RecMusLoose)
			{
				if(mu->IsoNearAll()/mu->IsoCentralAll() > 50.) {continue;}
				_hists2d["trmu_eta_pt_all"]->Fill(mu->Eta(), mu->Pt(), nlweight);
				if(Mu50.GetResult(mu) == 1)
				{
					_hists2d["trmu_eta_pt_passnoniso"]->Fill(mu->Eta(), mu->Pt(), nlweight);
				}
				if(Mu50.GetResult(mu) == 1 || MuIso.GetResult(mu) == 1 || MuTkIso.GetResult(mu) == 1)
				{
					_hists2d["trmu_eta_pt_pass"]->Fill(mu->Eta(), mu->Pt(), nlweight);
				}
			}
		}
	}
	int trres = trgres["ElJet"].GetResult();
	if(abs(trres) == 1)
	{
		int fels = 0;
		//if(nbjets >= 2)
		{
			for(OElectron* el : RecElsLoose)
			{
				if(el->IsoNearAll()/el->IsoCentralAll() > 50.) {continue;}
				_hists2d["trel_eta_pt_all"]->Fill(el->Eta(), el->Pt(), nlweight);
				if(ElJet_El.GetResult(el) == 1)
				{
					_hists2d["trel_eta_pt_passnoniso"]->Fill(el->Eta(), el->Pt(), nlweight);
					fels++;
				}
				if(ElJet_El.GetResult(el) == 1 || ElIso32.GetResult(el) == 1 || ElIso27.GetResult(el) == 1)
				{
					_hists2d["trel_eta_pt_pass"]->Fill(el->Eta(), el->Pt(), nlweight);
				}
				if(ElIso32.GetResult(el) == 1 || ElIso27.GetResult(el) == 1)
				{
					_hists2d["trel_eta_pt_passiso"]->Fill(el->Eta(), el->Pt(), nlweight);
				}
			}
		}		

		for(OElectron* el : RecEls)
		{

			_hists2d["treliso_eta_pt_all"]->Fill(el->Eta(), el->Pt(), nlweight);
			if(ElJet_El.GetResult(el) == 1)
			{
				_hists2d["treliso_eta_pt_passnoniso"]->Fill(el->Eta(), el->Pt(), nlweight);
				fels++;
			}
			if(ElJet_El.GetResult(el) == 1 || ElIso32.GetResult(el) == 1 || ElIso27.GetResult(el) == 1)
			{
				_hists2d["treliso_eta_pt_pass"]->Fill(el->Eta(), el->Pt(), nlweight);
			}

		}		
		if(fels != 0)
		{
			for(OJet* jet : RecAK4s)
			{
				_hists2d["trjet_eta_pt_all"]->Fill(jet->Eta(), jet->Pt(), nlweight);
				if(ElJet_Jet.GetResult(jet) == 1)
				{
					_hists2d["trjet_eta_pt_pass"]->Fill(jet->Eta(), jet->Pt(), nlweight);
				}
			}
		}
	}

}
