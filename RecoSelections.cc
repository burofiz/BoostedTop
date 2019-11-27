#include "BoostedTop.h"

void BoostedTop::SelectRECOParticles()
{
	if(ran_SelectRECOParticles_ == true){return;}
	ran_SelectRECOParticles_ = true;

	for(UInt_t i = 0 ; i < NumIOMuons() ; i++)
	{   
		OMuon mu(GetIOMuon(i));
		if(mu.Pt() > 15. && abs(mu.Eta()) < 2.4)
		{ 
			SMuons.push_back(mu);
			if(mu.ID(OMuon::TIGHT, true))
			{
				RecMusVeto.push_back(&SMuons.back());
			}

			if(mu.Pt() > 30. && mu.ID(OMuon::TIGHT, true))
			{   
				RecMus.push_back(&SMuons.back());
			}
			else if(((!ISOSIDEBAND_ && mu.Pt() > 50. && mu.IsoNearAll()/mu.IsoCentralAll() < 200.) || (ISOSIDEBAND_ && mu.Pt() > 30.)) && mu.ID(OMuon::TIGHT, false))
			{
				RecMusLoose.push_back(&SMuons.back());
			}
		}
	}

	bool hemaffected = false;
	for(UInt_t i = 0 ; i < NumIOElectrons() ; i++)
	{   
		OElectron el(GetIOElectron(i));
		if(el.Pt() > 15. && abs(el.Eta()) < 2.4)
		{
			//Excluding a phi-region in endcap due to HEM 15/16 issue
			if(UserTime() == 2018 && el.Eta() < -1.5 && el.Phi() > -1.5 && el.Phi() < -0.9)
			{
				if(IsData() && Run() >= 319077) {continue;}
				hemaffected = true;
			}

			SElectrons.push_back(el);
			if(el.ID(OElectron::TIGHT, true))
			{
				RecElsVeto.push_back(&SElectrons.back());
			}

			if(((UserTime() == 2018 && el.Pt() > 34.) || (UserTime() != 2018 && el.Pt() > 30.)) && el.ID(OElectron::TIGHT, true))
			{
				RecEls.push_back(&SElectrons.back());
			}
			else if(el.Pt() > 50. && el.ID(OElectron::TIGHT, false) && el.IsoNearAll()/el.IsoCentralAll() < 200.)
			{
				RecElsLoose.push_back(&SElectrons.back());
			}
		}
	}
	if(hemaffected) {weight*=0.3494;}

	for(UInt_t i = 0 ; i < NumIOPhotons() ; i++)
	{   
		OPhoton ph(GetIOPhoton(i));
		SPhotons.push_back(ph);
		if(ph.Pt() > 30. && abs(ph.Eta()) < 2.4 && ph.ID(OPhoton::MEDIUM))
		{
			RecPhs.push_back(&SPhotons.back());
		}
		if(ph.Pt() > 25. && abs(ph.Eta()) < 3.0 && ph.ID(OPhoton::LOOSE, false))
		{
			RecPrefirePhs.push_back(&SPhotons.back());
		}
	}


	TLorentzVector metcor;
	for(UInt_t i = 0 ; i < NumIOPFAK4Jets() ; i++)
	{   
		OJet jet(GetIOPFAK4Jet(i));

		if(IsMC())
		{
			double sf = jetscaler.GetRes(&jet, GenAllJets, Rho(), cuncjer_);
			sf *= jetscaler.GetScale(&jet, GenBJets, cuncjes_);
			metcor += jet.ApplySF(sf);
		}

		SAK4s.push_back(jet);
		if(jet.ID() && jet.Pt() > 30. && abs(jet.Eta()) < 2.4 && jet.Clean(RecMusVeto, RecElsVeto, RecPhs, 0.4))
		{
			RecAK4s.push_back(&SAK4s.back());
		}
		if(jet.Pt() > 25. && abs(jet.Eta()) < 3.0 && jet.Clean(RecMusVeto, {}, RecPrefirePhs, 0.4))
		{
			RecPrefireAK4s.push_back(&SAK4s.back());
		}
	}

	int mettype = (UserTime() == 2017 ? 1 : 0);
	TVector3 metv(GetIOMET(mettype).px()-metcor.Px(), GetIOMET(mettype).py()-metcor.Py(), 0.);
	met.SetPxPyPzE(metv.Px(), metv.Py(), 0., metv.Mag());


}

void BoostedTop::SelectRECOParticles2023()
{
	if(ran_SelectRECOParticles_ == true){return;}
	ran_SelectRECOParticles_ = true;


	for(UInt_t i = 0 ; i < NumIOMuons() ; i++)
	{   
		OMuon mu(GetIOMuon(i));
		if(mu.Pt() > 15. && abs(mu.Eta()) < 2.8)
		{   
			SMuons.push_back(mu);
			if((DELPHES && mu.IsoDelphes() < 0.15) || mu.ID(OMuon::TIGHT, true))
			{
				RecMusVeto.push_back(&SMuons.back());
			}

			if(ttgen.IsComplete() && ttgen.LLep()->DeltaR(mu) < 0.1)
			{
				if(DELPHES) {_hists1d["muon_iso"]->Fill(mu.IsoDelphes());}
				else{_hists1d["muon_iso"]->Fill(mu.IsoDB());}
			}

			if(mu.Pt() > 30. && ((DELPHES && mu.IsoDelphes() < 0.15) || mu.ID(OMuon::TIGHT, true)))
			{   
				RecMus.push_back(&SMuons.back());
			}
			else if(mu.Pt() > 60. && mu.ID(OMuon::TIGHT, false))
			{
				RecMusLoose.push_back(&SMuons.back());
			}
		}
	}

	for(UInt_t i = 0 ; i < NumIOElectrons() ; i++)
	{   
		OElectron el(GetIOElectron(i));
		if(el.Pt() > 15. && abs(el.Eta()) < 2.8)
		{
			SElectrons.push_back(el);
			if((DELPHES && el.IsoDelphes() < 0.15)|| el.ID(OElectron::TIGHT, true))
			{
				RecElsVeto.push_back(&SElectrons.back());
			}

			if(ttgen.IsComplete() && ttgen.LLep()->DeltaR(el) < 0.1)
			{
				if(DELPHES){_hists1d["el_iso"]->Fill(el.IsoDelphes());}
				else{_hists1d["el_iso"]->Fill(el.IsoDB());}
			}

			if(el.Pt() > 30. && ((DELPHES && el.IsoDelphes() < 0.15) || el.ID(OElectron::TIGHT, true)))
			{
				RecEls.push_back(&SElectrons.back());
			}
			else if(el.Pt() > 60. && el.ID(OElectron::LOOSE, false))
			{
				RecElsLoose.push_back(&SElectrons.back());
			}
		}
	}

	TLorentzVector metcor;

	for(UInt_t i = 0 ; i < NumIOPFAK4Jets() ; i++)
	{   
		OJet jet(GetIOPFAK4Jet(i));
		if(IsMC())
		{
			double sf = delphesjec.GetScale(&jet, GenBJets, cuncjes_);
			metcor += jet.ApplySF(sf);
		}
		if(jet.Pt() < 30. || abs(jet.Eta()) > 4.0 || !jet.Clean(RecMusVeto, RecElsVeto, RecPhs, 0.4)) {continue;}
		SAK4s.push_back(jet);
		RecAK4s.push_back(&SAK4s.back());
	}
	TVector3 metv(GetIOMET(0).px()-metcor.Px(), GetIOMET(0).py()-metcor.Py(), 0.);
	met.SetPxPyPzE(metv.Px(), metv.Py(), 0., metv.Mag());

}
