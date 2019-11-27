#ifndef TTEVENT
#define TTEVENT
#include <TLorentzVector.h>

#include <NeutrinoSolver.h>
#include <Analyse.h>

#include "TopJets.h"

class TTEvent
{
	private:
		bool recal_ = true;
		bool thadset_ = false;
		TLorentzVector* thad_ = nullptr;
		TLorentzVector* bhad_ = nullptr;
		TLorentzVector* jahad_ = nullptr;
		TLorentzVector* jbhad_ = nullptr;
		bool tlepset_ = false;
		TLorentzVector* tlep_ = nullptr;
		TLorentzVector* blep_ = nullptr;
		bool nlepset_ = false;
		TLorentzVector* nlep_ = nullptr;
		TLorentzVector* met_ = nullptr;
		TLorentzVector* llep_ = nullptr;
		TLorentzVector tt_;
		TLorentzVector t_cms_;
		TLorentzVector tb_cms_;
		vector<const TLorentzVector*> addjets;
		vector<const HadTopJet*> addhtjets;
		int lpdgid_ = 0;
		double dnu_ = numeric_limits<double>::max();
		double prob_ = numeric_limits<double>::max();
		double probhad_ = numeric_limits<double>::max();
		double problep_ = numeric_limits<double>::max();
		double ht_ = 0;
		double cts_ = 0;
		double bstar_ = 0;
		double bbarstar_ = 0;
		double deltabbstar_ = 0;

		void copy(const TTEvent& ttev)
		{
			recal_ = true;
			thadset_ = ttev.thadset_;
			thad_ = ttev.thad_;
			bhad_ = ttev.bhad_;
			jahad_ = ttev.jahad_;
			jbhad_ = ttev.jbhad_;
			tlepset_ = ttev.tlepset_;
			tlep_ = ttev.tlep_;
			blep_ = ttev.blep_;
			nlepset_ = ttev.nlepset_;
			nlep_ = ttev.nlep_;
			met_ = ttev.met_;
			llep_ = ttev.llep_;
			tt_ = ttev.tt_;
			t_cms_ = ttev.t_cms_;
			tb_cms_ = ttev.tb_cms_;
			lpdgid_ = ttev.lpdgid_;
			dnu_ = ttev.dnu_;
			prob_ = ttev.prob_;
			probhad_ = ttev.probhad_;
			problep_ = ttev.problep_;
			ht_ = ttev.ht_;
			cts_ = ttev.cts_;
			bstar_ = ttev.bstar_;
			bbarstar_ = ttev.bbarstar_;
			deltabbstar_ = ttev.deltabbstar_;
			addjets = ttev.addjets;
			if(!tlepset_) tlep_ = nullptr;
			if(!thadset_) thad_ = nullptr;
			if(!nlepset_) nlep_ = nullptr;
		}
	public:
		TTEvent() = default;
		TTEvent(const TTEvent& ttev)
		{
			copy(ttev);
		}

		TTEvent& operator=(const TTEvent& ttev)
		{
			copy(ttev);
			return *this;
		}

		void Reset()
		{
			if(!thadset_ && thad_ != nullptr) {delete thad_;}
			if(!tlepset_ && tlep_ != nullptr) {delete tlep_;}
			if(!nlepset_ && nlep_ != nullptr) {delete nlep_;}
			recal_ = true;
			thadset_ = false;
			thad_ = nullptr;
			bhad_ = nullptr;
			jahad_ = nullptr;
			jbhad_ = nullptr;
			tlepset_ = false;
			tlep_ = nullptr;
			blep_ = nullptr;
			nlepset_ = false;
			nlep_ = nullptr;
			met_ = nullptr;
			llep_ = nullptr;
			lpdgid_ = 0;
			dnu_ = numeric_limits<double>::max();
			prob_ = 0.;
			probhad_ = 0.;
			problep_ = 0.;
			ht_ = 0.;
			cts_ = 0.;
			
			addjets.clear(); 
		}
		~TTEvent()
		{
			if(!thadset_ && thad_ != nullptr) {delete thad_;}
			if(!tlepset_ && tlep_ != nullptr) {delete tlep_;}
			if(!nlepset_ && nlep_ != nullptr) {delete nlep_;}
		}
		void THad(TLorentzVector* thad){recal_= true; thad == nullptr ? thadset_ = false : thadset_ = true; thad_ = thad;}
		bool THadSet() const {return thadset_;}
		void BHad(TLorentzVector* bhad){recal_= true; bhad_ = bhad;}
		void JaHad(TLorentzVector* jahad){recal_= true; jahad_ = jahad;}
		void JbHad(TLorentzVector* jbhad){recal_= true; jbhad_ = jbhad;}
		void TLep(TLorentzVector* tlep){recal_= true; tlep == nullptr ? tlepset_ = false : tlepset_ = true; tlep_ = tlep;}
		bool TLepSet() const {return tlepset_;}
		void BLep(TLorentzVector* blep){recal_= true; blep_ = blep;}
		void LLep(TLorentzVector* llep, int lpdgid){recal_= true; llep_ = llep; lpdgid_ = lpdgid;}
		void MET(TLorentzVector* met){recal_= true; met_ = met;}
		void NLep(TLorentzVector* nlep){recal_= true; nlepset_ = true; nlep_ = nlep;}
		bool NLepSet() const {return nlepset_;}
		template<typename T> void SetAdditionalJets(const vector<T*>& jets)
		{   
			addjets.clear();
			ht_ = 0.;
			for(const T* jet : jets)
			{   
				if(IsJetIn(jet) == -1)
				{
					addjets.push_back(jet);
					ht_ += jet->Pt();
				}
			}
			sort(addjets.begin(), addjets.end(), [](const TLorentzVector* A, const TLorentzVector* B){return A->Pt() > B->Pt();});  
		}
		void SetAdditionalHTJets(const vector<HadTopJet*>& jets)
		{   
			addhtjets.clear();
			for(const HadTopJet* jet : jets)
			{   
				if(IsCompleteBL())
				{
					if(jet->CloseToMember(TLep())) {continue;}
				}
				if(IsCompleteBH())
				{
					if(jet->DeltaR(*THad()) < 0.4) {continue;}
				}
				if(IsCompleteRL())
				{
					if(jet->CloseToMember(BLep())) {continue;}
				}
				if(IsCompleteRH())
				{
					if(jet->CloseToMember(BHad())) {continue;}
					if(jet->CloseToMember(WJPtmax())) {continue;}
					if(jet->CloseToMember(WJPtmin())) {continue;}
				}
				addhtjets.push_back(jet);
			}
			sort(addhtjets.begin(), addhtjets.end(), [](const HadTopJet* A, const HadTopJet* B){return A->NNRes() > B->NNRes();});  
		}

		const vector<const TLorentzVector*>& AddJets() const {return addjets;}
        const TLorentzVector* GetJet(size_t n) const
        {   
            if(n == 0) return BLep();
            if(n == 1) return BHad();
            if(n == 2) return WJPtmax();
            if(n == 3) return WJPtmin();
            if(n > 3 && n-4 < addjets.size()) return addjets[n-4];
            
            return nullptr;
        }

        const TLorentzVector* GetJet(const TLorentzVector* jet) const
        {   
            int nj = IsJetIn(jet);
            return(GetJet(nj));
        }

        const vector<const HadTopJet*>& AddHTJets() const {return addhtjets;}
        int NAddJets() const {return addjets.size();}
        int NJets() const {return 4+NAddJets();}

		void CalculateTLep()
		{
			if(tlepset_ || !IsCompleteRL()) return;


			if(!nlepset_)
			{
				if(nlep_ != nullptr)
				{
					delete nlep_;	
				}
				NeutrinoSolver ns(llep_, blep_);
				nlep_ = new TLorentzVector(ns.GetBest(met_->Px(), met_->Py(), 1., 1., 0., dnu_));
			}
			if(tlep_ != nullptr)
			{
				delete tlep_;	
			}
			tlep_ = new TLorentzVector(*blep_ + *nlep_ + *llep_);

		}

		void CalculateTHad()
		{
			if(thadset_ || !IsCompleteRH()) return;

			if(thad_ != nullptr)
			{
				delete thad_;	
			}
			thad_ = new TLorentzVector(*bhad_ + *jahad_ + *jbhad_);
			if(GLAN->IsData() && GLAN->UserTime() == 2017)
			{
				double s = 0.99925+0.0091582/(1+exp(-0.0239413*(thad_->Pt()-184.856)));
				*thad_ *= 1/s;
			}
			//else if(GLAN->IsData() && GLAN->UserTime() == 2018)
			//{
			//	double pt = min(500., thad_->Pt());
			//	double s = 0.99745 + 3.88166E-5*pt;
			//	*thad_ *= 1/s;
			//}

		}

		int Calculate()
		{
			if(!recal_) return -1;
			recal_= false;
			CalculateTLep();
			CalculateTHad();
			
			if(tlep_ != nullptr && thad_ != nullptr)
			{
				tt_ = (*tlep_ + *thad_);
				t_cms_ = *T();
				t_cms_.Boost(-1.*TT()->BoostVector());
				tb_cms_ = *TBar();
				tb_cms_.Boost(-1.*TT()->BoostVector());
				cts_ = TT()->Vect().Dot(t_cms_.Vect())/TT()->P()/t_cms_.P();

				if(TB() != nullptr && TBarB() != nullptr)
				{
					TLorentzVector tb(*TB());
					tb.Boost(-1.*TT()->BoostVector());
					TVector3 tbtrans = tb.Vect() - t_cms_.Vect().Dot(tb.Vect())/t_cms_.P()/t_cms_.P()*t_cms_.Vect();

					TLorentzVector tbarb(*TBarB());
					tbarb.Boost(-1.*TT()->BoostVector());
					TVector3 tbarbtrans = tbarb.Vect() - t_cms_.Vect().Dot(tbarb.Vect())/t_cms_.P()/t_cms_.P()*t_cms_.Vect();

					TVector3 zaxis = tt_.Vect().Cross(t_cms_.Vect());
					zaxis = zaxis * (1./zaxis.Mag());

					bstar_ = zaxis.Dot(tbtrans)/tbtrans.Mag();
					bbarstar_ = zaxis.Dot(tbarbtrans)/tbarbtrans.Mag();
					deltabbstar_ = tbtrans.Dot(tbarbtrans)/tbtrans.Mag()/tbarbtrans.Mag();

				}

				//TLorentzVector ba(0,0,1,1);
				//TLorentzVector bb(0,0,-1,1);
				//ba.Boost(-1.*TT()->BoostVector());
				//bb.Boost(-1.*TT()->BoostVector());
				//ba = ba*(1./ba.P());
				//bb = bb*(1./bb.P());
				//
				//TVector3 bv(-1.*TT()->BoostVector());
				//cout << bv.Px() << " " << bv.Py() << " " <<bv.Pz() << endl;
				//cout << ba.Px() << " " << ba.Py() << " " <<ba.Pz() << " " <<ba.E() << endl;
				//cout << bb.Px() << " " << bb.Py() << " " <<bb.Pz() << " " <<bb.E() << endl;
				//TLorentzVector babb(ba+bb);
				//babb = babb*(1./babb.P());
				//cout << babb.Px() << " " << babb.Py() << " " <<babb.Pz() << " " <<babb.E() << endl;
				

				return 2;
			}
			return 0;
		}

		void Prob(double probhad, double problep) {probhad_ = -1.*log(probhad); problep_ = -1.*log(problep); prob_ = probhad_+problep_;}
		double Prob() const {return prob_;}
		double ProbHad() const {return probhad_;}
		double ProbLep() const {return problep_;}

		TLorentzVector* T() const {return LPDGID() < 0 ? TLep() : THad();}
		TLorentzVector* TB() const {return LPDGID() < 0 ? BLep() : BHad();}
		TLorentzVector* TBar() const {return LPDGID() < 0 ? THad() : TLep();}
		TLorentzVector* TBarB() const {return LPDGID() < 0 ? BHad() : BLep();}
		TLorentzVector* THad() const {return thad_;}
		TLorentzVector* THard() const {return THad()->Pt() > TLep()->Pt() ? THad() : TLep();}
		TLorentzVector* TSoft() const {return THad()->Pt() > TLep()->Pt() ? TLep() : THad();}
		TLorentzVector* BHad() const {return bhad_;}
		TLorentzVector* JaHad() const {return jahad_;}
		TLorentzVector* JbHad() const {return jbhad_;}
		TLorentzVector* WJPtmax() const {return jbhad_->Pt() > jahad_->Pt() ? jbhad_ : jahad_;}
		TLorentzVector* WJPtmin() const {return jbhad_->Pt() > jahad_->Pt() ? jahad_ : jbhad_;}
		TLorentzVector* TLep() const {return tlep_;}
		TLorentzVector* BLep() const {return blep_;}
		TLorentzVector* LLep() const {return llep_;}
		TLorentzVector* NLep() const {return nlep_;}
		TLorentzVector* MET() const {return met_;}
		const TLorentzVector* TT() const {return &tt_;}
		double Dnu() const {return sqrt(dnu_);}
		int LPDGID() const {return lpdgid_;}
		int IsJetIn(const TLorentzVector* jet) const
		{
			if(IsCompleteBL())
			{
				if(jet->DeltaR(*TLep()) < 0.6) return -2;
			}
			if(IsCompleteBH())
			{
				HadTopJet* th = dynamic_cast<HadTopJet*>(THad());
				if(th != nullptr)
				{
					if(th->CloseToMember(jet)) {return -3;}	
				}
				else
				{
					if(jet->DeltaR(*THad()) < 1.2) {return -3;}
				}
			}
			if(IsCompleteRL())
			{
				if(jet->DeltaR(*BLep()) < 0.4) return 0;
			}
			if(IsCompleteRH())
			{
				if(jet->DeltaR(*BHad()) < 0.4) return 1;
				if(jet->DeltaR(*WJPtmax()) < 0.4) return 2;
				if(jet->DeltaR(*WJPtmin()) < 0.4) return 3;
			}
			for(size_t i = 0 ; i < addjets.size() ; ++i)
			{
				if(jet->DeltaR(*addjets[i]) < 0.4) return i+4;
			}

			return -1;
		}
		double DRmaxHad() 	
		{
			return max({bhad_->DeltaR(*jahad_), bhad_->DeltaR(*jbhad_), jahad_->DeltaR(*jbhad_)});
		}
		double DRminHad() 	
		{
			return min({bhad_->DeltaR(*jahad_), bhad_->DeltaR(*jbhad_), jahad_->DeltaR(*jbhad_)});
		}
		//bool BoostedHad() 	
		//{
		//	return DRmaxHad() < 1.2;
		//}
		//bool ResolvedHad() 	
		//{
		//	return DRminHad() > 0.4;
		//}
		//bool BoostedLep() 	
		//{
		//	return llep_->DeltaR(*blep_) < 0.6;
		//}
		bool IsCompleteRL() const 
		{
			if((met_ == nullptr && nlep_ == nullptr) || llep_ == nullptr || blep_ == nullptr) return false;
			return true;
		}
		bool IsCompleteRH() const 
		{
			if( bhad_ == nullptr || jahad_ == nullptr || jbhad_ == nullptr) return false;
			if(bhad_ == jahad_) return false;
			if(bhad_ == jbhad_) return false;
			if(jbhad_ == jahad_) return false;
			return true;
		}
		bool IsCompleteRES() const {return IsCompleteRL() && IsCompleteRH() && blep_ != bhad_ && blep_ != jahad_ && blep_ != jbhad_;}
		//bool IsCompleteRES() const {return IsCompleteRL() && IsCompleteRH();}
		bool IsComplete() const {return IsCompleteBHBL() || IsCompleteRES() || (IsCompleteRL() && IsCompleteBH()) || (IsCompleteRH() && IsCompleteBL());}
		bool IsCompleteBH() const {return THadSet();}
		bool IsCompleteBL() const {return TLepSet();}
		bool IsCompleteBHBL() const {return IsCompleteBH() && IsCompleteBL();}
		bool IsCompleteHad() const {return THadSet() || IsCompleteRH();}
		bool IsCompleteLep() const {return TLepSet() || IsCompleteRL();}

		bool IsTHadCorrectKin(const TTEvent& ttev) const
		{
			if(!IsCompleteRH() || !ttev.IsCompleteRH())
			{
				double rm = 0.4;
				return THad()->DeltaR(*ttev.THad()) < rm;
			}
			else
			{
				double rm = 0.4;
				if(BHad()->DeltaR(*ttev.BHad()) > rm) {return false;}
				if(JaHad()->DeltaR(*ttev.JaHad()) < rm && JbHad()->DeltaR(*ttev.JbHad()) < rm) {return true;}
				if(JaHad()->DeltaR(*ttev.JbHad()) < rm && JbHad()->DeltaR(*ttev.JaHad()) < rm) {return true;}
				return false;
			}
		}

		bool IsTLepCorrectKin(const TTEvent& ttev) const
		{
			if(!IsCompleteRL() || !ttev.IsCompleteRL())
			{
				double rm = 0.4;
				return TLep()->DeltaR(*ttev.TLep()) < rm;
			}
			else
			{
				double rm = 0.4;
				if(BLep()->DeltaR(*ttev.BLep()) > rm) {return false;}
				if(LLep()->DeltaR(*ttev.LLep()) > 0.2) {return false;}
				return true;
			}
		}

		bool IsCorrectKin(const TTEvent& ttev) const
		{
			return IsTLepCorrectKin(ttev) && IsTHadCorrectKin(ttev); 
		}

		double Ht() const {return ht_;}

		double CosThetaCMS() const{return cts_;}
		double BPhiCMS() const{return bstar_;}
		double BBarPhiCMS() const{return bbarstar_;}
		double DeltaBBPhiCMS() const{return deltabbstar_;}

//		bool IsTHadCorrect(const TTEvent& ttev) const
//		{
//			if(thadset_ && ttev.thadset_ && thad_ == ttev.thad_) {return true;}
//			if(!thadset_ && !ttev.thadset_ && bhad_ == ttev.bhad_ && ((jahad_ == ttev.jahad_ && jbhad_ == ttev.jbhad_) || (jbhad_ == ttev.jahad_ && jahad_ == ttev.jbhad_))) {return true;}
//
//			return false;
//		}
//		bool IsTLepCorrect(const TTEvent& ttev) const
//		{
//			if(tlepset_ && ttev.tlepset_ && tlep_ == ttev.tlep_) {return true;}
//			if(!tlepset_ && !ttev.tlepset_ && blep_ == ttev.blep_ && llep_ == ttev.llep_) {return true;}
//
//			return false;
//		}
//		bool IsCorrect(const TTEvent& ttev)
//		{
//			if(IsTHadCorrect(ttev) && IsTLepCorrect(ttev)) return true;
//			return false;
//		}
};

#endif
