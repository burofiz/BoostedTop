#ifndef TOPJETS_H
#define TOPJETS_H

#include <NeutrinoSolver.h>
#include <OJet.h>

#include <TMatrixD.h>
#include <TVectorD.h>
#include <TRandom3.h>

class LepTopJet : public TLorentzVector
{
	private:
		OJet* jet_ = nullptr;
		TLorentzVector* met_ = nullptr;
		TLorentzVector* lep_ = nullptr;
		int lepPDGID_;
		double nnres_ = -1.;
		bool injet_ = false;
		TLorentzVector purejet_;
		TLorentzVector neutrino_;
		double dnu_;
		bool istop_ = false;
	public:
		LepTopJet(OJet* jet, TLorentzVector* met, TLorentzVector* lep, int lepPDGID, bool injet) : jet_(jet), met_(met), lep_(lep), lepPDGID_(lepPDGID), injet_(injet), purejet_(injet ? *jet-*lep : *jet)
	{

		TLorentzVector reco(purejet_ + *lep_);

		NeutrinoSolver ns(lep_, &purejet_);
		neutrino_ = ns.GetBest(met_->Px(), met_->Py(), 1., 1., 0., dnu_);
		if(dnu_ >= 0.)
		{
			reco += neutrino_;
			reco *= 0.97;
		}
	//	else
	//	{
	//		//neutrino_.SetXYZM(met_->Px(), met_->Py(), met_->Pt()/reco.Pt()*reco.Pz(), 0.);
	//		double s = 172.5/reco.M();
	//		reco *= s;
	//	}
		SetPxPyPzE(reco.Px(), reco.Py(), reco.Pz(), reco.E());
		//if(injet)
		//{
		//	SetPxPyPzE(tmp.Px() + jet_->Px(), tmp.Py() + jet_->Py(), tmp.Pz() + jet_->Pz(), tmp.Mag() + jet_->E());
		//}
		//else
		//{
		//	lepjet_ += *lep;
		//	SetPxPyPzE(tmp.Px() + lepjet_.Px(), tmp.Py() + lepjet_.Py(), tmp.Pz() + lepjet_.Pz(), tmp.Mag() + lepjet_.E());

		//}
	}
		TLorentzVector* Lep() const {return lep_;}
		TLorentzVector* MET() const {return met_;}
		double Dnu() const {return sqrt(dnu_);}
		const TLorentzVector& PureJet() const {return purejet_;}
		const TLorentzVector& Neutrino() const {return neutrino_;}
		OJet* Jet() const {return jet_;}
		int LepPDGID() const {return lepPDGID_;}
		double NNRes() const {return nnres_;}
		void NNRes(double nnres)  {nnres_ = nnres;}
		bool LepInJet() const {return injet_;}
		void IsTop(bool istop) {istop_ = istop;}
		bool IsTop() const {return istop_;}
		TLorentzVector LepJet() const {return PureJet()+*Lep();}
		template <class LEPTYPE> void SetTree(float* topdata)
		{
			LEPTYPE* lep = dynamic_cast<LEPTYPE*>(Lep());
			topdata[0] = LepJet().M()/200;
			topdata[1] = lep->Pt()/LepJet().Pt();
			topdata[2] = (LepJet() - *lep).M()/LepJet().M();
			topdata[3] = lep->IsoFarAll()/lep->IsoCentralAll()/0.2;
			topdata[4] = lep->IsoNearAll()/lep->IsoCentralAll()/400;
			//topdata[5] = (*jet_ + *met_ + *lep_).Mt()/200;
		}

};

class HadTopJet : public IOMergedSubJet, public TLorentzVector
{
	private:
		vector<const TLorentzVector*> jets_;
		double rescale_ = 1.;
		TLorentzVector JetSum_;
		TLorentzVector W_;
		TLorentzVector T_;
		double kmin = -1.;
		double nnres_ = -1.;
		vector<double> m2_;
		vector<double> m3_;
		double hardpairs_ = -1.;
		double hardjets_ = -1.;
		int matchedtpartons_ = 0;
		TLorentzVector matchedpt_;
		int matchedtbarpartons_ = 0;
		TLorentzVector matchedptbar_;
		bool needmergedjetst_ = true;
		bool needmergedjetstbar_ = true;
		bool ispshad_ = false;
		float btag_ = 0.;
		TVectorD sphv_;
		HadTopJet* genhtjet_ = nullptr;
		TLorentzVector PUSubJets(size_t n) const {return rescale_*SJet(SubJets(n));}
	public:
		static bool PS_MODE;
		template <class JETTYPE> HadTopJet(const IOMergedSubJet& htjet, const vector<JETTYPE*>& jets) : IOMergedSubJet(htjet), TLorentzVector(px(), py(), pz(), e()), JetSum_(px(), py(), pz(), e()), sphv_(3)
	{


		for(size_t mj = 0 ; mj < Num_MemberJets() ; ++mj)
		{
			JETTYPE* mjet = jets[MemberJets(mj)];
			jets_.push_back(mjet);
			needmergedjetst_ = needmergedjetst_ && mjet->matchedtpartons > 0;
			needmergedjetstbar_ = needmergedjetstbar_ && mjet->matchedtbarpartons > 0;
			matchedtpartons_ += mjet->matchedtpartons;
			matchedpt_ += mjet->ptpartons;
			matchedtbarpartons_ += mjet->matchedtbarpartons;
			matchedptbar_ += mjet->ptbarpartons;
		}
		Init();
	}

		void Init();
		void SetTree(float* tophaddata);
		void IsPS(bool ispshad) {ispshad_ = ispshad;}
		bool IsPS() const {return ispshad_;}
		bool IsCleanHadTopJet() const {return (PS_MODE && IsPS()) || (!PS_MODE && ((needmergedjetst_ && MatchedTPartons() == 3 && MatchedTBarPartons() == 0) || (needmergedjetstbar_ && MatchedTPartons() == 0 && MatchedTBarPartons() == 3)));}
		bool IsHadTopJet() const {return (PS_MODE && IsPS()) || (!PS_MODE && ((MatchedTPartons() == 3 && MatchedTBarPartons() == 0) || (MatchedTPartons() == 0 && MatchedTBarPartons() == 3)));}
		bool IsCleanBrokenHadTopJet() const {return !PS_MODE && ((needmergedjetst_ && MatchedTPartons() == 2 && MatchedTBarPartons() == 0) || (needmergedjetstbar_ && MatchedTPartons() == 0 && MatchedTBarPartons() == 2));}
		bool IsBrokenHadTopJet() const {return !PS_MODE && (MatchedTPartons() == 2 || MatchedTBarPartons() == 2);}
		bool IsSignalJet() const {return IsHadTopJet() || IsBrokenHadTopJet();}

		bool Overlap(const HadTopJet& b);
		bool  CloseToMember(const TLorentzVector* b, double dr = 0.8) const;
		int MatchedTPartons() const {return matchedtpartons_;}	
		const TLorentzVector& MatchedMomT() const {return matchedpt_;}	
		int MatchedTBarPartons() const {return matchedtbarpartons_;}
		const TLorentzVector& MatchedMomTBar() const {return matchedptbar_;}	
		bool MatchPartons() const {return matchedtbarpartons_ > 1  || matchedtpartons_ > 1;}
		int NumMatchPartons() const {return max(matchedtbarpartons_, matchedtpartons_);}
		bool IsTJet() const {return matchedtpartons_ == 3;}
		bool IsTBarJet() const {return matchedtbarpartons_ == 3;}
		const TLorentzVector& Jets() const {return JetSum_;}
		const TLorentzVector& W() const {return W_;}
		const TLorentzVector& T() const {return T_;}
		double NNRes() const {return nnres_;}
		void NNRes(double nnres)  {nnres_ = nnres;}
		bool NNCut() const {return NNRes() > 0.5;};
		const vector<double>& M3s() const {return m3_;}
		const vector<double>& M2s() const {return m2_;}
		double HardPairs() const {return hardpairs_;}
		double HardJets() const {return hardjets_;}
		double Sphericity() const {return 3./2*(sphv_(1)+sphv_(2));}
		double Aplanarity() const {return 3./2*sphv_(2);}
		double Frac(int n) {TLorentzVector tlv(PUSubJets(n));  tlv.Boost(BoostVector()); return tlv.Pt()/Pt();}
		void BTag(float btag) {btag_ = btag;}
		float BTag() const {return btag_;}
		HadTopJet* MatchedGenJet() const {return genhtjet_;};

		double ScaleRes(double chadtopscale_, double chadtopres_, HadTopJet* genhtjet);
};

#endif
