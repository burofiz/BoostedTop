#include "TopJets.h"


void HadTopJet::Init()
{
	m2_.clear();
	m3_.clear();
	hardjets_ = 0;

	for(size_t nja = 0 ; nja < Num_SubJets() ; ++nja)
	{
		if(PUSubJets(nja).E() > 15.) {hardjets_++;}
		for(size_t njb = 0 ; njb < nja ; ++njb)
		{   
			TLorentzVector tlv2 = PUSubJets(nja) + PUSubJets(njb);
			m2_.push_back(tlv2.M());
			for(size_t njc = 0 ; njc < njb ; ++njc)
			{
				TLorentzVector tlv3 = tlv2 + PUSubJets(njc);
				m3_.push_back(tlv3.M());
			}

			for(size_t njc = 0 ; njc < Num_SubJets() ; ++njc)
			{   
				if(nja == njc || njb == njc) {continue;}
				TLorentzVector tlv3 = tlv2 + PUSubJets(njc);
				double k = pow(tlv3.M()-172.5, 2) + pow(tlv2.M()-81.4, 2);
				if(kmin == -1 || k < kmin)
				{
					kmin = k;
					W_ = tlv2;
					T_ = tlv3;
				}	
			}
		}
	}
	//TVector3 bv(JetSum_.BoostVector());
	//SetPxPyPzE(T_.Px(), T_.Py(), T_.Pz(), T_.E());
	//Boost(bv);
	sort(m2_.begin(), m2_.end(), greater<double>());
	sort(m3_.begin(), m3_.end(), greater<double>());
	hardpairs_ = count_if(m2_.begin(), m2_.end(), [](const double& val){return val > 40.;});

	TMatrixD sph(3,3);
	sph.Zero();
	double sump=0.;
	for(size_t nja = 0 ; nja < Num_SubJets() ; ++nja)
	{
		TLorentzVector tlv(PUSubJets(nja));
		sph(0, 0) += tlv.Px()*tlv.Px()/tlv.P();
		sph(0, 1) += tlv.Px()*tlv.Py()/tlv.P();
		sph(0, 2) += tlv.Px()*tlv.Pz()/tlv.P();
		sph(1, 0) += tlv.Py()*tlv.Px()/tlv.P();
		sph(1, 1) += tlv.Py()*tlv.Py()/tlv.P();
		sph(1, 2) += tlv.Py()*tlv.Pz()/tlv.P();
		sph(2, 0) += tlv.Pz()*tlv.Px()/tlv.P();
		sph(2, 1) += tlv.Pz()*tlv.Py()/tlv.P();
		sph(2, 2) += tlv.Pz()*tlv.Pz()/tlv.P();
		sump += tlv.P();
	}
	sph *=1./sump;
	sph.EigenVectors(sphv_);
	//cout << sphv_(0) << " " << sphv_(1) << " " << sphv_(2) << endl;
}

void HadTopJet::SetTree(float* tophaddata)
{
	tophaddata[0] = m3_[0]/100;
	tophaddata[1] = Sphericity();
	tophaddata[2] = Num_SubJets() > 3 ? m3_[1]/100 : 0.;
	tophaddata[3] = Num_SubJets() > 3 ? m3_[1]/m3_[0] : 0.;
	tophaddata[4] = m2_[0]/100;
	tophaddata[5] = m2_[1]/100;
	tophaddata[6] = m2_[2]/100;
	tophaddata[7] = HardPairs()/10.;
	tophaddata[8] = Taus(1)/Taus(0);
	tophaddata[9] = Taus(2)/Taus(1);
	tophaddata[10] = Num_Taus() > 3 ? Taus(3)/Taus(2) : 0.;
	tophaddata[11] = Num_Taus() > 4 ? Taus(4)/Taus(3) : 0.;
	tophaddata[12] = PUSubJets(0).E()/100.;
	tophaddata[13] = PUSubJets(1).E()/100.;
	tophaddata[14] = PUSubJets(2).E()/100.;
	//tophaddata[15] = Aplanarity();
	tophaddata[15] = Num_SubJets() > 3 ? PUSubJets(3).E()/100. : 0.;
	tophaddata[16] = m3_[1]/JetSum_.M();
	TVector3 cross(PUSubJets(0).Vect().Cross(PUSubJets(1).Vect()));
	cross *= 1./cross.Mag();
	tophaddata[17] = abs(cross * PUSubJets(2).Vect())/PUSubJets(2).P();
	tophaddata[18] = Frac(0);
	tophaddata[19] = Frac(1);
	tophaddata[20] = Frac(2);

}		

bool HadTopJet::Overlap(const HadTopJet& b)
{
	for(size_t mj = 0 ; mj < Num_MemberJets() ; ++mj)
	{
		for(size_t mjb = 0 ; mjb < b.Num_MemberJets() ; ++mjb)
		{
			if(MemberJets(mj) == b.MemberJets(mjb)) {return true;}
		}
	}
	return false;
}
bool  HadTopJet::CloseToMember(const TLorentzVector* b, double dr) const
{
	for(const TLorentzVector* mj : jets_)
	{
		if(b->DeltaR(*mj) < dr) {return true;}
	}
	return false;
}


double HadTopJet::ScaleRes(double chadtopscale_, double chadtopres_, HadTopJet* genhtjet)
{
	genhtjet_ = genhtjet;
	rescale_ = chadtopscale_;
	if(genhtjet != nullptr)
	{
		rescale_ = (genhtjet->E() + chadtopres_*(chadtopscale_*E() - genhtjet->E()))/E();
	}
	SetPxPyPzE(rescale_*Px(), rescale_*Py(), rescale_*Pz(), rescale_*E());		
	JetSum_*=rescale_;
	Init();
	return rescale_;
}

bool HadTopJet::PS_MODE = false;

