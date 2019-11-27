#ifndef NNJET_H
#define NNJET_H

#include <helper.h>
#include <NN.h>
#include <OJet.h>

#include <map>
#include <iostream>

using namespace std;

//class NNJet
//{
//	private:
//		int insize_ = 14;
//		map<Bin, NeuralNet> NNs_;
//		Layer layin_;
//
//
//	public:
//		NNJet(const vector<double>& etabins, const vector<string>& NNfiles) : layin_(insize_)
//	{
//		for(size_t b = 0 ; b < etabins.size()-1 ; ++b)
//		{
//			Bin eta(etabins[b], etabins[b+1]);
//			NNs_[eta] = NeuralNet(NNfiles[b]);	
//		}
//	}
//
//
//		TLorentzVector getscale(OJet& jet, double rho, double vertices)	
//		{
//			if(abs(jet.Eta()) >= 2.5) {return TLorentzVector(0.,0.,0.,0.);}
//
//			double s = 1000.;
//			layin_.output()(0, 0) = jet.ChargedHadronEnergy()/s;
//			layin_.output()(1, 0) = jet.ChargedHadronEnergy()/max(1,jet.NumChargedHadrons())/s*50;
//			layin_.output()(2, 0) = (float)jet.NumChargedHadrons()/50;
//			layin_.output()(3, 0) = jet.NeutralHadronEnergy()/s;
//			layin_.output()(4, 0) = jet.NeutralHadronEnergy()/max(1,jet.NumNeutralHadrons())/s*50;
//			layin_.output()(5, 0) = (float)jet.NumNeutralHadrons()/50;
//			layin_.output()(6, 0) = jet.PhotonEnergy()/s;
//			layin_.output()(7, 0) = jet.PhotonEnergy()/max(1,jet.NumPhotons())/s*50;
//			layin_.output()(8, 0) = (float)jet.NumPhotons()/50;
//			layin_.output()(9, 0) = jet.MuonEnergy()/s;
//			layin_.output()(10, 0) = jet.ElectronEnergy()/s;
//			layin_.output()(11, 0) = rho/100;
//			layin_.output()(12, 0) = vertices/100.;
//			layin_.output()(13, 0) = jet.Area();
//
//			TMatrixD nnoutcut(NNs_[jet.Eta()].apply(layin_));
//			double r = nnoutcut(0,0)+1.;
////cout << "U: " << jet.RawE() << " " << r <<  " " << jet.EnergyCorrection() << " " << jet.E() << endl;
//			double scale = r*jet.EnergyCorrection();
//			return jet.ApplySF(scale);
//		}
//
//
//
//};


//class JetPlots
//{
//	private:
//		map<Bin, TH2D*> Hres_;
//		map<Bin, TH2D*> Hscale_;
//	public:
//		void Init(const string& prefix)
//		{
//			vector<double> etabins = {-2.5, -2.0, -1.5, -1.0, -0.5, 0., 0.5, 1.0, 1.5, 2.0, 2.5};
//
//			for(size_t b = 0 ; b < etabins.size()-1 ; ++b)
//			{
//				Bin eta(etabins[b], etabins[b+1]);
//				string name = prefix + "_res_" + to_string(b);
//				Hres_[eta] = new TH2D(name.c_str(), name.c_str(), 500, 0, 1000, 100, -1., 1.);	
//			}
//			for(size_t b = 0 ; b < etabins.size()-1 ; ++b)
//			{
//				Bin eta(etabins[b], etabins[b+1]);
//				string name = prefix + "_scale_" + to_string(b);
//				Hscale_[eta] = new TH2D(name.c_str(), name.c_str(), 500, 0, 1000, 100, -1., 1.);	
//			}
//
//		}
//
//		void Fill(OJet* reco, TLorentzVector* gen, double weight)
//		{
//			if(gen != nullptr){Hres_[reco->Eta()]->Fill(gen->E(), (reco->E()-gen->E())/gen->E(), weight);}
//			Hscale_[reco->Eta()]->Fill(reco->E(), reco->GetSF()-1., weight);
////			cout << "P: " << reco->RawE() << " " << reco->E() << " " << reco->GetSF() << endl;
//		}
//
//};







#endif

