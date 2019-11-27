#ifndef THRESIDUALCORRECTOR_H
#define THRESIDUALCORRECTOR_H
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TLorentzVector.h>

#include <iostream>

using namespace std;


class ResidualsPT
{
	private:
		TFile* tf_ = nullptr;
		TGraphAsymmErrors* dc_ = nullptr; 
		TGraphAsymmErrors* dm_ = nullptr; 
		TGraphAsymmErrors* dp_ = nullptr; 
		double ptmin_, ptmax_;
	public:

		void init(const string& filename, int beta)
		{
			TDirectory* cdir = gDirectory;

			tf_ = TFile::Open(filename.c_str());
			dc_ = dynamic_cast<TGraphAsymmErrors*>(tf_->Get(("scale_c_"+to_string(beta)).c_str()));
			dm_ = dynamic_cast<TGraphAsymmErrors*>(tf_->Get(("scale_m"+to_string(beta)).c_str()));
			dp_ = dynamic_cast<TGraphAsymmErrors*>(tf_->Get(("scale_p"+to_string(beta)).c_str()));
			double dummy;
			dc_->GetPoint(0, ptmin_, dummy);
			dc_->GetPoint(dc_->GetN()-1, ptmax_, dummy);

			cdir->cd();

		}

		double GetScale(TLorentzVector* jet, double sigma = 0.)
		{
			double pt = jet->Pt();
			pt = max(ptmin_, pt);
			pt = min(ptmax_, pt);

			double sc = dc_->Eval(pt);
			double se = 0.;
			if(sigma < 0)
			{
				se = dm_->Eval(pt);
			}
			else
			{
				se = dp_->Eval(pt);
			}

			return sc + abs(sigma)*(se - sc);
		}

		double GetError(TLorentzVector* jet, double sigma = 0.)
		{
			double pt = jet->Pt();
			pt = max(ptmin_, pt);
			pt = min(ptmax_, pt);

			double sc = dc_->Eval(pt);
			double se = 0.;
			if(sigma < 0)
			{
				se = dm_->Eval(pt);
			}
			else
			{
				se = dp_->Eval(pt);
			}

			return (sc+abs(sigma)*(se - sc))/sc;
		}

};

class THResidualCorrector
{
	private:
		vector<ResidualsPT> etabins;

		size_t GetEtaBin(TLorentzVector* jet)
		{
			if(abs(jet->Eta()) < -1.) {return 0;}
			else if(abs(jet->Eta()) < 0.) {return 1;}
			else if(abs(jet->Eta()) < 1.) {return 2;}
			else {return 3;}
		}
	public:
		THResidualCorrector() : etabins(4) {}

		void init(const string& filename)
		{
			for(int beta : {0,1,2,3})
			{
				etabins[beta].init(filename, beta);
			}
		}

		double GetScale(TLorentzVector* jet, double sigma = 0.) {return etabins[GetEtaBin(jet)].GetScale(jet, sigma);}

		double GetError(TLorentzVector* jet, double sigma = 0.) {return etabins[GetEtaBin(jet)].GetError(jet, sigma);}

};

#endif
