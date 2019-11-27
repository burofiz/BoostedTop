#ifndef DELPHESJEC_H
#define DELPHESJEC_H

#include <TFile.h>
#include <TH1D.h>
#include <TDirectory.h>


#include <OJet.h>

#include <iostream>
#include <vector>

using namespace std;


class DelphesJEC
{
private:
	TFile* tf_ = nullptr;
	TH1D* hbj_ = nullptr;
	TH1D* hlj_ = nullptr;


public:

	void Init(const string& filename)
	{
		TDirectory* cdir = gDirectory;
		tf_ = TFile::Open(filename.c_str());
		hbj_ = dynamic_cast<TH1D*>(tf_->Get("TOTAL_BJES_AntiKt4EMTopo_YR2018"));
		hlj_ = dynamic_cast<TH1D*>(tf_->Get("TOTAL_DIJET_AntiKt4EMTopo_YR2018"));

		cdir->cd();
	}

template <typename T>
	double GetScale(const OJet* jet, const vector<T*> genbjets, double sigmascale = 0.)
	{
		if(sigmascale == 0.) {return 1.;}
		bool BJET = (find_if(genbjets.begin(), genbjets.end(), [&](TLorentzVector* bp){return jet->DeltaR(*bp) < 0.3;}) != genbjets.end());
		double scale = sigmascale;
		if(BJET)
		{
			scale *= hbj_->Interpolate(min(jet->Pt(), 3000.));
		}
		else
		{
			scale *= hlj_->Interpolate(min(jet->Pt(), 3000.));
		}
		
		return scale+1.;	

	}


};






#endif
