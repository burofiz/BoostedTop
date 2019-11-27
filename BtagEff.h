#ifndef BTAGEFF_H
#define BTAGEFF_H
#include <TTree.h>
#include <helper.h>

class TTEvent;
class OJet;

class BtagEff
{
	private:
		TH1DCollection _hists1d;
		TTree* btagtree;
		Float_t j[17];
		Float_t prob;
		Float_t problep;
		Float_t probhad;
		Float_t probnu;
		Float_t met;
		Float_t weight;
		Float_t nvtx;
		Int_t run;
		Int_t typ;
		Int_t test;
		double btagselection;
		double bkgcutmin_, bkgcutmax_;
		void PrepareFill(OJet* jet);
	public:
		BtagEff();
		void Init(double btagsel, double bkgcutmin, double bkgcutmax);
		bool Fill(TTEvent& per, TTEvent& genper,double theweight = 1.);
};

#endif
