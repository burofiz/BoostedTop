#ifndef PREFIREPROB_H
#define PREFIREPROB_H
#include <TDirectory.h>
#include <TFile.h>
#include <TH2D.h>

#include <helper.h>
#include <OPhoton.h>
#include <OJet.h>
#include <RLObject.h>

using namespace std;

class L1Prefire : public RLObject
{
	private:
		map<Bin, double> ps_EG;
		map<Bin, double> ps_isoEG;
		map<Bin, double> ps_EGer;
		map<Bin, double> ps_isoEGer;
		vector< pair<int, int> > in_EG;
		vector< pair<int, int> > in_isoEG;
		vector< pair<int, int> > in_EGer;
		vector< pair<int, int>> in_isoEGer;

	public:
		void virtual  UpdateRun() override;
		void virtual  UpdateLB() override;

		double FireProb(const L1Object& l1eg);


};


class PrefireProb
{
	private:

		vector<double> ps_EG;
		vector<double> ps_isoEG;
		vector<double> ps_EGer;
		vector<double> ps_isoEGer;

		TFile* tf = nullptr;
		TH2D* prefiremap = nullptr;
		TH2D* prefiremapph = nullptr;
		TH2D* prefiremapjet = nullptr;

		ocontainer<OPhoton> Photons;
		ocontainer<OJet> Jets;

		void Select();

	public:

		void Init(const string& prefiremapfile, double sigma = 0.);

		~PrefireProb();

		double getProb() const;
		double FireProb(const L1Object& l1eg) const;

};

#endif
