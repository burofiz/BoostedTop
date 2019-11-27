#ifndef LEPEFFCORRECTIO_H
#define LEPEFFCORRECTIO_H
#include <TFile.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <OElectron.h>
#include <OMuon.h>
#include <OJet.h>

#include <string>
#include <vector>
#include <unordered_map>


using namespace std;

class LepEffScales
{
	private:
		string type_;
		TFile* efffile = nullptr;
		TH1D* heta = nullptr;
		vector<TH1D*> effhists_mc;
		vector<TH1D*> effhists_da;
		vector<double> results;

	public:

		LepEffScales() : results(3) {}
		void init(const string& filename, const string& type);

		const vector<double>& getvals(TLorentzVector* vec);

};

class LepEffCorrection
{
	private:
		unordered_map<string, LepEffScales> corrections_;

	public:

		void init(const string& mufile, const string& elfile);
		void init2(const string& mufile, const string& elfile);

		double correctionEvent(vector<OElectron*> els, vector<OMuon*> mus, vector<OJet*> jets, double sigmael = 0., double sigmamu = 0.);

		double correctionEvent2(vector<OElectron*> els, vector<OElectron*> looseels, vector<OMuon*> mus, vector<OMuon*> loosemus, vector<OJet*> jets, double sigmael, double sigmamu);

};


#endif
