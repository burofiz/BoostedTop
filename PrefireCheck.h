#ifndef PREFIRECHE_H
#define PREFIRECHE_H
#include <unordered_map>
#include <unordered_set>

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>

using namespace std;


class PrefireCheck
{
	private:
		unordered_map<Long64_t, unordered_set<Long64_t> > evmap_;

		Long64_t ID(Long64_t run, Long64_t lumi)
		{
			return lumi*10E7+run;
		}
	public:

		void Init(const string& EventFile)
		{
			TDirectory* cdir = gDirectory;
			TFile* evfile = TFile::Open(EventFile.c_str());
			TTree* evtr = dynamic_cast<TTree*>(evfile->Get("tree"));

			Long64_t run;
			Long64_t lumi;
			Long64_t event;

			evtr->SetBranchAddress("run", &run);
			evtr->SetBranchAddress("lumi", &lumi);
			evtr->SetBranchAddress("event", &event);


			for(int n = 0 ; n < evtr->GetEntries() ; ++n)
			{
				evtr->GetEntry(n);
				Long64_t id = ID(run, lumi);
				evmap_[id].insert(event);
			}

			cdir->cd();
		}

		bool Unprefired(Long64_t run, Long64_t lumi, Long64_t event)
		{
			Long64_t id = ID(run, lumi);
			auto checkid = evmap_.find(id);
			if(checkid == evmap_.end()) {return false;}

			return checkid->second.find(event) != checkid->second.end();

		}

};

#endif
