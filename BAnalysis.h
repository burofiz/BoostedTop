#ifndef H_BANALYSIS
#define H_BANALYSIS
#include <helper.h>
#include <Analyse.h>
#include <OMuon.h>
#include <GenSelectedParticle.h>

class TTEvent;
class BoostedTop;

class BAnalysis
{
	private:
		BoostedTop* AN = nullptr;
		TH1DCollection hists1d_;
		TH2DCollection hists2d_;

		ocontainer<OMuon> SMuons;
		vector<OMuon*> RecMusAllPt;
		ocontainer<GenSelectedParticle> SGenParticles;
		vector<GenSelectedParticle*> GenMusBTop;
		vector<GenSelectedParticle*> GenMusBnoTop;
		vector<GenSelectedParticle*> GenBsTop;


		bool init_ = false;
	public:

		void Init();
		BAnalysis();


		void RunEvent(const TTEvent& ttrecoev);//, const TTEvent& ttgenev);


};











#endif


