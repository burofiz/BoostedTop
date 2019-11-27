#ifndef TTPLOTS_H
#define TTPLOTS_H

#include <helper.h>

class TTEvent;
class BoostedTop;

class Plots
{
	public:
		string _prefix;
		TH1DCollection _hists1d;
		TH2DCollection _hists2d;
		Plots(const string& prefix);

		virtual void Init() = 0;
};

class TTPlots : public Plots
{
	private:
		vector<string> leptypes_ = {"mu", "el"};
	public:
		TTPlots(const string& prefix);

		virtual void Init();
		void Fill(const string& leptype, TTEvent& ev, double weight);
		void Fill(TTEvent& ev, double weight);
};

class HadTopJet;

class HadJetHists : public Plots
{
	private:
	public:
		HadJetHists(const string& prefix);
		virtual void Init();
		void Fill(HadTopJet* tj, double weight);
};

class LepTopJet;

class LepJetHists : public Plots
{
	private:
	public:
		LepJetHists(const string& prefix);
		virtual void Init();
		void Fill(LepTopJet* tj, double weight);
};

#endif
