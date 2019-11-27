#ifndef UNFOLDING_H
#define UNFOLDING_H

#include <helper.h>
#include <TTBarResponse1D.h>
#include <TTBarResponse2D.h>

#include <TH1D.h>

#include <PDFHists.h>

class TTEvent;
class BoostedTop;

class Unfolding
{
private:
	BoostedTop* AN = nullptr;

	TTBarResponse1D unfolding1d_;
	TTBarResponse2D unfolding2d_;
	TH2DCollection plot2d_;
	PDFHists<TH1D> pdf1d_;
	bool PDFUNCS_ = false;
public:
	Unfolding(const string& dir) : unfolding1d_(dir), unfolding2d_(dir) {}

	void Init(bool PDFUNCS = false);

	void FillTruth(const TTEvent& evt, double weight);
	void FillTruthReco(const TTEvent& truth_evt, const TTEvent& reco_evt, const string& category, double weight);
	void FillReco(const TTEvent& evt, const string& category, double weight);

	void FillRes(const TTEvent& truth_evt, const TTEvent& reco_evt, const string& category, double weight);

	float GetBin2D(const string& dist, double va, double vb) 
	{
		return unfolding2d_.GetBin(dist, va, vb)-0.5;
	}

};



#endif
