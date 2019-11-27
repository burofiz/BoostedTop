#include "BAnalysis.h"
#include "BoostedTop.h"


double PtRel(const TLorentzVector* mu, const TLorentzVector* jet)
{
	return mu->Perp(jet->Vect())/mu->P();
}

BAnalysis::BAnalysis() : 
	AN(dynamic_cast<BoostedTop*>(GLAN)),
	hists1d_("b1d"),
	hists2d_("b2d")

{}
double_t dist(const TVector2& v0, const TVector2& v1){
	return std::sqrt(std::pow(v0.Px()-v1.Px(),2) + std::pow(v0.DeltaPhi(v1),2));
}

double_t leng(const TVector2& v0){
	return std::sqrt(std::pow(v0.Px(),2) + std::pow(v0.Py(),2));
}
double_t AngleDiff(const TVector2& pv, const TVector2& dv){//angle between two TVector2's, always between 0 and pi
	double_t dProduct = pv.Px() * dv.Px() + pv.Py() * dv.Py();
	double_t angle = std::acos(dProduct/(leng(pv)*leng(dv)));
	return angle;
}
/*double_t SignedAnDiff(const TVector2& pv, const TVector2& dv){//angle between two TVector2's, signed, from -pi to pi
															  //determinant/dotproduct = tangent
	double_t dProduct = pv.Px() * dv.Px() + pv.Py() * dv.Py();
	double_t determinant = pv.Px() * dv.Px() - pv.Py() * dv.Py();
	return TMath::ATan2(dProduct, determinant);
}*/
void VTranslate(TVector2* v, const TVector2& t){//translates vector v by vector t
	v->Set(v->Px()+t.Px(), v->Py()+t.Py());
}
bool PhiFilter(double_t phi){//returns true iff pi-0.8 < phi < -pi+0.8
	if((Pi()-0.8>phi)&&(phi>0.8-Pi())){
		return true;
	}else{
		return false;
	}
}

double_t SignedAnDiff2(const TVector2& v0, const TVector2& v1){//Signed angle from v0 to v1 (meaning: rotate v0 by angle and get v1) 
	double_t angle = atan2(v1.Py(), v1.Px()) - atan2(v0.Py(),v0.Px());
	//if (angle < 0) { angle += 2 * TMath::Pi(); }
	if (angle > Pi())        { angle -= 2 * Pi(); }
	else if (angle <= -Pi()) { angle += 2 * Pi(); }
	return angle;
}
void BAnalysis::Init()
{
	if(init_) {return;} init_ = true;

	hists1d_.AddHist("muBT_pt", 200, 0., 200);
	hists1d_.AddHist("muBT_ptrel", 200, 0., 1);
	hists1d_.AddHist("muBKG_pt", 200, 0., 200);
	hists1d_.AddHist("muBKG_ptrel", 200, 0., 1);
	
	hists1d_.AddHist("BBjetb0", 100, -Pi(), Pi());
	hists1d_.AddHist("BBpull", 100, -Pi(), Pi());//Pull vectors angle with difference vector and each other
	hists1d_.AddHist("BBpullb0", 100, -Pi(), Pi());//Pull vectors angle with difference vector and each other
	hists1d_.AddHist("BBpullb1", 100, -Pi(), Pi());
	hists1d_.AddHist("BBPullbs", 100, 0, Pi());
	hists1d_.AddHist("BBMass", 100, 300, 1000);
	hists1d_.AddHist("BBbeamb0", 100, -Pi(), Pi());//Beam and jet
	hists1d_.AddHist("BBbeamb1", 100, -Pi(), Pi());
	hists2d_.AddHist("BBPullM1", 100, 300., 1000, 100, 0, Pi());//Pull vs Mass
	hists2d_.AddHist("BBPullM2", 100, 300., 1000, 100, 0, Pi());
	hists2d_.AddHist("BBpull2D", 100, 0, Pi(), 100, 0, Pi());
	hists2d_.AddHist("BBEtaPhi3", 100, -0.1, 0.1, 100, -0.1, 0.1); //Eta-Phi plots for the rotated systems
	hists2d_.AddHist("BBEtaPhi4", 100, -0.1, 0.1, 100, -0.1, 0.1);
	
	hists1d_.AddHist("BBpullbv", 100, -Pi(), Pi());//For angle between pullb0 and v0
	
	hists1d_.AddHist("WWpull", 100, -Pi(), Pi());//Same just for W jets 
	hists1d_.AddHist("WWjetv1", 100, -Pi(), Pi());
	hists1d_.AddHist("WWpullv0", 100, 0, Pi());//Same just for W jets
	hists1d_.AddHist("WWpullv1", 100, 0, Pi());
	hists2d_.AddHist("WWpullD", 100, -Pi(), Pi(), 100, -Pi(), Pi());
	hists1d_.AddHist("WWPullPull", 100, 0, 2*Pi());
	hists2d_.AddHist("WWEtaPhi3", 100, -0.1, 0.1, 100, -0.1, 0.1); //Eta-Phi plots for the rotated systems
	hists2d_.AddHist("WWEtaPhi4", 100, -0.1, 0.1, 100, -0.1, 0.1);
	hists1d_.AddHist("WWbeamv0", 100, 0, 2*Pi());//Pull vector and beam
	hists1d_.AddHist("WWbeamv1", 100, 0, 2*Pi());
	
	int bins[] = {0, 340, 400, 430, 460, 500, 600, 700, 800};
	int i=1;
	for(i=1; i<9; i++){//Store W-plots for different masses
		hists1d_.AddHist("WWbeamv0_"+std::to_string(bins[i]), 50, -Pi(), Pi());// simply creates histograms whose names contain the values in "bins at the end"
		hists1d_.AddHist("WWbeamv1_"+std::to_string(bins[i]), 50, -Pi(), Pi());
		hists1d_.AddHist("WWPullPull_"+std::to_string(bins[i]), 50, -Pi(), Pi());
		hists1d_.AddHist("WWpullv0_"+std::to_string(bins[i]), 50, -Pi(), Pi());
		hists1d_.AddHist("WWpullv1_"+std::to_string(bins[i]), 50, -Pi(), Pi());
		hists2d_.AddHist("WWpullD_"+std::to_string(bins[i]), 100, -Pi(), Pi(), 100, -Pi(), Pi());
		hists2d_.AddHist("WWEtaPhi3_"+std::to_string(bins[i]), 100, -0.1, 0.1, 100, -0.1, 0.1); 
		hists2d_.AddHist("WWEtaPhi4_"+std::to_string(bins[i]), 100, -0.1, 0.1, 100, -0.1, 0.1);
	}
	i=1;
	for(i=1; i<9; i++){//Generate B-plots for different masses
		hists1d_.AddHist("BBpullb0_"+std::to_string(bins[i]), 50, 0, 2*Pi());// simply creates histograms whose names contain the values in "bins at the end"
		hists1d_.AddHist("BBpullb1_"+std::to_string(bins[i]), 50, 0, 2*Pi());
		hists1d_.AddHist("BBPullbs_"+std::to_string(bins[i]), 50, -Pi(), 2*Pi());
		hists2d_.AddHist("BBpull2D_"+std::to_string(bins[i]), 50, 0, Pi(), 50, 0, Pi());
		hists1d_.AddHist("BBbeamb0_"+std::to_string(bins[i]), 50, -Pi(), Pi());//Beam and jet
		hists1d_.AddHist("BBbeamb1_"+std::to_string(bins[i]), 50, -Pi(), Pi());
		hists1d_.AddHist("BBpullbv_"+std::to_string(bins[i]), 50, -Pi(), Pi());
		hists2d_.AddHist("BBEtaPhi3_"+std::to_string(bins[i]), 50, -0.1, 0.1, 50, -0.1, 0.1); //Eta-Phi plots for the rotated systems
		hists2d_.AddHist("BBEtaPhi4_"+std::to_string(bins[i]), 50, -0.1, 0.1, 50, -0.1, 0.1);
	}
	
	/*
	 * Same as above, but with the data from the generated events
	 */
	//B jets
	/*
	hists1d_.AddHist("BBpullb0G", 200, 0, 2*Pi());//Pull vectors angle with difference vector and each other
	hists1d_.AddHist("BBpullb1G", 200, 0, 2*Pi());
	hists1d_.AddHist("BBPullbsG", 200, -Pi(), 2*Pi());
	hists1d_.AddHist("BBbeamb0G", 200, -Pi(), Pi());//Beam and jet
	hists1d_.AddHist("BBbeamb1G", 200, -Pi(), Pi());
	hists1d_.AddHist("BBpullbvG", 200, -Pi(), Pi());//For angle between pullb0 and v0
	hists2d_.AddHist("BBpull2DG", 200, 0, Pi(), 200, 0, Pi());
	hists2d_.AddHist("BBEtaPhi3G", 200, -0.1, 0.1, 200, -0.1, 0.1); //Eta-Phi plots for the rotated systems
	hists2d_.AddHist("BBEtaPhi4G", 200, -0.1, 0.1, 200, -0.1, 0.1);
	
	
	
	//W jets
	hists1d_.AddHist("WWpullv0G", 200, -Pi(), Pi());
	hists1d_.AddHist("WWpullv1G", 200, -Pi(), Pi());
	hists2d_.AddHist("WWpullDG", 200, -Pi(), Pi(), 200, -Pi(), Pi());
	hists1d_.AddHist("WWPullPullG", 200, -Pi(), Pi());
	hists1d_.AddHist("WWbeamv0G", 200, -Pi(), Pi());
	hists1d_.AddHist("WWbeamv1G", 200, -Pi(), Pi());
	hists2d_.AddHist("WWEtaPhi3G", 200, -0.1, 0.1, 200, -0.1, 0.1); //Eta-Phi plots for the rotated systems
	hists2d_.AddHist("WWEtaPhi4G", 200, -0.1, 0.1, 200, -0.1, 0.1);
	*/
	/*
	 * End of generated part
	 */
}

void BAnalysis::RunEvent(const TTEvent& ttrecoev)//{//, const TTEvent& ttgenev)
{
	SGenParticles.reset(AN->NumSelectedGenParticles());
	GenMusBTop.clear();
	GenMusBnoTop.clear();
	GenBsTop.clear();	
	/*
		-----------------------------------------------------------------------------
		Generated W jets
	*/
	/*
	
	if(ttgenev.AddJets().size()==0){//If there are no extraneous jets
		TVector2 diffVG(v0Gen.Px() - v1Gen.Px(),  v0Gen.Py() - v1Gen.Py());//Difference vector between the two jets
		OJet* wjmaxG = dynamic_cast<OJet*>(ttgenev.WJPtmax());//Casts to jet type
		OJet* wjminG = dynamic_cast<OJet*>(ttgenev.WJPtmin());
		//double_t asdfa = wjmaxG->Pull(1);
		TVector2 pullV0G(wjmaxG->Pull(0), wjmaxG->Pull(1)); //creates pull vectors
		TVector2 pullV1G(wjminG->Pull(0), wjminG->Pull(1)); 
		double_t angleV0G = SignedAnDiff2(pullV0G, diffVG);//computes angle between pull and difference vectors
		double_t angleV1G = SignedAnDiff2(pullV1G, diffVG);
		double_t anglePPG = AngleDiff(pullV0G, pullV1G);
		wjmaxG->Delete(); wjminG->Delete();
		if(PhiFilter(pullV0G.Py())&&PhiFilter(pullV1G.Py())){
			hists1d_["WWbeamv0G"]->Fill(TMath::ACos(pullV0G.Px()/TMath::Sqrt(pullV0G.Px()*pullV0G.Px()+pullV0G.Py()*pullV0G.Py())), AN->weight);
			hists1d_["WWbeamv1G"]->Fill(TMath::ACos(pullV1G.Px()/TMath::Sqrt(pullV1G.Px()*pullV1G.Px()+pullV1G.Py()*pullV1G.Py())), AN->weight);
			hists1d_["WWPullPullG"]->Fill(anglePPG, AN->weight);
			hists1d_["WWpullv0G"]->Fill(TMath::Abs(angleV0G), AN->weight);
			hists1d_["WWpullv1G"]->Fill(TMath::Abs(angleV1G), AN->weight);	
			hists2d_["WWpullDG"]->Fill(angleV0G, angleV1G, AN->weight);
			
			double_t angleVVG = std::atan2(diffVG.Py(), diffVG.Px());
			pullV0G = pullV0G.Rotate(-1*angleVVG);
			hists2d_["WWEtaPhi3G"]->Fill(pullV0G.Px(), pullV0G.Py(), AN->weight);
			pullV0G = pullV0G.Rotate(2*angleVVG);
			hists2d_["WWEtaPhi4G"]->Fill(pullV0G.Px(), pullV0G.Py(), AN->weight);
		}
	}*/
	
	/*
	-------------------------------------------------------------------------------
	//	Reconstructed W jet   //
	//						  //
	*/
	TVector2 v0(ttrecoev.WJPtmax()->Eta(), ttrecoev.WJPtmax()->Phi());
	TVector2 v1(ttrecoev.WJPtmin()->Eta(), ttrecoev.WJPtmin()->Phi());
	//std::cout<<"BA"<<endl;
	//TVector2 v0Gen(ttgenev.WJPtmax()->Eta(), ttgenev.WJPtmax()->Phi());
	//TVector2 v1Gen(ttgenev.WJPtmin()->Eta(), ttgenev.WJPtmin()->Phi());
	
	//if((ttrecoev.AddJets().size()==ttgenev.AddJets().size()) &&(ttrecoev.AddJets().size()==0)){// no extra jets
	
	if(ttrecoev.AddJets().size()==0){
		std::cout<<"aasd"<<endl;
		TVector2 diffV(v0.Px() - v1.Px(),  v0.Py() - v1.Py());//Difference vector between the two jets
	///	std::cout << "V0.x = " << v0.Px() << ", V0.y = " << v0.Py() << endl;
		//std::cout << "V1.x = " << v1.Px() << ", V1.y = " << v1.Py() << endl;
		//std::cout << "V01.x = " << v0.Px()-v1.Px() << ", V01.y = " << v0.Py()-v1.Py() << endl;
		OJet* wjmax = dynamic_cast<OJet*>(ttrecoev.WJPtmax());//Casts to jet type
		OJet* wjmin = dynamic_cast<OJet*>(ttrecoev.WJPtmin());
		TVector2 pullV0(wjmax->Pull(0), wjmax->Pull(1)); //creates pull vectors
		TVector2 pullV1(wjmin->Pull(0), wjmin->Pull(1)); 
		if(PhiFilter(pullV0.Py())&&PhiFilter(pullV1.Py())){
			double_t angleV0 = SignedAnDiff2(pullV0, diffV);//computes angle between pull and difference vectors
			double_t angleV1 = SignedAnDiff2(pullV1, diffV);
			double_t anglePP = AngleDiff(pullV0, pullV1);
			wjmax->Delete(); wjmin->Delete();
			
			hists1d_["WWpull"]->Fill(std::atan2(pullV0.Py(), pullV0.Px()));
			hists1d_["WWjetv1"]->Fill(TMath::ATan2(diffV.Py(), diffV.Px()));
			hists1d_["WWbeamv0"]->Fill(TMath::ACos(pullV0.Px()/TMath::Sqrt(pullV0.Px()*pullV0.Px()+pullV0.Py()*pullV0.Py())), AN->weight);
			hists1d_["WWbeamv1"]->Fill(TMath::ACos(pullV1.Px()/TMath::Sqrt(pullV1.Px()*pullV1.Px()+pullV1.Py()*pullV1.Py())), AN->weight);
			hists1d_["WWPullPull"]->Fill(anglePP, AN->weight);
			hists1d_["WWpullv0"]->Fill(TMath::Abs(angleV0), AN->weight);
			hists1d_["WWpullv1"]->Fill(TMath::Abs(angleV1), AN->weight);	
			hists2d_["WWpullD"]->Fill(angleV0, angleV1, AN->weight);
			double_t angleVV = std::atan2(diffV.Py(), diffV.Px());
			//std::cout << "diffV.x = " << diffV.Px() << ", diffV. y = " << diffV.Py() << ", Angle = " << TMath::ATan2(diffV.Py(), diffV.Px()) << endl;
		//	std::cout<<"---------------------"<<endl;
			pullV0 = pullV0.Rotate(-1*angleVV);
			pullV1 = pullV1.Rotate(-1*angleVV);
			diffV = diffV.Rotate(-1*angleVV);
			hists2d_["WWEtaPhi3"]->Fill(pullV0.Px(), pullV0.Py(), AN->weight);
			hists2d_["WWEtaPhi4"]->Fill(pullV1.Px(), pullV1.Py(), AN->weight);

			//This part is same as above, but bins by mass
			double_t mass = ttrecoev.TT()->M();
			int bins[] = {0, 340, 400, 430, 460, 500, 600, 700, 800};
			int i=1;
			for(i=1; i<8 && mass>bins[i]; i++){} //returns the index of the right bin.  Skips 0, because there are no negative masses.
			hists1d_["WWbeamv0_"+std::to_string(bins[i])]->Fill(TMath::ACos(pullV0.Px()/TMath::Sqrt(pullV0.Px()*pullV0.Px()+pullV0.Py()*pullV0.Py())), AN->weight);
			hists1d_["WWbeamv1_"+std::to_string(bins[i])]->Fill(TMath::ACos(pullV1.Px()/TMath::Sqrt(pullV1.Px()*pullV1.Px()+pullV1.Py()*pullV1.Py())), AN->weight);
			hists1d_["WWPullPull_"+std::to_string(bins[i])]->Fill(anglePP, AN->weight);
			hists1d_["WWpullv0_"+std::to_string(bins[i])]->Fill(TMath::Abs(angleV0), AN->weight);
			hists1d_["WWpullv1_"+std::to_string(bins[i])]->Fill(TMath::Abs(angleV1), AN->weight);	
			hists2d_["WWpullD_"+std::to_string(bins[i])]->Fill(angleV0, angleV1, AN->weight);
			hists2d_["WWEtaPhi3_"+std::to_string(bins[i])]->Fill(pullV0.Px(), pullV0.Py(), AN->weight);
			hists2d_["WWEtaPhi4_"+std::to_string(bins[i])]->Fill(pullV1.Px(), pullV1.Py(), AN->weight);
		}

	}
	/*
		-----------------------------------------------------------------------------
		B-Jet: Generated
	*/
	
	TLorentzVector*  BHad = ttrecoev.BHad();
	TLorentzVector*  BLep = ttrecoev.BLep();
	
	TVector2 b0(BHad->Eta(), BHad->Phi());
	TVector2 b1(BLep->Eta(), BLep->Phi());
	//TVector2 b0Gen(ttgenev.BHad()->Eta(), ttgenev.BHad()->Phi());//Generated values, considered "true";
	//TVector2 b1Gen(ttgenev.BLep()->Eta(), ttgenev.BLep()->Phi());
	
	
	/*if(ttgenev.AddJets().size()==0){
		TVector2 diffBG(b0Gen.Px() - b1Gen.Px(),  b0Gen.Py() - b1Gen.Py());//Difference vector between the two jets
		OJet* bhadG = dynamic_cast<OJet*>(BHad);//Casts to jet type
		OJet* blepG = dynamic_cast<OJet*>(BLep);
		TVector2 pullB0G(bhadG->Pull(0), bhadG->Pull(1)); //creates pull vectors
		TVector2 pullB1G(blepG->Pull(0), blepG->Pull(1)); 
		double_t angleB0G = SignedAnDiff2(pullB0G, diffBG);//computes angle between pull and difference vectors
		double_t angleB1G = SignedAnDiff2(pullB1G, diffBG);
		double_t anglePPBG = AngleDiff(pullB0G, pullB1G);
		bhadG->Delete(); blepG->Delete();
			
		if(PhiFilter(pullB0G.Py())&&PhiFilter(pullB1G.Py())){
			hists1d_["BBMassG"]->Fill(ttgenev.TT()->M(), AN->weight);
			hists1d_["BBpullb0G"]->Fill(angleB0G, AN->weight);
			hists1d_["BBpullb1G"]->Fill(angleB1G, AN->weight);
			hists1d_["BBPullbsG"]->Fill(anglePPBG, AN->weight);	
			hists2d_["BBpull2DG"]->Fill(angleB0G, angleB1G, AN->weight);
			hists1d_["BBbeamb0G"]->Fill(TMath::ACos(pullB0G.Px()/TMath::Sqrt(pullB0G.Px()*pullB0G.Px()+pullB0G.Py()*pullB0G.Py())), AN->weight);
			hists1d_["BBbeamb1G"]->Fill(TMath::ACos(pullB1G.Px()/TMath::Sqrt(pullB1G.Px()*pullB1G.Px()+pullB1G.Py()*pullB1G.Py())), AN->weight);
				
			//Do the same but with two random vectors. No correlation expected
			TVector2 v0G(ttgenev.WJPtmax()->Eta(), ttgenev.WJPtmax()->Phi());
			TVector2 diffBVG(b0Gen.Px()-v0G.Px(),  b0Gen.Py()-v0G.Py());
			double_t angleBV0G = AngleDiff(pullB0G, diffBVG);
			hists1d_["BBpullbvG"]->Fill(angleBV0G, AN->weight);
	
			double_t angleG = std::atan2(diffBG.Py(), diffBG.Px());
			pullB0G = pullB0G.Rotate(-1*angleG);
			//pullB1 = pullB1.Rotate(-1*angleG);
			hists2d_["BBEtaPhi3G"]->Fill(pullB0G.Px(), pullB0G.Py(), AN->weight);
			pullB0G = pullB0G.Rotate(2*angleG);
			hists2d_["BBEtaPhi4G"]->Fill(pullB1G.Px(), pullB1G.Py(), AN->weight);
		}
	}*/
	
	
	
	/*
		-----------------------------------------------------------------------------
		B-Jet: reconstructed
	*/
	if(ttrecoev.AddJets().size()==0){
		std::cout << "Accepted" << endl;
		TVector2 diffB(b0.Px() - b1.Px(),  b0.Py() - b1.Py());//Difference vector between the two jets
		OJet* bhad = dynamic_cast<OJet*>(BHad);//Casts to jet type
		OJet* blep = dynamic_cast<OJet*>(BLep);
		TVector2 pullB0(bhad->Pull(0), bhad->Pull(1)); //creates pull vectors
		TVector2 pullB1(blep->Pull(0), blep->Pull(1)); 
		double_t angleB0 = SignedAnDiff2(pullB0, diffB);//computes angle between pull and difference vectors
		double_t angleB1 = SignedAnDiff2(pullB1, diffB);
		double_t anglePPB = AngleDiff(pullB0, pullB1);
		bhad->Delete(); blep->Delete();
		
		if(PhiFilter(pullB0.Py())&&PhiFilter(pullB1.Py())){
			double_t mass = ttrecoev.TT()->M();
			hists1d_["BBjetb0"]->Fill(std::atan2(diffB.Py(), diffB.Px()));
			hists1d_["BBpull"]->Fill(std::atan2(pullB0.Py(), pullB0.Px()));
			hists1d_["BBMass"]->Fill(mass, AN->weight);
			hists2d_["BBPullM1"]->Fill(ttrecoev.TT()->M(), angleB0, AN->weight);	
			hists2d_["BBPullM2"]->Fill(ttrecoev.TT()->M(), angleB1, AN->weight);
			hists1d_["BBpullb0"]->Fill(TMath::Abs(angleB0), AN->weight);
			hists1d_["BBpullb1"]->Fill(TMath::Abs(angleB1), AN->weight);
			hists1d_["BBPullbs"]->Fill(anglePPB, AN->weight);	
			hists2d_["BBpull2D"]->Fill(angleB0, angleB1, AN->weight);
			hists1d_["BBbeamb0"]->Fill(TMath::ACos(pullB0.Px()/TMath::Sqrt(pullB0.Px()*pullB0.Px()+pullB0.Py()*pullB0.Py())), AN->weight);
			hists1d_["BBbeamb1"]->Fill(TMath::ACos(pullB1.Px()/TMath::Sqrt(pullB1.Px()*pullB1.Px()+pullB1.Py()*pullB1.Py())), AN->weight);
			
			//Do the same but with two random vectors. No correlation expected
			TVector2 v0(ttrecoev.WJPtmax()->Eta(), ttrecoev.WJPtmax()->Phi());
			TVector2 diffBV(b0.Px()-v0.Px(),  b0.Py()-v0.Py());
			double_t angleBV0 = AngleDiff(pullB0, diffBV);
			hists1d_["BBpullbv"]->Fill(angleBV0, AN->weight);
			
			double_t angle = std::atan2(diffB.Py(), diffB.Px());
			pullB0 = pullB0.Rotate(angle);
			pullB1 = pullB1.Rotate(angle);
			hists2d_["BBEtaPhi3"]->Fill(pullB0.Px(), pullB0.Py(), AN->weight);
			hists2d_["BBEtaPhi4"]->Fill(pullB1.Px(), pullB1.Py(), AN->weight);
			
			int bins[] = {0, 340, 400, 430, 460, 500, 600, 700, 800};
			int i=1;
			for(i=1; i<8 && mass>bins[i]; i++){} //returns the index of the right bin. For example mass=410 would return i=3;hists1d_["BBpullb0_"+std::to_string(bins[i])]->Fill(angleB0, AN->weight);
			hists1d_["BBpullb0_"+std::to_string(bins[i])]->Fill(angleB0, AN->weight);
			hists1d_["BBpullb1_"+std::to_string(bins[i])]->Fill(angleB1, AN->weight);
			hists1d_["BBPullbs_"+std::to_string(bins[i])]->Fill(anglePPB, AN->weight);	
			hists2d_["BBpull2D_"+std::to_string(bins[i])]->Fill(angleB0, angleB1, AN->weight);
			hists1d_["BBbeamb0_"+std::to_string(bins[i])]->Fill(TMath::ACos(pullB0.Px()/TMath::Sqrt(pullB0.Px()*pullB0.Px()+pullB0.Py()*pullB0.Py())), AN->weight);
			hists1d_["BBbeamb1_"+std::to_string(bins[i])]->Fill(TMath::ACos(pullB1.Px()/TMath::Sqrt(pullB1.Px()*pullB1.Px()+pullB1.Py()*pullB1.Py())), AN->weight);
			hists1d_["BBpullbv_"+std::to_string(bins[i])]->Fill(angleBV0, AN->weight);
			hists2d_["BBEtaPhi3_"+std::to_string(bins[i])]->Fill(pullB0.Px(), pullB0.Py(), AN->weight);
			hists2d_["BBEtaPhi4_"+std::to_string(bins[i])]->Fill(pullB1.Px(), pullB1.Py(), AN->weight);
			
		}
	}
	BHad->Delete(); BLep->Delete();
	
	/*
		------------------------------------------------------------------------------------
	*/
	
	SMuons.reset(AN->NumIOMuons());
	RecMusAllPt.clear();

	for(UInt_t i = 0 ; i < AN->NumSelectedGenParticles() ; i++)
	{   
		GenSelectedParticle gpar(AN->GetSelectedGenParticle(i));
		SGenParticles.push_back(gpar);
		if(abs(gpar.PDGID()) == 13 && abs(gpar.Eta()) < 2.4)
		{ 
			if(gpar.FromB() && gpar.FromTop())
			{
				GenMusBTop.push_back(&SGenParticles.back());
			}
			else if(gpar.FromB())
			{
				GenMusBnoTop.push_back(&SGenParticles.back());
			}
		}
		else if(abs(gpar.PDGID()) != 13 && gpar.FromTop() && !gpar.FromW())
		{
			GenBsTop.push_back(&SGenParticles.back());

		}
	}
//cout << "topbs: " << GenBsTop.size() << endl;

	for(UInt_t i = 0 ; i < AN->NumIOMuons() ; i++)
	{   
		OMuon mu(AN->GetIOMuon(i));
		if(mu.Pt() > 1. && abs(mu.Eta()) < 2.4)
		{ 
			SMuons.push_back(mu);
			if(mu.ID(OMuon::TIGHT, false) && mu.DeltaR(*ttrecoev.LLep()) > 0.0001 && (mu.DeltaR(*ttrecoev.BHad()) < 0.6 || mu.DeltaR(*ttrecoev.BLep()) < 0.6))
			{
				RecMusAllPt.push_back(&SMuons.back());
			}
		}
	}

	for(const OMuon* mu : RecMusAllPt)
	{
		const TLorentzVector* jetclose = mu->DeltaR(*ttrecoev.BHad()) < mu->DeltaR(*ttrecoev.BLep()) ? ttrecoev.BHad() : ttrecoev.BLep();
		auto itgmu =find_if(GenMusBTop.begin(), GenMusBTop.end(), [&](const TLorentzVector* gmu){return mu->DeltaR(*gmu) < 0.1;});
		if(itgmu != GenMusBTop.end())
		{
			hists1d_["muBT_pt"]->Fill(mu->Pt(), AN->weight);
			hists1d_["muBT_ptrel"]->Fill(PtRel(mu, jetclose), AN->weight);
		}
		else
		{
			hists1d_["muBKG_pt"]->Fill(mu->Pt(), AN->weight);
			hists1d_["muBKG_ptrel"]->Fill(PtRel(mu, jetclose), AN->weight);
		}
	}

}
		
		/*
		std::cout << "PullV0.x = " << pullV0.Px() << ", PullV0.y = " << pullV0.Py() <<endl;
		std::cout << "diffV.x = " << diffV.Px() << ", diffV.y = " << diffV.Py() <<endl;
		std::cout << "AngleDx = " << angleVV << endl;
		std::cout << "AnglePx = " << std::atan2(pullV0.Py(), pullV0.Px()) << endl;
		std::cout << "diffV.x = " << diffV.Px() << ", diffV.y = " << diffV.Py() <<endl;
		std::cout << "PullV0.x = " << pullV0.Px() << ", PullV0.y = " << pullV0.Py() <<endl;
		std::cout << "AngleDx = " << angleVV << endl;
		std::cout << "AnglePx = " << std::atan2(pullV0.Py(), pullV0.Px()) << endl;
		std::cout << "--------" << endl;
		*/
		
		/*
		//VTranslate(&v1, -1*v0);
		//VTranslate(&pullV1, -1*v0);
		//VTranslate(&pullV0, -1*v0);
		//VTranslate(&b0, -1*b0);
		//std::cout << "b0.x " << b0.Px() << " b0.y = " << b0.Py() <<endl;
		//v0.Set(0., 0.);
		//b1 = b1.Rotate(angle);
		//std::cout << "PullB0.x " << pullB1.Px() << " PullB0.y = " << pullB1.Py() <<endl;
		//std::cout << "Angle = " << angle << endl;
		//std::cout << "PullV0.y" << pullV0.Px() << " PullV0.y " << pullV0.Py() << "End" << endl;
		//std::cout << "PullB0.y" << pullB1.Px() << " PullB0.y " << pullB1.Py() <<endl;
		double_t angleV01 = SignedAnDiff2(pullV0, diffV);
		double_t angleV11 = SignedAnDiff2(pullV1, diffV);
		double_t anglePP1 = SignedAnDiff2(pullV0, pullV1);
		
		std::cout << "angleV0 = " << angleV0 << ", angleV1 = " << angleV1 << ", anglePP = " << anglePP <<endl;
		std::cout << "angleV01 = " << angleV01 << ", angleV11 = " << angleV11 << ", anglePP1 = " << anglePP1 <<endl;
		std::cout << "diffV.x = " << diffV.Px() << ", diffV. y = " << diffV.Py() << ", Angle = " << TMath::ATan2(diffV.Py(), diffV.Py()) << endl;
		std::cout << "pullV0.x = " << pullV0.Px() << ", pullV0. y = " << pullV0.Py() << endl;
		std::cout << "pullV1.x = " << pullV1.Px() << ", pullV1. y = " << pullV1.Py() << endl;
		  
		std::cout << "Angle diff = " <<100*((TMath::Abs(angleV01)+TMath::Abs(angleV11))-Abs(anglePP1))/Abs(anglePP1) << endl;
		std::cout << "Angle diff = " <<100*((TMath::Abs(angleV0)+TMath::Abs(angleV1))-Abs(anglePP))/Abs(anglePP) << endl;
		std::cout << "--------" << endl;
		/ *std::cout << "cos = " << pullV0.Px()/TMath::Sqrt(pullV0.Px()*pullV0.Px()+pullV0.Py()*pullV0.Py()) << "Angle = " << TMath::ACos(pullV0.Px()/TMath::Sqrt(pullV0.Px()*pullV0.Px()+pullV0.Py()*pullV0.Py())) <<endl;
		std::cout << "cos2 = " << pullV0.Py()/TMath::Sqrt(pullV0.Px()*pullV0.Px()+pullV0.Py()*pullV0.Py()) << "Angle = " << TMath::ACos(pullV0.Py()/TMath::Sqrt(pullV0.Px()*pullV0.Px()+pullV0.Py()*pullV0.Py())) <<endl;
		std::cout << "--------" << endl;*/
		
