#include "hist.h"//define all the histograms
#include "PRINT4VECTOR.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
using namespace std;
using namespace erhic;

#define MASS_MUON  0.1056

double sPN_bins[]={0.,1.0,2.0,3.0,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.7,5.0,5.5,6.0,7.0,8.0,9.0,10.0,12.0,15.0};
int sPN_nBins = sizeof(sPN_bins)/sizeof(sPN_bins[0]) -1;

TH1D* nk = new TH1D("nk",";k (GeV) ;n(k) GeV^{-3}",200,0,1);

TLorentzRotation BoostToHCM(TLorentzVector const &eBeam_lab,
                            TLorentzVector const &pBeam_lab,
                            TLorentzVector const &eScat_lab) {
   TLorentzVector q_lab=eBeam_lab - eScat_lab;
   TLorentzVector p_plus_q=pBeam_lab + q_lab;
   // boost to HCM
   TLorentzRotation boost=TLorentzRotation(-1.0*p_plus_q.BoostVector());
   TLorentzVector pBoost=boost*pBeam_lab;
   TVector3 axis=pBoost.BoostVector();
   // rotate away y-coordinate
   boost.RotateZ(-axis.Phi());
   // rotate away x-coordinate
   boost.RotateY(M_PI-axis.Theta());
   return boost;
}

void eC_src(const int nEvents = 40000, TString filename=""){

	TChain *tree = new TChain("EICTree");
	tree->Add("/eicdata/eic0003/ztu/BeAGLE_Packages/"+filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		double pxf = event->pxf;
		double pyf = event->pyf;
		double pzf = event->pzf;

		double k = sqrt(pxf*pxf+pyf*pyf+pzf*pzf);
		nk->Fill( k );

	}

	// for(int j = 0; j < nk->GetNbinsX(); j++){
	// 	double value = nk->GetBinContent(j+1);
	// 	double width = nk->GetBinWidth(j+1);
	// 	double bincenter = nk->GetBinCenter(j+1);
	// 	double error = nk->GetBinError(j+1);

	// 	nk->SetBinContent(j+1,  (value/(width*bincenter*bincenter)) );
	// 	nk->SetBinError(j+1, (error/(width*bincenter*bincenter)) );
	// }

	nk->SetMarkerStyle(20);
	TFile output("../rootfiles/"+filename+"_src_Beagle.root","RECREATE");
	nk->Write();




}