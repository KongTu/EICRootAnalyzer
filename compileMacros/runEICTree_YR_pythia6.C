#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include <eicsmear/erhic/EventBase.h>
#include <eicsmear/erhic/EventMC.h>
#include <eicsmear/erhic/EventPythia.h>
#include <eicsmear/erhic/Particle.h>
#include <eicsmear/erhic/ParticleMC.h>
#include <eicsmear/erhic/Pid.h>

#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TBranchElement.h"

#define PI            3.1415926

#define MASS_MUON     0.1056
#define MASS_ELECTRON 0.000511
#define MASS_JPSI 	  3.09688
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_DEUTERON 1.8756129
#define MASS_TRITON   2.7937167208086358
#define MASS_HE3      2.7937167208086358
#define MASS_ALPHA    3.7249556277448477
#define MASS_LI6      5.5874334416172715
#define MASS_C12      11.174866883234543
#define MASS_CA40     37.249556277448477
#define MASS_XE131    121.99229680864376
#define MASS_AU197    183.45406466643374
#define MASS_PB208    193.69769264273208

using namespace std;
using namespace erhic;


void runEICTree_YR_pythia6(const TString filename="ep_test_noycut_R=200_PARJ72", const int nEvents = 1e5){

	TChain *tree = new TChain("EICTree");
	tree->Add( filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	TFile* output = new TFile("output.root", "RECREATE");

	TH1D* hist_mult = new TH1D("hist_mult","",100,0,100);
	
	TH1D* hist_eta_all = new TH1D("hist_eta_all","",50,-5,8);
	TH1D* hist_eta_charged = new TH1D("hist_eta_charged","",50,-5,8);
	TH1D* hist_eta_neutral = new TH1D("hist_eta_neutral","",50,-5,8);

	TH1D* hist_etatheta_all = new TH1D("hist_etatheta_all","",50,-5,8);
	TH1D* hist_etatheta_charged = new TH1D("hist_etatheta_charged","",50,-5,8);
	TH1D* hist_etatheta_neutral = new TH1D("hist_etatheta_neutral","",50,-5,8);


	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;

		//event information:
		double trueQ2 = event->GetTrueQ2();
		double trueW2 = event->GetTrueW2();
		double trueX = event->GetTrueX();
		double trueY = event->GetTrueY();
		double trueNu = event->nu;
		double s_hat = event->GetHardS();
		double t_hat = event->t_hat;
		double u_hat = event->GetHardU();
		double photon_flux = event->GetPhotonFlux();
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();

		//event cuts - nothing.
		// if( trueQ2 < 1. ) continue;
		// if( trueY > 0.95 || trueY < 0.01 ) continue;

		//particle loop
		int multiplicity = 0;
		for(int j(0); j < nParticles; ++j ) {

			const erhic::ParticleMC* particle = event->GetTrack(j);

			int pdg = particle->GetPdgCode();
			int status = particle->GetStatus();
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
			double pt = particle->GetPt();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			double theta = particle->GetTheta(); 
			
			TLorentzVector p = particle->Get4Vector();
			int charge = particle->eA->charge;
			int NoBAM = particle->eA->NoBam;

			if( status != 1 ) continue;
			if( NoBAM != 0 && NoBAM != 2 ) continue;
			if( index == 1 || index == 2 || index == 3 ) continue;
			if( pdg == 12 || pdg == 14 || pdg == 16 ) continue;
			// if( pt<0.2 ) continue;
			multiplicity++;

			hist_eta_all->Fill( eta );
			if(charge != 0) hist_eta_charged->Fill( eta );
			if(charge == 0) hist_eta_neutral->Fill( eta );

			hist_etatheta_all->Fill( eta );
			if(charge != 0) hist_etatheta_charged->Fill( eta );
			if(charge == 0) hist_etatheta_neutral->Fill( eta );

		} // end of particle loop

		hist_mult->Fill( multiplicity );
		//fill histograms
	}

	output->Write();
	output->Close();


}
