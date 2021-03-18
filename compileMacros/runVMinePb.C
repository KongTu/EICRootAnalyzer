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


void runVMinePb(const TString filename="eA_TEST", const int nEvents = 40000){

	TChain *tree = new TChain("EICTree");
	tree->Add( filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	TFile* output = new TFile("output.root","RECREATE");
	TH1D* h_trueT = new TH1D("h_trueT",";-t (GeV^{2})", 100,0,1);
	TH2D* h_thetaVsMom[3];
	for(int k=0;k<3;k++){
		h_thetaVsMom[k] = new TH2D(Form("h_thetaVsMom_%d",k),";p (GeV);#theta (mrad)",2500,0,250,100,0,1500);
	}

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		int struck_nucleon = event->nucleon;
		double MASS_NUCLEON = MASS_PROTON;
		if( struck_nucleon==2112 ) MASS_NUCLEON = MASS_NEUTRON;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep+MASS_ELECTRON*MASS_ELECTRON));
		TLorentzVector p_beam(0.,0.,pztarg,sqrt(pztarg*pztarg+MASS_NUCLEON*MASS_NUCLEON));
		TLorentzVector e_scattered(0.,0.,0.,0.);

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
		
		double impact_parameter = event->b;
		double Tb = event->Thickness;
		double distance = event->d1st;
		int N_nevap = event->Nnevap;
		int N_pevap = event->Npevap;

		//event cuts
		if( event_process != 91 && event_process != 93 ) continue;
		if( trueQ2 < 1. || trueQ2 > 100. ) continue;
		if( trueY > 0.95 || trueY < 0.01 ) continue;

		//do analysis, or fill historgrams for event levels

		//particle loop
		bool hasJpsi = false;
		vector< double> angle_neutron, angle_proton, angle_photon;
		vector< double> momentum_neutron, momentum_proton, momentum_photon;
		for(int j(0); j < nParticles; ++j ) {

			const erhic::ParticleMC* particle = event->GetTrack(j);

			int pdg = particle->GetPdgCode();
			int status = particle->GetStatus();
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
			double pt = particle->GetPt();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			double rap = particle->GetRapidity();
			double mass = particle->GetM();
			double theta = particle->GetTheta(); 
			theta = theta*1000.0; //change to mrad;
			double mom = particle->GetP();

			//only stable particles or j/psi.
			if( status != 1 && pdg != 443 ) continue;

			if( pdg == 443 ){ //jpsi
				hasJpsi = true;
			}
			if( pdg == 2112 ){ // neutrons
				angle_neutron.push_back( theta );
				momentum_neutron.push_back( mom );
			}
			if( pdg == 2212 ){ // proton
				angle_proton.push_back( theta );
				momentum_proton.push_back( mom );
			}
			if( pdg == 22 ){ // photon
				angle_photon.push_back( theta );
				momentum_photon.push_back( mom );
			}
			//do analysis track-by-track

		} // end of particle loop

		if( hasJpsi ) {

			h_trueT->Fill( -t_hat );
			for(unsigned ipart=0;ipart<angle_neutron.size();ipart++){
				h_thetaVsMom[0]->Fill(momentum_neutron[ipart],angle_neutron[ipart]);
			}
			for(unsigned ipart=0;ipart<angle_proton.size();ipart++){
				h_thetaVsMom[1]->Fill(momentum_proton[ipart],angle_proton[ipart]);
			}
			for(unsigned ipart=0;ipart<angle_photon.size();ipart++){
				h_thetaVsMom[2]->Fill(momentum_photon[ipart],angle_photon[ipart]);
			}
		}
		//fill histograms
	}

	output->Write();
	output->Close();


}
