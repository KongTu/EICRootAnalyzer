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


void runInclusJpsi_BeAGLE(const TString filename="eA_TEST", const int nEvents = 2e5, TString system="eAu"){

	TChain *tree = new TChain("EICTree");
	tree->Add( filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	TFile* output=new TFile("../rootfiles/InclusJpsi_BeAGLE_"+system+".root","RECREATE");
	//histograms
	TH1D* h_W = new TH1D("h_W",";W (GeV)",3000,0,3000);
	TH1D* h_Jpsi_rapidity = new TH1D("h_Jpsi_rapidity",";y",100,-4,4);
	TH1D* h_Jpsi_pt = new TH1D("h_Jpsi_pt",";pt",100,0,10);
	TH1D* h_Jpsi_process = new TH1D("h_Jpsi_process",";process",200,0,200);
	TH1D* h_Jpsi_status = new TH1D("h_Jpsi_status",";status",100,0,100);

	TH1D* h_highPt_rapidity = new TH1D("h_highPt_rapidity",";y",100,-4,4);
	TH1D* h_highPt_pt = new TH1D("h_highPt_pt",";pt",100,0,30);
	TH1D* h_highPt_process = new TH1D("h_highPt_process",";process",200,0,200);

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		double Atarg = event->Atarg;
		double pztarg_total = pztarg*Atarg;
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
		if(event_process>90 && event_process<96) continue;

		h_W->Fill( sqrt(trueW2) );
		
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
			int charge= particle->eA->charge;
			if( pdg == 443 ) {
				h_Jpsi_rapidity->Fill( rap );
				h_Jpsi_pt->Fill( pt );
				h_Jpsi_process->Fill( event_process );
				h_Jpsi_status->Fill( status );
			}
			if( status!= 1) continue;
			if( pt < 3. ) continue;
			h_highPt_rapidity->Fill( rap );
			h_highPt_pt->Fill( pt );
			h_highPt_process->Fill( event_process );

		} // end of particle loop
	
	}

	output->Write();
	output->Close();

}
