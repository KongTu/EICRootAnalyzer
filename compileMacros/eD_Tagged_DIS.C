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

#define MASS_JPSI 	  3.09688
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
// #define MASS_DEUTERON 1.8756129
#define MASS_DEUTERON 1.8751019071673038
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

#define MASS_MUON  0.1056


void eD_Tagged_DIS(const int nEvents = 40000, TString filename="eD_dis_Tagged_highQ2_1M"){


	//input from BeAGLE root files
	TFile * output = new TFile("../rootfiles/eD_Tagged_DIS_Beagle.root","recreate");
	
	TChain *tree = new TChain("EICTree");
	tree->Add("/gpfs02/eic/DATA/KongTu/forDouglasJLab/"+filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	//all constants
	double totalXSection   = .0000450463252; //mb
	double nEventsTotal        = 500257.0;
	double Lint = nEventsTotal/totalXSection; // mb^{-1}
	double alpha2 = TMath::Power((1./137),2);
	double twopi = 2*PI;
	double mbToGeV_m2 = 2.56819;

	TH1D* h_nk = new TH1D("h_nk","h_nk",100,0,2);
	TH1D* h_HERA_Q2_10_13 = new TH1D("h_HERA_Q2_10_13","h_HERA_Q2_10_13",1000,0.00001,0.1);

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);

		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		double pznucl = event->pznucl;
		double Atarg = event->Atarg;
		double pztarg_total = pztarg*Atarg;

		double pxf = event->pxf;
		double pyf = event->pyf;
		double pzf = event->pzf;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep+0.00051*0.00051));
		TLorentzVector d_beam(0.,0.,pztarg_total,sqrt(pztarg_total*pztarg_total+MASS_DEUTERON*MASS_DEUTERON));
		TLorentzVector e_scattered(0.,0.,0.,0.);
		
		//boost vector for lab <--> d rest frame
		TVector3 b = d_beam.BoostVector();

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
		int struck_nucleon = event->nucleon;
		double nk_event = sqrt(pxf*pxf+pyf*pyf+pzf*pzf);

		//HERA inclusive cross section
		double event_weight = 1.;
		double Yc = 1.-TMath::Power((1-trueY),2);
		event_weight = (TMath::Power(trueQ2,2)*trueX) / (twopi*alpha2*Yc);
		event_weight = (event_weight*mbToGeV_m2) / Lint;
		double bin_width = h_HERA_Q2_10_13->GetBinWidth(1);
		event_weight = event_weight / bin_width;

		//event process and kinematic phase space
		if( event_process != 99 ) continue;
		if( trueQ2 < 10.  || trueQ2 > 13. ) continue;
		if( trueY > 0.95  || trueY < 0.01 ) continue;

		//try HERA inclusive cross section:
		h_HERA_Q2_10_13->Fill( trueX, event_weight );

		if( trueX > 0.009 || trueX < 0.007 ) continue;

		h_nk->Fill( nk_event );

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
			TLorentzVector ppart = particle->Get4Vector();

			if( index == 3 ) {
				e_scattered.SetPtEtaPhiM(pt,eta,phi,0.00051);
				// e_scattered = ppart;
			}
			if( status != 1 ) continue;


		} // end of particle loop


	}

	output->Write();
	output->Close();

}