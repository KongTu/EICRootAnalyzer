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


void eD_Tagged_DIS(const int nEvents = 40000, TString filename="Output_input_temp_91"){


	//input from BeAGLE root files
	TFile * output = new TFile("../rootfiles/eD_Tagged_DIS_Beagle.root","recreate");
	
	TChain *tree = new TChain("EICTree");
	// tree->Add("/gpfs02/eic/ztu/Analysis/BeAGLE/eD_Tagged_DIS/18x110_Q2_10_100/eD_Tagged_DIS_100M_batch_1/"+filename+".root" );
	tree->Add("/gpfs02/eic/ztu/Analysis/BeAGLE/eD_Tagged_DIS/18x110_Q2_1_10/eD_Tagged_DIS_100M_batch_5/"+filename+".root" );

	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	//all constants
	double totalXSection   = 0.00056464908244711964;; //mb
	double nEventsTotal        = 250084.0;
	double Lint = nEventsTotal/totalXSection; // mb^{-1}
	double alpha2 = TMath::Power((1./137),2);
	double twopi = 2*PI;
	double mbToGeV_m2 = 2.56819;
	double Q2binwidth = 3.0-2.0;

	//alex's xbj binning
	double xBinsArray[] = {0.0001, 0.0002, 0.0004, 0.0007, 0.001, 0.002, 0.004, 0.007, 0.01, 0.02, 0.04, 0.07, 0.1};
	double xBinsWidth[12];
	for(int bin=0;bin<12;bin++){
		xBinsWidth[bin] = xBinsArray[bin+1]-xBinsArray[bin];
	}
	TH1D* h_xbj[12];
	for(int bin=0;bin<12;bin++){
		h_xbj[bin] = new TH1D(Form("h_xbj_%d",bin),Form("h_xbj_%d",bin),1,0,1);
	}
	//
	TH1D* h_HERA_Q2_10_14 = new TH1D("h_HERA_Q2_10_14","h_HERA_Q2_10_14",12,xBinsArray);
	TH1D* h_HERA_Q2_10_14_alpha_1 = new TH1D("h_HERA_Q2_10_14_alpha_1","h_HERA_Q2_10_14_alpha_1",100,0.00001,0.1);
	TH1D* h_HERA_Q2_10_14_alpha_2 = new TH1D("h_HERA_Q2_10_14_alpha_2","h_HERA_Q2_10_14_alpha_2",100,0.00001,0.1);
	TH1D* h_alpha_spec = new TH1D("h_alpha_spec","h_alpha_spec",100,0,2);
	TH1D* h_nk = new TH1D("h_nk","h_nk",100,0,2);
	double bin_width = h_HERA_Q2_10_14->GetBinWidth(1);
	
	double alpha_binning[161];
	for(int ibin=0;ibin<161;ibin++){
		alpha_binning[ibin] = 0.4+ibin*0.01;
	}	
	TH1D* h_HERA_Q2_10_14_x007_009_alpha[160];
	TH1D* h_alpha_spec_everybin[160];
	for(int ibin=0;ibin<160;ibin++){
	 	h_HERA_Q2_10_14_x007_009_alpha[ibin] = new TH1D(Form("h_HERA_Q2_10_14_x007_009_alpha_%d",ibin),Form("h_HERA_Q2_10_14_x007_009_alpha_%d",ibin),100,0,0.15);
		h_alpha_spec_everybin[ibin] = new TH1D(Form("h_alpha_spec_everybin_%d",ibin),Form("h_alpha_spec_everybin_%d",ibin),100,0,2);
	}

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		if( (i%10000)==0 ) cout << "#Events = "<< i << endl;
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
		double photon_flux = event->GetPhotonFlux();
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();
		int struck_nucleon = event->nucleon;
		double nk_event = sqrt(pxf*pxf+pyf*pyf+pzf*pzf);
		h_nk->Fill( nk_event );//sanity check for my wavefunction;
		TLorentzVector spectator_4vect_irf;
		double Espec = 0.;
		if( struck_nucleon == 2212 ){
			Espec = sqrt(nk_event*nk_event+MASS_NEUTRON*MASS_NEUTRON);
		}
		else{
			Espec = sqrt(nk_event*nk_event+MASS_PROTON*MASS_PROTON);
		}
		spectator_4vect_irf.SetPxPyPzE(-pxf,-pyf,-pzf,Espec);

		//event process and kinematic phase space
		if( struck_nucleon != 2212 ) continue; //proton only
		if( event_process != 99 ) continue;
		if( trueQ2 < 2.  || trueQ2 > 3. ) continue;
		if( trueY > 0.95  || trueY < 0.01 ) continue;

		cout << "Event #"<<i<< "with xbj = " << trueX << endl;
		cout << "Event #"<<i<< "with Q2 = " << trueQ2 << endl;
		cout << "Event #"<<i<< "with y = " << trueY << endl;
		cout << "Event #"<<i<< "with struck_nucleon = " << struck_nucleon << endl;

		//HERA inclusive cross section
		double event_weight = 1.;
		double Yc = 1. + TMath::Power((1-trueY),2);
		for(int j(0); j < nParticles; ++j ) {
			const erhic::ParticleMC* particle = event->GetTrack(j);
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
			double pt = particle->GetPt();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			if( index == 3 ) {
				e_scattered.SetPtEtaPhiM(pt,eta,phi,0.00051);
				// e_scattered = ppart;
			}
		}
		TLorentzVector qbeam = e_beam - e_scattered;

		// double xd = trueQ2 / (2*d_beam.Dot(qbeam));
		// double gamma2 = (4.*TMath::Power(MASS_DEUTERON,2)*TMath::Power(xd,2)) / trueQ2;
		// double epsilon = (1. - trueY - gamma2*TMath::Power(trueY/2.,2)) / (1. - trueY + TMath::Power(trueY,2)/2. + gamma2*TMath::Power(trueY/2.,2) );
		// double compare = TMath::Power( trueY, 2) / (1. - epsilon);
		//test for two different flux factor
		// cout << "compare ~ " << compare << "   Yc ~ " << Yc << endl;
		for(int bin=0;bin<12;bin++){
			if(trueX>xBinsArray[bin]&&trueX<xBinsArray[bin+1]){
				bin_width=xBinsWidth[bin];
				h_xbj[bin]->Fill(trueX);
			}
		}
		event_weight = (TMath::Power(trueQ2,2)*trueX) / (twopi*alpha2*Yc);
		event_weight = event_weight * (mbToGeV_m2)/(Lint*bin_width*Q2binwidth);
		//fill HERA inclusive cross section for Q2(10,13) GeV**2:
		h_HERA_Q2_10_14->Fill( trueX, event_weight );
		//x bin [0.007,0.009]
		//alpha below is wrong.
		double Pplus = (spectator_4vect_irf.E() + spectator_4vect_irf.Pz()) / sqrt(2);
		double PdPlus = MASS_DEUTERON / sqrt(2);
		double alpha_spec = 2*Pplus / PdPlus;

		if(alpha_spec < 0.9){
			h_HERA_Q2_10_14_alpha_1->Fill( trueX, event_weight );
		}
		if(alpha_spec > 1.1 ){
			h_HERA_Q2_10_14_alpha_2->Fill( trueX, event_weight );
		}
		// double pt2 = pxf*pxf+pyf*pyf;
		// double alpha_spec_binwidth = -1; // will have to be rewritten by 20 alpha bins
		// double xbinwidth = (0.003-0.002);
		// double pt2binwidth = h_HERA_Q2_10_14_x007_009_alpha[0]->GetBinWidth(1);
		// h_alpha_spec->Fill( alpha_spec );
		// if( trueX > 0.003 || trueX < 0.002 ) continue;
		
		// int alpha_bin_index = 0;
		// for(int ibin=0;ibin<160;ibin++){
		// 	if( alpha_spec>alpha_binning[ibin] && alpha_spec<alpha_binning[ibin+1] ){
		// 		alpha_bin_index = ibin;
		// 		alpha_spec_binwidth = alpha_binning[ibin+1] - alpha_binning[ibin];
		// 	}
		// }
		// double event_weight_alphaPt2 = alpha_spec*(16.*TMath::Power(PI,1)*(TMath::Power(trueQ2,2)*trueX)) / (alpha2*Yc);
		// event_weight_alphaPt2 = event_weight_alphaPt2 * (mbToGeV_m2/(Lint*Q2binwidth*xbinwidth*pt2binwidth*alpha_spec_binwidth));
		// //filling all alpha bins
		// h_HERA_Q2_10_14_x007_009_alpha[alpha_bin_index]->Fill(pt2, event_weight_alphaPt2 );
		// h_alpha_spec_everybin[alpha_bin_index]->Fill( alpha_spec );

	}

	output->Write();
	output->Close();


}