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
#include "TLorentzRotation.h"
#include "TBranchElement.h"

#define PI            3.1415926

#define MASS_JPSI 	  3.09688
#define MASS_PROTON   0.93827208816
#define MASS_NEUTRON  0.93956542052
#define MASS_NUCLEON  0.93891875 //.93891875
#define MASS_DEUTERON 1.8756129 // 1.8756134 (most precise so far)
// #define MASS_DEUTERON 1.8751019071673038 (not precise enough)!
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

void eD_TaggedEMC_DIS( const int nEvents = 1e6 ){


	//input from BeAGLE root files
	TFile * output = new TFile("../rootfiles/eD_TaggedEMC_DIS_Beagle.root","recreate");
	
	TChain *tree = new TChain("EICTree");
	tree->Add("/gpfs02/eic/ztu/Analysis/BeAGLE/eD_Tagged_DIS/5x41_Q2_10/output/Output_eD_5x41_DIS_1M.root" );

	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	//all constants
	double factorInLumi = nEvents / 1e6 ;
	double totalXSection   = 1.8501226230870374E-002;; //mb
	double nEventsTotal    = nEvents*factorInLumi;
	double Lint = nEventsTotal/totalXSection; // mb^{-1}
	double alpha2 = TMath::Power((1./137),2);
	double twopi = 2*PI;
	double mbToGeV_m2 = 2.56819;
	double Q2binwidth = 15-10.;
	//alex's xbj binning
	double xBinsArray[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7};
	double xBinsWidth[6];
	for(int bin=0;bin<6;bin++){
		xBinsWidth[bin] = xBinsArray[bin+1]-xBinsArray[bin];
	}
	const int nPt2=100; 
	TH1D* h_xbj[6];
	TH1D* h_pt2[6][nPt2];
	for(int bin=0;bin<6;bin++){
		h_xbj[bin] = new TH1D(Form("h_xbj_%d",bin),Form("h_xbj_%d",bin),1,0,1);
		for(int ipt=0;ipt<nPt2;ipt++){
			h_pt2[bin][ipt] = new TH1D(Form("h_pt2_%d_%d",bin,ipt),Form("h_pt2_%d_%d",bin,ipt),1,0,1);
			h_pt2[bin][ipt]->Sumw2();
		}
	}
	//
	TH1D* h_HERA_Q2_10_15 = new TH1D("h_HERA_Q2_10_15","h_HERA_Q2_10_15",6,xBinsArray);
	h_HERA_Q2_10_15->Sumw2();
	
	TH1D* h_alpha_spec = new TH1D("h_alpha_spec","h_alpha_spec",100,0,2);
	TH1D* h_nk = new TH1D("h_nk","h_nk",100,0,2);
	double bin_width = h_HERA_Q2_10_15->GetBinWidth(1);
	
	double alpha_binning[6]={0.75,0.85,0.95,1.05,1.15,1.25};
	TH1D* h_HERA_Q2_10_15_x_alpha[6][5];
	TH1D* h_alpha_spec_everybin[5];
	for(int ibin=0;ibin<5;ibin++){
		h_alpha_spec_everybin[ibin] = new TH1D(Form("h_alpha_spec_everybin_%d",ibin),Form("h_alpha_spec_everybin_%d",ibin),100,0,2);
		for(int jbin=0;jbin<6;jbin++){
		 	h_HERA_Q2_10_15_x_alpha[jbin][ibin] = new TH1D(Form("h_HERA_Q2_10_15_x_alpha_%d_%d",jbin,ibin),Form("h_HERA_Q2_10_15_x_alpha_%d_%d",jbin,ibin),nPt2,0,0.15);
			h_HERA_Q2_10_15_x_alpha[jbin][ibin]->Sumw2();
		}
	}

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		if( (i%10000)==0 ) cout << "#Events = "<< i << endl;
		// cout << "#Events = "<< i << endl;
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		double pznucl = event->pznucl;
		double Atarg = event->Atarg;
		double pztarg_total = pztarg*Atarg;
		double EMC_weight = event->User3;

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
		if( struck_nucleon == 2212 ){Espec = sqrt(nk_event*nk_event+MASS_NEUTRON*MASS_NEUTRON);}
		else{Espec = sqrt(nk_event*nk_event+MASS_PROTON*MASS_PROTON);}
		
		//event process and kinematic phase space
		if( struck_nucleon != 2212 ) continue; //proton only
		if( trueQ2 < 10.  || trueQ2 > 15. ) continue;
		if( trueY > 0.95  || trueY < 0.01 ) continue;

		//HERA inclusive cross section
		double event_weight = 1.;
		double Yc = 1. + TMath::Power((1-trueY),2);
		for(int j(0); j < nParticles; ++j ) {
			const erhic::ParticleMC* particle = event->GetTrack(j);
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
			double pt = particle->GetPt();
			int status = particle->GetStatus();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			int pdg = particle->GetPdgCode();
			int orig = particle->GetParentIndex();
			if( index == 3 ) {
				e_scattered.SetPtEtaPhiM(pt,eta,phi,0.00051);
				// e_scattered = ppart;
			}
			if( status!= 1 ) continue;
		}

		TLorentzVector qbeam = e_beam - e_scattered;
		spectator_4vect_irf.SetPxPyPzE(-pxf,-pyf,-pzf,Espec);
		TLorentzVector trueSpect_lab = spectator_4vect_irf;
		for(int bin=0;bin<6;bin++){
			if(trueX>xBinsArray[bin]&&trueX<xBinsArray[bin+1]){
				bin_width=xBinsWidth[bin];
				h_xbj[bin]->Fill(trueX);
			}
		}
		event_weight = (TMath::Power(trueQ2,2)*trueX) / (twopi*alpha2*Yc);
		event_weight = event_weight * (mbToGeV_m2)/(Lint*bin_width*Q2binwidth);
		
		//below a test for one xbj bin
		trueSpect_lab.SetPtEtaPhiM(trueSpect_lab.Pt(),trueSpect_lab.Eta(), trueSpect_lab.Phi(),MASS_NUCLEON);
		double Pplus = (trueSpect_lab.E() + trueSpect_lab.Pz()) / sqrt(2);
		double PdPlus = MASS_DEUTERON / sqrt(2);
		double alpha_spec = 2*Pplus / PdPlus;
		double pt2 = trueSpect_lab.Pt()*trueSpect_lab.Pt();
		
		//fill inclusive.
		h_HERA_Q2_10_15->Fill( trueX, event_weight );

		double pt2weight = 1.0;
		double alpha_spec_binwidth = -1; // will have to be rewritten by alpha defined bins
		double pt2binwidth = h_HERA_Q2_10_15_x_alpha[0][0]->GetBinWidth(1);
		h_alpha_spec->Fill( alpha_spec );
		int x_bin_index = -1;
		for(int ix=0;ix<6;ix++){
			if( trueX>xBinsArray[ix]&&trueX<xBinsArray[ix+1] ){
				x_bin_index=ix;
			}
		}
		if( x_bin_index < 0 ) continue;
		int alpha_bin_index = 0;
		for(int ibin=0;ibin<5;ibin++){
			if( alpha_spec>alpha_binning[ibin] && alpha_spec<alpha_binning[ibin+1] ){
				alpha_bin_index = ibin;
				alpha_spec_binwidth = alpha_binning[ibin+1] - alpha_binning[ibin];
			}
		}
		//to calculate <pt2> per bin
		for(int ipt=0;ipt<nPt2;ipt++){
			if(pt2 > ipt*pt2binwidth && pt2 <(ipt+1)*pt2binwidth){
				h_pt2[x_bin_index][ipt]->Fill(pt2,pt2weight);
			}
		}
		//8*PI is correct, NOT 16*PI.
		double event_weight_alphaPt2 = alpha_spec*(8.*TMath::Power(PI,1)*(TMath::Power(trueQ2,2)*trueX)) / (alpha2*Yc);
		event_weight_alphaPt2 = event_weight_alphaPt2 * (mbToGeV_m2/(Lint*Q2binwidth*bin_width*pt2binwidth*alpha_spec_binwidth));
		//filling all alpha bins
		h_HERA_Q2_10_15_x_alpha[x_bin_index][alpha_bin_index]->Fill(pt2, event_weight_alphaPt2*pt2weight );
		h_alpha_spec_everybin[alpha_bin_index]->Fill( alpha_spec );

	}

	output->Write();
	output->Close();


}