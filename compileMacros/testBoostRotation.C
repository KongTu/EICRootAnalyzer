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


TLorentzRotation BoostToHCM(TLorentzVector const &eBeam_lab,
                            TLorentzVector const &pBeam_lab,
                            TLorentzVector const &eScat_lab) {
   TLorentzVector q_lab=eBeam_lab - eScat_lab;
   TLorentzVector p_plus_q=pBeam_lab + q_lab;
   // boost to HCM
   TLorentzRotation boost=TLorentzRotation(-1.0*p_plus_q.BoostVector());
   TLorentzVector pBoost=boost*pBeam_lab;
   TVector3 axis=pBoost.BoostVector();

   // rotate away x-coordinate
   boost.RotateY(M_PI-axis.Theta());

   TLorentzVector pBoost_escat=boost*eScat_lab;
   TVector3 axis_escat=pBoost_escat.BoostVector();

   // rotate away y-coordinate
   boost.RotateZ(-axis_escat.Phi());
   

   return boost;
}

TLorentzRotation BoostToHCM_base(TLorentzVector const &eBeam_lab,
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

void testBoostRotation(const int nEvents = 40000){


	//input from BeAGLE root files
	TFile * output = new TFile("../rootfiles/test_boostRotation.root","recreate");
	
	TChain *tree = new TChain("EICTree");
	tree->Add("/gpfs02/eic/DATA/BeAGLE/eAu/DIS/18x110_Q2_1_100/EICTree/18x110_Q2_1_100_batch_1_10M/*.root" );
	// tree->Add("/gpfs02/eic/ztu/Analysis/BeAGLE/eD_Tagged_DIS/18x110_Q2_10_100_noINC/eD_Tagged_DIS_100M_batch_2/*.root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	TH1D* h_nk = new TH1D("h_nk",";nk",100,0,1);
	
	TH1D* h_ptStar = new TH1D("h_ptStar",";ptStar",100,0,20);
	TH1D* h_ptStar2 = new TH1D("h_ptStar2",";ptStar",100,0,20);
	TH1D* h_ptStar_after = new TH1D("h_ptStar_after",";ptStar",100,0,20);
	TH1D* h_pt = new TH1D("h_pt",";pt",100,0,20);
	
	TH1D* h_etaStar = new TH1D("h_etaStar",";etaStar",100,-10,10);
	TH1D* h_etaStar2 = new TH1D("h_etaStar2",";etaStar",100,-10,10);
	TH1D* h_eta = new TH1D("h_eta",";eta",100,-10,10);
	
	TH1D* h_phiStar = new TH1D("h_phiStar",";phiStar",100,-PI,PI);
	TH1D* h_phiStar2 = new TH1D("h_phiStar2",";phiStar",100,-PI,PI);
	TH1D* h_phi = new TH1D("h_phi",";phi",100,-PI,PI);

	TH1D* h_zhad = new TH1D("h_zhad",";z",100,0,1);
	TH2D* h_zhadVsPtStar = new TH2D("h_zhadVsPtStar",";z;ptStar",100,0,1,100,0,20);
	TH2D* h_phiStar2D = new TH2D("h_phiStar2D",";phiStar;phiStar_new",100,-PI,PI,100,-PI,PI);


	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);

		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		double pznucl = event->pznucl;
		double Atarg = event->Atarg;
		double pztarg_total = pztarg*Atarg;
		double pznucl_total = pznucl*Atarg;

		double pxf = event->pxf;
		double pyf = event->pyf;
		double pzf = event->pzf;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep+0.00051*0.00051));
		TLorentzVector d_beam(0.,0.,pztarg,sqrt(pztarg*pztarg+0.93891*0.93891));
		TLorentzVector e_scattered(0.,0.,0.,0.);

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
		
		//event process and kinematic phase space
		if( event_process != 99 ) continue;
		if( trueQ2 < 20.||trueQ2 > 100. ) continue;
		if( trueY > 0.95  || trueY < 0.01 ) continue;
				
		//HERA inclusive cross section
		vector<TLorentzVector> list_of_particles;
		for(int j(0); j < nParticles; ++j ) {
			const erhic::ParticleMC* particle = event->GetTrack(j);
			int pdg = particle->GetPdgCode();
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
			double pt = particle->GetPt();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			int status = particle->GetStatus();
			double mass = particle->GetM();
			int charge = particle->eA->charge;
			int NoBAM = particle->eA->NoBam;
			TLorentzVector ppart = particle->Get4Vector();
			if( index == 3 ) {
				e_scattered.SetPtEtaPhiM(pt,eta,phi,0.00051);
			}
			if( status!=1 ) continue;
			if( particle->GetParentIndex()==3 )continue;
			if(pt<0.15||fabs(eta)>1.6||charge==0) continue;

			list_of_particles.push_back(particle->Get4Vector());
		}
		list_of_particles.push_back(e_scattered);

		TLorentzRotation boost_HCM = BoostToHCM(e_beam,d_beam,e_scattered);
		TLorentzRotation boost_HCM_base = BoostToHCM_base(e_beam,d_beam,e_scattered);

		for(unsigned j=0;j<list_of_particles.size();j++){
			
			TLorentzVector hstar =  boost_HCM*list_of_particles[j];
			TLorentzVector hstar2 =  boost_HCM_base*list_of_particles[j];
			
			if(j==list_of_particles.size()-1){
				cout << "e' pt = " << list_of_particles[j].Pt() << endl;
				cout << "e' eta = " << list_of_particles[j].Eta() << endl;
				cout << "e' phi = " << list_of_particles[j].Phi() << endl;
				cout << "e' ptStar = " << hstar.Pt() << endl;
				cout << "e' etaStar = " << hstar.Eta() << endl;
				cout << "e' phiStar = " << hstar.Phi() << endl;
			}
			else{
				double zhad_value = d_beam.Dot(list_of_particles[j]) / d_beam.Dot(e_beam-e_scattered);
				TLorentzVector part_new;
				part_new.SetPtEtaPhiM(list_of_particles[j].Pt(),list_of_particles[j].Eta(),list_of_particles[j].Phi(),0.13957);
				zhad_value = d_beam.Dot(part_new) / d_beam.Dot(e_beam-e_scattered);

				h_zhad->Fill(zhad_value);

				h_ptStar->Fill(hstar.Pt());
				h_etaStar->Fill(hstar.Eta());
				h_phiStar->Fill(hstar.Phi());

				h_ptStar2->Fill(hstar2.Pt());
				h_etaStar2->Fill(hstar2.Eta());
				h_phiStar2->Fill(hstar2.Phi());

				h_phiStar2D->Fill(hstar.Phi(),hstar2.Phi());

				h_zhadVsPtStar->Fill(zhad_value,hstar.Pt());
				if( zhad_value > 0.2 ) h_ptStar_after->Fill(hstar.Pt());

				h_pt->Fill(list_of_particles[j].Pt());
				h_eta->Fill(list_of_particles[j].Eta());
				h_phi->Fill(list_of_particles[j].Phi());
			}
		}
		
	}

	output->Write();
	output->Close();


}