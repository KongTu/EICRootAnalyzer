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
#define MASS_MUON  	  0.1056

using namespace std;
using namespace erhic;


TLorentzRotation RotateToLab(TLorentzVector const &eBeam_lab,
                            TLorentzVector const &pBeam_lab,
                            TLorentzVector const &eScat_lab) {
	
	TLorentzVector q_lab=eBeam_lab - eScat_lab;
	TLorentzVector q_irf=q_lab;
	TLorentzVector eScat_irf=eScat_lab;
	TVector3 pBoost=pBeam_lab.BoostVector();
	q_irf.Boost(-pBoost);
	eScat_irf.Boost(-pBoost);

	TLorentzRotation l;
	double angleTheta = q_irf.Theta();
	double anglePhi = eScat_irf.Phi();
	l.RotateY( angleTheta );
	l.RotateZ( anglePhi );

	return l;
}

void eD_Tagged_DIS_debug(const int nEvents = 40000){

	//input from BeAGLE root files
	TFile * output = new TFile("../rootfiles/eD_Tagged_DIS_Beagle_background.root","recreate");
	TChain *tree = new TChain("EICTree");
	tree->Add("/gpfs02/eic/ztu/Analysis/BeAGLE/eD_Tagged_DIS/18x110_Q2_10_100_noINC/eD_Tagged_DIS_1M_batch_1/*.root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

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
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();
		int struck_nucleon = event->nucleon;
		double nk_event = sqrt(pxf*pxf+pyf*pyf+pzf*pzf);
		double Espec = 0.;
		TLorentzVector trueSpect;
		if( struck_nucleon == 2212 ){
			Espec = sqrt(nk_event*nk_event+MASS_NEUTRON*MASS_NEUTRON);
		}
		else{
			Espec = sqrt(nk_event*nk_event+MASS_PROTON*MASS_PROTON);
		}
		trueSpect.SetPxPyPzE(-pxf,-pyf,-pzf,Espec);
		
		//event process and kinematic phase space
		if( event_process != 99 ) continue;
		if( trueQ2 < 10.  || trueQ2 > 13. ) continue;
		if( trueY > 0.95  || trueY < 0.01 ) continue;
				
		vector< TLorentzVector> saveListOfNucleons;
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
			TLorentzVector ppart = particle->Get4Vector();
			if( index == 3 ) {
				e_scattered.SetPtEtaPhiM(pt,eta,phi,0.00051);
			}
			if( status!=1 ) continue;
			if( TMath::Abs(pdg) == 2112 || TMath::Abs(pdg) == 2212 ) saveListOfNucleons.push_back( ppart );
		}

		TLorentzRotation rotateVector=RotateToLab(e_beam, d_beam, e_scattered);
		TLorentzVector trueSpect_lab = rotateVector*trueSpect;//rotation only
		trueSpect_lab.Boost(b);//longitudinal boost without rotation
		
		cout << "Event #" << i << " is a struck nucleon with ID: " << struck_nucleon << endl;
		cout << "True spectator pt " << trueSpect_lab.Pt() << " eta " << trueSpect_lab.Eta() << " phi " << trueSpect_lab.Phi() << " mass " << trueSpect_lab.M() << " total p " << trueSpect_lab.P() << endl;
		for(unsigned icand=0; icand<saveListOfNucleons.size(); icand++){
			cout << "candidate " << icand << " mass " << saveListOfNucleons[icand].M() 
			<< " pt " << saveListOfNucleons[icand].Pt() << " eta " << saveListOfNucleons[icand].Eta()  << " phi " << saveListOfNucleons[icand].Phi()  << " total p " << saveListOfNucleons[icand].P() <<endl;
		}
		saveListOfNucleons.clear();
	
	}

	output->Write();
	output->Close();


}