#include "hist.h"//define all the histograms
#include "PRINT4VECTOR.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include <sstream>
#include <string>

using namespace std;
using namespace erhic;

#define MASS_MUON  0.1056

void eD_photo_main(const int nEvents = 40000, TString filename=""){

	//just naming in the output file, only show ZDC parameters. 
	TString settings = "eD_photo";

	//input from BeAGLE root files
	TFile * output = new TFile("../rootfiles/"+filename+"_"+settings+"_main_Beagle.root","recreate");
	//histograms to be saved
	TH1D* t_truth = new TH1D("t_truth",";-t'(GeV)",100,0,2);
	TH1D* t_ZDC = new TH1D("t_ZDC",";-t'(GeV)",100,0,2);
	//define TTree
	TChain *tree = new TChain("EICTree");
	tree->Add("/gpfs02/eic/ztu/BeAGLE/BeAGLE_devK_SRC/"+filename+".root" );
	
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
		double trueNu = event->nu;
		double s_hat = event->GetHardS();
		double t_hat = event->t_hat;
		double u_hat = event->GetHardU();
		double photon_flux = event->GetPhotonFlux();
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();
		int struck_nucleon = event->nucleon;
		double nk_event = sqrt(pxf*pxf+pyf*pyf+pzf*pzf);
		
		if( event_process != 91 && event_process != 93 ) continue;
		if( trueQ2 > 1. ) continue;
		if( trueY > 0.95 || trueY < 0.05 ) continue;

		t_truth->Fill( -t_hat );
		
		int nParticles_ZDC = 0;
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
			if( theta < 5 && (pdg==2112 || pdg==22 || pdg==111)) nParticles_ZDC++;
		} // end of particle loop

		if( nParticles_ZDC>0 ) t_ZDC->Fill( -t_hat );

	}

	output->Write();
	output->Close();

}