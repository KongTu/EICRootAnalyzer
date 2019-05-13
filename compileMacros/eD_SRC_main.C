#include "hist.h"//define all the histograms
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
using namespace std;
using namespace erhic;

#define MASS_MUON  0.1056

TH1D* sPN = new TH1D("sPN","sPN",1000,2,14);
TH1D* h_trk = new TH1D("h_trk","h_trk",4000,0,4000);

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

void eD_SRC_main(const int nEvents = 40000){

	TChain *tree = new TChain("EICTree");
	tree->Add("/eicdata/eic0003/ztu/BeAGLE_devK/eD_src_1.root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		double Atarg = event->Atarg;
		double pztarg_total = pztarg*Atarg;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep));//neglecting e mass
		TLorentzVector d_beam(0.,0.,pztarg_total,sqrt(pztarg_total*pztarg_total+MASS_DEUTERON*MASS_DEUTERON));
		TLorentzVector e_scattered(0.,0.,0.,0.);

		double boostv =  d_beam.Pz()/d_beam.E();
		TVector3 b;

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
		
		if( event_process != 91 ) continue;
		if( trueQ2 < 1. ) continue;
		if( trueY > 0.85 || trueY < 0.05 ) continue;

		int nParticles_process = 0;
		TLorentzVector p_4vect, n_4vect;
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

			if( index == 3 ) {
				e_scattered.SetPtEtaPhiM(pt,eta,phi,mass);
			}
			if( status != 1 ) continue;
			
			TLorentzVector ppart = particle->Get4Vector();
			ppart.Boost(0,0,-boostv);
			ppart.Boost(b);

			if(pdg == 2212) p_4vect = ppart;
			if(pdg == 2112) n_4vect = ppart;

			nParticles_process++;

		} // end of particle loop

		h_trk->Fill( nParticles_process );
		sPN->Fill( (p_4vect+n_4vect).Mag() );

	}

	TFile output("../rootfiles/eD_SRC_main_Beagle.root","RECREATE");
	sPN->Write();
	h_trk->Write();



}