#include "hist.h"//define all the histograms
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
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
   // rotate away y-coordinate
   boost.RotateZ(-axis.Phi());
   // rotate away x-coordinate
   boost.RotateY(M_PI-axis.Theta());
   return boost;
}

void makeComparisonE665(const int nEvents = 40000){

	TChain *tree = new TChain("EICTree");
	tree->Add("/gpfs/mnt/gpfs02/eic/wanchang/BeAGLE/muXe/muXe_490x0_Q2_1_100_y_0.1_0.85_tau_7_Shd3_trigcut_US0_40k.root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);
	// EventPythia* event(NULL);
	// tree->SetBranchAddress("event", &event );

	//Using the TBranchElement is a hack to access the BeAGLE information.       
	

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		cout << "pzlep " << pzlep << endl;
		//event information:
		// double trueQ2 = event->GetTrueQ2();
		// double trueW2 = event->GetTrueW2();
		// double trueX = event->GetTrueX();
		// double trueY = event->GetTrueY();
		// double trueNu = event->nu;
		// double s_hat = event->GetHardS();
		// double t_hat = event->t_hat;
		// double u_hat = event->GetHardU();
		// double photon_flux = event->GetPhotonFlux();
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();
		int struck_nucleon = event->nucleon;
		if( event_process != 91 ) continue;

		int nParticles_process = 0;
		for(int j(0); j < nParticles; ++j ) {

			const erhic::ParticleMC* particle = event->GetTrack(j);

			int pdg = particle->GetPdgCode();
			int status = particle->GetStatus();
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
			double pt = particle->GetPt();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			double mass = particle->GetM();
			double theta = particle->GetTheta(); 
			theta = theta*1000.0; //change to mrad;
			double mom = particle->GetP();



			nParticles_process++;

		} // end of particle loop

		cout << "nParticles " << nParticles_process << endl;
	}



}