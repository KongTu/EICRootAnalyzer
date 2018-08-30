#include "hist.h"

using namespace erhic;
using namespace std;


TH2D* thetaNeutronVsthetaProton = new TH2D("thetaNeutronVsthetaProton",";#theta_{proton};#theta_{neutron}",100,0,0.01,100,0,0.01);
TH1D* deltaPhiLAB = new TH1D("deltaPhiLAB",";#phi_{n}#minus#phi_{p}",100,-2*PI-0.3,2*PI+0.3);
TH1D* deltaPhiION = new TH1D("deltaPhiION",";#phi_{n}#minus#phi_{p}",100,-2*PI-0.3,2*PI+0.3);

void plotTheta(int nEvents, TString inputFilename){

	TChain *tree = new TChain("EICTree");
	tree->Add("../../EICTree/eD_Jpsidiffnodecay_EICTree/eD_18x135_Q2_1_10_y_0.01_0.95_tau_7_noquench_kt=ptfrag=0.32_Shd1_ShdFac=1.32_Jpsidiffnodecay_test40k_"+inputFilename+".root" ); // Wild cards are allowed e.g. tree.Add("*.root" );

	EventPythia* event(NULL);// = new EventPythia;

	// EventBase* event(NULL);
	// EventBeagle* event_beagle(NULL);

	tree->SetBranchAddress("event", &event ); // Note &event, not event.

	//Using the TBranchElement is a hack to access the BeAGLE information.       
	TBranchElement* branch_atarg = (TBranchElement*) tree->GetBranch("Atarg");
	TBranchElement* branch_pz = (TBranchElement*) tree->GetBranch("pztarg");
	TBranchElement* branch_pzlep = (TBranchElement*) tree->GetBranch("pzlep");
	TBranchElement* branch_pxf = (TBranchElement*) tree->GetBranch("pxf");
	TBranchElement* branch_pyf = (TBranchElement*) tree->GetBranch("pyf");
	TBranchElement* branch_pzf = (TBranchElement*) tree->GetBranch("pzf");

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);

		//event information:
		double trueQ2 = event->GetTrueQ2();
		double trueW2 = event->GetTrueW2();
		double trueX = event->GetTrueX();
		double trueY = event->GetTrueY();
		double trueNu = event->GetTrueNu();
		double s_hat = event->GetHardS();
		double t_hat = event->t_hat;
		double u_hat = event->GetHardU();
		double photon_flux = event->GetPhotonFlux();
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();
		int struck_nucleon = event->nucleon;
		
		int nParticles_process = 0;

		TLorentzVector particle_4mom_proton;
		TLorentzVector particle_4mom_neutron;
		TLorentzVector particle_4mom_jpsi;
		TLorentzVector particle_4mom_photon;

		double pxf = branch_pxf->GetValue(0,0);
		double pyf = branch_pyf->GetValue(0,0);
		double pzf = branch_pzf->GetValue(0,0);
		double pF = pxf*pxf + pyf*pyf + pzf*pzf;
		

		if( pF < 0.3025 || pF > 0.36 ) continue;
		if( event_process != 91 ) continue;
		if( fabs(t_hat) > 0.1 ) continue;
		if( struck_nucleon != 2112 ) continue;

		for(int j(0); j < nParticles; ++j ) {

			const erhic::ParticleMC* particle = event->GetTrack(j);

			int pdg = particle->GetPdgCode();
			int status = particle->GetStatus();
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.

			if( index == 4 ){ //get gamma 4-momentum:

				particle_4mom_photon = particle->Get4Vector();
			}
			if( status != 1 ) continue; //only stable final-state particles 
			if( pdg == 443 ){//Jpsi

				particle_4mom_jpsi = particle->Get4Vector();
			}
			if( pdg == 2212 ){//proton

				particle_4mom_proton = particle->Get4Vector();
			}
			if( pdg == 2112 ){//neutron

				particle_4mom_neutron = particle->Get4Vector();
			}

			nParticles_process++;

		} // end of particle loop

	 
		thetaNeutronVsthetaProton->Fill( particle_4mom_proton.Theta(), particle_4mom_neutron.Theta() );
		deltaEtadeltaPhi->Fill(particle_4mom_neutron.Eta() -  particle_4mom_proton.Eta(), particle_4mom_neutron.Phi() -  particle_4mom_proton.Phi());
		deltaPhiLAB->Fill( particle_4mom_neutron.Phi() -  particle_4mom_proton.Phi() );
	
	}



	TString outfilename;
	outfilename = "_kinematics_eD.root";

   	TFile output("../rootfiles/"+inputFilename+outfilename,"RECREATE");

   	thetaNeutronVsthetaProton->Write();
   	deltaEtadeltaPhi->Write();
   	deltaPhiLAB->Write();




}