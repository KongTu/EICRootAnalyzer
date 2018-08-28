#include "hist.h"//define all the histograms

using namespace std;
using namespace erhic;

void run_SRCkicks(int nEvents, bool doKick, TString inputFilename){

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
      double t_hat = event->GetHardT();
      double u_hat = event->GetHardU();
      double photon_flux = event->GetPhotonFlux();
      int event_process = event->GetProcess();
      
      int nParticles = event->GetNTracks();
      //event t_hat

      // We now know the number of particles in the event, so loop over
      // the particles:
      int nParticles_process = 0;

      TLorentzVector particle_4mom;
      TLorentzVector t,k;
      TLorentzVector p3,p4,p5;

      TLorentzVector particle_4mom_proton_bKick;
      TLorentzVector particle_4mom_neutron_bKick;
      TLorentzVector particle_4mom_jpsi_bKick;

      TLorentzVector particle_4mom_proton;
      TLorentzVector particle_4mom_neutron;
      TLorentzVector particle_4mom_jpsi;

      TLorentzVector particle_4mom_photon;
      TLorentzVector particle_4mom_electron_prime;

      //if( event_process != 91 ) continue;
      
      /*E-M Conservation*/
      double pztarg_1 = 135.290727;
      double pztarg_2 = 135.103537;
      double pz_total = pztarg_1+pztarg_2;
      double total_energy = sqrt(pz_total*pz_total + MASS_DEUTERON*MASS_DEUTERON);
      //electron, neglect electron mass
      double pz_lepton = branch_pzlep->GetValue(0,0);
      double electron_mass = 0.00051;
      double total_lep_energy = sqrt(pz_lepton*pz_lepton + electron_mass*electron_mass);

      TLorentzVector total4Mom_deuteron(0., 0., pz_total, total_energy);
      TLorentzVector total4Mom_electron(0., 0., pz_lepton, total_lep_energy);

      TLorentzVector total4Mom_outgoing(0.,0.,0.,0.);
      TLorentzVector total4Mom_incoming = total4Mom_deuteron + total4Mom_electron;
      /*end*/

      cout << "total energy: " << total4Mom_incoming.E() << endl;

   } // for

}