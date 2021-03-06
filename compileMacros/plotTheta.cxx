#include "hist.h"
#include "PRINT4VECTOR.h"

using namespace erhic;
using namespace std;

TH2D* thetaNeutronVsthetaProton = new TH2D("thetaNeutronVsthetaProton",";#theta_{proton};#theta_{neutron}",100,0,0.01,100,0,0.01);
TH1D* deltaPhiLAB = new TH1D("deltaPhiLAB",";#phi_{n}#minus#phi_{p}",100,-2*PI-0.3,2*PI+0.3);
TH1D* deltaPhiION = new TH1D("deltaPhiION",";#phi_{n}#minus#phi_{p}",100,-2*PI-0.3,2*PI+0.3);

TH1D* px_neutron = new TH1D("px_neutron",";px",1000,-10,10);
TH1D* py_neutron = new TH1D("py_neutron",";py",1000,-10,10);
TH1D* pz_neutron = new TH1D("pz_neutron",";pz",1000,-10,10);

TH1D* px_proton = new TH1D("px_proton",";px",1000,-10,10);
TH1D* py_proton = new TH1D("py_proton",";py",1000,-10,10);
TH1D* pz_proton = new TH1D("pz_proton",";pz",1000,-10,10);

TH2D* pxVspxF_nucleon = new TH2D("pxVspxF_nucleon",";px;pxf",1000,-1,1,1000,-1,1);
TH2D* pyVspyF_nucleon = new TH2D("pyVspyF_nucleon",";py;pyf",1000,-1,1,1000,-1,1);
TH2D* pzVspzF_nucleon = new TH2D("pzVspzF_nucleon",";pz;pzf",1000,-1,1,1000,-1,1);

TH2D* relativePxPy_nucleon = new TH2D("relativePxPy_nucleon",";px;py",1000,-3,3,1000,-3,3);
TH1D* relativePz_nucleon = new TH1D("relativePz_nucleon",";pz",1000,-3,3);
TH1D* sNN_nucleon = new TH1D("sNN_nucleon",";s_{NN}",1000,0,15);

void plotTheta(int nEvents, TString inputFilename){

	TChain *tree = new TChain("EICTree");
	// tree->Add("../../EICTree/eD_Jpsidiffnodecay_EICTree/eD_18x135_Q2_1_10_y_0.01_0.95_tau_7_noquench_kt=ptfrag=0.32_Shd1_ShdFac=1.32_Jpsidiffnodecay_test40k_"+inputFilename+".root" ); // Wild cards are allowed e.g. tree.Add("*.root" );
	tree->Add("/eicdata/eic0003/ztu/EICTree/eD_FSI/"+inputFilename+".root" ); // Wild cards are allowed e.g. tree.Add("*.root" );

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
		double pF2 = pxf*pxf + pyf*pyf + pzf*pzf;
		
		/*hard-coded cuts*/
		//if( pF2 < 0.3025 || pF2 > 0.36 ) continue;
		if( event_process != 91 ) continue;
		//if( fabs(t_hat) > 0.1 ) continue;
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
	
		//Deuteron
		double pztarg_1 = 135;
		double pztarg_2 = 135;

		double Atarg = branch_atarg->GetValue(0,0);
		double pz_total = pztarg_1+pztarg_2;
		double total_energy = sqrt(pz_total*pz_total + MASS_DEUTERON*MASS_DEUTERON);

		TLorentzVector total4Mom_deuteron(0., 0., pz_total, total_energy);

		/* lorentz boost incoming particle*/
		double gamma_ion = total_energy/MASS_DEUTERON;
		double bz = pz_total/(gamma_ion*MASS_DEUTERON);

		TVector3 b;

		particle_4mom_proton.Boost(0,0,-bz);
		particle_4mom_proton.Boost(b);

		particle_4mom_neutron.Boost(0,0,-bz);
		particle_4mom_neutron.Boost(b);

		particle_4mom_photon.Boost(0,0,-bz);
		particle_4mom_photon.Boost(b);
      
		/*begin hand rotation*/
		TVector3 proton_v3 = particle_4mom_proton.Vect();
		TVector3 neutron_v3 = particle_4mom_neutron.Vect();
		TVector3 photon_v3 = particle_4mom_photon.Vect();

		double aa = proton_v3.Angle(photon_v3);
		double bb = neutron_v3.Angle(photon_v3);

		/*
		Use the angle between proton and photon, the projection (cos(aa)) is the z component of new proton
		then force the py = 0, and px > 0
		*/
		double mag2 = proton_v3.Mag2();
		double proton_pz = proton_v3.Mag()*TMath::Cos(aa);
		double proton_py = 0.0;
		double proton_px = sqrt(mag2 - proton_pz*proton_pz - proton_py*proton_py);

		/*
		Build new proton 4Vector
		*/
		TVector3 proton_v3_new(proton_px, proton_py, proton_pz);
		TLorentzVector particle_4mom_proton_new;
		particle_4mom_proton_new.SetVectM(proton_v3_new, MASS_PROTON);
		
		px_proton->Fill( proton_px );
		py_proton->Fill( proton_py );
		pz_proton->Fill( proton_pz );
		
		/*
		Proton cross photon gives the negative y axis direction, then make it unit vector;
		Use the y-axis and photon (z) to obtain the x unit vector;
		*/
		TVector3 y_axis = proton_v3.Cross(photon_v3);//x CROSS z = -y
		TVector3 y_unit = -y_axis.Unit(); // flip sign
		TVector3 x_axis = y_unit.Cross(photon_v3);// y CROSS z = x
		TVector3 x_unit = x_axis.Unit();
		
		/*Build 3Vector of new neutron*/
		mag2 = neutron_v3.Mag2();
		double neutron_pz = neutron_v3.Mag()*TMath::Cos(bb);
		double cc = neutron_v3.Angle(y_unit);
		double neutron_py = neutron_v3.Mag()*TMath::Cos(cc);//the neutron 3Vector's projection on y unit vector.
		double neutron_px = sqrt(mag2 - neutron_pz*neutron_pz - neutron_py*neutron_py);
		
		if( neutron_v3.Dot(x_unit) > 0 ) neutron_px = neutron_px;
		else if( neutron_v3.Dot(x_unit) < 0 ) neutron_px = -neutron_px;
		else cout << "wrong angle" << endl;

		/*Build new neutron 4Vector*/
		TVector3 neutron_v3_new(neutron_px, neutron_py, neutron_pz);
		TLorentzVector particle_4mom_neutron_new;
		particle_4mom_neutron_new.SetVectM(neutron_v3_new, MASS_NEUTRON);

		px_neutron->Fill( neutron_px );
		py_neutron->Fill( neutron_py );
		pz_neutron->Fill( neutron_pz );
		
		PhiDist_proton->Fill( particle_4mom_proton_new.Phi() );
		PhiDist_neutron->Fill( particle_4mom_neutron_new.Phi() );

 		deltaPhiION->Fill( particle_4mom_neutron_new.Phi() -  particle_4mom_proton_new.Phi() );

 		/*correlation between spectator momentum and pF*/

 		pxVspxF_nucleon->Fill( particle_4mom_proton_new.Px(), -pxf);
 		pyVspyF_nucleon->Fill( particle_4mom_proton_new.Py(), -pyf);
 		pzVspzF_nucleon->Fill( particle_4mom_proton_new.Pz(), -pzf);

 		relativePz_nucleon->Fill( particle_4mom_proton_new.Pz() - particle_4mom_neutron_new.Pz() );
 		relativePxPy_nucleon->Fill( particle_4mom_proton_new.Px() - particle_4mom_neutron_new.Px(), particle_4mom_proton_new.Py() - particle_4mom_neutron_new.Py() );
	
		TLorentzVector s = particle_4mom_proton_new + particle_4mom_neutron_new;
		double s_NN = s.Mag2();
		sNN_nucleon->Fill( s_NN );
	}

	TString outfilename;
	outfilename = "_kinematics_eD.root";

   	TFile output("../rootfiles/"+inputFilename+outfilename,"RECREATE");

   	thetaNeutronVsthetaProton->Write();
   	deltaEtadeltaPhi->Write();
   	deltaPhiLAB->Write();
   	deltaPhiION->Write();

   	px_neutron->Write();
   	py_neutron->Write();
   	pz_neutron->Write();

   	px_proton->Write();
   	py_proton->Write();
   	pz_proton->Write();

   	PhiDist_neutron->Write();
   	PhiDist_proton->Write();

   	pxVspxF_nucleon->Write();
   	pyVspyF_nucleon->Write();
   	pzVspzF_nucleon->Write();

   	relativePz_nucleon->Write();
   	relativePxPy_nucleon->Write();
   	sNN_nucleon->Write();



}