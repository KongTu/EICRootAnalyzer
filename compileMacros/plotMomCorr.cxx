#include "hist.h"
#include "PRINT4VECTOR.h"

using namespace erhic;
using namespace std;

TH2D* thetaNeutronVsthetaProton = new TH2D("thetaNeutronVsthetaProton",";#theta_{proton};#theta_{neutron}",100,0,0.01,100,0,0.01);
TH1D* deltaPhiLAB = new TH1D("deltaPhiLAB",";#phi_{n}#minus#phi_{p}",100,-2*PI-0.3,2*PI+0.3);
TH1D* deltaPhiION = new TH1D("deltaPhiION",";#phi_{n}#minus#phi_{p}",100,-2*PI-0.3,2*PI+0.3);

TH1D* px_electron = new TH1D("px_electron",";px",1000,-100,100);
TH1D* py_electron = new TH1D("py_electron",";py",10,-10,10);
TH1D* pz_electron = new TH1D("pz_electron",";pz",1000,-10000,10000);

TH1D* px_neutron = new TH1D("px_neutron",";px",500,-5,5);
TH1D* py_neutron = new TH1D("py_neutron",";py",500,-5,5);
TH1D* pz_neutron = new TH1D("pz_neutron",";pz",500,-5,5);

TH1D* px_proton = new TH1D("px_proton",";px",500,-5,5);
TH1D* py_proton = new TH1D("py_proton",";py",500,-5,5);
TH1D* pz_proton = new TH1D("pz_proton",";pz",500,-5,5);

TH2D* pxFVspx_nucleon = new TH2D("pxFVspx_nucleon",";px;pxf",1000,-2,2,1000,-2,2);
TH2D* pyFVspy_nucleon = new TH2D("pyFVspy_nucleon",";py;pyf",1000,-2,2,1000,-2,2);
TH2D* pzFVspz_nucleon = new TH2D("pzFVspz_nucleon",";pz;pzf",1000,-2,2,1000,-2,2);

TH2D* pxFVspx_spectator = new TH2D("pxFVspx_spectator",";px;-pxf",1000,-2,2,1000,-2,2);
TH2D* pyFVspy_spectator = new TH2D("pyFVspy_spectator",";py;-pyf",1000,-2,2,1000,-2,2);
TH2D* pzFVspz_spectator = new TH2D("pzFVspz_spectator",";pz;-pzf",1000,-2,2,1000,-2,2);

TH1D* pt_spectator = new TH1D("pt_spectator",";px",500,0,5);
TH1D* phi_spectator = new TH1D("phi_spectator",";py",500,-2*PI,2*PI);

TH1D* pt_nucleon = new TH1D("pt_nucleon",";px",500,0,5);
TH1D* phi_nucleon = new TH1D("phi_nucleon",";py",500,-2*PI,2*PI);


void plotMomCorr(int nEvents, TString inputFilename, double pFmin_, double pFmax_ ){

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
		TLorentzVector particle_4mom_electron;

		double pxf = branch_pxf->GetValue(0,0);
		double pyf = branch_pyf->GetValue(0,0);
		double pzf = branch_pzf->GetValue(0,0);
		double pF = sqrt(pxf*pxf + pyf*pyf + pzf*pzf);
		
		/*hard-coded cuts*/
		if( event_process != 91 ) continue;
		if( pF < pFmin_ || pF > pFmax_ ) continue;
		if( struck_nucleon != 2112 ) continue;

		for(int j(0); j < nParticles; ++j ) {

			const erhic::ParticleMC* particle = event->GetTrack(j);

			int pdg = particle->GetPdgCode();
			int status = particle->GetStatus();
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.

			if( index == 3){//get scatter e'

				particle_4mom_electron = particle->Get4Vector();
			}
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

		particle_4mom_electron.Boost(0,0,-bz);
		particle_4mom_electron.Boost(b);

		particle_4mom_proton.Boost(0,0,-bz);
		particle_4mom_proton.Boost(b);

		particle_4mom_neutron.Boost(0,0,-bz);
		particle_4mom_neutron.Boost(b);

		particle_4mom_photon.Boost(0,0,-bz);
		particle_4mom_photon.Boost(b);
      
		/*begin hand rotation*/
		TVector3 electron_v3 = particle_4mom_electron.Vect();
		TVector3 proton_v3 = particle_4mom_proton.Vect();
		TVector3 neutron_v3 = particle_4mom_neutron.Vect();
		TVector3 photon_v3 = particle_4mom_photon.Vect();

		double aa = electron_v3.Angle(photon_v3);
		double bb = proton_v3.Angle(photon_v3);
		double cc = neutron_v3.Angle(photon_v3);

		/*
		Use the angle between electron and photon, the projection (cos(aa)) is the z component of new electron
		then force the py = 0, and px > 0
		*/
		double mag2 = electron_v3.Mag2();
		double electron_pz = electron_v3.Mag()*TMath::Cos(aa);
		double electron_py = 0.0;
		double electron_px = sqrt(mag2 - electron_pz*electron_pz - electron_py*electron_py);

		px_electron->Fill( electron_px );
		py_electron->Fill( electron_py );
		pz_electron->Fill( electron_pz );

		/*
		Build new electron 4Vector
		*/
		TVector3 electron_v3_new(electron_px, electron_py, electron_pz);
		TLorentzVector particle_4mom_electron_new;
		particle_4mom_electron_new.SetVectM(electron_v3_new, 0.00051);

		/*
		Elctron cross photon gives the negative y axis direction, then make it unit vector;
		Use the y-axis and photon (z) to obtain the x unit vector;
		*/
		TVector3 y_axis = electron_v3.Cross(photon_v3);//x CROSS z = -y
		TVector3 y_unit = -y_axis.Unit(); // flip sign
		TVector3 x_axis = y_unit.Cross(photon_v3);// y CROSS z = x
		TVector3 x_unit = x_axis.Unit();

		/*Build 3Vector of new proton*/
		double proton_pz = proton_v3.Mag()*TMath::Cos(bb);
		double dd = proton_v3.Angle(y_unit);
		double proton_py = proton_v3.Mag()*TMath::Cos(dd);//the proton 3Vector's projection on y unit vector.
		dd = proton_v3.Angle(x_unit);
		double proton_px = proton_v3.Mag()*TMath::Cos(dd);//similar to px
		
		px_proton->Fill( proton_px );
		py_proton->Fill( proton_py );
		pz_proton->Fill( proton_pz );

		/*Build new proton 4Vector*/
		TVector3 proton_v3_new(proton_px, proton_py, proton_pz);
		TLorentzVector particle_4mom_proton_new;
		particle_4mom_proton_new.SetVectM(proton_v3_new, MASS_PROTON);

		/*Build 3Vector of new neutron*/
		double neutron_pz = neutron_v3.Mag()*TMath::Cos(cc);
		dd = neutron_v3.Angle(y_unit);
		double neutron_py = neutron_v3.Mag()*TMath::Cos(dd);//the neutron 3Vector's projection on y unit vector.
		dd = neutron_v3.Angle(x_unit);
		double neutron_px = neutron_v3.Mag()*TMath::Cos(dd);//similar to px
		
		px_neutron->Fill( neutron_px );
		py_neutron->Fill( neutron_py );
		pz_neutron->Fill( neutron_pz );

		/*Build new neutron 4Vector*/
		TVector3 neutron_v3_new(neutron_px, neutron_py, neutron_pz);
		TLorentzVector particle_4mom_neutron_new;
		particle_4mom_neutron_new.SetVectM(neutron_v3_new, MASS_NEUTRON);

		//done with ratation;
		PhiDist_proton->Fill( particle_4mom_proton_new.Phi() );
		PhiDist_neutron->Fill( particle_4mom_neutron_new.Phi() );

 		deltaPhiION->Fill( particle_4mom_neutron_new.Phi() -  particle_4mom_proton_new.Phi() );

 		/*correlation between spectator momentum and pF*/
 		pxFVspx_spectator->Fill( particle_4mom_proton_new.Px(), -pxf);
 		pyFVspy_spectator->Fill( particle_4mom_proton_new.Py(), -pyf);
 		pzFVspz_spectator->Fill( particle_4mom_proton_new.Pz(), -pzf);

 		/*correlation between struck momentum and pF*/
 		pxFVspx_nucleon->Fill( particle_4mom_neutron_new.Px(), pxf);
 		pyFVspy_nucleon->Fill( particle_4mom_neutron_new.Py(), pyf);
 		pzFVspz_nucleon->Fill( particle_4mom_neutron_new.Pz(), pzf);

 		pt_spectator->Fill( particle_4mom_proton_new.Pt() );
 		phi_spectator->Fill( particle_4mom_proton_new.Phi() );

		pt_nucleon->Fill( particle_4mom_neutron_new.Pt() );
 		phi_nucleon->Fill( particle_4mom_neutron_new.Phi() );


	}

	TString outfilename;
	outfilename = "_kinematics_eD.root";

	std::stringstream ss,nn;
	ss << pFmin_;
	std::string pFmin_string;
	ss >> pFmin_string;
	nn << pFmax_;
	std::string pFmax_string;
	nn >> pFmax_string;

   	TFile output("../rootfiles/"+inputFilename+"_"+pFmin_string+"_"+pFmax_string+outfilename,"RECREATE");

   	thetaNeutronVsthetaProton->Write();
   	deltaEtadeltaPhi->Write();
   	deltaPhiLAB->Write();
   	deltaPhiION->Write();

   	px_electron->Write();
   	py_electron->Write();
   	pz_electron->Write();

   	px_neutron->Write();
   	py_neutron->Write();
   	pz_neutron->Write();

   	px_proton->Write();
   	py_proton->Write();
   	pz_proton->Write();

   	PhiDist_neutron->Write();
   	PhiDist_proton->Write();

   	pxFVspx_nucleon->Write();
   	pyFVspy_nucleon->Write();
   	pzFVspz_nucleon->Write();

   	pxFVspx_spectator->Write();
   	pyFVspy_spectator->Write();
   	pzFVspz_spectator->Write();

   	pt_spectator->Write();
   	phi_spectator->Write();

	pt_nucleon->Write();
   	phi_nucleon->Write();


}