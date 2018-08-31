#include "hist.h"//define all the histograms
#include "PRINT4VECTOR.h"
using namespace std;
using namespace erhic;

void run_SRCkicks(int nEvents, bool doKick, int CASE, TString inputFilename){

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
		double trueNu = event->nu;
		double s_hat = event->GetHardS();
		double t_hat = event->t_hat;
		double u_hat = event->GetHardU();
		double photon_flux = event->GetPhotonFlux();
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();
		int struck_nucleon = event->nucleon;

		cout << "trueQ2 " << trueQ2 << endl;
		cout << "trueNu " << trueNu << endl;
		
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

		if( event_process != 91 ) continue;

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

			statusHist.Fill( status ); 

			if( pdg == 22 ){
				TLorentzVector p = particle->Get4Vector();
				PRINT4VECTOR(p, true);
			}

			if( index == 4 ){ //get gamma 4-momentum:

			particle_4mom_photon = particle->Get4Vector();
			PRINT4VECTOR(particle_4mom_photon,true);
			}
			if( index == 3 ){
			particle_4mom_electron_prime = particle->Get4Vector();
			PRINT4VECTOR(particle_4mom_electron_prime,true);
			}
			if( status != 1 ) continue; //only stable final-state particles 
			if( pdg == 443 ){//Jpsi

			particle_4mom_jpsi_bKick = particle->Get4Vector();
			particle_4mom_jpsi = particle->Get4Vector();
			PRINT4VECTOR(particle_4mom_jpsi,true);
			}
			if( pdg == 2212 ){//proton

				//SetPtEtaPhiM(pt,eta,phi,mass);//this won't work if there is no pT 
			particle_4mom_proton_bKick = particle->Get4Vector();
			particle_4mom_proton = particle->Get4Vector();
			PRINT4VECTOR(particle_4mom_proton,true);
			cout << "Mass proton: " << particle->GetM() << endl;

			}
			if( pdg == 2112 ){//neutron

			particle_4mom_neutron_bKick = particle->Get4Vector();
			particle_4mom_neutron = particle->Get4Vector();
			PRINT4VECTOR(particle_4mom_neutron,true);
			cout << "Mass neutron: " << particle->GetM() << endl;

			}

			nParticles_process++;

		} // end of particle loop

		if( doKick ){ 

			t = particle_4mom_neutron_bKick + particle_4mom_proton_bKick + particle_4mom_jpsi_bKick;
			
			double kick_px = 0.;
			double kick_py = 0.;
			double kick_pz = 0.;

			TF1 *fa_pt = new TF1("fa_pt","[0]*TMath::Abs(TMath::Exp([1]*TMath::Abs(x)))",0,2);
			fa_pt->SetParameter(0,1);
			fa_pt->SetParameter(1,-3);
			double kick_pt = fa_pt->GetRandom();

			TF1 *num = new TF1("num","[0]*1",-1,1);
			num->SetParameter(0,1);
			double prob = num->GetRandom();

			TF1 *phiran = new TF1("phiran","[0]*1",-PI,PI);
			phiran->SetParameter(0,1);
			double phi_kick = phiran->GetRandom();

			if( prob > 0 ){
				kick_px = sqrt(kick_pt*kick_pt/(1+TMath::Tan(phi_kick)*TMath::Tan(phi_kick)));
			}
			else{
				kick_px = -sqrt(kick_pt*kick_pt/(1+TMath::Tan(phi_kick)*TMath::Tan(phi_kick)));
			}
	
			kick_py = kick_px*TMath::Tan(phi_kick);

			TF1 *fa_pz = new TF1("fa_pz","[0]*TMath::Abs(TMath::Exp([1]*TMath::Abs(x)))",-2,2);
			fa_pz->SetParameter(0,1);
			fa_pz->SetParameter(1,-1);
			
			kick_pz = fa_pz->GetRandom();

			px_dist->Fill( kick_px );
			py_dist->Fill( kick_py );
			pz_dist->Fill( kick_pz );

			pt_dist->Fill( kick_pt );
			phi_dist->Fill( phi_kick );

			//proton 3 momentum:
			double p_px = particle_4mom_proton_bKick.Px();
			double p_py = particle_4mom_proton_bKick.Py();
			double p_pz = particle_4mom_proton_bKick.Pz();
			double p_E = sqrt(p_px*p_px + p_py*p_py + p_pz*p_pz + MASS_PROTON*MASS_PROTON);

			//neutron 3 momentum:
			double n_px = particle_4mom_neutron_bKick.Px();
			double n_py = particle_4mom_neutron_bKick.Py();
			double n_pz = particle_4mom_neutron_bKick.Pz(); 
			double n_E = sqrt(n_px*n_px + n_py*n_py + n_pz*n_pz + MASS_NEUTRON*MASS_NEUTRON);

			//Jpsi 3 momentum:
			double j_px = particle_4mom_jpsi_bKick.Px();
			double j_py = particle_4mom_jpsi_bKick.Py();
			double j_pz = particle_4mom_jpsi_bKick.Pz();
			double j_E = sqrt(j_px*j_px + j_py*j_py + j_pz*j_pz + MASS_JPSI*MASS_JPSI);

			double p_py_prime = p_py; 
			double n_py_prime = n_py; 
			double j_py_prime = j_py; 

			double p_px_prime = p_px; 
			double n_px_prime = n_px; 
			double j_px_prime = j_px; 

			double p_pz_prime = p_pz; 
			double n_pz_prime = n_pz; 
			double j_pz_prime = j_pz; 

			if( CASE == 1 ){
				
				double E_min = 1000.0;
				double comp_min = 0.;
				double delta_min = 0.;
				double kappa_min = 0.;

				int i_min = 0;
				int j_min = 0;
				int k_min = 0;

				double comp_init = -1;
				double delta_init = -5;
				double kappa_init = -5;

				const int iteration_1 = 20;
				const int iteration_2 = 10;

				double comp[iteration_1];
				double delta[iteration_2];
				double kappa[iteration_2];

				for(int jter = 0; jter < iteration_1; jter++){

				 double temp = comp_init+0.1*jter;
				 comp[jter] = temp;  
				}

				for(int jter = 0; jter < iteration_2; jter++){

				 double temp = delta_init+1.*jter;
				 delta[jter] = temp;

				 temp = kappa_init+1.*jter;
				 kappa[jter] = temp;
				}

				for(int iter = 0; iter < iteration_2; iter++){//delta
				 for(int jter = 0; jter < iteration_1; jter++){//comp
				    for(int kter = 0; kter < iteration_2; kter++){//kappa
					
					if( struck_nucleon == 2212 ){

						p_py_prime = p_py + kick_py;
					    n_py_prime = n_py - kick_py + delta[iter];
					    j_py_prime = j_py - delta[iter];

					    p_px_prime = p_px + kick_px; 
					    n_px_prime = n_px - kick_px + kappa[kter];
					    j_px_prime = j_px - kappa[kter];

					    p_pz_prime = p_pz + kick_pz;
					    n_pz_prime = n_pz - kick_pz + comp[jter];
					    j_pz_prime = j_pz - comp[jter];

					}
					else{

						p_py_prime = p_py - kick_py + delta[iter];
					    n_py_prime = n_py + kick_py ;
					    j_py_prime = j_py - delta[iter];

					    p_px_prime = p_px - kick_px + kappa[kter]; 
					    n_px_prime = n_px + kick_px;
					    j_px_prime = j_px - kappa[kter];

					    p_pz_prime = p_pz - kick_pz + comp[jter];
					    n_pz_prime = n_pz + kick_pz;
					    j_pz_prime = j_pz - comp[jter];
					}

				    double p_E_prime = sqrt(p_px_prime*p_px_prime + p_py_prime*p_py_prime + p_pz_prime*p_pz_prime + MASS_PROTON*MASS_PROTON);
				    double n_E_prime = sqrt(n_px_prime*n_px_prime + n_py_prime*n_py_prime + n_pz_prime*n_pz_prime + MASS_NEUTRON*MASS_NEUTRON);
				    double j_E_prime = sqrt(j_px_prime*j_px_prime + j_py_prime*j_py_prime + j_pz_prime*j_pz_prime + MASS_JPSI*MASS_JPSI);

				    p3.SetPxPyPzE(p_px_prime,p_py_prime,p_pz_prime,p_E_prime);
				    p4.SetPxPyPzE(n_px_prime,n_py_prime,n_pz_prime,n_E_prime);
				    p5.SetPxPyPzE(j_px_prime,j_py_prime,j_pz_prime,j_E_prime);

				    k = p3+p4+p5;

				    double E_DIFF = t.E() - k.E();
				    double pz_DIFF = t.Pz() - k.Pz();

				       if( fabs(E_DIFF) < fabs(E_min) ) {

				          E_min = E_DIFF;
				          comp_min = comp[jter];
				          delta_min = delta[iter];
				          kappa_min = kappa[kter];

				          i_min = iter;
				          j_min = jter;
				          k_min = kter;  

					          particle_4mom_proton = p3;
					          particle_4mom_neutron = p4;
					          particle_4mom_jpsi = p5;
				       }

				    }//loop1
				}//loop2
			}//loop3
			  
			if( i_min == 0 || j_min == 0 || k_min == 0 || i_min == 9 || j_min == 19 || k_min == 9 ) continue;//hit the boundary continue;
			// cout << "iter: " << i_min << " jter: " << j_min << " kter: " << k_min << endl;
			// cout << "E diff: " << E_DIFF <<  " comp: " << comp_min << " delta: " << delta_min << " kappa: " << kappa_min << endl;
	
		}
		if( CASE == 2 ){

			double aa[1000];
			double bb[1000];
			double cc;

			for(int jter = 0; jter < 1000; jter++){

				 double temp = 0.001*jter;
				 aa[jter] = temp;  
				 bb[jter] = temp;
			}

			for(int iter = 0; iter < 1000; iter++){
				for(int jter = 0; jter < 1000; jter++){

					cc = ((1-aa[iter])*p_py + (1-bb[jter])*n_py + (bb[jter]-aa[iter])*kick_py)/j_py + 1.;
					p_py_prime = aa[iter]*(p_py + kick_py);
				    n_py_prime = bb[jter]*(n_py - kick_py);
				    j_py_prime = cc*j_py;

				    cc = ((1-aa[iter])*p_px + (1-bb[jter])*n_px + (bb[jter]-aa[iter])*kick_px)/j_px + 1.;
				    p_px_prime = aa[iter]*(p_px + kick_px); 
				    n_px_prime = bb[jter]*(n_px - kick_px);
				    j_px_prime = cc*j_px;

				    cc = ((1-aa[iter])*p_pz + (1-bb[jter])*n_pz + (bb[jter]-aa[iter])*kick_pz)/j_pz + 1.;
				    p_pz_prime = aa[iter]*(p_pz + kick_pz);
				    n_pz_prime = bb[jter]*(n_pz - kick_pz);
				    j_pz_prime = cc*j_pz;

				    double p_E_prime = sqrt(p_px_prime*p_px_prime + p_py_prime*p_py_prime + p_pz_prime*p_pz_prime + MASS_PROTON*MASS_PROTON);
				    double n_E_prime = sqrt(n_px_prime*n_px_prime + n_py_prime*n_py_prime + n_pz_prime*n_pz_prime + MASS_NEUTRON*MASS_NEUTRON);
				    double j_E_prime = sqrt(j_px_prime*j_px_prime + j_py_prime*j_py_prime + j_pz_prime*j_pz_prime + MASS_JPSI*MASS_JPSI);

				    p3.SetPxPyPzE(p_px_prime,p_py_prime,p_pz_prime,p_E_prime);
				    p4.SetPxPyPzE(n_px_prime,n_py_prime,n_pz_prime,n_E_prime);
				    p5.SetPxPyPzE(j_px_prime,j_py_prime,j_pz_prime,j_E_prime);

				    k = p3+p4+p5;

				    double E_DIFF = t.E() - k.E();
				    double px_DIFF = t.Px() - k.Px();
				    double py_DIFF = t.Py() - k.Py();
				    double pz_DIFF = t.Pz() - k.Pz();

				    if(E_DIFF < 0.7) cout << "E_DIFF = " << E_DIFF << endl;
					
				}
			}	


		}
		
	}//end of kick

		//fill histograms:
		total4Mom_outgoing = particle_4mom_proton + particle_4mom_neutron + particle_4mom_jpsi + particle_4mom_electron_prime;

		/*fill histograms*/
		energy_corr->Fill(total4Mom_incoming.E() - total4Mom_outgoing.E());
		//Jpsi:
		PtDist_Jpsi->Fill( particle_4mom_jpsi.Pt() );
		EtaDist_Jpsi->Fill( particle_4mom_jpsi.Eta() );
		PhiDist_Jpsi->Fill( particle_4mom_jpsi.Phi() );

		//proton
		PtDist_proton->Fill( particle_4mom_proton.Pt() );
		EtaDist_proton->Fill( particle_4mom_proton.Eta() );
		PhiDist_proton->Fill( particle_4mom_proton.Phi() );
		AngleVsMom_proton->Fill(particle_4mom_proton.P(), particle_4mom_proton.Theta()*1000.);

		//neutron:
		PtDist_neutron->Fill( particle_4mom_neutron.Pt() );
		EtaDist_neutron->Fill( particle_4mom_neutron.Eta() );
		PhiDist_neutron->Fill( particle_4mom_neutron.Phi() );
		AngleVsMom_neutron->Fill(particle_4mom_neutron.P(), particle_4mom_neutron.Theta()*1000.);

		//delta eta delta phi:
		deltaEtadeltaPhi->Fill( particle_4mom_proton.Eta()-particle_4mom_neutron.Eta(), particle_4mom_proton.Phi()-particle_4mom_neutron.Phi());

		//t_hat
		T_dist->Fill( t_hat );

		//small t, namely the momentum transfer to the struck nucleon (proton)
		TLorentzVector t1_proton = particle_4mom_proton_bKick - total4Mom_deuteron;//(p'-p)
		double t_proton_squared = t1_proton.Mag2();

		t1_dist->Fill( t_proton_squared );

		TLorentzVector t2_proton = particle_4mom_proton - total4Mom_deuteron;//(p'-p)
		t_proton_squared = t2_proton.Mag2();

		t2_dist->Fill( t_proton_squared );

		//T, momentum transfer from photon to Jpsi
		TLorentzVector T_Jpsi = particle_4mom_jpsi - particle_4mom_photon;//delta
		double T_Jpsi_squared = T_Jpsi.Mag2();

		T_Jpsi_dist->Fill( T_Jpsi_squared );

		particle_4mom = particle_4mom_proton + particle_4mom_neutron;

		double sNN = particle_4mom.Mag2();//center of mass energy squared

		E_CM->Fill( sqrt(sNN) );
		//end COM
		//if( fabs(t_proton_squared)  > 0.5 && fabs(t_proton_squared) < 5.0 && fabs(T_Jpsi_squared) < 0.5 ) 
		sNN_dist->Fill( sNN );

		//t vs T
		tVsT->Fill(T_Jpsi_squared, t_proton_squared);



   	} // end of event loop

	TString outfilename;
	if( doKick ) outfilename = "_SRCkicks_eD_kick.root";
	else outfilename = "_SRCkicks_eD_nokick.root";

   	TFile output("../rootfiles/"+inputFilename+outfilename,"RECREATE");

	T_dist->Write();//T_distribution in the selected range
	T_Jpsi_dist->Write();
	t1_dist->Write();//t_distribution in the selected range
	t2_dist->Write();//t_distribution in the selected range

	PtDist_Jpsi->Write();
	EtaDist_Jpsi->Write();
	PhiDist_Jpsi->Write();

	PtDist_proton->Write();
	EtaDist_proton->Write();
	PhiDist_proton->Write();

	PtDist_neutron->Write();
	EtaDist_neutron->Write();
	PhiDist_neutron->Write();

	AngleVsMom_proton->Write();
	AngleVsMom_neutron->Write();

	tVsT->Write();
	deltaEtadeltaPhi->Write();
	sNN_dist->Write();
	energy_corr->Write();
	px_dist->Write();
	py_dist->Write();
	pz_dist->Write();
	pt_dist->Write();
	phi_dist->Write();

}