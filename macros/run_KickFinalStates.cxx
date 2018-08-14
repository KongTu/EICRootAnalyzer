#include "hist.h"//define all the histograms

TH1D* t1_dist = new TH1D("t1_dist",";t1", 200,-5,5);
TH1D* t2_dist = new TH1D("t2_dist",";t2", 200,-5,5);
TH1D* energy_corr = new TH1D("energy_corr",";E_{in} - E_{out}",600,-30,30);

void run_KickFinalStates( int nEvents, bool doKick, TString inputFilename ) {
   

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

	TLorentzVector total4Mom_electron;
	TLorentzVector total4Mom_deuteron;
	TLorentzVector total4Mom_outgoing;
	TLorentzVector total4Mom_incoming;

	for(int i(0); i < nEvents; ++i ) {
      
      // Read the next entry from the tree.
      tree->GetEntry(i);

      //event information:
      double trueQ2 = event->GetTrueQ2();
      double trueW2 = event->GetTrueW2();
         W2->Fill(trueW2);

      double trueX = event->GetTrueX();
      double trueY = event->GetTrueY();
      double trueNu = event->GetTrueNu();
      double s_hat = event->GetHardS();
      double t_hat = event->GetHardT();
      double u_hat = event->GetHardU();
      double photon_flux = event->GetPhotonFlux();
         photonFlux->Fill( photon_flux );

      int event_process = event->GetProcess();

      //Deuteron 4 momentum
      double pztarg = branch_pz->GetValue(0,0);
      double Atarg = branch_atarg->GetValue(0,0);
      double pz_total = pztarg*Atarg;
      double total_energy = sqrt(pz_total*pz_total + MASS_DEUTERON*MASS_DEUTERON);
      
      //electron 4 momentum
      double pz_lepton = branch_pzlep->GetValue(0,0);
      //double electron_mass = 0.00051;
      double total_lep_energy = sqrt(pz_lepton*pz_lepton+0.00051*0.00051);

      total4Mom_deuteron.SetPxPyPzE(0., 0., pz_total, total_energy);
      total4Mom_electron.SetPxPyPzE(0., 0., pz_lepton, total_lep_energy);

      total4Mom_outgoing.SetPxPyPzE(0.,0.,0.,0.);
      total4Mom_incoming = total4Mom_deuteron + total4Mom_electron;

      Q2VsX->Fill(trueX, trueQ2);
      W2VsFlux->Fill(photon_flux, trueW2);

      // The event contains a vector (array) of particles.
      int nParticles = event->GetNTracks();
      //event t_hat

      // We now know the number of particles in the event, so loop over
      // the particles:

      int nParticles_process = 0;

      TLorentzVector particle_4mom;

      TLorentzVector particle_4mom_proton_bKick;
      TLorentzVector particle_4mom_neutron_bKick;

      TLorentzVector particle_4mom_proton;
      TLorentzVector particle_4mom_neutron;

      TLorentzVector particle_4mom_photon;
      TLorentzVector particle_4mom_Jpsi;
      //TLorentzVector particle_4mom_electron;

      if( event_process != 91 ) continue;

      for(int j(0); j < nParticles; ++j ) {
         
         const erhic::ParticleMC* particle = event->GetTrack(j);
    
         int pdg = particle->GetPdgCode();
         int status = particle->GetStatus();
         int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
         double pt = particle->GetPt();
         double eta = particle->GetEta();
         double phi = particle->GetPhi();
         double theta = particle->GetTheta(); 
         theta = theta*1000.0; //change to mrad;
         double mom = particle->GetP();

         statusHist.Fill( status ); 

         if( index == 3 ){ //get scattered electron

         	//particle_4mom_electron = particle->PxPyPzE();

         }
         if( index == 4 ){ //get gamma 4-momentum:

            particle_4mom_photon = particle->PxPyPzE(); 
         }

         if( status != 1 ) continue; //only stable final-state particles 

            if( pdg == 443 ){//Jpsi

               PtDist_Jpsi->Fill( pt );
               EtaDist_Jpsi->Fill( eta );
               PhiDist_Jpsi->Fill( phi );

               double pt_jpsi = pt;//store Jpsi pt
               
               particle_4mom_Jpsi = particle->PxPyPzE();

            }
            if( pdg == 2212 ){//proton

            	particle_4mom_proton_bKick = particle->PxPyPzE();
            	particle_4mom_proton = particle->PxPyPzE();

				if( doKick ){

					TF1 *fa = new TF1("fa","[0]*TMath::Exp([1]*x)",0,3);
					fa->SetParameter(0,1);
					fa->SetParameter(1,-5);

					double kick = fa->GetRandom();

					double p_eta = eta;
					double p_phi = phi;
					double p_E = particle->GetE();
					double p_M = particle->GetM();
					double p_pT = pt;
					double p_pz = particle->GetPz();

					p_pT += kick;
					p_pz = sqrt(p_E*p_E - p_M*p_M - p_pT*p_pT);
					//p_E = sqrt(p_E*p_E + kick*kick);
					p_eta = TMath::ATanH(p_pz/sqrt(p_pz*p_pz+p_pT*p_pT));

					particle_4mom_proton.SetPtEtaPhiE(p_pT, p_eta, p_phi, p_E);

				}

				double new_pt = particle_4mom_proton.Pt();
				double new_eta = particle_4mom_proton.Eta();
				double new_phi = particle_4mom_proton.Phi();
				double new_theta = particle_4mom_proton.Theta();
				new_theta = new_theta*1000;
				double new_mom = particle_4mom_proton.P();

				PtDist_proton->Fill( new_pt );
				EtaDist_proton->Fill( new_eta );
				PhiDist_proton->Fill( new_phi );

				PtVsEta_proton->Fill(new_eta, new_pt);
				AngleVsMom_proton->Fill(new_mom, new_theta);

				double theta_proton = new_theta;//store proton angle
				double pt_proton = new_pt;//store proton pt

            }
            if( pdg == 2112 ){//neutron

				particle_4mom_neutron_bKick = particle->PxPyPzE();
				particle_4mom_neutron = particle->PxPyPzE();

            }

            nParticles_process++;

	} // end of particle loop

	if( nParticles_process != 4 ) continue;

	TLorentzVector bKick_PN;
	if( doKick ){

		bKick_PN = particle_4mom_proton_bKick + particle_4mom_neutron_bKick;
		particle_4mom_neutron = bKick_PN - particle_4mom_proton;//modify neutron kinematics.
	}

	// total4Mom_outgoing += particle_4mom_neutron;
	// total4Mom_outgoing += particle_4mom_proton;
	// total4Mom_outgoing += particle_4mom_Jpsi;
	// total4Mom_outgoing += particle_4mom_electron;

	//refill neutron kinematics:
	PtDist_neutron->Fill( particle_4mom_neutron.Pt() );
	EtaDist_neutron->Fill( particle_4mom_neutron.Eta() );
	PhiDist_neutron->Fill( particle_4mom_neutron.Phi() );
	
	double theta_neutron = particle_4mom_neutron.Theta();//store neutron angle
	theta_neutron = fabs(theta_neutron*1000);

	PtVsEta_neutron->Fill(particle_4mom_neutron.Eta(), particle_4mom_neutron.Pt());
	AngleVsMom_neutron->Fill(particle_4mom_neutron.P(), theta_neutron);

	//t_hat
	T_dist->Fill( t_hat );

	//small t, namely the momentum transfer to the struck nucleon (proton)
	TLorentzVector t1_proton = particle_4mom_proton_bKick - total4Mom_deuteron;//(p'-p)
	double t_proton_squared = t1_proton.Mag2();

	t1_dist->Fill( t_proton_squared );

	TLorentzVector t2_proton = particle_4mom_proton - total4Mom_deuteron;//(p'-p)
	double t_proton_squared = t2_proton.Mag2();

	t2_dist->Fill( t_proton_squared );

	//T, momentum transfer from photon to Jpsi
	TLorentzVector T_Jpsi = particle_4mom_Jpsi - particle_4mom_photon;//delta
	double T_Jpsi_squared = T_Jpsi.Mag2();

	T_Jpsi_dist->Fill( T_Jpsi_squared );

	//compute COM s_NN of proton and neutron system:
	double E_NN = particle_4mom_proton.E() + particle_4mom_neutron.E();

	TVector3 p_proton = particle_4mom_proton.Vect();
	TVector3 p_neutron = particle_4mom_neutron.Vect();

	TVector3 total_NN = p_proton+p_neutron;
	TVector3 V_cm = (1/E_NN)*total_NN;

	double bx = V_cm.X();
	double by = V_cm.Y();
	double bz = V_cm.Z();

	TVector3 b;
	particle_4mom_proton.Boost(-bx,-by,-bz);
	particle_4mom_proton.Boost(b);

	particle_4mom_neutron.Boost(-bx,-by,-bz);
	particle_4mom_neutron.Boost(b);

	particle_4mom = particle_4mom_proton + particle_4mom_neutron;

	double sNN = particle_4mom.E()*particle_4mom.E();//center of mass energy squared
	 
	E_CM->Fill( sqrt(sNN) );

	//end COM

	//hadron angle vs their center of mass energy
	AngleVssNN_proton->Fill( sNN, theta_proton);
	AngleVssNN_neutron->Fill( sNN, theta_neutron);

	//Qsquared vs Jpsi transverse momentum
	Q2VsJpsi->Fill( pt_jpsi, trueQ2);

	//Wsqured vs Jpsi transverse momentum
	W2VsJpsi->Fill( pt_jpsi, trueW2);

	//t vs pt^2-Q^2
	T_hatVsPt2->Fill( pt_jpsi*pt_jpsi-trueQ2, T_Jpsi_squared );//should this be t_hat

	//T and t vs Jpsi pt
	TvsPt->Fill( pt_jpsi, T_Jpsi_squared );
	tProtonVsPt->Fill( pt_jpsi, t_proton_squared );

	//T and t vs center of mass energy
	ThatVssNN->Fill(sNN, T_Jpsi_squared);
	tdisVssNN->Fill(sNN, t_proton_squared);

	//t vs T
	tVsT->Fill(T_Jpsi_squared, t_proton_squared);

	//center of mass energy vs Jpsi pt
	sNNvsPt->Fill( pt_jpsi, sNN );

	//proton pt vs Jpsi pt
	PtVsPt_protonVsJpsi->Fill(pt_jpsi, pt_proton);

	//multiplicity distribution
	Ntrk_process_all->Fill(nParticles);
	Ntrk_process->Fill(nParticles_process);

	//double energy_diff = total4Mom_incoming.E() - total4Mom_outgoing.E();
	double energy_diff = 0;
	energy_corr->Fill( energy_diff );



   } // for



   TString outfilename;
   if( doKick ) outfilename = "_Jpsinodecay_KickFinalStates_eD_kick.root";
   else outfilename = "_Jpsinodecay_KickFinalStates_eD_nokick.root";

   TFile output("../rootfiles/"+inputFilename+outfilename,"RECREATE");
   
   //TH1D
   Ntrk_process_all->Write();
   Ntrk_process->Write();
   statusHist.Write();

   E_CM->Write();
   W2->Write();
   photonFlux->Write();
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

   //TH2D
   Q2VsX->Write();
   W2VsFlux->Write();
  
   T_hatVsPt2->Write();
   TvsPt->Write();
   tProtonVsPt->Write();

   ThatVssNN->Write();
   tdisVssNN->Write();
   tVsT->Write();

   sNNvsPt->Write();
   Q2VsJpsi->Write();
   W2VsJpsi->Write();

   PtVsEta_proton->Write();
   PtVsEta_neutron->Write();

   PtVsPt_protonVsJpsi->Write();

   AngleVsMom_proton->Write();
   AngleVsMom_neutron->Write();

   AngleVssNN_proton->Write();
   AngleVssNN_neutron->Write();

   energy_corr->Write();


}