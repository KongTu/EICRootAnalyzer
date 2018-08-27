#include "hist.h"//define all the histograms

TH1D* t1_dist = new TH1D("t1_dist",";t1", 200,-5,5);
TH1D* t2_dist = new TH1D("t2_dist",";t2", 200,-5,5);
TH1D* energy_corr = new TH1D("energy_corr",";E_{in} - E_{out}",1000,-0.05,0.05);

TH1D* sNN_dist = new TH1D("sNN_dist","s_{_{NN}} ",300,0,30);
TH2D* deltaEtadeltaPhi = new TH2D("deltaEtadeltaPhi",";#eta;#phi",200,-20,20,30,-7,7);
TH1D* px_dist = new TH1D("px_dist",";px",1000,-10,10);
TH1D* py_dist = new TH1D("py_dist",";py",1000,-10,10);

void run_KickFinalStates_reduced( int nEvents, bool doKick, TString inputFilename ) {
   
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

         if( index == 4 ){ //get gamma 4-momentum:

            particle_4mom_photon = particle->Get4Vector(); 
         }
         if( index == 3 ){
            particle_4mom_electron_prime = particle->Get4Vector();
         }
         if( status != 1 ) continue; //only stable final-state particles 
         if( pdg == 443 ){//Jpsi
           
            particle_4mom_jpsi_bKick = particle->Get4Vector();
            particle_4mom_jpsi = particle->Get4Vector();
         }
         if( pdg == 2212 ){//proton

         	//SetPtEtaPhiM(pt,eta,phi,mass);//this won't work if there is no pT 
            particle_4mom_proton_bKick = particle->Get4Vector();
            particle_4mom_proton = particle->Get4Vector();
         }
         if( pdg == 2112 ){//neutron
            
            particle_4mom_neutron_bKick = particle->Get4Vector();
            particle_4mom_neutron = particle->Get4Vector();
         }

         nParticles_process++;

	} // end of particle loop

	if( nParticles_process != 4 ) continue;
   
   if( doKick ){ 

      vector<TLorentzVector> afterKick = kickit(particle_4mom_neutron_bKick,particle_4mom_proton_bKick,particle_4mom_jpsi_bKick);
      particle_4mom_proton = afterKick[0];
      particle_4mom_neutron = afterKick[1];
      particle_4mom_jpsi = afterKick[2];

   }//end of kick

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
	double t_proton_squared = t2_proton.Mag2();

	t2_dist->Fill( t_proton_squared );

	//T, momentum transfer from photon to Jpsi
	TLorentzVector T_Jpsi = particle_4mom_jpsi - particle_4mom_photon;//delta
	double T_Jpsi_squared = T_Jpsi.Mag2();

	particle_4mom = particle_4mom_proton + particle_4mom_neutron;

	double sNN = particle_4mom.Mag2();//center of mass energy squared
	 
	E_CM->Fill( sqrt(sNN) );
	//end COM
	//if( fabs(t_proton_squared)  > 0.5 && fabs(t_proton_squared) < 5.0 && fabs(T_Jpsi_squared) < 0.5 ) 
      sNN_dist->Fill( sNN );

	//t vs T
	tVsT->Fill(T_Jpsi_squared, t_proton_squared);

   } // for


   TString outfilename;
   if( doKick ) outfilename = "_Jpsinodecay_KickFinalStates_eD_kick.root";
   else outfilename = "_Jpsinodecay_KickFinalStates_eD_nokick.root";

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


}