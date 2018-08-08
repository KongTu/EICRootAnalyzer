#include "hist.h"//define all the histograms

void run_EvtParticlePlotter( int nEvents, bool doBoost, TString inputFilename ) {
   
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

   // Now we can do some analysis...
   // Loop over events:
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
         T_dist->Fill( t_hat );

      double u_hat = event->GetHardU();
      double photon_flux = event->GetPhotonFlux();
         photonFlux->Fill( photon_flux );

      int event_process = event->GetProcess();

      //initial proton in Deuteron
      double pztarg = branch_pz->GetValue(0,0);
      double Atarg = branch_atarg->GetValue(0,0);
      double pz_total = pztarg*(Atarg-1.0);//proton longitudinal momentum
      double total_energy = sqrt(pz_total*pz_total + MASS_PROTON*MASS_PROTON);

      TLorentzVector total4Mom_iProton(0., 0., pz_total, total_energy);

      Q2VsX->Fill(trueX, trueQ2);
      W2VsFlux->Fill(photon_flux, trueW2);

      // The event contains a vector (array) of particles.
      int nParticles = event->GetNTracks();
      //event t_hat

      // We now know the number of particles in the event, so loop over
      // the particles:

      int nParticles_process = 0;

      TLorentzVector particle_4mom;
      TLorentzVector particle_4mom_proton;
      TLorentzVector particle_4mom_neutron;

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


         if( status != 1 ) continue; //only stable final-state particles 

         // cout << "pdg " << pdg << endl;
         // cout << "status " << status << endl;
         // cout << "index " << index << endl;
         // cout << "pt " << pt << endl;

            if( pdg == 443 ){//Jpsi

               PtDist_Jpsi->Fill( pt );
               EtaDist_Jpsi->Fill( eta );
               PhiDist_Jpsi->Fill( phi );

               double pt_jpsi = pt;//store Jpsi pt

            }
            if( pdg == 2212 ){//proton

               PtDist_proton->Fill( pt );
               EtaDist_proton->Fill( eta );
               PhiDist_proton->Fill( phi );

               PtVsEta_proton->Fill(eta, pt);
               AngleVsMom_proton->Fill(mom, theta);
               
               double theta_proton = theta;//store proton angle
               double pt_proton = pt;//store proton pt
               
               particle_4mom_proton = particle->PxPyPzE();//store proton 4vector

            }
            if( pdg == 2112 ){//neutron

               PtDist_neutron->Fill( pt );
               EtaDist_neutron->Fill( eta );
               PhiDist_neutron->Fill( phi );

               PtVsEta_neutron->Fill(eta, pt);
               AngleVsMom_neutron->Fill(mom, theta);
               
               double theta_neutron = theta;//store neutron angle

               particle_4mom_neutron = particle->PxPyPzE();

            }

            nParticles_process++;

      } // end of particle loop

      if( nParticles_process < 4 ) continue; 
      
      //small t, namely the momentum transfer to the struck nucleon (proton)
      TLorentzVector t_proton = particle_4mom_proton - total4Mom_iProton;//(p'-p)
      double t_proton_squared = t_proton.Mag2();

      if( t_proton_squared > 0.0 ) continue;
      
      // if( t_proton_squared > 0.0 ){

      //    cout << "event number " << i << endl;
      //    cout << "event process " << event_process << endl;
      //    cout << "nParticles_process " << nParticles_process << endl;
      //    cout << "proton energy " << particle_4mom_proton.E() << endl;
      //    cout << "neutron energy " << particle_4mom_neutron.E() << endl;
      // }

      t_dist->Fill( t_proton_squared );

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
      T_hatVsPt2->Fill( pt_jpsi*pt_jpsi-trueQ2, t_hat );
      
      //T and t vs Jpsi pt
      TvsPt->Fill( pt_jpsi, t_hat );
      tProtonVsPt->Fill( pt_jpsi, t_proton_squared );

      //T and t vs center of mass energy
      ThatVssNN->Fill(sNN, t_hat);
      tdisVssNN->Fill(sNN, t_proton_squared);
      
      //t vs T
      tVsT->Fill(t_hat, t_proton_squared);

      //center of mass energy vs Jpsi pt
      sNNvsPt->Fill( pt_jpsi, sNN );
      
      //proton pt vs Jpsi pt
      PtVsPt_protonVsJpsi->Fill(pt_jpsi, pt_proton);

      //multiplicity distribution
      Ntrk_process_all->Fill(nParticles);
      Ntrk_process->Fill(nParticles_process);

   } // for

   TString outfilename;
   if( doBoost ) outfilename = "_Jpsinodecay_EvtParticlePlotter_eD_ionframe.root";
   else outfilename = "_Jpsinodecay_EvtParticlePlotter_eD.root";

   TFile output("../rootfiles/"+inputFilename+outfilename,"RECREATE");
   
   //TH1D
   Ntrk_process_all->Write();
   Ntrk_process->Write();
   statusHist.Write();

   E_CM->Write();
   W2->Write();
   photonFlux->Write();
   T_dist->Write();//T_distribution in the selected range
   t_dist->Write();//t_distribution in the selected range

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
  

}

