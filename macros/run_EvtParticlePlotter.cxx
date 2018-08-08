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
      double trueX = event->GetTrueX();
      double trueY = event->GetTrueY();
      double trueNu = event->GetTrueNu();
      double s_hat = event->GetHardS();
      double t_hat = event->GetHardT();
      double u_hat = event->GetHardU();
      int event_process = event->GetProcess();

      //initial proton in Deuteron
      double pztarg = branch_pz->GetValue(0,0);
      double Atarg = branch_atarg->GetValue(0,0);
      double pz_total = pztarg*(Atarg-1.0);//proton longitudinal momentum
      double P_mass = 0.93827;//1.8624778138724238 for D, 193.69769264273208 for Pb
      double total_energy = sqrt(pz_total*pz_total + P_mass*P_mass);

      TLorentzVector total4Mom_iProton(0., 0., pz_total, total_energy);

      Q2VsX->Fill(trueX, trueQ2);

      // The event contains a vector (array) of particles.
      int nParticles = event->GetNTracks();
      //event t_hat

      // We now know the number of particles in the event, so loop over
      // the particles:

      int nParticles_process_91 = 0;
      int nParticles_process_93 = 0;

      TLorentzVector particle_4mom;
      TLorentzVector particle_4mom_proton;
      TLorentzVector particle_4mom_neutron;

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
         
         if( event_process == 91 ){

            PtDist_process_91->Fill( pt );
            EtaDist_process_91->Fill( eta );
            PhiDist_process_91->Fill( phi );

            if( pdg == 443 ){//Jpsi

               PtDist_process_91_Jpsi->Fill( pt );
               EtaDist_process_91_Jpsi->Fill( eta );
               PhiDist_process_91_Jpsi->Fill( phi );

               double pt_91_jpsi = particle->GetPt();

               if( fabs(eta) < 2.0 && fabs(t_hat) < 1.0 ){
                  
                  double pt_91_jpsi_differential = particle->GetPt();
               }

            }
            if( pdg == 2212 ){//proton

               PtDist_process_91_proton->Fill( pt );
               EtaDist_process_91_proton->Fill( eta );
               PhiDist_process_91_proton->Fill( phi );

               PtVsEta_process_91_proton->Fill(eta, pt);
               AngleVsMom_process_91_proton->Fill(mom, theta);
               
               double theta_91_proton = theta;//for later use
               double pt_91_proton = particle->GetPt();
               particle_4mom_proton = particle->PxPyPzE();

            }
            if( pdg == 2112 ){//neutron

               PtDist_process_91_neutron->Fill( pt );
               EtaDist_process_91_neutron->Fill( eta );
               PhiDist_process_91_neutron->Fill( phi );

               PtVsEta_process_91_neutron->Fill(eta, pt);
               AngleVsMom_process_91_neutron->Fill(mom, theta);
               double theta_91_neutron = theta;//for later use

               particle_4mom_neutron = particle->PxPyPzE();

            }

            nParticles_process_91++;

         }
         if( event_process == 93 ){

            PtDist_process_93->Fill( pt );
            EtaDist_process_93->Fill( eta );
            PhiDist_process_93->Fill( phi );

            if( pdg == 443 ){//Jpsi

               PtDist_process_93_Jpsi->Fill( pt );
               EtaDist_process_93_Jpsi->Fill( eta );
               PhiDist_process_93_Jpsi->Fill( phi );

               double pt_93_jpsi = particle->GetPt();
            }
            if( pdg == 2212 ){//proton

               PtDist_process_93_proton->Fill( pt );
               EtaDist_process_93_proton->Fill( eta );
               PhiDist_process_93_proton->Fill( phi );

               PtVsEta_process_93_proton->Fill(eta, pt);
               AngleVsMom_process_93_proton->Fill(mom, theta);

               double theta_93_proton = theta;
               double pt_93_proton = particle->GetPt();
            }
            if( pdg == 2112 ){//neutron

               PtDist_process_93_neutron->Fill( pt );
               EtaDist_process_93_neutron->Fill( eta );
               PhiDist_process_93_neutron->Fill( phi );

               PtVsEta_process_93_neutron->Fill(eta, pt);
               AngleVsMom_process_93_neutron->Fill(mom, theta);
               double theta_93_neutron = theta;//for later use
            }

            nParticles_process_93++;

         } 

      } // for

      //small t, namely the momentum transfer to the struck nucleon (proton)
      TLorentzVector t_proton = particle_4mom_proton - total4Mom_iProton;//(p'-p)
      double t_proton_squared = t_proton.Mag2();

         t_proton_dist->Fill( t_proton_squared );

      //compute COM s_NN of proton and neutron system:
      double E_NN = particle_4mom_proton.E() + particle_4mom_neutron.E();
      double Pt_proton = particle_4mom_proton.Pt();

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

      double sNN = particle_4mom.E();//center of mass energy
         
         E_CM->Fill( sNN );
      
      //end COM

      //hadron angle vs their center of mass energy
      AngleVssNN_process_91_proton->Fill( sNN, theta_91_proton);
      AngleVssNN_process_93_proton->Fill( sNN, theta_93_proton);
      AngleVssNN_process_91_neutron->Fill( sNN, theta_91_neutron);
      AngleVssNN_process_93_neutron->Fill( sNN, theta_93_neutron);
      
      //Qsquared vs Jpsi transverse momentum
      Q2VsJpsi_91->Fill( pt_91_jpsi, trueQ2);
      Q2VsJpsi_93->Fill( pt_93_jpsi, trueQ2);
      
      //t vs pt^2-Q^2
      T_hatVsPt2->Fill( pt_91_jpsi*pt_91_jpsi-trueQ2, t_hat );
      
      //t vs Jpsi pt
      TvsPt_91->Fill( pt_91_jpsi, t_hat );
      TvsPt_93->Fill( pt_93_jpsi, t_hat );

      //sNN vs Jpsi pt
      if( fabs(t_hat) < 1.0 && fabs(t_proton_squared) > 1.0 ) {

         T_dist->Fill( t_hat );
         t_dist->Fill( t_proton_squared );
         sNNvsPt_91->Fill( pt_91_jpsi_differential, sNN*sNN );
      }
      //Jpsi pt vs proton pt
      PtVsPt_process_91_protonVsJpsi->Fill(pt_91_jpsi, pt_91_proton);
      PtVsPt_process_93_protonVsJpsi->Fill(pt_93_jpsi, pt_93_proton);//this may have more than one proton in one event

      //multiplicity distribution
      Ntrk_process_all->Fill(nParticles);
      Ntrk_process_91->Fill(nParticles_process_91);
      Ntrk_process_93->Fill(nParticles_process_93);

   } // for

   TString outfilename;
   if( doBoost ) outfilename = "_Jpsinodecay_EvtParticlePlotter_eD_ionframe.root";
   else outfilename = "_Jpsinodecay_EvtParticlePlotter_eD.root";

   TFile output("../rootfiles/"+inputFilename+outfilename,"RECREATE");
   
   Q2VsX->Write();
   
   E_CM->Write();
   T_dist->Write();//T_distribution in the selected range
   t_dist->Write();//t_distribution in the selected range
   t_proton_dist->Write();//t_distribution for proton for entire range

   T_hatVsPt2->Write();
   TvsPt_91->Write();
   TvsPt_93->Write();
   sNNvsPt_91->Write();

   Q2VsJpsi_91->Write();
   Q2VsJpsi_93->Write();

   statusHist.Write();

   Ntrk_process_all->Write();
   Ntrk_process_91->Write();
   Ntrk_process_93->Write();

   PtDist_process_91->Write();
   EtaDist_process_91->Write();
   PhiDist_process_91->Write();

   PtDist_process_91_Jpsi->Write();
   EtaDist_process_91_Jpsi->Write();
   PhiDist_process_91_Jpsi->Write();

   PtDist_process_91_proton->Write();
   EtaDist_process_91_proton->Write();
   PhiDist_process_91_proton->Write();

   PtDist_process_91_neutron->Write();
   EtaDist_process_91_neutron->Write();
   PhiDist_process_91_neutron->Write();

   PtDist_process_93->Write();
   EtaDist_process_93->Write();
   PhiDist_process_93->Write();

   PtDist_process_93_Jpsi->Write();
   EtaDist_process_93_Jpsi->Write();
   PhiDist_process_93_Jpsi->Write();

   PtDist_process_93_proton->Write();
   EtaDist_process_93_proton->Write();
   PhiDist_process_93_proton->Write();

   PtDist_process_93_neutron->Write();
   EtaDist_process_93_neutron->Write();
   PhiDist_process_93_neutron->Write();

   PtVsEta_process_91_proton->Write();
   PtVsEta_process_93_proton->Write();
   PtVsEta_process_91_neutron->Write();
   PtVsEta_process_93_neutron->Write();

   PtVsPt_process_91_protonVsJpsi->Write();
   PtVsPt_process_93_protonVsJpsi->Write();

   AngleVsMom_process_91_proton->Write();
   AngleVsMom_process_93_proton->Write();
   AngleVsMom_process_91_neutron->Write();
   AngleVsMom_process_93_neutron->Write();

   AngleVssNN_process_91_proton->Write();
   AngleVssNN_process_93_proton->Write();
   AngleVssNN_process_91_neutron->Write();
   AngleVssNN_process_93_neutron->Write();
  

}

