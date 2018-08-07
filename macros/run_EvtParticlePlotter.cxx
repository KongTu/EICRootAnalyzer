
// read.cxx
//
// Created by TB on 6/13/11.
// Copyright 2011 BNL. All rights reserved.
//
// Example of how to read a file produced by BuildTree for a simple analysis.
// To run, in ROOT do:
// root [0] .L /path/to/read.cxx
// root [1] read("myInputFile.root", 10000 )
#include "hist.h"

void run_EvtParticlePlotter( int nEvents, bool doBoost, TString inputFilename ) {
   
   // If the analysis solely uses TTree::Draw statements, you don't need to load
   // the shared library. You will receive warnings such as
   // Warning in <TClass::TClass>: no dictionary for class Particle is available
   // but these can be ignored. However if you wish to work with the event
   // objects themselves, the shared library must be loaded:
   // gSystem->Load("/afs/rhic.bnl.gov/eic/MACROS/BuildTree/lib/BuildTree.so" ); // To use the master version
   
   //gSystem->Load("../lib/BuildTree.so" ); // To use your own compiled version
   
   // The TTrees are named EICTree.
   // Create a TChain for trees with this name.
   TChain *tree = new TChain("EICTree");
   
   // Add the file(s) we want to analyse to the chain.
   // We could add multiple files if we wanted.
   
   tree->Add("../../EICTree/eD_Jpsidiffnodecay_EICTree/eD_18x135_Q2_1_10_y_0.01_0.95_tau_7_noquench_kt=ptfrag=0.32_Shd1_ShdFac=1.32_Jpsidiffnodecay_test40k_"+inputFilename+".root" ); // Wild cards are allowed e.g. tree.Add("*.root" );
   
   // tree.Add(/path/to/otherFileNames ); // etc... 
   
   // Create an object to store the current event from the tree.
   // This is how we access the values in the tree.
   // If you want to use generator-specific values, then
   // the event type should match the type in the TTree. Valid types are
   // EventPythia, EventPepsi, EventRapgap, EventDjangoh, EventMilou.
   // If you only need shared quantities like x, Q2 and the particle list
   // you can use EventBase and the macro will be general for any Monte Carlo.
   
   EventPythia* event(NULL);// = new EventPythia;
   
   //EventBase* event(NULL);
   // EventBeagle* event_beagle(NULL);
   // Now associate the contents of the branch with the buffer.
   // The events are stored in a branch named event:
   tree->SetBranchAddress("event", &event ); // Note &event, not event.
   // tree.SetBranchAddress("event", &event_beagle ); // Note &event, not event.
      
   TBranchElement* branch_atarg = (TBranchElement*) tree->GetBranch("Atarg");
   TBranchElement* branch_pz = (TBranchElement*) tree->GetBranch("pztarg");
   TBranchElement* branch_pzlep = (TBranchElement*) tree->GetBranch("pzlep");
   TBranchElement* branch_pxf = (TBranchElement*) tree->GetBranch("pxf");
   TBranchElement* branch_pyf = (TBranchElement*) tree->GetBranch("pyf");
   TBranchElement* branch_pzf = (TBranchElement*) tree->GetBranch("pzf");

   // Histograms for our analysis.
   TH1D* Ntrk_process_91 = new TH1D("Ntrk_process_91",";Ntrk_process_91", 100, 0, 100);
   TH1D* Ntrk_process_93 = new TH1D("Ntrk_process_93",";Ntrk_process_93", 100, 0, 100);
   TH1D* Ntrk_process_all = new TH1D("Ntrk_process_all",";Ntrk_process_all", 100, 0, 100);
   TH1D statusHist("statusHist", "status distribution  ", 50, 0, 50 );
   
   TH1D* Ntrk = new TH1D("Ntrk",";Ntrk", 100, 0, 100);
   TH2D* pTvsThat = new TH2D("pTvsThat",";pT;t_hat", 1000,0,10,1000,-10,10);
   //more in hist.h
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
    
      // Let's just select charged pions for this example:
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

         if( status != 1 ) continue;
         if( event_process == 91 ){

            PtDist_process_91->Fill( pt );
            EtaDist_process_91->Fill( eta );
            PhiDist_process_91->Fill( phi );

            if( pdg == 443 ){//Jpsi

               PtDist_process_91_Jpsi->Fill( pt );
               EtaDist_process_91_Jpsi->Fill( eta );
               PhiDist_process_91_Jpsi->Fill( phi );

               double pt_91_jpsi = particle->GetPt();
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

      //compute center of mass energy of proton and neutron system:
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

      AngleVssNN_process_91_proton->Fill( sNN, theta_91_proton);
      AngleVssNN_process_93_proton->Fill( sNN, theta_93_proton);
      AngleVssNN_process_91_neutron->Fill( sNN, theta_91_neutron);
      AngleVssNN_process_93_neutron->Fill( sNN, theta_93_neutron);
      
      Q2VsJpsi_91->Fill( pt_91_jpsi, trueQ2);
      Q2VsJpsi_93->Fill( pt_93_jpsi, trueQ2);
      T_hat->Fill( pt_91_jpsi*pt_91_jpsi-trueQ2, t_hat );
      TvsPt->Fill( pt_91_jpsi, t_hat );
      
      PtVsPt_process_91_protonVsJpsi->Fill(pt_91_jpsi, pt_91_proton);
      PtVsPt_process_93_protonVsJpsi->Fill(pt_93_jpsi, pt_93_proton);//this may have more than one proton in one event

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
   T_hat->Write();
   TvsPt->Write();

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

