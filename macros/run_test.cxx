
// read.cxx
//
// Created by TB on 6/13/11.
// Copyright 2011 BNL. All rights reserved.
//
// Example of how to read a file produced by BuildTree for a simple analysis.
// To run, in ROOT do:
// root [0] .L /path/to/read.cxx
// root [1] read("myInputFile.root", 10000 )

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TBranchElement.h"
#include "TChain.h"

void run_test(int nEvents ) {
   
   // If the analysis solely uses TTree::Draw statements, you don't need to load
   // the shared library. You will receive warnings such as
   // Warning in <TClass::TClass>: no dictionary for class Particle is available
   // but these can be ignored. However if you wish to work with the event
   // objects themselves, the shared library must be loaded:
   gSystem->Load("/afs/rhic.bnl.gov/eic/MACROS/BuildTree/lib/BuildTree.so" ); // To use the master version
   
   //gSystem->Load("../lib/BuildTree.so" ); // To use your own compiled version
   
   // The TTrees are named EICTree.
   // Create a TChain for trees with this name.
   TChain *tree = new TChain("EICTree");
   
   // Add the file(s) we want to analyse to the chain.
   // We could add multiple files if we wanted.
   tree->Add("../../EICTree/eD_Jpsidiffnodecay_EICTree/eD_18x135_Q2_1_10_y_0.01_0.95_tau_7_noquench_kt=ptfrag=0.32_Shd1_ShdFac=1.32_Jpsidiffnodecay_test40k_wrongpf.root" ); // Wild cards are allowed e.g. tree.Add("*.root" );
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
   EventBeagle* event_beagle(NULL);
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

   // Now we can do some analysis...
   
   // Histograms for our analysis.
   TH1D ptHist("ptHist", "pT of all particles", 500, 0.0, 10 );
   TH1D statusHist("statusHist", "status distribution  ", 50, 0, 50 );
   
   TH1D* Ntrk = new TH1D("Ntrk",";Ntrk", 100, 0, 100);
   TH2D* pTvsThat = new TH2D("pTvsThat",";pT;t_hat", 1000,0,10,1000,-10,10);
   
   TH1D* energy_corr = new TH1D("energy_corr",";E_{in} - E_{out}",600,-30,30);
   TH1D* px_corr = new TH1D("px_corr","; px_{in} - px_{out}", 600,-30,30);
   TH1D* py_corr = new TH1D("py_corr","; py_{in} - py_{out}", 600,-30,30);
   TH1D* pz_corr = new TH1D("pz_corr","; pz_{in} - pz_{out}", 600,-30,30);

   TH1D* energy_corr_ionframe = new TH1D("energy_corr_ionframe",";E_{in} - E_{out}",600,-3000,3000);
   TH1D* pz_corr_ionframe = new TH1D("pz_corr_ionframe","; pz_{in} - pz_{out}", 600,-3000,3000);

   TH2D* energyVsQ2_2Dcorr = new TH2D("energyVsQ2_2Dcorr",";E_{in} - E_{out};Q2",600,-30,30, 200,0,10);
   TH2D* energyVsW2_2Dcorr = new TH2D("energyVsW2_2Dcorr",";E_{in} - E_{out};W2",600,-30,30, 2000,0,20000);
   TH2D* energyVsX_2Dcorr = new TH2D("energyVsX_2Dcorr",";E_{in} - E_{out};X",600,-30,30, 1000,0,0.5);
   TH2D* energyVsY_2Dcorr = new TH2D("energyVsY_2Dcorr",";E_{in} - E_{out};Y",600,-30,30, 1000,0,1.0);
   TH2D* energyVsNu_2Dcorr = new TH2D("energyVsNu_2Dcorr",";E_{in} - E_{out};Nu",600,-30,30, 200, 0,10000);
   TH2D* energyVsPf_2Dcorr = new TH2D("energyVsPf_2Dcorr",";E_{in} - E_{out};pf",600,-30,30, 200, -10,10);
   TH2D* energyVsPtf_2Dcorr = new TH2D("energyVsPtf_2Dcorr",";E_{in} - E_{out};ptf",600,-30,30, 200, -10,10);
   TH2D* energyVsProcess_2Dcorr = new TH2D("energyVsProcess_2Dcorr",";E_{in} - E_{out};process",600,-30,30, 10, 90,100);


   // Loop over events:
   for(int i(0); i < nEvents; ++i ) {
      
      // Read the next entry from the tree.
      tree->GetEntry(i);

      //fermi momentum in the ion rest frame with gamma* direction as z
      double pxf = branch_pxf->GetValue(0,0);
      double pyf = branch_pyf->GetValue(0,0);
      double pzf = branch_pzf->GetValue(0,0);
      double ptf = sqrt(pxf*pxf + pyf*pyf);
      double pf = sqrt(pxf*pxf + pyf*pyf + pzf*pzf);

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

      //Deuteron
      double pztarg = branch_pz->GetValue(0,0);
      double Atarg = branch_atarg->GetValue(0,0);
      double pz_total = pztarg*Atarg;
      double D_mass = 1.8755;//1.8755
      double total_energy = sqrt(pz_total*pz_total + D_mass*D_mass);

      //electron, neglect electron mass
      double pz_lepton = branch_pzlep->GetValue(0,0);
      double electron_mass = 0.00055;
      double total_lep_energy = sqrt(pz_lepton*pz_lepton + electron_mass*electron_mass);

      //Lab frame:
      TLorentzVector total4Mom_deuteron(0., 0., pz_total, total_energy);
      TLorentzVector total4Mom_electron(0., 0., pz_lepton, total_lep_energy);

      cout << "lab pz" << total4Mom_electron.Pz() << endl;
      total4Mom_electron.Boost(0.0,0.0,pz_total);
      //cout << "boost pz " << total4Mom_electron.Pz() << endl;

      //Lab frame
      TLorentzVector total4Mom_outgoing(0.,0.,0.,0.);
      TLorentzVector total4Mom_incoming = total4Mom_deuteron + total4Mom_electron;

      // We now know the number of particles in the event, so loop over
      // the particles:

     // The event contains a vector (array) of particles.
      int nParticles = event->GetNTracks();
 
      for(int j(0); j < nParticles; ++j ) {
         const erhic::ParticleMC* particle = event->GetTrack(j);
	 
	 // Let's just select charged pions for this example:
         int pdg = particle->GetPdgCode();
         int status = particle->GetStatus();
         int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
         
         // cout << "----- check if this is exchanged photon ------- " << endl;
         // cout << "index = " << index << " pdg = " << pdg << endl;
         // TLorentzVector photon_4mom = particle->PxPyPzE();
         // cout << "gamma px = " << photon_4mom.Px() << endl;
         // cout << "gamma py = " << photon_4mom.Py() << endl;
         // cout << "gamma pz = " << photon_4mom.Pz() << endl;

         
         if( status == 1 ){
            
            TLorentzVector particle_4mom = particle->PxPyPzE();//lab frame
            total4Mom_outgoing += particle_4mom;
         }
            
	      ptHist.Fill(particle->GetPt());
         statusHist.Fill( status ); 

      } // for

      double particle_pt = sqrt(total4Mom_outgoing.Px()*total4Mom_outgoing.Px() + total4Mom_outgoing.Py()*total4Mom_outgoing.Py());
      
      //1D histogram:
      pTvsThat->Fill( particle_pt, t_hat );
      Ntrk->Fill(nParticles);

      double energy_diff = total4Mom_incoming.E() - total4Mom_outgoing.E();
      energy_corr->Fill( energy_diff );
      px_corr->Fill( total4Mom_incoming.Px() - total4Mom_outgoing.Px() );
      py_corr->Fill( total4Mom_incoming.Py() - total4Mom_outgoing.Py() );
      pz_corr->Fill( total4Mom_incoming.Pz() - total4Mom_outgoing.Pz() );

      //2D histogram:
      energyVsQ2_2Dcorr->Fill(energy_diff, trueQ2);
      energyVsW2_2Dcorr->Fill(energy_diff, trueW2);
      energyVsX_2Dcorr->Fill(energy_diff, trueX);
      energyVsY_2Dcorr->Fill(energy_diff, trueY);
      energyVsNu_2Dcorr->Fill(energy_diff, trueNu);
      energyVsPf_2Dcorr->Fill(energy_diff, pf);
      energyVsPtf_2Dcorr->Fill(energy_diff, ptf);
      energyVsProcess_2Dcorr->Fill(energy_diff, event_process);

   } // for

   TFile output("wrongpf_JpsiNodecay_eD.root","RECREATE");
   
   ptHist.Write();    
   statusHist.Write();
   pTvsThat->Write();
   Ntrk->Write();
   energy_corr->Write();
   px_corr->Write();
   py_corr->Write();
   pz_corr->Write();

   energy_corr_ionframe->Write();
   pz_corr_ionframe->Write();

   energyVsQ2_2Dcorr->Write();
   energyVsW2_2Dcorr->Write();
   energyVsX_2Dcorr->Write();
   energyVsY_2Dcorr->Write();
   energyVsNu_2Dcorr->Write();
   energyVsPf_2Dcorr->Write();
   energyVsPtf_2Dcorr->Write();
   energyVsProcess_2Dcorr->Write();

}

