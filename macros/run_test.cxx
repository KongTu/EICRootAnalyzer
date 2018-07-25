
// read.cxx
//
// Created by TB on 6/13/11.
// Copyright 2011 BNL. All rights reserved.
//
// Example of how to read a file produced by BuildTree for a simple analysis.
// To run, in ROOT do:
// root [0] .L /path/to/read.cxx
// root [1] read("myInputFile.root", 10000 )

void run_test(TString inFileNames, int nEvents ) {
   
   // If the analysis solely uses TTree::Draw statements, you don't need to load
   // the shared library. You will receive warnings such as
   // Warning in <TClass::TClass>: no dictionary for class Particle is available
   // but these can be ignored. However if you wish to work with the event
   // objects themselves, the shared library must be loaded:
   gSystem->Load("/afs/rhic.bnl.gov/eic/MACROS/BuildTree/lib/BuildTree.so" ); // To use the master version
   
   //gSystem->Load("../lib/BuildTree.so" ); // To use your own compiled version
   
   // The TTrees are named EICTree.
   // Create a TChain for trees with this name.
   TChain tree("EICTree" );
   
   // Add the file(s) we want to analyse to the chain.
   // We could add multiple files if we wanted.
   tree.Add(inFileNames ); // Wild cards are allowed e.g. tree.Add("*.root" );
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
   
   // Now associate the contents of the branch with the buffer.
   // The events are stored in a branch named event:
   tree.SetBranchAddress("event", &event ); // Note &event, not event.
      
   TBranchElement* branch = (TBranchElement*) tree.GetBranch("Atarg");
   TBranchElement* branch_pz = (TBranchElement*) tree.GetBranch("pztarg");
   TBranchElement* branch_pzlep = (TBranchElement*) tree.GetBranch("pzlep");


   // Now we can do some analysis...
   
   // We record the largest particle pT we find here:
   //double highestPt(-1. );
   
   // Histograms for our analysis.
   TH1D ptHist("ptHist", "pT of charged pions", 500, 0.0, 10 );
   TH1D statusHist("statusHist", "status distribution  ", 50, 0, 50 );
   
   TH1D* Ntrk = new TH1D("Ntrk",";Ntrk", 100, 0, 100);
   TH2D* pTvsThat = new TH2D("pTvsThat",";pT;t_hat", 1000,0,10,1000,-10,10);
   TH2D* energy_corr = new TH2D("energy_corr",";incoming;outgoing", 100,0,1000,100,0,1000);
   TH2D* px_corr = new TH2D("px_corr",";incoming;outgoing", 100,0,1000,100,0,1000);
   TH2D* py_corr = new TH2D("py_corr",";incoming;outgoing", 100,0,1000,100,0,1000);
   TH2D* pz_corr = new TH2D("pz_corr",";incoming;outgoing", 100,0,1000,100,0,1000);

   // Loop over events:
   for(int i(0); i < nEvents; ++i ) {
      
      // Read the next entry from the tree.
      tree.GetEntry(i);

      //Deuteron
      double pztarg = branch_pz->GetValue(0,0);
      double Atarg = branch->GetValue(0,0);
      double pz_total = pztarg*Atarg;
      double D_mass = 1.8755;
      double total_energy = sqrt(pz_total*pz_total + D_mass*D_mass);

      //electron, neglect electron mass
      int pz_lepton = branch_pzlep->GetValue(0,0);
      double total_lep_energy = sqrt(pz_lepton*pz_lepton);


      TLorentzVector total4Mom_deuteron(total_energy, 0., 0., pz_total);
      TLorentzVector total4Mom_electron(total_lep_energy, 0., 0., pz_lepton);

      // The event contains a vector (array) of particles.
      int nParticles = event->GetNTracks();
      //event t_hat
      double t_hat = event->GetHardT();
      
      TLorentzVector total4Mom_outgoing(0.,0.,0.,0.);
      TLorentzVector total4Mom_incoming = total4Mom_deuteron + total4Mom_electron;



      cout << "--------- event " << i << "------- " << endl;

      cout << "energy " << total4Mom_incoming.E() << endl;
      cout << "pz " << total4Mom_incoming.Pz() << endl;


      // We now know the number of particles in the event, so loop over
      // the particles:
      for(int j(0); j < nParticles; ++j ) {
         const erhic::ParticleMC* particle = event->GetTrack(j);
	 
	 // Let's just select charged pions for this example:
         int pdg = particle->GetPdgCode();
         int status = particle->GetStatus();
         int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.

         TLorentzVector particle_4mom = particle->PxPyPzE();

         if( status == 1 ){
            total4Mom_outgoing += particle_4mom;
         }
            
	      ptHist.Fill(particle->GetPt());
         statusHist.Fill( status ); 

      } // for

      cout << "---------------- " << endl;


      double particle_pt = sqrt(total4Mom_outgoing.Px()*total4Mom_outgoing.Px() + total4Mom_outgoing.Py()*total4Mom_outgoing.Py());
      pTvsThat->Fill( particle_pt, t_hat );
      Ntrk->Fill(nParticles);
      energy_corr->Fill( total4Mom_incoming.E(), total4Mom_outgoing.E() );
      px_corr->Fill( total4Mom_incoming.Px(), total4Mom_outgoing.Px() );
      py_corr->Fill( total4Mom_incoming.Py(), total4Mom_outgoing.Py() );
      pz_corr->Fill( total4Mom_incoming.Pz(), total4Mom_outgoing.Pz() );

      	
   } // for

   TFile output("test.root","RECREATE");
   ptHist.Write();    
   statusHist.Write();
   pTvsThat->Write();
   Ntrk->Write();
   energy_corr->Write();
   px_corr->Write();
   py_corr->Write();
   pz_corr->Write();

}

