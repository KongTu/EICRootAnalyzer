
// read.cxx
//
// Created by TB on 6/13/11.
// Copyright 2011 BNL. All rights reserved.
//
// Example of how to read a file produced by BuildTree for a simple analysis.
// To run, in ROOT do:
// root [0] .L /path/to/read.cxx
// root [1] read("myInputFile.root", 10000 )

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

   // Now we can do some analysis...
   
   // Histograms for our analysis.
   TH1D* Ntrk_process_91 = new TH1D("Ntrk_process_91",";Ntrk_process_91", 100, 0, 100);
   TH1D* Ntrk_process_93 = new TH1D("Ntrk_process_93",";Ntrk_process_93", 100, 0, 100);
   TH1D* Ntrk_process_all = new TH1D("Ntrk_process_all",";Ntrk_process_all", 100, 0, 100);
   TH1D statusHist("statusHist", "status distribution  ", 50, 0, 50 );
   
   TH1D* Ntrk = new TH1D("Ntrk",";Ntrk", 100, 0, 100);
   TH2D* pTvsThat = new TH2D("pTvsThat",";pT;t_hat", 1000,0,10,1000,-10,10);
   
   TH1D* PtDist_process_91 = new TH1D("PtDist_process_91",";PtDist_process_91", 200, 0,10);
   TH1D* PhiDist_process_91 = new TH1D("PhiDist_process_91",";PhiDist_process_91", 2000, 0,10);
   TH1D* EtaDist_process_91 = new TH1D("EtaDist_process_91",";EtaDist_process_91", 2000, -10,10);

   TH1D* PtDist_process_93 = new TH1D("PtDist_process_93",";PtDist_process_93", 200, 0,10);
   TH1D* PhiDist_process_93 = new TH1D("PhiDist_process_93",";PhiDist_process_93", 2000, 0,10);
   TH1D* EtaDist_process_93 = new TH1D("EtaDist_process_93",";EtaDist_process_93", 2000, -10,10);

   TH1D* PtDist_process_91_Jpsi = new TH1D("PtDist_process_91_Jpsi",";PtDist_process_91_Jpsi", 200, 0,10);
   TH1D* PhiDist_process_91_Jpsi = new TH1D("PhiDist_process_91_Jpsi",";PhiDist_process_91_Jpsi", 2000, 0,10);
   TH1D* EtaDist_process_91_Jpsi = new TH1D("EtaDist_process_91_Jpsi",";EtaDist_process_91_Jpsi", 2000, -10,10);

   TH1D* PtDist_process_93_Jpsi = new TH1D("PtDist_process_93_Jpsi",";PtDist_process_93_Jpsi", 200, 0,10);
   TH1D* PhiDist_process_93_Jpsi = new TH1D("PhiDist_process_93_Jpsi",";PhiDist_process_93_Jpsi", 2000, 0,10);
   TH1D* EtaDist_process_93_Jpsi = new TH1D("EtaDist_process_93_Jpsi",";EtaDist_process_93_Jpsi", 2000, -10,10);

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


      // The event contains a vector (array) of particles.
      int nParticles = event->GetNTracks();
      //event t_hat

      // We now know the number of particles in the event, so loop over
      // the particles:

      int nParticles_process_91 = 0;
      int nParticles_process_93 = 0;

      for(int j(0); j < nParticles; ++j ) {
         
         const erhic::ParticleMC* particle = event->GetTrack(j);
    
      // Let's just select charged pions for this example:
         int pdg = particle->GetPdgCode();
         int status = particle->GetStatus();
         int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
         double pt = particle->GetPt();
         double eta = particle->GetEta();
         double phi = particle->GetPhi();

         statusHist.Fill( status ); 

         if( status != 1 ) continue;
         if( event_process == 91 ){

            PtDist_process_91->Fill( pt );
            EtaDist_process_91->Fill( eta );
            PhiDist_process_91->Fill( phi );

            if( pdg == 443 ){

               PtDist_process_91_Jpsi->Fill( pt );
               EtaDist_process_91_Jpsi->Fill( eta );
               PhiDist_process_91_Jpsi->Fill( phi );

            }
            nParticles_process_91++;
         }
         if( event_process == 93 ){

            PtDist_process_93->Fill( pt );
            EtaDist_process_93->Fill( eta );
            PhiDist_process_93->Fill( phi );

            if( pdg == 443 ){

               PtDist_process_93_Jpsi->Fill( pt );
               EtaDist_process_93_Jpsi->Fill( eta );
               PhiDist_process_93_Jpsi->Fill( phi );

            }
            nParticles_process_93++;
         } 


      } // for

      Ntrk_process_all->Fill(nParticles);
      Ntrk_process_91->Fill(nParticles_process_91);
      Ntrk_process_93->Fill(nParticles_process_93);

   } // for

   TString outfilename;
   if( doBoost ) outfilename = "_Jpsinodecay_EvtParticlePlotter_eD_ionframe.root";
   else outfilename = "_Jpsinodecay_EvtParticlePlotter_eD.root";

   TFile output("../rootfiles/"+inputFilename+outfilename,"RECREATE");
   
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

   PtDist_process_93->Write();
   EtaDist_process_93->Write();
   PhiDist_process_93->Write();

   PtDist_process_93_Jpsi->Write();
   EtaDist_process_93_Jpsi->Write();
   PhiDist_process_93_Jpsi->Write();

}

