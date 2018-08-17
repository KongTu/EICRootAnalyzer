#define PI            3.1415926

#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_DEUTERON 1.87783999999999995
#define MASS_TRITON   2.7937167208086358
#define MASS_HE3      2.7937167208086358
#define MASS_ALPHA    3.7249556277448477
#define MASS_LI6      5.5874334416172715
#define MASS_C12      11.174866883234543
#define MASS_CA40     37.249556277448477
#define MASS_XE131    121.99229680864376
#define MASS_AU197    183.45406466643374
#define MASS_PB208    193.69769264273208

void run_DMPJet3( int nEvents, bool doBoost, TString inputFilename ) {
   
   TChain *tree = new TChain("EICTree");
   
   tree->Add("../../EICTree/eD_DPMjet3/eD_18x135_y_0.01_0.95_Q2_1_20_1kevents.root" ); // Wild cards are allowed e.g. tree.Add("*.root" );

   EventDpmjet* event(NULL);// = new EventPythia;
   //EventBase* event(NULL);
   //EventBeagle* event_beagle(NULL);

   tree->SetBranchAddress("event", &event ); // Note &event, not event.
   // tree.SetBranchAddress("event", &event_beagle ); // Note &event, not event.
   
   TBranchElement* branch_process1 = (TBranchElement*) tree->GetBranch("process1");
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

   TH2D* energyVsQ2_2Dcorr = new TH2D("energyVsQ2_2Dcorr",";E_{in} - E_{out};Q2",600,-30,30, 200,0,10);
   TH2D* energyVsW2_2Dcorr = new TH2D("energyVsW2_2Dcorr",";E_{in} - E_{out};W2",600,-30,30, 2000,0,20000);
   TH2D* energyVsX_2Dcorr = new TH2D("energyVsX_2Dcorr",";E_{in} - E_{out};X",600,-30,30, 1000,0,0.5);
   TH2D* energyVsY_2Dcorr = new TH2D("energyVsY_2Dcorr",";E_{in} - E_{out};Y",600,-30,30, 1000,0,1.0);
   TH2D* energyVsNu_2Dcorr = new TH2D("energyVsNu_2Dcorr",";E_{in} - E_{out};Nu",600,-30,30, 200, 0,10000);
   TH2D* energyVsPf_2Dcorr = new TH2D("energyVsPf_2Dcorr",";E_{in} - E_{out};pf",600,-30,30, 2000, -1,1);
   TH2D* energyVsPtf_2Dcorr = new TH2D("energyVsPtf_2Dcorr",";E_{in} - E_{out};ptf",600,-30,30, 2000, -1,1);
   TH2D* energyVsProcess_2Dcorr = new TH2D("energyVsProcess_2Dcorr",";E_{in} - E_{out};process",600,-30,30, 10, 90,100);


   // Loop over events:
   for(int i(0); i < nEvents; ++i ) {
      
      // Read the next entry from the tree.
      tree->GetEntry(i);

      //fermi momentum in the ion rest frame with gamma* direction as z
      // double pxf = branch_pxf->GetValue(0,0);
      // double pyf = branch_pyf->GetValue(0,0);
      // double pzf = branch_pzf->GetValue(0,0);
      // double ptf = sqrt(pxf*pxf + pyf*pyf);
      // double pf = sqrt(pxf*pxf + pyf*pyf + pzf*pzf);

      //event information:
      double trueQ2 = event->GetQ2();
      double trueW2 = event->GetW2();
      double trueX = event->GetX();
      double trueY = event->GetY();
      double trueNu = event->GetNu();
    
      int event_process = branch_process1->GetValue(0,0);

      cout << "process " << event_process << endl;
      if( event_process != 5 ) continue;

      //Deuteron hard code for DMPJet
      double pztarg = 135.0;
      double Atarg = 2.0;
      double pz_total = pztarg*Atarg;
      double total_energy = sqrt(pz_total*pz_total + MASS_DEUTERON*MASS_DEUTERON);

      //electron, neglect electron mass
      double pz_lepton = -17.99482;
      double electron_mass = 0.00051;
      double total_lep_energy = 17.99481964;//sqrt(pz_lepton*pz_lepton + electron_mass*electron_mass);

      TLorentzVector total4Mom_deuteron(0., 0., pz_total, total_energy);
      TLorentzVector total4Mom_electron(0., 0., pz_lepton, total_lep_energy);

      /* lorentz boost incoming particle*/
      double gamma_ion = total_energy/MASS_DEUTERON;
      double bz = pz_total/(gamma_ion*MASS_DEUTERON);
      // double gamma_ion = sqrt(135.*135. + MASS_NEUTRON*MASS_NEUTRON)/MASS_NEUTRON;
      // double bz = 135./(gamma_ion*MASS_NEUTRON);

      cout << "gamma factor: " << gamma_ion << endl;

      TVector3 b;

      if( doBoost ){
         total4Mom_electron.Boost(0,0,-bz);
         total4Mom_electron.Boost(b);

         total4Mom_deuteron.Boost(0,0,-bz);
         total4Mom_deuteron.Boost(b);
      }
      //end here

      // The event contains a vector (array) of particles.
      int nParticles = event->GetNTracks();
      //event t_hat

      // cout << "----------- event " << i << "--------------" << endl;
      // cout << "electron pz " << total4Mom_electron.Pz() << endl;
      // cout << "deuteron pz " << total4Mom_deuteron.Pz() << endl;


      TLorentzVector total4Mom_outgoing(0.,0.,0.,0.);
      TLorentzVector total4Mom_incoming = total4Mom_deuteron + total4Mom_electron;

      TLorentzVector test_proton(-0.01865,-0.0203,0.00038,0.93867);
      test_proton.Boost(0,0,bz);
      test_proton.Boost(b);

      cout << "fermi energy in lab frame " << test_proton.E() << endl;
      cout << "fermi pz in lab frame " << test_proton.Pz() << endl; 

      // cout << "incoming pz " << total4Mom_incoming.Pz() << endl;
      // cout << "incoming E " << total4Mom_incoming.E() << endl;

      // We now know the number of particles in the event, so loop over
      // the particles:
      for(int j(0); j < nParticles; ++j ) {
         const erhic::ParticleMC* particle = event->GetTrack(j);
    
    // Let's just select charged pions for this example:
         int pdg = particle->GetPdgCode();
         int status = particle->GetStatus();
         int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
 
         TLorentzVector particle_4mom = particle->PxPyPzE();
         
         if(doBoost){
               particle_4mom.Boost(0,0,-bz);
               particle_4mom.Boost(b);

         }

         if( index == 1 ) {

            cout << "pdg " << pdg << endl;
            cout << "status" << status << endl;
            cout << "index " << index << endl;
            cout << "mass " << particle->GetM() << endl;
            cout << std::setprecision(10) << "particle_4mom px " << particle_4mom.Px() << endl;
            cout << std::setprecision(10) << "particle_4mom py " << particle_4mom.Py() << endl;
            cout << std::setprecision(10) << "particle_4mom pz " << particle_4mom.Pz() << endl;
            cout << std::setprecision(10) << "particle_4mom E " << particle_4mom.E() << endl;
         }


         // if( pdg == 2212 || pdg == 2112 ){
         //    cout << "pdg " << pdg << endl;
         //    cout << "status" << status << endl;
         //    cout << "index " << index << endl;
         //    cout << "mass " << particle->GetM() << endl;
         //    cout << "particle_4mom px " << particle_4mom.Px() << endl;
         //    cout << "particle_4mom py " << particle_4mom.Py() << endl;
         //    cout << "particle_4mom pz " << particle_4mom.Pz() << endl;
         //    cout << "particle_4mom E " << particle_4mom.E() << endl;
         // }
         if( status == 1 ){
            
            cout << "pdg " << pdg << endl;
            cout << "status" << status << endl;
            cout << "index " << index << endl;
            cout << "mass " << particle->GetM() << endl;
            cout << "particle_4mom px " << particle_4mom.Px() << endl;
            cout << "particle_4mom py " << particle_4mom.Py() << endl;
            cout << "particle_4mom pz " << particle_4mom.Pz() << endl;
            cout << "particle_4mom E " << particle_4mom.E() << endl;

            total4Mom_outgoing += particle_4mom;   
         }
            
         ptHist.Fill(particle->GetPt());
         statusHist.Fill( status ); 

      } // for

      // cout << "outgoing pz " << total4Mom_outgoing.Pz() << endl;
      // cout << "outgoing E " << total4Mom_outgoing.E() << endl;


      double particle_pt = sqrt(total4Mom_outgoing.Px()*total4Mom_outgoing.Px() + total4Mom_outgoing.Py()*total4Mom_outgoing.Py());
      
      //1D histogram:
      Ntrk->Fill(nParticles);

      double energy_diff = total4Mom_incoming.E() - total4Mom_outgoing.E();
      energy_corr->Fill( energy_diff );
      px_corr->Fill( total4Mom_incoming.Px() - total4Mom_outgoing.Px() );
      py_corr->Fill( total4Mom_incoming.Py() - total4Mom_outgoing.Py() );
      pz_corr->Fill( total4Mom_incoming.Pz() - total4Mom_outgoing.Pz() );

      cout << "e diff " << energy_diff << endl;
      cout << "pz diff " << total4Mom_incoming.Pz() - total4Mom_outgoing.Pz() << endl;
      //2D histogram:
      energyVsQ2_2Dcorr->Fill(energy_diff, trueQ2);
      energyVsW2_2Dcorr->Fill(energy_diff, trueW2);
      energyVsX_2Dcorr->Fill(energy_diff, trueX);
      energyVsY_2Dcorr->Fill(energy_diff, trueY);
      energyVsNu_2Dcorr->Fill(energy_diff, trueNu);


   } // for

   TString outfilename;
   if( doBoost ) outfilename = "_DMPJet_eD_ionframe.root";
   else outfilename = "_DMPJet_eD.root";

   TFile output("../rootfiles/"+inputFilename+outfilename,"RECREATE");
   
   ptHist.Write();    
   statusHist.Write();
   pTvsThat->Write();
   Ntrk->Write();
   energy_corr->Write();
   px_corr->Write();
   py_corr->Write();
   pz_corr->Write();

   energyVsQ2_2Dcorr->Write();
   energyVsW2_2Dcorr->Write();
   energyVsX_2Dcorr->Write();
   energyVsY_2Dcorr->Write();
   energyVsNu_2Dcorr->Write();
   energyVsPf_2Dcorr->Write();
   energyVsPtf_2Dcorr->Write();
   energyVsProcess_2Dcorr->Write();

}

