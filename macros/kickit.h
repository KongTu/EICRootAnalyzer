vector<TLorentzVector> kickit(TLorentzVector particle_4mom_neutron_bKick, TLorentzVector particle_4mom_proton_bKick, TLorentzVector particle_4mom_jpsi_bKick){

   TLorentzVector t,k;
   TLorentzVector p3,p4,p5;
   TLorentzVector particle_4mom_proton,particle_4mom_neutron,particle_4mom_jpsi;

   t = particle_4mom_neutron_bKick + particle_4mom_proton_bKick + particle_4mom_jpsi_bKick;

   TF1 *fa = new TF1("fa","[0]*TMath::Abs(TMath::Exp([1]*x))",0,10);
   fa->SetParameter(0,1);
   fa->SetParameter(1,-3);

   double kick_px = 0.;
   double kick_py = fa->GetRandom();

   TF1 *phiran = new TF1("phiran","[0]*1",-PI,PI);
   phiran->SetParameter(0,1);
   double phi_kick = phiran->GetRandom();

   if( phi_kick > 0 ){
      kick_px = kick_py/TMath::Tan(phi_kick);
   }
   else if( phi_kick < 0 ) {
      kick_py = -kick_py;
      kick_px = kick_py/TMath::Tan(phi_kick);
   }

   px_dist->Fill( kick_px );
   py_dist->Fill( kick_py );

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

   double E_min = 1.0;
   double comp_min = 0.;
   double delta_min = 0.;
   double kappa_min = 0.;

   int i_min = 0;
   int j_min = 0;
   int k_min = 0;

   double comp_init = -50;
   double delta_init = -5;
   double kappa_init = -5;

   const int iteration_1 = 100;
   const int iteration_2 = 10;

   double comp[iteration_1];
   double delta[iteration_2];
   double kappa[iteration_2];

   for(int jter = 0; jter < iteration_1; jter++){

      double temp = comp_init+1.*jter;
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

         double p_py_prime = p_py + kick_py;
         double n_py_prime = n_py - kick_py + delta[iter];
         double j_py_prime = j_py - delta[iter];

         double p_px_prime = p_px + kick_px; 
         double n_px_prime = n_px - kick_px + kappa[kter];
         double j_px_prime = j_px - kappa[kter];

         double p_pz_prime = p_pz + comp[jter];
         double n_pz_prime = n_pz - comp[jter];
         double j_pz_prime = j_pz;

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

         }
      }
   }

   vector<TLorentzVector> ParticleCollection;
   ParticleCollection.push_back( particle_4mom_proton );
   ParticleCollection.push_back( particle_4mom_neutron );
   ParticleCollection.push_back( particle_4mom_jpsi );
   
   if( i_min == 0 || j_min == 0 || k_min == 0 || i_min == 9 || j_min == 99 || k_min == 9 ) return;//hit the boundary continue;

   return ParticleCollection;
   // cout << "iter: " << i_min << " jter: " << j_min << " kter: " << k_min << endl;
   // cout << "E diff: " << E_min <<  " comp: " << comp_min << " delta: " << delta_min << " kappa: " << kappa_min << endl;
}