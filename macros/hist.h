#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"

#define PI            3.1415926

#define MASS_JPSI 	  3.09688
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_DEUTERON 1.8624778138724238
#define MASS_TRITON   2.7937167208086358
#define MASS_HE3      2.7937167208086358
#define MASS_ALPHA    3.7249556277448477
#define MASS_LI6      5.5874334416172715
#define MASS_C12      11.174866883234543
#define MASS_CA40     37.249556277448477
#define MASS_XE131    121.99229680864376
#define MASS_AU197    183.45406466643374
#define MASS_PB208    193.69769264273208
// Histograms for our analysis.

//TH2D DIS kinematics:
TH2D* Q2VsX = new TH2D("Q2VsX",";x;Q^{2}",10000,0.00001,1,2000,0,200);
TH2D* W2VsFlux = new TH2D("W2VsFlux",";#Phi;W^{2}",5000,0,0.1,1000,0,10000);

TH2D* T_hatVsPt2 = new TH2D("T_hatVsPt2",";p^{2}_{T}-Q^{2} (GeV);T", 500,0,50,200,-5,0);
TH2D* TvsPt = new TH2D("TvsPt",";p_{T} (GeV);T", 100,0,10,200,-5,0);
TH2D* tProtonVsPt = new TH2D("tProtonVsPt",";p_{T} (GeV);T", 100,0,10,200,-5,0);

TH2D* ThatVssNN = new TH2D("ThatVssNN",";s_{_{NN}} (GeV^{2});T" ,200,0,20,200,-5,0);
TH2D* tdisVssNN = new TH2D("tdisVssNN",";s_{_{NN}} (GeV^{2});t" ,200,0,20,400,-5,5);
TH2D* tVsT = new TH2D("tVsT",";T;t" ,200,-5,0,400,-5,5);

TH2D* sNNvsPt = new TH2D("sNNvsPt",";p_{T} (GeV);s_{_{NN}} (GeV^{2})", 100,0,10,200,0,20);
TH2D* Q2VsJpsi = new TH2D("Q2VsJpsi",";p_{T};Q^{2}",100,0,10,2000,0,200);
TH2D* W2VsJpsi = new TH2D("W2VsJpsi",";p_{T};W^{2}",100,0,10,1000,0,10000);

TH2D* PtVsEta_proton = new TH2D("PtVsEta_proton",";#eta;p_{T} (GeV)", 2000, -20,20, 200, 0,10);
TH2D* PtVsEta_neutron = new TH2D("PtVsEta_neutron",";#eta;p_{T} (GeV)", 2000, -20,20, 200, 0,10);

TH2D* PtVsPt_protonVsJpsi = new TH2D("PtVsPt_protonVsJpsi",";p_{T} (GeV);p_{T} (GeV)", 200, 0,10, 200, 0,10);

TH2D* AngleVsMom_proton = new TH2D("AngleVsMom_proton",";p (GeV);#theta (mrad)",2500,0,250,300,0,30);
TH2D* AngleVsMom_neutron = new TH2D("AngleVsMom_neutron",";p (GeV);#theta (mrad)",2500,0,250,300,0,30);

TH2D* AngleVssNN_proton = new TH2D("AngleVssNN_proton",";s_{_{NN}} (GeV^{2});#theta (mrad)",2500,0,10,300,0,30);
TH2D* AngleVssNN_neutron = new TH2D("AngleVssNN_neutron",";s_{_{NN}} (GeV^{2});#theta (mrad)",2500,0,10,300,0,30);

//TH1D event variables:
TH1D* Ntrk_process_all = new TH1D("Ntrk_process_all",";Ntrk_process_all", 100, 0, 100);
TH1D* Ntrk_process = new TH1D("Ntrk_process",";Ntrk_process", 100, 0, 100);
TH1D statusHist("statusHist", "status distribution  ", 50, 0, 50 );

TH1D* E_CM = new TH1D("E_CM",";#sqrt{s_{_{NN}}} (GeV)", 200,0,10);
TH1D* W2 = new TH1D("W2",";W^{2} (GeV^{2})", 1000,0,10000);
TH1D* photonFlux = new TH1D("photonFlux",";#Phi", 5000,0,0.1);
TH1D* T_dist = new TH1D("T_dist",";T", 200,-5,0);
TH1D* T_Jpsi_dist = new TH1D("T_Jpsi_dist",";T", 200,-5,0);
TH1D* t_dist = new TH1D("t_dist",";t", 200,-5,5);

//TH1D:
TH1D* PtDist_Jpsi = new TH1D("PtDist_Jpsi",";PtDist_Jpsi", 200, 0,10);
TH1D* PhiDist_Jpsi = new TH1D("PhiDist_Jpsi",";PhiDist_Jpsi", 200, -10,10);
TH1D* EtaDist_Jpsi = new TH1D("EtaDist_Jpsi",";EtaDist_Jpsi", 2000, -20,20);

TH1D* PtDist_proton = new TH1D("PtDist_proton",";PtDist_proton", 200, 0,10);
TH1D* PhiDist_proton = new TH1D("PhiDist_proton",";PhiDist_proton", 200, -10,10);
TH1D* EtaDist_proton = new TH1D("EtaDist_proton",";EtaDist_proton", 2000, -20,20);

TH1D* PtDist_neutron = new TH1D("PtDist_neutron",";PtDist_neutron", 200, 0,10);
TH1D* PhiDist_neutron = new TH1D("PhiDist_neutron",";PhiDist_neutron", 200, -10,10);
TH1D* EtaDist_neutron = new TH1D("EtaDist_neutron",";EtaDist_neutron", 2000, -20,20);

TH1D* px_dist = new TH1D("px_dist",";px",1000,-10,10);
TH1D* py_dist = new TH1D("py_dist",";py",1000,-10,10);

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
   
   //if( i_min == 0 || j_min == 0 || k_min == 0 || i_min == 9 || j_min == 99 || k_min == 9 ) return ;//hit the boundary continue;

   return ParticleCollection;
   // cout << "iter: " << i_min << " jter: " << j_min << " kter: " << k_min << endl;
   // cout << "E diff: " << E_min <<  " comp: " << comp_min << " delta: " << delta_min << " kappa: " << kappa_min << endl;
}