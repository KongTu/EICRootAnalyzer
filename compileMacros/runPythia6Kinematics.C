#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include <eicsmear/erhic/EventBase.h>
#include <eicsmear/erhic/EventMC.h>
#include <eicsmear/erhic/EventPythia.h>
#include <eicsmear/erhic/Particle.h>
#include <eicsmear/erhic/ParticleMC.h>
#include <eicsmear/erhic/Pid.h>

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
#include "TLorentzRotation.h"
#include "TBranchElement.h"
#include "TRandom.h"

#define PI            3.1415926

#define MASS_MUON     0.1056
#define MASS_PION     0.13957
#define MASS_ELECTRON 0.000511
#define MASS_JPSI 	  3.09688
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_DEUTERON 1.8756129
#define MASS_TRITON   2.7937167208086358
#define MASS_HE3      2.7937167208086358
#define MASS_ALPHA    3.7249556277448477
#define MASS_LI6      5.5874334416172715
#define MASS_C12      11.174866883234543
#define MASS_CA40     37.249556277448477
#define MASS_XE131    121.99229680864376
#define MASS_AU197    183.45406466643374
#define MASS_PB208    193.69769264273208

using namespace std;
using namespace erhic;

double ME = 0.00051;

TLorentzRotation BoostToHCM(TLorentzVector const &eBeam_lab,
                            TLorentzVector const &pBeam_lab,
                            TLorentzVector const &eScat_lab) 
{
   TLorentzVector q_lab=eBeam_lab - eScat_lab;
   TLorentzVector p_plus_q=pBeam_lab + q_lab;
   // boost to HCM
   TLorentzRotation boost=TLorentzRotation(-1.0*p_plus_q.BoostVector());
   TLorentzVector pBoost=boost*pBeam_lab;
   TVector3 axis=pBoost.BoostVector();
   // rotate away y-coordinate
   boost.RotateZ(-axis.Phi());
   // rotate away x-coordinate
   boost.RotateY(M_PI-axis.Theta());

   TLorentzVector pBoost_escat=boost*eScat_lab;
   TVector3 axis_escat=pBoost_escat.BoostVector();

   // rotate away y-coordinate
   boost.RotateZ(-axis_escat.Phi());

   return boost;
}

TLorentzRotation BoostToHCM_es(TLorentzVector const &eBeam_lab,
                               TLorentzVector const &pBeam_lab,
                               TLorentzVector const &eScat_lab, 
                               double Q2_es, 
                               double y_es) {

   double escat_lab_es_E = (Q2_es)/(4.*eBeam_lab.E()) + eBeam_lab.E()*(1.-y_es);
   double b_par = 4.*eBeam_lab.E()*eBeam_lab.E()*(1.-y_es)/(Q2_es);
   double escat_lab_es_theta = TMath::ACos((1.-b_par)/(1.+b_par));
   
   double escat_lab_es_pz = sqrt(escat_lab_es_E*escat_lab_es_E - ME*ME)*TMath::Cos(escat_lab_es_theta);
   double escat_lab_es_pt = sqrt(escat_lab_es_E*escat_lab_es_E - ME*ME - escat_lab_es_pz*escat_lab_es_pz);
   double escat_lab_es_eta = -TMath::Log(TMath::Tan(escat_lab_es_theta/2.));
   double phi_elec = eScat_lab.Phi();

   TLorentzVector eScat_lab_ES;
   eScat_lab_ES.SetPtEtaPhiE(escat_lab_es_pt, escat_lab_es_eta, phi_elec, escat_lab_es_E);

   //same as before
   TLorentzVector q_lab=eBeam_lab - eScat_lab_ES;
   TLorentzVector p_plus_q=pBeam_lab + q_lab;
   // boost to HCM
   TLorentzRotation boost=TLorentzRotation(-1.0*p_plus_q.BoostVector());
   TLorentzVector pBoost=boost*pBeam_lab;
   TVector3 axis=pBoost.BoostVector();
   // rotate away y-coordinate
   boost.RotateZ(-axis.Phi());
   // rotate away x-coordinate
   boost.RotateY(M_PI-axis.Theta());

   TLorentzVector pBoost_escat=boost*eScat_lab;
   TVector3 axis_escat=pBoost_escat.BoostVector();

   // rotate away y-coordinate
   boost.RotateZ(-axis_escat.Phi());

   return boost;

}

//boost to HCM frame with kinematics from da/e method 
TLorentzRotation BoostToHCM_da(TLorentzVector const &eBeam_lab,
                               TLorentzVector const &pBeam_lab,
                               TLorentzVector const &eScat_lab, 
                               double Q2_da) {

   double E_da = Q2_da / (2 * eBeam_lab.E() * (1 + TMath::Cos( eScat_lab.Theta() )) );
   double escat_lab_da_pz = sqrt(E_da*E_da - ME*ME)*TMath::Cos(eScat_lab.Theta());
   double escat_lab_da_pt = sqrt(E_da*E_da - ME*ME - escat_lab_da_pz*escat_lab_da_pz);
   double escat_lab_da_eta = -TMath::Log(TMath::Tan(eScat_lab.Theta()/2.));
   TLorentzVector eScat_lab_DA;
   eScat_lab_DA.SetPtEtaPhiE(escat_lab_da_pt, escat_lab_da_eta, eScat_lab.Phi(), E_da);

   //same as before
   TLorentzVector q_lab=eBeam_lab - eScat_lab_DA;
   TLorentzVector p_plus_q=pBeam_lab + q_lab;

   // boost to HCM
   TLorentzRotation boost=TLorentzRotation(-1.0*p_plus_q.BoostVector());
   TLorentzVector pBoost=boost*pBeam_lab;
   TVector3 axis=pBoost.BoostVector();
   // rotate away y-coordinate
   boost.RotateZ(-axis.Phi());
   // rotate away x-coordinate
   boost.RotateY(M_PI-axis.Theta());

   TLorentzVector pBoost_escat=boost*eScat_lab;
   TVector3 axis_escat=pBoost_escat.BoostVector();

   // rotate away y-coordinate
   boost.RotateZ(-axis_escat.Phi());

   return boost;

}

TLorentzVector smearParticle( TLorentzVector part){

	bool isElectron=false;
	if( TMath::Abs(part.M()-ME) < 1e-5 ) isElectron=true;

	double resolution= 0.5/sqrt(part.E());
	resolution = sqrt(resolution*resolution + 0.02*0.02);//hadron
	if(isElectron){
		resolution = 0.12/sqrt(part.E());
		resolution = sqrt(resolution*resolution + 0.01*0.01);//electron
	}

	double phi = part.Phi();
	double dE = gRandom->Gaus(0,resolution);
	double E_new = part.E()*(1.+dE);
	double angular_res = 0.001;
	double dTheta = gRandom->Gaus(0,angular_res);
	double theta_new = part.Theta()*(1.+dTheta);

	part.SetE(E_new);
	double P_new = sqrt(part.E()*part.E() - part.M()*part.M());
	double Pz_new = P_new*TMath::Cos(theta_new);
	double Px_new = P_new*TMath::Sin(theta_new)*TMath::Cos(phi);
	double Py_new = P_new*TMath::Sin(theta_new)*TMath::Sin(phi);

	TLorentzVector part_new;
	part_new.SetPxPyPzE(Px_new,Py_new,Pz_new,E_new);

	return part_new;

}


void runPythia6Kinematics(const int nEvents = 1e5){

	//input from BeAGLE root files
	TFile * output = new TFile("../rootfiles/ep_kinematics.root","recreate");

	TChain *tree = new TChain("EICTree");
	tree->Add( "/gpfs02/eic/ztu/BeAGLE/BeAGLE_ep_2022-01-21/ep.root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	 // Histograms.
    TH1D* h_Epz = new TH1D("h_Epz","E-p_{z} (GeV)",100,0,70);

    TH1D* h_Q2_truth = new TH1D("h_Q2_truth",";Q^{2}_{truth}",200,5,1e4);
    TH1D* h_Q2_e = new TH1D("h_Q2_e",";Q^{2}_{electron}",200,5,1e4);
    TH1D* h_Q2_es = new TH1D("h_Q2_es",";Q^{2}_{e-#Sigma}",200,5,1e4);
    TH1D* h_Q2_da = new TH1D("h_Q2_da",";Q^{2}_{DA}",200,5,1e4);

    TH1D* h_y_truth = new TH1D("h_y_truth",";y_{truth}",200,0,1);
    TH1D* h_y_e = new TH1D("h_y_e",";y_{electron}",200,0,1);
    TH1D* h_y_es = new TH1D("h_y_es",";y_{e-#Sigma}",200,0,1);
    TH1D* h_y_da = new TH1D("h_y_da",";y_{DA}",200,0,1);

    TH1D* h_x_truth = new TH1D("h_x_truth",";x_{truth}",1000,0,1);
    TH1D* h_x_e = new TH1D("h_x_e",";x_{electron}",1000,0,1);
    TH1D* h_x_es = new TH1D("h_x_es",";x_{e-#Sigma}",1000,0,1);
    TH1D* h_x_da = new TH1D("h_x_da",";x_{DA}",1000,0,1);

    TH1D* h_phiStar_e = new TH1D("h_phiStar_e",";phiStar_{electron}",100,-3.14,3.14);
    TH1D* h_phiStar_es = new TH1D("h_phiStar_es",";phiStar_{e-#Sigma}",100,-3.14,3.14);
    TH1D* h_phiStar_da = new TH1D("h_phiStar_da",";phiStar_{DA}",100,-3.14,3.14);

    TH1D* h_Q2_es_res = new TH1D("h_Q2_es_res",";Q^{2}_{e} - Q^{2}_{es} / Q^{2}_{e}",100,-1,1);
    TH1D* h_Q2_da_res = new TH1D("h_Q2_da_res",";Q^{2}_{e} - Q^{2}_{da} / Q^{2}_{e}",100,-1,1);

    TH1D* h_y_es_res = new TH1D("h_y_es_res",";y^{2}_{e} - y^{2}_{es} / y^{2}_{e}",100,-1,1);
    TH1D* h_y_da_res = new TH1D("h_y_da_res",";y^{2}_{e} - y^{2}_{da} / y^{2}_{e}",100,-1,1);

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		if( (i%10000)==0 ) cout << "#Events = "<< i << endl;

		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		int struck_nucleon = event->nucleon;
		double MASS_NUCLEON = MASS_PROTON;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep+MASS_ELECTRON*MASS_ELECTRON));
		TLorentzVector p_beam(0.,0.,pztarg,sqrt(pztarg*pztarg+MASS_NUCLEON*MASS_NUCLEON));
		TLorentzVector scat_e(0.,0.,0.,0.);
		double s = (e_beam+p_beam).Mag2();

		//event information:
		double trueQ2 = event->GetTrueQ2();
		double trueX = event->GetTrueX();
		double trueY = event->GetTrueY();

		int nParticles = event->GetNTracks();

		//event cuts
		// if( event_process != 99 ) continue;
		if( trueQ2 < 180. || trueQ2 > 7000. ) continue;
		if( trueY > 0.8 || trueY < 0.2 ) continue;

		h_Q2_truth->Fill( trueQ2 );
		h_y_truth->Fill( trueY );
		h_x_truth->Fill( trueX );

		//particle loop

		/*
        Sigma e and Sigma h
        */
		TLorentzVector hfs(0,0,0,0);
		TLorentzVector part4v(0,0,0,0);
		for(int j(0); j < nParticles; ++j ) {
			const erhic::ParticleMC* particle = event->GetTrack(j);
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
			double pt = particle->GetPt();
			int status = particle->GetStatus();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			double theta = particle->GetTheta()*TMath::RadToDeg();
			int pdg = particle->GetPdgCode();
			int orig = particle->GetParentIndex();
			part4v = particle->Get4Vector();

			if( index == 3 ) {
				scat_e=particle->Get4Vector();
				scat_e = smearParticle(scat_e);
			}
			if( status!= 1 ) continue;
			if( (part4v-scat_e).P()<1e-4 ) continue;
			if( theta > 174 || theta < 4) continue;//LAr+SpaCal acceptance at H1.
			part4v = smearParticle( part4v );
			hfs += part4v;
		} // end of particle loop
		TLorentzVector q_beam = e_beam-scat_e;

		double sigma_e = scat_e.E() - scat_e.Pz();
        double sigma_h = hfs.E() - hfs.Pz();
        double EPz = sigma_e + sigma_h;
        h_Epz->Fill( EPz );

		/*
        electron method
        */

        double Q2_e    = -q_beam.Mag2();
        double x_e     = Q2_e / (2. * p_beam * q_beam);
        double y_e     = (p_beam * q_beam) / (p_beam * e_beam);
      
        // Fill kinematics histograms.
        h_Q2_e->Fill( Q2_e );
        h_y_e->Fill( y_e );
        h_x_e->Fill( x_e );

        // if(EPz<35.||EPz>70.) continue;
        /*
        e-sigma method
        */
        double Q2_es = Q2_e;
        double y_sigma = sigma_h / EPz;
        double Q2_sigma =  TMath::Power(scat_e.Pt(),2) / (1. - y_sigma);
        double x_sigma = Q2_sigma / (s*y_sigma);
        double x_es = x_sigma;
        double y_es = Q2_es / (s*x_es);

        h_Q2_es->Fill( Q2_es );
        h_y_es->Fill( y_es );
        h_x_es->Fill( x_es );

        /*
        double-angle method
        */
        double theta_h = TMath::ATan( sigma_h / hfs.Pt() )*2.;
        double y_da = TMath::Tan( theta_h / 2.) / (TMath::Tan(scat_e.Theta()/2.) + TMath::Tan( theta_h / 2.));
        double Q2_da = 4.*TMath::Power(e_beam.E(),2)* (1./TMath::Tan(scat_e.Theta()/2.) ) / (TMath::Tan(scat_e.Theta()/2.) + TMath::Tan( theta_h / 2.));
        double x_da = Q2_da / (s*y_da);

        h_Q2_da->Fill( Q2_da );
        h_y_da->Fill( y_da );
        h_x_da->Fill( x_da );

        h_Q2_es_res->Fill( (Q2_e - Q2_es)/Q2_e );
        h_Q2_da_res->Fill( (Q2_e - Q2_da)/Q2_e );

        h_y_es_res->Fill( (y_e - y_es)/y_e );
        h_y_da_res->Fill( (y_e - y_da)/y_e );

        /*
        Boost to HCM frame for phiStar
        */
        // Loop over final particles in the event.
        TLorentzRotation use_e_to_rotate;
        TLorentzRotation use_es_to_rotate;
        TLorentzRotation use_da_to_rotate;

        TLorentzVector part_e;
        TLorentzVector part_es;
        TLorentzVector part_da;

        for(int i(0); i < nParticles; ++i ) {
			const erhic::ParticleMC* particle = event->GetTrack(i);
			double pt = particle->GetPt();
			int status = particle->GetStatus();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			double mass = particle->GetM();
			part4v = particle->Get4Vector();
			if(status!=1) continue;
			part4v = smearParticle( part4v );
        	if( (part4v-scat_e).P() < 1e-4 ) continue; //skip e'
            if(part4v.Pt() < 0.15 || TMath::Abs(part4v.Eta())>1.75) continue;
            double zhad = p_beam.Dot(part4v) / p_beam.Dot(q_beam);
            if(zhad < 0.2) continue;

            use_e_to_rotate = BoostToHCM(e_beam,p_beam,scat_e);
            use_es_to_rotate = BoostToHCM_es(e_beam,p_beam,scat_e,Q2_es,y_es);
            use_da_to_rotate = BoostToHCM_da(e_beam,p_beam,scat_e,Q2_da);

            part_e = use_e_to_rotate*part4v;
            part_es = use_es_to_rotate*part4v;
            part_da = use_da_to_rotate*part4v;

            h_phiStar_e->Fill( part_e.Phi() );
            h_phiStar_es->Fill( part_es.Phi() );
            h_phiStar_da->Fill( part_da.Phi() );
            
        }

	}

	output->Write();
	output->Close();


}
