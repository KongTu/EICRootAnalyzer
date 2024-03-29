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
#include "TBranchElement.h"

#define PI            3.1415926

#define MASS_MUON     0.1056
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


void runUPC_BeAGLE(const TString filename="eA_TEST", const int nEvents = 40000, bool reweight_=true, const TString collider="RHIC", TString system="eA"){

	TChain *tree = new TChain("EICTree");
	tree->Add( filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	TFile* input = new TFile("eicToUPC_"+collider+".root","READ");
	TH1D* eicToUPC=(TH1D*) input->Get("h_photonWeight");

	TFile* output=new TFile("../rootfiles/UPC_BeAGLE_"+collider+"_"+system+".root","RECREATE");
	
	TH1D* h_trueQ2[2],*h_trueNu[2],*h_trueW[2],*h_trueX[2], *h_Ntrk[2],
	*h_Nevap[2],*h_Tb[2],*h_b[2],*h_d[2],
	*h_charged_eta[2],*h_charged_pt[2],
	*h_all_eta[2],*h_all_pt[2];
	TH1D* h_measQ2[2],*h_measW[2], *h_measX[2];
	TH2D* h_photonFluxVsNu[2],*h_trueWvsNevap[2], *h_Wsmear[2];

	for(int j=0;j<2;j++){
		h_trueQ2[j] = new TH1D(Form("h_trueQ2_%d",j),";Q^{2} (GeV^{2})",100,1e-4,1);
		h_trueNu[j] = new TH1D(Form("h_trueNu_%d",j),";#nu (GeV)",180,1e-3,18);
		h_trueW[j] = new TH1D(Form("h_trueW_%d",j),";W (GeV)",500,0,500);
		h_trueX[j] = new TH1D(Form("h_trueX_%d",j),";x",500,1e-5,1e-1);
		h_Ntrk[j] = new TH1D(Form("h_Ntrk_%d",j),";N",100,-0.5,99.5);

		h_Nevap[j] = new TH1D(Form("h_Nevap_%d",j),";N_{neutron}",60,-0.5,59.5);
		h_Tb[j] = new TH1D(Form("h_Tb_%d",j),";T_{b}",60,0,16);
		h_b[j] = new TH1D(Form("h_b_%d",j),";b",60,0,10);
		h_d[j] = new TH1D(Form("h_d_%d",j),";d",60,0,16);

		h_charged_eta[j] = new TH1D(Form("h_charged_eta_%d",j),";#eta",100,-3,10);
		h_charged_pt[j] = new TH1D(Form("h_charged_pt_%d",j),";p_{T}",100,0,10);

		h_all_eta[j] = new TH1D(Form("h_all_eta_%d",j),";#eta",100,-3,10);
		h_all_pt[j] = new TH1D(Form("h_all_pt_%d",j),";p_{T}",100,0,10);
		
		h_measQ2[j] = new TH1D(Form("h_measQ2_%d",j),";Q^{2} (GeV^{2})",100,1e-4,1);
		h_measW[j] = new TH1D(Form("h_measW_%d",j),";W (GeV)",500,0,500);
		h_measX[j] = new TH1D(Form("h_measX_%d",j),";x",500,1e-5,1e-1);

		h_photonFluxVsNu[j] = new TH2D(Form("h_photonFluxVsNu_%d",j),";#nu (GeV);flux",180,1e-3,18,100,0,0.3);
		h_trueWvsNevap[j] = new TH2D(Form("h_trueWvsNevap_%d",j),";trueW;Nevap",500,0,500,60,-0.5,59.5);
		h_Wsmear[j] = new TH2D(Form("h_Wsmear_%d",j),";trueW;measW",500,0,500,500,0,500);
	}
	

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		if( i%10000 == 0 ) cout << "Number of events = " << i << endl;
		
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		double Atarg = event->Atarg;
		double pztarg_total = pztarg*Atarg;
		int struck_nucleon = event->nucleon;
		double MASS_NUCLEON = MASS_PROTON;
		if( struck_nucleon==2112 ) MASS_NUCLEON = MASS_NEUTRON;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep+MASS_ELECTRON*MASS_ELECTRON));
		TLorentzVector p_beam(0.,0.,pztarg,sqrt(pztarg*pztarg+MASS_NUCLEON*MASS_NUCLEON));
		TLorentzVector e_scattered(0.,0.,0.,0.);
		TLorentzVector Au_beam(0.,0.,-pztarg_total,sqrt(pztarg_total*pztarg_total+MASS_AU197*MASS_AU197));
		TLorentzVector photon_emitter=Au_beam;
		if(system=="ep") photon_emitter.SetPxPyPzE(0.,0.,-pztarg,sqrt(pztarg*pztarg+MASS_NUCLEON*MASS_NUCLEON));


		//event information:
		double trueQ2 = event->GetTrueQ2();
		double trueW2 = event->GetTrueW2();
		double trueX = event->GetTrueX();
		double trueY = event->GetTrueY();
		double trueNu = event->nu;
		double s_hat = event->GetHardS();
		double t_hat = event->t_hat;
		double u_hat = event->GetHardU();
		double photon_flux = event->GetPhotonFlux();
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();
		
		double impact_parameter = event->b;
		double Tb = event->Thickness;
		double distance = event->d1st;
		int N_nevap = event->Nnevap;
		int N_pevap = event->Npevap;

		const erhic::ParticleMC* particle_escat = event->GetTrack(2);
		e_scattered=particle_escat->Get4Vector();
		TLorentzVector q=(e_beam-e_scattered);
		double phot_energy=q.E();
		double weight=eicToUPC->GetBinContent(eicToUPC->FindBin(phot_energy));
		if(!reweight_) weight=1.0;
		if(phot_energy>11.&&collider=="RHIC") weight=0.;
		//Event histograms
		h_trueQ2[0]->Fill( trueQ2, 1);
		h_trueW[0]->Fill(sqrt(trueW2), 1);
		h_trueNu[0]->Fill( q.E(), 1); 
		h_trueX[0]->Fill( trueX, 1); 
		h_Nevap[0]->Fill(N_nevap, 1);
		h_Tb[0]->Fill(Tb, 1);
		h_b[0]->Fill(impact_parameter, 1);
		h_d[0]->Fill(distance, 1);
		h_photonFluxVsNu[0]->Fill(q.E(), photon_flux,1);
		h_trueWvsNevap[0]->Fill(sqrt(trueW2),N_nevap,1);

		h_trueQ2[1]->Fill( trueQ2, weight);
		h_trueW[1]->Fill(sqrt(trueW2), weight);
		h_trueNu[1]->Fill( q.E(), weight);
		h_trueX[1]->Fill( trueX, weight); 
		h_Nevap[1]->Fill(N_nevap, weight);
		h_Tb[1]->Fill(Tb, weight);
		h_b[1]->Fill(impact_parameter, weight);
		h_d[1]->Fill(distance, weight);
		h_trueWvsNevap[1]->Fill(sqrt(trueW2),N_nevap,weight);

		//cut on energy
		if(trueNu<2.5) continue;

		TLorentzVector hfs(0,0,0,0);
		//particle loop
		int Ntrk=0;
		for(int j(0); j < nParticles; ++j ) {

			const erhic::ParticleMC* particle = event->GetTrack(j);

			int pdg = particle->GetPdgCode();
			int status = particle->GetStatus();
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
			double pt = particle->GetPt();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			double rap = particle->GetRapidity();
			double mass = particle->GetM();
			double theta = particle->GetTheta(); 
			theta = theta*1000.0; //change to mrad;
			double mom = particle->GetP();
			int charge= particle->eA->charge;
			if( status!= 1) continue;
			if(TMath::Abs(particle->Get4Vector().E()-e_scattered.E())<1e-1) continue;//no scat e
			if((eta>-1.5&&eta<1.5&&pt>0.2)){//||(eta>2.5&&eta<4.0)){//STAR forward upgrade acceptance
				hfs+=particle->Get4Vector();
				if(charge!=0)Ntrk++;
			}
			if(pt<0.2) continue;
			h_all_pt[0]->Fill(pt, 1.);
			h_all_eta[0]->Fill(eta, 1.);
			h_all_pt[1]->Fill(pt, weight);
			h_all_eta[1]->Fill(eta, weight);
		
			if( charge==0 ) continue;
			//charged particles
			h_charged_pt[0]->Fill(pt, 1.);
			h_charged_eta[0]->Fill(eta, 1.);
			h_charged_pt[1]->Fill(pt, weight);
			h_charged_eta[1]->Fill(eta, weight);

		} // end of particle loop
		
		h_Ntrk[0]->Fill(Ntrk, 1);
		h_Ntrk[1]->Fill(Ntrk, weight);

		//Hadron only method.
		double sigma_had=hfs.E()-hfs.Pz();
		double measY=sigma_had/(2*photon_emitter.E());
		double measQ2=hfs.Pt()*hfs.Pt()/(1-measY);
		double s=(photon_emitter+p_beam).Mag2();
		double measW=sqrt(s*measY - measQ2 - MASS_NUCLEON*MASS_NUCLEON);
		double measX=measQ2/(s*measY);
	
		h_measQ2[0]->Fill(measQ2, 1);
		h_measW[0]->Fill(measW, 1);
		h_measX[0]->Fill(measX, 1);

		h_measQ2[1]->Fill(measQ2, weight);
		h_measW[1]->Fill(measW, weight);
		h_measX[1]->Fill(measX, weight);

		h_Wsmear[0]->Fill(sqrt(trueW2),measW,1);
		h_Wsmear[1]->Fill(sqrt(trueW2),measW,weight);

	}

	TH1D* h_photonFlux_nu = new TH1D("h_photonFlux_nu",";#nu (GeV)",180,1e-3,18);
	for(int ibin=0;ibin<h_photonFlux_nu->GetNbinsX();ibin++){
		TH1D* TEMP=new TH1D(Form("temp_%d",ibin),";",100,0,0.3);
		double flux=h_photonFluxVsNu[0]->ProjectionY(Form("temp_%d",ibin),ibin+1,ibin+1)->GetMean();
		h_photonFlux_nu->SetBinContent(ibin+1, flux);
		TEMP->Delete();
	}

	output->Write();
	output->Close();

}
