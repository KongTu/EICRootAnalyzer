#include "hist.h"//define all the histograms
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
using namespace std;
using namespace erhic;

#define MASS_MUON  0.1056
#define MASS_ELECTRON  0.00051
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957

TH1D* dNdetaStar = new TH1D("dNdetaStar","dNdetaStar",100,-10,10);
TH1D* dNdetaStar_p = new TH1D("dNdetaStar_p","dNdetaStar_p",100,-10,10);
TH1D* dNdetaStar_m = new TH1D("dNdetaStar_m","dNdetaStar_m",100,-10,10);
TH1D* dNdeta = new TH1D("dNdeta","dNdeta",100,-10,10);
TH1D* h_trk = new TH1D("h_trk","h_trk",4000,0,4000);

TH1D* h_neutE = new TH1D("h_neutE","E (GeV)", 100,0,10000);
TH1D* h_nNeutrons = new TH1D("h_nNeutrons","h_nNeutrons", 100,0,100);

TH2D* h_NpevapVsNnevap = new TH2D("h_NpevapVsNnevap",";Nnevap;Npevap",40,0,40,10,0,10);
TH2D* h_NnevapVsb = new TH2D("h_NnevapVsb",";Nnevap;b",40,0,40,100,0,15);
TH2D* h_NnevapVsTb = new TH2D("h_NnevapVsTb",";Nnevap;Tb",40,0,40,100,0,15);
TH2D* h_bVsTb = new TH2D("h_bVsTb",";b (fm);Tb (fm)", 100,0,10,100,0,15);

TH2D* h_neutEVsb = new TH2D("h_neutEVsb",";neutE;b",100,0,10000,100,0,15);
TH2D* h_neutEVsTb = new TH2D("h_neutEVsTb",";neutE;Tb",100,0,10000,100,0,15);

TH1D* h_particleE = new TH1D("h_particleE","E (GeV)", 100,0,100);
TH2D* h_particleEVsb = new TH2D("h_particleEVsb",";particleE;b",100,0,100,100,0,15);
TH2D* h_particleEVsTb = new TH2D("h_particleEVsTb",";particleE;Tb",100,0,100,100,0,15);

TLorentzRotation BoostToHCM(TLorentzVector const &eBeam_lab,
                            TLorentzVector const &pBeam_lab,
                            TLorentzVector const &eScat_lab) {
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
   return boost;
}

void makeBeAGLE_centrality(const int nEvents = 40000){

	TChain *tree = new TChain("EICTree");
	tree->Add("/eicdata/eic0003/ztu/BeAGLE_Centrality/ePb.root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		int struck_nucleon = event->nucleon;
		double MASS_NUCLEON = MASS_PROTON;
		if( struck_nucleon==2112 ) MASS_NUCLEON = MASS_NEUTRON;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep+MASS_ELECTRON*MASS_ELECTRON));
		TLorentzVector p_beam(0.,0.,pztarg,sqrt(pztarg*pztarg+MASS_NUCLEON*MASS_NUCLEON));
		TLorentzVector e_scattered(0.,0.,0.,0.);

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
		int N_nevap = event->Nnevap;
		int N_pevap = event->Npevap;

		if( event_process != 99 ) continue;
		if( trueQ2 < 1. || trueQ2 > 20. ) continue;
		if( trueY > 0.95 || trueY < 0.1 ) continue;

		h_NpevapVsNnevap->Fill(N_nevap,N_pevap);
		h_NnevapVsb->Fill(N_nevap,impact_parameter);
		h_NnevapVsTb->Fill(N_nevap,Tb);
		h_bVsTb->Fill(Tb,impact_parameter);

		int nParticles_process = 0;
		int nNeutrons = 0;
		double sumNeutronEnergy = 0.;
		double sumChargeParticleEnergy = 0.;
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

			if( index == 3 ) {
				e_scattered.SetPtEtaPhiM(pt,eta,phi,mass);
			}
			if( status != 1 ) continue;
			
			//ZDC neutron selection
			if( pdg == 2112 ){
				if( TMath::Abs(theta) < 4.0 ){
					nNeutrons++;
					sumNeutronEnergy += sqrt(mom*mom + mass*mass);
				}
			}
			//central detector total energy
			if( fabs(pdg) == 211 || fabs(pdg) == 321 || fabs(pdg) == 2212 ){
				if( pt > 0.2 && fabs(eta) < 4.0 ){
					sumChargeParticleEnergy += sqrt(mom*mom + mass*mass);
				}
			}

			if( mom < 0.2 || mom > 10. ) continue;
			if( fabs(pdg) != 211 && fabs(pdg) != 321 && fabs(pdg) != 2212 ) continue;
			
			TLorentzRotation boost_MC_HCM = BoostToHCM(e_beam,p_beam,e_scattered);	
			TLorentzVector h = particle->Get4Vector();
			TLorentzVector hStar = boost_MC_HCM*h;

			dNdetaStar->Fill( hStar.Rapidity() );
			dNdeta->Fill( eta );
			nParticles_process++;

			if(	pdg > 0 ) dNdetaStar_p->Fill( hStar.Rapidity() );
			if(	pdg < 0 ) dNdetaStar_m->Fill( hStar.Rapidity() );

		} // end of particle loop

		h_neutEVsb->Fill( sumNeutronEnergy, impact_parameter );
		h_neutEVsTb->Fill( sumNeutronEnergy, Tb );
		h_neutE->Fill( sumNeutronEnergy );
		h_nNeutrons->Fill( nNeutrons );
		
		h_particleE->Fill( sumChargeParticleEnergy );
		h_particleEVsb->Fill( sumChargeParticleEnergy, impact_parameter );
		h_particleEVsTb->Fill( sumChargeParticleEnergy, Tb );


		h_trk->Fill( nParticles_process );

	}

	TFile output("../rootfiles/ePb_18x135_centrality.root","RECREATE");
	dNdetaStar->Write();
	dNdetaStar_p->Write();
	dNdetaStar_m->Write();
	dNdeta->Write();
	h_trk->Write();
	h_neutE->Write();
	h_particleE->Write();
	h_nNeutrons->Write();
	h_NpevapVsNnevap->Write();
	h_NnevapVsb->Write();
	h_NnevapVsTb->Write();
	h_bVsTb->Write();
	h_neutEVsb->Write();
	h_neutEVsTb->Write();
	h_particleEVsb->Write();
	h_particleEVsTb->Write();


}