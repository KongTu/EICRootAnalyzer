#include "hist.h"//define all the histograms
#include "PRINT4VECTOR.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
using namespace std;
using namespace erhic;

#define MASS_MUON  0.1056

TH1D* that = new TH1D("that","that",200,0,10);
TH1D* tjpsi = new TH1D("tjpsi","tjpsi",200,0,10);
TH1D* sPN = new TH1D("sPN","sPN",200,0,14);
TH1D* sPN_4pt2 = new TH1D("sPN_4pt2","sPN_4pt2",200,0,14);
TH1D* sPN_Jpsi = new TH1D("sPN_Jpsi","sPN_Jpsi",200,0,14);
TH1D* sPN_Jpsi_fix = new TH1D("sPN_Jpsi_fix","sPN_Jpsi_fix",200,0,14);
TH1D* sPN_Jpsi_fix_oneTagged = new TH1D("sPN_Jpsi_fix_oneTagged","sPN_Jpsi_fix_oneTagged",200,0,14);
TH1D* h_trk = new TH1D("h_trk","h_trk",50,0,50);

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

void eD_SRC_main(const int nEvents = 40000, TString filename=""){

	TChain *tree = new TChain("EICTree");
	tree->Add("/eicdata/eic0003/ztu/BeAGLE_devK/"+filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		double Atarg = event->Atarg;
		double pztarg_total = pztarg*Atarg;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep));//neglecting e mass
		TLorentzVector d_beam(0.,0.,pztarg_total,sqrt(pztarg_total*pztarg_total+MASS_DEUTERON*MASS_DEUTERON));
		TLorentzVector e_scattered(0.,0.,0.,0.);

		double boostv =  d_beam.Pz()/d_beam.E();
		TVector3 b;

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
		int struck_nucleon = event->nucleon;
		
		if( event_process != 91 ) continue;
		if( trueQ2 < 1. ) continue;
		if( trueY > 0.85 || trueY < 0.05 ) continue;
		bool struckproton = false;
		if( struck_nucleon == 2212 ) struckproton = true;

		that->Fill( fabs(t_hat) );

		int nParticles_process = 0;
		TLorentzVector p_4vect, n_4vect,j_4vect,q;
		TLorentzVector p_4vect_irf, n_4vect_irf,j_4vect_irf,q_irf;
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
			//photon 4vector
			q = e_beam-e_scattered;
			q_irf = q;
			
			TLorentzVector ppart = particle->Get4Vector();
			if(pdg == 443 ) j_4vect = ppart;
			if(pdg == 2212) p_4vect = ppart;
			if(pdg == 2112) n_4vect = ppart;

			ppart.Boost(0,0,-boostv);
			ppart.Boost(b);

			if(pdg == 443 ) j_4vect_irf = ppart;
			if(pdg == 2212) p_4vect_irf = ppart;
			if(pdg == 2112) n_4vect_irf = ppart; 

			q_irf.Boost(0,0,-boostv);
			q_irf.Boost(b);

			nParticles_process++;

		} // end of particle loop

		double pt2 = j_4vect.Pt()*j_4vect.Pt();
		tjpsi->Fill( pt2-trueQ2 );
		h_trk->Fill( nParticles_process );
		
		if( struckproton ){
			TLorentzVector n_partner_4vect_irf;
			n_partner_4vect_irf.SetPxPyPzE(-n_4vect_irf.Px(), -n_4vect_irf.Py(), -n_4vect_irf.Pz(), sqrt(n_4vect_irf.P()*n_4vect_irf.P()+MASS_PROTON*MASS_PROTON) );
			//intrinsic Beagle n(k) distribution -> SNN distribution
			sPN->Fill( (n_partner_4vect_irf+n_4vect_irf).Mag2() );
			//approximation by spectator nucleon pt in the lab frame
			sPN_4pt2->Fill( 4*n_4vect.Pt()*n_4vect.Pt() );
		} 
		else{
			TLorentzVector p_partner_4vect_irf;
			p_partner_4vect_irf.SetPxPyPzE(-p_4vect_irf.Px(), -p_4vect_irf.Py(), -p_4vect_irf.Pz(), sqrt(p_4vect_irf.P()*p_4vect_irf.P()+MASS_NEUTRON*MASS_NEUTRON) );
			sPN->Fill( (p_partner_4vect_irf+p_4vect_irf).Mag2() );
			sPN_4pt2->Fill( 4*p_4vect.Pt()*p_4vect.Pt() );
		}

		//inclusive J/psi measurement, convolution of exp and intrinsic n(k)
		sPN_Jpsi->Fill( (p_4vect_irf+n_4vect_irf).Mag2() );
		//remove dependence of Q2 and Jpsi production
		sPN_Jpsi_fix->Fill( (p_4vect_irf+n_4vect_irf+j_4vect_irf-q_irf).Mag2() );

		TLorentzVector n_4vect_cal; 
		n_4vect_cal = q_irf - j_4vect_irf - p_4vect_irf;
		//this one doesn't work, because we loose the information of how relative pn moves.
		//therefore by definition, this is going to be zero. 
		sPN_Jpsi_fix_oneTagged->Fill( (p_4vect_irf+n_4vect_cal+j_4vect_irf-q_irf).Mag2() );

		cout << "test ~ " << (p_4vect_irf+n_4vect_irf+j_4vect_irf-q_irf).Pz() << endl;
	}

	TFile output("../rootfiles/eD_SRC_main_Beagle.root","RECREATE");
	that->Write();
	tjpsi->Write();
	sPN->Write();
	sPN_4pt2->Write();
	sPN_Jpsi->Write();
	sPN_Jpsi_fix->Write();
	sPN_Jpsi_fix_oneTagged->Write();
	h_trk->Write();



}