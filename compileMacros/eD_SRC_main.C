#include "hist.h"//define all the histograms
#include "PRINT4VECTOR.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
using namespace std;
using namespace erhic;

#define MASS_MUON  0.1056

double sPN_bins[]={0.,1.0,2.0,3.0,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.7,5.0,5.5,6.0,7.0,8.0,9.0,10.0,12.0,15.0};
int sPN_nBins = sizeof(sPN_bins)/sizeof(sPN_bins[0]) -1;

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

void eD_SRC_main(const int nEvents = 40000, TString filename="", const bool doSmear_ = false){

	TFile * output = new TFile("../rootfiles/eD_SRC_main_Beagle.root","recreate");
	TH1D* h_trk = new TH1D("h_trk","h_trk",50,0,50);
	TH1D* that = new TH1D("that","that",200,0,10);
	TH1D* tjpsi = new TH1D("tjpsi","tjpsi",200,0,10);
	TH2D* nRes = new TH2D("nRes","",60,-30,30,20,-0.1,0.1);
	TH1D* nucleon_t = new TH1D("nucleon_t","nucleon_t",200,0,10);
	TH2D* sPN_t = new TH2D("sPN_t",";t;s",200,0,10,sPN_nBins,sPN_bins);
	TH1D* sPN = new TH1D("sPN","sPN",sPN_nBins,sPN_bins);
	TH1D* sPN_4pt2 = new TH1D("sPN_4pt2","sPN_4pt2",200,0,10);
	TH1D* sPN_Jpsi = new TH1D("sPN_Jpsi","sPN_Jpsi",sPN_nBins,sPN_bins);
	TH1D* sPN_Jpsi_fix = new TH1D("sPN_Jpsi_fix","sPN_Jpsi_fix",sPN_nBins,sPN_bins);
	TH1D* sPN_Jpsi_fix_noMass = new TH1D("sPN_Jpsi_fix_noMass","sPN_Jpsi_fix_noMass",sPN_nBins,sPN_bins);

	TH1D* nk_truth = new TH1D("nk_truth","k (GeV/c)", 100,0,3);
	TH1D* nk_spectator = new TH1D("nk_spectator",";k (GeV/c)", 100, 0,3);
	TH1D* nk_spectator_pt = new TH1D("nk_spectator_pt",";k (GeV/c)", 100, 0,3);
	TH1D* nk_allfinalstate = new TH1D("nk_allfinalstate",";k (GeV/c)", 100, 0,3);
	TH1D* d_k = new TH1D("d_k","d_k",100,-1,1);

	TChain *tree = new TChain("EICTree");
	tree->Add("/eicdata/eic0003/ztu/BeAGLE_devK/"+filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	double energy_resolution = 0.5;//50%
	TF1* smear_e = new TF1("smear_e","gaus(0)",-30,30);
	smear_e->SetParameter(0,1);
	smear_e->SetParameter(1,0);
	smear_e->SetParameter(2, sqrt( (energy_resolution/sqrt(135.))*(energy_resolution/sqrt(135.)) + 0.05*0.05 )*135. );

	TF1* smear_theta = new TF1("smear_theta","gaus(0)",-0.001,0.001);
	smear_theta->SetParameter(0,1);
	smear_theta->SetParameter(1,0);
	smear_theta->SetParameter(2,0.00025);

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		double Atarg = event->Atarg;
		double pztarg_total = pztarg*Atarg;
		double pxf = event->pxf;
		double pyf = event->pyf;
		double pzf = event->pzf;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep));//neglecting e mass
		TLorentzVector d_beam(0.,0.,pztarg_total,sqrt(pztarg_total*pztarg_total+MASS_DEUTERON*MASS_DEUTERON));
		TLorentzVector e_scattered(0.,0.,0.,0.);

		smear_e->SetParameter(2, sqrt( (energy_resolution/sqrt(pztarg))*(energy_resolution/sqrt(pztarg)) + 0.05*0.05 )*pztarg );

		TVector3 b = d_beam.BoostVector();

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

		nk_truth->Fill( sqrt(pxf*pxf+pyf*pyf+pzf*pzf) );

		that->Fill( fabs(t_hat) );

		int nParticles_process = 0;
		TLorentzVector n_4vect_unsmear;
		TLorentzVector p_4vect, n_4vect,j_4vect,q;
		TLorentzVector p_4vect_irf, n_4vect_irf,j_4vect_irf,q_irf,d_beam_irf;

		d_beam_irf = d_beam;

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
			if(pdg == 2112) {n_4vect = ppart; n_4vect_unsmear = n_4vect;}

			/*
			- do energy and scattering angle smearing 
			- together with acceptance cuts
			*/
			
			if(doSmear_){
				double angle = ppart.Theta();
				//acceptance cuts for proton
				if( pdg == 2212 ){
					if( (angle > 0.005 && angle < 0.007) || angle > 0.022 ) p_4vect.SetPxPyPzE(0.,0.,0.,0.);
				}
				//acceptance cuts for neutron
				if( pdg == 2112 ){
					if( angle > 0.004 ) {
						n_4vect.SetPxPyPzE(0.,0.,0.,0.);
					}
					else{
						//smearing neutron
						double E_n = ppart.E();
						double delta_E = smear_e->GetRandom();
						E_n = E_n + delta_E;
						double delta_Theta = 0.;
						angle = angle + delta_Theta;
						double Pz_n2 = (E_n*E_n - MASS_NEUTRON*MASS_NEUTRON)/(1+TMath::Sin(angle)*TMath::Sin(angle));
						double Pz_n = sqrt(Pz_n2);
						double Pt_n2 = (E_n*E_n - MASS_NEUTRON*MASS_NEUTRON - Pz_n2);
						double Pt_n = sqrt(Pt_n2);
						double Px_n = Pt_n*TMath::Cos(ppart.Phi());
						double Py_n = Pt_n*TMath::Sin(ppart.Phi());

						n_4vect.SetPxPyPzE(Px_n, Py_n, Pz_n, E_n);
					}
				}
			}

			if(pdg == 443 ) j_4vect_irf = j_4vect;
			if(pdg == 2212) p_4vect_irf = p_4vect;
			if(pdg == 2112) n_4vect_irf = n_4vect; 

		
			nParticles_process++;

		} // end of particle loop

		if( p_4vect.E() == 0 || n_4vect.E() == 0 ) continue;

		//boost
		j_4vect_irf.Boost(-b);
		p_4vect_irf.Boost(-b);
		n_4vect_irf.Boost(-b);
		q_irf.Boost(-b);
		d_beam_irf.Boost(-b);

		double pt2 = j_4vect.Pt()*j_4vect.Pt();
		tjpsi->Fill( pt2 );
		h_trk->Fill( nParticles_process );
		nRes->Fill( n_4vect_unsmear.E()-n_4vect.E(), n_4vect_unsmear.Theta()-n_4vect.Theta() );
		
		if( struckproton ){
			TLorentzVector n_partner_4vect_irf;
			n_partner_4vect_irf.SetPxPyPzE(-n_4vect_irf.Px(), -n_4vect_irf.Py(), -n_4vect_irf.Pz(), sqrt(n_4vect_irf.P()*n_4vect_irf.P()+MASS_PROTON*MASS_PROTON) );
			//intrinsic Beagle n(k) distribution -> SNN distribution
			sPN->Fill( (n_partner_4vect_irf+n_4vect_irf).Mag2() );
			//approximation by spectator nucleon pt in the lab frame
			sPN_4pt2->Fill( 4*n_4vect.Pt()*n_4vect.Pt() );

			//use spectator only:
			nk_spectator->Fill( n_4vect_irf.P() );
			//use pt of the lab frame particle to determine k
			TLorentzVector n;
			n.SetPtEtaPhiM( n_4vect.Pt(), n_4vect.Eta(), n_4vect.Phi(), MASS_NEUTRON);
			n.Boost(-b);
			nk_spectator_pt->Fill( n.P() );
			
		} 
		else{
			TLorentzVector p_partner_4vect_irf;
			p_partner_4vect_irf.SetPxPyPzE(-p_4vect_irf.Px(), -p_4vect_irf.Py(), -p_4vect_irf.Pz(), sqrt(p_4vect_irf.P()*p_4vect_irf.P()+MASS_NEUTRON*MASS_NEUTRON) );
			sPN->Fill( (p_partner_4vect_irf+p_4vect_irf).Mag2() );
			sPN_4pt2->Fill( 4*p_4vect.Pt()*p_4vect.Pt() );

			nk_spectator->Fill( p_4vect_irf.P() );
			//use pt of the lab frame particle to determine k
			TLorentzVector p;
			p.SetPtEtaPhiM( p_4vect.Pt(), p_4vect.Eta(), p_4vect.Phi(), MASS_PROTON);
			p.Boost(-b);
			nk_spectator_pt->Fill( p.P() );

		}

		//inclusive J/psi measurement, convolution of exp and intrinsic n(k)
		sPN_Jpsi->Fill( (p_4vect_irf+n_4vect_irf).Mag2() );
		//remove dependence of Q2 and Jpsi production
		sPN_Jpsi_fix->Fill( (p_4vect_irf+n_4vect_irf+j_4vect_irf-q_irf).Mag2() );
		nucleon_t->Fill( (p_4vect_irf+n_4vect_irf - d_beam_irf).Mag2() );
		sPN_t->Fill((p_4vect_irf+n_4vect_irf - d_beam_irf).Mag2(), (p_4vect_irf+n_4vect_irf+j_4vect_irf-q_irf).Mag2());
		//use all final state particles:
		TLorentzVector pn = p_4vect_irf+n_4vect_irf+j_4vect_irf-q_irf;
		double Epn = pn.E();
		double k = sqrt( Epn*Epn/4. - MASS_NEUTRON*MASS_NEUTRON );//use proton mass to simplify
		nk_allfinalstate->Fill( k );
		d_k->Fill( pn.P() ); //cross check with net zero momentum in the IRF.

	}

	output->Write();
	output->Close();




}