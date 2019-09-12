#include "hist.h"//define all the histograms
#include "PRINT4VECTOR.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include <sstream>
#include <string>

using namespace std;
using namespace erhic;

#define MASS_MUON  0.1056

double sPN_bins[]={0.,1.0,2.0,3.0,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.7,5.0,5.5,6.0,7.0,8.0,9.0,10.0,12.0,15.0};
int sPN_nBins = sizeof(sPN_bins)/sizeof(sPN_bins[0]) -1;

double nk_bins[]={0.,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
					0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,
					0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,
					0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,
					0.42,0.44,0.46,0.48,0.5,
					0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1.0};

int nk_nBins = sizeof(nk_bins)/sizeof(nk_bins[0]) -1;

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

//mathematica one of the two solutions are physical.
Double_t getCorrJz(Double_t qzkz, Double_t numn, Double_t jx, Double_t jy, Double_t px, Double_t py, Double_t Mp){

	double Md = MASS_DEUTERON;
	double Mj = MASS_JPSI;

	double finalJz = (qzkz*(TMath::Power(jx,2) + TMath::Power(jy,2) + TMath::Power(Mj,2) - TMath::Power(Mp,2) + TMath::Power(Md + numn,2) - 
        TMath::Power(px,2) - TMath::Power(py,2) - TMath::Power(qzkz,2)) - 
     sqrt(TMath::Power(Md + numn,2)*(TMath::Power(jx,4) + TMath::Power(jy,4) + TMath::Power(Md,4) - 
         2*TMath::Power(Md,2)*TMath::Power(Mj,2) + TMath::Power(Mj,4) - 2*TMath::Power(Md,2)*TMath::Power(Mp,2) - 
         2*TMath::Power(Mj,2)*TMath::Power(Mp,2) + TMath::Power(Mp,4) + 4*TMath::Power(Md,3)*numn - 
         4*Md*TMath::Power(Mj,2)*numn - 4*Md*TMath::Power(Mp,2)*numn + 
         6*TMath::Power(Md,2)*TMath::Power(numn,2) - 2*TMath::Power(Mj,2)*TMath::Power(numn,2) - 
         2*TMath::Power(Mp,2)*TMath::Power(numn,2) + 4*Md*TMath::Power(numn,3) + TMath::Power(numn,4) - 
         2*TMath::Power(Md,2)*TMath::Power(px,2) - 2*TMath::Power(Mj,2)*TMath::Power(px,2) + 
         2*TMath::Power(Mp,2)*TMath::Power(px,2) - 4*Md*numn*TMath::Power(px,2) - 
         2*TMath::Power(numn,2)*TMath::Power(px,2) + TMath::Power(px,4) - 2*TMath::Power(Md,2)*TMath::Power(py,2) - 
         2*TMath::Power(Mj,2)*TMath::Power(py,2) + 2*TMath::Power(Mp,2)*TMath::Power(py,2) - 
         4*Md*numn*TMath::Power(py,2) - 2*TMath::Power(numn,2)*TMath::Power(py,2) + 
         2*TMath::Power(px,2)*TMath::Power(py,2) + TMath::Power(py,4) + 
         2*(TMath::Power(Mj,2) + TMath::Power(Mp,2) - TMath::Power(Md + numn,2) + TMath::Power(px,2) + 
            TMath::Power(py,2))*TMath::Power(qzkz,2) + TMath::Power(qzkz,4) - 
         2*TMath::Power(jy,2)*(-TMath::Power(Mj,2) + TMath::Power(Mp,2) + TMath::Power(Md + numn,2) + 
            TMath::Power(px,2) + TMath::Power(py,2) - TMath::Power(qzkz,2)) + 
         2*TMath::Power(jx,2)*(TMath::Power(jy,2) + TMath::Power(Mj,2) - TMath::Power(Mp,2) - 
            TMath::Power(Md + numn,2) - TMath::Power(px,2) - TMath::Power(py,2) + TMath::Power(qzkz,2)))))/
   (2.*(Md + numn - qzkz)*(Md + numn + qzkz));

   return finalJz;
}

Double_t getCorrPz(Double_t qzkz, Double_t numn, Double_t jx, Double_t jy, Double_t px, Double_t py, Double_t Mp){

	double Md = MASS_DEUTERON;
	double Mj = MASS_JPSI;

	double finalPz = (qzkz*(-TMath::Power(jx,2) - TMath::Power(jy,2) - TMath::Power(Mj,2) + TMath::Power(Mp,2) + TMath::Power(Md + numn,2) + 
        TMath::Power(px,2) + TMath::Power(py,2) - TMath::Power(qzkz,2)) + 
     sqrt(TMath::Power(Md + numn,2)*(TMath::Power(jx,4) + TMath::Power(jy,4) + TMath::Power(Md,4) - 
         2*TMath::Power(Md,2)*TMath::Power(Mj,2) + TMath::Power(Mj,4) - 2*TMath::Power(Md,2)*TMath::Power(Mp,2) - 
         2*TMath::Power(Mj,2)*TMath::Power(Mp,2) + TMath::Power(Mp,4) + 4*TMath::Power(Md,3)*numn - 
         4*Md*TMath::Power(Mj,2)*numn - 4*Md*TMath::Power(Mp,2)*numn + 
         6*TMath::Power(Md,2)*TMath::Power(numn,2) - 2*TMath::Power(Mj,2)*TMath::Power(numn,2) - 
         2*TMath::Power(Mp,2)*TMath::Power(numn,2) + 4*Md*TMath::Power(numn,3) + TMath::Power(numn,4) - 
         2*TMath::Power(Md,2)*TMath::Power(px,2) - 2*TMath::Power(Mj,2)*TMath::Power(px,2) + 
         2*TMath::Power(Mp,2)*TMath::Power(px,2) - 4*Md*numn*TMath::Power(px,2) - 
         2*TMath::Power(numn,2)*TMath::Power(px,2) + TMath::Power(px,4) - 2*TMath::Power(Md,2)*TMath::Power(py,2) - 
         2*TMath::Power(Mj,2)*TMath::Power(py,2) + 2*TMath::Power(Mp,2)*TMath::Power(py,2) - 
         4*Md*numn*TMath::Power(py,2) - 2*TMath::Power(numn,2)*TMath::Power(py,2) + 
         2*TMath::Power(px,2)*TMath::Power(py,2) + TMath::Power(py,4) + 
         2*(TMath::Power(Mj,2) + TMath::Power(Mp,2) - TMath::Power(Md + numn,2) + TMath::Power(px,2) + 
            TMath::Power(py,2))*TMath::Power(qzkz,2) + TMath::Power(qzkz,4) - 
         2*TMath::Power(jy,2)*(-TMath::Power(Mj,2) + TMath::Power(Mp,2) + TMath::Power(Md + numn,2) + 
            TMath::Power(px,2) + TMath::Power(py,2) - TMath::Power(qzkz,2)) + 
         2*TMath::Power(jx,2)*(TMath::Power(jy,2) + TMath::Power(Mj,2) - TMath::Power(Mp,2) - 
            TMath::Power(Md + numn,2) - TMath::Power(px,2) - TMath::Power(py,2) + TMath::Power(qzkz,2)))))/
   (2.*(Md + numn - qzkz)*(Md + numn + qzkz));

   return finalPz;
}

//solution 2 using lightcone kinematics
Double_t getCorrPzLF(Double_t Ennz, Double_t Ennz2, Double_t nuqzmd, Double_t nuqzmd2, Double_t jx, Double_t jy, Double_t px, Double_t py, Double_t Mp){

	double Mj = MASS_JPSI;

	double finalPz = ((-Ennz + Ennz2 + nuqzmd - nuqzmd2)*
      (-TMath::Power(jx,2) - TMath::Power(jy,2) - TMath::Power(Mj,2) + TMath::Power(Mp,2) + 
        (Ennz - nuqzmd)*(Ennz2 - nuqzmd2) + TMath::Power(px,2) + TMath::Power(py,2)) + 
     sqrt(TMath::Power(Ennz + Ennz2 - nuqzmd - nuqzmd2,2)*
       (TMath::Power(jx,4) + TMath::Power(jy,4) + 2*TMath::Power(jy,2)*TMath::Power(Mj,2) + TMath::Power(Mj,4) - 
         2*TMath::Power(jy,2)*TMath::Power(Mp,2) - 2*TMath::Power(Mj,2)*TMath::Power(Mp,2) + TMath::Power(Mp,4) + 
         2*Ennz2*TMath::Power(jy,2)*nuqzmd + 2*Ennz2*TMath::Power(Mj,2)*nuqzmd + 
         2*Ennz2*TMath::Power(Mp,2)*nuqzmd + TMath::Power(Ennz2,2)*TMath::Power(nuqzmd,2) + 
         TMath::Power(Ennz,2)*TMath::Power(Ennz2 - nuqzmd2,2) - 2*TMath::Power(jy,2)*nuqzmd*nuqzmd2 - 
         2*TMath::Power(Mj,2)*nuqzmd*nuqzmd2 - 2*TMath::Power(Mp,2)*nuqzmd*nuqzmd2 - 
         2*Ennz2*TMath::Power(nuqzmd,2)*nuqzmd2 + TMath::Power(nuqzmd,2)*TMath::Power(nuqzmd2,2) - 
         2*TMath::Power(jy,2)*TMath::Power(px,2) - 2*TMath::Power(Mj,2)*TMath::Power(px,2) + 
         2*TMath::Power(Mp,2)*TMath::Power(px,2) + 2*Ennz2*nuqzmd*TMath::Power(px,2) - 
         2*nuqzmd*nuqzmd2*TMath::Power(px,2) + TMath::Power(px,4) + 
         2*(-TMath::Power(jy,2) - TMath::Power(Mj,2) + TMath::Power(Mp,2) + Ennz2*nuqzmd - 
            nuqzmd*nuqzmd2 + TMath::Power(px,2))*TMath::Power(py,2) + TMath::Power(py,4) + 
         2*TMath::Power(jx,2)*(TMath::Power(jy,2) + TMath::Power(Mj,2) - TMath::Power(Mp,2) + Ennz2*nuqzmd - 
            nuqzmd*nuqzmd2 - TMath::Power(px,2) - TMath::Power(py,2)) - 
         2*Ennz*(Ennz2 - nuqzmd2)*(TMath::Power(jx,2) + TMath::Power(jy,2) + TMath::Power(Mj,2) + 
            TMath::Power(Mp,2) + Ennz2*nuqzmd - nuqzmd*nuqzmd2 + TMath::Power(px,2) + TMath::Power(py,2)))
       ))/(4.*(Ennz - nuqzmd)*(Ennz2 - nuqzmd2));

	return finalPz;
}

Double_t getCorrJzLF(Double_t Ennz, Double_t Ennz2, Double_t nuqzmd, Double_t nuqzmd2, Double_t jx, Double_t jy, Double_t px, Double_t py, Double_t Mp){

	double Mj = MASS_JPSI;

	double finalJz = -((Ennz - Ennz2 - nuqzmd + nuqzmd2)*
       (TMath::Power(jx,2) + TMath::Power(jy,2) + TMath::Power(Mj,2) - TMath::Power(Mp,2) + 
         (Ennz - nuqzmd)*(Ennz2 - nuqzmd2) - TMath::Power(px,2) - TMath::Power(py,2)) + 
      sqrt(TMath::Power(Ennz + Ennz2 - nuqzmd - nuqzmd2,2)*
        (TMath::Power(jx,4) + TMath::Power(jy,4) + 2*TMath::Power(jy,2)*TMath::Power(Mj,2) + TMath::Power(Mj,4) - 
          2*TMath::Power(jy,2)*TMath::Power(Mp,2) - 2*TMath::Power(Mj,2)*TMath::Power(Mp,2) + TMath::Power(Mp,4) + 
          2*Ennz2*TMath::Power(jy,2)*nuqzmd + 2*Ennz2*TMath::Power(Mj,2)*nuqzmd + 
          2*Ennz2*TMath::Power(Mp,2)*nuqzmd + TMath::Power(Ennz2,2)*TMath::Power(nuqzmd,2) + 
          TMath::Power(Ennz,2)*TMath::Power(Ennz2 - nuqzmd2,2) - 2*TMath::Power(jy,2)*nuqzmd*nuqzmd2 - 
          2*TMath::Power(Mj,2)*nuqzmd*nuqzmd2 - 2*TMath::Power(Mp,2)*nuqzmd*nuqzmd2 - 
          2*Ennz2*TMath::Power(nuqzmd,2)*nuqzmd2 + TMath::Power(nuqzmd,2)*TMath::Power(nuqzmd2,2) - 
          2*TMath::Power(jy,2)*TMath::Power(px,2) - 2*TMath::Power(Mj,2)*TMath::Power(px,2) + 
          2*TMath::Power(Mp,2)*TMath::Power(px,2) + 2*Ennz2*nuqzmd*TMath::Power(px,2) - 
          2*nuqzmd*nuqzmd2*TMath::Power(px,2) + TMath::Power(px,4) + 
          2*(-TMath::Power(jy,2) - TMath::Power(Mj,2) + TMath::Power(Mp,2) + Ennz2*nuqzmd - 
             nuqzmd*nuqzmd2 + TMath::Power(px,2))*TMath::Power(py,2) + TMath::Power(py,4) + 
          2*TMath::Power(jx,2)*(TMath::Power(jy,2) + TMath::Power(Mj,2) - TMath::Power(Mp,2) + Ennz2*nuqzmd - 
             nuqzmd*nuqzmd2 - TMath::Power(px,2) - TMath::Power(py,2)) - 
          2*Ennz*(Ennz2 - nuqzmd2)*
           (TMath::Power(jx,2) + TMath::Power(jy,2) + TMath::Power(Mj,2) + TMath::Power(Mp,2) + Ennz2*nuqzmd - 
             nuqzmd*nuqzmd2 + TMath::Power(px,2) + TMath::Power(py,2)))))/
   (4.*(Ennz - nuqzmd)*(Ennz2 - nuqzmd2));

	return finalJz;
}

void eD_SRC_main(const int nEvents = 40000, TString filename="", const bool doSmear_ = false, const bool doAcceptance_ = false, const double rZDC = 1.){

	std::ostringstream os;
	os << (int) doSmear_;
	os << (int) doAcceptance_;
	os << "_ZDC_" << (double) rZDC;
	std::string str = os.str();
	TString settings = (TString) str;

	TFile * output = new TFile("../rootfiles/"+filename+"_"+settings+"_main_Beagle.root","recreate");
	
	TH1D* h_trk = new TH1D("h_trk","h_trk",50,0,50);
	TH1D* that = new TH1D("that","that",200,0,10);
	TH1D* tjpsi = new TH1D("tjpsi","tjpsi",200,0,10);
	TH2D* nRes = new TH2D("nRes","",60,-30,30,20,-0.1,0.1);
	TH1D* nucleon_t = new TH1D("nucleon_t","nucleon_t",200,0,10);
	TH2D* sPN_t = new TH2D("sPN_t",";t;s",200,0,10,sPN_nBins,sPN_bins);
	TH2D* sPN_k = new TH2D("sPN_k",";t;s",200,0,1,sPN_nBins,sPN_bins);
	TH1D* sPN = new TH1D("sPN","sPN",sPN_nBins,sPN_bins);
	TH1D* sPN_4pt2 = new TH1D("sPN_4pt2","sPN_4pt2",200,0,10);
	TH1D* sPN_Jpsi_fix = new TH1D("sPN_Jpsi_fix","sPN_Jpsi_fix",sPN_nBins,sPN_bins);

	TH1D* nk_truth = new TH1D("nk_truth","k (GeV/c)", nk_nBins, nk_bins);
	TH1D* nk_spectator = new TH1D("nk_spectator",";k (GeV/c)", nk_nBins, nk_bins);
	TH2D* EvsPz = new TH2D("EvsPz",";pz;E",500,-0.01,0.01,500,-0.01,0.01);
	TH2D* EvsPzFix = new TH2D("EvsPzFix",";pz;E",500,-0.01,0.01,500,-0.01,0.01);
	TH1D* Pp_old = new TH1D("Pp_old","",500,0,5);
	TH1D* Pp_new = new TH1D("Pp_new","",500,0,5);
	TH1D* Pp_new1 = new TH1D("Pp_new1","",500,0,5);
	TH1D* Pp_new2 = new TH1D("Pp_new2","",500,0,5);


	TChain *tree = new TChain("EICTree");
	tree->Add("/eicdata/eic0003/ztu/BeAGLE_devK/"+filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	double energy_resolution = rZDC;//50%
	TF1* smear_e = new TF1("smear_e","gaus(0)",-30,30);
	smear_e->SetParameter(0,1);
	smear_e->SetParameter(1,0);
	smear_e->SetParameter(2, sqrt( (energy_resolution/sqrt(135.))*(energy_resolution/sqrt(135.)) + 0.05*0.05 )*135. );

	TF1* smear_theta = new TF1("smear_theta","gaus(0)",-0.001,0.001);
	smear_theta->SetParameter(0,1);
	smear_theta->SetParameter(1,0);
	smear_theta->SetParameter(2,0.00025);//assume 8x8 0.2m x 0.2m ZDC 27meter away from IR at eRHIC. 
										//resolution is smallest distance/sqrt(12) ~ 0.007

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

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep+0.00051*0.00051));
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
		double nk_event = sqrt(pxf*pxf+pyf*pyf+pzf*pzf);
		
		if( event_process != 91 ) continue;
		if( trueQ2 < 1. ) continue;
		if( trueY > 0.85 || trueY < 0.05 ) continue;
		bool struckproton = false;
		if( struck_nucleon == 2212 ) struckproton = true;

		that->Fill( fabs(t_hat) );

		int nParticles_process = 0;
		TLorentzVector n_4vect_unsmear;
		TLorentzVector p_4vect, n_4vect,j_4vect,q;
		TLorentzVector p_4vect_irf, n_4vect_irf,j_4vect_irf,q_irf,d_beam_irf;
		TLorentzVector jnew,pnew,nnew;
		TLorentzVector lfjnew,lfpnew,lfnnew;
		TLorentzVector jnew2,pnew2;

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
			TLorentzVector ppart = particle->Get4Vector();

			if( index == 3 ) {
				e_scattered.SetPtEtaPhiM(pt,eta,phi,0.00051);
				// e_scattered = ppart;
			}
			if( status != 1 ) continue;

			//photon 4vector
			q = e_beam-e_scattered;
			q_irf = q;
			
			
			if(pdg == 443 ) j_4vect = ppart;
			if(pdg == 2212) p_4vect = ppart;
			if(pdg == 2112) {n_4vect = ppart; n_4vect_unsmear = n_4vect;}

			/*
			- do energy and scattering angle smearing 
			- together with acceptance cuts
			*/
			
			if( doAcceptance_ ){
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
				}
				if( doSmear_ ){
					//smearing neutron
					double E_n = ppart.E();
					double delta_E = smear_e->GetRandom();
					E_n = E_n + delta_E;
					double delta_Theta = smear_theta->GetRandom();
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

			if(pdg == 443 ) j_4vect_irf = j_4vect;
			if(pdg == 2212) p_4vect_irf = p_4vect;
			if(pdg == 2112) n_4vect_irf = n_4vect; 

		
			nParticles_process++;

		} // end of particle loop

		nk_truth->Fill( sqrt(pxf*pxf+pyf*pyf+pzf*pzf) );

		if( p_4vect.E() == 0 || n_4vect.E() == 0 ) continue;

		cout << " Event ~ " << i << endl;

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
		

		TLorentzVector struck_4vect_irf, spectator_4vect_irf;
		double struck_mass = MASS_PROTON;
		double spectator_mass = MASS_NEUTRON;
		if( struckproton ) {
			struck_4vect_irf = p_4vect_irf;
			spectator_4vect_irf = n_4vect_irf;
			struck_mass = MASS_PROTON;
			spectator_mass = MASS_NEUTRON;
		}
		else{
			struck_4vect_irf = n_4vect_irf;
			spectator_4vect_irf = p_4vect_irf;
			struck_mass = MASS_NEUTRON;
			spectator_mass = MASS_PROTON;
		}

		//use spectator only:
		nk_spectator->Fill( spectator_4vect_irf.P() );
		
		/*
		fixing deuteron momentum nonconservation:
		*/

		TLorentzVector testp = q_irf+d_beam_irf-j_4vect_irf-struck_4vect_irf-spectator_4vect_irf;

		//approach 1
		double qzkz = q_irf.Pz() - (spectator_4vect_irf.Pz());//qz-kz
		double numn = q_irf.E() - spectator_4vect_irf.E();//sqrt( MASS_NEUTRON*MASS_NEUTRON + pxf*pxf+pyf*pyf+pzf*pzf )
		double jx = j_4vect_irf.Px();
		double jy = j_4vect_irf.Py();
		double px = struck_4vect_irf.Px();
		double py = struck_4vect_irf.Py();

		double jz = getCorrJz(qzkz,numn,jx,jy,px,py,struck_mass);
		double pz = getCorrPz(qzkz,numn,jx,jy,px,py,struck_mass);

		double px_new = struck_4vect_irf.Px();
		double py_new = struck_4vect_irf.Py();
		double pz_new = pz;
		pnew.SetPxPyPzE(px_new,py_new,pz_new, sqrt( struck_mass*struck_mass + px_new*px_new + py_new*py_new + pz_new*pz_new));
		
		double jx_new = j_4vect_irf.Px();
		double jy_new = j_4vect_irf.Py();
		double jz_new = jz;
		jnew.SetPxPyPzE(jx_new,jy_new,jz_new, sqrt( MASS_JPSI*MASS_JPSI + jx_new*jx_new + jy_new*jy_new + jz_new*jz_new));
	
		//approach 2
		double Ennz = spectator_4vect_irf.E() + spectator_4vect_irf.Pz();
		double Ennz2 = spectator_4vect_irf.E() - spectator_4vect_irf.Pz();
		double nuqzmd = q_irf.E()+q_irf.Pz()+MASS_DEUTERON;
		double nuqzmd2 = q_irf.E()-q_irf.Pz()+MASS_DEUTERON;
		jx = j_4vect_irf.Px()+spectator_4vect_irf.Px();
		jy = j_4vect_irf.Py()+spectator_4vect_irf.Py();
		px = struck_4vect_irf.Px()-spectator_4vect_irf.Px();
		py = struck_4vect_irf.Py()-spectator_4vect_irf.Py();

		double lfpz = getCorrPzLF(Ennz,Ennz2,nuqzmd,nuqzmd2,jx,jy,px,py,struck_mass);
		double lfjz = getCorrJzLF(Ennz,Ennz2,nuqzmd,nuqzmd2,jx,jy,px,py,struck_mass);

		px_new = struck_4vect_irf.Px()-spectator_4vect_irf.Px();
		py_new = struck_4vect_irf.Py()-spectator_4vect_irf.Py();
		pz_new = lfpz;
		lfpnew.SetPxPyPzE(px_new,py_new,pz_new, sqrt( struck_mass*struck_mass + px_new*px_new + py_new*py_new + pz_new*pz_new));
		
		jx_new = j_4vect_irf.Px()+spectator_4vect_irf.Px();
		jy_new = j_4vect_irf.Py()+spectator_4vect_irf.Py();
		jz_new = lfjz;
		lfjnew.SetPxPyPzE(jx_new,jy_new,jz_new, sqrt( MASS_JPSI*MASS_JPSI + jx_new*jx_new + jy_new*jy_new + jz_new*jz_new));

		TLorentzVector testnew2 = q_irf+d_beam_irf-lfjnew-lfpnew-spectator_4vect_irf;
		cout << "check momentum conservation approach 2, total change q+d-j-p'-n' should be 0 now: " << endl;
		PRINT4VECTOR(testnew2,1);

		EvsPz->Fill(testp.Pz(), testp.E());
		EvsPzFix->Fill(testnew2.Pz(), testnew2.E());

		//approach 3

		/*
		- Start trying off-shell intermediate conditions
		- Assume same proton and neutron mass. 
		*/
		double kmag = spectator_4vect_irf.P();
		double MnuclOff = sqrt(0.25*MASS_DEUTERON*MASS_DEUTERON - kmag*kmag);
		double PpOff = sqrt( struck_mass*struck_mass - MnuclOff*MnuclOff + struck_4vect_irf.P()*struck_4vect_irf.P() );
		double Ptheta = struck_4vect_irf.Theta();
		double Pnewpt = PpOff*TMath::Sin(Ptheta);
		TLorentzVector Poff4vector;Poff4vector.SetPtEtaPhiM(Pnewpt,struck_4vect_irf.Eta(),struck_4vect_irf.Phi(),MnuclOff);
		TLorentzVector Pon4vectorNew; 
		double Ppx = Poff4vector.Px() - spectator_4vect_irf.Px();
		double Ppy = Poff4vector.Py() - spectator_4vect_irf.Py();
		double Ppz = Poff4vector.Pz() - spectator_4vect_irf.Pz();
		Pon4vectorNew.SetPxPyPzE(Ppx,Ppy,Ppz,sqrt(Ppx*Ppx+Ppy*Ppy+Ppz*Ppz+struck_mass*struck_mass) );

		qzkz = q_irf.Pz() - (spectator_4vect_irf.Pz());//qz-kz
		numn = q_irf.E() - spectator_4vect_irf.E();//sqrt( MASS_NEUTRON*MASS_NEUTRON + pxf*pxf+pyf*pyf+pzf*pzf )
		jx = j_4vect_irf.Px()+spectator_4vect_irf.Px()-(Poff4vector.Px()-struck_4vect_irf.Px());
		jy = j_4vect_irf.Py()+spectator_4vect_irf.Py()-(Poff4vector.Py()-struck_4vect_irf.Py());
		px = Pon4vectorNew.Px();
		py = Pon4vectorNew.Py();

		jz = getCorrJz(qzkz,numn,jx,jy,px,py,struck_mass);
		pz = getCorrPz(qzkz,numn,jx,jy,px,py,struck_mass);

		px_new = px;
		py_new = py;
		pz_new = pz;
		pnew2.SetPxPyPzE(px_new,py_new,pz_new, sqrt( struck_mass*struck_mass + px_new*px_new + py_new*py_new + pz_new*pz_new));
		
		jx_new = jx;
		jy_new = jy;
		jz_new = jz;
		jnew2.SetPxPyPzE(jx_new,jy_new,jz_new, sqrt( MASS_JPSI*MASS_JPSI + jx_new*jx_new + jy_new*jy_new + jz_new*jz_new));

		TLorentzVector testnew3 = q_irf+d_beam_irf-jnew2-pnew2-spectator_4vect_irf;
		cout << "check momentum conservation approach 3, total change q+d-j-p'-n' should be 0 now: " << endl;
		PRINT4VECTOR(testnew3,1);

		cout << "Let's compare different kinematics method:" << endl;
		cout << "proton old"<<endl;
		PRINT4VECTOR(struck_4vect_irf,1);
		cout << "proton new"<<endl;
		PRINT4VECTOR(pnew,1);
		cout << "proton new lf"<<endl;
		PRINT4VECTOR(lfpnew,1);
		cout << "proton new 2"<<endl;
		PRINT4VECTOR(pnew2,1);
		cout << "jpsi old"<<endl;
		PRINT4VECTOR(j_4vect_irf,1);
		cout << "jpsi new"<<endl;
		PRINT4VECTOR(jnew,1);
		cout << "jpsi new lf"<<endl;
		PRINT4VECTOR(lfjnew,1);
		cout << "jpsi new 2"<<endl;
		PRINT4VECTOR(jnew2,1);

		Pp_old->Fill( struck_4vect_irf.P() );
		Pp_new->Fill( pnew.P() );
		Pp_new1->Fill( lfpnew.P() );
		Pp_new2->Fill( pnew2.P() );
		
		sPN_Jpsi_fix->Fill( (pnew2+spectator_4vect_irf+jnew2-q_irf).Mag2() );
		nucleon_t->Fill( (pnew2+spectator_4vect_irf - d_beam_irf).Mag2() );
		sPN_t->Fill((pnew2+spectator_4vect_irf - d_beam_irf).Mag2(), (pnew2+spectator_4vect_irf+jnew2-q_irf).Mag2());
		sPN_k->Fill(nk_event, (pnew2+spectator_4vect_irf+jnew2-q_irf).Mag2());

		/*This wouldn't work because of the offshell mass*/
		// double Epn = pn.E();
		// double EpnRed2 = Epn*Epn - MASS_NEUTRON*MASS_NEUTRON - MASS_PROTON*MASS_PROTON; 
		// double k = sqrt( Epn*Epn/4. - MASS_NEUTRON*MASS_NEUTRON );//use proton mass to simplify
		
	}

	output->Write();
	output->Close();




}