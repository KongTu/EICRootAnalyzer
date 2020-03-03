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

//only for ZDC acceptance
double acceptanceGlobal = 0.005;

// solutions for momentum non-conservations, 
// won't be needed for the next version of BeAGLE
// mathematica gives only one of the two solutions are physical.

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

//for both neutron and proton acceptances
bool passDetector(TLorentzVector p, TVector3 b){

	/*
	- do acceptance cuts -
	*/
	
	bool pass = true;

	p.Boost(b);//boost to lab frame;

	bool isNeutron = true;
	if( p.M() < MASS_PROTON+0.0001 ) isNeutron = false;

	if( isNeutron ){
		if( p.Theta() > acceptanceGlobal ) pass = false;
	}
	else{
		if( (p.Theta() > 0.005 && p.Theta() < 0.007) || p.Theta() > 0.022 ) pass = false;
	}

	return pass;
}

//for neutron energy/position resolutions;
TLorentzVector afterNeutronDetector(TLorentzVector p, TVector3 b, TF1*smear_e_zdc, TF1*smear_theta_zdc){

	TLorentzVector pafter;
	if( p.E() == 0. ) {
		return p;
	}
	//boost to lab frame;
	p.Boost(b);

	//smearing neutron
	double E_n = p.E();
	E_n = E_n*(1+smear_e_zdc->GetRandom());
	double delta_Theta = smear_theta_zdc->GetRandom();
	double angle = p.Theta() + delta_Theta;
	double Pp = sqrt(E_n*E_n - MASS_NEUTRON*MASS_NEUTRON);
	double Pz_n = Pp*TMath::Cos(angle);
	double Px_n = Pp*TMath::Sin(angle)*TMath::Cos(p.Phi());
	double Py_n = Pp*TMath::Sin(angle)*TMath::Sin(p.Phi());

	pafter.SetPxPyPzE(Px_n, Py_n, Pz_n, E_n);

	//boost back to IRF;
	pafter.Boost(-b);
	return pafter;
}
//for proton pt resolutions;
TLorentzVector afterProtonDetector(TLorentzVector p, TVector3 b,TF1*smear_pt_proton){

	TLorentzVector pafter;
	if( p.E() == 0. ) {
		return p;
	}
	//boost to lab frame;
	p.Boost(b);

	double pt = p.Pt();
	double eta = p.Eta();
	double phi = p.Phi();
	double Mass = p.M();
	double Pp = p.P();
	double angle = p.Theta();

	pt = pt*(1+smear_pt_proton->GetRandom());
	pafter.SetPtEtaPhiM(pt,eta,phi,MASS_PROTON);

	//boost back to IRF;
	pafter.Boost(-b);
	return pafter;
}

void eD_SRC_main(const int nEvents = 40000, TString filename="", const int hitNucleon_ = 0, const bool doSmear_ = false, const bool doAcceptance_ = false, const double rZDC = 0.5, const double acceptance=0.005, const double ptreso_ = 0.03){

	//just naming in the output file, only show ZDC parameters. 
	acceptanceGlobal = acceptance;
	std::ostringstream os;
	os << "hitNucleon_" << (int) hitNucleon_;
	os << "_dosmear_" <<(int) doSmear_;
	os << "_doaccept_" << (int) doAcceptance_;
	os << "_ZDCreso_" << (double) rZDC;
	os << "_ZDCaccept_" << (double) acceptance;
	os << "_RPreso_" << (double) ptreso_;
	std::string str = os.str();
	TString settings = (TString) str;

	//not so important now
	TFile* input = new TFile("./inputSd.root","READ");
	TH1D* h_spectral_pt_input = (TH1D*) input->Get("h_spectral_pt");
	h_spectral_pt_input->Scale(1./h_spectral_pt_input->Integral());

	//input from BeAGLE root files
	TFile * output = new TFile("../rootfiles/"+filename+"_"+settings+"_main_Beagle.root","recreate");
	
	//histograms to be saved
	TH1D* nk_truth = new TH1D("nk_truth","k (GeV/c)", nk_nBins, nk_bins);
	TH1D* nk_truth_uniformbins = new TH1D("nk_truth_uniformbins","k (GeV/c)", 200,0,1.4);
	TH2D* h_ThetaVsEnergy_Spectator = new TH2D("h_ThetaVsEnergy_Spectator",";E_{spectator} (GeV);#theta",300,0,200,200,0,100);
	TH2D* h_ThetaVsEnergy_Struck = new TH2D("h_ThetaVsEnergy_Struck",";E_{struck} (GeV);#theta",300,0,200,200,0,100);
	TH1D* sPN = new TH1D(Form("sPN"),"sPN",sPN_nBins,sPN_bins);
	TH1D* Pt_struck = new TH1D("Pt_struck",";p_{T} (GeV)",200,0,1.4);
	TH1D* Pz_struck = new TH1D("Pz_struck",";p_{z} (GeV)",200,-1,1);
	TH1D* Pp_struck = new TH1D("Pp_struck",";p (GeV)",200,0,1.4);
	TH1D* Pt_spectator = new TH1D("Pt_spectator",";p_{T} (GeV)",200,0,1.4);
	TH1D* Pz_spectator = new TH1D("Pz_spectator",";p_{z} (GeV)",200,-1,1);
	TH1D* Pp_spectator = new TH1D("Pp_spectator",";p (GeV)",200,0,1.4);
	TH1D* Pt_VM = new TH1D("Pt_VM",";p_{T} (GeV)",200,0,10);
	TH1D* Pt2_VM = new TH1D("Pt2_VM",";p^{2}_{T} (GeV)",200,0,5);
	TH1D* Pz_VM = new TH1D("Pz_VM",";p_{z} (GeV)",200,-20,20);
	TH1D* Pp_VM = new TH1D("Pp_VM",";p (GeV)",200,0,10);
	TH1D* alpha_spectator = new TH1D("alpha_spectator",";#alpha_{spec}",100,0,2);
	TH1D* ttprime = new TH1D("ttprime",";-t'(GeV)",100,0,2);
	TH1D* t_eej = new TH1D("t_eej",";-t'(GeV)",100,0,2);
	TH1D* t_nprimeprime = new TH1D("t_nprimeprime",";-t'(GeV)",100,0,2);
	TH1D* t_truth = new TH1D("t_truth",";-t'(GeV)",100,0,2);
	TH2D* t_compare = new TH2D("t_compare",";-t'(GeV)",100,0,2,100,0,2);
	TH2D* h_ttprime_alpha = new TH2D("h_ttprime_alpha",";#alpha_{p};-t'",200,0,2,1000,0,1);
	TH2D* h_dNdAlphadPt2 = new TH2D("h_dNdAlphadPt2",";#alpha_{p};p_{T} (GeV/c)'",500,0,2,1000,0,1);
	TH2D* h_ThetaRprimePm = new TH2D("h_ThetaRprimePm",";#theta_{r'};p_{m} (GeV/c)",200,0,PI,200,0,1.4);
	TH1D* h_spectral_pt = new TH1D("h_spectral_pt",";p_{T} (GeV/c)",500,0,1);
	TH1D* h_spectralAtPole = new TH1D("h_spectralAtPole",";-t' (GeV)^{2}",500,0,0.5);

	TChain *tree = new TChain("EICTree");
	tree->Add("/gpfs02/eic/ztu/BeAGLE/BeAGLE_devK_SRC/"+filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	//ZDC for neutron
	double energy_resolution = rZDC;//default 50%
	double energy_resolution_constant_term = 0.05; //default 5%
	double beam_momentum = 135.; // 135 GeV for Deuteron default
	TF1* smear_e_zdc = new TF1("smear_e_zdc","gaus(0)",-5,5);
	smear_e_zdc->SetParameter(0,1);
	smear_e_zdc->SetParameter(1,0);
	smear_e_zdc->SetParameter(2, sqrt( TMath::Power((energy_resolution/sqrt(beam_momentum)),2) 
		+ TMath::Power(energy_resolution_constant_term,2)) );//giving resolution in percent
	//resolution terms adding in quadrature. 

	TF1* smear_theta_zdc = new TF1("smear_theta_zdc","gaus(0)",-0.001,0.001);
	smear_theta_zdc->SetParameter(0,1);
	smear_theta_zdc->SetParameter(1,0);
	double angle_reso = 3e-6;
	smear_theta_zdc->SetParameter(2,angle_reso);//absolute angles smearing
	//1cm position resolution --> 3 microRad resolution
	//Yuji's Letter of Intent in EIC R&D proposal

	//RP,B0,Ext Sensor for proton
	TF1* smear_pt_proton = new TF1("smear_pt_proton","gaus(0)",-10,10);
	smear_pt_proton->SetParameter(0,1);
	smear_pt_proton->SetParameter(1,0);
	smear_pt_proton->SetParameter(2,ptreso_ );//3% default resolution dpt/pt 

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
		
		//overwrite with the correct beam energy:
		beam_momentum = pztarg;
		if( i = 0 ){//only set it once
			smear_e_zdc->SetParameter(2, sqrt( TMath::Power((energy_resolution/sqrt(beam_momentum)),2) 
		+ TMath::Power(energy_resolution_constant_term,2)) );//giving resolution in percent
		}
		
		
		//boost vector for lab <--> d rest frame
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

		// use hitNucleon_ to choose only hit proton/neutron or mixing
		bool struckproton = false;
		if( struck_nucleon == 2212 ) struckproton = true;
		if( hitNucleon_ == 0){
			if(!struckproton) continue;
		}else if( hitNucleon_ == 1 ){
			if(struckproton) continue;
		}else{
			//otherwise it's mixing of both.
		}
		
		int nParticles_process = 0;
		TLorentzVector p_4vect, n_4vect,j_4vect,q;
		TLorentzVector p_4vect_irf, n_4vect_irf,j_4vect_irf,q_irf,d_beam_irf;
		TLorentzVector jnew,pnew;
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
			q = e_beam - e_scattered;
			q_irf = q;
			
			if(pdg == 443 ) j_4vect = ppart;//jpsi
			if(pdg == 2212) p_4vect = ppart;//proton
			if(pdg == 2112) {n_4vect = ppart;}//neutron

			//prepare for boost later
			if(pdg == 443 ) j_4vect_irf = j_4vect;
			if(pdg == 2212) p_4vect_irf = p_4vect;
			if(pdg == 2112) n_4vect_irf = n_4vect; 

			nParticles_process++;

		} // end of particle loop

		if( nParticles_process != 4 ) continue;
		
		//fill n(k) or dN/dk distribution, but averaged over all direction
		//LFKine tells us pzf is not symmetric in the lab frame
		nk_truth->Fill( nk_event );
		nk_truth_uniformbins->Fill( nk_event );

		//boost to d rest frame
		j_4vect_irf.Boost(-b);
		p_4vect_irf.Boost(-b);
		n_4vect_irf.Boost(-b);
		q_irf.Boost(-b);
		d_beam_irf.Boost(-b);

		//assign who's struck and who's spectator
		TLorentzVector struck_4vect_irf, spectator_4vect_irf;
		TLorentzVector struck_4vect, spectator_4vect;
		double struck_mass = MASS_PROTON;
		double spectator_mass = MASS_NEUTRON;
		if( struckproton ) {
			struck_4vect_irf = p_4vect_irf;
			spectator_4vect_irf = n_4vect_irf;
			struck_mass = MASS_PROTON;
			spectator_mass = MASS_NEUTRON;
			spectator_4vect = n_4vect;
			struck_4vect = p_4vect;
		}
		else{
			struck_4vect_irf = n_4vect_irf;
			spectator_4vect_irf = p_4vect_irf;
			struck_mass = MASS_NEUTRON;
			spectator_mass = MASS_PROTON;
			spectator_4vect = p_4vect;
			struck_4vect = n_4vect;
		}
		
		/*
		fixing deuteron momentum nonconservation:
		*/

		// approach 1
		double qzkz = q_irf.Pz() - (spectator_4vect_irf.Pz());//qz-kz
		double numn = q_irf.E() - spectator_4vect_irf.E();//sqrt( MASS_NEUTRON*MASS_NEUTRON + pxf*pxf+pyf*pyf+pzf*pzf )
		double jx = j_4vect_irf.Px();
		double jy = j_4vect_irf.Py();
		double px = struck_4vect_irf.Px();
		double py = struck_4vect_irf.Py();

		double jz = getCorrJz(qzkz,numn,jx,jy,px,py,struck_mass);
		double pz = getCorrPz(qzkz,numn,jx,jy,px,py,struck_mass);

		double px_new = px;
		double py_new = py;
		double pz_new = pz;
		pnew.SetPxPyPzE(px_new,py_new,pz_new, sqrt( struck_mass*struck_mass + px_new*px_new + py_new*py_new + pz_new*pz_new) );
		
		double jx_new = jx;
		double jy_new = jy;
		double jz_new = jz;
		jnew.SetPxPyPzE(jx_new,jy_new,jz_new, sqrt( MASS_JPSI*MASS_JPSI + jx_new*jx_new + jy_new*jy_new + jz_new*jz_new) );


		/*
		- Start trying off-shell intermediate conditions
		- Assume same proton and neutron mass. 
		*/

		//approach 2 is now removed, can be added back anytime. 


		//accpetance and smearing before filling histograms:
		if( doAcceptance_ ) {
			if( !passDetector(pnew,b) ) pnew.SetPxPyPzE(0,0,0,0);
			if( !passDetector(spectator_4vect_irf,b) ) spectator_4vect_irf.SetPxPyPzE(0,0,0,0);
		}
		if( doSmear_ ){
			if( struckproton ){
				pnew = afterProtonDetector(pnew,b,smear_pt_proton); 
				spectator_4vect_irf = afterNeutronDetector(spectator_4vect_irf,b,smear_e_zdc,smear_theta_zdc);
			}
			else{
				spectator_4vect_irf = afterProtonDetector(spectator_4vect_irf,b,smear_pt_proton); 
				pnew = afterNeutronDetector(pnew,b,smear_e_zdc,smear_theta_zdc);
			}
		}

		//Only both proton and neutron in acceptance are kept
		if( pnew.E() == 0 || spectator_4vect_irf.E() == 0 ) continue;

		spectator_4vect_irf.Boost(b);
		spectator_4vect = spectator_4vect_irf;
		spectator_4vect_irf.Boost(-b);
		//fill lab frame theta vs Energy for struck and spectator
		h_ThetaVsEnergy_Spectator->Fill(spectator_4vect.E(), spectator_4vect.Theta()*1000. );
		TLorentzVector pnew_lab;
		pnew.Boost(b);
		pnew_lab = pnew;
		pnew.Boost(-b);
		h_ThetaVsEnergy_Struck->Fill(pnew_lab.E(), pnew_lab.Theta()*1000. );

		//filling alpha of spectator
		double Pplus = (spectator_4vect_irf.E() + spectator_4vect_irf.Pz()) / sqrt(2);
		double PdPlus = MASS_DEUTERON / sqrt(2);
		double alpha_spec = 2*Pplus / PdPlus;
		double alpha_stru = 2. - alpha_spec;
		alpha_spectator->Fill( alpha_spec );

		//filling t' distribution
		double tt = (spectator_4vect_irf - d_beam_irf).Mag2();
		tt = tt - TMath::Power(pnew.M(),2);
		ttprime->Fill( -tt );
		h_ttprime_alpha->Fill( alpha_spec, -tt );

		//filling t distribution 
		// 0) first to flip the virtual photon 4 vector to be e'-e
		q_irf.Boost(b);
		q_irf = -q_irf;
		q_irf.Boost(-b);
		// 1) (e'-e+Jpsi)**2
		double t1_uppervtx = (q_irf + jnew).Mag2();
		t_eej->Fill( -t1_uppervtx );
		// 2) (p - (n''))**2
		// use LF kinematics to calculate the struck nucleon pz, E before interactions.
		Double_t E_bInt = (alpha_stru*MASS_DEUTERON)/4. + (spectator_4vect_irf.Px()*spectator_4vect_irf.Px()+
			spectator_4vect_irf.Py()*spectator_4vect_irf.Py()+struck_mass*struck_mass)/(alpha_stru*MASS_DEUTERON);
		Double_t Pz_bInt = -(alpha_stru*MASS_DEUTERON)/4. + (spectator_4vect_irf.Px()*spectator_4vect_irf.Px()+
			spectator_4vect_irf.Py()*spectator_4vect_irf.Py()+struck_mass*struck_mass)/(alpha_stru*MASS_DEUTERON);
		//new 4 vector for struck nucleon before interaction;
		TLorentzVector n_primeprime;
		n_primeprime.SetPxPyPzE(-spectator_4vect_irf.Px(),-spectator_4vect_irf.Py(),
		Pz_bInt,E_bInt);
		double t2_uppervtx = (pnew - n_primeprime).Mag2();
		t_nprimeprime->Fill( -t2_uppervtx );

		// the difference might be the primary interaction still
		// doesn't know about the fermi momentum
		t_compare->Fill(-t_hat, -t2_uppervtx);
		
		//true t? Not sure what the true t is in eD.
		t_truth->Fill( -t_hat );

		//spectral function
		if(alpha_spec > 0 && spectator_4vect_irf.Pt() > 0. ) {
			h_dNdAlphadPt2->Fill( alpha_spec, spectator_4vect_irf.Pt(), 1./(2*PI*spectator_4vect_irf.Pt()) );
			if( alpha_spec >= 0.99 && alpha_spec < 1.01 ) {
				h_spectral_pt->Fill(spectator_4vect_irf.Pt()*spectator_4vect_irf.Pt() );
				double MASS_NUCLEON = (MASS_NEUTRON + MASS_PROTON) / 2.;
				double epsilon = 2*MASS_NUCLEON - MASS_DEUTERON;
				double a2 = MASS_NUCLEON*epsilon - epsilon*epsilon/4.;
				double Ra = 4*sqrt(MASS_NUCLEON*MASS_NUCLEON-a2)*TMath::Gamma(2-alpha_spec)*TMath::Gamma(2-alpha_spec);
				double Sd = h_spectral_pt_input->GetBinContent( h_spectral_pt_input->FindBin( -tt ) );
				h_spectralAtPole->Fill( -tt, Sd*(tt*tt)/Ra );
			}
		}

		//angle between photon and spectator in d rest frame
		double angle = spectator_4vect_irf.Angle(q_irf.Vect());
		h_ThetaRprimePm->Fill( angle, spectator_4vect_irf.P() );

		TLorentzVector pn_final;
		if( pnew.E() != 0. && spectator_4vect_irf.E() != 0.){
			pn_final = pnew+spectator_4vect_irf;
			sPN->Fill( pn_final.Mag2() );		
		}

		//struck nucleon 3 momentum:
		Pt_struck->Fill( pnew.Pt() );
		Pz_struck->Fill( pnew.Pz() );
		Pp_struck->Fill( pnew.P() );

		//spectator nucleon 3 momentum:
		Pt_spectator->Fill( spectator_4vect_irf.Pt() );
		Pz_spectator->Fill( spectator_4vect_irf.Pz() );
		Pp_spectator->Fill( spectator_4vect_irf.P() );

		//Jpsi VM 3 momentum in lab frame:
		jnew.Boost(b);
		Pt_VM->Fill( jnew.Pt() );
		Pt2_VM->Fill( jnew.Pt()*jnew.Pt() );
		Pz_VM->Fill( jnew.Pz() );
		Pp_VM->Fill( jnew.P() );
	}

	output->Write();
	output->Close();

}