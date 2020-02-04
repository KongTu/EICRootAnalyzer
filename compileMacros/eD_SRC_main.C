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

double acceptanceGlobal = 0.004;

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

vector<double> getPspa(TLorentzVector p){

	vector< double> temp;
	if( p.E() == 0. ){
		temp.push_back(-999.);
		temp.push_back(-999.);
	}
	else{
		double zdcip = 28.8;
		double dp_struck = zdcip*TMath::Tan(p.Theta());
		double P_sx = dp_struck*TMath::Cos(p.Phi());
		double P_sy = dp_struck*TMath::Sin(p.Phi());
		temp.push_back(P_sx);
		temp.push_back(P_sy);
	}

	return temp;
}

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

//only for ZDC for now
TLorentzVector afterDetector(TLorentzVector p, TVector3 b, TF1*smear_e, TF1*smear_theta){

	bool isNeutron = true;
	TLorentzVector pafter;

	if( p.E() == 0. ) {
		return p;
	}
	//boost to lab frame;
	p.Boost(b);

	if( p.M() < MASS_PROTON+0.0001 ) isNeutron = false;
	if( !isNeutron ) {
		pafter = p;
	}
	else{
		//smearing neutron
		double E_n = p.E();
		double delta_E = smear_e->GetRandom();
		E_n = E_n + delta_E;
		double delta_Theta = smear_theta->GetRandom();
		double angle = p.Theta() + delta_Theta;
		double Pz_n2 = (E_n*E_n - MASS_NEUTRON*MASS_NEUTRON)/(1+TMath::Sin(angle)*TMath::Sin(angle));
		double Pz_n = sqrt(Pz_n2);
		double Pt_n2 = (E_n*E_n - MASS_NEUTRON*MASS_NEUTRON - Pz_n2);
		double Pt_n = sqrt(Pt_n2);
		double Px_n = Pt_n*TMath::Cos(p.Phi());
		double Py_n = Pt_n*TMath::Sin(p.Phi());

		pafter.SetPxPyPzE(Px_n, Py_n, Pz_n, E_n);
	}

	//boost back to IRF;
	pafter.Boost(-b);
	return pafter;
}

void eD_SRC_main(const int nEvents = 40000, TString filename="", const bool doSmear_ = false, const bool doAcceptance_ = false, const double rZDC = 1., const double acceptance=0.004){

	acceptanceGlobal = acceptance;
	std::ostringstream os;
	os << (int) doSmear_;
	os << (int) doAcceptance_;
	os << "_ZDC_" << (double) rZDC;
	os << "_" << (double) acceptance;
	std::string str = os.str();
	TString settings = (TString) str;

	TFile * output = new TFile("../rootfiles/"+filename+"_"+settings+"_main_Beagle.root","recreate");
		
	TH1D* nk_truth = new TH1D("nk_truth","k (GeV/c)", nk_nBins, nk_bins);
	TH1D* nk_truth_uniformbins = new TH1D("nk_truth_uniformbins","k (GeV/c)", 200,0,1.4);
	TH2D* h_ThetaVsEnergy_Spectator = new TH2D("h_ThetaVsEnergy_Spectator",";E_{spectator} (GeV);#theta",300,0,200,200,0,100);
	TH2D* h_ThetaVsEnergy_Struck = new TH2D("h_ThetaVsEnergy_Struck",";E_{spectator} (GeV);#theta",300,0,200,200,0,100);
	TH1D* sPN = new TH1D(Form("sPN"),"sPN",sPN_nBins,sPN_bins);
	TH1D* Pt_struck = new TH1D("Pt_struck",";p_{T} (GeV)",200,0,1);
	TH1D* Pz_struck = new TH1D("Pz_struck",";p_{z} (GeV)",200,-1,1);
	TH1D* Pp_struck = new TH1D("Pp_struck",";p (GeV)",200,0,1.4);
	TH1D* Pt_spectator = new TH1D("Pt_spectator",";p_{T} (GeV)",200,0,1);
	TH1D* Pz_spectator = new TH1D("Pz_spectator",";p_{z} (GeV)",200,-1,1);
	TH1D* Pp_spectator = new TH1D("Pp_spectator",";p (GeV)",200,0,1.4);
	TH2D* spa_struck = new TH2D("spa_struck",";x(m);y(m)",200,-1,1,200,-1,1);
	TH2D* spa_spectator = new TH2D("spa_spectator",";x(m);y(m)",200,-1,1,200,-1,1);
	TH1D* alpha_spectator = new TH1D("alpha_spectator",";#alpha_{spec}",100,0,2);
	TH1D* ttprime = new TH1D("ttprime",";-t'(GeV)",100,0,2);
	TH2D* h_ttprime_alpha = new TH2D("h_ttprime_alpha",";#alpha_{p};-t'",200,0,2,1000,0,0.1);
	TH2D* h_dNdAlphadPt2 = new TH2D("h_dNdAlphadPt2",";#alpha_{p};p_{T} (GeV/c)'",500,0,2,500,0,2);

	TChain *tree = new TChain("EICTree");
	tree->Add("/gpfs02/eic/ztu/BeAGLE/BeAGLE_devK_SRC/"+filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	double energy_resolution = rZDC;//50%
	double energy_resolution_constant_term = 0.05; //5%
	double beam_momentum = 110.; // 110 GeV for Deuteron now
	TF1* smear_e = new TF1("smear_e","gaus(0)",-30,30);
	smear_e->SetParameter(0,1);
	smear_e->SetParameter(1,0);
	smear_e->SetParameter(2, sqrt( TMath::Power((energy_resolution/sqrt(beam_momentum)),2) 
		+ TMath::Power(energy_resolution_constant_term,2))*beam_momentum );
	//resolution adding in quadrature. 

	TF1* smear_theta = new TF1("smear_theta","gaus(0)",-0.001,0.001);
	smear_theta->SetParameter(0,1);
	smear_theta->SetParameter(1,0);
	double angle_reso = 3e-6;
	smear_theta->SetParameter(2,angle_reso);
	//1cm position resolution --> 3 microRad resolution
	//Yuji's Letter of Intent in EIC R&D proposal

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
		bool struckproton = false;
		if( struck_nucleon == 2212 ) struckproton = true;

		int nParticles_process = 0;
		TLorentzVector n_4vect_unsmear;
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
			q = e_beam-e_scattered;
			q_irf = q;
			
			if(pdg == 443 ) j_4vect = ppart;//jpsi
			if(pdg == 2212) p_4vect = ppart;//proton
			if(pdg == 2112) {n_4vect = ppart; n_4vect_unsmear = n_4vect;}//neutron

			//prepare for boost later
			if(pdg == 443 ) j_4vect_irf = j_4vect;
			if(pdg == 2212) p_4vect_irf = p_4vect;
			if(pdg == 2112) n_4vect_irf = n_4vect; 

			nParticles_process++;

		} // end of particle loop

		//fill n(k) or dN/dk distribution, but averaged over all direction
		//LFKine tells us pzf is not symmetric in the lab frame
		nk_truth->Fill( nk_event );
		nk_truth_uniformbins->Fill( nk_event );

		//just a protection
		if( p_4vect.E() == 0 || n_4vect.E() == 0 ) continue;

		//boost
		j_4vect_irf.Boost(-b);
		p_4vect_irf.Boost(-b);
		n_4vect_irf.Boost(-b);
		q_irf.Boost(-b);
		d_beam_irf.Boost(-b);

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
		
		

		//filling histograms:
		if( doAcceptance_ ) {
			if( !passDetector(pnew,b) ) pnew.SetPxPyPzE(0,0,0,0);
			if( !passDetector(spectator_4vect_irf,b) ) spectator_4vect_irf.SetPxPyPzE(0,0,0,0);
		}
		if( doSmear_ ){
			pnew = afterDetector(pnew,b,smear_e,smear_theta); //dummy for now, only smear ZDC neutron
			spectator_4vect_irf = afterDetector(spectator_4vect_irf,b,smear_e,smear_theta);
		}

		//fill lab frame theta vs Energy for struck and spectator
		h_ThetaVsEnergy_Spectator->Fill(spectator_4vect.Pz(), spectator_4vect.Theta()*1000. );
		TLorentzVector pnew_lab;
		pnew.Boost(b);
		pnew_lab = pnew;
		pnew.Boost(-b);
		h_ThetaVsEnergy_Struck->Fill(pnew_lab.Pz(), pnew_lab.Theta()*1000. );

		//filling alpha of spectator
		double Pplus = (spectator_4vect_irf.E() + spectator_4vect_irf.Pz()) / sqrt(2);
		double PdPlus = MASS_DEUTERON / sqrt(2);
		double alpha_spec = 2*Pplus / PdPlus;
		alpha_spectator->Fill( alpha_spec );

		//filling t' distribution
		double tt = (spectator_4vect_irf - d_beam_irf).Mag2();
		tt = tt - TMath::Power(spectator_4vect_irf.M(),2);
		ttprime->Fill( -tt );
		h_ttprime_alpha->Fill( alpha_spec, -tt );

		//spectral function
		h_dNdAlphadPt2->Fill( alpha_spec, spectator_4vect_irf.Pt(), 1./(2*PI*spectator_4vect_irf.Pt()) );


		TLorentzVector pn_final;
		if( pnew.E() != 0. && spectator_4vect_irf.E() != 0.){
			pn_final = pnew+spectator_4vect_irf;
			sPN->Fill( pn_final.Mag2() );		
		}
		//filling histograms:
		Pt_struck->Fill( pnew.Pt() );
		Pz_struck->Fill( pnew.Pz() );
		Pp_struck->Fill( pnew.P() );

		//use spectator only:
		Pt_spectator->Fill( spectator_4vect_irf.Pt() );
		Pz_spectator->Fill( spectator_4vect_irf.Pz() );
		Pp_spectator->Fill( spectator_4vect_irf.P() );

		//spatial distributions
		vector< double> pos;
		if( !struckproton ){
			//in case has smearing
			spectator_4vect_irf.Boost(b);
			pnew.Boost(b);
			pos.clear(); pos = getPspa(spectator_4vect_irf);
			spa_struck->Fill(pos[0],pos[1]);
			pos.clear(); pos = getPspa(pnew);
			spa_spectator->Fill(pos[0],pos[1]);
		}
		else{
			spectator_4vect_irf.Boost(b);
			pnew.Boost(b);
			pos.clear(); pos = getPspa(pnew);
			spa_struck->Fill(pos[0],pos[1]);
			pos.clear(); pos = getPspa(spectator_4vect_irf);
			spa_spectator->Fill(pos[0],pos[1]);
		}

	}

	for(int ibin=0;ibin<Pp_spectator->GetNbinsX();ibin++){
		Pp_spectator->SetBinContent(ibin+1, Pp_spectator->GetBinContent(ibin+1)*(1./TMath::Power(Pp_spectator->GetBinCenter(ibin+1),2)));
	}

	output->Write();
	output->Close();

}