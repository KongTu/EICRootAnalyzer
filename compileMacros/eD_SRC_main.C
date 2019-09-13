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
		if( p.Theta() > 0.008 ) pass = false;
	}
	else{
		if( (p.Theta() > 0.005 && p.Theta() < 0.007) || p.Theta() > 0.022 ) pass = false;
	}

	return pass;
}

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
	TH1D* nk_truth = new TH1D("nk_truth","k (GeV/c)", nk_nBins, nk_bins);
	TH2D* EvsPz = new TH2D("EvsPz",";pz;E",500,-0.01,0.01,500,-0.01,0.01);
	TH2D* EvsPzFix = new TH2D("EvsPzFix",";pz;E",500,-0.01,0.01,500,-0.01,0.01);
	
	TH1D* Pp_mag[5];
	TH2D* P_spa[5];
	TH1D* nucleon_t[5]; 
	TH1D* sPN[5]; 
	TH1D* sPN_Fpt2[5]; 
	TH2D* sPN_t[5]; 
	TH2D* sPN_k[5]; 
	TH2D* sPN_Fpt2_k[5];
	for(int i=0;i<5;i++){
		Pp_mag[i] = new TH1D(Form("Pp_mag_%d",i),";P (GeV/c)",500,0,5);
		P_spa[i] = new TH2D(Form("P_spa_%d",i),";x(m);y(m)",200,-1,1,200,-1,1);
		nucleon_t[i] = new TH1D(Form("nucleon_t_%d",i),"t (GeV^{2})",200,-10,10);
		sPN[i] = new TH1D(Form("sPN_%d",i),"sPN",sPN_nBins,sPN_bins);
		sPN_Fpt2[i] = new TH1D(Form("sPN_Fpt2_%d",i),"sPN_Fpt2",200,0,10);
		sPN_t[i] = new TH2D(Form("sPN_t_%d",i),";t;s",200,-10,10,sPN_nBins,sPN_bins);
		sPN_k[i] = new TH2D(Form("sPN_k_%d",i),";k;s",200,0,1,sPN_nBins,sPN_bins);
		sPN_Fpt2_k[i] = new TH2D(Form("sPN_Fpt2_k_%d",i),";k;4p^{2}_{T}",200,0,1,200,0,10);
	}
	
	TH1D* Np_mag[2];
	TH2D* N_spa[2];
	TH1D* nk_spectator[2];
	for(int i=0;i<2;i++){
		Np_mag[i] = new TH1D(Form("Np_mag_%d",i),";P (GeV/c)",500,0,5);
		N_spa[i] = new TH2D(Form("N_spa_%d",i),";x(m);y(m)",200,-1,1,200,-1,1);
		nk_spectator[i] = new TH1D(Form("nk_spectator_%d",i),";k (GeV/c)", nk_nBins, nk_bins);
	}


	TChain *tree = new TChain("EICTree");
	tree->Add("/eicdata/eic0003/ztu/BeAGLE_devK/"+filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	double energy_resolution = rZDC;//50%
	TF1* smear_e = new TF1("smear_e","gaus(0)",-30,30);
	smear_e->SetParameter(0,1);
	smear_e->SetParameter(1,0);
	smear_e->SetParameter(2, sqrt( (energy_resolution/sqrt(135.))*(energy_resolution/sqrt(135.)) + 0.02*0.02 )*135. );

	TF1* smear_theta = new TF1("smear_theta","gaus(0)",-0.001,0.001);
	smear_theta->SetParameter(0,1);
	smear_theta->SetParameter(1,0);
	double dis_reso = 0.1/sqrt(135.0);
	double angle_reso = TMath::ATan(dis_reso/28.8);
	smear_theta->SetParameter(2,angle_reso);//assume 28.8m away from IP and 10cm/sqrt(E) resolution

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
		TLorentzVector jnew,pnew;
		TLorentzVector jnew1,pnew1;
		TLorentzVector jnew2,pnew2;
		TLorentzVector jnew3,pnew3,nnew3;

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

			if(pdg == 443 ) j_4vect_irf = j_4vect;
			if(pdg == 2212) p_4vect_irf = p_4vect;
			if(pdg == 2112) n_4vect_irf = n_4vect; 

		
			nParticles_process++;

		} // end of particle loop

		nk_truth->Fill( sqrt(pxf*pxf+pyf*pyf+pzf*pzf) );

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
		pnew1.SetPxPyPzE(px_new,py_new,pz_new, sqrt( struck_mass*struck_mass + px_new*px_new + py_new*py_new + pz_new*pz_new));
		
		jx_new = j_4vect_irf.Px()+spectator_4vect_irf.Px();
		jy_new = j_4vect_irf.Py()+spectator_4vect_irf.Py();
		jz_new = lfjz;
		jnew1.SetPxPyPzE(jx_new,jy_new,jz_new, sqrt( MASS_JPSI*MASS_JPSI + jx_new*jx_new + jy_new*jy_new + jz_new*jz_new));

		TLorentzVector testnew2 = q_irf+d_beam_irf-jnew1-pnew1-spectator_4vect_irf;
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


		// approach 4 with touching the spectator
		// final state struck nucleon + k = 3 vectors for the "kick"
		// evenly distribute it to both nucleons. 
		// should consider the off shell mass for spectator as well, the kick is on the 3 vectors when still off shell
		
		TVector3 kick(pnew2.Px()-(-spectator_4vect_irf.Px()), pnew2.Py()-(-spectator_4vect_irf.Py()), pnew2.Pz()-(-spectator_4vect_irf.Pz()) );
		double kick_x = kick.Px();
		double kick_y = kick.Py();
		double kick_z = kick.Pz();
		
		kmag = spectator_4vect_irf.P();
		MnuclOff = sqrt(0.25*MASS_DEUTERON*MASS_DEUTERON - kmag*kmag);
		PpOff = sqrt( spectator_mass*spectator_mass - MnuclOff*MnuclOff + spectator_4vect_irf.P()*spectator_4vect_irf.P() );
		Ptheta = spectator_4vect_irf.Theta();
		Pnewpt = PpOff*TMath::Sin(Ptheta);
		
		TLorentzVector Noff4vector; Noff4vector.SetPtEtaPhiM(Pnewpt,spectator_4vect_irf.Eta(),spectator_4vect_irf.Phi(),MnuclOff);
		double spectator_Px = Noff4vector.Px() + kick_x/2.0;
		double spectator_Py = Noff4vector.Py() + kick_y/2.0;
		double spectator_Pz = Noff4vector.Pz() + kick_z/2.0;
		double spectator_E = sqrt(spectator_Px*spectator_Px+spectator_Py*spectator_Py+spectator_Pz*spectator_Pz+spectator_mass*spectator_mass);
		//now putting spectator back on shell with a kick.
		nnew3.SetPxPyPzE(spectator_Px,spectator_Py,spectator_Pz,spectator_E);
		
		qzkz = q_irf.Pz() - nnew3.Pz();
		numn = q_irf.E() - nnew3.E();
		jx = j_4vect_irf.Px()+spectator_4vect_irf.Px()-(Poff4vector.Px()-struck_4vect_irf.Px())-(Noff4vector.Px()-spectator_4vect_irf.Px());
		jy = j_4vect_irf.Py()+spectator_4vect_irf.Py()-(Poff4vector.Py()-struck_4vect_irf.Py())-(Noff4vector.Py()-spectator_4vect_irf.Py());
		px = Pon4vectorNew.Px()-(kick_x/2.0);
		py = Pon4vectorNew.Py()-(kick_y/2.0);

		jz = getCorrJz(qzkz,numn,jx,jy,px,py,struck_mass);
		pz = getCorrPz(qzkz,numn,jx,jy,px,py,struck_mass);

		px_new = px;
		py_new = py;
		pz_new = pz;
		pnew3.SetPxPyPzE(px_new,py_new,pz_new, sqrt( struck_mass*struck_mass + px_new*px_new + py_new*py_new + pz_new*pz_new));
		
		jx_new = jx;
		jy_new = jy;
		jz_new = jz;
		jnew3.SetPxPyPzE(jx_new,jy_new,jz_new, sqrt( MASS_JPSI*MASS_JPSI + jx_new*jx_new + jy_new*jy_new + jz_new*jz_new));


		//filling histograms:
		if( doAcceptance_ ) {
			if( !passDetector(struck_4vect_irf,b) ) struck_4vect_irf.SetPxPyPzE(0,0,0,0);
			if( !passDetector(pnew,b) ) pnew.SetPxPyPzE(0,0,0,0);
			if( !passDetector(pnew1,b) ) pnew1.SetPxPyPzE(0,0,0,0);
			if( !passDetector(pnew2,b) ) pnew2.SetPxPyPzE(0,0,0,0);
			if( !passDetector(pnew3,b) ) pnew3.SetPxPyPzE(0,0,0,0);
			if( !passDetector(spectator_4vect_irf,b) ) spectator_4vect_irf.SetPxPyPzE(0,0,0,0);
			if( !passDetector(nnew3,b) ) nnew3.SetPxPyPzE(0,0,0,0);
		}
		if( doSmear_ ){
			struck_4vect_irf = afterDetector(struck_4vect_irf,b,smear_e,smear_theta);
			pnew = afterDetector(pnew,b,smear_e,smear_theta);
			pnew1 = afterDetector(pnew1,b,smear_e,smear_theta);
			pnew2 = afterDetector(pnew2,b,smear_e,smear_theta);
			pnew3 = afterDetector(pnew3,b,smear_e,smear_theta);
			spectator_4vect_irf = afterDetector(spectator_4vect_irf,b,smear_e,smear_theta);
			nnew3 = afterDetector(nnew3,b,smear_e,smear_theta);
		}

		//filling histograms:
		Pp_mag[0]->Fill( struck_4vect_irf.P() );
		Pp_mag[1]->Fill( pnew.P() );
		Pp_mag[2]->Fill( pnew1.P() );
		Pp_mag[3]->Fill( pnew2.P() );
		Pp_mag[4]->Fill( pnew3.P() );

		Np_mag[0]->Fill( spectator_4vect_irf.P() );
		Np_mag[1]->Fill( nnew3.P() );


		TLorentzVector pn_final;double Fpt2=0.0;
		if( struck_4vect_irf.E() != 0. && spectator_4vect_irf.E() != 0.){
			pn_final = struck_4vect_irf+spectator_4vect_irf;
			struck_4vect_irf.Boost(b);spectator_4vect_irf.Boost(b);
			Fpt2 = (struck_4vect_irf.Pt()+spectator_4vect_irf.Pt())*(struck_4vect_irf.Pt()+spectator_4vect_irf.Pt());
			struck_4vect_irf.Boost(-b);spectator_4vect_irf.Boost(-b);

			//default BeAGLE:
			nucleon_t[0]->Fill( (pn_final - d_beam_irf).Mag2() );
			sPN[0]->Fill( pn_final.Mag2() );
			sPN_Fpt2[0]->Fill( Fpt2 );//4*spectator pt**2
			sPN_t[0]->Fill((pn_final - d_beam_irf).Mag2(), pn_final.Mag2() );
			sPN_k[0]->Fill(nk_event, pn_final.Mag2());
			sPN_Fpt2_k[0]->Fill(nk_event, Fpt2 );
		}
		if( pnew.E() != 0. && spectator_4vect_irf.E() != 0.){
			//approach 1:
			pn_final = pnew+spectator_4vect_irf;
			pnew.Boost(b);spectator_4vect_irf.Boost(b);
			Fpt2 = (pnew.Pt()+spectator_4vect_irf.Pt())*(pnew.Pt()+spectator_4vect_irf.Pt());
			pnew.Boost(-b);spectator_4vect_irf.Boost(-b);
			
			nucleon_t[1]->Fill( (pn_final - d_beam_irf).Mag2() );
			sPN[1]->Fill( pn_final.Mag2() );
			sPN_Fpt2[1]->Fill( Fpt2 );//4*spectator pt**2
			sPN_t[1]->Fill((pn_final - d_beam_irf).Mag2(), pn_final.Mag2() );
			sPN_k[1]->Fill(nk_event, pn_final.Mag2());
			sPN_Fpt2_k[1]->Fill(nk_event, Fpt2 );
		}
		if( pnew1.E() != 0. && spectator_4vect_irf.E() != 0.){
			//approach 2:
			pn_final = pnew1+spectator_4vect_irf;
			pnew1.Boost(b);spectator_4vect_irf.Boost(b);
			Fpt2 = (pnew1.Pt()+spectator_4vect_irf.Pt())*(pnew1.Pt()+spectator_4vect_irf.Pt());
			pnew1.Boost(-b);spectator_4vect_irf.Boost(-b);		

			nucleon_t[2]->Fill( (pn_final - d_beam_irf).Mag2() );
			sPN[2]->Fill( pn_final.Mag2() );
			sPN_Fpt2[2]->Fill( Fpt2 );//4*spectator pt**2
			sPN_t[2]->Fill((pn_final - d_beam_irf).Mag2(), pn_final.Mag2() );
			sPN_k[2]->Fill(nk_event, pn_final.Mag2());
			sPN_Fpt2_k[2]->Fill(nk_event, Fpt2 );
		}
		if( pnew2.E() != 0. && spectator_4vect_irf.E() != 0.){
			//approach 3:
			pn_final = pnew2+spectator_4vect_irf;
			pnew2.Boost(b);spectator_4vect_irf.Boost(b);
			Fpt2 = (pnew2.Pt()+spectator_4vect_irf.Pt())*(pnew2.Pt()+spectator_4vect_irf.Pt());
			pnew2.Boost(-b);spectator_4vect_irf.Boost(-b);		

			nucleon_t[3]->Fill( (pn_final - d_beam_irf).Mag2() );
			sPN[3]->Fill( pn_final.Mag2() );
			sPN_Fpt2[3]->Fill( Fpt2 );//4*spectator pt**2
			sPN_t[3]->Fill((pn_final - d_beam_irf).Mag2(), pn_final.Mag2() );
			sPN_k[3]->Fill(nk_event, pn_final.Mag2());
			sPN_Fpt2_k[3]->Fill(nk_event, Fpt2 );
		}
		if( pnew3.E() != 0. && nnew3.E() != 0.){
			//approach 3:
			pn_final = pnew3+nnew3;
			pnew3.Boost(b);nnew3.Boost(b);
			Fpt2 = (pnew3.Pt()+nnew3.Pt())*(pnew3.Pt()+nnew3.Pt());
			pnew3.Boost(-b);nnew3.Boost(-b);		

			nucleon_t[4]->Fill( (pn_final - d_beam_irf).Mag2() );
			sPN[4]->Fill( pn_final.Mag2() );
			sPN_Fpt2[4]->Fill( Fpt2 );//4*spectator pt**2
			sPN_t[4]->Fill((pn_final - d_beam_irf).Mag2(), pn_final.Mag2() );
			sPN_k[4]->Fill(nk_event, pn_final.Mag2());
			sPN_Fpt2_k[4]->Fill(nk_event, Fpt2 );
		}

		//use spectator only:
		if( spectator_4vect_irf.E() != 0. ) nk_spectator[0]->Fill( spectator_4vect_irf.P() );
		if( nnew3.E() != 0. ) nk_spectator[1]->Fill( nnew3.P() );

		//spatial distributions, first boost back in lab frame:
		
		if(struck_4vect_irf.E() != 0.) struck_4vect_irf.Boost(b);
		if(pnew.E() != 0.) pnew.Boost(b);
		if(pnew1.E() != 0.) pnew1.Boost(b);
		if(pnew2.E() != 0.) pnew2.Boost(b);
		if(pnew3.E() != 0.) pnew3.Boost(b);
		if(spectator_4vect_irf.E() != 0.) spectator_4vect_irf.Boost(b);
		if(nnew3.E() != 0.) nnew3.Boost(b);

		if( !struckproton ){
			vector< double> pos = getPspa(struck_4vect_irf);
			P_spa[0]->Fill(pos[0],pos[1]);
			pos.clear(); pos = getPspa(pnew);
			P_spa[1]->Fill(pos[0],pos[1]);
			pos.clear(); pos = getPspa(pnew1);
			P_spa[2]->Fill(pos[0],pos[1]);
			pos.clear(); pos = getPspa(pnew2);
			P_spa[3]->Fill(pos[0],pos[1]);
			pos.clear(); pos = getPspa(pnew3);
			P_spa[4]->Fill(pos[0],pos[1]);
		}
		else{
			vector< double> pos;
			pos.clear(); pos = getPspa(spectator_4vect_irf);
			N_spa[0]->Fill(pos[0],pos[1]);
			pos.clear(); pos = getPspa(nnew3);
			N_spa[1]->Fill(pos[0],pos[1]);
		}

	}

	output->Write();
	output->Close();




}