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

#define PI            3.1415926

#define MASS_JPSI 	  3.09688
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_DEUTERON 1.8756129
// #define MASS_DEUTERON 1.8751019071673038
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

#define MASS_MUON  0.1056

int findSpectator(TVector3 p, int charge=-99){
	
	int candidate=0;
	bool couldBeNeutron = true;
	if(charge!=0) couldBeNeutron = false;
	if(couldBeNeutron&&(p.Theta()*1e3<4.5))candidate=1;
	if(!couldBeNeutron&&(p.Theta()*1e3<20)) candidate=2;

	return candidate;
}
int isMatch(TLorentzVector trueSpect, TLorentzVector taggedSpect){
	if(TMath::Abs(trueSpect.Pt()-taggedSpect.Pt())<3e-3 
		&& TMath::Abs(trueSpect.Eta()-taggedSpect.Eta())<3e-1 ){
			// && TMath::Abs(trueSpect.M()-taggedSpect.M())<1e-4 ){
		return 1;
	}
	else{
		return 0;
	}
}
TLorentzRotation RotateToLab(TLorentzVector const &eBeam_lab,
                            TLorentzVector const &pBeam_lab,
                            TLorentzVector const &eScat_lab) {
	
	TLorentzVector q_lab=eBeam_lab - eScat_lab;
	TLorentzVector q_irf=q_lab;
	TLorentzVector eScat_irf=eScat_lab;
	TVector3 pBoost=pBeam_lab.BoostVector();
	q_irf.Boost(-pBoost);
	eScat_irf.Boost(-pBoost);

	TLorentzRotation l;
	double angleTheta = q_irf.Theta();
	double anglePhi = eScat_irf.Phi();
	l.RotateY( angleTheta );
	l.RotateZ( PI-anglePhi );

	return l;

}

void eD_Tagged_DIS_background(const int nEvents = 40000, double HFSaccept=4.0, bool cutPtBal_=false, TString filename="Output_input_temp_91"){


	//input from BeAGLE root files
	TFile * output = new TFile("../rootfiles/eD_Tagged_DIS_Beagle_background.root","recreate");
	
	TChain *tree = new TChain("EICTree");
	// tree->Add("/gpfs02/eic/ztu/Analysis/BeAGLE/eD_Tagged_DIS/18x110_Q2_10_100/eD_Tagged_DIS_100M_batch_1/"+filename+".root" );
	tree->Add("/gpfs02/eic/ztu/Analysis/BeAGLE/eD_Tagged_DIS/18x110_Q2_10_100_noINC/eD_Tagged_DIS_1M_batch_1/*.root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);

	//all constants
	double totalXSection   = .0000450463252; //mb
	double nEventsTotal        = 500257.0;
	double Lint = nEventsTotal/totalXSection; // mb^{-1}
	double alpha2 = TMath::Power((1./137),2);
	double twopi = 2*PI;
	double mbToGeV_m2 = 2.56819;
	double Q2binwidth = 13.0-10.0;

	TH2D* h_taggingEfficiency_pt2 = new TH2D("h_taggingEfficiency_pt2",";p^{2}_{T,tagged}(GeV^{2});p^{2}_{T,truth}(GeV^{2})", 100, 0, 0.15, 100, 0, 0.15);
	
	TH1D* h_taggingEfficiency = new TH1D("h_taggingEfficiency","",3,-1,2);
	TH1D* h_taggingEfficiency_step2 = new TH1D("h_taggingEfficiency_step2","",3,-1,2);
	
	TH2D* h_ptBalance = new TH2D("h_ptBalance",";pt_{hfsQ};pt_{spec}", 100, 0, 2, 100, 0, 2);
	TH1D* h_ptBalance1D = new TH1D("h_ptBalance1D",";#Delta pt", 100,-1,1);
	
	TH1D* h_beforeTagging = new TH1D("h_beforeTagging","; pt^{2}", 100,0,0.2);
	TH1D* h_afterTagging = new TH1D("h_afterTagging","; pt^{2}", 100,0,0.2);
	TH1D* h_allTagging = new TH1D("h_allTagging","; pt^{2}", 100,0,0.2);
	TH1D* h_wrongTagging = new TH1D("h_wrongTagging","; pt^{2}", 100,0,0.2);

	TH1D* h_HERA_Q2_10_13 = new TH1D("h_HERA_Q2_10_13","h_HERA_Q2_10_13",100,0.00001,0.1);
	TH1D* h_alpha_spec = new TH1D("h_alpha_spec","h_alpha_spec",100,0,2);
	TH1D* h_nk = new TH1D("h_nk","h_nk",100,0,2);
	double bin_width = h_HERA_Q2_10_13->GetBinWidth(1);
	
	double alpha_binning[161];
	for(int ibin=0;ibin<161;ibin++){
		alpha_binning[ibin] = 0.4+ibin*0.01;
	}	
	TH1D* h_HERA_Q2_10_13_x007_009_alpha[160];
	TH1D* h_alpha_spec_everybin[160];
	for(int ibin=0;ibin<160;ibin++){
	 	h_HERA_Q2_10_13_x007_009_alpha[ibin] = new TH1D(Form("h_HERA_Q2_10_13_x007_009_alpha_%d",ibin),Form("h_HERA_Q2_10_13_x007_009_alpha_%d",ibin),100,0,0.15);
		h_alpha_spec_everybin[ibin] = new TH1D(Form("h_alpha_spec_everybin_%d",ibin),Form("h_alpha_spec_everybin_%d",ibin),100,0,2);
	}


	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);

		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		double pznucl = event->pznucl;
		double Atarg = event->Atarg;
		double pztarg_total = pztarg*Atarg;
		double pznucl_total = pznucl*Atarg;

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
		double photon_flux = event->GetPhotonFlux();
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();
		int struck_nucleon = event->nucleon;
		double nk_event = sqrt(pxf*pxf+pyf*pyf+pzf*pzf);
		double Espec = 0.;
		TLorentzVector trueSpect;
		if( struck_nucleon == 2212 ){
			Espec = sqrt(nk_event*nk_event+MASS_NEUTRON*MASS_NEUTRON);
		}
		else{
			Espec = sqrt(nk_event*nk_event+MASS_PROTON*MASS_PROTON);
		}
		trueSpect.SetPxPyPzE(-pxf,-pyf,-pzf,Espec);
		
		h_nk->Fill( nk_event );//sanity check for my wavefunction;
		
		//event process and kinematic phase space
		if( event_process != 99 ) continue;
		if( trueQ2 < 10.  || trueQ2 > 13. ) continue;
		if( trueY > 0.95  || trueY < 0.01 ) continue;
				
		//HERA inclusive cross section
		double event_weight = 1.;
		double Yc = 1. + TMath::Power((1-trueY),2);
		double Emax=-1.;
		double etaMax=-1;
		int bestCandidate=-1;
		TVector3 bestCandidateVector(-1,-1,-1);
		TLorentzVector hfsCand(0.,0.,0.,0.);
		vector< TLorentzVector> saveListOfNucleons;
		for(int j(0); j < nParticles; ++j ) {
			const erhic::ParticleMC* particle = event->GetTrack(j);
			int pdg = particle->GetPdgCode();
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
			double pt = particle->GetPt();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			int status = particle->GetStatus();
			double mass = particle->GetM();
			int charge = particle->eA->charge;
			int NoBAM = particle->eA->NoBam;
			TLorentzVector ppart = particle->Get4Vector();
			if( index == 3 ) {
				e_scattered.SetPtEtaPhiM(pt,eta,phi,0.00051);
			}
			if( status!=1 ) continue;
			if( TMath::Abs(pdg) == 2112 || TMath::Abs(pdg) == 2212 ) saveListOfNucleons.push_back( ppart );
			// saveListOfNucleons.push_back( ppart );
			TLorentzVector part4pion; part4pion.SetPtEtaPhiM(pt,eta,phi,0.13957);//assume pions
		    //sum over HFS excluding elec' within main detector acceptance;
		    if(!(isMatch(ppart,e_scattered)) && TMath::Abs(part4pion.Eta())<HFSaccept ) hfsCand += part4pion;
		    // if(!(isMatch(ppart,e_scattered)) && !(isMatch(part4pion,trueSpect)) ) hfsCand += part4pion;
			TVector3 part; part.SetPtEtaPhi(pt, eta, phi);
			if(TMath::Abs(pdg)==211) part.SetTheta(part.Theta()*6.725);//rigidity change
			if(TMath::Abs(pdg)==321) part.SetTheta(part.Theta()*1.901);//rigidity change
			int spec_cand = findSpectator(part, charge);
			if( spec_cand ){
				if(part.Eta()>etaMax && part.Mag()>Emax ) {
					etaMax=part.Eta();
					Emax=part.Mag();
					bestCandidate=spec_cand;
					bestCandidateVector=part;
				}
			}
		}
		TLorentzRotation rotateVector=RotateToLab(e_beam, d_beam, e_scattered);
		TLorentzVector trueSpect_lab = rotateVector*trueSpect;//rotation only
		trueSpect_lab.Boost(b);//longitudinal boost without rotation
		trueSpect.Boost(b);

		// TLorentzVector qbeam_IRF=(e_beam - e_scattered);
		// TLorentzVector e_scattered_IRF = e_scattered;
		// e_scattered.Boost(-b);
		// qbeam_IRF.Boost(-b);
		// double angleTheta = qbeam_IRF.Theta();
		// double anglePhi = e_scattered.Phi();
		// rotateVector.RotateY( angleTheta );
		// rotateVector.RotateZ( PI-anglePhi );
		// TLorentzVector trueSpect_lab = rotateVector*trueSpect;
		// trueSpect_lab.Boost(b);
		// trueSpect.Boost(b);
		
		cout << "before rotaton pt " << trueSpect.Pt() << " mass " << trueSpect.M() << " eta " << trueSpect.Eta() << " phi " << trueSpect.Phi() << " total p " << trueSpect.P() << endl; 
		cout << "after rotaton pt " << trueSpect_lab.Pt() << " mass " << trueSpect_lab.M() << " eta " << trueSpect_lab.Eta() << " phi " << trueSpect_lab.Phi() << " total p " << trueSpect_lab.P() << endl; 
		for(unsigned icand=0; icand<saveListOfNucleons.size(); icand++){
			cout << "candidate " << icand << " mass " << saveListOfNucleons[icand].M() 
			<< " pt " << saveListOfNucleons[icand].Pt() << " eta " << saveListOfNucleons[icand].Eta()  << " total p " << saveListOfNucleons[icand].P() <<endl;
		}
		//don't touch below
		
		h_beforeTagging->Fill( TMath::Power(trueSpect.Pt(),2) );
		//virtual photon
		TLorentzVector qbeam = e_beam - e_scattered;
		//initialize spectator 4vect
		TLorentzVector spectator_4vect_irf;
		if(bestCandidate<0) continue;
		if(bestCandidate==1) {
			spectator_4vect_irf.SetPtEtaPhiM(bestCandidateVector.Pt(), bestCandidateVector.Eta(), bestCandidateVector.Phi(), MASS_NEUTRON);
		}if(bestCandidate==2){
			spectator_4vect_irf.SetPtEtaPhiM(bestCandidateVector.Pt(), bestCandidateVector.Eta(), bestCandidateVector.Phi(), MASS_PROTON);
		}
		h_taggingEfficiency->Fill(isMatch(trueSpect, spectator_4vect_irf));
		//pt balance 2D and 1D
		h_ptBalance->Fill( (qbeam-hfsCand).Pt(), spectator_4vect_irf.Pt() );
		h_ptBalance1D->Fill( (qbeam-hfsCand).Pt() - spectator_4vect_irf.Pt() );
		// cut pt Bal
		if( cutPtBal_ ) {
			if( ((qbeam-hfsCand).Pt()-spectator_4vect_irf.Pt())>0.01 ) continue;
		}
		//algo step 1 eff.
		h_allTagging->Fill( TMath::Power(spectator_4vect_irf.Pt(),2) );
		if( !isMatch(trueSpect, spectator_4vect_irf) ) h_wrongTagging->Fill( TMath::Power(spectator_4vect_irf.Pt(),2) );
		if( isMatch(trueSpect, spectator_4vect_irf) ) h_afterTagging->Fill( TMath::Power(trueSpect.Pt(),2) );
		//if turn on cut on pt balance variable.
		h_taggingEfficiency_step2->Fill(isMatch(trueSpect, spectator_4vect_irf));
		h_taggingEfficiency_pt2->Fill( TMath::Power(spectator_4vect_irf.Pt(),2), TMath::Power(trueSpect.Pt(),2) );
		// if( !isMatch(trueSpect, spectator_4vect_irf) ){
		// 	cout << "start~" << i << " struck " << struck_nucleon << endl;
		// 	cout << "true spectator pt " << trueSpect.Pt() << " eta " << trueSpect.Eta() << " mass " << trueSpect.M() << " total p " << trueSpect.P() << endl;
		// 	cout << "tagged spectator pt " << spectator_4vect_irf.Pt() << " eta " << spectator_4vect_irf.Eta() << " mass " << spectator_4vect_irf.M() << " total p " << spectator_4vect_irf.P() << endl;
		// 	cout << "is matched " << isMatch(trueSpect, spectator_4vect_irf) << endl;
		// 	for(unsigned icand=0; icand<saveListOfNucleons.size(); icand++){
		// 		cout << "candidate " << icand << " mass " << saveListOfNucleons[icand].M() 
		// 		<< " pt " << saveListOfNucleons[icand].Pt() << " eta " << saveListOfNucleons[icand].Eta()  << " total p " << saveListOfNucleons[icand].P() <<endl;
		// 	}
		// }
		saveListOfNucleons.clear();
		//boost back to IRF, continue analysis on cross sections
		spectator_4vect_irf.Boost(-b);
		double xd = trueQ2 / (2*d_beam.Dot(qbeam));
		double gamma2 = (4.*TMath::Power(MASS_DEUTERON,2)*TMath::Power(xd,2)) / trueQ2;
		double epsilon = (1. - trueY - gamma2*TMath::Power(trueY/2.,2)) / (1. - trueY + TMath::Power(trueY,2)/2. + gamma2*TMath::Power(trueY/2.,2) );
		double compare = TMath::Power( trueY, 2) / (1. - epsilon);
	
		event_weight = (TMath::Power(trueQ2,2)*trueX) / (twopi*alpha2*Yc);
		event_weight = event_weight * (mbToGeV_m2)/(Lint*bin_width*Q2binwidth);
		//fill HERA inclusive cross section for Q2(10,13) GeV**2:
		h_HERA_Q2_10_13->Fill( trueX, event_weight );
		//x bin [0.007,0.009]
		double Pplus = (spectator_4vect_irf.E() + spectator_4vect_irf.Pz()) / sqrt(2);
		double PdPlus = MASS_DEUTERON / sqrt(2);
		double alpha_spec = 2*Pplus / PdPlus;
		double pt2 = TMath::Power(spectator_4vect_irf.Pt(),2);
		double alpha_spec_binwidth = -1; // will have to be rewritten by 20 alpha bins
		double xbinwidth = (0.003-0.002);
		double pt2binwidth = h_HERA_Q2_10_13_x007_009_alpha[0]->GetBinWidth(1);
		h_alpha_spec->Fill( alpha_spec );
		if( trueX > 0.003 || trueX < 0.002 ) continue;
		
		int alpha_bin_index = 0;
		for(int ibin=0;ibin<160;ibin++){
			if( alpha_spec>alpha_binning[ibin] && alpha_spec<alpha_binning[ibin+1] ){
				alpha_bin_index = ibin;
				alpha_spec_binwidth = alpha_binning[ibin+1] - alpha_binning[ibin];
			}
		}
		double event_weight_alphaPt2 = alpha_spec*(16.*TMath::Power(PI,1)*(TMath::Power(trueQ2,2)*trueX)) / (alpha2*Yc);
		event_weight_alphaPt2 = event_weight_alphaPt2 * (mbToGeV_m2/(Lint*Q2binwidth*xbinwidth*pt2binwidth*alpha_spec_binwidth));
		//filling all alpha bins
		h_HERA_Q2_10_13_x007_009_alpha[alpha_bin_index]->Fill(pt2, event_weight_alphaPt2 );
		h_alpha_spec_everybin[alpha_bin_index]->Fill( alpha_spec );

	}

	output->Write();
	output->Close();


}