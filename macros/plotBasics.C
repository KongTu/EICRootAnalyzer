#include "RiceStyle.h"

using namespace std;

void plotBasics(){

	TFile* file = new TFile("wrongpf_JpsiNodecay_eP.root");
	TH1D* pz_corr = (TH1D*) file->Get("pz_corr");
	TH1D* energy_corr = (TH1D*) file->Get("energy_corr");

	TCanvas* c1 = new TCanvas("c1","",1,1,800,400);
	c1->Divide(2,1);
	c1->cd(1);
	pz_corr->SetTitle("wrongpf_eP");
	//pz_corr->GetXaxis()->SetRangeUser(240,270);
	pz_corr->Draw("");
	
	c1->cd(2);
	energy_corr->SetTitle("wrongpf_eP");
	//energy_corr->GetXaxis()->SetRangeUser(280,320);
	energy_corr->Draw("");

	c1->Print("pz_corr_wrongpf_eP.pdf");

}