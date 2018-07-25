#include "RiceStyle.h"

using namespace std;

void plotBasics(){

	TFile* file = new TFile("test_fixpf.root");
	TH2D* pz_corr = (TH2D*) file->Get("pz_corr");
	TH2D* energy_corr = (TH2D*) file->Get("energy_corr");

	TCanvas* c1 = new TCanvas("c1","",1,1,800,400);
	c1->Divide(2,1);
	c1->cd(1);
	pz_corr->SetTitle("fixpf");
	pz_corr->GetXaxis()->SetRangeUser(240,270);
	pz_corr->GetYaxis()->SetRangeUser(240,270);
	pz_corr->Draw("colz");
	
	c1->cd(2);
	energy_corr->SetTitle("fixpf");
	energy_corr->GetXaxis()->SetRangeUser(280,320);
	energy_corr->GetYaxis()->SetRangeUser(280,320);
	energy_corr->Draw("colz");

	c1->Print("pz_corr.pdf");

}