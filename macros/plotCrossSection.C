#include "RiceStyle.h"

using namespace std;

void plotCrossSection(){

	gStyle->SetErrorX(0);

	TFile* file = new TFile("../rootfiles/fixpf_Jpsinodecay_EvtParticlePlotter_eD.root");

	TH2D* TvsPt = (TH2D*) file->Get("TvsPt");

	double N_Jpsi[10];
	double bin_number[]={1,40,41,80,81,120,121,160,161,200};
	double bin_hist[] = {0.0,1.0,2.0,3.0,4.0,5.0};
	TH1D* dsigdT = new TH1D("dsigdT",";|T|(GeV/c^{2})", 5, bin_hist);
	for(int ibin = 0; ibin < 5; ibin++){

		double yield = TvsPt->ProjectionX(Form("pt_distribution_%d", ibin), bin_number[2*ibin], bin_number[2*ibin+1])->GetEntries();
		
		dsigdT->SetBinContent(5-ibin, yield);
		dsigdT->SetBinError(5-ibin, sqrt(yield));
	}

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetLogy(1);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	//gStyle->SetPadBorderMode(0.1);
	//gStyle->SetOptTitle(0);

	TH1D* base2 = makeHist("base2", "", "|T| (GeV^{2})", "dN/dT (GeV^{-2})", 100,0,5,kBlack);
	base2->GetYaxis()->SetRangeUser(0.1, 100000);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.1,1.25);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.6);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.6);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.6);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.6);
	base2->GetXaxis()->SetNdivisions(4,6,0);
	base2->GetYaxis()->SetNdivisions(4,6,0);

	base2->Draw();

	dsigdT->SetMarkerStyle(20);
	dsigdT->Draw("Psame");
	
}