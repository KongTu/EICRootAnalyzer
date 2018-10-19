#include "RiceStyle.h"

using namespace std;

/*
cross section:

dsigma/dtdTdy vs X

- fill X distribution given the events you are interested in, with selection on t, T, y and etc...
- normalize by the range of such selections (of course we should know the average value too)
- normalize by integrated luminosity

*/


void plotCrossSection(){

	gStyle->SetErrorX(0);

	TFile* file = new TFile("../rootfiles/fixpf_Jpsinodecay_KickFinalStates_eD_nokick.root");
	TFile* file2 = new TFile("../rootfiles/fixpf_Jpsinodecay_KickFinalStates_eD_kick.root");

	TH2D* TvsPt = (TH2D*) file->Get("TvsPt");

	double N_Jpsi[10];
	double bin_number[]={1,40,41,80,81,120,121,160,161,200};
	double bin_hist[] = {0.0,1.0,2.0,3.0,4.0,5.0};
	TH1D* dsigdT = new TH1D("dsigdT",";|T|(GeV/c^{2})", 5, bin_hist);
	
	for(int ibin = 0; ibin < 5; ibin++){

		double raw_yield = TvsPt->ProjectionX(Form("pt_distribution_%d", ibin), bin_number[2*ibin], bin_number[2*ibin+1])->GetEntries();
		double Tweighted_yield = raw_yield/(dsigdT->GetBinCenter(5-ibin));

		dsigdT->SetBinContent(5-ibin,  Tweighted_yield );
		dsigdT->SetBinError(5-ibin, sqrt(raw_yield)/(dsigdT->GetBinCenter(5-ibin)));

	}


	TH2D* sNNvsPt = (TH2D*) file->Get("sNNvsPt");
	TH2D* sNNvsPt2 = (TH2D*) file2->Get("sNNvsPt");

	double bin_number_sNN[] = {1,10,11,20,21,30,31,40,41,50,51,60,61,70,71,80,81,90,91,100};
	double bin_hist_sNN[] = {0,1,2,3,4,5,6,7,8,9,10};
	TH1D* dNdTdy = new TH1D("dNdTdy",";s_{_{NN}}(GeV^{2})", 10, bin_hist_sNN);
	TH1D* dNdTdy2 = new TH1D("dNdTdy2",";s_{_{NN}}(GeV^{2})", 10, bin_hist_sNN);

	for(int ibin = 0; ibin < 10; ibin++){

		double raw_yield = sNNvsPt->ProjectionX(Form("pt_distribution_%d", ibin), bin_number_sNN[2*ibin], bin_number_sNN[2*ibin+1])->GetEntries();
		double differential_yield = raw_yield; //4 is the total unit of rapidity
		
		//cout << "raw_yield " << differential_yield << endl;

		dNdTdy->SetBinContent(ibin+1, differential_yield);
		dNdTdy->SetBinError(ibin+1, sqrt(raw_yield) );

		double raw_yield = sNNvsPt2->ProjectionX(Form("pt_distribution2_%d", ibin), bin_number_sNN[2*ibin], bin_number_sNN[2*ibin+1])->GetEntries();
		double differential_yield = raw_yield; //4 is the total unit of rapidity
		
		//cout << "raw_yield " << differential_yield << endl;

		dNdTdy2->SetBinContent(ibin+1, differential_yield);
		dNdTdy2->SetBinError(ibin+1, sqrt(raw_yield) );
	}

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetLogy(1);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	//gStyle->SetPadBorderMode(0.1);
	//gStyle->SetOptTitle(0);

	TH1D* base2 = makeHist("base2", "", "|T_{hard}| (GeV^{2})", "dN/dT (GeV^{-2})", 100,0,5,kBlack);
	base2->GetYaxis()->SetRangeUser(0.01, 500000);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.1,1.25);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.6);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.6);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.6);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.6);
	base2->GetXaxis()->SetNdivisions(4,6,0);
	base2->GetYaxis()->SetNdivisions(6,6,0);

	base2->Draw();

	dsigdT->SetMarkerStyle(20);
	dsigdT->Draw("Psame");

	TLatex* r42 = new TLatex(0.61, 0.85, "#gamma+D -> J/#psi+p+n");
    r42->SetNDC();
    r42->SetTextSize(23);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);

    TLatex* r43 = new TLatex(0.76,0.91, "BeAGLE");
    r43->SetNDC();
    r43->SetTextSize(0.04);
    
    TLatex* r44 = new TLatex(0.71,0.91, "");
    r44->SetNDC();
    r44->SetTextSize(22);
    r44->SetTextFont(53);

    r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");


    TCanvas* c2 = new TCanvas("c2","c2",600,600);
	gPad->SetLogy(1);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	//gStyle->SetPadBorderMode(0.1);
	//gStyle->SetOptTitle(0);

	TH1D* base3 = makeHist("base3", "", "s_{_{NN}} (GeV^{2})", "dN/dTdtdy (GeV^{-4})", 100,0,10,kBlack);
	base3->GetYaxis()->SetRangeUser(0.01, 5000000);
	base3->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base3,1.1,1.25);

	base3->GetYaxis()->SetTitleOffset(1.3);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.6);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.6);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.6);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.6);
	base3->GetXaxis()->SetNdivisions(4,6,0);
	base3->GetYaxis()->SetNdivisions(6,6,0);


	TF1* myFunc = new TF1("myfunc","[0]*TMath::Power(x,-14)",3,7);

	//dNdTdy->Fit("myfunc");
		base3->Draw();

	dNdTdy->SetMarkerStyle(20);
	dNdTdy->Draw("Psame");

	dNdTdy2->SetMarkerStyle(24);
	dNdTdy2->Draw("Psame");

	r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
	
}