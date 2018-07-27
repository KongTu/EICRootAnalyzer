#include "RiceStyle.h"

using namespace std;

void plotBasics(){

	TFile* file[20];
	file[0] = new TFile("../rootfiles/wrongpf_JpsiNodecay_eD.root");
	file[1] = new TFile("../rootfiles/wrongpf_JpsiNodecay_eD_ionframe.root");
	file[2] = new TFile("../rootfiles/fixpf_JpsiNodecay_eD.root");
	file[3] = new TFile("../rootfiles/fixpf_JpsiNodecay_eD_ionframe.root");
	file[4] = new TFile("../rootfiles/highpf_JpsiNodecay_eD.root");
	file[5] = new TFile("../rootfiles/highpf_JpsiNodecay_eD_ionframe.root");
	file[6] = new TFile("../rootfiles/zeropf_JpsiNodecay_eD.root");
	file[7] = new TFile("../rootfiles/zeropf_JpsiNodecay_eD_ionframe.root");
	
	TH1D* pz_corr[8];
	TH1D* energy_corr[8];
	
	TH2D* energyVsQ2_2Dcorr[8];
	TH2D* energyVsW2_2Dcorr[8];
	TH2D* energyVsX_2Dcorr[8];
	TH2D* energyVsY_2Dcorr[8];
	TH2D* energyVsNu_2Dcorr[8];
	TH2D* energyVsPf_2Dcorr[8];
	TH2D* energyVsPtf_2Dcorr[8];
	TH2D* energyVsProcess_2Dcorr[8];

	for(int i = 0; i < 8; i++){

		pz_corr[i] = (TH1D*) file[i]->Get("pz_corr");
		energy_corr[i] = (TH1D*) file[i]->Get("energy_corr");
		energyVsQ2_2Dcorr[i] = (TH2D*) file[i]->Get("energyVsQ2_2Dcorr");
		energyVsW2_2Dcorr[i] = (TH2D*) file[i]->Get("energyVsW2_2Dcorr");
		energyVsX_2Dcorr[i] = (TH2D*) file[i]->Get("energyVsX_2Dcorr");
		energyVsY_2Dcorr[i] = (TH2D*) file[i]->Get("energyVsY_2Dcorr");
		energyVsNu_2Dcorr[i] = (TH2D*) file[i]->Get("energyVsNu_2Dcorr");
		energyVsPf_2Dcorr[i] = (TH2D*) file[i]->Get("energyVsPf_2Dcorr");
		energyVsPtf_2Dcorr[i] = (TH2D*) file[i]->Get("energyVsPtf_2Dcorr");
		energyVsProcess_2Dcorr[i] = (TH2D*) file[i]->Get("energyVsProcess_2Dcorr");
	}

	TCanvas* c1 = new TCanvas("c1","",1,1,700,700);
	c1->Divide(3,3);
	c1->cd(1);
	pz_corr[0]->SetTitle("p_{z,in} - p_{z,out}");
	//pz_corr->GetXaxis()->SetRangeUser(240,270);
	pz_corr[0]->Draw("");
	
	c1->cd(2);
	energy_corr[0]->SetTitle("E_{in} - E_{out}");
	//energy_corr->GetXaxis()->SetRangeUser(280,320);
	energy_corr[0]->Draw("");

	c1->cd(3);
	gPad->SetLogz();
	energyVsQ2_2Dcorr[0]->SetTitle("Q2 vs E_{in}-E_{out}");
	//energyVsQ2_2Dcorr->GetXaxis()->SetRangeUser(280,320);
	energyVsQ2_2Dcorr[0]->Draw("colz");

	c1->cd(4);
	gPad->SetLogz();
	energyVsW2_2Dcorr[0]->SetTitle("W2 vs E_{in}-E_{out}");
	//energyVsW2_2Dcorr->GetXaxis()->SetRangeUser(280,320);
	energyVsW2_2Dcorr[0]->Draw("colz");

	c1->cd(5);
	gPad->SetLogz();
	energyVsX_2Dcorr[0]->SetTitle("x vs E_{in}-E_{out}");
	//energyVsX_2Dcorr->GetXaxis()->SetRangeUser(280,320);
	energyVsX_2Dcorr[0]->Draw("colz");

	c1->cd(6);
	gPad->SetLogz();
	energyVsY_2Dcorr[0]->SetTitle("y vs E_{in}-E_{out}");
	//energyVsY_2Dcorr->GetXaxis()->SetRangeUser(280,320);
	energyVsY_2Dcorr[0]->Draw("colz");

	c1->cd(7);
	gPad->SetLogz();
	energyVsPf_2Dcorr[0]->SetTitle("pf vs E_{in}-E_{out}");
	energyVsPf_2Dcorr->GetYaxis()->SetRangeUser(-0.3,0.3);
	energyVsPf_2Dcorr[0]->Draw("colz");

	c1->cd(8);
	gPad->SetLogz();
	energyVsPtf_2Dcorr[0]->SetTitle("ptf vs E_{in}-E_{out}");
	energyVsPtf_2Dcorr->GetYaxis()->SetRangeUser(-0.3,0.3);
	energyVsPtf_2Dcorr[0]->Draw("colz");

	c1->cd(9);
	gPad->SetLogz();
	energyVsProcess_2Dcorr[0]->SetTitle("process vs E_{in}-E_{out}");
	//energyVsProcess_2Dcorr->GetXaxis()->SetRangeUser(280,320);
	energyVsProcess_2Dcorr[0]->Draw("colz");

	c1->Print("test.pdf");

}