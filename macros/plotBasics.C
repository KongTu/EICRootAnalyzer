#include "RiceStyle.h"

using namespace std;

void plotBasics(){

	TString name[8];
	name[0] = "wrongpf_JpsiNodecay_eD";
	name[1] = "wrongpf_JpsiNodecay_eD_ionframe";
	name[2] = "fixpf_JpsiNodecay_eD";
	name[3] = "fixpf_JpsiNodecay_eD_ionframe";
	name[4] = "highpf_JpsiNodecay_eD";
	name[5] = "highpf_JpsiNodecay_eD_ionframe";
	name[6] = "zeropf_JpsiNodecay_eD";
	name[7] = "zeropf_JpsiNodecay_eD_ionframe";

	TFile* file[20];
	for(int i = 0; i < 8; i++){file[i] = new TFile("../rootfiles/"+name[i]+".root");}

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

	TCanvas* c1[8];

	for(int i = 0; i < 8; i++){

		c1[i] = new TCanvas("c1","",1,1,700,700);
		c1[i]->Divide(3,3);
		c1[i]->cd(1);
		pz_corr[i]->SetTitle("p_{z,in} - p_{z,out}");
		//pz_corr->GetXaxis()->SetRangeUser(240,270);
		pz_corr[i]->Draw("");
		
		c1[i]->cd(2);
		energy_corr[i]->SetTitle("E_{in} - E_{out}");
		//energy_corr->GetXaxis()->SetRangeUser(280,320);
		energy_corr[i]->Draw("");

		c1[i]->cd(3);
		gPad->SetLogz();
		energyVsQ2_2Dcorr[i]->SetTitle("Q2 vs E_{in}-E_{out}");
		//energyVsQ2_2Dcorr->GetXaxis()->SetRangeUser(280,320);
		energyVsQ2_2Dcorr[i]->Draw("colz");

		c1[i]->cd(4);
		gPad->SetLogz();
		energyVsW2_2Dcorr[i]->SetTitle("W2 vs E_{in}-E_{out}");
		//energyVsW2_2Dcorr->GetXaxis()->SetRangeUser(280,320);
		energyVsW2_2Dcorr[i]->Draw("colz");

		c1[i]->cd(5);
		gPad->SetLogz();
		energyVsX_2Dcorr[i]->SetTitle("x vs E_{in}-E_{out}");
		//energyVsX_2Dcorr->GetXaxis()->SetRangeUser(280,320);
		energyVsX_2Dcorr[i]->Draw("colz");

		c1[i]->cd(6);
		gPad->SetLogz();
		energyVsY_2Dcorr[i]->SetTitle("y vs E_{in}-E_{out}");
		//energyVsY_2Dcorr->GetXaxis()->SetRangeUser(280,320);
		energyVsY_2Dcorr[i]->Draw("colz");

		c1[i]->cd(7);
		gPad->SetLogz();
		energyVsPf_2Dcorr[i]->SetTitle("pf vs E_{in}-E_{out}");
		energyVsPf_2Dcorr[i]->GetYaxis()->SetRangeUser(-0.1,0.3);
		energyVsPf_2Dcorr[i]->Draw("colz");

		c1[i]->cd(8);
		gPad->SetLogz();
		energyVsPtf_2Dcorr[i]->SetTitle("ptf vs E_{in}-E_{out}");
		energyVsPtf_2Dcorr[i]->GetYaxis()->SetRangeUser(-0.1,0.3);
		energyVsPtf_2Dcorr[i]->Draw("colz");

		c1[i]->cd(9);
		gPad->SetLogz();
		energyVsProcess_2Dcorr[i]->SetTitle("process vs E_{in}-E_{out}");
		//energyVsProcess_2Dcorr->GetXaxis()->SetRangeUser(280,320);
		energyVsProcess_2Dcorr[i]->Draw("colz");

		c1[i]->Print(name[i]+".pdf");

	} 

	

}