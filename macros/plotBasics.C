#include "RiceStyle.h"

using namespace std;

void plotBasics(){

	TString name[20];
	name[0] = "wrongpf_JpsiNodecay_eD";
	name[1] = "wrongpf_JpsiNodecay_eD_ionframe";
	name[2] = "fixpf_JpsiNodecay_eD";
	name[3] = "fixpf_JpsiNodecay_eD_ionframe";
	name[4] = "highpf_JpsiNodecay_eD";
	name[5] = "highpf_JpsiNodecay_eD_ionframe";
	name[6] = "zeropf_JpsiNodecay_eD";
	name[7] = "zeropf_JpsiNodecay_eD_ionframe";
	name[8] = "fixpf_JpsiNodecay_ePb";
	name[9] = "fixpf_JpsiNodecay_ePb_ionframe";
	name[10] = "wrongpf_Jpsilept_eD";
	name[11] = "wrongpf_Jpsilept_eD_ionframe";
	name[12] = "fixpf_Jpsilept_eD";
	name[13] = "fixpf_Jpsilept_eD_ionframe";
	name[14] = "highpf_Jpsilept_eD";
	name[15] = "highpf_Jpsilept_eD_ionframe";

	TFile* file[20];
	for(int i = 0; i < 16; i++){file[i] = new TFile("../rootfiles/"+name[i]+".root");}

	TH1D* pz_corr[20];
	TH1D* energy_corr[20];
	
	TH2D* energyVsQ2_2Dcorr[20];
	TH2D* energyVsW2_2Dcorr[20];
	TH2D* energyVsX_2Dcorr[20];
	TH2D* energyVsY_2Dcorr[20];
	TH2D* energyVsNu_2Dcorr[20];
	TH2D* energyVsPf_2Dcorr[20];
	TH2D* energyVsPtf_2Dcorr[20];
	TH2D* energyVsProcess_2Dcorr[20];

	for(int i = 0; i < 16; i++){

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

	TCanvas* c1[20];

	for(int i = 0; i < 16; i++){

		c1[i] = new TCanvas("c1",name[i],1,1,700,700);
		c1[i]->Divide(3,3);
		c1[i]->cd(1);
		pz_corr[i]->SetTitle("p_{z,in} - p_{z,out}");
		pz_corr->GetXaxis()->SetRangeUser(-5,5);
		pz_corr[i]->Draw("");
		
		c1[i]->cd(2);
		energy_corr[i]->SetTitle("E_{in} - E_{out}");
		energy_corr->GetXaxis()->SetRangeUser(-5,5);
		energy_corr[i]->Draw("");

		c1[i]->cd(3);
		gPad->SetLogz();
		energyVsQ2_2Dcorr[i]->SetTitle("Q2 vs E_{in}-E_{out}");
		energyVsQ2_2Dcorr->GetXaxis()->SetRangeUser(-5,5);
		energyVsQ2_2Dcorr[i]->Draw("colz");

		c1[i]->cd(4);
		gPad->SetLogz();
		energyVsW2_2Dcorr[i]->SetTitle("W2 vs E_{in}-E_{out}");
		energyVsW2_2Dcorr->GetXaxis()->SetRangeUser(-5,5);
		energyVsW2_2Dcorr[i]->Draw("colz");

		c1[i]->cd(5);
		gPad->SetLogz();
		energyVsX_2Dcorr[i]->SetTitle("x vs E_{in}-E_{out}");
		energyVsX_2Dcorr->GetXaxis()->SetRangeUser(-5,5);
		energyVsX_2Dcorr[i]->Draw("colz");

		c1[i]->cd(6);
		gPad->SetLogz();
		energyVsY_2Dcorr[i]->SetTitle("y vs E_{in}-E_{out}");
		energyVsY_2Dcorr->GetXaxis()->SetRangeUser(-5,5);
		energyVsY_2Dcorr[i]->Draw("colz");

		c1[i]->cd(7);
		gPad->SetLogz();
		energyVsPf_2Dcorr[i]->SetTitle("pf vs E_{in}-E_{out}");
		energyVsPf_2Dcorr->GetXaxis()->SetRangeUser(-5,5);
		energyVsPf_2Dcorr[i]->GetYaxis()->SetRangeUser(-0.1,0.3);
		energyVsPf_2Dcorr[i]->Draw("colz");

		c1[i]->cd(8);
		gPad->SetLogz();
		energyVsPtf_2Dcorr[i]->SetTitle("ptf vs E_{in}-E_{out}");
		energyVsPtf_2Dcorr[i]->GetYaxis()->SetRangeUser(-0.1,0.3);
		energyVsPtf_2Dcorr->GetXaxis()->SetRangeUser(-5,5);
		energyVsPtf_2Dcorr[i]->Draw("colz");

		c1[i]->cd(9);
		gPad->SetLogz();
		energyVsProcess_2Dcorr[i]->SetTitle("process vs E_{in}-E_{out}");
		energyVsProcess_2Dcorr->GetXaxis()->SetRangeUser(-5,5);
		energyVsProcess_2Dcorr[i]->Draw("colz");

		c1[i]->Print(name[i]+".pdf");

	} 

	

}