#include "RiceStyle.h"

using namespace std;

void plotBasics(){

	TString name[50];
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
	name[16] = "_DMPJet_eD";
	name[17] = "_DMPJet_eD_ionframe";
	
	name[18] = "highpf_JpsiNodecay_eT";
	name[19] = "highpf_JpsiNodecay_eT_ionframe";

	name[20] = "highpf_JpsiNodecay_eHe3";
	name[21] = "highpf_JpsiNodecay_eHe3_ionframe";

	name[22] = "highpf_JpsiNodecay_eAlpha";
	name[23] = "highpf_JpsiNodecay_eAlpha_ionframe";

	name[24] = "highpf_JpsiNodecay_eLi";
	name[25] = "highpf_JpsiNodecay_eLi_ionframe";

	name[26] = "highpf_JpsiNodecay_eC";
	name[27] = "highpf_JpsiNodecay_eC_ionframe";

	name[28] = "highpf_JpsiNodecay_eCa";
	name[29] = "highpf_JpsiNodecay_eCa_ionframe";

	name[30] = "highpf_JpsiNodecay_ePb";
	name[31] = "highpf_JpsiNodecay_ePb_ionframe";

	name[32] = "fixpf_JpsiGeneral_eD";
	name[33] = "fixpf_JpsiGeneral_eD_ionframe";
	
	name[34] = "highpf_JpsiGeneral_eD";
	name[35] = "highpf_JpsiGeneral_eD_ionframe";
	
	int start_i = 32;
	int end_i = 36;

	TFile* file[50];
	for(int i = start_i; i < end_i; i++){file[i] = new TFile("../rootfiles/"+name[i]+".root");}

	TH1D* pz_corr[50];
	TH1D* energy_corr[50];
	
	TH2D* energyVsQ2_2Dcorr[50];
	TH2D* energyVsW2_2Dcorr[50];
	TH2D* energyVsX_2Dcorr[50];
	TH2D* energyVsY_2Dcorr[50];
	TH2D* energyVsNu_2Dcorr[50];
	TH2D* energyVsPf_2Dcorr[50];
	TH2D* energyVsPtf_2Dcorr[50];
	TH2D* energyVsProcess_2Dcorr[50];

	for(int i = start_i; i < end_i; i++){

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

	TCanvas* c1[50];

	
	for(int i = start_i; i < end_i; i++){

		double x_range = 30.;
		if( i%2 != 0 ) {x_range = 1.0;}

		c1[i] = new TCanvas(Form("c1_%d",i),name[i],1,1,700,700);
		c1[i]->Divide(3,3);
		c1[i]->cd(1);
		gPad->SetLogy(1);
		pz_corr[i]->SetTitle("p_{z,in} - p_{z,out}");
		pz_corr[i]->GetXaxis()->SetRangeUser(-x_range,x_range);
		pz_corr[i]->Draw("");
		
		c1[i]->cd(2);
		gPad->SetLogy(1);
		energy_corr[i]->SetTitle("E_{in} - E_{out}");
		energy_corr[i]->GetXaxis()->SetRangeUser(-x_range,x_range);
		energy_corr[i]->Draw("");

		c1[i]->cd(3);
		gPad->SetLogz();
		energyVsQ2_2Dcorr[i]->SetTitle("Q2 vs E_{in}-E_{out}");
		energyVsQ2_2Dcorr[i]->GetXaxis()->SetRangeUser(-x_range,x_range);
		energyVsQ2_2Dcorr[i]->Draw("colz");

		c1[i]->cd(4);
		gPad->SetLogz();
		energyVsW2_2Dcorr[i]->SetTitle("W2 vs E_{in}-E_{out}");
		energyVsW2_2Dcorr[i]->GetXaxis()->SetRangeUser(-x_range,x_range);
		energyVsW2_2Dcorr[i]->Draw("colz");

		c1[i]->cd(5);
		gPad->SetLogz();
		energyVsX_2Dcorr[i]->SetTitle("x vs E_{in}-E_{out}");
		energyVsX_2Dcorr[i]->GetXaxis()->SetRangeUser(-x_range,x_range);
		energyVsX_2Dcorr[i]->Draw("colz");

		c1[i]->cd(6);
		gPad->SetLogz();
		energyVsY_2Dcorr[i]->SetTitle("y vs E_{in}-E_{out}");
		energyVsY_2Dcorr[i]->GetXaxis()->SetRangeUser(-x_range,x_range);
		energyVsY_2Dcorr[i]->Draw("colz");

		c1[i]->cd(7);
		gPad->SetLogz();
		energyVsPf_2Dcorr[i]->SetTitle("pf vs E_{in}-E_{out}");
		energyVsPf_2Dcorr[i]->GetXaxis()->SetRangeUser(-x_range,x_range);
		energyVsPf_2Dcorr[i]->GetYaxis()->SetRangeUser(-0.1,0.8);
		energyVsPf_2Dcorr[i]->Draw("colz");

		c1[i]->cd(8);
		gPad->SetLogz();
		energyVsPtf_2Dcorr[i]->SetTitle("ptf vs E_{in}-E_{out}");
		energyVsPtf_2Dcorr[i]->GetYaxis()->SetRangeUser(-0.1,0.8);
		energyVsPtf_2Dcorr[i]->GetXaxis()->SetRangeUser(-x_range,x_range);
		energyVsPtf_2Dcorr[i]->Draw("colz");

		c1[i]->cd(9);
		gPad->SetLogz();
		energyVsProcess_2Dcorr[i]->SetTitle("process vs E_{in}-E_{out}");
		energyVsProcess_2Dcorr[i]->GetXaxis()->SetRangeUser(-x_range,x_range);
		energyVsProcess_2Dcorr[i]->Draw("colz");

		c1[i]->Print("../figures/momentum_conservation/JpsiGeneral/"+name[i]+".pdf");

	} 

	

}