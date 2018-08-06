#include <iostream>
#include <vector>
#include <sstream>
#include <string>

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

//TH2D DIS kinematics:
TH2D* Q2VsX = new TH2D("Q2VsX",";x;Q^{2}",10000,0.00001,1,2000,0,200);

//TH2D:
TH2D* PtVsEta_process_91_proton = new TH2D("PtVsEta_process_91_proton",";#eta;p_{T} (GeV)", 2000, -20,20, 200, 0,10);
TH2D* PtVsEta_process_91_neutron = new TH2D("PtVsEta_process_91_neutron",";#eta;p_{T} (GeV)", 2000, -20,20, 200, 0,10);

TH2D* PtVsEta_process_93_proton = new TH2D("PtVsEta_process_93_proton",";#eta;p_{T} (GeV)", 2000, -20,20, 200, 0,10);
TH2D* PtVsEta_process_93_neutron = new TH2D("PtVsEta_process_93_neutron",";#eta;p_{T} (GeV)", 2000, -20,20, 200, 0,10);

TH2D* PtVsPt_process_91_protonVsJpsi = new TH2D("PtVsPt_process_91_protonVsJpsi",";p_{T} (GeV);p_{T} (GeV)", 200, 0,10, 200, 0,10);
TH2D* PtVsPt_process_93_protonVsJpsi = new TH2D("PtVsPt_process_93_protonVsJpsi",";p_{T} (GeV);p_{T} (GeV)", 200, 0,10, 200, 0,10);

TH2D* AngleVsMom_process_91_proton = new TH2D("AngleVsMom_process_91_proton",";p (GeV);#theta (mrad)",2500,0,250,300,0,0.03);
TH2D* AngleVsMom_process_91_neutron = new TH2D("AngleVsMom_process_91_neutron",";p (GeV);#theta (mrad)",2500,0,250,300,0,0.03);

TH2D* AngleVsMom_process_93_proton = new TH2D("AngleVsMom_process_93_proton",";p (GeV);#theta (mrad)",2500,0,250,300,0,0.03);
TH2D* AngleVsMom_process_93_neutron = new TH2D("AngleVsMom_process_93_neutron",";p (GeV);#theta (mrad)",2500,0,250,300,0,0.03);


//TH1D:
TH1D* PtDist_process_91 = new TH1D("PtDist_process_91",";PtDist_process_91", 200, 0,10);
TH1D* PhiDist_process_91 = new TH1D("PhiDist_process_91",";PhiDist_process_91", 200, 0,10);
TH1D* EtaDist_process_91 = new TH1D("EtaDist_process_91",";EtaDist_process_91", 2000, -20,20);

TH1D* PtDist_process_93 = new TH1D("PtDist_process_93",";PtDist_process_93", 200, 0,10);
TH1D* PhiDist_process_93 = new TH1D("PhiDist_process_93",";PhiDist_process_93", 200, 0,10);
TH1D* EtaDist_process_93 = new TH1D("EtaDist_process_93",";EtaDist_process_93", 2000, -20,20);

TH1D* PtDist_process_91_Jpsi = new TH1D("PtDist_process_91_Jpsi",";PtDist_process_91_Jpsi", 200, 0,10);
TH1D* PhiDist_process_91_Jpsi = new TH1D("PhiDist_process_91_Jpsi",";PhiDist_process_91_Jpsi", 200, 0,10);
TH1D* EtaDist_process_91_Jpsi = new TH1D("EtaDist_process_91_Jpsi",";EtaDist_process_91_Jpsi", 2000, -20,20);

TH1D* PtDist_process_93_Jpsi = new TH1D("PtDist_process_93_Jpsi",";PtDist_process_93_Jpsi", 200, 0,10);
TH1D* PhiDist_process_93_Jpsi = new TH1D("PhiDist_process_93_Jpsi",";PhiDist_process_93_Jpsi", 200, 0,10);
TH1D* EtaDist_process_93_Jpsi = new TH1D("EtaDist_process_93_Jpsi",";EtaDist_process_93_Jpsi", 2000, -20,20);

TH1D* PtDist_process_91_proton = new TH1D("PtDist_process_91_proton",";PtDist_process_91_proton", 200, 0,10);
TH1D* PhiDist_process_91_proton = new TH1D("PhiDist_process_91_proton",";PhiDist_process_91_proton", 200, 0,10);
TH1D* EtaDist_process_91_proton = new TH1D("EtaDist_process_91_proton",";EtaDist_process_91_proton", 2000, -20,20);

TH1D* PtDist_process_93_proton = new TH1D("PtDist_process_93_proton",";PtDist_process_93_proton", 200, 0,10);
TH1D* PhiDist_process_93_proton = new TH1D("PhiDist_process_93_proton",";PhiDist_process_93_proton", 200, 0,10);
TH1D* EtaDist_process_93_proton = new TH1D("EtaDist_process_93_proton",";EtaDist_process_93_proton", 2000, -20,20);

TH1D* PtDist_process_91_neutron = new TH1D("PtDist_process_91_neutron",";PtDist_process_91_neutron", 200, 0,10);
TH1D* PhiDist_process_91_neutron = new TH1D("PhiDist_process_91_neutron",";PhiDist_process_91_neutron", 200, 0,10);
TH1D* EtaDist_process_91_neutron = new TH1D("EtaDist_process_91_neutron",";EtaDist_process_91_neutron", 2000, -20,20);

TH1D* PtDist_process_93_neutron = new TH1D("PtDist_process_93_neutron",";PtDist_process_93_neutron", 200, 0,10);
TH1D* PhiDist_process_93_neutron = new TH1D("PhiDist_process_93_neutron",";PhiDist_process_93_neutron", 200, 0,10);
TH1D* EtaDist_process_93_neutron = new TH1D("EtaDist_process_93_neutron",";EtaDist_process_93_neutron", 2000, -20,20);
