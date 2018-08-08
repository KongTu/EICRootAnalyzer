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

#define PI            3.1415926

#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_DEUTERON 1.8624778138724238
#define MASS_TRITON   2.7937167208086358
#define MASS_HE3      2.7937167208086358
#define MASS_ALPHA    3.7249556277448477
#define MASS_LI6      5.5874334416172715
#define MASS_C12      11.174866883234543
#define MASS_CA40     37.249556277448477
#define MASS_XE131    121.99229680864376
#define MASS_AU197    183.45406466643374
#define MASS_PB208    193.69769264273208
// Histograms for our analysis.
TH1D* Ntrk_process_91 = new TH1D("Ntrk_process_91",";Ntrk_process_91", 100, 0, 100);
TH1D* Ntrk_process_93 = new TH1D("Ntrk_process_93",";Ntrk_process_93", 100, 0, 100);
TH1D* Ntrk_process_all = new TH1D("Ntrk_process_all",";Ntrk_process_all", 100, 0, 100);
TH1D statusHist("statusHist", "status distribution  ", 50, 0, 50 );

TH1D* Ntrk = new TH1D("Ntrk",";Ntrk", 100, 0, 100);
TH2D* pTvsThat = new TH2D("pTvsThat",";pT;t_hat", 1000,0,10,1000,-10,10);

//TH2D DIS kinematics:
TH2D* Q2VsX = new TH2D("Q2VsX",";x;Q^{2}",10000,0.00001,1,2000,0,200);
TH2D* W2VsFlux = new TH2D("W2VsFlux",";#Phi;W^{2}",5000,0,0.1,1000,0,10000);

//TH2D:
TH2D* PtVsEta_process_91_proton = new TH2D("PtVsEta_process_91_proton",";#eta;p_{T} (GeV)", 2000, -20,20, 200, 0,10);
TH2D* PtVsEta_process_91_neutron = new TH2D("PtVsEta_process_91_neutron",";#eta;p_{T} (GeV)", 2000, -20,20, 200, 0,10);

TH2D* PtVsEta_process_93_proton = new TH2D("PtVsEta_process_93_proton",";#eta;p_{T} (GeV)", 2000, -20,20, 200, 0,10);
TH2D* PtVsEta_process_93_neutron = new TH2D("PtVsEta_process_93_neutron",";#eta;p_{T} (GeV)", 2000, -20,20, 200, 0,10);

TH2D* PtVsPt_process_91_protonVsJpsi = new TH2D("PtVsPt_process_91_protonVsJpsi",";p_{T} (GeV);p_{T} (GeV)", 200, 0,10, 200, 0,10);
TH2D* PtVsPt_process_93_protonVsJpsi = new TH2D("PtVsPt_process_93_protonVsJpsi",";p_{T} (GeV);p_{T} (GeV)", 200, 0,10, 200, 0,10);

TH2D* AngleVsMom_process_91_proton = new TH2D("AngleVsMom_process_91_proton",";p (GeV);#theta (mrad)",2500,0,250,300,0,30);
TH2D* AngleVsMom_process_91_neutron = new TH2D("AngleVsMom_process_91_neutron",";p (GeV);#theta (mrad)",2500,0,250,300,0,30);
TH2D* AngleVsMom_process_93_proton = new TH2D("AngleVsMom_process_93_proton",";p (GeV);#theta (mrad)",2500,0,250,300,0,30);
TH2D* AngleVsMom_process_93_neutron = new TH2D("AngleVsMom_process_93_neutron",";p (GeV);#theta (mrad)",2500,0,250,300,0,30);

TH2D* AngleVssNN_process_91_proton = new TH2D("AngleVssNN_process_91_proton",";s_{_{NN}} (GeV^{2});#theta (mrad)",2500,0,10,300,0,30);
TH2D* AngleVssNN_process_93_proton = new TH2D("AngleVssNN_process_93_proton",";s_{_{NN}} (GeV^{2});#theta (mrad)",2500,0,10,300,0,30);
TH2D* AngleVssNN_process_91_neutron = new TH2D("AngleVssNN_process_91_neutron",";s_{_{NN}} (GeV^{2});#theta (mrad)",2500,0,10,300,0,30);
TH2D* AngleVssNN_process_93_neutron = new TH2D("AngleVssNN_process_93_neutron",";s_{_{NN}} (GeV^{2});#theta (mrad)",2500,0,10,300,0,30);

TH2D* Q2VsJpsi_91 = new TH2D("Q2VsJpsi_91",";p_{T};Q^{2}",100,0,10,2000,0,200);
TH2D* Q2VsJpsi_93 = new TH2D("Q2VsJpsi_93",";p_{T};Q^{2}",100,0,10,2000,0,200);

TH2D* W2VsJpsi_91 = new TH2D("W2VsJpsi_91",";p_{T};W^{2}",100,0,10,1000,0,10000);
TH2D* W2VsJpsi_93 = new TH2D("W2VsJpsi_93",";p_{T};W^{2}",100,0,10,1000,0,10000);

TH2D* T_hatVsPt2 = new TH2D("T_hatVsPt2",";p^{2}_{T}-Q^{2} (GeV);T", 500,0,50,200,-5,0);

TH2D* TvsPt_91 = new TH2D("TvsPt_91",";p_{T} (GeV);T", 100,0,10,200,-5,0);
TH2D* TvsPt_93 = new TH2D("TvsPt_93",";p_{T} (GeV);T", 100,0,10,200,-5,0);

TH2D* tProtonVsPt_91 = new TH2D("tProtonVsPt_91",";p_{T} (GeV);T", 100,0,10,200,-5,0);
TH2D* tProtonVsPt_93 = new TH2D("tProtonVsPt_93",";p_{T} (GeV);T", 100,0,10,200,-5,0);

TH2D* ThatVssNN = new TH2D("ThatVssNN",";s_{_{NN}} (GeV^{2});T" ,200,0,20,200,-5,0);
TH2D* tdisVssNN = new TH2D("tdisVssNN",";s_{_{NN}} (GeV^{2});t" ,200,0,20,400,-5,5);

TH2D* tVsT = new TH2D("tVsT",";T;t" ,200,-5,0,400,-5,5);

TH2D* sNNvsPt_91 = new TH2D("sNNvsPt_91",";p_{T} (GeV);s_{_{NN}} (GeV^{2})", 100,0,10,200,0,20);

//TH1D event variables:
TH1D* E_CM = new TH1D("E_CM",";#sqrt{s_{_{NN}}} (GeV)", 200,0,10);
TH1D* T_dist = new TH1D("T_dist",";T", 200,-5,0);
TH1D* t_dist = new TH1D("t_dist",";t", 200,-5,5);

TH1D* t_proton_dist = new TH1D("t_proton_dist",";t",200,-5,5);
TH1D* W2 = new TH1D("W2",";W^{2} (GeV^{2})", 1000,0,10000);
TH1D* photonFlux = new TH1D("photonFlux",";#Phi", 5000,0,0.1);

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
