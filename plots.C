//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

Int_t plots(){

  TString fileN;
  TString objectName;
  fileN = "./hist.root";

  TFile *f1 = new TFile(fileN.Data());

  //TH1D *h1_1 = (TH1D*)f1->Get("h1_theta_deg");
  //TH1D *h1_1 = (TH1D*)f1->Get("h1_phi_deg");
  //TH1D *h1_1 = (TH1D*)f1->Get("h1_distToEarth");
  //TH1D *h1_1 = (TH1D*)f1->Get("h1_distToTerzina");
  //TH1D *h1_1 = (TH1D*)f1->Get("h1_angleTrzinaTrk_deg");
  TH1D *h1_1 = (TH1D*)f1->Get("h1_nphotons_per_m2");
  TH2D *h2_1 = (TH2D*)f1->Get("h2_distToEarth_vs_theta_norm");
    
  TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);

  gPad->SetGridx();
  gPad->SetGridy();
  
  h1_1->SetLineColor(kBlack);
  h1_1->SetLineWidth(3.0);

  //h1_1->SetTitle("");
  //h1_1->Draw();
  h2_1->SetTitle("");
  h2_1->Draw("ZCOLOR");
  h2_1->GetXaxis()->SetTitle("Theta, deg");
  h2_1->GetYaxis()->SetTitle("Distance to the surface, km");

  //h1_1->GetXaxis()->SetTitle("Theta, deg");
  //h1_1->GetXaxis()->SetTitle("Phi, deg");
  //h1_1->GetXaxis()->SetTitle("Distance to the surface, km");
  //h1_1->GetXaxis()->SetTitle("Distance to terzina, km");
  //h1_1->GetXaxis()->SetTitle("Angle between terzina and proton track, deg");
  //h1_1->GetXaxis()->SetTitle("Photon density, 1/m^2");
  //h1_1->GetXaxis()->SetRangeUser(50.0,90.0);
  //h1_1->GetXaxis()->SetRangeUser(250.0,290.0);

  return 0;
}
