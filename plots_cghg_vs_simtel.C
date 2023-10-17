//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

Int_t plots_cghg_vs_simtel(){

  TString fileN01;
  TString fileN02;
  fileN01 = "./hist.root";
  fileN02 = "../pyeventio_example/hist_fast_proton_nsb_1x_33837ev.root";
  TFile *f01 = new TFile(fileN01.Data());
  TFile *f02 = new TFile(fileN02.Data());
  //
  TH1D *h1_01_x = (TH1D*)f01->Get("h1_x1_int_cut");
  TH1D *h1_01_y = (TH1D*)f01->Get("h1_y1_int_cut");
  TH1D *h1_01_az = (TH1D*)f01->Get("h1_azimuth_deg");
  TH1D *h1_01_al = (TH1D*)f01->Get("h1_altitude_deg");
  //  
  TH1D *h1_02_x = (TH1D*)f02->Get("h1_xcore_km");
  TH1D *h1_02_y = (TH1D*)f02->Get("h1_ycore_km");
  TH1D *h1_02_az = (TH1D*)f02->Get("h1_azimuth_deg");
  TH1D *h1_02_al = (TH1D*)f02->Get("h1_altitude_deg");
  //  
  h1_01_x->SetLineColor(kBlack);
  h1_01_x->SetLineWidth(3.0);
  h1_02_x->SetLineColor(kRed);
  h1_02_x->SetLineWidth(3.0);
  //
  h1_01_y->SetLineColor(kBlack);
  h1_01_y->SetLineWidth(3.0);
  h1_02_y->SetLineColor(kRed);
  h1_02_y->SetLineWidth(3.0);
  //  
  h1_01_az->SetLineColor(kBlack);
  h1_01_az->SetLineWidth(3.0);
  h1_02_az->SetLineColor(kRed);
  h1_02_az->SetLineWidth(3.0);
  //  
  h1_01_al->SetLineColor(kBlack);
  h1_01_al->SetLineWidth(3.0);
  h1_02_al->SetLineColor(kRed);
  h1_02_al->SetLineWidth(3.0);
  //
  //
  TLegend *leg1 = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg1->AddEntry(h1_01_x, "Cosmique gamma hadron generator", "apl");
  leg1->AddEntry(h1_02_x, "Simtel (with at least one p.e.)", "apl");
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  gPad->SetGridx();
  gPad->SetGridy();
  //
  c1->SetRightMargin(0.01);
  c1->SetLeftMargin(0.1);
  c1->SetTopMargin(0.01);
  c1->SetBottomMargin(0.1);
  //  
  h1_01_x->SetTitle("");
  h1_02_x->SetTitle("");
  h1_02_x->Draw();
  h1_02_x->GetXaxis()->SetTitle("Core x, km");
  h1_01_x->Draw("sames");
  leg1->Draw("same");
  //
  //    
  TCanvas *c2 = new TCanvas("c2","c2",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  gPad->SetGridx();
  gPad->SetGridy();
  //
  c2->SetRightMargin(0.01);
  c2->SetLeftMargin(0.1);
  c2->SetTopMargin(0.01);
  c2->SetBottomMargin(0.1);
  //
  h1_01_y->SetTitle("");
  h1_02_y->SetTitle("");
  h1_02_y->Draw();
  h1_02_y->GetXaxis()->SetTitle("Core y, km");
  h1_01_y->Draw("sames");
  leg1->Draw("same");
  //
  //    
  TCanvas *c3 = new TCanvas("c3","c3",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  gPad->SetGridx();
  gPad->SetGridy();
  //
  c3->SetRightMargin(0.01);
  c3->SetLeftMargin(0.1);
  c3->SetTopMargin(0.01);
  c3->SetBottomMargin(0.1);
  //  
  h1_01_az->SetTitle("");
  h1_02_az->SetTitle("");
  h1_02_az->Draw();
  h1_02_az->GetXaxis()->SetTitle("Azimuth, deg");
  h1_01_az->Draw("sames");
  leg1->Draw("same");
  //
  //    
  TCanvas *c4 = new TCanvas("c4","c4",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  gPad->SetGridx();
  gPad->SetGridy();
  //
  c4->SetRightMargin(0.01);
  c4->SetLeftMargin(0.1);
  c4->SetTopMargin(0.01);
  c4->SetBottomMargin(0.1);
  //  
  h1_01_al->SetTitle("");
  h1_02_al->SetTitle("");
  h1_02_al->Draw();
  h1_02_al->GetXaxis()->SetTitle("Altitude, deg");
  h1_01_al->Draw("sames");
  //
  leg1->Draw("same");
  //
  //
  //
  //h1_01->GetXaxis()->SetTitle("Theta, deg");
  //h1_01->GetYaxis()->SetTitle("Distance to the surface, km");
  //
  //h1_1->GetXaxis()->SetTitle("Theta, deg");
  //h1_1->GetXaxis()->SetTitle("Phi, deg");
  //h1_1->GetXaxis()->SetTitle("Distance to the surface, km");
  //h1_1->GetXaxis()->SetTitle("Distance to terzina, km");
  //h1_1->GetXaxis()->SetTitle("Angle between terzina and proton track, deg");
  //h1_1->GetXaxis()->SetTitle("Photon density, 1/m^2");
  //h1_1->GetXaxis()->SetRangeUser(50.0,90.0);
  //h1_1->GetXaxis()->SetRangeUser(250.0,290.0);
  //
  return 0;
}
