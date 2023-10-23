//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

void norm_hist( TH1D *h1, TH1D *h1_n);

Int_t plots_cghg_vs_simtel(){

  TString fileN01;
  TString fileN02;
  //fileN01 = "./hist.root";
  fileN01 = "./hist1km.root";
  //fileN02 = "../pyeventio_example/hist_fast_proton_nsb_1x_33837ev.root";
  fileN02 = "../pyeventio_example/hist_fast_proton_nsb_1x.root";
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
  //
  TH1D *h1_01_x_n = new TH1D();
  TH1D *h1_02_x_n = new TH1D();
  TH1D *h1_01_y_n = new TH1D();
  TH1D *h1_02_y_n = new TH1D();
  TH1D *h1_01_az_n = new TH1D();
  TH1D *h1_02_az_n = new TH1D();
  TH1D *h1_01_al_n = new TH1D();
  TH1D *h1_02_al_n = new TH1D();
  h1_01_x_n->SetNameTitle("h1_01_x_n","h1_01_x_n");
  h1_02_x_n->SetNameTitle("h1_02_x_n","h1_02_x_n");
  h1_01_y_n->SetNameTitle("h1_01_y_n","h1_01_y_n");
  h1_02_y_n->SetNameTitle("h1_02_y_n","h1_02_y_n");
  h1_01_az_n->SetNameTitle("h1_01_az_n","h1_01_az_n");
  h1_02_az_n->SetNameTitle("h1_02_az_n","h1_02_az_n");
  h1_01_al_n->SetNameTitle("h1_01_al_n","h1_01_al_n");
  h1_02_al_n->SetNameTitle("h1_02_al_n","h1_02_al_n");
  norm_hist(h1_01_x,h1_01_x_n);
  norm_hist(h1_02_x,h1_02_x_n);
  norm_hist(h1_01_y,h1_01_y_n);
  norm_hist(h1_02_y,h1_02_y_n);
  norm_hist(h1_01_az,h1_01_az_n);
  norm_hist(h1_02_az,h1_02_az_n);
  norm_hist(h1_01_al,h1_01_al_n);
  norm_hist(h1_02_al,h1_02_al_n);
  //
  //  
  h1_01_x->SetLineColor(kBlack);
  h1_01_x->SetLineWidth(3.0);
  h1_02_x->SetLineColor(kRed);
  h1_02_x->SetMarkerColor(kRed);
  h1_02_x->SetLineWidth(3.0);
  //  
  h1_01_x_n->SetLineColor(kBlack);
  h1_01_x_n->SetLineWidth(3.0);
  h1_02_x_n->SetLineColor(kRed);
  h1_02_x_n->SetMarkerColor(kRed);
  h1_02_x_n->SetLineWidth(3.0);
  //
  h1_01_y->SetLineColor(kBlack);
  h1_01_y->SetLineWidth(3.0);
  h1_02_y->SetLineColor(kRed);
  h1_02_y->SetMarkerColor(kRed);
  h1_02_y->SetLineWidth(3.0);
  //
  h1_01_y_n->SetLineColor(kBlack);
  h1_01_y_n->SetLineWidth(3.0);
  h1_02_y_n->SetLineColor(kRed);
  h1_02_y_n->SetMarkerColor(kRed);
  h1_02_y_n->SetLineWidth(3.0);
  //  
  h1_01_az->SetLineColor(kBlack);
  h1_01_az->SetLineWidth(3.0);
  h1_02_az->SetLineColor(kRed);
  h1_02_az->SetMarkerColor(kRed);
  h1_02_az->SetLineWidth(3.0);
  //  
  h1_01_az_n->SetLineColor(kBlack);
  h1_01_az_n->SetLineWidth(3.0);
  h1_02_az_n->SetLineColor(kRed);
  h1_02_az_n->SetMarkerColor(kRed);
  h1_02_az_n->SetLineWidth(3.0);
  //  
  h1_01_al->SetLineColor(kBlack);
  h1_01_al->SetLineWidth(3.0);
  h1_02_al->SetLineColor(kRed);
  h1_02_al->SetMarkerColor(kRed);
  h1_02_al->SetLineWidth(3.0);
  //  
  h1_01_al_n->SetLineColor(kBlack);
  h1_01_al_n->SetLineWidth(3.0);
  h1_02_al_n->SetLineColor(kRed);
  h1_02_al_n->SetMarkerColor(kRed);
  h1_02_al_n->SetLineWidth(3.0);
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
  h1_01_x_n->SetTitle("");
  h1_02_x_n->SetTitle("");
  h1_02_x_n->Draw();
  h1_02_x_n->GetXaxis()->SetTitle("Core x, km");
  h1_01_x_n->Draw("sames");
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
  h1_01_y_n->SetTitle("");
  h1_02_y_n->SetTitle("");
  h1_02_y_n->Draw();
  h1_02_y_n->GetXaxis()->SetTitle("Core y, km");
  h1_01_y_n->Draw("sames");
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
  h1_01_az_n->SetTitle("");
  h1_02_az_n->SetTitle("");
  h1_02_az_n->Draw();
  h1_02_az_n->GetXaxis()->SetTitle("Azimuth, deg");
  h1_01_az_n->Draw("sames");
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
  h1_01_al_n->SetTitle("");
  h1_02_al_n->SetTitle("");
  h1_02_al_n->Draw();
  h1_02_al_n->GetXaxis()->SetTitle("Altitude, deg");
  h1_01_al_n->Draw("sames");
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

void norm_hist( TH1D *h1, TH1D *h1_n){
  Int_t nBins = h1->GetNbinsX();
  Double_t xl = h1->GetBinLowEdge(1);
  Double_t xr = h1->GetBinLowEdge(nBins) + h1->GetBinWidth(nBins);
  Double_t norm = h1->GetEntries();
  cout<<"nBins "<<nBins<<endl
      <<"xl    "<<xl<<endl
      <<"xr    "<<xr<<endl;
  h1_n->SetBins(nBins, xl, xr);
  if(norm>0.0){
    for(Int_t i = 1;i<=nBins;i++){
      h1_n->SetBinContent(i,h1->GetBinContent(i)/norm);
    }
  }
}
