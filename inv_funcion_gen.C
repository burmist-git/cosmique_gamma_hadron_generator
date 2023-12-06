//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

Double_t get_EE_dist_invf(TRandom3 *rnd, Double_t val_Emin, Double_t val_Emax);

Double_t e_pdf(Double_t *x, Double_t *par){
  Double_t A  = par[0];
  Double_t B  = par[1];
  Double_t C  = par[2];
  Double_t val = A/(TMath::Power(x[0],B)) + C;
  return val;
}

Int_t inv_funcion_gen(){
  //
  Int_t nn = 10000000;
  //
  Double_t val_Emin = 10.0;
  Double_t val_Emax = 1000;
  //
  TRandom3 *rnd = new TRandom3(1231231);  
  //
  TH1D *h1 = new TH1D("h1","h1", 10000, val_Emin, val_Emax);

  for(Int_t j = 0;j<5;j++){
    cout<<j<<endl;
    for(Int_t i = 0;i<nn;i++)
      h1->Fill(get_EE_dist_invf(rnd,val_Emin,val_Emax));
  }
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //
  h1->SetLineColor(kBlack);
  h1->SetLineWidth(3.0);
  //
  h1->SetTitle("");
  h1->Draw();
  //
  const Int_t npar = 3;
  Double_t inParameters[npar];
  Double_t e_min  = val_Emin;
  Double_t e_max  = val_Emax;
  inParameters[1] = 2.0;
  inParameters[2] = 0.0;
  inParameters[0] = 13000.0;
  //
  TF1 *f_e_pdf = new TF1( "e_pdf", e_pdf, val_Emin, val_Emax, npar);
  f_e_pdf->SetParameters(inParameters);
  f_e_pdf->SetParName(0, "A");
  f_e_pdf->SetParName(1, "B");
  f_e_pdf->SetParName(2, "C");
  f_e_pdf->FixParameter(2,inParameters[2]);
  f_e_pdf->FixParameter(1,inParameters[1]);
  //
  h1->Fit("e_pdf","","",val_Emin, val_Emax);
  //
  return 0;
}

Double_t get_EE_dist_invf(TRandom3 *rnd,
			  Double_t val_Emin, Double_t val_Emax){
  Double_t C = 1.0/(1.0/val_Emin - 1.0/val_Emax);
  return 1.0/(1.0/val_Emin - rnd->Uniform()/C);
}

