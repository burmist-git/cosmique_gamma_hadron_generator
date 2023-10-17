//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

Double_t get_intensity_of_primary_nucleons( Double_t ekin, Int_t N);
void get_DAMPE(TGraphErrors *gr, TString name_dat_file);

Double_t getFlux(TF1 *f_cpf, Double_t el, Double_t er);

void get_flux_p_PDG(TGraph *gr, TString name_dat_file);

Double_t cpf(Double_t *x, Double_t *par){
  Double_t A  = par[0];
  Double_t B  = par[1];
  Double_t C  = 0.0;
  Double_t val = A/(TMath::Power(x[0],B)) + C;
  return val;
}

Int_t intensity_of_primary_nucleons(){
  //
  TString fileN;
  //
  const Int_t n = 10000;
  Double_t I[n];
  Double_t etot[n];
  //
  Double_t etot_min = 1;      //in GeV (5 GeV)
  Double_t etot_max = 1000000; //in GeV (100 TeV)
  //
  Double_t N = 1;
  //
  TGraph *gr_proton = new TGraph();
  gr_proton->SetNameTitle("gr_proton","gr_proton");
  TGraphErrors *gr_DAMPE = new TGraphErrors();
  gr_DAMPE->SetNameTitle("gr_DAMPE","gr_DAMPE");
  TGraphErrors *gr_DAMPE_He = new TGraphErrors();
  gr_DAMPE_He->SetNameTitle("gr_DAMPE_He","gr_DAMPE_He");
  TGraphErrors *gr_diff = new TGraphErrors();
  gr_diff->SetNameTitle("gr_diff","gr_diff");
  TGraphErrors *gr_diff_PDG = new TGraphErrors();
  gr_diff_PDG->SetNameTitle("gr_diff_PDG","gr_diff_PDG");
  TGraphErrors *gr_diff_PDG_par = new TGraphErrors();
  gr_diff_PDG_par->SetNameTitle("gr_diff_PDG_par","gr_diff_PDG_par");
  //
  TGraph *gr_flux_p_PDG = new TGraph();
  gr_flux_p_PDG->SetNameTitle("gr_flux_p_PDG","gr_flux_p_PDG");  
  //
  TGraphErrors *gr_diff_DAMPE_fit_p_PDG = new TGraphErrors();
  gr_diff_DAMPE_fit_p_PDG->SetNameTitle("gr_diff_DAMPE_fit_p_PDG","gr_diff_DAMPE_fit_p_PDG");
  //
  for(Int_t i = 0;i<n;i++){
    etot[i] = etot_min + (etot_max - etot_min)/(n - 1)*i;
    I[i] = get_intensity_of_primary_nucleons(etot[i], N);
    gr_proton->SetPoint(gr_proton->GetN(),etot[i],I[i]);
  }
  //cout<<"10      "<<gr_proton->Eval(10)<<endl
  //  <<"100     "<<gr_proton->Eval(100)<<endl
  //  <<"1000    "<<gr_proton->Eval(1000)<<endl
  //  <<"10000   "<<gr_proton->Eval(10000)<<endl
  //  <<"100000  "<<gr_proton->Eval(100000)<<endl;
  get_DAMPE(gr_DAMPE,"DAMPE.dat");
  get_DAMPE(gr_DAMPE_He,"DAMPE_He.dat");
  get_flux_p_PDG(gr_flux_p_PDG, "flux_P.csv");  
  //
  Double_t ee;
  Double_t ff;
  Double_t ff_d;
  for(Int_t i = 0;i<gr_DAMPE->GetN();i++){
    gr_DAMPE->GetPoint(i,ee,ff_d);
    ff = get_intensity_of_primary_nucleons(ee, N);
    gr_diff->SetPoint(i,ee,(ff - ff_d)/ff_d);
    gr_diff_PDG->SetPoint(i,ee,(gr_flux_p_PDG->Eval(ee) - ff_d)/ff_d);
  }
  for(Int_t i = 0;i<gr_flux_p_PDG->GetN();i++){
    gr_flux_p_PDG->GetPoint(i,ee,ff_d);
    ff = get_intensity_of_primary_nucleons(ee, N);
    gr_diff_PDG_par->SetPoint(i,ee,(ff - ff_d)/ff_d);
  }
  //  
  TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,1200,600);
  c1->Divide(2,1);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  c1->cd(1);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  gPad->SetLogy();
  //
  //gr_proton->Draw("APL");
  //gr_DAMPE->Draw("APL");
  //
  TMultiGraph *mg = new TMultiGraph();
  gr_proton->SetMarkerStyle(1);
  gr_proton->SetLineWidth(2);
  gr_proton->SetMarkerColor(kBlack);
  gr_DAMPE->SetMarkerStyle(20);
  gr_DAMPE->SetMarkerColor(kRed+2);
  gr_DAMPE->SetLineWidth(2);
  gr_DAMPE->SetLineColor(kRed+2);
  gr_DAMPE_He->SetMarkerStyle(20);
  gr_DAMPE_He->SetMarkerColor(kBlue+2);
  gr_DAMPE_He->SetLineWidth(2);
  gr_DAMPE_He->SetLineColor(kBlue+2);
  //
  gr_flux_p_PDG->SetLineWidth(2);
  gr_flux_p_PDG->SetMarkerColor(kMagenta+2);
  gr_flux_p_PDG->SetLineColor(kMagenta+2);
  //
  mg->Add(gr_proton);
  mg->Add(gr_DAMPE);
  mg->Add(gr_flux_p_PDG);
  //mg->Add(gr_DAMPE_He);
  //
  mg->Draw("APL");
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr_DAMPE, "DAMPE Collaboration,Sci. Adv.2019;5:eaax3793", "apl");
  leg->AddEntry(gr_proton, "PDG parametrisation: 1.8*10^{4}/E^{2.7}", "apl");
  leg->AddEntry(gr_flux_p_PDG, "PDG data", "apl");
  leg->Draw();
  
  //
  mg->GetXaxis()->SetTitle("E, GeV");
  mg->GetYaxis()->SetTitle("dN/dE, m^{-2} s^{-1} sr^{-1} GeV^{-1}");
  //
  c1->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  //
  gr_diff->SetMarkerStyle(1);
  gr_diff->SetLineWidth(2);
  gr_diff->SetMarkerColor(kBlack);
  gr_diff->SetTitle("");
  //
  gr_diff_PDG->SetMarkerStyle(1);
  gr_diff_PDG->SetLineWidth(2);
  gr_diff_PDG->SetMarkerColor(kBlue+2);
  gr_diff_PDG->SetLineColor(kBlue+2);
  gr_diff_PDG->SetTitle("");
  //
  gr_diff_PDG_par->SetMarkerStyle(1);
  gr_diff_PDG_par->SetLineWidth(2);
  gr_diff_PDG_par->SetMarkerColor(kRed+2);
  gr_diff_PDG_par->SetLineColor(kRed+2);
  gr_diff_PDG_par->SetTitle("");
  //
  TMultiGraph *mg2 = new TMultiGraph();
  mg2->Add(gr_diff);
  mg2->Add(gr_diff_PDG);
  mg2->Add(gr_diff_PDG_par);
  mg2->Draw("APL");
  //
  mg2->GetYaxis()->SetTitle("Relative difference");
  mg2->GetXaxis()->SetTitle("E, GeV");
  //
  TLegend *leg2 = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg2->AddEntry(gr_diff, "DAMPE vs. PDG par.", "apl");
  leg2->AddEntry(gr_diff_PDG, "DAMPE vs. PDG data", "apl");
  leg2->AddEntry(gr_diff_PDG_par, "PDG par. vs. PDG data", "apl");
  leg2->Draw();
  //
  //

  TCanvas *c2 = new TCanvas("c2",fileN.Data(),10,10,1200,600);
  c2->Divide(2,1);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  c2->cd(1);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  gPad->SetLogy();
  //

  //Fit
  //
  const Int_t npar = 2;
  Double_t inParameters[npar];
  Double_t e_min  = etot_min;
  Double_t e_max  = etot_max;
  //Double_t f_max  = 7.0*1.0e-15;
  //Double_t f_min  = 1.0*1.0e-27;
  //inParameters[1] = 3.0;
  //inParameters[2] = 0.0;
  //inParameters[0] = f_max*TMath::Power(e_min,inParameters[1]);
  //inParameters[3] = 0.0;
  //
  //inParameters[0] = 1.8e+04;
  //inParameters[1] = 2.70;
  //inParameters[2] = 0.0;
  //
  //inParameters[0] = 8.98902e+03;
  //inParameters[1] = 2.66195e+00;
  //inParameters[2] = -1.91206e-10;
  //inParameters[2] = 0.0;
  inParameters[0] = 9.67854e+03;
  //inParameters[1] = 2.67400e+00;
  inParameters[1] = 2674.0/1000.0;
  //
  TF1 *f_cpf = new TF1( "f_cpf", cpf, e_min, e_max, npar);
  f_cpf->SetParameters(inParameters);
  f_cpf->SetParName(0, "A");
  f_cpf->SetParName(1, "B");
  //f_cpf->SetParName(2, "C");
  //f_cpf->FixParameter(0,inParameters[0]);
  //f_cpf->FixParameter(1,inParameters[1]);
  //f_cpf->FixParameter(2,inParameters[2]);
  gr_DAMPE->Fit("f_cpf","","",e_min, e_max);
  ///////////////////
  gr_DAMPE->SetTitle("Flux(E) = A/E^{B}");
  //gr_DAMPE->Draw("APL");
  //gr_DAMPE->GetXaxis()->SetTitle("E, GeV");
  //gr_DAMPE->GetYaxis()->SetTitle("dN/dE, m^{-2} s^{-1} sr^{-1} GeV^{-1}");
  //gr_flux_p_PDG->Draw("same APL");
  //
  TMultiGraph *mg3 = new TMultiGraph();
  mg3->SetTitle("Flux(E) = A/E^{B}");
  mg3->Add(gr_DAMPE);
  mg3->Add(gr_flux_p_PDG);
  mg3->Draw("APL");
  //
  mg3->GetXaxis()->SetTitle("E, GeV");
  mg3->GetYaxis()->SetTitle("dN/dE, m^{-2} s^{-1} sr^{-1} GeV^{-1}");
  //
  TLegend *leg3 = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg3->AddEntry(gr_DAMPE, "DAMPE", "apl");
  leg3->AddEntry(gr_flux_p_PDG, "PDG", "apl");
  leg3->Draw();

  c2->cd(2);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  //gPad->SetLogy();
  //
  for(Int_t i = 0;i<gr_flux_p_PDG->GetN();i++){
    gr_flux_p_PDG->GetPoint(i,ee,ff_d);
    ff = f_cpf->Eval(ee);
    gr_diff_DAMPE_fit_p_PDG->SetPoint(i,ee,TMath::Abs((ff - ff_d)/ff_d));
  }
  //
  gr_diff_DAMPE_fit_p_PDG->SetMarkerStyle(1);
  gr_diff_DAMPE_fit_p_PDG->SetLineWidth(2);
  gr_diff_DAMPE_fit_p_PDG->SetMarkerColor(kBlack);
  gr_diff_DAMPE_fit_p_PDG->SetTitle("");
  gr_diff_DAMPE_fit_p_PDG->Draw();
  gr_diff_DAMPE_fit_p_PDG->GetYaxis()->SetTitle("Relative difference");
  gr_diff_DAMPE_fit_p_PDG->GetXaxis()->SetTitle("E, GeV");
  
  //
  //
  //  
  //for(i = 0;i<nChannels;i++){
  //gr_Arr[i]->SetLineColor(colorArr[i]);
  //gr_Arr[i]->SetLineWidth(3.0);
  //gr_Arr[i]->SetMarkerColor(colorArr[i]);
  //gr_Arr[i]->SetMarkerStyle(markerArr[i]);
  //mg->Add(gr_Arr[i]);
  //}
  //  
  //h1_1->SetLineColor(kBlack);
  //h1_1->SetLineWidth(3.0);

  //h1_1->SetTitle("");
  //h1_1->Draw();
  //h2_1->SetTitle("");
  //h2_1->Draw("ZCOLOR");
  //h2_1->GetXaxis()->SetTitle("Theta, deg");
  //h2_1->GetYaxis()->SetTitle("Distance to the surface, km");

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

Double_t get_intensity_of_primary_nucleons( Double_t e, Int_t N){
  return 1.8*10000.0/TMath::Power(e/N,2.7);
}

//E     Emin  Emax     F    s_stat s_ana s_had (GeV-1_m-2_s-1_sr-1)
//49.8  39.8  63.1     2.97 0.00 0.14 0.20  1.0e-1
void get_DAMPE(TGraphErrors *gr, TString name_dat_file){
  //
  Double_t E;
  Double_t Emin;
  Double_t Emax;
  Double_t F;
  Double_t s_stat;
  Double_t s_ana;
  Double_t s_had;
  Double_t frac;
  Double_t Flux;
  Double_t Flux_err;
  //
  Int_t np = 0;
  //  
  ifstream indata;
  string mot;
  indata.open(name_dat_file.Data());
  assert(indata.is_open());  
  for(Int_t i = 0;i<8;i++)
    indata  >> mot;
  while(indata>>E>>Emin>>Emax>>F>>s_stat>>s_ana>>s_had>>frac){
    //cout<<E<<endl;
    Flux = F*frac;
    Flux_err = TMath::Sqrt(s_stat*s_stat*frac*frac + s_ana*s_ana*frac*frac + s_had*s_had*frac*frac);    
    gr->SetPoint(np,E,Flux);
    gr->SetPointError(np,(Emax - Emin)/2.0,Flux_err);
    np++;
    //cout<<mot<<endl;
  }
  indata.close();
}

Double_t getFlux(TF1 *f_cpf, Double_t el, Double_t er){
  //
  Int_t n = 10000000;
  Double_t dx = (er - el)/n;
  Double_t x;
  Double_t integral_flux = 0.0;
  for(Int_t i = 0;i<n;i++){
    x = el + dx/2.0 + dx*i;
    integral_flux = integral_flux + dx*f_cpf->Eval(x);
    //cout<<f_cpf->Eval(x)<<endl;
  }
  //cout<<"integral_flux = "<<integral_flux<<endl;
  return integral_flux;
}

void get_flux_p_PDG(TGraph *gr, TString name_dat_file){
  //
  Double_t Ekin;
  Double_t F;
  //  
  ifstream indata;
  indata.open(name_dat_file.Data());
  assert(indata.is_open());  
  while(indata>>Ekin>>F){
    //cout<<E<<endl;
    gr->SetPoint(gr->GetN(),Ekin+1,F);
  }
  indata.close();
}
