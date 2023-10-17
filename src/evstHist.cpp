//my
#include "evstHist.hh"

//c, c++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <vector>

//root
#include <TVector2.h>
#include <TPolyLine.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TText.h>
#include <TMath.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TCrown.h>
#include <TArc.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TPad.h>
#include <TString.h>
#include <TFile.h>
#include <TAxis.h>
#include <TVector2.h>
#include <TImage.h>

using namespace std;

evstHist::evstHist(): evstHist("","")
{
}

evstHist::evstHist(const char* name, const char* title,
		   Double_t val_Emin, Double_t val_Emax, Int_t val_N_bins_E,
		   Double_t val_Thetamin, Double_t val_Thetamax, Int_t val_N_bins_t) : TH2Poly()
{    
  //
  SetName(name);
  SetTitle(title);
  //
  //_Emin = 1.0;     //1   GeV
  //_Emax = 100000;  //100 TeV
  //_N_bins_E  = 25; //
  _Emin = val_Emin;
  _Emax = val_Emax;
  _N_bins_E  = val_N_bins_E;
  //
  //
  _Eminarr = new Double_t[_N_bins_E];
  _Emaxarr = new Double_t[_N_bins_E];
  Double_t *E_bins = new Double_t[_N_bins_E+1];
  //
  //
  Double_t bin_log;
  Double_t E_log_l;
  Double_t E_log_r;
  bin_log = (TMath::Log10(_Emax) - TMath::Log10(_Emin))/_N_bins_E;
  for(Int_t i = 0;i<_N_bins_E; i++){
    E_log_l = TMath::Log10(_Emin) + bin_log*i;
    E_log_r = TMath::Log10(_Emin) + bin_log*(i+1);
    _Eminarr[i] = TMath::Power(10,E_log_l);
    _Emaxarr[i] = TMath::Power(10,E_log_r); 
    E_bins[i] = _Eminarr[i];
    cout<<E_bins[i]<<endl;
  }
  E_bins[_N_bins_E] = _Emax;
  _h1_E = new TH1D("_h1_E","_h1_E", _N_bins_E, E_bins);
  //  
  //_Thetamin = 0.0;  //deg
  //_Thetamax = 10.0; //deg
  //_N_bins_t = 10;   //
  //
  _Thetamin = val_Thetamin;
  _Thetamax = val_Thetamax;
  _N_bins_t = val_N_bins_t;
  //
  _h1_theta = new TH1D("_h1_theta","_h1_theta", _N_bins_t, _Thetamin, _Thetamax);
  //
  _Thetaminarr = new Double_t[_N_bins_t];
  _Thetamaxarr = new Double_t[_N_bins_t];
  for(Int_t i = 0;i<_N_bins_t; i++){
    _Thetaminarr[i] = (_Thetamax - _Thetamin)/_N_bins_t*i;
    _Thetamaxarr[i] = (_Thetamax - _Thetamin)/_N_bins_t*(i+1);
  }
  //
  Int_t n = 5;  
  Double_t* x;
  Double_t* y;
  x = new Double_t [n];
  y = new Double_t [n];
  //
  for(Int_t i_E = 0;i_E<_N_bins_E;i_E++){
    for(Int_t i_t = 0;i_t<_N_bins_t;i_t++){
      //
      //we go clockwise
      x[0] = _Thetaminarr[i_t];
      y[0] = _Eminarr[i_E];
      //
      x[1] = x[0];
      y[1] = _Emaxarr[i_E];
      //
      x[2] = _Thetamaxarr[i_t];
      y[2] = y[1];
      //
      x[3] = x[2];
      y[3] = y[0];
      //
      x[4] = x[0];
      y[4] = y[0];
      //
      AddBin(n,x,y);
    }
  }  
}

evstHist::~evstHist(){
}

void evstHist::test(){
  for(Int_t i = 0;i<GetNcells();i++)
    SetBinContent(i+1,i+1);
};

void evstHist::Draw_hist(TString fileName){

  //kDeepSea=51,          kGreyScale=52,    kDarkBodyRadiator=53,
  //kBlueYellow= 54,      kRainBow=55,      kInvertedDarkBodyRadiator=56,
  //kBird=57,             kCubehelix=58,    kGreenRedViolet=59,
  //kBlueRedYellow=60,    kOcean=61,        kColorPrintableOnGrey=62,
  //kAlpine=63,           kAquamarine=64,   kArmy=65,
  //kAtlantic=66,         kAurora=67,       kAvocado=68,
  //kBeach=69,            kBlackBody=70,    kBlueGreenYellow=71,
  //kBrownCyan=72,        kCMYK=73,         kCandy=74,
  //kCherry=75,           kCoffee=76,       kDarkRainBow=77,
  //kDarkTerrain=78,      kFall=79,         kFruitPunch=80,
  //kFuchsia=81,          kGreyYellow=82,   kGreenBrownTerrain=83,
  //kGreenPink=84,        kIsland=85,       kLake=86,
  //kLightTemperature=87, kLightTerrain=88, kMint=89,
  //kNeon=90,             kPastel=91,       kPearl=92,
  //kPigeon=93,           kPlum=94,         kRedBlue=95,
  //kRose=96,             kRust=97,         kSandyTerrain=98,
  //kSienna=99,           kSolar=100,       kSouthWest=101,
  //kStarryNight=102,     kSunset=103,      kTemperatureMap=104,
  //kThermometer=105,     kValentine=106,   kVisibleSpectrum=107,
  //kWaterMelon=108,      kCool=109,        kCopper=110,
  //kGistEarth=111,       kViridis=112,     kCividis=113

  gStyle->SetPalette(kInvertedDarkBodyRadiator);

  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->SetRightMargin(0.12);
  c1->SetLeftMargin(0.12);
  c1->SetTopMargin(0.1);
  c1->SetBottomMargin(0.15);
  c1->SetLogy();
  c1->SetLogz();
  
  TH2F *frame = new TH2F("h2","h2", 40, _Thetamin, _Thetamax, 40, _Emin, _Emax);
  frame->SetTitle("");
  //frame->GetXaxis()->SetTitle("x (mm)");
  //frame->GetYaxis()->SetTitle("y (mm)");
  //frame->GetXaxis()->CenterTitle();
  //frame->GetYaxis()->CenterTitle();
  //frame->GetYaxis()->SetTitleOffset(1.5);
  frame->SetStats(kFALSE);
  frame->Draw();

  //SetMaximum(50);
  //SetMinimum(0);

  SetMarkerSize(0.6);
  
  //Draw("same TEXT");
  Draw("same ZCOLOR TEXT");
  //Draw("same ZCOLOR");
  //c1->SaveAs("Draw_map_test.pdf");
  c1->SaveAs(fileName.Data());
}


