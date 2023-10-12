//my
#include "cpv.hh"

//root
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TVector3.h>
#include <TGraph2D.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <bits/stdc++.h>

using namespace std;

void cpv::Loop(TString histOut){
  TGraph2D *gr2D_int1 = new TGraph2D();
  gr2D_int1->SetNameTitle("gr2D_int1", "gr2D_int1");
  TGraph2D *gr2D_sphere = new TGraph2D();
  gr2D_sphere->SetNameTitle("gr2D_sphere","gr2D_sphere");
  TGraph2D *gr2D_gen = new TGraph2D();
  gr2D_gen->SetNameTitle("gr2D_gen","gr2D_gen");
  //
  TH1D *h1_theta = new TH1D("h1_theta","h1_theta",1000,0.0,2*TMath::Pi());
  TH1D *h1_phi = new TH1D("h1_phi","h1_phi",1000,0.0,2*TMath::Pi());
  TH1D *h1_theta_deg = new TH1D("h1_theta_deg","h1_theta_deg",10000,0.0,181.0);
  TH1D *h1_phi_deg = new TH1D("h1_phi_deg","h1_phi_deg",1000,0.0,360.0);
  //
  Double_t theta_deg;
  Double_t phi_deg;
  //
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if(jentry%100000 == 0)
      cout<<jentry<<endl;
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //
    gr2D_int1->SetPoint(gr2D_int1->GetN(),x1_int,y1_int,z1_int);
    //
    gr2D_gen->SetPoint(gr2D_gen->GetN(),x0,y0,z0);
    //
    theta_deg = theta*180.0/TMath::Pi();
    phi_deg = phi*180.0/TMath::Pi();
    //
    h1_theta->Fill(theta);
    h1_phi->Fill(phi);
    //
    h1_theta_deg->Fill(theta_deg);
    h1_phi_deg->Fill(phi_deg);
  }
  //
  Double_t thetaE;
  Double_t costhetaE;
  Double_t phiE;
  TVector3 v1;
  TRandom3 *rnd = new TRandom3(21341);
  for(Int_t i = 0;i<1000;i++){
    costhetaE = rnd->Uniform(-1.0,1.0);
    thetaE = TMath::ACos(costhetaE);
    phiE = rnd->Uniform(0.0,TMath::TwoPi());
    v1.SetMagThetaPhi(6371,thetaE,phiE);
    gr2D_sphere->SetPoint(gr2D_sphere->GetN(),v1.x(),v1.y(),v1.z());
  }

  
  /*
  //
  TGraph2D *gr2D_sphere = new TGraph2D();
  gr2D_sphere->SetTitle("gr2D_sphere");
  gr2D_sphere->SetName("gr2D_sphere");
  //
  TGraph2D *gr2D_pint = new TGraph2D();
  gr2D_pint->SetTitle("gr2D_pint");
  gr2D_pint->SetName("gr2D_pint");
  //
  TGraph2D *gr2D_p0 = new TGraph2D();
  gr2D_p0->SetTitle("gr2D_p0");
  gr2D_p0->SetName("gr2D_p0");
  //
  TGraph2D *gr2D_pe0 = new TGraph2D();
  gr2D_pe0->SetTitle("gr2D_pe0");
  gr2D_pe0->SetName("gr2D_pe0");
  //
  //
  TH1D *h1_theta_deg = new TH1D("h1_theta_deg","h1_theta_deg",1000,0.0,180.0);
  TH1D *h1_phi_deg = new TH1D("h1_phi_deg","h1_phi_deg",1000,-10.0,370);
  TH1D *h1_distToEarth = new TH1D("h1_distToEarth","h1_distToEarth",1000,0.0,65.0);
  TH1D *h1_distToTerzina = new TH1D("h1_distToTerzina","h1_distToTerzina",1000,0.0,45.0);
  TH1D *h1_angleTrzinaTrk_deg = new TH1D("h1_angleTrzinaTrk_deg","h1_angleTrzinaTrk_deg",1000,0.0,190.0);
  //
  TH1D *h1_heightAboveEarth = new TH1D("h1_heightAboveEarth","heightAboveEarth",1000,0.0,1000);
  //
  TH1D *h1_nphotons_per_m2 = new TH1D("h1_nphotons_per_m2","nphotons_per_m2",1000,0.0,1000);
  //
  TH2D *h2_distToEarth_vs_theta = new TH2D("h2_distToEarth_vs_theta","distToEarth vs theta",50,66.5,70,50,0.0,60.0);
  //
  TH2D *h2_distToEarth_vs_theta_norm = new TH2D("h2_distToEarth_vs_theta_norm","distToEarth vs theta norm",50,66.5,70,50,0.0,60.0);
  //
  TH2D *h2_distToEarth_vs_theta_00 = new TH2D("h2_distToEarth_vs_theta_00","distToEarth vs theta 00 phi<260",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_01 = new TH2D("h2_distToEarth_vs_theta_01","distToEarth vs theta 01 260<=phi<262",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_02 = new TH2D("h2_distToEarth_vs_theta_02","distToEarth vs theta 02 262<=phi<264",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_03 = new TH2D("h2_distToEarth_vs_theta_03","distToEarth vs theta 03 264<=phi<266",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_04 = new TH2D("h2_distToEarth_vs_theta_04","distToEarth vs theta 04 266<=phi<268",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_05 = new TH2D("h2_distToEarth_vs_theta_05","distToEarth vs theta 05 268<=phi<270",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_06 = new TH2D("h2_distToEarth_vs_theta_06","distToEarth vs theta 06 270<=phi<272",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_07 = new TH2D("h2_distToEarth_vs_theta_07","distToEarth vs theta 07 272<=phi<274",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_08 = new TH2D("h2_distToEarth_vs_theta_08","distToEarth vs theta 08 274<=phi<276",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_09 = new TH2D("h2_distToEarth_vs_theta_09","distToEarth vs theta 09 276<=phi<278",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_10 = new TH2D("h2_distToEarth_vs_theta_10","distToEarth vs theta 10 278<=phi<=280",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_11 = new TH2D("h2_distToEarth_vs_theta_11","distToEarth vs theta 11 phi>280",50,66.5,70,50,0.0,60.0);

  TH2D *h2_distToEarth_vs_theta_norm_00 = new TH2D("h2_distToEarth_vs_theta_norm_00","distToEarth vs theta norm 00 phi<260",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_norm_01 = new TH2D("h2_distToEarth_vs_theta_norm_01","distToEarth vs theta norm 01 260<=phi<262",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_norm_02 = new TH2D("h2_distToEarth_vs_theta_norm_02","distToEarth vs theta norm 02 262<=phi<264",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_norm_03 = new TH2D("h2_distToEarth_vs_theta_norm_03","distToEarth vs theta norm 03 264<=phi<266",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_norm_04 = new TH2D("h2_distToEarth_vs_theta_norm_04","distToEarth vs theta norm 04 266<=phi<268",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_norm_05 = new TH2D("h2_distToEarth_vs_theta_norm_05","distToEarth vs theta norm 05 268<=phi<270",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_norm_06 = new TH2D("h2_distToEarth_vs_theta_norm_06","distToEarth vs theta norm 06 270<=phi<272",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_norm_07 = new TH2D("h2_distToEarth_vs_theta_norm_07","distToEarth vs theta norm 07 272<=phi<274",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_norm_08 = new TH2D("h2_distToEarth_vs_theta_norm_08","distToEarth vs theta norm 08 274<=phi<276",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_norm_09 = new TH2D("h2_distToEarth_vs_theta_norm_09","distToEarth vs theta norm 09 276<=phi<278",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_norm_10 = new TH2D("h2_distToEarth_vs_theta_norm_10","distToEarth vs theta norm 10 278<=phi<=280",50,66.5,70,50,0.0,60.0);
  TH2D *h2_distToEarth_vs_theta_norm_11 = new TH2D("h2_distToEarth_vs_theta_norm_11","distToEarth vs theta norm 11 phi>280",50,66.5,70,50,0.0,60.0);

  //
  //
  Double_t phiDeg;
  //
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  Double_t heightAboveEarth = 0.0;//km
  const Double_t earthR = 6371;   // km
  Int_t eventIDmy = 0;
  double nphotons_per_m2 = 0.0;
  TString inShowerRootFileName; 
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  //for (Long64_t jentry=0; jentry<100;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if(jentry%1000 == 0)
      cout<<jentry<<endl;
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(angleTrzinaTrk>170.0/180.0*TMath::Pi()){
      gr2D_pint->SetPoint(gr2D_pint->GetN(),x_int,y_int,z_int);
      gr2D_p0->SetPoint(gr2D_p0->GetN(),x0,y0,z0);
      gr2D_pe0->SetPoint(gr2D_pe0->GetN(),xe0,ye0,ze0);
      heightAboveEarth = TMath::Sqrt(x_int*x_int + y_int*y_int + z_int*z_int) - earthR;
      //
      h1_theta_deg->Fill(theta*180.0/TMath::Pi());
      h1_phi_deg->Fill(phi*180.0/TMath::Pi());
      phiDeg=phi*180.0/TMath::Pi();
      h1_distToEarth->Fill(distToEarth);
      h1_distToTerzina->Fill(distToTerzina);
      h1_angleTrzinaTrk_deg->Fill(angleTrzinaTrk*180.0/TMath::Pi());
      h1_heightAboveEarth->Fill(heightAboveEarth);
      // 
      //old function for ini creation
      //createEASCherSim_ini(eventIDmy,theta*180.0/TMath::Pi(),heightAboveEarth);
      //create_cosmique_proton_generator_info(eventIDmy,theta*180.0/TMath::Pi(),phi*180.0/TMath::Pi(),heightAboveEarth);
      inShowerRootFileName = getShowerRootFileName(eventIDmy);
      //readEASCherSim( inShowerRootFileName, distToTerzina, nphotons_per_m2);
      readEASCherSimNewFormat( inShowerRootFileName, distToTerzina, nphotons_per_m2);
      printTrkInfo(eventIDmy,nphotons_per_m2);
      //cout<<"nphotons_per_m2 = "<<nphotons_per_m2<<endl;
      //
      h1_nphotons_per_m2->Fill(nphotons_per_m2);
      //
      h2_distToEarth_vs_theta->Fill(theta*180.0/TMath::Pi(),distToEarth);
      h2_distToEarth_vs_theta_norm->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      //
      if(phi*180.0/TMath::Pi()<260.0){
	h2_distToEarth_vs_theta_00->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_00->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      if(phi*180.0/TMath::Pi() >= 260.0 && phi*180.0/TMath::Pi() < 262.0){
	h2_distToEarth_vs_theta_01->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_01->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      if(phi*180.0/TMath::Pi() >= 262.0 && phi*180.0/TMath::Pi() < 264.0){
	h2_distToEarth_vs_theta_02->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_02->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      if(phi*180.0/TMath::Pi() >= 264.0 && phi*180.0/TMath::Pi() < 266.0){
	h2_distToEarth_vs_theta_03->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_03->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      if(phi*180.0/TMath::Pi() >= 266.0 && phi*180.0/TMath::Pi() < 268.0){
	h2_distToEarth_vs_theta_04->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_04->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      if(phi*180.0/TMath::Pi() >= 268.0 && phi*180.0/TMath::Pi() < 270.0){
	h2_distToEarth_vs_theta_05->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_05->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      if(phi*180.0/TMath::Pi() >= 270.0 && phi*180.0/TMath::Pi() < 272.0){
	h2_distToEarth_vs_theta_06->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_06->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      if(phi*180.0/TMath::Pi() >= 272.0 && phi*180.0/TMath::Pi() < 274.0){
	h2_distToEarth_vs_theta_07->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_07->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      if(phi*180.0/TMath::Pi() >= 274.0 && phi*180.0/TMath::Pi() < 276.0){
	h2_distToEarth_vs_theta_08->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_08->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      if(phi*180.0/TMath::Pi() >= 276.0 && phi*180.0/TMath::Pi() < 278.0){
	h2_distToEarth_vs_theta_09->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_09->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      if(phi*180.0/TMath::Pi() >= 278.0 && phi*180.0/TMath::Pi() <= 280.0){
	h2_distToEarth_vs_theta_10->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_10->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      if(phi*180.0/TMath::Pi() > 280.0){
	h2_distToEarth_vs_theta_11->Fill(theta*180.0/TMath::Pi(),distToEarth);
	h2_distToEarth_vs_theta_norm_11->Fill(theta*180.0/TMath::Pi(),distToEarth,nphotons_per_m2);
      }
      //
      eventIDmy++;
    }
  }
  //
  h2D2div(h2_distToEarth_vs_theta_norm,h2_distToEarth_vs_theta);
  h2D2div(h2_distToEarth_vs_theta_norm_00,h2_distToEarth_vs_theta_00);
  h2D2div(h2_distToEarth_vs_theta_norm_01,h2_distToEarth_vs_theta_01);
  h2D2div(h2_distToEarth_vs_theta_norm_02,h2_distToEarth_vs_theta_02);
  h2D2div(h2_distToEarth_vs_theta_norm_03,h2_distToEarth_vs_theta_03);
  h2D2div(h2_distToEarth_vs_theta_norm_04,h2_distToEarth_vs_theta_04);
  h2D2div(h2_distToEarth_vs_theta_norm_05,h2_distToEarth_vs_theta_05);
  h2D2div(h2_distToEarth_vs_theta_norm_06,h2_distToEarth_vs_theta_06);
  h2D2div(h2_distToEarth_vs_theta_norm_07,h2_distToEarth_vs_theta_07);
  h2D2div(h2_distToEarth_vs_theta_norm_08,h2_distToEarth_vs_theta_08);
  h2D2div(h2_distToEarth_vs_theta_norm_09,h2_distToEarth_vs_theta_09);
  h2D2div(h2_distToEarth_vs_theta_norm_10,h2_distToEarth_vs_theta_10);
  h2D2div(h2_distToEarth_vs_theta_norm_11,h2_distToEarth_vs_theta_11);
  //
  */
  TFile* rootFile = new TFile(histOut.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<histOut.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<histOut.Data()<<endl;
  //
  gr2D_int1->Write();
  gr2D_sphere->Write();
  gr2D_gen->Write();
  //
  h1_theta->Write();
  h1_phi->Write();
  //
  h1_theta_deg->Write();
  h1_phi_deg->Write();
  //
  //gr2D_sphere3->Write();
  /*
  gr2D_sphere->Write();
  gr2D_pint->Write();
  gr2D_p0->Write();
  gr2D_pe0->Write();
  //
  h1_theta_deg->Write();
  h1_phi_deg->Write();
  h1_distToEarth->Write();
  h1_distToTerzina->Write();
  h1_angleTrzinaTrk_deg->Write();
  //
  h1_heightAboveEarth->Write();
  //
  h1_nphotons_per_m2->Write();
  //
  h2_distToEarth_vs_theta->Write();
  h2_distToEarth_vs_theta_norm->Write();
  //
  h2_distToEarth_vs_theta_00->Write();
  h2_distToEarth_vs_theta_01->Write();
  h2_distToEarth_vs_theta_02->Write();
  h2_distToEarth_vs_theta_03->Write();
  h2_distToEarth_vs_theta_04->Write();
  h2_distToEarth_vs_theta_05->Write();
  h2_distToEarth_vs_theta_06->Write();
  h2_distToEarth_vs_theta_07->Write();
  h2_distToEarth_vs_theta_08->Write();
  h2_distToEarth_vs_theta_09->Write();
  h2_distToEarth_vs_theta_10->Write();
  h2_distToEarth_vs_theta_11->Write();
  //
  h2_distToEarth_vs_theta_norm_00->Write();
  h2_distToEarth_vs_theta_norm_01->Write();
  h2_distToEarth_vs_theta_norm_02->Write();
  h2_distToEarth_vs_theta_norm_03->Write();
  h2_distToEarth_vs_theta_norm_04->Write();
  h2_distToEarth_vs_theta_norm_05->Write();
  h2_distToEarth_vs_theta_norm_06->Write();
  h2_distToEarth_vs_theta_norm_07->Write();
  h2_distToEarth_vs_theta_norm_08->Write();
  h2_distToEarth_vs_theta_norm_09->Write();
  h2_distToEarth_vs_theta_norm_10->Write();
  h2_distToEarth_vs_theta_norm_11->Write();
  //
  */
  rootFile->Close();
}

/*
void cpv::createEASCherSim_ini( Int_t eventIDmy, Double_t theta_deg, Double_t heightAboveEarth_km){
  char buffer [8];
  sprintf(buffer ,"%07d", eventIDmy);
  //TString fileName = "../EASCherSim/ini/EASCherSim_";
  TString fileName;
  //
  TString iniDir = "../easchersim/ini/";
  Int_t nFilesPerDir = 200;
  Int_t dirNameID = eventIDmy/nFilesPerDir;
  char buffer2 [6];
  sprintf(buffer2 ,"%05d", dirNameID);
  TString fileNameIDstr =  buffer2;
  iniDir += fileNameIDstr;
  iniDir += "/";
  //
  TString fileNameID =  buffer;
  fileName = iniDir;
  fileName += "EASCherSim_";

  fileName += fileNameID;
  fileName += ".ini";
  cout<<fileName<<endl;
  //
  ofstream myfile;
  myfile.open (fileName.Data());
  myfile<<"[shower]"<<endl;
  myfile<<"energy = 1e+17"<<endl;
  myfile<<"angle = "<<theta_deg<<", viewing"<<endl;
  myfile<<"distFirstInt = 0.0, distance"<<endl;
  myfile<<endl;
  myfile<<"[detector]"<<endl;
  myfile<<"altitude = "<<heightAboveEarth_km<<endl;
  myfile<<endl;
  myfile<<"[atmosphere]"<<endl;
  myfile<<"scattering = yes"<<endl;
  myfile<<endl;
  myfile<<"[magField]"<<endl;
  myfile<<"useField = off"<<endl;
  myfile<<endl;
  myfile<<"[runSettings]"<<endl;
  myfile<<"useGreisenParametrization = no"<<endl;
  myfile<<endl;
  myfile<<"[output]"<<endl;
  myfile<<"plots = none"<<endl;
  myfile.close();
}

TString cpv::getShowerRootFileName(Int_t eventIDmy){
  // 
  //../easchersim/root/00009/EASCherSim_0001894.ini.npz.root
  //
  char buffer [8];
  sprintf(buffer ,"%07d", eventIDmy);
  //TString fileName = "../EASCherSim/ini/EASCherSim_";
  TString rootDir = "../easchersim/root/";
  Int_t nFilesPerDir = 200;
  Int_t dirNameID = eventIDmy/nFilesPerDir;
  char buffer2 [6];
  sprintf(buffer2 ,"%05d", dirNameID);
  TString fileNameIDstr =  buffer2;
  rootDir += fileNameIDstr;
  rootDir += "/";
  //
  TString fileNameID =  buffer;
  //
  TString fileName = rootDir;
  fileName += "EASCherSim_";
  //
  fileName += fileNameID;
  fileName += ".ini.npz.root";
  //cout<<fileName<<endl;
  //
  //TString fileName = "/home/dpncguest/home2/work/POEMMA/geant4/cosmique_proton_generator/cosmique_proton_generator.info";
  //ofstream myfile;
  //char buffer[6];
  //sprintf(buffer ,"%05d", eventIDmy);
  //TString fileNameID =  buffer;
  ////
  //TString fileNameini = "EASCherSim_";
  //fileNameini += fileNameID;
  //fileNameini += ".ini";
  ////
  //TString fileNameroot = "/home/dpncguest/home2/work/POEMMA/geant4/EASCherSim/root/";
  //fileNameroot += fileNameini;
  //fileNameroot += ".root";
  return fileName;
}

void cpv::printTrkInfo( Int_t eventIDmy, double nphotons_per_m2){
  //theta       phi 
  //1.17706 4.86938 
  //         x_int        y_int         z_int            xe0           ye0          ze0 
  //-8.87926639948 1.0402504604 6898.51118961 -390.672021698 2412.97338838 5884.0611177 
  //   distToEarth distToTerzina angleTrzinaTrk nphotons_per_m2
  //0.596380690828 9.28598761906 2.99660248562  0.0128641130477
  // 
  //../easchersim/root/00009/EASCherSim_0001894.ini.npz.root
  //
  char buffer [8];
  sprintf(buffer ,"%07d", eventIDmy);
  //TString fileName = "../EASCherSim/ini/EASCherSim_";
  TString rootDir = "../easchersim/root/";
  Int_t nFilesPerDir = 200;
  Int_t dirNameID = eventIDmy/nFilesPerDir;
  char buffer2 [6];
  sprintf(buffer2 ,"%05d", dirNameID);
  TString fileNameIDstr =  buffer2;
  rootDir += fileNameIDstr;
  rootDir += "/";
  //
  TString fileNameID =  buffer;
  //
  TString fileName = rootDir;
  fileName += "trkInfo_";
  //
  fileName += fileNameID;
  fileName += ".dat";
  //
  if(file_exists_test(getShowerRootFileName(eventIDmy))){
    //FILE *fp01;
    //fp01 = fopen(fileName.Data(), "w");
    //fprintf(fp01, " theta        m \n");
    //fprintf(fp01, "%10.7f %10.7f \n",x,y);
    //fclose(fp01);
    //
    ofstream myfile;
    myfile.open (fileName.Data());
    myfile<<"theta"<<" "
    	  <<"phi"<<" "
    	  <<"x_int"<<" "
    	  <<"y_int"<<" "
    	  <<"z_int"<<" ";
    myfile<<"xe0"<<" "
    	  <<"ye0"<<" "
    	  <<"ze0"<<" "
    	  <<"distToEarth"<<" "
    	  <<"distToTerzina"<<" "
	  <<"angleTrzinaTrk"<<" "
    	  <<"nphotons_per_m2"<<endl;
    myfile<<theta<<" "
    	  <<phi<<" "
    	  <<std::setprecision(12)<<x_int<<" "
    	  <<std::setprecision(12)<<y_int<<" "
    	  <<std::setprecision(12)<<z_int<<" ";
    myfile<<std::setprecision(12)<<xe0<<" "
    	  <<std::setprecision(12)<<ye0<<" "
    	  <<std::setprecision(12)<<ze0<<" "
    	  <<distToEarth<<" "
    	  <<distToTerzina<<" "
	  <<angleTrzinaTrk<<" "
    	  <<nphotons_per_m2<<endl;
    myfile.close();
  }
}

void cpv::create_cosmique_proton_generator_info(Int_t eventIDmy, Double_t theta_deg, Double_t phi_deg, Double_t heightAboveEarth_km){
  //
  // 00000 EASCherSim_00000.ini EASCherSim_00000.ini.root     15.1845       27.3905          519.776  67.9146  280.445 || ang = 67.9146 alt = 519.776
  // 00651 EASCherSim_00651.ini EASCherSim_00651.ini.root     5.90176       21.8164          533.125  67.4229  274.476 || ang = 67.4229 alt = 533.125
  // 01708 EASCherSim_01708.ini EASCherSim_01708.ini.root     9.08457       32.0732          500.556  68.2915  273.212 || ang = 68.2915 alt = 500.556
  //
  TString fileName = "/home/dpncguest/home2/work/POEMMA/geant4/cosmique_proton_generator/cosmique_proton_generator.info";
  ofstream myfile;
  char buffer[6];
  sprintf(buffer ,"%05d", eventIDmy);
  TString fileNameID =  buffer;
  //
  TString fileNameini = "EASCherSim_";
  fileNameini += fileNameID;
  fileNameini += ".ini";
  //
  TString fileNameroot = fileNameini;
  fileNameroot += ".root";
  //
  if(eventIDmy == 0){
    myfile.open (fileName.Data());
    myfile<<"iniFiles: /home/dpncguest/home2/work/POEMMA/geant4/EASCherSim/ini/"<<endl
	  <<"rootFiles: /home/dpncguest/home2/work/POEMMA/geant4/EASCherSim/root/"<<endl;
    myfile<<setw(10)<<"EVNTID"
	  <<setw(21)<<"iniFiles"
	  <<setw(26)<<"rootFiles"
	  <<setw(12)<<"distToEarth"
	  <<setw(14)<<"distToTerzina"
	  <<setw(17)<<"heightAboveEarth"
	  <<setw(9)<<"trkTheta"
	  <<setw(9)<<"trkPhi"
	  <<endl;
    myfile<<setw(10)<<fileNameID
	  <<setw(21)<<fileNameini
      	  <<setw(26)<<fileNameroot
	  <<setw(12)<<distToEarth
	  <<setw(14)<<distToTerzina
	  <<setw(17)<<heightAboveEarth_km
	  <<setw(9)<<theta_deg
	  <<setw(9)<<phi_deg
	  <<endl;
  }
  else{
    //
    myfile.open (fileName.Data(),std::ios_base::app);
    myfile<<setw(10)<<fileNameID
	  <<setw(21)<<fileNameini
      	  <<setw(26)<<fileNameroot
	  <<setw(12)<<distToEarth
	  <<setw(14)<<distToTerzina
	  <<setw(17)<<heightAboveEarth_km
	  <<setw(9)<<theta_deg
	  <<setw(9)<<phi_deg
	  <<endl;
  }
  myfile.close();
}
*/

/*
void cpv::readEASCherSim( TString inRootFileName, double dist,  double &nphotons_per_m2){
  TString blackList[7];
  blackList[0] = "/home/dpncguest/home2/work/POEMMA/geant4/EASCherSim/root/EASCherSim_00045.ini.root";
  blackList[1] = "/home/dpncguest/home2/work/POEMMA/geant4/EASCherSim/root/EASCherSim_00220.ini.root";
  blackList[2] = "/home/dpncguest/home2/work/POEMMA/geant4/EASCherSim/root/EASCherSim_00598.ini.root";
  blackList[3] = "/home/dpncguest/home2/work/POEMMA/geant4/EASCherSim/root/EASCherSim_01050.ini.root";
  blackList[4] = "/home/dpncguest/home2/work/POEMMA/geant4/EASCherSim/root/EASCherSim_01118.ini.root";
  blackList[5] = "/home/dpncguest/home2/work/POEMMA/geant4/EASCherSim/root/EASCherSim_01127.ini.root";
  blackList[6] = "/home/dpncguest/home2/work/POEMMA/geant4/EASCherSim/root/EASCherSim_01675.ini.root";
  if(inRootFileName != blackList[0] &&
     inRootFileName != blackList[1] &&
     inRootFileName != blackList[2] &&
     inRootFileName != blackList[3] &&
     inRootFileName != blackList[4] &&
     inRootFileName != blackList[5] &&
     inRootFileName != blackList[6]){
  TFile *fIn = new TFile(inRootFileName.Data());
  fIn->cd();
  //fIn->ls();
  //fIn->Print();

  //cout<<"ddd"<<endl;
  
  // -----------------------------------------------------------------------------
  // Shower parameter
  // -----------------------------------------------------------------------------
  //cout << endl;
  //cout << endl;
  //cout<<"ddd"<<endl;
  TTree *showerProp = (TTree*)fIn->Get("showerProp");
  //showerProp->ls();
  //showerProp->Print();
  //cout<<"eeee"<<endl;
  //
  Double_t zenith, energy, startalt;
  TBranch *energyBr = showerProp->GetBranch("energy");
  energyBr->SetAddress(&energy);
  //cout<<"ttt"<<endl;
  //
  TBranch *zenithBr = showerProp->GetBranch("zenith");
  zenithBr->SetAddress(&zenith);
  //cout<<"eeeeee"<<endl;
  //TBranch *startaltBr = showerProp->GetBranch("startalt");
  //startaltBr->SetAddress(&startalt);
  //cout<<"rreeeeee"<<endl;
  Int_t nEntries = showerProp->GetEntries();
  //cout<<"nEntries:"<<" "<<nEntries<<endl;
  showerProp->GetEntry(0);
  //cout<< "Zenith angle = " << zenith <<endl;
  //cout<< "E = " << energy << " PeV" << endl;
  //cout<< "Starting altitude = " << startalt << " km" << endl;
  // -----------------------------------------------------------------------------
  
  // -----------------------------------------------------------------------------
  // Cherenkov light parameter
  // -----------------------------------------------------------------------------
  //cout << endl;
  //cout << endl;
  TTree *cherPhProp = (TTree*)fIn->Get("cherPhProp");
  //cherPhProp->ls();
  //cherPhProp->Print();
  TH1D *wavelength = NULL;
  TH1D *distance = NULL;
  TH1D *time_offset = NULL;
  TH1D *angle_offset = NULL;
  // wavelength
  TBranch *wavelengthBr = cherPhProp->GetBranch("wavelength");
  wavelengthBr->SetAddress(&wavelength);
  cherPhProp->GetEntry(0);
  TCanvas *cwl = new TCanvas("cwl", "cwl", 800, 600);
  cwl->cd();
  wavelength->SetDirectory(0);
  //wavelength->GetXaxis()->SetTitle("wavelength [nm]");
  //wavelength->Draw();
  // distance
  TBranch *distanceBr = cherPhProp->GetBranch("distance");
  distanceBr->SetAddress(&distance);
  cherPhProp->GetEntry(0);
  TCanvas *cdst = new TCanvas("cdst", "cdst", 800, 600);
  cdst->cd();
  distance->SetDirectory(0);
  TH1D *hdistance = (TH1D *)(distance)->Clone();
  hdistance->SetTitle("distance");
  //hdistance->GetYaxis()->SetTitle("photon density [pht/m^{2}]");
  //hdistance->GetXaxis()->SetTitle("distance [km]");
  //hdistance->Draw();
  // time_offset
  TBranch *time_offsetBr = cherPhProp->GetBranch("time_offset");
  time_offsetBr->SetAddress(&time_offset);
  cherPhProp->GetEntry(0);
  TCanvas *ctoff = new TCanvas("ctoff", "ctoff", 800, 600);
  ctoff->cd();
  time_offset->SetDirectory(0);
  //time_offset->GetXaxis()->SetTitle("ns");
  //time_offset->Draw("histo");
  // angle_offset
  TBranch *angle_offsetBr = cherPhProp->GetBranch("angle_offset");
  angle_offsetBr->SetAddress(&angle_offset);
  cherPhProp->GetEntry(0);
  TCanvas *cangoff = new TCanvas("cangoff", "cangoff", 800, 600);
  cangoff->cd();
  angle_offset->SetDirectory(0);
  //angle_offset->GetXaxis()->SetTitle("deg");
  //angle_offset->Draw("histo");
  //----------------------------------------------------------------------------------------


  //----------------------------------------------------------------------------------------
  //double dist=40.0;
  // time distribution
  TTreeReader fReaderCherPh_time;
  TTreeReaderValue<std::vector<TH1D>> fDistTimeVec = {fReaderCherPh_time, "time_dist"};
  fReaderCherPh_time.SetTree("cherPhProp");
  fReaderCherPh_time.SetEntry(0);
  // get number of photons at aperture
  auto nBinDist = distance->FindBin(dist);

  nphotons_per_m2 = distance->GetBinContent((int)nBinDist);

  //cout<<"nphotons_per_m2 = "<<nphotons_per_m2<<endl;

  //cout<<inRootFileName<<endl;
  
  fIn->Close();
  }  
*/
  /*
  
  //cout<<"distance->GetNbinsX() = "<<distance->GetNbinsX()<<endl;
  //cout<<"distance->FindBin(7.0) = "<<distance->FindBin(7.0)<<endl;
  //cout<<"distance->FindBin(0.0) = "<<distance->FindBin(0.0)<<endl;
  //cout<<"distance->FindBin(15.0) = "<<distance->FindBin(15.0)<<endl;
  //cout<<"distance->FindBin(100.0) = "<<distance->FindBin(100.0)<<endl;
  //cout<<"distance->FindBin(200.0) = "<<distance->FindBin(200.0)<<endl;
  //cout<<"distance->GetBinCenter(1) = "<<distance->GetBinCenter(1)<<endl;
  //cout<<"distance->GetBinCenter(63) = "<<distance->GetBinCenter(63)<<endl;
  //cout<<"distance->GetBinCenter(400) = "<<distance->GetBinCenter(400)<<endl;
  //cout<<"distance->GetBinLowEdge(1) = "<<distance->GetBinLowEdge(1)<<endl;
  //cout<<"distance->GetBinWidth(1) = "<<distance->GetBinWidth(1)<<endl;

  //cout<<"nBinDist "<<nBinDist<<endl;
  //auto ApertureSize = 1.0; // in square meters
  //auto nPh = distance->GetBinContent(nBinDist) * ApertureSize;
  TH1D *timeHist = NULL;
  timeHist = (TH1D*)((*fDistTimeVec)[nBinDist-1]).Clone();
  //timeHist->GetYaxis()->SetTitle("photon density [pht/m^{2}]");
  timeHist->GetYaxis()->SetTitle("a.u.");
  timeHist->GetXaxis()->SetTitle("time [ns]");
  // Plot histogram
  TCanvas *ctmHst = new TCanvas("ctmHst", "ctmHst", 800, 600);
  ctmHst->cd();
  ctmHst->SetLogx();
  ctmHst->SetLogy();
  //timeHist->Draw();

  //----------------------------------------------------------------------------------------
  // Angle distribution
  TTreeReader fReaderCherPh_angle;
  TTreeReaderValue<std::vector<TH1D>> fDistAngleVec = {fReaderCherPh_angle, "angle_dist"};
  fReaderCherPh_angle.SetTree("cherPhProp");
  fReaderCherPh_angle.SetEntry(0);
  TH1D *angleHist = NULL;
  angleHist = (TH1D*)((*fDistAngleVec)[nBinDist-1]).Clone();
  angleHist->GetYaxis()->SetTitle("photon density [counts/m^{2}]");
  angleHist->GetXaxis()->SetTitle("angle, deg");
  //Plot histogram
  TCanvas *canglHst = new TCanvas("canglHst", "canglHst", 800, 600);
  canglHst->cd();
  canglHst->SetLogx();
  canglHst->SetLogy();
  //angleHist->Draw();

  */
  //assert(0);
//}

/*
void cpv::readEASCherSimNewFormat( TString inRootFileName, double dist,  double &nphotons_per_m2){
  if(file_exists_test(inRootFileName)){
    TFile *fIn = new TFile(inRootFileName.Data());
    fIn->cd();
    TH1D *distance = (TH1D*)fIn->Get("r");
    auto nBinDist = distance->FindBin(dist);
    nphotons_per_m2 = distance->GetBinContent((int)nBinDist);
    fIn->Close();
  }
  else {
    nphotons_per_m2 = 0.0;
  }
}

bool cpv::file_exists_test(TString inRootFileName) {
  if (FILE *file = fopen(inRootFileName.Data(), "r")) {
    fclose(file);
    return true;
  }else {
    return false;
  }  
  return false;
}
*/
