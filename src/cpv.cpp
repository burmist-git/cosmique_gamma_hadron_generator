//my
#include "cpv.hh"
#include "evstHist.hh"

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

void cpv::Loop(TString histOut,
	       Double_t y0_LST01_km, Double_t x0_LST01_km,
	       Int_t max_trg_ev, Double_t norm_theta_p_t_deg_max, Double_t norm_r_core_max){
  //
  Double_t theta_deg;
  Double_t phi_deg;
  Double_t theta_p_t;
  Double_t theta_p_t_deg;
  Double_t azimuth_deg;
  Double_t altitude_deg;
  Double_t r_from_tel;
  //
  TH1D *h1_theta_deg = new TH1D("h1_theta_deg","h1_theta_deg",10000,0.0,181.0);
  TH1D *h1_phi_deg   = new TH1D("h1_phi_deg","h1_phi_deg",10000,-360.0,360.0);
  TH1D *h1_theta_p_t_deg = new TH1D("h1_theta_p_t_deg","h1_theta_p_t_deg",1000,0.0,181);
  //
  TH1D *h1_azimuth_deg = new TH1D("h1_azimuth_deg","h1_azimuth_deg",400,0.0,360);
  TH1D *h1_altitude_deg = new TH1D("h1_altitude_deg","h1_altitude_deg",400,0.0,180);
  //
  TH1D *h1_x1_int = new TH1D("h1_x1_int","h1_x1_int",200,-1.8,1.8);
  TH1D *h1_y1_int = new TH1D("h1_y1_int","h1_y1_int",200,-1.8,1.8);
  //
  TH1D *h1_r_from_tel = new TH1D("h1_r_from_tel","h1_r_from_tel",1000,0.0,2.0);
  //
  Int_t n_trg_ev = 0;
  //
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  TVector3 v_det(1.0*TMath::Sin(20.0/180.0*TMath::Pi()),0.0,1.0*TMath::Cos(20.0/180.0*TMath::Pi()));
  ////////
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = LoadTree(jentry);
    if(jentry%100000000 == 0)
      cout<<jentry<<endl;
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //
    theta_deg = theta*180.0/TMath::Pi();
    phi_deg = conv_phi(phi)*180.0/TMath::Pi();
    //
    TVector3 v_prot(-vx,-vy,-vz);
    theta_p_t = TMath::ACos(v_prot.Dot(v_det)/v_prot.Mag()/v_det.Mag());
    theta_p_t_deg = theta_p_t*180/TMath::Pi();
    if(theta_p_t_deg<=norm_theta_p_t_deg_max){
      r_from_tel = TMath::Sqrt((x0_LST01_km - x1_int)*(x0_LST01_km - x1_int) + (y0_LST01_km - y1_int)*(y0_LST01_km - y1_int));
      if(r_from_tel<=norm_r_core_max){
	TVector3 v_prot_azimuth_alt;
	v_prot_azimuth_alt.SetMagThetaPhi(1.0,theta,conv_phi(phi));
	TVector3 v_prot_azimuth_alt_inv(-v_prot_azimuth_alt.x(),-v_prot_azimuth_alt.y(),-v_prot_azimuth_alt.z());
	azimuth_deg = conv_phi(phi)*180/TMath::Pi();
	altitude_deg = 90.0 - v_prot_azimuth_alt_inv.Theta()*180/TMath::Pi();
	//
	h1_theta_p_t_deg->Fill(theta_p_t_deg);
	h1_theta_deg->Fill(theta_deg);
	h1_phi_deg->Fill(phi_deg);
	//
	h1_azimuth_deg->Fill(azimuth_deg);
	h1_altitude_deg->Fill(altitude_deg);
	//
	h1_x1_int->Fill(x1_int);
	h1_y1_int->Fill(y1_int);
	//
	h1_r_from_tel->Fill(r_from_tel);
	//
	n_trg_ev++;
      }
    }
    if(n_trg_ev>=max_trg_ev)
      break;    
  }
  //
  //
  //
  TFile* rootFile = new TFile(histOut.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<histOut.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<histOut.Data()<<endl;
  //
  h1_theta_p_t_deg->Write();
  h1_theta_deg->Write();
  h1_phi_deg->Write();
  //
  h1_azimuth_deg->Write();
  h1_altitude_deg->Write();
  //
  h1_x1_int->Write();
  h1_y1_int->Write();
  //
  h1_r_from_tel->Write();
  //
  rootFile->Close();
}

void cpv::Loop(TString histOut,
	       Double_t ellipse_A_r_km, Double_t ellipse_B_r_km, Double_t theta_p_t_deg_max,
	       Double_t y0_LST01_km, Double_t x0_LST01_km,
	       Int_t max_trg_ev,
	       Double_t norm_n_simtel,
	       Double_t norm_theta_p_t_deg_max,
	       Double_t norm_Energy_min,
	       Double_t norm_Energy_max,
	       Double_t norm_r_core_max,
	       Double_t norm_simtel_Energy_min, Double_t norm_simtel_Energy_max, TString evH_E2_data_out){
  //
  if(theta_p_t_deg_max == 0.0 && norm_theta_p_t_deg_max == 0.0){
    theta_p_t_deg_max = 0.01;
    norm_theta_p_t_deg_max = 0.01;
  }
  //
  Double_t norm_n_simtel_cphg = 0.0;
  //
  TH1D *h1_theta_deg = new TH1D("h1_theta_deg","h1_theta_deg",10000,0.0,181.0);
  TH1D *h1_phi_deg   = new TH1D("h1_phi_deg","h1_phi_deg",10000,-360.0,360.0);
  TH1D *h1_theta_p_t_deg = new TH1D("h1_theta_p_t_deg","h1_theta_p_t_deg",1000,0.0,181);
  //
  TH1D *h1_azimuth_deg = new TH1D("h1_azimuth_deg","h1_azimuth_deg",400,0.0,360);
  TH1D *h1_altitude_deg = new TH1D("h1_altitude_deg","h1_altitude_deg",400,0.0,180);
  //
  TH1D *h1_x1_int = new TH1D("h1_x1_int","h1_x1_int",400,-1.8,1.8);
  TH1D *h1_y1_int = new TH1D("h1_y1_int","h1_y1_int",400,-1.8,1.8);
  //
  TH1D *h1_r_from_tel = new TH1D("h1_r_from_tel","h1_r_from_tel",1000,0.0,2.0);
  //
  Double_t val_Emin = 1.0;      // GeV
  Double_t val_Emax = 100000;   // GeV
  Int_t val_N_bins_E = 25;
  //
  Double_t val_Thetamin = 0.0;  // deg
  Double_t val_Thetamax = 10.0; // deg
  Int_t val_N_bins_t = 10;
  //
  //
  evstHist *evH = new evstHist("evH","evH",
			       val_Emin, val_Emax, val_N_bins_E,
			       val_Thetamin, val_Thetamax, val_N_bins_t);  
  //
  evstHist *evH_E2 = new evstHist("evH_E2","evH_E2",
				  val_Emin, val_Emax, val_N_bins_E,
				  val_Thetamin, val_Thetamax, val_N_bins_t);  
  //
  Double_t theta_deg;
  Double_t phi_deg;
  Double_t theta_p_t;
  Double_t theta_p_t_deg;
  Double_t azimuth_deg;
  Double_t altitude_deg;
  Double_t r_from_tel;
  //
  Int_t n_trg_ev = 0;
  //
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  TVector3 v_det(1.0*TMath::Sin(20.0/180.0*TMath::Pi()),0.0,1.0*TMath::Cos(20.0/180.0*TMath::Pi()));
  ////////
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if(jentry%100000000 == 0)
      cout<<jentry<<endl;
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //
    TVector2 v_core(x1_int,y1_int);
    Double_t r_core = TMath::Sqrt(x1_int*x1_int+y1_int*y1_int);
    Double_t r_ellipse = TMath::Sqrt(ellipse_A_r_km*TMath::Cos(v_core.Phi())*ellipse_A_r_km*TMath::Cos(v_core.Phi()) +
				     ellipse_B_r_km*TMath::Sin(v_core.Phi())*ellipse_B_r_km*TMath::Sin(v_core.Phi()));
    if(r_core<=r_ellipse){
      //
      theta_deg = theta*180.0/TMath::Pi();
      phi_deg = conv_phi(phi)*180.0/TMath::Pi();
      //
      TVector3 v_prot(-vx,-vy,-vz);
      theta_p_t = TMath::ACos(v_prot.Dot(v_det)/v_prot.Mag()/v_det.Mag());
      theta_p_t_deg = theta_p_t*180/TMath::Pi();
      //
      TVector3 v_prot_azimuth_alt;
      v_prot_azimuth_alt.SetMagThetaPhi(1.0,theta,conv_phi(phi));
      TVector3 v_prot_azimuth_alt_inv(-v_prot_azimuth_alt.x(),-v_prot_azimuth_alt.y(),-v_prot_azimuth_alt.z());
      azimuth_deg = conv_phi(phi)*180/TMath::Pi();
      altitude_deg = 90.0 - v_prot_azimuth_alt_inv.Theta()*180/TMath::Pi();
      ///////////
      ///////////
      if(theta_p_t_deg<=theta_p_t_deg_max){
	//
	r_from_tel = TMath::Sqrt((x0_LST01_km - x1_int)*(x0_LST01_km - x1_int) + (y0_LST01_km - y1_int)*(y0_LST01_km - y1_int));
	//
	h1_theta_p_t_deg->Fill(theta_p_t_deg);
	h1_theta_deg->Fill(theta_deg);
	h1_phi_deg->Fill(phi_deg);
	//
	h1_azimuth_deg->Fill(azimuth_deg);
	h1_altitude_deg->Fill(altitude_deg);
	//
	h1_x1_int->Fill(x1_int);
	h1_y1_int->Fill(y1_int);
	//
	h1_r_from_tel->Fill(r_from_tel);
	//
	if(theta_p_t_deg_max == 0.01 && norm_theta_p_t_deg_max == 0.01){
	  evH->Fill(theta_p_t_deg+0.001, 1.1);
	  evH->Fill_rcore(theta_p_t_deg+0.001, 1.1, r_from_tel*1000);
	}
	else{
	  evH->Fill(theta_p_t_deg, 1.1);
	  evH->Fill_rcore(theta_p_t_deg, 1.1, r_from_tel*1000);
	}
	//
	n_trg_ev++;
	//
	if(theta_p_t_deg<norm_theta_p_t_deg_max){
	  if(r_from_tel<=norm_r_core_max){
	    norm_n_simtel_cphg++;
	  }
	}
      }
      if(n_trg_ev>=max_trg_ev)
	break;
    }
  }
  //
  //
  //
  //
  Double_t energy_integral_val_Emin_val_Emax = get_tot_flux_simulation(val_Emin, val_Emax);
  Double_t energy_integral_simtel = get_tot_flux_simulation(norm_simtel_Energy_min, norm_simtel_Energy_max);
  Double_t energy_integral_simtel_cut = get_tot_flux_simulation(norm_Energy_min, norm_Energy_max);
  //
  Double_t tot_simtel_events = norm_n_simtel/(norm_n_simtel_cphg/evH->Integral())/(energy_integral_simtel_cut/energy_integral_simtel);
  Double_t tot_simtel_events_extended = tot_simtel_events/(energy_integral_simtel/energy_integral_val_Emin_val_Emax);
  //
  Double_t fluence_E2;
  //
  cout<<"evH->Integral()                 "<<evH->Integral()<<endl
      <<"norm_n_simtel                   "<<norm_n_simtel<<endl
      <<"norm_n_simtel_cphg              "<<norm_n_simtel_cphg<<endl
      <<"energy_integral_simtel          "<<energy_integral_simtel<<endl
      <<"energy_integral_simtel_cut      "<<energy_integral_simtel_cut<<endl
      <<"tot_simtel_events               "<<tot_simtel_events<<endl
      <<"tot_simtel_events_extended      "<<tot_simtel_events_extended<<endl;
  //
  Int_t i_cell;
  Double_t e_min_integral;
  Double_t e_max_integral;
  //
  Double_t th_min_integral;
  Double_t th_max_integral;
  //
  Double_t e_bin_center;
  Double_t th_bin_center;
  //
  Int_t core_hist_id;
  Double_t bin_cont_core;
  //
  for(Int_t i_theta = 0;i_theta<val_N_bins_t;i_theta++){
    for(Int_t i_E = 0;i_E<val_N_bins_E;i_E++){
      i_cell = (i_E)*val_N_bins_t+(i_theta+1);
      //evH->Fill_rcore(theta_p_t_deg, 1.1, r_from_tel*1000);
      //
      e_min_integral = evH->get_E_hist()->GetBinLowEdge(i_E+1);
      e_max_integral = evH->get_E_hist()->GetBinLowEdge(i_E+1) + evH->get_E_hist()->GetBinWidth(i_E+1);
      th_min_integral = evH->get_theta_hist()->GetBinLowEdge(i_theta+1);
      th_max_integral = evH->get_theta_hist()->GetBinLowEdge(i_theta+1) + evH->get_theta_hist()->GetBinWidth(i_theta+1);
      e_bin_center = (e_max_integral + e_min_integral)/2.0;
      th_bin_center = (th_max_integral + th_min_integral)/2.0;
      //
      if(energy_integral_val_Emin_val_Emax*evH->GetBinContent(i_theta+1)>0.0 && evH->Integral()>0.0)
	fluence_E2 = tot_simtel_events_extended*get_tot_flux_simulation(e_min_integral,e_max_integral)/energy_integral_val_Emin_val_Emax*evH->GetBinContent(i_theta+1)/evH->Integral();
      else
	fluence_E2 = 0.0;
      evH_E2->SetBinContent( i_cell, fluence_E2);
      core_hist_id=evH_E2->get_bin_ID(e_bin_center, th_bin_center)-1;
      //cout<<"core_hist_id = "<<core_hist_id<<endl;
      for(Int_t i_core=0;i_core<evH_E2->get_v_r().at(core_hist_id)->GetNbinsX();i_core++){
	if(evH->get_v_r().at(evH->get_bin_ID(1.1, th_bin_center)-1)->Integral()>0.0)
	  bin_cont_core = evH->get_v_r().at(evH->get_bin_ID(1.1, th_bin_center)-1)->GetBinContent(i_core+1)/evH->get_v_r().at(evH->get_bin_ID(1.1, th_bin_center)-1)->Integral();
	else
	  bin_cont_core = 0.0;
	bin_cont_core *= fluence_E2;
	//cout<<"i_core = "<<i_core<<endl
	//  <<"bin_cont_core = "<<bin_cont_core<<endl;
	evH_E2->get_v_r().at(core_hist_id)->SetBinContent(i_core+1, bin_cont_core);
      }
    }
  } 
  //
  //
  TFile* rootFile = new TFile(histOut.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<histOut.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<histOut.Data()<<endl;
  //
  for(unsigned int i = 0;i<evH_E2->get_v_r().size();i++)
    evH_E2->get_v_r().at(i)->Write();
  //
  h1_theta_deg->Write();
  h1_phi_deg->Write();
  h1_theta_p_t_deg->Write();
  //
  h1_azimuth_deg->Write();
  h1_altitude_deg->Write();
  //
  h1_x1_int->Write();
  h1_y1_int->Write();
  //
  h1_r_from_tel->Write();
  //
  evH->Draw_hist("")->Write();
  evH_E2->Draw_hist("")->Write();
  //
  evH_E2->DumpBinContent(evH_E2_data_out, true);
  //
  rootFile->Close();
}

void cpv::Loop(TString histOut){
  //
  //gamma diff
  //const Double_t ellipse_A_r_km = 1.1126;
  //const Double_t ellipse_B_r_km = 1.0000;
  //proton diff
  const Double_t ellipse_A_r_km = 1.7321;
  const Double_t ellipse_B_r_km = 1.5000;
  //
  //
  Double_t val_Emin = 1.0;    // GeV
  Double_t val_Emax = 100000; // GeV
  Int_t val_N_bins_E = 25;
  //
  Double_t val_Thetamin = 0.0;  //deg
  Double_t val_Thetamax = 10.0; //deg
  Int_t val_N_bins_t = 10;
  //
  TGraph2D *gr2D_int1 = new TGraph2D();
  gr2D_int1->SetNameTitle("gr2D_int1", "gr2D_int1");
  TGraph2D *gr2D_sphere = new TGraph2D();
  gr2D_sphere->SetNameTitle("gr2D_sphere","gr2D_sphere");
  TGraph2D *gr2D_gen = new TGraph2D();
  gr2D_gen->SetNameTitle("gr2D_gen","gr2D_gen");
  //
  TH1D *h1_theta = new TH1D("h1_theta","h1_theta",1000,0.0,2*TMath::Pi());
  TH1D *h1_phi   = new TH1D("h1_phi","h1_phi",1000,0.0,2*TMath::Pi());
  TH1D *h1_theta_deg = new TH1D("h1_theta_deg","h1_theta_deg",10000,0.0,181.0);
  TH1D *h1_phi_deg   = new TH1D("h1_phi_deg","h1_phi_deg",10000,-360.0,360.0);
  TH1D *h1_theta_p_t = new TH1D("h1_theta_p_t","h1_theta_p_t",1000,0.0,2*TMath::Pi());
  TH1D *h1_theta_p_t_deg = new TH1D("h1_theta_p_t_deg","h1_theta_p_t_deg",1000,0.0,181);
  //
  TH1D *h1_azimuth_deg = new TH1D("h1_azimuth_deg","h1_azimuth_deg",400,0.0,360);
  TH1D *h1_altitude_deg = new TH1D("h1_altitude_deg","h1_altitude_deg",400,0.0,180);
  //
  TH1D *h1_theta_deg_cut = new TH1D("h1_theta_deg_cut","h1_theta_deg_cut",10000,0.0,181.0);
  TH1D *h1_phi_deg_cut = new TH1D("h1_phi_deg_cut","h1_phi_deg_cut",10000,-360.0,360.0);
  TH1D *h1_theta_p_t_deg_cut = new TH1D("h1_theta_p_t_deg_cut","h1_theta_p_t_deg_cut",val_N_bins_t,val_Thetamin,val_Thetamax);
  //
  //
  Double_t theta_deg;
  Double_t phi_deg;
  Double_t theta_p_t;
  Double_t theta_p_t_deg;
  //
  Double_t azimuth_deg;
  Double_t altitude_deg;
  //
  TH1D *h1_x0 = new TH1D("h1_x0","h1_x0",400,-1800,1800);
  TH1D *h1_x0_cut = new TH1D("h1_x0_cut","h1_x0_cut",400,-1800,1800);
  TH1D *h1_y0 = new TH1D("h1_y0","h1_y0",400,-1800,1800);
  TH1D *h1_y0_cut = new TH1D("h1_y0_cut","h1_y0_cut",400,-1800,1800);
  TH1D *h1_z0 = new TH1D("h1_z0","h1_z0",4000,6000,7000);
  TH1D *h1_z0_cut = new TH1D("h1_z0_cut","h1_z0_cut",4000,6000,7000);
  //
  TH1D *h1_x1_int = new TH1D("h1_x1_int","h1_x1_int",400,-1.8,1.8);
  TH1D *h1_x1_int_cut = new TH1D("h1_x1_int_cut","h1_x1_int_cut",200,-1.8,1.8);
  TH1D *h1_y1_int = new TH1D("h1_y1_int","h1_y1_int",400,-1.8,1.8);
  TH1D *h1_y1_int_cut = new TH1D("h1_y1_int_cut","h1_y1_int_cut",200,-1.8,1.8);
  //
  TH2D *h2_y1_int_vs_x1_int = new TH2D("h2_y1_int_vs_x1_int","h2_y1_int_vs_x1_int", 400,-1.8,1.8, 400,-1.8,1.8);
  TH2D *h2_y1_int_vs_x1_int_cut = new TH2D("h2_y1_int_vs_x1_int_cut","h2_y1_int_vs_x1_int_cut", 400,-1.8,1.8, 400,-1.8,1.8);
  //
  //
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  TVector3 v_det(1.0*TMath::Sin(20.0/180.0*TMath::Pi()),0.0,1.0*TMath::Cos(20.0/180.0*TMath::Pi()));
  ////////
  //
  evstHist *evH = new evstHist("evH","evH",
			       val_Emin, val_Emax, val_N_bins_E,
			       val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_integral = new evstHist("evH_integral","evH_integral",
					val_Emin, val_Emax, val_N_bins_E,
					val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_simtel_all = new evstHist("evH_simtel_all","evH_simtel_all",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_integral_diff_protons = new evstHist("evH_integral_diff_protons","evH_integral_diff_protons",
						     val_Emin, val_Emax, val_N_bins_E,
						     val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_integral_diff_gammas = new evstHist("evH_integral_diff_gammas","evH_integral_diff_gammas",
						     val_Emin, val_Emax, val_N_bins_E,
						     val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_integral_diff_ratio = new evstHist("evH_integral_diff_ratio","evH_integral_diff_ratio",
						     val_Emin, val_Emax, val_N_bins_E,
						     val_Thetamin, val_Thetamax, val_N_bins_t);  
  //
  //evH->test();
  ////////
  //evH->Draw_hist("evH.pdf");
  ////////
  Int_t nnpoints = 0;
  Int_t ngr2_points = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if(jentry%10000000 == 0)
      cout<<jentry<<endl;
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //  
    TVector2 v_core(x1_int,y1_int);
    Double_t r_core = TMath::Sqrt(x1_int*x1_int+y1_int*y1_int);
    Double_t r_ellipse = TMath::Sqrt(ellipse_A_r_km*TMath::Cos(v_core.Phi())*ellipse_A_r_km*TMath::Cos(v_core.Phi()) +
				     ellipse_B_r_km*TMath::Sin(v_core.Phi())*ellipse_B_r_km*TMath::Sin(v_core.Phi()));
    if(r_core<=r_ellipse){
      //
      if(ngr2_points<10000){
	gr2D_int1->SetPoint(gr2D_int1->GetN(),x1_int,y1_int,z1_int);
	gr2D_gen->SetPoint(gr2D_gen->GetN(),x0,y0,z0);
	ngr2_points++;
      }
      //
      theta_deg = theta*180.0/TMath::Pi();
      phi_deg = conv_phi(phi)*180.0/TMath::Pi();
      //
      h1_theta->Fill(theta);
      h1_phi->Fill(conv_phi(phi));
      //
      h1_theta_deg->Fill(theta_deg);
      h1_phi_deg->Fill(phi_deg);
      //
      TVector3 v_prot(-vx,-vy,-vz);
      theta_p_t = TMath::ACos(v_prot.Dot(v_det)/v_prot.Mag()/v_det.Mag());
      theta_p_t_deg = theta_p_t*180/TMath::Pi();
      //
      h1_theta_p_t->Fill(theta_p_t);
      h1_theta_p_t_deg->Fill(theta_p_t_deg);
      //
      h1_x0->Fill(x0);
      h1_y0->Fill(y0);
      h1_z0->Fill(z0);
      //
      h1_x1_int->Fill(x1_int);
      h1_y1_int->Fill(y1_int);
      h2_y1_int_vs_x1_int->Fill(x1_int,y1_int);
      //
      TVector3 v_prot_azimuth_alt;
      v_prot_azimuth_alt.SetMagThetaPhi(1.0,theta,conv_phi(phi));
      TVector3 v_prot_azimuth_alt_inv(-v_prot_azimuth_alt.x(),-v_prot_azimuth_alt.y(),-v_prot_azimuth_alt.z());
      //azimuth_deg = v_prot_azimuth_alt_inv.Phi()*180/TMath::Pi();
      azimuth_deg = conv_phi(phi)*180/TMath::Pi();
      altitude_deg = 90.0 - v_prot_azimuth_alt_inv.Theta()*180/TMath::Pi();
      ///////////
      //
      Double_t x0_LST01 = -70.93/1000.0;
      Double_t y0_LST01 = -52.07/1000.0;
      Double_t r_from_tel = TMath::Sqrt((x0_LST01 - x1_int)*(x0_LST01 - x1_int) + (y0_LST01 - y1_int)*(y0_LST01 - y1_int));
      //
      //
      ///////////
      if(theta_p_t_deg<=10.0){
 	//if(theta_p_t_deg<3.0){
	//if(r_from_tel<=0.15){
	//if(theta_p_t_deg<2.0){
	//if(r_from_tel<=0.07){
	//if(nnpoints<4492){
	//for(Int_t iii = 0; iii<49; iii++){
	h1_theta_deg_cut->Fill(theta_deg);
	h1_phi_deg_cut->Fill(phi_deg);
	h1_x0_cut->Fill(x0);
	h1_y0_cut->Fill(y0);
	h1_z0_cut->Fill(z0);
	//
	h1_x1_int_cut->Fill(x1_int);
	h1_y1_int_cut->Fill(y1_int);
	//
	h2_y1_int_vs_x1_int_cut->Fill(x1_int,y1_int);
	//
	h1_theta_p_t_deg_cut->Fill(theta_p_t_deg);
	//
	h1_azimuth_deg->Fill(azimuth_deg);
	h1_altitude_deg->Fill(altitude_deg);
	//
	nnpoints++;
	//}
	//}
	//}
      }
    }
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
  //
  //
  Double_t nev_theta;
  //
  // sphere 200 km protons
  //Double_t Ntot = 10*1000000000.0*501.0;
  //const Double_t earthR = 6371000;                          // m
  //const Double_t generationH = 200000;                      // m
  //const Double_t generationH_from_c = earthR + generationH; // m
  //const Double_t theta_generation = 15.0/180.0*TMath::Pi();
  //const Double_t ellipse_A_r = 1732.1;
  //const Double_t ellipse_B_r = 1500.0;
  //const Double_t solid_angle = 4*TMath::Pi();
  //
  //corsica 120 km protons
  //Double_t Ntot = 1000000000.0*38.0;
  //const Double_t earthR = 6371000;                          // m
  //const Double_t generationH = 120000;                      // m
  //const Double_t generationH_from_c = earthR + generationH; // m
  //const Double_t theta_generation = 15.0/180.0*TMath::Pi(); // 33.6539 + 1.5
  //const Double_t ellipse_A_r = 1732.1;
  //const Double_t ellipse_B_r = 1500.0;
  //const Double_t solid_angle = 2*TMath::Pi()*(1.0-TMath::Cos(theta_generation));
  //
  //corsica 120 km protons
  Double_t Ntot = 1000000000.0*38.0;
  const Double_t earthR = 6371000;                          // m
  const Double_t generationH = 120000;                      // m
  const Double_t generationH_from_c = earthR + generationH; // m
  const Double_t theta_generation = 15.0/180.0*TMath::Pi(); // 33.6539 + 1.5
  const Double_t ellipse_A_r = ellipse_A_r_km*1000.0;
  const Double_t ellipse_B_r = ellipse_B_r_km*1000.0;
  const Double_t solid_angle = 2*TMath::Pi()*(1.0-TMath::Cos(theta_generation));
  //  
  const Double_t effective_area = TMath::Pi()*ellipse_A_r*ellipse_B_r;
  //Double_t surfaceTotal = get_surface_fluxs(generationH_from_c, theta_generation);
  Double_t surfaceTotal = TMath::Pi()*(33.6539+1.5)*1000*(33.6539+1.5)*1000;
  //
  cout<<"surfaceTotal   "<<surfaceTotal<<endl
      <<"Ntot           "<<Ntot<<endl
      <<"effective_area "<<effective_area<<endl;
  //
  //
  Double_t solid_angle_integral;
  Double_t tot_flux_dE_perSR_simulation;
  Double_t ntot_to_sim = 0.0;
  Double_t ntot_to_sim_6degless = 0.0;
  Double_t ntot_to_sim_6degmore = 0.0;
  //
  Int_t i_cell;
  Double_t e_min_integral;
  Double_t e_max_integral;
  Double_t theta_min_integral;
  Double_t theta_max_integral;
  //
  Double_t tot_flux_dE_perM2_perS_perSR;
  Double_t tot_flux_in_the_acceptance;
  Double_t acceptance;
  for(Int_t i_theta = 0;i_theta<val_N_bins_t;i_theta++){
    //cout<<evH->get_theta_hist()->GetBinLowEdge(i_theta+1)<<endl;
    for(Int_t i_E = 0;i_E<val_N_bins_E;i_E++){
      i_cell = (i_E)*val_N_bins_t+(i_theta+1);
      //
      //cout<<evH->get_E_hist()->GetBinLowEdge(i_E+1)<<" "
      //<<evH->get_E_hist()->GetBinLowEdge(i_E+1) + evH->get_E_hist()->GetBinWidth(i_E+1)<<endl;
      //cout<<evH->get_theta_hist()->GetBinLowEdge(i_theta+1)<<endl;
      //
      nev_theta = h1_theta_p_t_deg_cut->GetBinContent(i_theta+1);
      acceptance = nev_theta/Ntot;
      //
      e_min_integral = evH->get_E_hist()->GetBinLowEdge(i_E+1);
      e_max_integral = evH->get_E_hist()->GetBinLowEdge(i_E+1) + evH->get_E_hist()->GetBinWidth(i_E+1);
      cout<<e_min_integral<<" "<<e_max_integral<<endl;
      tot_flux_dE_perM2_perS_perSR = get_tot_flux(e_min_integral,e_max_integral);
      //tot_flux_dE_perM2_perS_perSR = get_tot_flux_gamma(e_min_integral,e_max_integral);
      tot_flux_in_the_acceptance = surfaceTotal*tot_flux_dE_perM2_perS_perSR*solid_angle*acceptance;
      //
      evH->SetBinContent(i_cell,tot_flux_in_the_acceptance);
      //
      theta_min_integral = evH->get_theta_hist()->GetBinLowEdge(i_theta+1);
      theta_max_integral = evH->get_theta_hist()->GetBinLowEdge(i_theta+1) + evH->get_theta_hist()->GetBinWidth(i_theta+1);
      solid_angle_integral = 2*TMath::Pi()*(TMath::Cos(theta_min_integral/180.0*TMath::Pi()) - TMath::Cos(theta_max_integral/180.0*TMath::Pi()));
      evH_integral->SetBinContent(i_cell,tot_flux_dE_perM2_perS_perSR*solid_angle_integral*effective_area);
      //
      //tot_flux_dE_perSR_simulation = get_tot_flux_simulation(e_min_integral,e_max_integral)*solid_angle_integral*effective_area*5.44902e+09/779133.0/0.099991;
      //ntot_to_sim = ntot_to_sim + tot_flux_dE_perSR_simulation;
      //evH_simtel_all->SetBinContent(i_cell,tot_flux_dE_perSR_simulation);
      //
      ///////////////
      ///////////////
      evH_integral_diff_protons->SetBinContent(i_cell,get_tot_flux(e_min_integral,e_max_integral)*solid_angle_integral);
      evH_integral_diff_gammas->SetBinContent(i_cell,get_tot_flux_gamma(e_min_integral,e_max_integral)*solid_angle_integral);
      ///////////////
      ///////////////
      //
      //*1.05433e+09/333646.42
      tot_flux_dE_perSR_simulation = get_tot_flux_simulation(e_min_integral,e_max_integral)*solid_angle_integral*effective_area*43821.9;
      ntot_to_sim = ntot_to_sim + tot_flux_dE_perSR_simulation;
      if(theta_max_integral<=6.0)
	ntot_to_sim_6degless += tot_flux_dE_perSR_simulation;
      else
	ntot_to_sim_6degmore += tot_flux_dE_perSR_simulation;
      //
      evH_simtel_all->SetBinContent(i_cell,tot_flux_dE_perSR_simulation);
    }
  }
  //for(Int_t i = 1;i<=evH->GetNcells();i++){
  //h1_theta_p_t_deg_cut
  //
  //}
  //evH->test();
  ////////
  evH_integral_diff_ratio->SetMinimum(1.0e-5);
  evH_integral_diff_ratio->SetMaximum(1.0e-4);
  evH_integral_diff_ratio->Divide(evH_integral_diff_gammas, evH_integral_diff_protons);
  evH_integral_diff_ratio->Draw_hist("evH_integral_diff_ratio.pdf");
  ////////
  evH->SetMaximum(1.0e+4);
  evH->SetMinimum(1.0e-6);
  evH->Draw_hist("evH_gamma_diff.pdf");
  //Proton rate, Hz
  evH_integral->SetMaximum(1.0e+4);
  evH_integral->SetMinimum(1.0e-6);
  evH_integral->Draw_hist("evH_integral_gamma_diff.pdf","Gamma bkg. galactic rate, Hz");
  evH_simtel_all->SetMaximum(1.0e+10);
  evH_simtel_all->SetMinimum(1.0e+2);
  evH_simtel_all->Draw_hist("evH_simtel_all_gamma_diff.pdf");
  //
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
  h1_theta_p_t->Write();
  h1_theta_p_t_deg->Write();
  //
  h1_theta_deg_cut->Write();
  h1_phi_deg_cut->Write();
  //
  h1_azimuth_deg->Write();
  h1_altitude_deg->Write();
  //
  h1_x0->Write();
  h1_y0->Write();
  h1_z0->Write();
  h1_x0_cut->Write();
  h1_y0_cut->Write();
  h1_z0_cut->Write();
  //
  h1_x1_int->Write();
  h1_y1_int->Write();
  h1_x1_int_cut->Write();
  h1_y1_int_cut->Write();
  //
  h2_y1_int_vs_x1_int->Write();
  h2_y1_int_vs_x1_int_cut->Write();
  //
  h1_theta_p_t_deg_cut->Write();
  //
  evH->Write();
  evH_simtel_all->Write();
  evH->Draw_hist("")->Write();
  evH_integral->Write();
  evH_integral->Draw_hist("")->Write();
  evH_simtel_all->Draw_hist("")->Write();
  //
  evH->DumpBinContent("flux.dat");
  evH_simtel_all->DumpBinContent("flux_simtel_all_gamma_diff.dat");
  //
  cout<<"ntot_to_sim                      "<<ntot_to_sim<<endl
      <<"evH_simtel_all->GetTotIntegral() "<<evH_simtel_all->GetTotIntegral()<<endl
      <<"                                 "<<evH_simtel_all->GetIntegral(val_Emin,val_Emax, val_Thetamin, val_Thetamax)<<endl
      <<"                                 "<<evH_simtel_all->GetIntegral(10,val_Emax, val_Thetamin, val_Thetamax)<<endl
      <<"                                 "<<evH_simtel_all->GetIntegral(10,val_Emax, val_Thetamin, val_Thetamax)/evH_simtel_all->GetIntegral(val_Emin,val_Emax, val_Thetamin, val_Thetamax)<<endl;
  //
  cout<<"get_tot_flux(50,3000.0)                  = "<<get_tot_flux(50,3000.0)<<endl
      <<"get_tot_flux_gamma(50.0,3000.0)          = "<<get_tot_flux_gamma(50.0,3000.0)<<endl
      <<"get_tot_flux_MAGIC_tel_crab(50.0,3000.0) = "<<get_tot_flux_MAGIC_tel_crab(1000000, 50.0, 3000.0)<<endl;
  //
  cout<<"ntot_to_sim_6degless = "<<ntot_to_sim_6degless<<endl
      <<"ntot_to_sim_6degmore = "<<ntot_to_sim_6degmore<<endl;
  //
  //
  solid_angle_integral = 2*TMath::Pi()*(TMath::Cos(0.0/180.0*TMath::Pi()) - TMath::Cos(6.0/180.0*TMath::Pi()));
  tot_flux_dE_perM2_perS_perSR = get_tot_flux_simulation(5.0,50000.0);
  Double_t not_norm_tot_fluence=tot_flux_dE_perM2_perS_perSR*solid_angle_integral*effective_area;
  Double_t norm_tot_fluence=1.05433e+09;
  Double_t correction_tot_fluence=norm_tot_fluence/not_norm_tot_fluence;  
  //
  cout<<"solid_angle_integral         "<<solid_angle_integral<<endl
      <<"tot_flux_dE_perM2_perS_perSR "<<tot_flux_dE_perM2_perS_perSR<<endl
      <<"not_norm_tot_fluence         "<<not_norm_tot_fluence<<endl
      <<"norm_tot_fluence             "<<norm_tot_fluence<<endl
      <<"correction_tot_fluence       "<<correction_tot_fluence<<endl;  
  //
  rootFile->Close();
}

Double_t cpv::get_tot_flux(Double_t e_min, Double_t e_max){
  return 5781.68*(1.0/TMath::Power(e_min,1.674)-1.0/TMath::Power(e_max,1.674));
}

Double_t cpv::get_tot_flux_gamma(Double_t e_min, Double_t e_max){
  return 1.5611068/TMath::Power(e_max,1.51752)-1.7976545/TMath::Power(e_max,1.52789) - 1.5611068/TMath::Power(e_min,1.51752) + 1.7976545/TMath::Power(e_min,1.52789);
}

Double_t cpv::function_log_parabola_fit_MAGIC_tel_crab_GeV(Double_t e){
  if(e<50.0)
    return 0.0;
  else if(e>3000.0)
    return 0.0;
  return 10000.0*((3.23)*1.0e-14)*TMath::Power((e/1.00000/1000.0),(-2.47 - 0.24*TMath::Log(e/1.00000/1000.0)));
}

//e_min = 5.0;
//e_max = 4000.0;
Double_t cpv::get_diff_flux_ele_pos_power_law(Double_t e){
  if(e<5.0)
    return 0.0;
  else if(e>4000.0)
    return 0.0;
  Double_t A1 = -1.86680e-09;
  Double_t B1 = 3.52962e+03;
  Double_t C1 = 2.08619e+00;
  Double_t A2 = 2.33684e-08;
  Double_t B2 = 1.80334e+03;
  Double_t C2 = 3.06795e+00;
  Double_t val = A1*(TMath::Power(e/B1,-C1)) + A2*(TMath::Power(e/B2,-C2));
  return val;
}

//e_min = 5.0;
//e_max = 4000.0;
Double_t cpv::get_tot_flux_ele_pos(Int_t nn, Double_t e_min, Double_t e_max){
  Double_t de = (e_max - e_min)/(nn-1);
  Double_t f;
  Double_t e;
  Double_t integral = 0.0;
  for(Int_t i = 0;i<nn;i++){
    e = de*i + e_min;
    f=get_diff_flux_ele_pos_power_law(e);
    integral += f*de;
  }
  return integral;
}

//e_min = 50.0;
//e_max = 3000.0;
Double_t cpv::get_tot_flux_MAGIC_tel_crab(Int_t nn, Double_t e_min, Double_t e_max){
  Double_t de = (e_max - e_min)/(nn-1);
  Double_t f;
  Double_t e;
  Double_t integral = 0.0;
  for(Int_t i = 0;i<nn;i++){
    e = de*i + e_min;
    f=function_log_parabola_fit_MAGIC_tel_crab_GeV(e);
    integral += f*de;
  }
  return integral;
}

Double_t cpv::get_tot_flux_simulation(Double_t e_min, Double_t e_max){
  return (1.0/TMath::Power(e_min,1.0)-1.0/TMath::Power(e_max,1.0));
}

Double_t cpv::get_surface_fluxs(Double_t r_in_m, Double_t theta){
  return 2.0*TMath::Pi()*r_in_m*r_in_m*(1-TMath::Cos(theta));
}

Double_t cpv::conv_phi(Double_t phiv){
  if(phiv<0)
    return 2.0*TMath::Pi()+phiv;
  return phiv;
}

void cpv::calculate_rate(Double_t ellipse_min_r, TString pdf_out, TString dat_out,
			 TString frame_title, Double_t min_f, Double_t max_f, TString part_type){
  //
  Double_t val_Emin = 1.0;      // GeV
  Double_t val_Emax = 100000;   // GeV
  Int_t val_N_bins_E = 25;
  //
  Double_t val_Thetamin = 0.0;  // deg
  Double_t val_Thetamax = 10.0; // deg
  Int_t val_N_bins_t = 10;
  //
  evstHist *evH_flux = new evstHist("evH_flux","evH_flux",
				    val_Emin, val_Emax, val_N_bins_E,
				    val_Thetamin, val_Thetamax, val_N_bins_t);
  Double_t area = TMath::Pi()*ellipse_min_r*ellipse_min_r;
  //
  Double_t solid_angle;
  //
  Int_t i_cell;
  Double_t e_min_integral;
  Double_t e_max_integral;
  Double_t theta_min_integral;
  Double_t theta_max_integral;
  //
  Double_t flux_dE_perM2_perS_perSR;
  Double_t flux;
  //
  for(Int_t i_theta = 0;i_theta<val_N_bins_t;i_theta++){
    for(Int_t i_E = 0;i_E<val_N_bins_E;i_E++){
      i_cell = (i_E)*val_N_bins_t+(i_theta+1);
      if(part_type == "gamma_crab")
	i_cell = (i_E)*val_N_bins_t+(0+1);
      //
      e_min_integral = evH_flux->get_E_hist()->GetBinLowEdge(i_E+1);
      e_max_integral = evH_flux->get_E_hist()->GetBinLowEdge(i_E+1) + evH_flux->get_E_hist()->GetBinWidth(i_E+1);
      //
      theta_min_integral = evH_flux->get_theta_hist()->GetBinLowEdge(i_theta+1);
      theta_max_integral = evH_flux->get_theta_hist()->GetBinLowEdge(i_theta+1) + evH_flux->get_theta_hist()->GetBinWidth(i_theta+1);
      //
      if(part_type == "proton_diff")
	flux_dE_perM2_perS_perSR = get_tot_flux(e_min_integral,e_max_integral);
      else if(part_type == "gamma_diff_galactic")
	flux_dE_perM2_perS_perSR = get_tot_flux_gamma(e_min_integral,e_max_integral);
      else if(part_type == "gamma_crab")
	flux_dE_perM2_perS_perSR = get_tot_flux_MAGIC_tel_crab( 1000000, e_min_integral, e_max_integral);
      else if(part_type == "ele_pos")
	flux_dE_perM2_perS_perSR = get_tot_flux_ele_pos( 1000000, e_min_integral, e_max_integral);
      //
      solid_angle = 2*TMath::Pi()*(TMath::Cos(theta_min_integral/180.0*TMath::Pi()) - TMath::Cos(theta_max_integral/180.0*TMath::Pi()));      
      if(part_type == "gamma_crab")
	solid_angle = 1.0;
      flux = flux_dE_perM2_perS_perSR*area*solid_angle;
      //
      evH_flux->SetBinContent(i_cell,flux);
      //
    }
  }
  //
  evH_flux->SetMinimum(min_f);
  evH_flux->SetMaximum(max_f);
  evH_flux->Draw_hist(pdf_out.Data(),frame_title.Data());
  evH_flux->DumpBinContent(dat_out.Data());
  cout<<part_type <<" tot flux : "<<evH_flux->Integral()<<endl;
}

void cpv::calculate_diff_proton_rate(){
  Double_t ellipse_min_r = 1500.0;
  TString pdf_out = "flux_diff_protons.pdf";
  TString dat_out = "flux_diff_protons.dat";
  TString frame_title = "Proton flux, Hz";
  Double_t min_f = 0.1;
  Double_t max_f = 4*1.0e+8;
  TString part_type = "proton_diff";
  //
  calculate_rate(ellipse_min_r, pdf_out, dat_out,
		 frame_title, min_f, max_f, part_type);
}

void cpv::calculate_diff_gamma_rate(){
  Double_t ellipse_min_r = 1000.0;
  TString pdf_out = "flux_gamma_diff_galactic.pdf";
  TString dat_out = "flux_gamma_diff_galactic.dat";
  TString frame_title = "Gamma diff galactic flux, Hz";
  Double_t min_f = 1.0e-6;
  Double_t max_f = 1.0e+4;
  TString part_type = "gamma_diff_galactic";
  //
  calculate_rate(ellipse_min_r, pdf_out, dat_out,
		 frame_title, min_f, max_f, part_type);
}

void cpv::calculate_diff_electron_rate(){
  Double_t ellipse_min_r = 1000.0;
  TString pdf_out = "flux_ele_pos_diff.pdf";
  TString dat_out = "flux_ele_pos_diff.dat";
  TString frame_title = "e+e- diff flux, Hz";
  Double_t min_f = 1.0e-6;
  Double_t max_f = 1.0e+4;
  TString part_type = "ele_pos";
  //
  calculate_rate(ellipse_min_r, pdf_out, dat_out,
		 frame_title, min_f, max_f, part_type);

}

void cpv::calculate_crab_gamma_rate(){
  Double_t ellipse_min_r = 800.0;
  TString pdf_out = "flux_gamma_crab.pdf";
  TString dat_out = "flux_gamma_crab.dat";
  TString frame_title = "Gamma crab flux, Hz";
  Double_t min_f = 1.0e-3;
  Double_t max_f = 1.0e+1;
  TString part_type = "gamma_crab";
  //
  calculate_rate(ellipse_min_r, pdf_out, dat_out,
		 frame_title, min_f, max_f, part_type);
}
