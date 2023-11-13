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

void cpv::Loop(TString histOut){
  //
  const Double_t ellipse_A_r_km = 1.1126;
  const Double_t ellipse_B_r_km = 1.0000;
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
  //evH->test();
  ////////
  //evH->Draw_hist("evH.pdf");
  ////////
  Int_t nnpoints = 0;
  Int_t ngr2_points = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if(jentry%1000000 == 0)
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
      if(theta_p_t_deg<10.0){
 	//if(theta_p_t_deg<3.0){
	//if(r_from_tel<=0.15){
	if(theta_p_t_deg<2.0){
	  if(r_from_tel<=0.07){
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
	  }
	}
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
  const Double_t ellipse_A_r = ellipse_A_r_km/1000.0;
  const Double_t ellipse_B_r = ellipse_B_r_km/1000.0;
  const Double_t solid_angle = 2*TMath::Pi()*(1.0-TMath::Cos(theta_generation));
  //  
  const Double_t effective_area = TMath::Pi()*ellipse_A_r*ellipse_B_r;
  //Double_t surfaceTotal = get_surface_fluxs(generationH_from_c, theta_generation);
  Double_t surfaceTotal = TMath::Pi()*(33.6539+1.5)*1000*(33.6539+1.5)*1000;
  cout<<"surfaceTotal   "<<surfaceTotal<<endl
      <<"Ntot           "<<Ntot<<endl
      <<"effective_area "<<effective_area<<endl;
  //
  //
  Double_t solid_angle_integral;
  Double_t tot_flux_dE_perSR_simulation;
  Double_t ntot_to_sim = 0.0;
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
      tot_flux_in_the_acceptance = surfaceTotal*tot_flux_dE_perM2_perS_perSR*solid_angle*acceptance;
      //
      evH->SetBinContent(i_cell,tot_flux_in_the_acceptance);
      //
      theta_min_integral = evH->get_theta_hist()->GetBinLowEdge(i_theta+1);
      theta_max_integral = evH->get_theta_hist()->GetBinLowEdge(i_theta+1) + evH->get_theta_hist()->GetBinWidth(i_theta+1);
      solid_angle_integral = 2*TMath::Pi()*(TMath::Cos(theta_min_integral/180.0*TMath::Pi()) - TMath::Cos(theta_max_integral/180.0*TMath::Pi()));
      evH_integral->SetBinContent(i_cell,tot_flux_dE_perM2_perS_perSR*solid_angle_integral*effective_area);
      //
      tot_flux_dE_perSR_simulation = get_tot_flux_simulation(e_min_integral,e_max_integral)*solid_angle_integral*effective_area*5.44902e+09/779133.0/0.099991;
      ntot_to_sim = ntot_to_sim + tot_flux_dE_perSR_simulation;
      evH_simtel_all->SetBinContent(i_cell,tot_flux_dE_perSR_simulation);
    }
  }
  //for(Int_t i = 1;i<=evH->GetNcells();i++){
  //h1_theta_p_t_deg_cut
  //
  //}
  //evH->test();
  ////////
  evH->SetMaximum(1.0e+9);
  evH->SetMinimum(1.0e-1);
  evH->Draw_hist("evH.pdf");
  evH_integral->Draw_hist("evH_integral.pdf");
  evH_simtel_all->Draw_hist("evH_simtel_all.pdf");
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
  evH->Draw_hist("")->Write();
  evH_integral->Write();
  evH_integral->Draw_hist("")->Write();
  evH_simtel_all->Draw_hist("")->Write();
  //
  evH->DumpBinContent("flux.dat");
  evH_simtel_all->DumpBinContent("flux_simtel_all.dat");
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
  rootFile->Close();
}

Double_t cpv::get_tot_flux(Double_t e_min, Double_t e_max){
  return 5781.68*(1.0/TMath::Power(e_min,1.674)-1.0/TMath::Power(e_max,1.674));
}

Double_t cpv::get_tot_flux_gamma(Double_t e_min, Double_t e_max){
  return 1.5611068/TMath::Power(e_max,1.51752)-1.7976545/TMath::Power(e_max,1.52789) - 1.5611068/TMath::Power(e_min,1.51752) + 1.7976545/TMath::Power(e_min,1.52789);
}

Double_t cpv::function_log_parabola_fit_MAGIC_tel_crab_GeV(Double_t e){
  return 10000.0*((3.23)*1.0e-14)*TMath::Power((e/1.00000/1000.0),(-2.47 - 0.24*TMath::Log(e/1.00000/1000.0)));
}

// e_min = 50.0;
// e_max = 3000.0;
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
