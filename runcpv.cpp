//my
#include "src/cpv.hh"

//root
#include "TROOT.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <time.h>

using namespace std;

int main(int argc, char *argv[]){
  if(argc == 4 && atoi(argv[1])==0){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    cpv a(rootFilesList);
    a.Loop(outRootFileF);
  }
  else if(argc == 4 && atoi(argv[1])==1){
    TString inRootFiles = argv[2];
    TString outRootFileF = argv[3];
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"inRootFiles   : "<<inRootFiles<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl;
    cpv a( inRootFiles, atoi(argv[1]));
    a.Loop(outRootFileF);
  }
  else if(argc == 2 && atoi(argv[1])==2){
    cpv a;
    a.calculate_diff_proton_rate();
    a.calculate_diff_gamma_rate();
    a.calculate_crab_gamma_rate();
    a.calculate_diff_electron_rate();
    //
    cout<<"crab      (800m)  = "<<a.get_tot_flux_MAGIC_tel_crab(1000000, 50.0, 3000.0)*TMath::Pi()*800.0*800.0<<endl;
    cout<<"prot      (1500m) = "<<a.get_tot_flux(50.0,3000.0)*2*TMath::Pi()*(TMath::Cos(0.0/180.0*TMath::Pi()) - TMath::Cos(0.5/180.0*TMath::Pi()))*TMath::Pi()*1500.0*1500.0<<endl;
    cout<<"gamma gal (1000m) = "<<a.get_tot_flux_gamma(50.0,3000.0)*2*TMath::Pi()*(TMath::Cos(0.0/180.0*TMath::Pi()) - TMath::Cos(0.5/180.0*TMath::Pi()))*TMath::Pi()*1000.0*1000.0<<endl;
    cout<<"ele pos   (1000m) = "<<a.get_tot_flux_ele_pos(1000000, 50.0, 3000.0)*2*TMath::Pi()*(TMath::Cos(0.0/180.0*TMath::Pi()) - TMath::Cos(0.5/180.0*TMath::Pi()))*TMath::Pi()*1000.0*1000.0<<endl;
    //
    cout<<"crab      (1000m) = "<<a.get_tot_flux_MAGIC_tel_crab(1000000, 50.0, 3000.0)*TMath::Pi()*1000.0*1000.0<<endl;
    cout<<"prot      (1000m) = "<<a.get_tot_flux(50.0,3000.0)*2*TMath::Pi()*(TMath::Cos(0.0/180.0*TMath::Pi()) - TMath::Cos(0.5/180.0*TMath::Pi()))*TMath::Pi()*1000.0*1000.0<<endl;
    cout<<"gamma gal (1000m) = "<<a.get_tot_flux_gamma(50.0,3000.0)*2*TMath::Pi()*(TMath::Cos(0.0/180.0*TMath::Pi()) - TMath::Cos(0.5/180.0*TMath::Pi()))*TMath::Pi()*1000.0*1000.0<<endl;
    cout<<"ele pos   (1000m) = "<<a.get_tot_flux_ele_pos(1000000, 50.0, 3000.0)*2*TMath::Pi()*(TMath::Cos(0.0/180.0*TMath::Pi()) - TMath::Cos(0.5/180.0*TMath::Pi()))*TMath::Pi()*1000.0*1000.0<<endl;
    //
  }
  else if(argc == 6 && atoi(argv[1])==3){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    TString part_type = argv[4];
    Int_t max_trg_ev = atoi(argv[5]);
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl
	<<"part_type     : "<<part_type<<endl
	<<"max_trg_ev    : "<<max_trg_ev<<endl;
    //
    Double_t ellipse_A_r_km;
    Double_t ellipse_B_r_km;
    Double_t theta_p_t_deg_max;
    Double_t x0_LST01_km;
    Double_t y0_LST01_km;
    //
    Double_t norm_n_simtel;
    Double_t norm_theta_p_t_deg_max;
    Double_t norm_Energy_min;
    Double_t norm_Energy_max;
    Double_t norm_r_core_max;
    //
    Double_t norm_simtel_Energy_min;
    Double_t norm_simtel_Energy_max;
    //
    TString evH_E2_data_out;
    //
    cpv a(rootFilesList);
    //
    if(part_type == "proton_diff"){
      //
      ellipse_A_r_km = 1.7321;
      ellipse_B_r_km = 1.5000;
      theta_p_t_deg_max = 10;
      //
      x0_LST01_km = -70.93/1000.0;
      y0_LST01_km = -52.07/1000.0;
      //
      norm_n_simtel = 10172;
      norm_theta_p_t_deg_max = 3.0;
      norm_Energy_min = 3.98107*1000.0;
      norm_Energy_max = 100.0*1000.0;
      norm_r_core_max = 150.0/1000.0;
      //
      norm_simtel_Energy_min = 10.0;
      norm_simtel_Energy_max = 100.0*1000.0;
      //
      evH_E2_data_out = "proton_diff_simtel.dat";
      //
      a.Loop(outRootFileF, ellipse_A_r_km, ellipse_B_r_km, theta_p_t_deg_max, x0_LST01_km, y0_LST01_km, max_trg_ev,
	     norm_n_simtel, norm_theta_p_t_deg_max, norm_Energy_min, norm_Energy_max, norm_r_core_max,
	     norm_simtel_Energy_min, norm_simtel_Energy_max, evH_E2_data_out);
    }
    else if(part_type == "gamma_diff"){
      //
      ellipse_A_r_km = 1.1126;
      ellipse_B_r_km = 1.0000;
      theta_p_t_deg_max = 6.0;
       //
      x0_LST01_km = -70.93/1000.0;
      y0_LST01_km = -52.07/1000.0;
      //
      norm_n_simtel = 12674;
      norm_theta_p_t_deg_max = 2.0;
      norm_Energy_min = 0.2*1000.0;
      norm_Energy_max = 50.0*1000.0;
      norm_r_core_max = 70.0/1000.0;
      //
      norm_simtel_Energy_min = 5.0;
      norm_simtel_Energy_max = 50.0*1000.0;
      //
      evH_E2_data_out = "gamma_diff_galactic_simtel.dat";
      //
      a.Loop(outRootFileF, ellipse_A_r_km, ellipse_B_r_km, theta_p_t_deg_max, x0_LST01_km, y0_LST01_km, max_trg_ev,
	     norm_n_simtel, norm_theta_p_t_deg_max, norm_Energy_min, norm_Energy_max, norm_r_core_max,
	     norm_simtel_Energy_min, norm_simtel_Energy_max, evH_E2_data_out);
    }
    else if(part_type == "ele_pos"){
      //
      ellipse_A_r_km = 1.1126;
      ellipse_B_r_km = 1.0000;
      theta_p_t_deg_max = 6.0;
       //
      x0_LST01_km = -70.93/1000.0;
      y0_LST01_km = -52.07/1000.0;
      //
      norm_n_simtel = 5704;
      norm_theta_p_t_deg_max = 2.0;
      norm_Energy_min = 0.2*1000.0;
      norm_Energy_max = 5.0*1000.0;
      norm_r_core_max = 70.0/1000.0;
      //
      norm_simtel_Energy_min = 5.0;
      norm_simtel_Energy_max = 5.0*1000.0;
      //
      evH_E2_data_out = "electron_simtel.dat";
      //
      a.Loop(outRootFileF, ellipse_A_r_km, ellipse_B_r_km, theta_p_t_deg_max, x0_LST01_km, y0_LST01_km, max_trg_ev,
	     norm_n_simtel, norm_theta_p_t_deg_max, norm_Energy_min, norm_Energy_max, norm_r_core_max,
	     norm_simtel_Energy_min, norm_simtel_Energy_max, evH_E2_data_out);
    }
  }
  else if(argc == 6 && atoi(argv[1])==4){
    TString rootFile = argv[2];
    TString outRootFileF = argv[3];
    TString part_type = argv[4];
    Int_t max_trg_ev = atoi(argv[5]);
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"rootFile     : "<<rootFile<<endl
	<<"outRootFileF : "<<outRootFileF<<endl
	<<"part_type    : "<<part_type<<endl
	<<"max_trg_ev   : "<<max_trg_ev<<endl;
    //
    Double_t ellipse_A_r_km;
    Double_t ellipse_B_r_km;
    Double_t theta_p_t_deg_max;
    Double_t x0_LST01_km;
    Double_t y0_LST01_km;
    //
    Double_t norm_n_simtel;
    Double_t norm_theta_p_t_deg_max;
    Double_t norm_Energy_min;
    Double_t norm_Energy_max;
    Double_t norm_r_core_max;
    //
    Double_t norm_simtel_Energy_min;
    Double_t norm_simtel_Energy_max;
    //
    TString evH_E2_data_out;
    //
    cpv a(rootFile,1);
    //
    if(part_type == "gamma"){
      //
      ellipse_A_r_km = 0.851;
      ellipse_B_r_km = 0.8;
      theta_p_t_deg_max = 0.0;
      //
      x0_LST01_km = -70.93/1000.0;
      y0_LST01_km = -52.07/1000.0;
      //
      norm_n_simtel = 15918;
      norm_theta_p_t_deg_max = 0.0;
      norm_Energy_min = 0.2*1000.0;
      norm_Energy_max = 50.0*1000.0;
      norm_r_core_max = 70.0/1000.0;
      //
      norm_simtel_Energy_min = 5.0;
      norm_simtel_Energy_max = 50.0*1000.0;
      //
      evH_E2_data_out = "gamma_on_axis_simtel.dat";
      //
      a.Loop(outRootFileF, ellipse_A_r_km, ellipse_B_r_km, theta_p_t_deg_max, x0_LST01_km, y0_LST01_km, max_trg_ev,
	     norm_n_simtel, norm_theta_p_t_deg_max, norm_Energy_min, norm_Energy_max, norm_r_core_max,
	     norm_simtel_Energy_min, norm_simtel_Energy_max, evH_E2_data_out);
    }
  }
  else if(argc == 6 && atoi(argv[1])==5){
    TString rootFilesList = argv[2];
    TString outRootFileF = argv[3];
    TString part_type = argv[4];
    Int_t max_trg_ev = atoi(argv[5]);
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"rootFilesList : "<<rootFilesList<<endl
	<<"outRootFileF  : "<<outRootFileF<<endl
	<<"part_type     : "<<part_type<<endl
	<<"max_trg_ev    : "<<max_trg_ev<<endl;
    //
    Double_t x0_LST01_km = -70.93/1000.0;
    Double_t y0_LST01_km = -52.07/1000.0;
    //
    Double_t norm_theta_p_t_deg_max;
    Double_t norm_r_core_max;
    //
    cpv a(rootFilesList);
    //
    if(part_type == "proton_diff"){
      //
      norm_theta_p_t_deg_max = 3.0;
      norm_r_core_max = 150.0/1000.0;
      //
      cout<<"norm_theta_p_t_deg_max "<<norm_theta_p_t_deg_max<<endl
	  <<"norm_r_core_max        "<<norm_r_core_max<<endl;
      //
      a.Loop(outRootFileF, y0_LST01_km, x0_LST01_km, max_trg_ev, norm_theta_p_t_deg_max, norm_r_core_max);
    }
    else if(part_type == "gamma_diff"){
      //
      norm_theta_p_t_deg_max = 2.0;
      norm_r_core_max = 70.0/1000.0;
      //
      cout<<"norm_theta_p_t_deg_max "<<norm_theta_p_t_deg_max<<endl
	  <<"norm_r_core_max        "<<norm_r_core_max<<endl;
      //
      a.Loop(outRootFileF, y0_LST01_km, x0_LST01_km, max_trg_ev, norm_theta_p_t_deg_max, norm_r_core_max);
    }
    else if(part_type == "ele_pos"){
      //
      norm_theta_p_t_deg_max = 2.0;
      norm_r_core_max = 70.0/1000.0;
      //
      cout<<"norm_theta_p_t_deg_max "<<norm_theta_p_t_deg_max<<endl
	  <<"norm_r_core_max        "<<norm_r_core_max<<endl;
      //
      a.Loop(outRootFileF, y0_LST01_km, x0_LST01_km, max_trg_ev, norm_theta_p_t_deg_max, norm_r_core_max);
    }
  }
  else if(argc == 6 && atoi(argv[1])==6){
    TString rootFile = argv[2];
    TString outRootFileF = argv[3];
    TString part_type = argv[4];
    Int_t max_trg_ev = atoi(argv[5]);
    cout<<"--> Parameter calculation from the WF <--"<<endl
	<<"rootFile     : "<<rootFile<<endl
	<<"outRootFileF : "<<outRootFileF<<endl
	<<"part_type    : "<<part_type<<endl
	<<"max_trg_ev   : "<<max_trg_ev<<endl;
    //
    Double_t x0_LST01_km = -70.93/1000.0;
    Double_t y0_LST01_km = -52.07/1000.0;
    //
    Double_t norm_theta_p_t_deg_max;
    Double_t norm_r_core_max;
    //
    cpv a(rootFile,1);
    //
    if(part_type == "gamma"){
      //
      norm_theta_p_t_deg_max = 2.0;
      norm_r_core_max = 70.0/1000.0;
      //
      cout<<"norm_theta_p_t_deg_max "<<norm_theta_p_t_deg_max<<endl
	  <<"norm_r_core_max        "<<norm_r_core_max<<endl;
      //
      a.Loop(outRootFileF, y0_LST01_km, x0_LST01_km, max_trg_ev, norm_theta_p_t_deg_max, norm_r_core_max);
    }
  }
  else{
    cout<<" --> ERROR in input arguments "<<endl
	<<" runID [1] = 0 (execution ID number)"<<endl
      	<<"       [2] - file with list of the root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 1 (execution ID number)"<<endl
      	<<"       [2] - in root files"<<endl
	<<"       [3] - name of root file with histograms"<<endl;
    cout<<" runID [1] = 2 (calculate diff. proton, diff. gamma, diff. electron and crab gamma rates)"<<endl;
    cout<<" runID [1] = 3 (execution ID number)"<<endl
      	<<"       [2] - in root file list"<<endl
	<<"       [3] - name of root file with histograms"<<endl
      	<<"       [4] - part type (proton_diff, gamma_diff, ele_pos)"<<endl
	<<"       [5] - max_trg_ev"<<endl;
    cout<<" runID [1] = 4 (execution ID number)"<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl
      	<<"       [4] - part type (gamma)"<<endl
	<<"       [5] - max_trg_ev"<<endl;
    cout<<" runID [1] = 5 (execution ID number) - cghg vs simtel"<<endl
      	<<"       [2] - in root file list"<<endl
	<<"       [3] - name of root file with histograms"<<endl
      	<<"       [4] - part type (proton_diff, gamma_diff, ele_pos)"<<endl
	<<"       [5] - max_trg_ev"<<endl;
    cout<<" runID [1] = 6 (execution ID number) - cghg vs simtel"<<endl
      	<<"       [2] - in root file"<<endl
	<<"       [3] - name of root file with histograms"<<endl
      	<<"       [4] - part type (gamma)"<<endl
	<<"       [5] - max_trg_ev"<<endl;

  }
  return 0;
}
