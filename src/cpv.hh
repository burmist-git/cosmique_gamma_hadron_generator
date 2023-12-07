#ifndef cpv_hh
#define cpv_hh

//My
#include "cpvbase.hh"

//root
#include <TROOT.h>
#include <TMath.h>

class TChain;
class TFile;
class TTree;
class TString;
class TBranch;

static const Double_t earthR = 6371;                           // km
static const Double_t earthD = earthR*2;                       // km
static const Double_t satelliteH = 525;                        // km
static const Double_t satelliteH_from_c = earthR + satelliteH; // km
static const Double_t track_speed = TMath::C()/1000;           // km/s

class cpv: public cpvbase {
public:
  cpv() : cpvbase()
  {
  }

  cpv(TString fileList) : cpvbase(fileList)
  {
  }

  cpv(TString file, Int_t key) : cpvbase(file, key)
  {
  }

  void Loop(TString histOut,
	    Double_t ellipse_A_r_km, Double_t ellipse_B_r_km, Double_t theta_p_t_deg_max,
	    Double_t y0_LST01_km, Double_t x0_LST01_km,
	    Int_t max_trg_ev,
	    Double_t norm_n_simtel, Double_t norm_theta_p_t_deg_max, Double_t norm_Energy_min, Double_t norm_Energy_max, Double_t norm_r_core_max,
	    Double_t norm_simtel_Energy_min, Double_t norm_simtel_Energy_max, TString evH_E2_data_out);

  void Loop(TString histOut,
	    Double_t y0_LST01_km, Double_t x0_LST01_km,
	    Int_t max_trg_ev, Double_t norm_theta_p_t_deg_max, Double_t norm_r_core_max);
  
  void Loop(TString histOut);
  void createEASCherSim_ini( Int_t eventIDmy, Double_t theta_deg, Double_t heightAboveEarth_km);
  void create_cosmique_proton_generator_info(Int_t eventIDmy, Double_t theta_deg, Double_t phi_deg, Double_t heightAboveEarth_km);
  void readEASCherSim( TString inRootFileName, double dist,  double &nphotons_per_m2);
  void readEASCherSimNewFormat( TString inRootFileName, double dist,  double &nphotons_per_m2);
  TString getShowerRootFileName(Int_t eventIDmy);
  bool file_exists_test(TString inRootFileName);
  void printTrkInfo( Int_t eventIDmy, double nphotons_per_m2);
  void saveTrkInfoAndPhotonDencity( Int_t eventIDmy, Double_t nphotons_per_m2);
  //
  Double_t conv_phi(Double_t phiv);
  //
  Double_t get_tot_flux(Double_t e_min, Double_t e_max);
  Double_t get_tot_flux_gamma(Double_t e_min, Double_t e_max);
  Double_t get_surface_fluxs(Double_t r_in_m, Double_t theta);
  Double_t get_tot_flux_simulation(Double_t e_min, Double_t e_max);
  Double_t function_log_parabola_fit_MAGIC_tel_crab_GeV(Double_t e);
  Double_t get_tot_flux_MAGIC_tel_crab(Int_t nn=1000000, Double_t e_min=50, Double_t e_max=3000);
  Double_t get_diff_flux_ele_pos_power_law(Double_t e);
  Double_t get_tot_flux_ele_pos(Int_t nn=1000000, Double_t e_min=5, Double_t e_max=4000);
  //
  void calculate_rate(Double_t ellipse_min_r = 1500.0,
		      TString pdf_out = "flux_diff_protons.pdf",
		      TString dat_out = "flux_diff_protons.dat",
		      TString frame_title = "Proton flux, Hz",
		      Double_t min_f = 0.1,
		      Double_t max_f = 4*1.0e+8,
		      TString part_type = "proton_diff");
  void calculate_diff_proton_rate();
  void calculate_diff_gamma_rate();
  void calculate_diff_electron_rate();
  void calculate_crab_gamma_rate();
  //
};

#endif
