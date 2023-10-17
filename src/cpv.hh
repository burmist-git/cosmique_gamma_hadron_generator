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
  cpv(TString fileList) : cpvbase(fileList)
  {
  }

  cpv(TString file, Int_t key) : cpvbase(file, key)
  {
  }

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
  Double_t get_tot_flux(Double_t e_min, Double_t e_max);
  Double_t get_surface_fluxs(Double_t r_in_m, Double_t theta);
};

#endif
