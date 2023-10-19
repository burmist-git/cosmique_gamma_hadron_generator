//root
#include "TROOT.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFree.h"
#include "TBranch.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <vector>

using namespace std;

////
const Double_t earthR = 6371;                             // km
const Double_t earthD = earthR*2;                         // km
const Double_t generationH = 1;                           // km
//const Double_t generationH = 200;                           // km
const Double_t generationH_from_c = earthR + generationH; // km
const Double_t track_speed = TMath::C()/1000;             // km/s

//
const Double_t cos_min_generation = TMath::Cos(3.0/180.0*TMath::Pi());
//const Double_t cos_min_generation = TMath::Cos(15.0/180.0*TMath::Pi());

const Double_t x0_LST01 = -70.93/1000.0;                  //km
const Double_t y0_LST01 = -52.07/1000.0;                  //km
const Double_t z0_LST01 = earthR;                         //km
////
const TVector3 LST_01_r0(x0_LST01,y0_LST01,z0_LST01);
const TVector3 LST_01_n0(0.0,TMath::Sin(20.0/180.0*TMath::Pi()),TMath::Cos(20.0/180.0*TMath::Pi()));
////
void genTrk(Double_t &x0, Double_t &y0, Double_t &z0, Double_t &vx, Double_t &vy, Double_t &vz, Double_t &theta, Double_t &phi, TRandom3 *rnd);
Int_t getIntersection(Double_t &x1_int, Double_t &y1_int, Double_t &z1_int, Double_t &t1_int,
		      Double_t &x2_int, Double_t &y2_int, Double_t &z2_int, Double_t &t2_int,
		      Double_t x0, Double_t y0, Double_t z0, Double_t vx, Double_t vy, Double_t vz);
Bool_t getLSTArea(Double_t x1_int, Double_t y1_int, Double_t ellipse_A, Double_t ellipse_B);
////
//Double_t getDistToEarth_andPosition(Double_t x0, Double_t y0, Double_t z0, Double_t vx, Double_t vy, Double_t vz, Double_t &xe0, Double_t &ye0, Double_t &ze0);
//Double_t getDistToTerzinaOfIntersection(Double_t x_int, Double_t y_int, Double_t z_int);
//Double_t getAngleBetweenTrzinaAndTrk(Double_t vx, Double_t vy, Double_t vz);
////
int main(int argc, char *argv[]){
  if(argc == 6 && atoi(argv[1])==0){
    //
    Int_t nEvents = atoi(argv[2]);
    TString outputRootFile = argv[3];
    Int_t randomSeed = atoi(argv[4]);
    Int_t statisticsMultiplyFactor = atoi(argv[5]);
    //
    cout<<"nEvents                  "<<nEvents<<endl
	<<"outputRootFile           "<<outputRootFile<<endl
      	<<"randomSeed               "<<randomSeed<<endl
      	<<"statisticsMultiplyFactor "<<statisticsMultiplyFactor<<endl;
    //
    Double_t ellipse_A = 1732.1/1000.0;
    Double_t ellipse_B = 1500.0/1000.0;
    //
    TRandom3 *rnd = new TRandom3(randomSeed);
    //
    Double_t x0, y0, z0;
    Double_t vx, vy, vz;
    Double_t theta;
    Double_t phi;
    Double_t x1_int, y1_int, z1_int, t1_int;
    Double_t x2_int, y2_int, z2_int, t2_int;
    //
    TFile *hfile = new TFile(outputRootFile.Data(), "RECREATE", "cosmique proton generator", 1);
    if (hfile->IsZombie()) {
      cout << "PROBLEM with the initialization of the output ROOT ntuple " 
	   << outputRootFile << ": check that the path is correct!!!"
	   << endl;
      assert(0);
    }
    TTree *tree = new TTree("T", "cosmique proton generator");
    hfile->SetCompressionLevel(2);
    tree->SetAutoSave(1000000);
    // Create new event
    TTree::SetBranchStyle(0);
    //
    tree->Branch("x0",&x0, "x0/D");
    tree->Branch("y0",&y0, "y0/D");
    tree->Branch("z0",&z0, "z0/D");
    //
    tree->Branch("vx",&vx, "vx/D");
    tree->Branch("vy",&vy, "vy/D");
    tree->Branch("vz",&vz, "vz/D");
    //
    tree->Branch("theta",&theta, "theta/D");
    tree->Branch("phi",&phi, "phi/D");
    //
    tree->Branch("x1_int",&x1_int, "x1_int/D");
    tree->Branch("y1_int",&y1_int, "y1_int/D");
    tree->Branch("z1_int",&z1_int, "z1_int/D");
    tree->Branch("t1_int",&t1_int, "t1_int/D");
    //
    tree->Branch("x2_int",&x2_int, "x2_int/D");
    tree->Branch("y2_int",&y2_int, "y2_int/D");
    tree->Branch("z2_int",&z2_int, "z2_int/D");
    tree->Branch("t2_int",&t2_int, "t2_int/D");
    //
    for(Int_t j = 0;j<statisticsMultiplyFactor;j++){
      for(Int_t i = 0;i<nEvents;i++){
	genTrk( x0, y0, z0, vx, vy, vz, theta, phi,rnd);
	if(getIntersection( x1_int, y1_int, z1_int, t1_int,
			    x2_int, y2_int, z2_int, t2_int,
			    x0, y0, z0, vx, vy, vz)>1){
	  if(t1_int>0.0 && t2_int > 0.0){
	    if(getLSTArea(x1_int, y1_int, ellipse_A, ellipse_B)){
	      tree->Fill();
	    }
	  }
	}
      }
    }
    //
    hfile = tree->GetCurrentFile();
    hfile->Write();
    hfile->Close();
  }
  else{
    //------------------------------------------------
    cout<<" --> ERROR in input arguments "<<endl
	<<" runID [1] = 0 (execution ID number)"<<endl
      	<<"       [2] - n events"<<endl
	<<"       [3] - output root file"<<endl
	<<"       [4] - random seed"<<endl
	<<"       [5] - statistics multiply factor"<<endl;
  }
  return 0;
}

/*
Double_t getDistToTerzinaOfIntersection(Double_t x_int, Double_t y_int, Double_t z_int){
  TVector3 v_int(x_int,y_int,z_int);
  TVector3 v = terzina_r0 - v_int;
  return v.Mag();
}
*/

void genTrk(Double_t &x0, Double_t &y0, Double_t &z0, Double_t &vx, Double_t &vy, Double_t &vz, Double_t &theta, Double_t &phi, TRandom3 *rnd){
  //
  Double_t theta_r0 = TMath::ACos(rnd->Uniform(cos_min_generation,1.0));
  Double_t phi_r0   = rnd->Uniform(0.0,TMath::TwoPi());
  theta  = TMath::ACos(rnd->Uniform(-1.0,1.0));
  phi    = rnd->Uniform(0.0,TMath::TwoPi());
  TVector3 v_r0;
  TVector3 v_v;
  v_r0.SetMagThetaPhi(generationH_from_c,theta_r0,phi_r0);
  v_v.SetMagThetaPhi(track_speed,theta,phi);
  x0 = v_r0.x();
  y0 = v_r0.y();
  z0 = v_r0.z();
  vx = v_v.x();
  vy = v_v.y();
  vz = v_v.z();
}

Bool_t getLSTArea(Double_t x1_int, Double_t y1_int, Double_t ellipse_A, Double_t ellipse_B){
  //h2_ycore_vs_xcore_ring->Fill(,_anaConf.*TMath::Sin(angleSim));
  TVector2 v_v(x1_int, y1_int);
  Double_t xe_max = ellipse_A*TMath::Cos(v_v.Phi());
  Double_t ye_max = ellipse_B*TMath::Sin(v_v.Phi());
  Double_t re_max = TMath::Sqrt(xe_max*xe_max + ye_max*ye_max); 
  if(re_max>=v_v.Mod())
    return true;
  return false;
}
  

Int_t getIntersection(Double_t &x1_int, Double_t &y1_int, Double_t &z1_int, Double_t &t1_int,
		      Double_t &x2_int, Double_t &y2_int, Double_t &z2_int, Double_t &t2_int,
		      Double_t x0, Double_t y0, Double_t z0, Double_t vx, Double_t vy, Double_t vz){
  //
  x1_int = -999.0;
  y1_int = -999.0;
  z1_int = -999.0;
  t1_int = -999.0;
  x2_int = -999.0;
  y2_int = -999.0;
  z2_int = -999.0;
  t2_int = -999.0;
  //
  TVector3 v_r0(x0,y0,z0);
  TVector3 v_v(vx,vy,vz);
  //  
  Double_t det2 = 4*(v_v.Dot(v_r0)*v_v.Dot(v_r0)-v_v.Mag2()*(v_r0.Mag2() - earthR*earthR));
  if(det2 == 0.0){
    Double_t t = -v_v.Dot(v_r0)/v_v.Mag2();
    t1_int = t;
    TVector3 v;
    v = v_r0 + v_v*t;
    x1_int = v.x();
    y1_int = v.y();
    z1_int = v.z();
    return 1;
  }
  else if(det2 < 0.0){
    return 0;
  }
  else{
    Double_t det = TMath::Sqrt(det2);
    Double_t t1 = (-2*v_v.Dot(v_r0) - det)/2/v_v.Mag2();
    Double_t t2 = (-2*v_v.Dot(v_r0) + det)/2/v_v.Mag2();
    if(t1<t2){
      t1_int = t1;
      t2_int = t2;
    }
    else{
      t1_int = t2;
      t2_int = t1;
    }
    //
    TVector3 v1;
    TVector3 v2;
    //
    v1 = v_r0 + v_v*t1_int;
    x1_int = v1.x();
    y1_int = v1.y();
    z1_int = v1.z();
    //
    v2 = v_r0 + v_v*t2_int;
    x2_int = v2.x();
    y2_int = v2.y();
    z2_int = v2.z();
    return 2;
  }
}

/*
Double_t getDistToEarth_andPosition(Double_t x0, Double_t y0, Double_t z0, Double_t vx, Double_t vy, Double_t vz, Double_t &xe0, Double_t &ye0, Double_t &ze0){
  //
  TVector3 v_r0(x0,y0,z0);
  TVector3 v_v(vx,vy,vz);
  //
  Double_t div = v_v.Dot(v_v);
  if(div == 0.0)
    assert(0);
  Double_t t = -v_r0.Dot(v_v)/div;
  TVector3 v;
  v = v_r0 + v_v*t;
  //
  xe0 = v.x();
  ye0 = v.y();
  ze0 = v.z();
  //
  return (v.Mag()-earthR);
}

Double_t getAngleBetweenTrzinaAndTrk(Double_t vx, Double_t vy, Double_t vz){
  TVector3 v_v(vx,vy,vz);
  return TMath::ACos(v_v.Dot(terzina_n0)/v_v.Mag()/terzina_n0.Mag());
}
*/
