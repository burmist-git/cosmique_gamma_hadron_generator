Int_t plots_cpv(){
  
  TString fileN;
  TString objectName;
  //fileN = "./hist.root";
  fileN = "./hist120km_cut_28427ev.root";
  TFile *f1 = new TFile(fileN.Data());

  TGraph2D *gr2D_sphere = (TGraph2D*)f1->Get("gr2D_sphere");
  TGraph2D *gr2D_int1 = (TGraph2D*)f1->Get("gr2D_int1");
  TGraph2D *gr2D_gen = (TGraph2D*)f1->Get("gr2D_gen");

  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800);
  gr2D_sphere->SetMarkerColor(kBlack);
  gr2D_sphere->SetMarkerStyle(7);
  gr2D_int1->SetMarkerColor(kRed);
  gr2D_int1->SetMarkerStyle(20);
  gr2D_gen->SetMarkerColor(kBlue);
  gr2D_gen->SetMarkerStyle(7);

  gr2D_sphere->SetTitle("");
  gr2D_sphere->Draw("P");
  gr2D_int1->Draw("sameP");
  gr2D_gen->Draw("sameP");

  gr2D_sphere->GetXaxis()->SetTitle("x");
  gr2D_sphere->GetYaxis()->SetTitle("y");
  gr2D_sphere->GetZaxis()->SetTitle("z");

  return 0;
}
