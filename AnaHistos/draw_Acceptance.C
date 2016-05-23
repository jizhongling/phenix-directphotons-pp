void draw_Acceptance()
{
  TFile *fsig = new TFile("Acceptance-signal.root");
  TH1 *h_phsig = (TH1*)fsig->Get("h_photon");

  TFile *ftot = new TFile("Acceptance-total.root");
  TH1 *h_phtot = (TH1*)ftot->Get("h_photon");

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_phsig, h_phtot);
  gr->GetXaxis()->SetRangeUser(0., 30.);
  gr->GetYaxis()->SetRangeUser(0., 0.1);
  gr->Draw("AP");

  c->Print("Acceptance.pdf");
}
