void draw_Acceptance()
{
  TFile *fsig = new TFile("Acceptance-signal.root");
  TH1 *h_phsig = (TH1*)fsig->Get("h_photon");

  TFile *ftot = new TFile("Acceptance-total.root");
  TH1 *h_phtot = (TH1*)ftot->Get("h_photon");

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_phsig, h_phtot);
  gr->SetTitle("Acceptance");
  gr->GetXaxis()->SetTitle("p_{T} [GeV]");
  gr->GetYaxis()->SetTitle("Acceptance");
  gr->GetXaxis()->SetRangeUser(0., 30.);
  gr->GetYaxis()->SetRangeUser(0., 0.3);
  gr->Draw("AP");

  c->Print("Acceptance.pdf");

  Double_t nsig = h_phsig->Integral(2,30);
  Double_t ntot = h_phtot->Integral(2,30);
  cout << "Acceptancebar = " << nsig/ntot << endl;
}
