void draw_BR()
{
  TFile *f = new TFile("particles.root");
  TH2 *h2_particles = (TH2*)f->Get("h2_particles");
  TH1* h_others = h2_particles->ProjectionX("h_others", 2, 3);
  TH1* h_pion = h2_particles->ProjectionX("h_pion", 1, 1);

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_others, h_pion);
  gr->GetXaxis()->SetRangeUser(0., 30.);
  gr->GetYaxis()->SetRangeUser(0., 0.5);
  gr->Draw("AP");

  c->Print("BR.pdf");

  Double_t RatioAbar = h_others->Integral(11,30) / h_pion->Integral(11,30);
  cout << "RatioAbar=" << RatioAbar << endl;
}
