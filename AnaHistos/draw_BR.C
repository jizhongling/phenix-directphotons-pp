void draw_BR()
{
  TFile *f = new TFile("particles.root");
  TH2 *h2_particles = (TH2*)f->Get("h2_particles");
  TH1* h_others = h2_particles->ProjectionX("h_others", 2, 3);
  TH1* h_pion = h2_particles->ProjectionX("h_pion", 1, 1);

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_others, h_pion);
  gr->SetTitle("Branching ratio");
  gr->GetXaxis()->SetTitle("p_{T} [GeV]");
  gr->GetYaxis()->SetTitle("#frac{#eta+#omega}{#pi^{0}}");
  gr->GetXaxis()->SetRangeUser(0., 30.);
  gr->GetYaxis()->SetRangeUser(0., 0.5);
  gr->Draw("AP");

  c->Print("BR.pdf");

  Int_t n = gr->GetN();
  Double_t *x = gr->GetX();
  Double_t *y = gr->GetY();
  for(Int_t i=0; i<n; i++)
    cout << x[i] << ",";
  cout << endl;
  for(Int_t i=0; i<n; i++)
    cout << y[i] << ",";
  cout << endl;
  Double_t BRbar = h_others->Integral(11,30) / h_pion->Integral(11,30);
  cout << "BRbar = " << BRbar << endl;
}
