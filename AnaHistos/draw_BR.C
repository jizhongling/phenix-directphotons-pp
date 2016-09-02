void draw_BR()
{
  gROOT->ProcessLine(".L ReadGraph.C");

  TFile *f = new TFile("particles.root");
  TH2 *h2_particles = (TH2*)f->Get("h2_particles");
  TH1* h_others = h2_particles->ProjectionX("h_others", 2, 3);
  TH1* h_pion = h2_particles->ProjectionX("h_pion", 1, 1);

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_others, h_pion);
  gr->SetTitle("Two photons acceptance");
  gr->GetXaxis()->SetTitle("p_{T} [GeV]");
  gr->GetYaxis()->SetTitle("#frac{#eta+#omega}{#pi^{0}}");
  gr->GetXaxis()->SetRangeUser(0., 12.);
  gr->GetYaxis()->SetRangeUser(0., 0.5);
  gr->Draw("AP");

  c->Print("BR.pdf");

  Double_t gx[30], gy[30], egy[30];
  ReadGraph(gr, gx, gy, egy);

  Double_t nothers = h_others->Integral(11,30);
  Double_t npion = h_pion->Integral(11,30);
  Double_t BRbar = nothers / npion;
  Double_t eBRbar = BRbar * sqrt( 1./nothers + 1./npion ); 
  for(Int_t ipt=10; ipt<30; ipt++)
  {
    gy[ipt] = BRbar;
    egy[ipt] = eBRbar;
  }

  TGraphErrors *grout = new TGraphErrors(30, gx, gy, 0, egy);
  grout->SetName("gr_0");
  TFile *fout = new TFile("BR.root", "RECREATE");
  grout->Write();
  fout->Close();
}
