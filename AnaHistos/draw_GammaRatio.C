#include "DivideFunctions.h"

void draw_GammaRatio()
{
  TFile *f = new TFile("data/GammaRatio-histo.root");
  TH2 *h2_particles = (TH2*)f->Get("h2_particles");
  TH1* h_others = h2_particles->ProjectionX("h_others", 2, 3);
  TH1* h_pion = h2_particles->ProjectionX("h_pion", 1, 1);

  mc();
  mcd();

  TGraphErrors *gr = DivideHisto(h_others, h_pion);
  gr->SetTitle("Photon Ratio");
  aset(gr, "p_{T} [GeV]","#frac{#eta+#omega}{#pi^{0}}", 0.,12., 0.,0.5);
  style(gr, 20, 1);
  gr->Draw("AP");

  c0->Print("plots/GammaRatio.pdf");
}
