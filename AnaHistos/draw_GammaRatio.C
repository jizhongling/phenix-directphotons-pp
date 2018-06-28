#include "DivideFunctions.h"

void draw_GammaRatio()
{
  TFile *f = new TFile("data/GammaRatio-histo.root");
  TH3 *h3_particles = (TH3*)f->Get("h3_particles");
  h3_particles->GetZaxis()->SetRange(1,7);
  TH1* h_others = h3_particles->ProjectionX("h_others", 3,5);
  TH1* h_pion = h3_particles->ProjectionX("h_pion", 2,2);

  mc();
  mcd();

  TGraphErrors *gr = DivideHisto(h_others, h_pion);
  gr->SetTitle("Photon Ratio");
  aset(gr, "p_{T} [GeV]","#frac{#eta+#omega+#eta'}{#pi^{0}}", 2.,30., 0.,0.5);
  style(gr, 20, 1);
  gr->Draw("AP");

  c0->Print("plots/GammaRatio.pdf");
}
