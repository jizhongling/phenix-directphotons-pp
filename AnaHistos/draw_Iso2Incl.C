#include "DivideFunctions.h"

void draw_Iso2Incl()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo.root");
  TH1 *h_photon = (TH1*)f->Get("h_photon_eta025");
  TH1 *h_isophoton = (TH1*)f->Get("h_isophoton_eta025");

  mc();
  mcd();

  TGraphErrors *gr = DivideHisto(h_isophoton, h_photon);
  aset(gr, "p_{T} [GeV]","isophoton/inclphoton", 5.,30.);
  style(gr, 20, 1);
  gr->Draw("AP");

  c0->Print("plots/Iso2Incl.pdf");
}
