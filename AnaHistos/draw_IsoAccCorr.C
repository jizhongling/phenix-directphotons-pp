#include "DivideFunctions.h"

void draw_IsoAccCorr()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-photon.root");
  THnSparse *hn_photon = (THnSparse*)f->Get("hn_hadron");
  hn_photon->GetAxis(3)->SetRange(1,2);  // prompt photon
  hn_photon->GetAxis(4)->SetRange(2,2);  // InAcc

  hn_photon->GetAxis(6)->SetRange(1,2);  // econe_acc[2]
  hn_photon->GetAxis(5)->SetRange(2,2);  // econe_all
  TH1 *h_isoall = hn_photon->Projection(0);  // pt_truth
  hn_photon->GetAxis(5)->SetRange(1,2);  // econe_all
  hn_photon->GetAxis(6)->SetRange(2,2);  // econe_acc[2]
  TH1 *h_isoacc = hn_photon->Projection(0);  // pt_truth

  mc();
  mcd();
  TGraphErrors *gr = DivideHisto(h_isoall, h_isoacc);
  aset(gr, "p_{T} (GeV/c)","IsoAll/IsoAcc", 5.,30., 0.8,1.);
  style(gr, 20, 1);
  gr->Draw("AP");
  c0->Print("plots/IsoAll2IsoAcc.pdf");
}
