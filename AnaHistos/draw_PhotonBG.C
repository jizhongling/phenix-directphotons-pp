#include "DivideFunctions.h"

void draw_PhotonBG()
{
  TFile *f_pythia = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-chargedpion.root");

  THnSparse *hn_pion = (THnSparse*)f_pythia->Get("hn_pion");
  hn_pion->GetAxis(7)->SetRange(1,1);
  hn_pion->GetAxis(6)->SetRange(3,3);
  TH1 *h_pi0 = hn_pion->Projection(1);
  h_pi0->SetName("h_pi0");

  THnSparse *hn_photon = (THnSparse*)f_pythia->Get("hn_photon");
  hn_photon->GetAxis(5)->SetRange(1,1);
  hn_photon->GetAxis(4)->SetRange(3,3);
  TH1 *h_photon = hn_photon->Projection(1);
  h_photon->SetName("h_photon");

  THnSparse *hn_hadron = (THnSparse*)f_pythia->Get("hn_hadron");
  hn_hadron->GetAxis(3)->SetRange(2,2);
  hn_hadron->GetAxis(4)->SetRange(3,6);
  TH1 *h_hadron = hn_hadron->Projection(1);
  h_hadron->SetName("h_hadron");

  mc();
  mcd();

  TGraphErrors *gr_photon = DivideHisto(h_photon, h_pi0);
  TGraphErrors *gr_hadron = DivideHisto(h_hadron, h_pi0);
  aset(gr_photon, "p_{T} [GeV]","#gamma/#pi^{0}", 5.,30., 1.,3.2);
  style(gr_photon, 20, 1);
  style(gr_hadron, 21, 2);
  gr_photon->Draw("AP");
  gr_hadron->Draw("P");

  c0->Print("plots/PhotonBG.pdf");
}
