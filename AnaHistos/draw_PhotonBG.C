#include "DivideFunctions.h"

void draw_PhotonBG()
{
  TFile *f_pythia = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-minbias.root");

  THnSparse *hn_pion = (THnSparse*)f_pythia->Get("hn_pion");
  hn_pion->GetAxis(7)->SetRange(1,1);  // isys
  hn_pion->GetAxis(6)->SetRange(3,3);  // ival
  TH1 *h_pi0 = hn_pion->Projection(1);  // reco pT
  h_pi0->SetName("h_pi0");

  THnSparse *hn_hadron = (THnSparse*)f_pythia->Get("hn_hadron");
  hn_hadron->GetAxis(4)->SetRange(2,2);  // InAcc
  hn_hadron->GetAxis(3)->SetRange(6,6);  // other hadron
  TH1 *h_hadron = hn_hadron->Projection(1);  // reco pT
  h_hadron->SetName("h_hadron");

  mc();
  mcd();

  TGraphErrors *gr_hadron = DivideHisto(h_hadron, h_pi0);
  aset(gr_hadron, "p_{T} [GeV/c]","#gamma_{BG}/#pi^{0}", 5.,30., 0.,5e-3);
  style(gr_hadron, 20, 1);
  gr_hadron->Draw("AP");

  c0->Print("plots/PhotonBG.pdf");
}
