#include "DivideFunctions.h"

void draw_ChargedPionRatio()
{
  TFile *f_pythia = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-minbias.root");

  THnSparse *hn_pion = (THnSparse*)f_pythia->Get("hn_pion");
  hn_pion->GetAxis(7)->SetRange(1,1);
  hn_pion->GetAxis(6)->SetRange(3,3);
  hn_pion->GetAxis(3)->SetRange(1,6);

  THnSparse *hn_photon = (THnSparse*)f_pythia->Get("hn_photon");
  hn_photon->GetAxis(5)->SetRange(5,5);
  hn_photon->GetAxis(4)->SetRange(3,3);
  hn_photon->GetAxis(2)->SetRange(1,6);

  mc();
  mcd();

  for(int iso=0; iso<2; iso++)
  {
    hn_pion->GetAxis(4)->SetRange(1+iso,2);
    TH1 *h_pi0 = hn_pion->Projection(1);
    h_pi0->SetName("h_pi0");

    hn_photon->GetAxis(3)->SetRange(1+iso,2);
    TH1 *h_chpi = hn_photon->Projection(1);
    h_chpi->SetName("h_chpi");

    TGraphErrors *gr_ratio = DivideHisto(h_chpi, h_pi0);
    aset(gr_ratio, "p_{T} [GeV/c]","#pi^{#pm}/#pi^{0}", 5.,30., 0.,0.2);
    style(gr_ratio, 20+iso, 1+iso);
    gr_ratio->Draw(iso==0?"AP":"P");

    delete h_pi0;
    delete h_chpi;
  }

  c0->Print("plots/ChargedPionRatio.pdf");
}
