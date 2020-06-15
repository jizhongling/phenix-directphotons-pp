#include "DivideFunctions.h"

void draw_JetsRatio()
{
  TFile *f_pythia = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-GenPH-histo-jets.root");
  THnSparse *hn_photon = (THnSparse*)f_pythia->Get("hn_photon");
  hn_photon->GetAxis(4)->SetRange(3,3);

  TFile *f_pisa = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/HadronResponse-histo-jets.root");
  THnSparse *hn_1photon = (THnSparse*)f_pisa->Get("hn_1photon");
  hn_1photon->GetAxis(3)->SetRange(2,2);

  mc();
  mcd();

  const int ngroup = 2;
  for(int iso=0; iso<2; iso++)
  {
    hn_photon->GetAxis(3)->SetRange(1+iso,2);
    TH1 *h_photon = hn_photon->Projection(1);
    h_photon->SetName("h_photon");
    h_photon->Rebin(ngroup);

    hn_1photon->GetAxis(2)->SetRange(1+iso,2);
    TH1 *h_jets = hn_1photon->Projection(0);
    h_jets->SetName("h_jets");
    h_jets->Rebin(ngroup);

    TGraphErrors *gr_jets = DivideHisto(h_jets, h_photon);
    aset(gr_jets, "p_{T} [GeV]","N_{jets}/N_{#gamma}", 5.,30., 0.5,1.1);
    style(gr_jets, 20+iso, 1+iso);
    gr_jets->Draw(iso==0?"AP":"P");

    delete h_photon;
    delete h_jets;
  }

  c0->Print("plots/JetsRatio.pdf");
}
