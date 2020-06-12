#include "DivideFunctions.h"

void draw_ChargedPionRatio()
{
  const char *pname[2] = {"PbSc", "PbGl"};
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  TFile *f_pythia = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-chargedpion.root");
  THnSparse *hn_photon = (THnSparse*)f_pythia->Get("hn_photon");
  hn_photon->GetAxis(4)->SetRange(3,3);

  mc();
  mcd();
  
  for(int iso=0; iso<2; iso++)
    for(int part=0; part<1; part++)
    {
      hn_photon->GetAxis(3)->SetRange(1+iso,2);
      hn_photon->GetAxis(2)->SetRange(secl[part],sech[part]);
      hn_photon->GetAxis(5)->SetRange(5,5);
      TH1 *h_pi0 = hn_photon->Projection(1);
      h_pi0->SetName("h_pi0");
      h_pi0->Scale(1./3.);  // filled each ival 3 times
      hn_photon->GetAxis(5)->SetRange(6,6);
      TH1 *h_chpi = hn_photon->Projection(1);
      h_chpi->SetName("h_chpi");

      TGraphErrors *gr_ratio = DivideHisto(h_chpi, h_pi0);
      aset(gr_ratio, "p_{T} [GeV]","#pi^{#pm}/#pi^{0}", 5.,30., 0.,0.2);
      style(gr_ratio, 20+4*iso+part, 1+part);
      gr_ratio->Draw(iso==0&&part==0?"AP":"P");

      delete h_pi0;
      delete h_chpi;
    }

  c0->Print("plots/ChargedPionRatio.pdf");
}
