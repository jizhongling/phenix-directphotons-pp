#include "DivideFunctions.h"

void draw_MissingRatio()
{
  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-warn-histo.root");
  THnSparse *hn_missing = (THnSparse*)f->Get("hn_missing");

  mc();
  mcd();

  for(Int_t part=0; part<3; part++)
  {
    hn_missing->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_missing->GetAxis(3)->SetRange(3,3);
    TH1 *h_2photon = hn_missing->Projection(1);
    hn_missing->GetAxis(3)->SetRange(2,2);
    TH1 *h_1photon = hn_missing->Projection(1);

    TGraphErrors *gr_miss = DivideHisto(h_1photon, h_2photon);
    aset(gr_miss, "p_{T} [GeV]","R", 0.,30., 0.,2.);
    style(gr_miss, 20+part, 1+part);
    if(part==0)
      gr_miss->Draw("AP");
    else
      gr_miss->Draw("P");

    delete h_2photon;
    delete h_1photon;
  }

  c0->Print("plots/MissingRatio.pdf");
}
