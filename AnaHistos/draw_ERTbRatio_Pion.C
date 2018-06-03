#include "FitMinv.h"

void draw_ERTbRatio_Pion()
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  TAxis *axis_sec = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);
  TAxis *axis_cut = hn_pion->GetAxis(3);
  TAxis *axis_type = hn_pion->GetAxis(4);
  
  mc(0, 2,2);

  for(Int_t part=0; part<2; part++)
  {
    Double_t npion_ertc, enpion_ertc;
    Double_t npion_ertb, enpion_ertb;
    TH1 *h_minv;

    axis_sec->SetRange(secl[part],sech[part]);
    axis_pt->SetRange(21+part*2,30);
    axis_cut->SetRange(4,4);

    mcd(0, part+1);
    axis_type->SetRange(3,3);
    h_minv = hn_pion->Projection(2);
    h_minv->SetTitle( Form("ERT4x4c Part %d",part) );
    FitMinv(h_minv, npion_ertc, enpion_ertc, kTRUE, 0.10,0.17);
    delete h_minv;

    mcd(0, part+3);
    axis_type->SetRange(2,2);
    h_minv = hn_pion->Projection(2);
    h_minv->SetTitle( Form("ERT4x4b Part %d",part) );
    FitMinv(h_minv, npion_ertb, enpion_ertb, kTRUE, 0.10,0.17);
    delete h_minv;

    Double_t ratio = npion_ertc / npion_ertb;
    Double_t eratio = ratio * sqrt( pow(enpion_ertc/npion_ertc,2.) + pow(enpion_ertb/npion_ertb,2.) );

    cout << "Part " << part << ", ratio = " << ratio << ", eratio = " << eratio << endl;
  }

  c0->Print("plots/ERTbRatio-pion.pdf");
}
