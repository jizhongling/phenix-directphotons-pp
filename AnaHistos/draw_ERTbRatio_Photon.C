#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"

void draw_ERTbRatio_Photon()
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  QueryTree *qt_ratio = new QueryTree("data/ERTbRatio-photon.root", "RECREATE");

  QueryTree *qt_misscorr = new QueryTree("data/MissCorr.root");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  TAxis *axis_1sec = hn_1photon->GetAxis(0);
  TAxis *axis_1pt = hn_1photon->GetAxis(1);
  TAxis *axis_1cut = hn_1photon->GetAxis(3);
  TAxis *axis_1type = hn_1photon->GetAxis(4);

  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");
  TAxis *axis_sec = hn_2photon->GetAxis(0);
  TAxis *axis_pt = hn_2photon->GetAxis(1);
  TAxis *axis_minv = hn_2photon->GetAxis(2);
  TAxis *axis_cut = hn_2photon->GetAxis(4);
  TAxis *axis_type = hn_2photon->GetAxis(5);

  for(int part=0; part<3; part++)
  {
    mc(part, 4,3);
    mc(part+3, 4,3);
  }

  for(int part=0; part<3; part++)
    for(int ipt=20; ipt<30; ipt++)
    {
      double nphoton_ertc, npion_ertc, enpion_ertc;
      double nphoton_ertb, npion_ertb, enpion_ertb;
      TH1 *h_photon, *h_minv;

      axis_1pt->SetRange(ipt+1,ipt+1);
      axis_1cut->SetRange(4,4);

      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);
      axis_cut->SetRange(4,4);

      axis_1type->SetRange(3,3);
      h_photon = hn_1photon->Projection(0);
      nphoton_ertc = h_photon->Integral(secl[part],sech[part]);
      delete h_photon;
      mcd(part, ipt+1);
      axis_type->SetRange(3,3);
      h_minv = hn_2photon->Projection(2);
      h_minv->Rebin(10);
      h_minv->Scale(0.5);
      h_minv->SetTitle( Form("ERT4x4c Part %d",part) );
      FitMinv(h_minv, npion_ertc, enpion_ertc, true, 0.10,0.17);
      delete h_minv;

      axis_1type->SetRange(2,2);
      h_photon = hn_1photon->Projection(0);
      nphoton_ertb = h_photon->Integral(secl[part],sech[part]);
      delete h_photon;
      mcd(part+3, ipt+1);
      axis_type->SetRange(2,2);
      h_minv = hn_2photon->Projection(2);
      h_minv->Scale(0.5);
      h_minv->SetTitle( Form("ERT4x4b Part %d",part) );
      FitMinv(h_minv, npion_ertb, enpion_ertb, true, 0.10,0.17);
      delete h_minv;

      double xpt, MissCorr, eMissCorr;
      qt_misscorr->Query(ipt, part, xpt, MissCorr, eMissCorr);
      MissCorr = 1.;

      double ndir_ertc = nphoton_ertc - MissCorr*npion_ertc*2.;
      double endir_ertc = sqrt( nphoton_ertc + pow(eMissCorr*npion_ertc*2.,2) + pow(MissCorr*enpion_ertc*2.,2) );
      double ndir_ertb = nphoton_ertb - MissCorr*npion_ertb*2.;
      double endir_ertb = sqrt( nphoton_ertb + pow(eMissCorr*npion_ertb*2.,2) + pow(MissCorr*enpion_ertb*2.,2) );
      double ratio = ndir_ertc / ndir_ertb;
      double eratio = ratio * sqrt( pow(endir_ertc/ndir_ertc,2) + pow(endir_ertb/ndir_ertb,2) );

      if( TMath::Finit(ratio+eratio) && eratio > 0.)
        qt_ratio->Fill(ipt, part, xpt, ratio, eratio);
    }

  mc(6, 3,1);

  for(int part=0; part<3; part++)
  {
    TGraphErrors *gr = qt_ratio->Graph(part);
    mcd(6, part+1);
    gr->Set(igp[part]);
    aset(gr, "p_{T} [GeV]");
    style(gr, part+20, part+1);
    gr->Draw("AP");
    gr->Fit("pol0", "Q","", 10.,30.);
  }

  c6->Print("plots/ERTbRatio-photon.pdf");

  qt_ratio->Write();
  for(int part=0; part<6; part++)
    mcw( part, Form("ERT%c-part%d",99-part/3,part%3) );
  qt_ratio->Close();
}
