#include "GlobalVars.h"
#include "ReadGraph.h"
#include "FitMinv.h"

void draw_ERTbRatio_Photon()
{
  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};

  TGraphErrors *gr[3];
  Int_t igp[3] = {};
  for(Int_t part=0; part<3; part++)
  {
    gr[part] = new TGraphErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
  }

  Double_t xMiss[3][npT] = {}, Miss[3][npT] = {}, eMiss[3][npT] = {};
  for(Int_t part=0; part<3; part++)
  {
    ReadGraph<TGraphAsymmErrors>("data/HadronRatio.root", part, xMiss[part], Miss[part], eMiss[part]);
    mc(part, 4,3);
    mc(part+3, 4,3);
  }

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

  for(Int_t part=0; part<3; part++)
    for(Int_t ipt=20; ipt<30; ipt++)
    {
      Double_t nphoton_ertc, npion_ertc, enpion_ertc;
      Double_t nphoton_ertb, npion_ertb, enpion_ertb;
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
      FitMinv(h_minv, npion_ertc, enpion_ertc, kTRUE, 0.10,0.17);
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
      FitMinv(h_minv, npion_ertb, enpion_ertb, kTRUE, 0.10,0.17);
      delete h_minv;

      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Int_t ipMiss = Get_ipt(xMiss[part], xx);
      Double_t aMiss = 1.;//Miss[part][ipMiss];
      Double_t eaMiss = eMiss[part][ipMiss];

      Double_t ndir_ertc = nphoton_ertc - aMiss*npion_ertc*2.;
      Double_t endir_ertc = sqrt( nphoton_ertc + pow(eaMiss*npion_ertc*2.,2.) + pow(aMiss*enpion_ertc*2.,2.) );
      Double_t ndir_ertb = nphoton_ertb - aMiss*npion_ertb*2.;
      Double_t endir_ertb = sqrt( nphoton_ertb + pow(eaMiss*npion_ertb*2.,2.) + pow(aMiss*enpion_ertb*2.,2.) );
      Double_t ratio = ndir_ertc / ndir_ertb;
      Double_t eratio = ratio * sqrt( pow(endir_ertc/ndir_ertc,2.) + pow(endir_ertb/ndir_ertb,2.) );

      if(eratio > 0.)
      {
        gr[part]->SetPoint(igp[part], xx, ratio);
        gr[part]->SetPointError(igp[part], 0., eratio);
        igp[part]++;
      }
    }

  mc(6, 3,1);

  for(Int_t part=0; part<3; part++)
  {
    mcd(6, part+1);
    gr[part]->Set(igp[part]);
    aset(gr[part], "p_{T} [GeV]");
    style(gr[part], part+20, part+1);
    gr[part]->Draw("AP");
    gr[part]->Fit("pol0", "Q","", 10.,30.);
  }

  c6->Print("plots/ERTbRatio-photon.pdf");

  TFile *f_out = new TFile("data/ERTbRatio-photon.root", "RECREATE");
  for(Int_t part=0; part<6; part++)
  {
    mcw( part, Form("ERT%c-part%d",99-part/3,part%3) );
    if(part<3)
      gr[part]->Write();
  }
  f_out->Close();
}
