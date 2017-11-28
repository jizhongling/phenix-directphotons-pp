#include "GlobalVars.h"
#include "FitMinv.h"
#include "GetEfficiency.h"

void draw_ProbEff()
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  TGraphAsymmErrors *gr[2];
  Int_t igp[2] = {};
  for(Int_t part=0; part<2; part++)
  {
    gr[part] = new TGraphAsymmErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
    for(Int_t ic=0; ic<2; ic++)
      mc(part*2+ic, 6,5);
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  TAxis *axis_sec = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);
  TAxis *axis_cut = hn_pion->GetAxis(3);
  TAxis *axis_type = hn_pion->GetAxis(4);

  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      axis_type->SetRange(3,3);
      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);
      TH1 *h_minv;

      Double_t nt, ent;
      mcd(part*2, ipt+1);
      axis_cut->SetRange(2,2);
      h_minv = hn_pion->Projection(2);
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, nt, ent);
      if( ipt > 22 )
        nt = h_minv->Integral(12,16);
      delete h_minv;

      Double_t np, enp;
      mcd(part*2+1, ipt+1);
      axis_cut->SetRange(4,4);
      h_minv = hn_pion->Projection(2);
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, np, enp);
      if( ipt > 22 )
        np = h_minv->Integral(12,16);
      delete h_minv;

      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t yy, eyyl, eyyh;
      if( !GetEfficiency(nt,np, yy,eyyl,eyyh) )
      {
        eyyl = yy * sqrt( pow(ent/nt,2.) + pow(enp/np,2.) );
        eyyh = 0.;
      }
      if( yy >= 0. && eyyl >= 0. && eyyl < TMath::Infinity() )
      {
        gr[part]->SetPoint(igp[part], xx, yy);
        gr[part]->SetPointError(igp[part], 0.,0., eyyl,eyyh);
        igp[part]++;
      }
    }

  TFile *f_sasha = new TFile("data/sasha-prob.root");

  mc(5, 2,1);

  gr[0]->SetTitle("Prob Eff for PbSc");
  gr[1]->SetTitle("Prob Eff for PbGl");

  for(Int_t part=0; part<2; part++)
  {
    mcd(5, part+1);
    gr[part]->Set(igp[part]);
    aset(gr[part], "p_{T} [GeV]","Prob Eff", 0.,30., 0.8,1.1);
    style(gr[part], part+20, part+1);
    gr[part]->Draw("AP");
    TGraphErrors *gr_sasha = (TGraphErrors*)f_sasha->Get( Form("gr_%d",part) );
    gr_sasha->Draw("C");
  }

  TFile *f_out = new TFile("data/ProbEff.root", "RECREATE");
  for(Int_t part=0; part<2; part++)
  {
    gr[part]->Write();
    mcw( part*2, Form("part%d-total",part) );
    mcw( part*2+1, Form("part%d-passed",part) );
  }
  f_out->Close();

  c5->Print("plots/ProbEff.pdf");
}
