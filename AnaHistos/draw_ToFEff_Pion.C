#include "GlobalVars.h"
#include "FitMinv.h"
#include "GetEfficiency.h"

void draw_ToFEff_Pion()
{
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};
  const char *name[2] = {"PbSc", "PbGl"};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  TAxis *axis_sec = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);
  TAxis *axis_cut = hn_pion->GetAxis(3);
  TAxis *axis_type = hn_pion->GetAxis(4);

  TGraphAsymmErrors *gr[2];
  int igr[2] = {};
  for(int part=0; part<2; part++)
    gr[part] =  new TGraphAsymmErrors(npT);

  mc(0, 6,5);
  mc(1, 6,5);
  mc(2, 2,1);

  for(int part=0; part<2; part++)
    for(int ipt=0; ipt<npT; ipt++)
    {
      axis_type->SetRange(3,3);
      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);
      TH1 *h_minv;

      mcd(0, ipt+1);
      double nt, ent;
      axis_cut->SetRange(3,3);
      h_minv = hn_pion->Projection(2); 
      h_minv->Rebin(10);
      h_minv->SetTitle(Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]));
      FitMinv(h_minv, nt, ent);
      delete h_minv;

      mcd(1, ipt+1);
      double np, enp;
      axis_cut->SetRange(4,4);
      h_minv = hn_pion->Projection(2); 
      h_minv->Rebin(10);
      h_minv->SetTitle(Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]));
      FitMinv(h_minv, np, enp);
      delete h_minv;

      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double yy, eyyl, eyyh;
      if( !GetEfficiency(nt,np, yy,eyyl,eyyh) )
      {
        eyyl = yy * sqrt( pow(ent/nt,2.) + pow(enp/np,2.) );
        eyyh = 0.;
      }
      if( yy >= 0. && eyyl >= 0. && eyyl < TMath::Infinity() )
      {
        gr[part]->SetPoint(igr[part], xx, yy);
        gr[part]->SetPointError(igr[part], 0.,0., eyyl,eyyh);
        igr[part]++;
      }
    }

  for(int part=0; part<2; part++)
  {
    mcd(2, part+1);
    gr[part]->SetTitle( Form("ToF efficeincy for %s", name[part]) );
    aset(gr[part], "p_{T} [GeV]", "Eff", 0.,30., 0.,1.1);
    style(gr[part], 24, kRed);
    gr[part]->Draw("APE");
    gr[part]->Fit("pol0", "Q","", 9.,30.);

    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr[part]->FindObject("stats");
    st->SetY1NDC(0.6);
    st->SetY2NDC(0.8);
  }

  TFile *f_out = new TFile("data/ToFEff-fit.root", "RECREATE");
  mcw(0, "PbSc.pdf");
  mcw(1, "PbGl.pdf");
  f_out->Close();

  c2->Print("plots/ToFEff-pion.pdf");
}
