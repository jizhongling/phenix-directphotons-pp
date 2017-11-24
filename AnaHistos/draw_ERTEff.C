#include "GlobalVars.h"
#include "FitMinv.h"
#include "GetEfficiency.h"

void draw_ERTEff()
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  TGraphAsymmErrors *gr[2];
  Int_t igp[2] = {};
  for(Int_t part=0; part<2; part++)
  {
    gr[part] = new TGraphAsymmErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
    mc(part, 6,5);
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  THnSparse *hn_trig = (THnSparse*)f->Get("hn_ert_pion");
  TAxis *axis_sec = hn_trig->GetAxis(0);
  TAxis *axis_pt = hn_trig->GetAxis(1);
  TAxis *axis_minv = hn_trig->GetAxis(2);
  TAxis *axis_cond = hn_trig->GetAxis(3);

  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);
      TH1 *h_minv;

      Double_t nt, ent;
      mcd(part, ipt+1);
      axis_cond->SetRange(1,1);
      h_minv = hn_trig->Projection(2);
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, nt, ent);
      delete h_minv;

      Double_t np, enp;
      mcd(part, ipt+1);
      axis_cond->SetRange(4,4);
      h_minv = hn_trig->Projection(2);
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, np, enp);
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

  mc(2);
  mcd(2);

  for(Int_t part=0; part<2; part++)
  {
    gr[part]->SetName(Form("gr_%d",part));
    gr[part]->SetTitle("ERT_4x4c trigger efficeincy for #pi^{0}");
    aset(gr[part], "p_{T} [GeV]","Eff", 0.,30., 0.,1.1);
    style(gr[part], part+20, part+1);
    if(part==0)
      gr[part]->Draw("APE");
    else
      gr[part]->Draw("PE");
    TGraph *gr_sasha =  new TGraph( Form("sasha-trig-part%d.txt",part) );
    gr_sasha->Draw("C");
  }

  legi(0, 0.7,0.2,0.9,0.4);
  leg0->AddEntry(gr[0], "PbSc", "LPE");
  leg0->AddEntry(gr[1], "PbGl", "LPE");
  leg0->Draw();

  c2->Print("ERTEff.pdf");

  TFile *f_out = new TFile("data/ERTEff.root", "RECREATE");
  for(Int_t part=0; part<2; part++)
    gr[part]->Write();
  f_out->Close();
}
