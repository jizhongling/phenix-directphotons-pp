#include "GlobalVars.h"
#include "FitMinv.h"
#include "GetEfficiency.h"

void draw_Acceptance_Pion()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  SetWeight();

  TGraphAsymmErrors *gr[3];
  int igp[3] = {};
  for(int part=0; part<3; part++)
  {
    gr[part] = new TGraphAsymmErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
    mc(part, 6,5);
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo.root");

  TH1 *h_total = (TH1*)f->Get("h_pion");

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  TAxis *axis_pt0 = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);
  TAxis *axis_sec = hn_pion->GetAxis(3);
  TAxis *axis_peak = hn_pion->GetAxis(4);

  for(int part=0; part<3; part++)
    for(int ipt=0; ipt<npT; ipt++)
    {
      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double ww = cross_pi0->Eval(xx);

      TH1 *h_tt = (TH1*)h_total->Clone();
      h_tt->Scale(1./ww);
      double nt = h_tt->GetBinContent(ipt+1);
      double ent = sqrt(nt);
      delete h_tt;

      double np, enp;
      mcd(part, ipt+1);
      //axis_peak->SetRange(3,3);  // Require two peaks
      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);
      TH1 *h_minv = hn_pion->Projection(2);
      h_minv->Rebin(10);
      h_minv->Scale(1./ww);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, np, enp);
      // Do not use fit result, use all invariant mass instead
      np = h_minv->Integral(0,-1);
      enp = sqrt(np);
      delete h_minv;

      double yy, eyyl, eyyh;
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


  mc(3);
  mcd(3);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);

  for(int part=0; part<3; part++)
  {
    gr[part]->Set(igp[part]);
    gr[part]->SetTitle("#pi^{0} acceptance");
    aset(gr[part], "p_{T} [GeV]","acceptance", 0.,30., 0.,0.12);
    style(gr[part], part+20, part+1);
    if(part==0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");
    leg0->AddEntry(gr[part], pname[part], "P");
    TGraph *gr_sasha =  new TGraph( Form("data/sasha-acc-part%d.txt",part) );
    gr_sasha->Draw("C");
  }

  leg0->Draw();
  c3->Print("plots/Acceptance-pion.pdf");

  TFile *f_out = new TFile("data/Acceptance-pion.root", "RECREATE");
  for(int part=0; part<3; part++)
  {
    mcw( part, Form("part%d",part) );
    gr[part]->Write();
  }
  f_out->Close();
}
