#include "GlobalVars.h"
#include "QueryTree.h"

void draw_HadronResponse()
{
  const char *pname[2] = {"PbSc", "PbGl"};
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  QueryTree *qt_res = new QueryTree("data/HadronResponse.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/HadronResponse-histo-photon.root");
  TFile *f_hadron = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/HadronResponse-hadron-histo.root");
  TFile *f_em = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/HadronResponse-em-histo.root");

  THnSparse *hn_cluster = (THnSparse*)f->Get("hn_cluster");
  THnSparse *hn_dc = (THnSparse*)f->Get("hn_dc");
  THnSparse *hn_emcal[2];  // hadron:0 ; em:1
  hn_emcal[0] = (THnSparse*)f_hadron->Get("hn_emcal");
  hn_emcal[1] = (THnSparse*)f_em->Get("hn_emcal");

  for(int em=0; em<2; em++)
    for(int part=0; part<2; part++)
    {
      int ic = part + 2*em;
      mc(ic, 6,5);
      for(int ipt=0; ipt<npT; ipt++)
      {
        mcd(ic, ipt+1);
        hn_emcal[em]->GetAxis(2)->SetRange(secl[part],sech[part]);  // sector
        hn_emcal[em]->GetAxis(0)->SetRange(ipt+1,ipt+1);  // pt

        hn_emcal[em]->GetAxis(3)->SetRange(1,1);  // truth
        TH1 *h_mc = hn_emcal[em]->Projection(1);
        h_mc->SetName("h_mc");
        hn_emcal[em]->GetAxis(3)->SetRange(2,2);  // reco
        TH1 *h_reco = hn_emcal[em]->Projection(1);
        h_reco->SetName("h_reco");

        const double max = h_mc->GetMaximum();
        h_mc->SetTitle( Form("p_{T}: %.1f-%.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
        aset(h_mc, "Energy (GeV)","", 0.,2.*pTbin[ipt], 0.,1.1*max);
        style(h_mc, 20, 1);
        style(h_reco, 21, 2);
        h_mc->DrawCopy("HISTO");
        h_reco->DrawCopy("HISTO SAME");
        delete h_mc;
        delete h_reco;
      } // ipt
    } // part, em

  mc(6, 2,1);
  for(int part=0; part<2; part++)
  {
    mc(4+part, 3,3);
    const int step = 4;
    int ipad = 1;
    for(int ipt=0; ipt<npT; ipt+=step)
    {
      int ilow = ipt;
      int ihigh = ipt+step < npT ? ipt+step : npT;
      hn_cluster->GetAxis(4)->SetRange(secl[part],sech[part]);  // sector
      hn_cluster->GetAxis(0)->SetRange(ilow+1,ihigh);  // pt

      TH1 *h_mc = hn_cluster->Projection(2);
      TH1 *h_reco = hn_cluster->Projection(3);
      const double max = h_mc->GetMaximum();

      double st_mc[4], st_reco[4];
      h_mc->GetXaxis()->SetRange(3,400);
      h_reco->GetXaxis()->SetRange(3,400);
      h_mc->GetStats(st_mc);
      h_reco->GetStats(st_reco);
      double xpt = ( pTbin[ilow] + pTbin[ihigh] ) / 2.;
      double sum_mc = st_mc[2];
      double erel_sum_mc = sqrt(st_mc[1]) / st_mc[0];
      double sum_reco = st_reco[2];
      double erel_sum_reco = sqrt(st_reco[1]) / st_reco[0];
      double rsum = sum_reco / sum_mc;
      double ersum = 1. * sqrt( erel_sum_mc*erel_sum_mc + erel_sum_reco*erel_sum_reco );
      if( rsum > 0.1 && rsum < 10. )
        qt_res->Fill(ipt/4, part, xpt, rsum, ersum);

      mcd(4+part, ipad);
      gPad->SetLogy();
      h_mc->SetTitle( Form("p_{T}: %.1f-%.1f GeV, reco/mc=%.2f",pTbin[ilow],pTbin[ihigh],rsum) );
      aset(h_mc, "Energy (GeV)","", 0.,40., 1e-8,2.*max);
      style(h_mc, 20, 1);
      style(h_reco, 21, 2);
      h_mc->DrawCopy();
      h_reco->DrawCopy("SAME");
      delete h_mc;
      delete h_reco;
      ipad++;
    } // ipt

    mcd(6, 1+part);
    TGraphErrors *gr = qt_res->Graph(part);
    gr->SetTitle(pname[part]);
    aset(gr, "p_{T} (GeV/c)","sum_reco/sum_mc", 0.,30., 0.,2.);
    style(gr, 20, 1);
    gr->Draw();
    gr->Fit("pol0", "","", 3.,30.);
  } // part

  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      int ic = 7 + ns + 2*we;
      mc(ic, 6,5);
      for(int ipt=0; ipt<npT; ipt++)
      {
        mcd(ic, ipt+1);
        hn_dc->GetAxis(2)->SetRange(ns+1,ns+1);  // north or south
        hn_dc->GetAxis(3)->SetRange(we+1,we+1);  // west or east
        hn_dc->GetAxis(0)->SetRange(ipt+1,ipt+1);  // pt

        hn_dc->GetAxis(4)->SetRange(1,1);  // mc
        TH1 *h_mc = hn_dc->Projection(1);
        h_mc->SetName("h_mc");
        hn_dc->GetAxis(4)->SetRange(2,2);  // reco
        TH1 *h_reco = hn_dc->Projection(1);
        h_reco->SetName("h_reco");

        const double max = h_reco->GetMaximum();
        h_mc->SetTitle( Form("p_{T}: %.1f-%.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
        aset(h_mc, "Energy (GeV)","", 0.,2.*pTbin[ipt], 0.,1.1*max);
        style(h_mc, 20, 1);
        style(h_reco, 21, 2);
        h_mc->DrawCopy("HISTO");
        h_reco->DrawCopy("HISTO SAME");
        delete h_mc;
        delete h_reco;
      } // ipt
    } // we, ns

  qt_res->Write();
  for(int em=0; em<2; em++)
    for(int part=0; part<2; part++)
    {
      int ic = part + 2*em;
      mcw( ic, Form("em%d-%s",em,pname[part]) );
    }
  for(int part=0; part<2; part++)
  {
    int ic = 4 + part;
    mcw( ic, Form("cluster-%s",pname[part]) );
  }
  mcw( 6, "rsum" );
  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      int ic = 7 + ns + 2*we;
      mcw( ic, Form("dc-ns%d-we%d",ns,we) );
    }
  qt_res->Close();
}
