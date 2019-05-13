#include "GlobalVars.h"

void draw_HadronResponse()
{
  const char *pname[2] = {"PbSc", "PbGl"};
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/HadronResponse-histo.root");
  THnSparse *hn_emcal = (THnSparse*)f->Get("hn_emcal");
  THnSparse *hn_cluster = (THnSparse*)f->Get("hn_cluster");
  THnSparse *hn_dc = (THnSparse*)f->Get("hn_dc");

  for(int em=0; em<2; em++)
    for(int part=0; part<2; part++)
    {
      int ic = part + 2*em;
      mc(ic, 6,5);
      for(int ipt=0; ipt<npT; ipt++)
      {
        mcd(ic, ipt+1);
        hn_emcal->GetAxis(5)->SetRange(em+1,em+1);  // is em shower?
        hn_emcal->GetAxis(4)->SetRange(secl[part],sech[part]);  // sector
        hn_emcal->GetAxis(0)->SetRange(ipt+1,ipt+1);  // pt

        TH1 *h_mc = hn_emcal->Projection(2);
        TH1 *h_reco = hn_emcal->Projection(3);
        const double max = h_reco->GetMaximum();
        h_mc->SetTitle( Form("p_{T}: %.1f-%.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
        aset(h_mc, "Energy [GeV]","", 0.,2.*pTbin[ipt], 0.,1.1*max);
        style(h_mc, 20, 1);
        style(h_reco, 21, 2);
        h_mc->DrawCopy();
        h_reco->DrawCopy("SAME");
        delete h_mc;
        delete h_reco;
      } // ipt
    } // part, em

  for(int part=0; part<2; part++)
  {
    mc(4+part);
    int ipt = 4;
    //for(int ipt=0; ipt<npT; ipt++)
    {
      mcd(4+part, ipt+1);
      hn_cluster->GetAxis(5)->SetRange(2,2);  // is photon shower?
      hn_cluster->GetAxis(4)->SetRange(secl[part],sech[part]);  // sector
      hn_cluster->GetAxis(0)->SetRange(ipt+1,30);  // pt

      TH1 *h_mc = hn_cluster->Projection(2);
      TH1 *h_reco = hn_cluster->Projection(3);
      const double max = h_mc->GetMaximum();
      h_mc->SetTitle( Form("p_{T}: %.1f-%.1f GeV",pTbin[ipt],30.) );
      aset(h_mc, "Energy [GeV]","", 0.,40., 0.,1.1*max);
      style(h_mc, 20, 1);
      style(h_reco, 21, 2);
      h_mc->DrawCopy();
      h_reco->DrawCopy("SAME");
      h_mc->GetXaxis()->SetRange(30,400);
      h_reco->GetXaxis()->SetRange(30,400);
      cout << "part " << part << ", E_reco/E_mc = " << h_reco->GetMean()/h_mc->GetMean() << endl;
      delete h_mc;
      delete h_reco;
    } // ipt
  } // part

  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      int ic = 6 + ns + 2*we;
      mc(ic, 6,5);
      for(int ipt=0; ipt<npT; ipt++)
      {
        mcd(ic, ipt+1);
        hn_dc->GetAxis(4)->SetRange(ns+1,ns+1);  // north or south
        hn_dc->GetAxis(5)->SetRange(we+1,we+1);  // west or east
        hn_dc->GetAxis(0)->SetRange(ipt+1,ipt+1);  // pt

        TH1 *h_mc = hn_dc->Projection(2);
        TH1 *h_reco = hn_dc->Projection(3);
        const double max = h_reco->GetMaximum();
        h_mc->SetTitle( Form("p_{T}: %.1f-%.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
        aset(h_mc, "Energy [GeV]","", 0.,2.*pTbin[ipt], 0.,1.1*max);
        style(h_mc, 20, 1);
        style(h_reco, 21, 2);
        h_mc->DrawCopy();
        h_reco->DrawCopy("SAME");
        delete h_mc;
        delete h_reco;
      } // ipt
    } // we, ns

  TFile *f_out = new TFile("data/HadronResponse.root", "RECREATE");
  for(int em=0; em<2; em++)
    for(int part=0; part<2; part++)
    {
      int ic = part + 2*em;
      mcw( ic, Form("em%d-%s",em,pname[part]) );
    }
  for(int part=0; part<2; part++)
    mcw( 4+part, Form("cluster-%s",pname[part]) );
  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      int ic = 6 + ns + 2*we;
      mcw( ic, Form("dc-ns%d-we%d",ns,we) );
    }
  f_out->Close();
}
