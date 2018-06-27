#include "GlobalVars.h"

void draw_YieldCmpSample()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  TGraph *gr[3];
  int igp[3] = {};
  for(int part=0; part<3; part++)
    gr[part] =  new TGraph(25);

  const int runnumber = 392294;
  TFile *f_histo = new TFile( Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/13568/data/PhotonHistos-%d.root",runnumber) );
  TFile *f_node = new TFile( Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber) );

  THnSparse *hn_histo = (THnSparse*)f_histo->Get("hn_bbc_pion");
  TAxis *axis_histo_sec = hn_histo->GetAxis(0);
  TAxis *axis_histo_pt = hn_histo->GetAxis(1);
  TAxis *axis_histo_bbc = hn_histo->GetAxis(3);
  //hn_histo->Reset();
  //for(int isoboth=0; isoboth<2; isoboth++)
  //  for(int isopair=0; isopair<2; isopair++)
  //  {
  //    int cut = 3;
  //    int evtype = 2;
  //    int bbc10cm = 1;
  //    int ih = 2*4*3*2*isoboth + 4*3*2*isopair + 3*2*cut + 2*evtype + bbc10cm;
  //    THnSparse *hn_tmp = (THnSparse*)f_histo->Get(Form("hn_2photon_%d",ih));
  //    hn_histo->Add(hn_tmp);
  //    delete hn_tmp;
  //  }

  THnSparse *hn_node = (THnSparse*)f_node->Get("hn_bbc_pion");
  TAxis *axis_node_sec = hn_node->GetAxis(0);
  TAxis *axis_node_pt = hn_node->GetAxis(1);
  TAxis *axis_node_bbc = hn_node->GetAxis(3);

  TH3 *h3_histo = (TH3*)f_histo->Get("h3_minv");
  TH3 *h3_node = (TH3*)f_node->Get("h3_minv");

  for(int part=0; part<3; part++)
    for(int ipt=0; ipt<30; ipt++)
    {
      TH1 *h_minv;

      axis_histo_bbc->SetRange(2,2);
      axis_histo_sec->SetRange(secl[part],sech[part]);
      axis_histo_pt->SetRange(ipt+1,ipt+1);
      h_minv = hn_histo->Projection(2);
      //h_minv = h3_histo->ProjectionZ("h_minv", secl[part],sech[part], ipt+1,ipt+1);
      double npion_histo = h_minv->Integral(120,160);
      delete h_minv;

      axis_node_bbc->SetRange(2,2);
      axis_node_sec->SetRange(secl[part],sech[part]);
      axis_node_pt->SetRange(ipt+1,ipt+1);
      h_minv = hn_node->Projection(2);
      //h_minv = h3_node->ProjectionZ("h_minv", secl[part],sech[part], ipt+1,ipt+1);
      double npion_node = h_minv->Integral(120,160);
      delete h_minv;

      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double ratio = npion_histo / npion_node;
      gr[part]->SetPoint(igp[part], xx, ratio);
      igp[part]++;
    }

  mc();
  mcd();
  legi();

  for(int part=0; part<3; part++)
  {
    gr[part]->Set(igp[part]);
    aset(gr[part], "p_{T} [GeV]","#frac{histo}{node}", 0.,30., 0.,2.);
    style(gr[part], part+24, part+1);
    if(part == 0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");
    leg0->AddEntry(gr[part], pname[part], "P");
  }
  leg0->Draw();

  c0->Print("plots/YieldCmpSample.pdf");
}
