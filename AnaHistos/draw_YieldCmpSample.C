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

  const int runnumber = 393180;
  //TFile *f_histo = new TFile( Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/13597/data/PhotonHistos-%d.root",runnumber) );
  TFile *f_histo = new TFile( Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root",runnumber) );
  TFile *f_node = new TFile( Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-total.root",runnumber) );

  // h_histo[part]
  TH2 *h2_histo[3];
  TH2 *h2_histo_t = (TH2*)f_histo->Get("h2_2photon_0");
  h2_histo_t->Reset();
  int bbc10cm = 1;
  int evtype = 2;
  int cut = 3;
  for(int part=0; part<3; part++)
  {
    h2_histo[part] = (TH2*)h2_histo_t->Clone(Form("h2_histo_%d",part));
    for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
      for(int pattern=0; pattern<3; pattern++)
        for(int isoboth=0; isoboth<2; isoboth++)
          for(int isopair=0; isopair<2; isopair++)
          {
            int ih = sector + 8*pattern + 3*8*isoboth + 2*3*8*isopair + 2*2*3*8*cut + 4*2*2*3*8*evtype + 3*4*2*2*3*8*bbc10cm;
            TH2 *h2_tmp = (TH2*)f_histo->Get(Form("h2_2photon_%d",ih));
            h2_histo[part]->Add(h2_tmp);
            delete h2_tmp;
          }
  }

  THnSparse *hn_node = (THnSparse*)f_node->Get("hn_2photon");
  TAxis *axis_node_sec = hn_node->GetAxis(0);
  TAxis *axis_node_pt = hn_node->GetAxis(1);
  TAxis *axis_node_cut = hn_node->GetAxis(4);
  TAxis *axis_node_type = hn_node->GetAxis(5);
  axis_node_type->SetRange(3,3);
  axis_node_cut->SetRange(4,4);

  for(int part=0; part<3; part++)
    for(int ipt=0; ipt<30; ipt++)
    {
      TH1 *h_minv;

      h_minv = h2_histo[part]->ProjectionY("h_minv", ipt+1,ipt+1);
      double npion_histo = h_minv->Integral(120,160);
      delete h_minv;

      axis_node_sec->SetRange(secl[part],sech[part]);
      axis_node_pt->SetRange(ipt+1,ipt+1);
      h_minv = hn_node->Projection(2);
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
    aset(gr[part], "p_{T} [GeV]","#frac{histo}{node}", 0.,30., 0.9999,1.0001);
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
