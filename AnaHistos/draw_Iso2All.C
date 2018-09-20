#include "GlobalVars.h"
#include "ReadGraph.h"
#include "FitMinv.h"
#include "Chi2Fit.h"

void draw_Iso2All()
{
  const double PI = TMath::Pi();

  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  TGraphErrors *gr[2];
  int igp[2] = {};
  for(int iph=0; iph<2; iph++)
  {
    gr[iph] = new TGraphErrors(npT);
    gr[iph]->SetName(Form("gr_%d",iph));
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[isolated][evtype][part]
  TH1 *h_1photon[2][3][3];
  TH2 *h2_2photon[3][3];
  TH2 *h2_isoboth[3][3];
  TH2 *h2_isopair[3][3];

  int bbc10cm = 1;
  int tof = 1;
  int prob = 1;
  int ival = 3;

  TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
  h_1photon_t->Reset();
  for(int isolated=0; isolated<2; isolated++)
    for(int evtype=1; evtype<3; evtype++)
      for(int part=0; part<3; part++)
      {
        h_1photon[isolated][evtype][part] = (TH1*)h_1photon_t->Clone(Form("h_1photon_iso%d_type%d_part%d",isolated,evtype,part));
        for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
          for(int pattern=0; pattern<3; pattern++)
          {
            int ih = sector + 8*pattern + 3*8*isolated + 2*3*8*tof + 2*2*3*8*prob + 2*2*2*3*8*evtype + 3*2*2*2*3*8*bbc10cm + 2*3*2*2*2*3*8*ival;
            TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
            h_1photon[isolated][evtype][part]->Add(h_tmp);
            delete h_tmp;
          }
        if(isolated==1)
          h_1photon[0][evtype][part]->Add(h_1photon[1][evtype][part]);
      }

  TH2 *h2_2photon_t = (TH2*)f->Get("h2_2photon_0");
  h2_2photon_t->Reset();

  for(int evtype=1; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h2_2photon[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_2photon_type%d_part%d",evtype,part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int pattern=0; pattern<3; pattern++)
          for(int isopair=0; isopair<2; isopair++)
            for(int isoboth=0; isoboth<2; isoboth++)
            {
              int ih = sector + 8*pattern + 3*8*isoboth + 2*3*8*isopair + 2*2*3*8*tof + 2*2*2*3*8*prob + 2*2*2*2*3*8*evtype + 3*2*2*2*2*3*8*bbc10cm + 2*3*2*2*2*2*3*8*ival;
              TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
              h2_2photon[evtype][part]->Add(h2_tmp);
              delete h2_tmp;
            }
    }

  for(int evtype=1; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h2_isoboth[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_isoboth_type%d_part%d",evtype,part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int pattern=0; pattern<3; pattern++)
          for(int isopair=0; isopair<2; isopair++)
          {
            int isoboth = 1;
            int ih = sector + 8*pattern + 3*8*isoboth + 2*3*8*isopair + 2*2*3*8*tof + 2*2*2*3*8*prob + 2*2*2*2*3*8*evtype + 3*2*2*2*2*3*8*bbc10cm + 2*3*2*2*2*2*3*8*ival;
            TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
            h2_isoboth[evtype][part]->Add(h2_tmp);
            delete h2_tmp;
          }
    }

  for(int evtype=1; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h2_isopair[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_isopair_type%d_part%d",evtype,part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int pattern=0; pattern<3; pattern++)
          for(int isoboth=0; isoboth<2; isoboth++)
          {
            int isopair = 1;
            int ih = sector + 8*pattern + 3*8*isoboth + 2*3*8*isopair + 2*2*3*8*tof + 2*2*2*3*8*prob + 2*2*2*2*3*8*evtype + 3*2*2*2*2*3*8*bbc10cm + 2*3*2*2*2*2*3*8*ival;
            TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
            h2_isopair[evtype][part]->Add(h2_tmp);
            delete h2_tmp;
          }
    }

  const double Pile[3] = {0.905, 0.905, 0.865};
  //const double IsoPile[3] = {1.27, 1.32, 1.23};
  //const double IsoPile[3] = {1.04, 1.07, 1.02};
  const double IsoPile[3] = {1.11, 1.11, 1.08};
  const double ePile = 0.02;
  const double A = 0.28;
  const double eA = 0.04;

  double xVeto[3][npT] = {}, Veto[3][npT] = {}, eVeto[3][npT] = {};
  double xMiss[3][npT] = {}, Miss[3][npT] = {}, eMiss[3][npT] = {};
  double xMerge[3][npT] = {}, Merge[3][npT] = {}, eMerge[3][npT] = {};
  double xBadPass[3][npT] = {}, BadPass[3][npT] = {}, eBadPass[3][npT] = {};

  double bck[2][npT] = {
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.10, 1.15, 1.20, 1.30,  1, 1, 1 },
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.08, 1.08, 1.11, 1.11, 1.11, 1.11, 1.11 }
  };

  double meff[2][npT] = {
    { 1, 1, 0.96, 0.97, 0.98, 0.985, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.985, 0.995, 0.995, 0.99, 0.98, 0.95, 1,1,1,1,1, },
    { 1, 1, 0.95, 0.97, 0.975, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.995, 0.995, 0.99, 0.99, 0.98, 0.98, 0.98, 0.98, 1, 1 }
  };

  for(int part=0; part<3; part++)
  {
    ReadGraph<TGraphAsymmErrors>("data/SelfVeto.root", part, xVeto[part], Veto[part], eVeto[part]);
    ReadGraph<TGraphErrors>("data/MissingRatio.root", part, xMiss[part], Miss[part], eMiss[part]);
    ReadGraph<TGraphAsymmErrors>("data/Merge-photon.root", part, xMerge[part], Merge[part], eMerge[part]);
    ReadGraph<TGraphErrors>("data/MergePassRate.root", part/2, xBadPass[part], BadPass[part], eBadPass[part]);
  }

  for(int ipt=2; ipt<npT; ipt++)
  {
    double xx, yy[2][3], eyy[2][3];

    for(int part=0; part<3; part++)
    {
      int evtype = 2;
      if(ipt < 22)  // <14GeV use ERT_4x4c
        evtype = 2;
      else  // >14GeV use ERT_4x4b
        evtype = 1;

      double nphoton = h_1photon[0][evtype][part]->GetBinContent(ipt+1);
      double nisophoton = h_1photon[1][evtype][part]->GetBinContent(ipt+1);

      TH1 *h_minv;

      double npion = 1., enpion = 1.;
      h_minv = (TH1*)h2_2photon[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, npion, enpion, true, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, npion, enpion, true, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, npion, enpion, false, 0.10,0.17);
      npion /= bck[part/2][ipt] * meff[part/2][ipt];
      delete h_minv;

      double nisoboth = 1., enisoboth = 1.;
      h_minv = (TH1*)h2_isoboth[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      // don't subtract background
      FitMinv(h_minv, nisoboth, enisoboth, false, 0.10,0.17);
      nisoboth /= 1.1;
      delete h_minv;

      double nisopair = 1., enisopair = 1.;
      h_minv = (TH1*)h2_isopair[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, nisopair, enisopair, true, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, nisopair, enisopair, true, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, nisopair, enisopair, false, 0.10,0.17);
      nisopair /= bck[part/2][ipt] * meff[part/2][ipt];
      delete h_minv;

      xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      int ipVeto = Get_ipt(xVeto[part], xx);
      int ipMiss = Get_ipt(xMiss[part], xx);
      double aMiss = Miss[part][ipMiss];
      double eaMiss = eMiss[part][ipMiss];
      int ipMerge = Get_ipt(xMerge[part], xx);
      double aMerge = Merge[part][ipMerge];
      double eaMerge = eMerge[part][ipMerge];
      int ipBadPass = Get_ipt(xBadPass[part], xx);
      double aBadPass = BadPass[part][ipBadPass];
      double eaBadPass = eBadPass[part][ipBadPass];
      if(ipt<23)
      {
        aBadPass = 0.;
        eaBadPass = 0.;
      }
      double aMissPass = aMiss + aMerge * aBadPass;
      double eaMissPass = sqrt( eaMiss*eaMiss + pow(eaMerge*aBadPass,2) + pow(aMerge*eaBadPass,2) );
      double aMissAll = aMiss + aMerge * 2.;
      double eaMissAll = sqrt( eaMiss*eaMiss + eaMerge*eaMerge*4. );

      double npi0 = npion + aMissPass * npion;
      double enpi0 = sqrt( pow(eaMissPass*npion,2) + pow((1.+aMissPass)*enpion,2) );
      double nisopi0 = nisoboth + aMissPass * nisopair;
      double enisopi0 = sqrt( enisoboth*enisoboth + pow(eaMissPass*nisopair,2) + pow(aMissPass*enisopair,2) );
      //yy[0][part] = ( nisopi0 * IsoPile[part] ) / ( npi0 * Pile[part] );
      yy[0][part] = nisopair / npion;
      eyy[0][part] = yy[0][part] * sqrt( pow(enisopi0/nisopi0,2) + pow(enpi0/npi0,2) + pow(ePile/Pile[part],2) + pow(ePile/IsoPile[part],2) );

      double ndir = nphoton - aMissPass * npion - A * ( 1. + aMissAll ) * npion;
      double endir = sqrt( nphoton + pow((1.+A)*(1.+aMissPass)*enpion,2) + pow((1.+A)*eaMissPass*npion,2) );
      double nisodir = nisophoton - ( nisoboth + aMissPass * nisopair ) - A * Veto[part][ipVeto] * ( 1. + aMissAll ) * nisopair;
      double enisodir = sqrt( nisophoton + enisoboth*enisoboth + pow(eaMissPass*nisopair,2) + pow(aMissPass*enisopair,2) );
      if( enisodir != enisodir ) enisodir = sqrt( nisophoton + enisoboth*enisoboth );
      yy[1][part] = ( nisodir * IsoPile[part] ) / ( ndir * Pile[part] );
      eyy[1][part] = yy[1][part] * sqrt( pow(endir/ndir,2) + pow(enisodir/nisodir,2) + pow(ePile/Pile[part],2) + pow(ePile/IsoPile[part],2) );
    } // part

    for(int iph=0; iph<2; iph++)
    {
      double ybar, eybar;
      Chi2Fit(3, yy[iph], eyy[iph], ybar, eybar);
      if( ybar > 0. && eybar > 0. && eybar < TMath::Infinity() )
      {
        gr[iph]->SetPoint(igp[iph], xx, ybar);
        gr[iph]->SetPointError(igp[iph], 0., eybar);
        igp[iph]++;
      }
    } // iph
  } // ipt

  mc();
  mcd();
  legi(0, 0.2,0.7,0.4,0.9);
  for(int iph=0; iph<2; iph++)
  {
    gr[iph]->Set(igp[iph]);
    aset(gr[iph], "p_{T} [GeV]", "Iso/All", 6.,30., 0.,1.4);
    style(gr[iph], iph+20, iph+1);
    if(iph==0)
      gr[iph]->Draw("AP");
    else
      gr[iph]->Draw("P");
  }
  leg0->AddEntry(gr[0], "#pi^{0}", "P");
  leg0->AddEntry(gr[1], "#gamma_{dir}", "P");
  leg0->Draw();
  c0->Print("plots/Iso2All-pair-notrack.pdf");
}
