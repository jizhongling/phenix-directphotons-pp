#include "GlobalVars.h"
#include "ReadGraph.h"
#include "FitMinv.h"
#include "Chi2Fit.h"

void draw_CrossSection_IsoPhoton()
{
  const double PI = TMath::Pi();

  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  // function for pT weight for direct photon
  cross_ph = new TF1("cross_ph", "x**(-[1]-[2]*log(x/[0]))*(1-(x/[0])**2)**[3]*[4]", 0, 30);
  cross_ph->SetParameters(255., 5.98, 0.273, 14.43, 1.);
  double ndata = 0.;
  double nfit = 0.;

  TGraphErrors *gr[4];  // PbScW, PbScE, PbGl, Combined
  int igp[4] = {};
  for(int part=0; part<4; part++)
  {
    gr[part] = new TGraphErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[evtype][part]
  TH1 *h_1photon[3][3];
  TH2 *h2_isoboth[3][3];
  TH2 *h2_isopair[3][3];

  int bbc10cm = 1;
  int tof = 1;
  int prob = 1;

  TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
  h_1photon_t->Reset();
  for(int evtype=1; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h_1photon[evtype][part] = (TH1*)h_1photon_t->Clone(Form("h_1photon_type%d_part%d",evtype,part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int pattern=0; pattern<3; pattern++)
        {
          int isolated = 1;
          int ih = sector + 8*pattern + 3*8*isolated + 2*3*8*tof + 2*2*3*8*prob + 2*2*2*3*8*evtype + 3*2*2*2*3*8*bbc10cm;
          TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
          h_1photon[evtype][part]->Add(h_tmp);
          delete h_tmp;
        }
    }

  TH2 *h2_2photon_t = (TH2*)f->Get("h2_2photon_0");
  h2_2photon_t->Reset();
  for(int evtype=1; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h2_isoboth[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_isoboth_type%d_part%d",evtype,part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int pattern=0; pattern<3; pattern++)
          for(int isopair=0; isopair<2; isopair++)
          {
            int isoboth = 1;
            int ih = sector + 8*pattern + 3*8*isoboth + 2*3*8*isopair + 2*2*3*8*tof + 2*2*2*3*8*prob + 2*2*2*2*3*8*evtype + 3*2*2*2*2*3*8*bbc10cm;
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
            int ih = sector + 8*pattern + 3*8*isoboth + 2*3*8*isopair + 2*2*3*8*tof + 2*2*2*3*8*prob + 2*2*2*2*3*8*evtype + 3*2*2*2*2*3*8*bbc10cm;
            TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
            h2_isopair[evtype][part]->Add(h2_tmp);
            delete h2_tmp;
          }
    }

  const double DeltaEta = 0.5;
  //const double NBBC =  3.59e11;  // from DAQ
  const double NBBC =  3.54e11;  // from rejection power
  const double XBBC = 32.51e9;
  const double eXBBC = 3.24e9;
  const double Pile[3] = {1.27, 1.32, 1.23};
  const double ePile = 0.02;
  const double TrigBBC = 0.91;
  const double eTrigBBC = 0.01;
  const double ToF[3] = {0.992, 0.992, 0.997};
  const double eToF[3] = {0.002, 0.002, 0.002};
  const double Conv[3] = {0.849, 0.959, 0.959};
  const double eConv[3] = {0.027, 0.023, 0.023};
  const double Norm[3] = {0.321, 0.314, 0.243};
  const double eNorm[3] = {0.005, 0.006, 0.005};
  const double A = 0.28;
  const double eA = 0.04;

  double xAcc[3][npT] = {}, Acc[3][npT] = {}, eAcc[3][npT] = {};
  double xTrigERT[3][npT] = {}, TrigERT[3][npT] = {}, eTrigERT[3][npT] = {};
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

  double Prob[2][npT] = {
    { 1, 1, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 },
    { 1, 1, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97 }
  };
  double eProb = 0.02;

  for(int part=0; part<3; part++)
  {
    ReadGraph<TGraphAsymmErrors>("data/Acceptance-photon.root", part, xAcc[part], Acc[part], eAcc[part]);
    ReadGraph<TGraphAsymmErrors>("data/ERTEff-photon.root", part/2, xTrigERT[part], TrigERT[part], eTrigERT[part]);
    ReadGraph<TGraphAsymmErrors>("data/SelfVeto.root", part, xVeto[part], Veto[part], eVeto[part]);
    ReadGraph<TGraphErrors>("data/MissingRatio.root", part, xMiss[part], Miss[part], eMiss[part]);
    ReadGraph<TGraphAsymmErrors>("data/Merge-photon.root", part, xMerge[part], Merge[part], eMerge[part]);
    ReadGraph<TGraphErrors>("data/MergePassRate.root", part/2, xBadPass[part], BadPass[part], eBadPass[part]);
    mc(part, 6,5);
    mc(part+3, 6,5);
  }
  TrigERT[0][0] = TrigERT[1][0] = 0.948;
  eTrigERT[0][0] = eTrigERT[1][0] = 0.004;
  TrigERT[2][0] = 0.651;
  eTrigERT[2][0] = 0.008;

  for(int ipt=2; ipt<npT; ipt++)
  {
    double xx, yy[3], eyy[3];

    for(int part=0; part<3; part++)
    {
      int evtype = 2;
      if(ipt < 22)  // <14GeV use ERT_4x4c
        evtype = 2;
      else  // >14GeV use ERT_4x4b
        evtype = 1;

      double nphoton = h_1photon[evtype][part]->GetBinContent(ipt+1);

      TH1 *h_minv;

      mcd(part, ipt+1);
      double nisoboth = 1., enisoboth = 1.;
      h_minv = (TH1*)h2_isoboth[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      // don't subtract background
      FitMinv(h_minv, nisoboth, enisoboth, false, 0.10,0.17);
      nisoboth /= 1.1;
      delete h_minv;

      mcd(part+3, ipt+1);
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
      int ipAcc = Get_ipt(xAcc[part], xx);
      int ipTrigERT = Get_ipt(xTrigERT[part], xx);
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
      if(ipt >= 20)
        ipTrigERT = 0;
      double aMissPass = aMiss + aMerge * aBadPass;
      double eaMissPass = sqrt( eaMiss*eaMiss + pow(eaMerge*aBadPass,2.) + pow(aMerge*eaBadPass,2.) );
      double aMissAll = aMiss + aMerge * 2.;
      double eaMissAll = sqrt( eaMiss*eaMiss + eaMerge*eaMerge*4. );
      double ndir = nphoton - ( nisoboth + aMissPass * nisopair ) - A * Veto[part][ipVeto] * ( 1. + aMissAll ) * nisopair;
      double endir = sqrt( nphoton + enisoboth*enisoboth + pow(eaMissPass*nisopair,2.) + pow(aMissPass*enisopair,2.) );
      if( endir != endir ) endir = sqrt( nphoton + enisoboth*enisoboth );
      if(ipt >= 22)  // >14GeV use ERT_4x4b
      {
        ndir *= Norm[part];
        endir *= Norm[part] * (1+eNorm[part]);
      }
      yy[part] = (XBBC/NBBC) / (2*PI*xx) / (pTbin[ipt+1]-pTbin[ipt]) / DeltaEta
        * ndir / Acc[part][ipAcc] / TrigERT[part][ipTrigERT] / Prob[part/2][ipt]
        / ToF[part] / Conv[part] / TrigBBC * Pile[part];
      eyy[part] = yy[part] * sqrt( pow(endir/ndir,2.)
          + pow(eAcc[part][ipAcc]/Acc[part][ipAcc],2.)
          + pow(eTrigERT[part][ipTrigERT]/TrigERT[part][ipTrigERT],2.)
          + pow(eProb/Prob[part/2][ipt],2.)
          + pow(eToF[part]/ToF[part],2.) + pow(eConv[part]/Conv[part],2.)
          //+ pow(eTrigBBC/TrigBBC,2.) + pow(ePile/Pile[part],2.) + pow(eXBBC/XBBC,2.)
          );
      if( yy[part] > 0. && eyy[part] > 0. && eyy[part] < TMath::Infinity() )
      {
        gr[part]->SetPoint(igp[part], xx, yy[part]);
        gr[part]->SetPointError(igp[part], 0., eyy[part]);
        igp[part]++;
      }
    } // part

    double ybar, eybar;
    Chi2Fit(3, yy, eyy, ybar, eybar);
    if( ybar > 0. && eybar > 0. && eybar < TMath::Infinity() )
    {
      gr[3]->SetPoint(igp[3], xx, ybar);
      gr[3]->SetPointError(igp[3], 0., eybar);
      igp[3]++;
    }
    if( ipt >= 12 )
    {
      ndata += ybar;
      nfit += cross_ph->Eval(xx);
    }
  } // ipt
  cross_ph->SetParameter(4, ndata/nfit);

  mc(6, 2,1);
  legi(0, 0.4,0.7,0.7,0.9);

  for(int part=0; part<4; part++)
  {
    gr[part]->Set(igp[part]);
    mcd(6, part/3+1);
    gPad->SetLogy();
    aset(gr[part], "p_{T} [GeV]", "Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{-3}]", 6.,30., 1e-1, 1e4);
    style(gr[part], part+20, part+1);
    if(part%3==0)
    {
      gr[part]->Draw("AP");
      cross_ph->DrawCopy("LSAME");
    }
    else
      gr[part]->Draw("P");
    if(part<3)
      leg0->AddEntry(gr[part], pname[part], "P");
  }

  gr[0]->SetTitle("Separated");
  gr[3]->SetTitle("Combined");
  mcd(6, 1);
  leg0->Draw();
  mcd(6, 2);
  c6->Print("plots/CrossSection-isophoton.pdf");

  TFile *f_out = new TFile("data/CrossSection-isophoton.root", "RECREATE");
  for(int part=0; part<4; part++)
  {
    if(part<3)
    {
      mcw( part, Form("Minv-isoboth-part%d",part) );
      mcw( part+3, Form("Minv-isopair-part%d",part) );
    }
    gr[part]->Write();
  }
  f_out->Close();
}
