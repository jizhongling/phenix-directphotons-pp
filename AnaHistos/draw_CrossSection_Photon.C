#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "Chi2Fit.h"

void draw_CrossSection_Photon()
{
  const double PI = TMath::Pi();

  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

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
  TH2 *h2_2photon[3][3];
  TH2 *h2_2photon2pt[3][3];

  int bbc10cm = 1;
  int ival = 1;

  TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
  h_1photon_t = (TH1*)h_1photon_t->Clone();
  h_1photon_t->Reset();
  for(int evtype=1; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h_1photon[evtype][part] = (TH1*)h_1photon_t->Clone(Form("h_1photon_type%d_part%d",evtype,part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int evenodd=0; evenodd<2; evenodd++)
          for(int pattern=0; pattern<3; pattern++)
            for(int isolated=0; isolated<2; isolated++)
            {
              int ih = sector + 8*evenodd + 2*8*pattern + 3*2*8*isolated + 2*3*2*8*evtype + 3*2*3*2*8*bbc10cm + 2*3*2*3*2*8*ival;
              TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
              h_1photon[evtype][part]->Add(h_tmp);
              delete h_tmp;
            }
    }

  TH2 *h2_2photon_t = (TH2*)f->Get("h2_2photon_0");
  h2_2photon_t = (TH2*)h2_2photon_t->Clone();
  h2_2photon_t->Reset();
  for(int evtype=1; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h2_2photon[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_2photon_type%d_part%d",evtype,part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int pattern=0; pattern<3; pattern++)
          for(int evenodd=0; evenodd<2; evenodd++)
            for(int isoboth=0; isoboth<2; isoboth++)
              for(int isopair=0; isopair<2; isopair++)
              {
                int ih = sector + 8*evenodd + 2*8*pattern + 3*2*8*isoboth + 2*3*2*8*isopair + 2*2*3*2*8*evtype + 3*2*2*3*2*8*bbc10cm + 2*3*2*2*3*2*8*ival;
                TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
                h2_2photon[evtype][part]->Add(h2_tmp);
                delete h2_tmp;
              }
    }

  for(int evtype=1; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h2_2photon2pt[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_2photon2pt_type%d_part%d",evtype,part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int pattern=0; pattern<3; pattern++)
          for(int evenodd=0; evenodd<2; evenodd++)
            for(int isoboth=0; isoboth<2; isoboth++)
              for(int isopair=0; isopair<2; isopair++)
              {
                int ih = sector + 8*evenodd + 2*8*pattern + 3*2*8*isoboth + 2*3*2*8*isopair + 2*2*3*2*8*evtype + 3*2*2*3*2*8*bbc10cm + 2*3*2*2*3*2*8*ival;
                TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon2pt_%d",ih));
                h2_2photon2pt[evtype][part]->Add(h2_tmp);
                delete h2_tmp;
              }
    }

  const double DeltaEta = 1.0;
  //const double NBBC =  3.59e11;  // from DAQ
  const double NBBC =  3.54e11;  // from rejection power
  const double XBBC = 32.51e9;
  const double eXBBC = 3.24e9;
  const double Pile[3] = {0.894, 0.886, 0.920};
  const double ePile = 0.01;
  const double TrigBBC = 0.91;
  const double eTrigBBC = 0.01;
  const double ToF[3] = {0.992, 0.992, 0.997};
  const double eToF[3] = {0.002, 0.002, 0.002};
  const double Conv[3] = {0.849, 0.959, 0.959};
  const double eConv[3] = {0.027, 0.023, 0.023};
  const double Norm[3] = {0.321, 0.314, 0.243};
  const double eNorm[3] = {0.005, 0.006, 0.005};
  const double A = 0.24;
  const double eA = 0.04;

  QueryTree *qt_acc = new QueryTree("data/Acceptance-photon.root");
  QueryTree *qt_ert = new QueryTree("data/ERTEff-photon.root");
  QueryTree *qt_miss = new QueryTree("data/MissingRatio.root");
  QueryTree *qt_mergecorr1 = new QueryTree("data/MergeCorr-1photon.root");
  QueryTree *qt_mergecorr2 = new QueryTree("data/MergeCorr-2photon.root");

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
    mc(part, 6,5);
    mc(part+3, 6,5);
  }

  for(int ipt=0; ipt<npT; ipt++)
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

      mcd(part, ipt+1);
      double n2photon = 1., en2photon = 1.;
      TH1 *h_minv = (TH1*)h2_2photon[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, n2photon, en2photon, true, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, n2photon, en2photon, true, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, n2photon, en2photon, false, 0.10,0.17);
      n2photon /= bck[part/2][ipt] * meff[part/2][ipt];
      delete h_minv;

      mcd(part+3, ipt+1);
      double n2photon2pt = 1., en2photon2pt = 1.;
      h_minv = (TH1*)h2_2photon2pt[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, n2photon2pt, en2photon2pt, true, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, n2photon2pt, en2photon2pt, true, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, n2photon2pt, en2photon2pt, false, 0.10,0.17);
      n2photon2pt /= bck[part/2][ipt] * meff[part/2][ipt];
      delete h_minv;

      double xx, Acc, eAcc, TrigERT, eTrigERT, Miss, eMiss, MergeCorr2, eMergeCorr2, MergeCorr1, eMergeCorr1;
      qt_acc->Query(ipt, part, xx, Acc, eAcc);
      qt_ert->Query(ipt, part/2, xx, TrigERT, eTrigERT);
      qt_miss->Query(ipt, part, xx, Miss, eMiss);
      qt_mergecorr1->Query(ipt, part, xx, MergeCorr1, eMergeCorr1);
      qt_mergecorr2->Query(ipt, part, xx, MergeCorr2, eMergeCorr2);

      if(ipt >= 20)
      {
        if(part<2)
        {
          TrigERT = 0.948;
          eTrigERT = 0.004;
        }
        else
        {
          TrigERT = 0.651;
          eTrigERT = 0.008;
        }
      }

      double ndir = nphoton - ( Miss + MergeCorr1 ) * n2photon - MergeCorr2 * n2photon2pt;
      double endir = sqrt( nphoton + pow((eMiss+eMergeCorr1)*n2photon,2) + pow((Miss+MergeCorr2)*en2photon,2) + pow(eMergeCorr2*n2photon2pt,2) + pow(MergeCorr2*en2photon2pt,2) );
      if(ipt >= 22)  // >14GeV use ERT_4x4b
      {
        ndir *= Norm[part];
        endir *= Norm[part] * (1+eNorm[part]);
      }

      yy[part] = (XBBC/NBBC) / (2*PI*xx) / (pTbin[ipt+1]-pTbin[ipt]) / DeltaEta
        * ndir / Acc / TrigERT / Prob[part/2][ipt]
        / ToF[part] / Conv[part] / TrigBBC * Pile[part];
      eyy[part] = yy[part] * sqrt( pow(endir/ndir,2)
          + pow(eAcc/Acc,2)
          + pow(eTrigERT/TrigERT,2)
          + pow(eProb/Prob[part/2][ipt],2)
          + pow(eToF[part]/ToF[part],2) + pow(eConv[part]/Conv[part],2)
          //+ pow(eTrigBBC/TrigBBC,2) + pow(ePile/Pile[part],2) + pow(eXBBC/XBBC,2)
          );
      if( TMath::Finite(yy[part] + eyy[part]) )
      {
        gr[part]->SetPoint(igp[part], xx, yy[part]);
        gr[part]->SetPointError(igp[part], 0., eyy[part]);
        igp[part]++;
      }
    } // part

    double ybar, eybar;
    Chi2Fit(3, yy, eyy, ybar, eybar);
    if( TMath::Finite(ybar + eybar) )
    {
      gr[3]->SetPoint(igp[3], xx, ybar);
      gr[3]->SetPointError(igp[3], 0., eybar);
      igp[3]++;
    }
  } // ipt

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
      gr[part]->Draw("AP");
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
  c6->Print("plots/CrossSection-photon.pdf");

  TFile *f_out = new TFile("data/CrossSection-photon.root", "RECREATE");
  for(int part=0; part<4; part++)
  {
    if(part<3)
    {
      mcw( part, Form("Minv-2photon-part%d",part) );
      mcw( part+3, Form("Minv-2photon2pt-part%d",part) );
    }
    gr[part]->Write();
  }
  f_out->Close();
}
