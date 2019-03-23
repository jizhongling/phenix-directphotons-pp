#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "Chi2Fit.h"

void draw_CrossSection_Photon()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  const double PI = TMath::Pi();
  const double DeltaEta = 1.0;
  //const double NBBC =  3.59e11;  // from DAQ
  const double NBBC =  3.54e11;  // from rejection power
  const double XBBC = 32.51e9;
  const double eXBBC = 3.24e9;
  const double Pile[3] = {0.890, 0.890, 0.916};
  const double ePile = 0.01;
  const double TrigBBC = 0.91;
  const double eTrigBBC = 0.01;
  const double ToF[3] = {0.992, 0.992, 0.997};
  const double eToF[3] = {0.002, 0.002, 0.002};
  const double Conv[3] = {0.849, 0.959, 0.959};
  const double eConv[3] = {0.027, 0.023, 0.023};
  const double Norm[3] = {0.320, 0.320, 0.243};
  const double eNorm[3] = {0.005, 0.006, 0.005};
  const double A = 0.24;
  const double eA = 0.04;

  const double Prob[2][npT] = {
    { 1, 1, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 },
    { 1, 1, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97 }
  };
  const double eProb = 0.02;

  QueryTree *qt_cross = new QueryTree("data/CrossSection-photon.root", "RECREATE");
  QueryTree *qt_asee = new QueryTree("data/ASee.root", "RECREATE");

  QueryTree *qt_acc = new QueryTree("data/Acceptance-photon.root");
  QueryTree *qt_ert = new QueryTree("data/ERTEff-photon.root");
  QueryTree *qt_miss = new QueryTree("data/MissingRatio.root");
  QueryTree *qt_miss_eta = new QueryTree("data/MissingRatio-eta.root");
  QueryTree *qt_merge1 = new QueryTree("data/Merge-1photon.root");
  QueryTree *qt_merge2 = new QueryTree("data/Merge-2photon.root");
  QueryTree *qt_badpass = new QueryTree("data/MergePassRate.root");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[evtype][part]
  TH1 *h_1photon[3][3];
  TH2 *h2_2photon[3][3];
  TH2 *h2_2photon2pt[3][3];

  int bbc10cm = 1;
  int tof = 1;
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
              int ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isolated + 8*2*3*2*evtype + 8*2*3*2*3*bbc10cm + 8*2*3*2*3*2*ival + 8*2*3*2*3*2*4*(tof-1);
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
                int ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isoboth + 8*2*3*2*isopair + 8*2*3*2*2*evtype + 8*2*3*2*2*3*bbc10cm + 8*2*3*2*2*3*2*ival + 8*2*3*2*2*3*2*4*(tof-1);
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
                int ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isoboth + 8*2*3*2*isopair + 8*2*3*2*2*evtype + 8*2*3*2*2*3*bbc10cm + 8*2*3*2*2*3*2*ival + 8*2*3*2*2*3*2*4*(tof-1);
                TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon2pt_%d",ih));
                h2_2photon2pt[evtype][part]->Add(h2_tmp);
                delete h2_tmp;
              }
    }

  for(int part=0; part<3; part++)
  {
    mc(part, 6,5);
    mc(part+3, 6,5);
  }

  for(int ipt=0; ipt<npT; ipt++)
  {
    double xpt, yy[3], eyy[3];

    for(int part=0; part<3; part++)
    {
      int evtype = 2;
      if(ipt < 22)  // <14GeV use ERT_4x4c
        evtype = 2;
      else  // >14GeV use ERT_4x4b
        evtype = 1;

      double nphoton = h_1photon[evtype][part]->GetBinContent(ipt+1);
      double enphoton = sqrt(nphoton);

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

      double xpt, Acc, eAcc, TrigERT, eTrigERT, Miss, eMiss, MissEta, eMissEta, Merge1, eMerge1, Merge2, eMerge2, BadPass, eBadPass;
      qt_acc->Query(ipt, part, xpt, Acc, eAcc);
      qt_ert->Query(ipt, part/2, xpt, TrigERT, eTrigERT);
      qt_miss->Query(ipt, part, xpt, Miss, eMiss);
      qt_miss_eta->Query(ipt, part, xpt, MissEta, eMissEta);
      qt_merge1->Query(ipt, part, xpt, Merge1, eMerge1);
      qt_merge2->Query(ipt, part, xpt, Merge2, eMerge2);
      qt_badpass->Query(ipt, part/2, xpt, BadPass, eBadPass);

      if(ipt >= 20)
      {
        if(part < 2)
        {
          TrigERT = 0.949;
          eTrigERT = 0.004;
        }
        else
        {
          TrigERT = 0.651;
          eTrigERT = 0.008;
        }
      }

      double ASee = A * (1.+MissEta)/(1.+2.*MissEta) * (1+2.*Miss+Merge1);
      double eASee = sqrt((pow(A,2)*pow(eMerge1,2)* pow(1. + Miss,2))/pow(1. + 2.*Miss,2) + (pow(eA,2)*pow(1. + Miss,2)* pow(1. + Merge1 + 2.*Miss,2))/ pow(1. + 2.*Miss,2) + pow(eMiss,2)*pow((2.*A*(1. + Miss))/ (1. + 2.*Miss) - (2.*A*(1. + Miss)* (1. + Merge1 + 2.*Miss))/ pow(1. + 2.*Miss,2) + (A*(1. + Merge1 + 2.*Miss))/(1. + 2.*Miss) ,2));
      qt_asee->Fill(ipt, part, xpt, ASee, eASee);

      double ndir = nphoton/Conv[part] - (1. + Miss + Merge1*Conv[part]*(1.-Conv[part]) + ASee) * n2photon/pow(Conv[part],2) - Merge2/2.*BadPass * n2photon2pt;
      double endir = sqrt(pow(enphoton,2)/pow(Conv[part],2) + 0.25*pow(BadPass,2)*pow(en2photon2pt,2)* pow(Merge2,2) + (pow(en2photon,2)* pow(1. + (1. - Conv[part])*Conv[part]*Merge1 + Miss + (A*(1. + Miss)* (1. + Merge1 + 2.*Miss))/ (1. + 2.*Miss),2))/pow(Conv[part],4) + (pow(eA,2)*pow(1. + Miss,2)* pow(1. + Merge1 + 2.*Miss,2)* pow(n2photon,2))/ (pow(Conv[part],4)*pow(1. + 2.*Miss,2)) + (pow(eMerge1,2)* pow((1. - Conv[part])*Conv[part] + (A*(1. + Miss))/(1. + 2.*Miss),2)* pow(n2photon,2))/pow(Conv[part],4) + (pow(eMiss,2)* pow(1 + (2.*A*(1. + Miss))/ (1. + 2.*Miss) - (2.*A*(1. + Miss)* (1. + Merge1 + 2.*Miss))/ pow(1. + 2.*Miss,2) + (A*(1. + Merge1 + 2.*Miss))/ (1. + 2.*Miss),2)*pow(n2photon,2))/ pow(Conv[part],4) + 0.25*pow(BadPass,2)*pow(eMerge2,2)* pow(n2photon2pt,2) + 0.25*pow(eBadPass,2)*pow(Merge2,2)* pow(n2photon2pt,2) + pow(eConv[part],2)*pow(-((((1. - Conv[part])* Merge1 - Conv[part]*Merge1)*n2photon)/ pow(Conv[part],2)) + (2*(1. + (1. - Conv[part])*Conv[part]*Merge1 + Miss + (A*(1. + Miss)* (1. + Merge1 + 2.*Miss))/ (1. + 2.*Miss))*n2photon)/ pow(Conv[part],3) - nphoton/pow(Conv[part],2),2));

      if(ipt >= 22)  // >14GeV use ERT_4x4b
      {
        ndir *= Norm[part];
        endir *= Norm[part] * (1+eNorm[part]);
      }

      yy[part] = (XBBC/NBBC) / (2*PI*xpt) / (pTbin[ipt+1]-pTbin[ipt]) / DeltaEta
        * ndir / Acc / TrigERT / Prob[part/2][ipt]
        / ToF[part] / TrigBBC * Pile[part];
      eyy[part] = yy[part] * sqrt( pow(endir/ndir,2)
          + pow(eAcc/Acc,2)
          + pow(eTrigERT/TrigERT,2)
          //+ pow(eProb/Prob[part/2][ipt],2)
          //+ pow(eToF[part]/ToF[part],2)
          //+ pow(eTrigBBC/TrigBBC,2) + pow(ePile/Pile[part],2) + pow(eXBBC/XBBC,2)
          );
      if( TMath::Finite(yy[part] + eyy[part]) )
        qt_cross->Fill(ipt, part, xpt, yy[part], eyy[part]);
    } // part

    double ybar, eybar;
    Chi2Fit(3, yy, eyy, ybar, eybar);
    if( TMath::Finite(ybar + eybar) )
      qt_cross->Fill(ipt, 3, xpt, ybar, eybar);
  } // ipt

  mc(6, 2,1);
  legi(0, 0.4,0.7,0.7,0.9);

  for(int part=0; part<4; part++)
  {
    TGraphErrors *gr = qt_cross->Graph(part);
    mcd(6, part/3+1);
    gPad->SetLogy();
    if(part == 0)
      gr->SetTitle("Separated");
    else if(part == 3)
      gr->SetTitle("Combined");
    aset(gr, "p_{T} [GeV]","Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{-3}]", 0.,30., 1e-1,1e6);
    style(gr, part+20, part+1);
    if(part%3==0)
      gr->Draw("AP");
    else
      gr->Draw("P");
    if(part<3)
      leg0->AddEntry(gr, pname[part], "P");
  }

  mcd(6, 1);
  leg0->Draw();
  mcd(6, 2);
  c6->Print("plots/CrossSection-photon.pdf");

  qt_cross->Write();
  for(int part=0; part<3; part++)
    mcw( part, Form("Minv-2photon-part%d",part) );
  qt_cross->Close();

  mc(7);
  mcd(7);
  for(int part=0; part<3; part++)
  {
    TGraphErrors *gr = qt_asee->Graph(part);
    aset(gr, "p_{T} [GeV]","A'", 6.,30.);
    style(gr, part+20, part+1);
    if(part==0)
      gr->Draw("AP");
    else
      gr->Draw("P");
  }
  leg0->Draw();
  c7->Print("plots/ASee.pdf");
  qt_asee->Save();

}
