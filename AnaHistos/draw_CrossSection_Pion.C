#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "Chi2Fit.h"

void draw_CrossSection_Pion()
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
  const double BR = 0.988;
  const double Pile[3] = {0.905, 0.905, 0.865};
  const double IsoPile[3] = {1.11, 1.11, 1.08};
  const double ePile = 0.01;
  const double TrigBBC = 0.91;
  const double eTrigBBC = 0.01;
  const double ToF[3] = {0.985, 0.985, 0.995};
  const double eToF[3] = {0.003, 0.003, 0.003};
  const double Conv[3] = {0.720, 0.919, 0.919};
  const double eConv[3] = {0.046, 0.044, 0.044};
  const double Norm[3] = {0.323, 0.323, 0.253};
  const double eNorm[3] = {0.005, 0.006, 0.004};

  const double xpt_shift[npT] =  { 0.25, 0.75, 1.215, 1.719, 2.223, 2.725, 3.228, 3.730, 4.231, 4.732, 5.234, 5.735, 6.237, 6.738, 7.238, 7.739, 8.240, 8.740, 9.241, 9.741, 10.88, 12.90, 14.91, 16.92, 18.93, 20.94, 22.94, 24.95, 27, 29 };

  const double Prob[2][npT] = {
    { 1, 1, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92 },
    { 1, 1, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.93, 0.94, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95 }
  };
  const double eProb = 0.02;

  const double acc_gl[npT] = { 0.01, 0.01, 0.0089, 0.0156, 0.0195, 0.0234, 0.0255, 0.0270, 0.0279, 0.0287, 0.0296, 0.0301, 0.0309, 0.0311, 0.0320, 0.0326, 0.0329, 0.0333, 0.0339, 0.0340, 0.0346, 0.0360, 0.0374, 0.0383, 0.0397, 0.0406, 0.0419, 0.0423, 0.0431, 0.0434, };
  const double emerge_gl[npT] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.00, 0.993, 0.959, 0.881, 0.751, 0.623, 0.481, 0.366, 0.261 };
  const double eff_22_gl[npT] =  { 0.000, 0.000, 0.000, 0.001, 0.002, 0.009, 0.029, 0.066, 0.128, 0.210, 0.289, 0.352, 0.433, 0.476, 0.517, 0.563, 0.600, 0.630, 0.650, 0.660, 0.676, 0.676, 0.676, 0.676, 0.676, 0.676, 0.676, 0.676, 0.676, 0.676, }; // From fit (only bins 8-10 GeV tuned)

  QueryTree *qt_cross = new QueryTree("data/CrossSection-pion.root", "RECREATE");

  QueryTree *qt_acc = new QueryTree("data/Acceptance-pion.root");
  QueryTree *qt_ert = new QueryTree("data/ERTEff-pion.root");
  QueryTree *qt_merge = new QueryTree("data/Merge.root");

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[evtype][part]
  TH2 *h2_pion[3][3];
  TH2 *h2_isopion[3][3];

  int bbc10cm = 1;
  int tof = 1;
  int prob = 1;
  int checkmap = 1;
  int ival = 1;

  TH2 *h2_pion_t = (TH2*)f->Get("h2_pion_0");
  h2_pion_t = (TH2*)h2_pion_t->Clone();
  h2_pion_t->Reset();

  for(int evtype=0; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h2_pion[evtype][part] = (TH2*)h2_pion_t->Clone(Form("h2_pion_type%d_part%d",evtype,part));
      h2_isopion[evtype][part] = (TH2*)h2_pion_t->Clone(Form("h2_isopion_type%d_part%d",evtype,part));

      for(int evenodd=0; evenodd<2; evenodd++)
        for(int pattern=0; pattern<3; pattern++)
          for(int isolated=0; isolated<2; isolated++)
          {
            int ih = part + 3*evenodd + 3*2*pattern + 3*2*3*evtype + 3*2*3*4*tof + 3*2*3*4*2*prob + 3*2*3*4*2*2*bbc10cm + 3*2*3*4*2*2*2*checkmap + 3*2*3*4*2*2*2*2*isolated + 3*2*3*4*2*2*2*2*2*ival;
            TH2 *h2_tmp = (TH2*)f->Get(Form("h2_pion_%d",ih));
            h2_pion[evtype][part]->Add(h2_tmp);
            if(isolated == 1)
              h2_isopion[evtype][part]->Add(h2_tmp);
            delete h2_tmp;
          } // isolated
    }

  for(int part=0; part<3; part++)
  {
    mc(part, 6,5);
    mc(part+3, 6,5);
  }

  for(int ipt=0; ipt<npT; ipt++)
  {
    double xpt, yy[3], eyy[3], yyiso[3], eyyiso[3];

    for(int part=0; part<3; part++)
    {
      int evtype = 2;
      if(ipt < 22)  // <14GeV use ERT_4x4c
        evtype = 2;
      else  // >14GeV use ERT_4x4b
        evtype = 1;

      TH1 *h_minv;

      mcd(part, ipt+1);
      double npion = 1., enpion = 1.;
      h_minv = (TH1*)h2_pion[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, npion, enpion, true, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, npion, enpion, true, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, npion, enpion, false, 0.10,0.17);
      //cout << "Part " << part << ", pT = " << (pTbin[ipt]+pTbin[ipt+1])/2. << ": " << npion << endl;
      if(ipt >= 22)  // >14GeV use ERT_4x4b
      {
        npion *= Norm[part];
        enpion *= Norm[part];
      }
      delete h_minv;

      mcd(part+3, ipt+1);
      double nisopion = 1., enisopion = 1.;
      h_minv = (TH1*)h2_isopion[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, nisopion, enisopion, true, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, nisopion, enisopion, true, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, nisopion, enisopion, false, 0.10,0.17);
      if(ipt >= 22)  // >14GeV use ERT_4x4b
      {
        nisopion *= Norm[part];
        enisopion *= Norm[part];
      }
      delete h_minv;

      double Acc, eAcc, TrigERT, eTrigERT, Merge, eMerge;
      qt_acc->Query(ipt, part, xpt, Acc, eAcc);
      qt_ert->Query(ipt, part/2, xpt, TrigERT, eTrigERT);
      qt_merge->Query(ipt, part/2, xpt, Merge, eMerge);

      if(ipt >= 20)
      {
        if(part < 2)
        {
          TrigERT = 0.952;
          eTrigERT = 0.003;
        }
        else
        {
          TrigERT = 0.661;
          eTrigERT = 0.008;
        }
      }

      xpt = xpt_shift[ipt];

      yy[part] = (XBBC/NBBC) / (2*PI*xpt) / (pTbin[ipt+1]-pTbin[ipt]) / DeltaEta
        * npion / BR / bck[part/2][ipt] / meff[part/2][ipt]
        / Acc / Merge / TrigERT / Prob[part/2][ipt]
        / ToF[part] / Conv[part] / TrigBBC * Pile[part];
      eyy[part] = yy[part] * sqrt( pow(enpion/npion,2)
          + pow(eAcc/Acc,2)
          + pow(eMerge/Merge,2)
          + pow(eTrigERT/TrigERT,2)
          + pow(eProb/Prob[part/2][ipt],2)
          + pow(eToF[part]/ToF[part],2) + pow(eConv[part]/Conv[part],2)
          //+ pow(eTrigBBC/TrigBBC,2) + pow(ePile/Pile[part],2) + pow(eXBBC/XBBC,2)
          );
      if( TMath::Finite(yy[part]+eyy[part]) )
        qt_cross->Fill(ipt, part, xpt, yy[part], eyy[part]);

      yyiso[part] = (XBBC/NBBC) / (2*PI*xpt) / (pTbin[ipt+1]-pTbin[ipt]) / DeltaEta
        * nisopion / BR / bck[part/2][ipt] / meff[part/2][ipt]
        / Acc / Merge / TrigERT / Prob[part/2][ipt]
        / ToF[part] / Conv[part] / TrigBBC * IsoPile[part];
      eyyiso[part] = yyiso[part] * sqrt( pow(enpion/npion,2)
          + pow(eAcc/Acc,2)
          + pow(eMerge/Merge,2)
          + pow(eTrigERT/TrigERT,2)
          + pow(eProb/Prob[part/2][ipt],2)
          + pow(eToF[part]/ToF[part],2) + pow(eConv[part]/Conv[part],2)
          //+ pow(eTrigBBC/TrigBBC,2) + pow(ePile/IsoPile[part],2) + pow(eXBBC/XBBC,2)
          );
      if( TMath::Finite(yyiso[part]+eyyiso[part]) )
        qt_cross->Fill(ipt, part+4, xpt, yyiso[part], eyyiso[part]);
    } // part

    double ybar, eybar;
    if(ipt < 25)
      Chi2Fit(3, yy, eyy, ybar, eybar);
    else
    {
      ybar = yy[2];
      eybar = eyy[2];
    }
    if( TMath::Finite(ybar+eybar) )
      qt_cross->Fill(ipt, part, xpt, ybar, eybar);

    double yisobar, eyisobar;
    Chi2Fit(3, yyiso, eyyiso, yisobar, eyisobar);
    if( TMath::Finite(yisobar+eyisobar) )
      qt_cross->Fill(ipt, 7, xpt, yisobar, eyisobar);
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
    aset(gr, "p_{T} [GeV]", "Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{-3}]", 6.,30., 1e-1,1e5);
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
  gr_sasha->Draw("L");
  c6->Print("plots/CrossSection-pion.pdf");

  qt_cross->Write();
  for(int part=0; part<3; part++)
  {
    mcw( part, Form("Minv-part%d",part) );
    mcw( part+3, Form("Minv-iso-part%d",part) );
  }
  qt_cross->Close();
}
