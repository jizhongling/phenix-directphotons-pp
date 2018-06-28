#include "GlobalVars.h"
#include "ReadGraph.h"
#include "FitMinv.h"
#include "Chi2Fit.h"

void draw_CrossSection_Pion()
{
  const double PI = TMath::Pi();

  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  TGraphErrors *gr[4];  // PbScW, PbScE, PbGl, Combined
  int igp[4] = {};
  for(int part=0; part<4; part++)
  {
    gr[part] = new TGraphErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  TAxis *axis_sec = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);
  TAxis *axis_cut = hn_pion->GetAxis(3);
  TAxis *axis_type = hn_pion->GetAxis(4);

  const double DeltaEta = 1.0;
  //const double NBBC =  3.59e11;  // from DAQ
  const double NBBC =  3.54e11;  // from rejection power
  const double XBBC = 32.51e9;
  const double eXBBC = 3.24e9;
  const double BR = 0.988;
  //const double Pile = 0.94;
  const double Pile[3] = {0.905, 0.905, 0.865};
  const double ePile = 0.01;
  const double TrigBBC = 0.91;
  const double eTrigBBC = 0.01;
  const double ToF[3] = {0.985, 0.985, 0.995};
  const double eToF[3] = {0.003, 0.003, 0.003};
  //const double Conv[3] = {0.784, 0.983, 0.983};
  //const double eConv[3] = {0.033, 0.017, 0.017};
  const double Conv[3] = {0.720, 0.919, 0.919};
  const double eConv[3] = {0.046, 0.044, 0.044};
  const double Norm[3] = {0.326, 0.326, 0.251};
  const double eNorm[3] = {0.001, 0.001, 0.004};

  double xAcc[3][npT] = {}, Acc[3][npT] = {}, eAcc[3][npT] = {};
  double xMerge[3][npT] = {}, Merge[3][npT] = {}, eMerge[3][npT] = {};
  double xTrigERT[3][npT] = {}, TrigERT[3][npT] = {}, eTrigERT[3][npT] = {};

  double xpt[npT] =  { 0.25, 0.75, 1.215, 1.719, 2.223, 2.725, 3.228, 3.730, 4.231, 4.732, 5.234, 5.735, 6.237, 6.738, 7.238, 7.739, 8.240, 8.740, 9.241, 9.741, 10.88, 12.90, 14.91, 16.92, 18.93, 20.94, 22.94, 24.95, 27, 29 };

  double bck[2][npT] = {
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.10, 1.15, 1.20, 1.30,  1, 1, 1 },
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.08, 1.08, 1.11, 1.11, 1.11, 1.11, 1.11 }
  };

  double meff[2][npT] = {
    { 1, 1, 0.96, 0.97, 0.98, 0.985, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.985, 0.995, 0.995, 0.99, 0.98, 0.95, 1,1,1,1,1, },
    { 1, 1, 0.95, 0.97, 0.975, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.995, 0.995, 0.99, 0.99, 0.98, 0.98, 0.98, 0.98, 1, 1 }
  };

  double Prob[2][npT] = {
    { 1, 1, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92 },
    { 1, 1, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.93, 0.94, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95 }
  };
  double eProb = 0.02;

  double acc_gl[npT] = { 0.01, 0.01, 0.0089, 0.0156, 0.0195, 0.0234, 0.0255, 0.0270, 0.0279, 0.0287, 0.0296, 0.0301, 0.0309, 0.0311, 0.0320, 0.0326, 0.0329, 0.0333, 0.0339, 0.0340, 0.0346, 0.0360, 0.0374, 0.0383, 0.0397, 0.0406, 0.0419, 0.0423, 0.0431, 0.0434, };
  double emerge_gl[npT] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.00, 0.993, 0.959, 0.881, 0.751, 0.623, 0.481, 0.366, 0.261 };
  double eff_22_gl[npT] =  { 0.000, 0.000, 0.000, 0.001, 0.002, 0.009, 0.029, 0.066, 0.128, 0.210, 0.289, 0.352, 0.433, 0.476, 0.517, 0.563, 0.600, 0.630, 0.650, 0.660, 0.676, 0.676, 0.676, 0.676, 0.676, 0.676, 0.676, 0.676, 0.676, 0.676, }; // From fit (only bins 8-10 GeV tuned)

  for(int part=0; part<3; part++)
  {
    ReadGraph<TGraphAsymmErrors>("data/Acceptance-pion.root", part, xAcc[part], Acc[part], eAcc[part]);
    ReadGraph<TGraphAsymmErrors>("data/Merge.root", part/2, xMerge[part], Merge[part], eMerge[part]);
    ReadGraph<TGraphAsymmErrors>("data/ERTEff-pion.root", part/2, xTrigERT[part], TrigERT[part], eTrigERT[part]);
    mc(part, 6,5);
  }
  TrigERT[0][0] = TrigERT[1][0] = 0.951;
  eTrigERT[0][0] = eTrigERT[1][0] = 0.004;
  TrigERT[2][0] = 0.666;
  eTrigERT[2][0] = 0.009;

  for(int ipt=0; ipt<npT; ipt++)
  {
    double xx, yy[3], eyy[3];

    for(int part=0; part<3; part++)
    {
      if(ipt < 22)  // <14GeV use ERT_4x4c
        axis_type->SetRange(3,3);
      else  // >14GeV use ERT_4x4b
        axis_type->SetRange(2,2);
      axis_cut->SetRange(4,4);
      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);

      mcd(part, ipt+1);
      double npion = 1., enpion = 1.;
      TH1 *h_minv = hn_pion->Projection(2);
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
        enpion *= Norm[part] * (1+eNorm[part]);
      }
      delete h_minv;


      //xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      xx = xpt[ipt];
      int ipAcc = Get_ipt(xAcc[part], xx);
      int ipMerge = Get_ipt(xMerge[part], xx);
      int ipTrigERT = Get_ipt(xTrigERT[part], xx);
      if(ipt >= 20)
        ipTrigERT = 0;
      if(part==3)
      {
        Acc[part][ipt] = acc_gl[ipt];
        ipAcc = ipt;
        Merge[part][ipt] = emerge_gl[ipt];
        ipMerge = ipt;
        TrigERT[part][ipt] = eff_22_gl[ipt];
        ipTrigERT = ipt;
      }
      yy[part] = (XBBC/NBBC) / (2*PI*xx) / (pTbin[ipt+1]-pTbin[ipt]) / DeltaEta
        * npion / BR / bck[part/2][ipt] / meff[part/2][ipt]
        / Acc[part][ipAcc] / Merge[part][ipMerge]
        / TrigERT[part][ipTrigERT] / Prob[part/2][ipt]
        / ToF[part] / Conv[part] / TrigBBC * Pile[part];
      eyy[part] = yy[part] * sqrt( pow(enpion/npion,2.)
          + pow(eAcc[part][ipAcc]/Acc[part][ipAcc],2.)
          + pow(eMerge[part][ipMerge]/Merge[part][ipMerge],2.)
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
    //if(ipt < 25)
    //  Chi2Fit(3, yy, eyy, ybar, eybar);
    //else
    {
      ybar = yy[2];
      eybar = eyy[2];
    }
    if( ybar > 0. && eybar > 0. && eybar < TMath::Infinity() )
    {
      gr[3]->SetPoint(igp[3], xx, ybar);
      gr[3]->SetPointError(igp[3], 0., eybar);
      igp[3]++;
    }
  } // ipt

  mc(3, 2,1);
  legi(0, 0.4,0.7,0.7,0.9);

  for(int part=0; part<4; part++)
  {
    gr[part]->Set(igp[part]);
    mcd(3, part/3+1);
    gPad->SetLogy();
    aset(gr[part], "p_{T} [GeV]", "Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{-3}]", 6.,30., 1e-1,1e5);
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
  mcd(3, 1);
  leg0->Draw();
  mcd(3, 2);
  gr_sasha->Draw("L");
  c3->Print("plots/CrossSection-pion.pdf");

  TFile *f_out = new TFile("data/CrossSection-pion.root", "RECREATE");
  for(int part=0; part<4; part++)
  {
    if(part<3)
      mcw( part, Form("Minv-part%d",part) );
    gr[part]->Write();
  }
  f_out->Close();
}
