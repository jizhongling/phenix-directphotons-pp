#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "Chi2Fit.h"
#include "Sasha-cross_sep-get_cross_gamma3.h"

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
  const double Conv[3] = {0.720, 0.943, 0.919};
  const double eConv[3] = {0.046, 0.044, 0.044};
  //const double Norm[3] = {0.328, 0.323, 0.251};  // Sasha's value
  const double Norm[3] = {0.322, 0.320, 0.252};
  const double eNorm[3] = {0.005, 0.006, 0.004};

  QueryTree *qt_cross = new QueryTree("data/CrossSection-pion.root", "RECREATE");

  QueryTree *qt_acc = new QueryTree("data/Acceptance-pion.root");
  QueryTree *qt_merge = new QueryTree("data/Merge.root");
  QueryTree *qt_ert = new QueryTree("data/ERTEff-pion.root");
  QueryTree *qt_pt = new QueryTree("data/PtShift.root");

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-Sasha.root");

  // h[iso][evtype][part]
  TH2 *h2_pion[2][3][3];

  int tof = 1;
  int prob = 1;
  int checkmap = 1;
  int ival = 1;

  TH2 *h2_pion_t = (TH2*)f->Get("h2_pion_0");
  h2_pion_t = (TH2*)h2_pion_t->Clone();
  h2_pion_t->Reset();

  for(int isolated=0; isolated<2; isolated++)
    for(int evtype=0; evtype<3; evtype++)
      for(int part=0; part<3; part++)
      {
        h2_pion[isolated][evtype][part] = (TH2*)h2_pion_t->Clone(Form("h2_%spion_type%d_part%d",isolated?"iso":"",evtype,part));
        int ih = part + 3*evtype + 3*3*tof + 3*3*2*prob + 3*3*2*2*checkmap + 3*3*2*2*2*isolated + 3*3*2*2*2*2*ival;
        TH2 *h2_tmp = (TH2*)f->Get(Form("h2_pion_%d",ih));
        h2_pion[isolated][evtype][part]->Add(h2_tmp);
        if(isolated == 1)
          h2_pion[0][evtype][part]->Add(h2_tmp);
      }

  for(int part=0; part<3; part++)
  {
    mc(part, 6,5);
    mc(part+3, 6,5);
  }

  for(int ipt=0; ipt<npT; ipt++)
  {
    double dummy, xpt, xsec[2][3], exsec[2][3],
           rsys[2][3], ersys[2][3];  // iso, part

    for(int part=0; part<3; part++)
    {
      int evtype = 2;
      if(ipt < 22)  // <14GeV use ERT_4x4c
        evtype = 2;
      else  // >14GeV use ERT_4x4b
        evtype = 1;

      double Acc, eAcc, Merge, eMerge, TrigERT, eTrigERT;
      qt_acc->Query(ipt, part, xpt, Acc, eAcc);
      qt_merge->Query(ipt, part/2, dummy, Merge, eMerge);
      qt_ert->Query(ipt, part/2, dummy, TrigERT, eTrigERT);

      if(xpt > 10.)
      {
        if(part < 2)
        {
          TrigERT = 0.960;
          eTrigERT = 0.002;
        }
        else
        {
          TrigERT = 0.673;
          eTrigERT = 0.005;
        }
      }

      // Use Sasha's values
      if(false)
      {
        Acc = acc[part][ipt];
        Merge = emerge[part/2][ipt];
        TrigERT = eff_22[part/2][ipt];
        if(ipt < 22)
          npion = np_22[part][ipt];
      }

      for(int iso=0; iso<2; iso++)
      {
        TH1 *h_minv;
        mcd(part, ipt+1);
        double npion = 1., enpion = 1.;
        h_minv = (TH1*)h2_pion[iso][evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
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

        qt_pt->Query(ipt, 0, dummy, xpt, dummy);
        xsec[iso][part] = (XBBC/NBBC) / (2*PI*xpt) / (pTbin[ipt+1]-pTbin[ipt]) / DeltaEta
          * npion / BR / bck[part/2][ipt] / meff[part/2][ipt]
          / Acc / Merge / TrigERT / Prob[part/2][ipt]
          / ToF[part] / Conv[part] / TrigBBC * Pile[part];
        exsec[iso][part] = xsec[iso][part]*enpion/npion;
        if( TMath::Finite(xsec[iso][part]+exsec[iso][part]) )
          qt_cross->Fill(ipt, part+4*iso, xpt, xsec[iso][part], exsec[iso][part]);

        const double eFactor = 1e-3;
        rsys[iso][part] = sqrt(
            pow(eAcc/Acc,2)
            + pow(eMerge/Merge,2)
            + pow(eTrigERT/TrigERT,2)
            + pow(eConv[part]/Conv[part],2)
            + pow(eProb/Prob[part/2][ipt],2)
            + pow(eToF[part]/ToF[part],2)
            + pow(ePile/IsoPile[part],2)
            + pow(eTrigBBC/TrigBBC,2)
            );
        ersys[iso][part] = rsys[iso][part]*eFactor*exsec[iso][part]/xsec[iso][part];
      } // iso
    } // part

    for(int iso=0; iso<2; iso++)
    {
      double xsecbar = xsec[iso][2];
      double exsecbar = exsec[iso][2];
      if(xpt < 20.)
        Chi2Fit(3, xsec[iso], exsec[iso], xsecbar, exsecbar);

      double rbar = rsys[iso][2];
      double erbar = ersys[iso][2];
      if(xpt < 20.)
        Chi2Fit(3, rsys[iso], ersys[iso], rbar, erbar);

      double etotal = sqrt(exsecbar*exsecbar + xsecbar*xsecbar*rbar*rbar);
      if( TMath::Finite(etotal) )
        qt_cross->Fill(ipt, 3+4*iso, xpt, xsecbar, exsecbar);
    }
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
    aset(gr, "p_{T} [GeV]", "Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{3}]", 6.,30., 1e-1,1e5);
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
