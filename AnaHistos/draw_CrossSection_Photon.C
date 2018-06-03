#include "GlobalVars.h"
#include "ReadGraph.h"
#include "FitMinv.h"
#include "Chi2Fit.h"

void draw_CrossSection_Photon()
{
  const Double_t PI = TMath::Pi();

  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};

  TGraphErrors *gr[4];  // PbScW, PbScE, PbGl, Combined
  Int_t igp[4] = {};
  for(Int_t part=0; part<4; part++)
  {
    gr[part] = new TGraphErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  TAxis *axis_1sec = hn_1photon->GetAxis(0);
  TAxis *axis_1pt = hn_1photon->GetAxis(1);
  TAxis *axis_1cut = hn_1photon->GetAxis(3);
  TAxis *axis_1type = hn_1photon->GetAxis(4);

  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");
  TAxis *axis_sec = hn_2photon->GetAxis(0);
  TAxis *axis_pt = hn_2photon->GetAxis(1);
  TAxis *axis_minv = hn_2photon->GetAxis(2);
  TAxis *axis_cut = hn_2photon->GetAxis(4);
  TAxis *axis_type = hn_2photon->GetAxis(5);

  const Double_t DeltaEta = 1.0;
  //const Double_t NBBC =  3.59e11;  // from DAQ
  const Double_t NBBC =  3.54e11;  // from rejection power
  const Double_t XBBC = 32.51e9;
  const Double_t eXBBC = 3.24e9;
  const Double_t Pile[3] = {0.891, 0.891, 0.878};
  const Double_t ePile = 0.01;
  const Double_t TrigBBC = 0.91;
  const Double_t eTrigBBC = 0.01;
  const Double_t ToF[3] = {0.992, 0.992, 0.997};
  const Double_t eToF[3] = {0.002, 0.002, 0.002};
  const Double_t Conv[3] = {0.849, 0.959, 0.959};
  const Double_t eConv[3] = {0.027, 0.023, 0.023};
  const Double_t Norm[3] = {0.321, 0.314, 0.243};
  const Double_t eNorm[3] = {0.005, 0.006, 0.005};

  Double_t xAcc[3][npT] = {}, Acc[3][npT] = {}, eAcc[3][npT] = {};
  Double_t xMiss[3][npT] = {}, Miss[3][npT] = {}, eMiss[3][npT] = {};
  Double_t xTrigERT[3][npT] = {}, TrigERT[3][npT] = {}, eTrigERT[3][npT] = {};

  Double_t bck[2][npT] = {
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.10, 1.15, 1.20, 1.30,  1, 1, 1 },
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.08, 1.08, 1.11, 1.11, 1.11, 1.11, 1.11 }
  };

  Double_t meff[2][npT] = {
    { 1, 1, 0.96, 0.97, 0.98, 0.985, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.985, 0.995, 0.995, 0.99, 0.98, 0.95, 1,1,1,1,1, },
    { 1, 1, 0.95, 0.97, 0.975, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.995, 0.995, 0.99, 0.99, 0.98, 0.98, 0.98, 0.98, 1, 1 }
  };

  Double_t Prob[2][npT] = {
    { 1, 1, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 },
    { 1, 1, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97 }
  };
  Double_t eProb = 0.02;

  for(Int_t part=0; part<3; part++)
  {
    ReadGraph<TGraphAsymmErrors>("data/Acceptance-photon.root", part, xAcc[part], Acc[part], eAcc[part]);
    ReadGraph<TGraphAsymmErrors>("data/HadronRatio.root", part, xMiss[part], Miss[part], eMiss[part]);
    ReadGraph<TGraphAsymmErrors>("data/ERTEff-photon.root", part/2, xTrigERT[part], TrigERT[part], eTrigERT[part]);
    mc(part, 6,5);
  }
  TrigERT[0][0] = TrigERT[1][0] = 0.948;
  eTrigERT[0][0] = eTrigERT[1][0] = 0.004;
  TrigERT[2][0] = 0.651;
  eTrigERT[2][0] = 0.008;

  for(Int_t ipt=0; ipt<npT; ipt++)
  {
    Double_t xx, yy[3], eyy[3];

    for(Int_t part=0; part<3; part++)
    {
      if(ipt < 22)  // <14GeV use ERT_4x4c
      {
        axis_1type->SetRange(3,3);
        axis_type->SetRange(3,3);
      }
      else  // >14GeV use ERT_4x4b
      {
        axis_1type->SetRange(2,2);
        axis_type->SetRange(2,2);
      }
      axis_1cut->SetRange(4,4);
      axis_1pt->SetRange(ipt+1,ipt+1);
      axis_cut->SetRange(4,4);
      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);

      TH1 *h_photon = hn_1photon->Projection(0);
      Double_t nphoton = h_photon->Integral(secl[part],sech[part]);
      delete h_photon;

      mcd(part, ipt+1);
      Double_t npion = 1., enpion = 1.;
      TH1 *h_minv = hn_2photon->Projection(2);
      h_minv->Scale(0.5);
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, npion, enpion, kTRUE, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, npion, enpion, kTRUE, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, npion, enpion, kFALSE, 0.10,0.17);
      //cout << "Part " << part << ", pT = " << (pTbin[ipt]+pTbin[ipt+1])/2. << ": " << npion << endl;
      npion /= bck[part/2][ipt] * meff[part/2][ipt];
      delete h_minv;

      xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Int_t ipAcc = Get_ipt(xAcc[part], xx);
      Int_t ipMiss = Get_ipt(xMiss[part], xx);
      Int_t ipTrigERT = Get_ipt(xTrigERT[part], xx);
      if(ipt >= 20)
        ipTrigERT = 0;
      Double_t ndir = nphoton - Miss[part][ipMiss]*npion*2.;
      Double_t endir = sqrt( nphoton + pow(eMiss[part][ipMiss]*npion*2.,2.) + pow(Miss[part][ipMiss]*enpion*2.,2.) );
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
          + pow(eMiss[part][ipMiss]/Miss[part][ipMiss],2.)
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

    Double_t ybar, eybar;
    if(ipt < 25)
      Chi2Fit(3, yy, eyy, ybar, eybar);
    else
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

  for(Int_t part=0; part<4; part++)
  {
    gr[part]->Set(igp[part]);
    mcd(3, part/3+1);
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
  mcd(3, 1);
  leg0->Draw();
  mcd(3, 2);
  c3->Print("plots/CrossSection-photon.pdf");

  TFile *f_out = new TFile("data/CrossSection-photon.root", "RECREATE");
  for(Int_t part=0; part<4; part++)
  {
    if(part<3)
      mcw( part, Form("Minv-part%d",part) );
    gr[part]->Write();
  }
  f_out->Close();
}
