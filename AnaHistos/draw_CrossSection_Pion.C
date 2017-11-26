#include "GlobalVars.h"
#include "ReadGraph.h"
#include "FitMinv.h"
#include "Chi2Fit.h"

void draw_CrossSection_Pion()
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

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  TAxis *axis_sec = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);
  TAxis *axis_cond = hn_pion->GetAxis(3);

  const Double_t DeltaEta = 1.0;
  const Double_t NBBC =  3.62e11;
  const Double_t XBBC = 32.51e9;
  const Double_t eXBBC = 3.24e9;
  const Double_t Pile = 0.94;
  const Double_t ePile = 0.01;
  const Double_t TrigBBC = 0.91;
  const Double_t eTrigBBC = 0.01;
  const Double_t ToF[3] = {0.984, 0.984, 0.997};
  const Double_t eToF[3] = {0.003, 0.003, 0.007};
  const Double_t Conv[3] = {0.760, 1., 1.};
  const Double_t eConv[3] = {0.033, 0., 0.};

  Double_t xAcc[3][npT] = {}, Acc[3][npT] = {}, eAcc[3][npT] = {};
  Double_t xMerge[3][npT] = {}, Merge[3][npT] = {}, eMerge[3][npT] = {};
  Double_t xTrigERT[3][npT] = {}, TrigERT[3][npT] = {}, eTrigERT[3][npT] = {};
  Double_t xProb[3][npT] = {}, Prob[3][npT] = {}, eProb[3][npT] = {};

  for(Int_t part=0; part<3; part++)
  {
    ReadGraph<TGraphAsymmErrors>("data/Acceptance.root", part, xAcc[part], Acc[part], eAcc[part]);
    ReadGraph<TGraphAsymmErrors>("data/Merge.root", part/2, xMerge[part], Merge[part], eMerge[part]);
    ReadGraph<TGraphAsymmErrors>("data/ERTEff.root", part/2, xTrigERT[part], TrigERT[part], eTrigERT[part]);
    ReadGraph<TGraphAsymmErrors>("data/ProbEff.root", part/2, xProb[part], Prob[part], eProb[part]);
    mc(part, 6,5);
  }

  for(Int_t ipt=0; ipt<npT; ipt++)
  {
    Double_t xx, yy[3], eyy[3];

    for(Int_t part=0; part<3; part++)
    {
      axis_cond->SetRange(2,2);
      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);

      mcd(part, ipt+1);
      Double_t npion = 1., enpion = 1.;
      TH1 *h_minv = hn_pion->Projection(2);
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T} %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      FitMinv(h_minv, npion, enpion);
      delete h_minv;

      xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Int_t ipAcc = Get_ipt(xAcc[part], xx);
      Int_t ipMerge = Get_ipt(xMerge[part], xx);
      Int_t ipTrigERT = Get_ipt(xTrigERT[part], xx);
      Int_t ipProb = Get_ipt(xProb[part], xx);
      yy[part] = npion * (XBBC/NBBC) / (2*PI*xx) / (pTbin[ipt+1]-pTbin[ipt]) / DeltaEta
        / Acc[part][ipAcc] / Merge[part][ipMerge]
        / TrigERT[part][ipTrigERT] / Prob[part][ipProb]
        / ToF[part] / Conv[part] / TrigBBC * Pile;
      eyy[part] = yy[part] * sqrt( pow(enpion/npion,2.)
          + pow(eAcc[part][ipAcc]/Acc[part][ipAcc],2.)
          + pow(eMerge[part][ipMerge]/Merge[part][ipMerge],2.)
          + pow(eTrigERT[part][ipTrigERT]/TrigERT[part][ipTrigERT],2.)
          + pow(eProb[part][ipProb]/Prob[part][ipProb],2.)
          + pow(eToF[part]/ToF[part],2.) + pow(eConv[part]/Conv[part],2.)
          //+ pow(eTrigBBC/TrigBBC,2.) + pow(ePile/Pile,2.) + pow(eXBBC/XBBC,2.)
          );
      if( yy[part] > 0. && eyy[part] > 0. && eyy[part] < TMath::Infinity() )
      {
        gr[part]->SetPoint(igp[part], xx, yy[part]);
        gr[part]->SetPointError(igp[part], 0., eyy[part]);
        igp[part]++;
      }
    } // part

    Double_t ybar, eybar;
    Chi2Fit(3, yy, eyy, ybar, eybar);
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
    aset(gr[part], "p_{T} [GeV]", "Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{-3}]", 6.,20., 1.,1e5);
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
  c3->Print("plots/CrossSection-pion.pdf");

  TFile *f_out = new TFile("data/CrossSection-pion.root", "RECREATE");
  for(Int_t part=0; part<4; part++)
  {
    if(part<3)
      mcw( part, Form("Minv-part%d",part) );
    gr[part]->Write();
  }
  f_out->Close();
}
