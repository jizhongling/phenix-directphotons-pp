#include "GlobalVars.h"
#include "FitMinv.h"

void GenerateGraph(TFile *f, TObjArray *Glist, Int_t part)
{
  const Double_t PI = TMath::Pi();

  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};

  TGraphErrors *gr = new TGraphErrors(npT);
  Int_t igp = 0;

  const Double_t NBBC =  3.62e11;
  const Double_t XBBC = 32.5e9;
  const Double_t eXBBC = 3.25e9;
  const Double_t TrigBBC = 0.91;
  const Double_t eTrigBBC = 0.01;
  const Double_t Conv[3] = {0.760, 0., 0.};
  const Double_t eConv[3] = {0.033, 0., 0.};
  Double_t Acc[30], TrigERT[30];
  Double_t eAcc[30], eTrigERT[30];

  ReadGraphAsymmErrors("ERTEff.root", part, gx, TrigERT, eTrigERT);
  ReadGraphErrors("Acceptance.root", part, gx, Acc, eAcc);

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");

  TAxis *axis_sec = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);
  TAxis *axis_cond = hn_pion->GetAxis(3);

  axis_cond->SetRange(2,2);
  axis_sec->SetRange(secl[part], sech[part]);

  mc(1, 6,5);

  Int_t ipad = 1;
  for(Int_t ipt=1; ipt<30; ipt++)
  {
    mcd(1, ipad++);

    Double_t npion, enpion;
    axis_pt->SetRange(ipt+1,ipt+1);
    TH1 *h_minv = hn_pion->Projection(2);
    h_minv->Rebin(10);
    h_minv->SetTitle( Form("p_{T} %3.1f-%3.1f", pTbin[ipt], pTbin[ipt+1]) );
    FitMinv(h_minv, npion, enpion);
    delete h_minv;

    xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
    yy = (XBBC/NBBC) * npion / (2*PI*gx[ipt]) / ((pTbin[ipt+1]-pTbin[ipt])*1.0)
      / Acc[ipt] / Conv[ipt] / TrigBBC / TrigERT[ipt];
    eyy = yy[ipt] * sqrt(
        pow(enpion/npion,2.) + pow(eTrigERT[ipt]/TrigERT[ipt],2.)
        + pow(eAcc[ipt]/Acc[ipt],2.) +  pow(eConv[ipt]/Conv[ipt],2.)
        // + pow(eXBBC/XBBC,2.) + pow(eTrigBBC/TrigBBC,2.)
        );
    if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
    {
      gr->SetPoint(igp, xx, yy);
      gr->SetPointError(igp, 0., eyy);
      igp++;
    }
  }

  Glist->AddAtAndExpand(gr,part);
  gr->SetName(Form("gr_%d",part));
  gr->SetTitle("#pi^{0} Cross Section");

  c1->Print(Form("CrossSection-part%d.pdf",part));
  delete c1;

  return;
}

void draw_CrossSection_Pion()
{
  gROOT->ProcessLine(".L ReadGraph.C");
  gROOT->ProcessLine(".L Chi2Fit.C");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");
  TObjArray *Glist = new TObjArray();
  TGraphErrors *gr;
  TLegend *leg;

  for(Int_t part=0; part<3; part++)
    GenerateGraph(f, Glist, part);

  for(Int_t im=0; im<1; im++)
  {
    for(Int_t i=2; i<3; i++)
    {
      Double_t gx[30], gy[30], egy[30];
      Double_t gyy[3][30], egyy[3][30];

      for(Int_t part=0; part<3; part++)
        ReadGraph((TGraphErrors*)Glist->At(36*im+12*trig+4*i+part), gx, gyy[part], egyy[part]);

      for(Int_t ipt=0; ipt<30; ipt++)
      {
        Double_t yy[3];
        Double_t eyy[3];
        for(Int_t part=0; part<3; part++)
        {
          yy[part] = gyy[part][ipt];
          eyy[part] = egyy[part][ipt];
        }
        Chi2Fit(2, yy, eyy, gy[ipt], egy[ipt]);
      }

      gr = new TGraphErrors(30, gx, gy, 0, egy);
      Glist->AddAtAndExpand(gr, 36*im+12*trig+4*i+3);
      gr->SetName(Form("gr_%d",36*im+12*trig+4*i+3));
      gr->SetTitle(((TGraphErrors*)Glist->At(36*im+12*trig+4*i))->GetTitle());
    }
  }

  mc(0, 3,2);
  legi(0, 0.4, 0.7, 0.7, 0.9);

  const char *legname[4] = {"PbScW", "PbScE", "PbGlE", "Combined"};
  for(Int_t im=0; im<1; im++)
  {
    for(Int_t i=2; i<3; i++)
      for(Int_t part=0; part<3; part++)
      {
        mcd(0, i+3*(part/3)+1);
        gr = (TGraphErrors*)Glist->At(36*im+12*trig+4*i+part);
        gr->GetXaxis()->SetTitle("p_{T} [GeV]");
        gr->GetXaxis()->SetRangeUser(0., 30.);
        if(i==0)
        {
          aset(gr, "p_{T} [GeV]", "Ratio", 0.,30., 0.,1.2);
        }
        else
        {
          gPad->SetLogy();
          aset(gr, "p_{T} [GeV]", "Ed^{3}#sigma/dp^{3} [pb*GeV^{-2}*c^{-3}]", 5.,20., 1.,1e6);
        }
        style(gr, part+20, part+1);
        if(part%3 == 0)
        {
          gr->Draw("AP");
          leg0->Draw();
        }
        else
        {
          gr->Draw("P");
        }
        leg0->AddEntry(gr, legname[part], "P");
      }
    c0->Print(Form("CrossSection-method-%d-ert%c.pdf",im,97+trig));
    c0->Clear("D");
  }

  TFile *fout = new TFile(Form("CrossSection-ert%c.root",97+trig), "RECREATE");
  Glist->Write();
  fout->Close();
}
