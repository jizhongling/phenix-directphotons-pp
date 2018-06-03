#include "GlobalVars.h"
#include "ReadGraph.h"

void draw_HadronRatio()
{
  TFile *f_out = new TFile("data/HadronRatio.root", "RECREATE");
  TGraphErrors *gr[3];
  Int_t igp[3] = {};

  const Double_t A = 0.22;
  const Double_t eA = 0.04;
  Double_t xMiss[3][npT] = {}, Miss[3][npT] = {}, eMiss[3][npT] = {};
  Double_t xMerge[3][npT] = {}, Merge[3][npT] = {}, eMerge[3][npT] = {};
  Double_t xBadPass[3][npT] = {}, BadPass[3][npT] = {}, eBadPass[3][npT] = {};

  mc();
  mcd();

  for(Int_t part=0; part<3; part++)
  {
    gr[part] = new TGraphErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));

    ReadGraph<TGraphErrors>("data/MissingRatio.root", part, xMiss[part], Miss[part], eMiss[part]);
    ReadGraph<TGraphAsymmErrors>("data/Merge-photon.root", part, xMerge[part], Merge[part], eMerge[part]);
    ReadGraph<TGraphErrors>("data/MergePassRate.root", part/2, xBadPass[part], BadPass[part], eBadPass[part]);

    for(Int_t ipt=0; ipt<30; ipt++)
    {
      Double_t xpT = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Int_t ipMiss = Get_ipt(xMiss[part], xpT);
      Double_t aMiss = Miss[part][ipMiss];
      Double_t eaMiss = eMiss[part][ipMiss];
      Int_t ipMerge = Get_ipt(xMerge[part], xpT);
      Double_t aMerge = Merge[part][ipMerge];
      Double_t eaMerge = eMerge[part][ipMerge];
      Int_t ipBadPass = Get_ipt(xBadPass[part], xpT);
      Double_t aBadPass = BadPass[part][ipBadPass];
      Double_t eaBadPass = eBadPass[part][ipBadPass];
      if(ipt<23)
      {
        aBadPass = 0.;
        eaBadPass = 0.;
      }

      Double_t aMissCorr = aMiss + aMerge * aBadPass;
      Double_t eaMissCorr = sqrt( eaMiss*eaMiss + pow(eaMerge*aBadPass,2.) + pow(eaBadPass*aMerge,2.) );
      Double_t aHadronR = (1. + aMissCorr) * (1. + A);
      Double_t eaHadronR = sqrt( pow(eaMissCorr*(1.+A),2.) + pow(eA*(1.+aMissCorr),2.) );

      gr[part]->SetPoint(igp[part], xpT, aHadronR);
      gr[part]->SetPointError(igp[part], 0., eaHadronR);
      igp[part]++;
    }

    gr[part]->Set(igp[part]);
    f_out->cd();
    gr[part]->Write();

    gr[part]->SetTitle("Total hadron correction");
    aset(gr[part], "p_{T} [GeV]","(1+R)*(1+A)", 0.,30., 1.,7.);
    style(gr[part], part+20, part+1);
    if(part==0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");
  }

  c0->Print("plots/HadronRatio.pdf");
  f_out->Close();
}
