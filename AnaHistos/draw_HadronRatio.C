#include "GlobalVars.h"
#include "ReadGraph.h"

void draw_HadronRatio()
{
  TFile *f_out = new TFile("data/HadronRatio.root", "RECREATE");
  TGraphErrors *gr[3];
  int igp[3] = {};

  const double A = 0.22;
  const double eA = 0.04;
  double xMiss[3][npT] = {}, Miss[3][npT] = {}, eMiss[3][npT] = {};
  double xMerge[3][npT] = {}, Merge[3][npT] = {}, eMerge[3][npT] = {};
  double xBadPass[3][npT] = {}, BadPass[3][npT] = {}, eBadPass[3][npT] = {};

  mc();
  mcd();

  for(int part=0; part<3; part++)
  {
    gr[part] = new TGraphErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));

    ReadGraph<TGraphErrors>("data/MissingRatio.root", part, xMiss[part], Miss[part], eMiss[part]);
    ReadGraph<TGraphAsymmErrors>("data/Merge-photon.root", part, xMerge[part], Merge[part], eMerge[part]);
    ReadGraph<TGraphErrors>("data/MergePassRate.root", part/2, xBadPass[part], BadPass[part], eBadPass[part]);

    for(int ipt=0; ipt<30; ipt++)
    {
      double xpT = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      int ipMiss = Get_ipt(xMiss[part], xpT);
      double aMiss = Miss[part][ipMiss];
      double eaMiss = eMiss[part][ipMiss];
      int ipMerge = Get_ipt(xMerge[part], xpT);
      double aMerge = Merge[part][ipMerge];
      double eaMerge = eMerge[part][ipMerge];
      int ipBadPass = Get_ipt(xBadPass[part], xpT);
      double aBadPass = BadPass[part][ipBadPass];
      double eaBadPass = eBadPass[part][ipBadPass];
      if(ipt<23)
      {
        aBadPass = 0.;
        eaBadPass = 0.;
      }

      double aMissCorr = aMiss + aMerge * aBadPass;
      double eaMissCorr = sqrt( eaMiss*eaMiss + pow(eaMerge*aBadPass,2.) + pow(eaBadPass*aMerge,2.) );
      double aHadronR = (1. + aMissCorr) * (1. + A);
      double eaHadronR = sqrt( pow(eaMissCorr*(1.+A),2.) + pow(eA*(1.+aMissCorr),2.) );

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
