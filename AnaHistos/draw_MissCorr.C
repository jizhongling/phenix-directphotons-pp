#include "GlobalVars.h"
#include "ReadGraph.h"

void draw_MissCorr()
{
  TFile *f_out = new TFile("data/MissCorr.root", "RECREATE");
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

      double aMissPass = aMiss + aMerge * aBadPass;
      double eaMissPass = sqrt( eaMiss*eaMiss + pow(eaMerge*aBadPass,2.) + pow(aMerge*eaBadPass,2.) );
      double aMissAll = aMiss + aMerge;
      double eaMissAll = sqrt( eaMiss*eaMiss + eaMerge*eaMerge );
      double aMissCorr = (1. + aMissPass) + A * (1. + aMissAll);
      double eaMissCorr = sqrt( pow(eA*(1.+aMissPass),2.) + pow((1.+A)*eaMissPass,2.) );

      gr[part]->SetPoint(igp[part], xpT, aMissCorr);
      gr[part]->SetPointError(igp[part], 0., eaMissCorr);
      igp[part]++;
    }

    gr[part]->Set(igp[part]);
    f_out->cd();
    gr[part]->Write();

    gr[part]->SetTitle("Missing ratio correction");
    aset(gr[part], "p_{T} [GeV]","R_{corr}", 0.,30., 0.,5.);
    style(gr[part], part+20, part+1);
    if(part==0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");
  }

  c0->Print("plots/MissCorr.pdf");
  f_out->Close();
}
