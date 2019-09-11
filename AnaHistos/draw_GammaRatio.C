#include "GlobalVars.h"
#include "DivideFunctions.h"

void draw_GammaRatio()
{
  const double BR[4] = {0.988, 0.394, 0.0445, 0.021};

  //TFile *f = new TFile("data/GammaRatio-histo.root");
  //TH3 *h3_particles = (TH3*)f->Get("h3_particles");
  //h3_particles->GetZaxis()->SetRange(1,7);
  //TH1* h_others = h3_particles->ProjectionX("h_others", 3,5);
  //TH1* h_pion = h3_particles->ProjectionX("h_pion", 2,2);
  //TGraphErrors *gr = DivideHisto(h_others, h_pion);

  TGraphErrors *gr[2];
  for(int iph=0; iph<2; iph++)
    gr[iph] = new TGraphErrors(npT);
  int igp[2] = {};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-GenPH-histo-minbias.root");
  TH1 *h_photo_eta050 = (TH1*)f->Get("h_photon_eta050");
  THnSparse *hn_hadron = (THnSparse*)f->Get("hn_hadron");
  TH1 *h_hadron[4];
  for(int id=0; id<4; id++)
  {
    hn_hadron->GetAxis(1)->SetRange(id+1,id+1);
    h_hadron[id] = hn_hadron->Projection(0);
  }

  for(int ipt=4; ipt<npT; ipt++)
  {
    double npion = h_hadron[0]->GetBinContent(ipt+1) * BR[0];
    double enpion = h_hadron[0]->GetBinError(ipt+1) * BR[0];
    double nothers = 0.;
    double enothers = 0.;
    for(int id=1; id<4; id++)
    {
      nothers += h_hadron[id]->GetBinContent(ipt+1) * BR[id];
      enothers += h_hadron[id]->GetBinError(ipt+1) * BR[id];
    }

    for(int iph=0; iph<2; iph++)
    {
      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      if(iph == 1)
        nothers = h_photon_eta050->GetBinContent(ipt+1) * 2.;
      double ratio = nothers / npion;
      double eratio = ratio * sqrt( pow(enothers/nothers,2) + pow(enpion/npion,2) );
      if( TMath::Finite(ratio+eratio) )
      {
        gr[iph]->SetPoint(igp[iph], xpt, ratio);
        gr[iph]->SetPointError(igp[iph], 0., eratio);
        igp[iph]++;
      }
    }
  }
  for(int iph=0; iph<2; iph++)
    gr[iph]->Set(igp[iph]);

  mc(0, 2,1);
  for(int iph=0; iph<2; iph++)
  {
    mcd(0, iph+1);
    gr[iph]->SetTitle("Photon Ratio");
    if(iph == 0)
      aset(gr[iph], "p_{T} [GeV]","#frac{#eta+#omega+#eta'}{#pi^{0}}", 5.1,30., 0.,0.5);
    else if(iph == 1)
      aset(gr[iph], "p_{T} [GeV]","#frac{prompt photons}{#pi^{0}}", 5.1,30., 0.,0.5);
    style(gr[iph], 20, 1);
    gr[iph]->Draw("AP");
  }
  c0->Print("plots/GammaRatio-minbias.pdf");
}
