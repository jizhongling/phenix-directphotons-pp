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

  TGraphErrors *gr = new TGraphErrors(npT);
  int igp = 0;

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-GenPH-histo-minbias.root");
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

    double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
    double ratio = nothers / npion;
    double eratio = ratio * sqrt( pow(enothers/nothers,2) + pow(enpion/npion,2) );
    if( TMath::Finite(ratio+eratio) )
    {
      gr->SetPoint(igp, xpt, ratio);
      gr->SetPointError(igp, 0., eratio);
      igp++;
    }
  }
  gr->Set(igp);

  mc();
  mcd();

  gr->SetTitle("Photon Ratio");
  aset(gr, "p_{T} [GeV]","#frac{#eta+#omega+#eta'}{#pi^{0}}", 5.1,30., 0.,0.5);
  style(gr, 20, 1);
  gr->Draw("AP");
  gr->Fit("pol0", "Q","", 6.1,30.);

  c0->Print("plots/GammaRatio-minbias.pdf");
}
