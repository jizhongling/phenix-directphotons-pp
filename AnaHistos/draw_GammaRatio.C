#include "GlobalVars.h"
#include "QueryTree.h"
#include "DivideFunctions.h"

void draw_GammaRatio()
{
  const double BR = 1.98;

  QueryTree *qt_A = new QueryTree("data/GammaRatio.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-syserr.root");

  THnSparse *hn_hadron = (THnSparse*)f->Get("hn_hadron");
  hn_hadron->GetAxis(1)->SetRange(1,2);  // |eta| < 0.5
  TH1 *h_hadron[5];  // isophoton, pion, eta, omega, eta prime
  for(int id=0; id<5; id++)
  {
    hn_hadron->GetAxis(2)->SetRange(id+2,id+2);
    h_hadron[id] = hn_hadron->Projection(0);
    h_hadron[id]->SetName( Form("h_hadron%d",id) );
  }

  for(int iph=0; iph<2; iph++)
    for(int ipt=4; ipt<npT; ipt++)
    {
      double npion = h_hadron[1]->GetBinContent(ipt+1);
      double enpion = h_hadron[1]->GetBinError(ipt+1);

      double nothers = 0.;
      double enothers = 0.;
      if(iph == 0)
      {
        for(int id=2; id<5; id++)
        {
          nothers += h_hadron[id]->GetBinContent(ipt+1);
          enothers += h_hadron[id]->GetBinError(ipt+1);
        }
      }
      else if(iph == 1)
      {
        nothers += h_hadron[0]->GetBinContent(ipt+1) * BR;
        enothers += h_hadron[0]->GetBinError(ipt+1) * BR;
      }

      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      double ratio = nothers / npion;
      double eratio = ratio * sqrt( pow(enothers/nothers,2) + pow(enpion/npion,2) );
      if( TMath::Finite(ratio+eratio) )
        qt_A->Fill(ipt, iph, xpt, ratio, eratio);
    } // iph, ipt

  mc(0, 2,1);
  for(int iph=0; iph<2; iph++)
  {
    mcd(0, iph+1);
    TGraphErrors *gr = qt_A->Graph(iph);
    if(iph == 0)
    {
      gr->SetTitle("#gamma Ratio");
      aset(gr, "p_{T} [GeV]","#frac{#eta+#omega+#eta'}{#pi^{0}}", 5.1,30., 0.3,0.5);
    }
    else if(iph == 1)
    {
      gr->SetTitle("#gamma/#pi^{0} Ratio");
      aset(gr, "p_{T} [GeV]","#frac{#gamma}{#pi^{0}}", 5.1,30., 0.,1.);
    }
    style(gr, 20, 1);
    gr->Draw("AP");
  }
  c0->Print("plots/GammaRatio-minbias.pdf");

  qt_A->Save();
}
