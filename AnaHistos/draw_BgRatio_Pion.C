#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "BgGPRMinv.h"

void draw_BgRatio_Pion()
{
  gSystem->Load("libGausProc.so");

  QueryTree *qt_rbg = new QueryTree("data/BgRatio-pion.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-Inseok.root");

  for(int icr=0; icr<2; icr++)
  {
    TH2 *h2_pion = (TH2*)f->Get(Form("h2_pion_pol_%d",icr));
    TH2 *h2_tmp = (TH2*)f->Get(Form("h2_pion_pol_%d",icr+2));
    h2_pion->Add(h2_tmp);

    for(int ipt=0; ipt<npT_pol; ipt++)
    {
      int ptbin_first = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt]);
      int ptbin_last = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt+1]) - 1;
      TH1 *h_minv = h2_pion->ProjectionY("h_minv", ptbin_first,ptbin_last);

      double npeak, enpeak, nbg, enbg;
      BgGPRMinv(h_minv, npeak, enpeak, nbg, enbg);

      double rbg = nbg/npeak;
      double erbg = rbg*sqrt(enpeak*enpeak/npeak/npeak + enbg*enbg/nbg/nbg);

      double xpt = (pTbin_pol[ipt] + pTbin_pol[ipt+1]) / 2.;
      qt_rbg->Fill(ipt, icr, xpt, rbg, erbg);
    } // ipt
  } // icr

  qt_rbg->Save();
}
