#include "GlobalVars.h"
#include "QueryTree.h"

void draw_YieldKEN2()
{
  const char *tname[2] = {"sig", "bg"};

  QueryTree *qt_ken2 = new QueryTree("data/YeildKEN2-pion.pdf", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  for(int ibg=0; ibg<2; ibg++)
    for(int icr=0; icr<2; icr++)
    {
      TH2 *h2_mul = (TH2*)f->Get( Form("h2_mul_pion_%s_%d", tname[ibg], icr) );
      for(int ipt=0; ipt<npT_pol; ipt++)
      {
        TH1 *h_mul = h2_mul->ProjectionY("h_mul", ipt+1,ipt+1);

        double stats[4];
        h_mul->GetStats(stats);
        double ken2 = stats[3]/stats[2];

        int part = ibg + 2*icr;
        double xpt = (pTbin_pol[ipt] + pTbin_pol[ipt+1]) / 2.;
        qt_ken2->Fill(ipt, part, xpt, ken2, 0.);

        delete h_mul;
      }
    }

  qt_ken2->Save();
}
