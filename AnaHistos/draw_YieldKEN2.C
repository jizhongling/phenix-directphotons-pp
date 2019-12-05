#include "GlobalVars.h"
#include "QueryTree.h"

void draw_YieldKEN2()
{
  QueryTree *qt_ken2 = new QueryTree("data/YieldKEN2-isophoton.pdf", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  for(int imul=0; imul<5; imul++)
    for(int beam=0; beam<3; beam++)
      for(int icr=0; icr<2; icr++)
        for(int checkmap=0; checkmap<2; checkmap++)
        {
          int ih = imul + 5*beam + 5*3*icr + 5*3*2*checkmap;
          TH2 *h2_mul = (TH2*)f->Get( Form("h2_mul_photon_%d", ih) );

          for(int ipt=0; ipt<npT_pol; ipt++)
          {
            TH1 *h_mul = h2_mul->ProjectionY("h_mul", ipt+1,ipt+1);

            double stats[4];
            h_mul->GetStats(stats);
            double ken2 = stats[3]/stats[2];

            double xpt = (pTbin_pol[ipt] + pTbin_pol[ipt+1]) / 2.;
            qt_ken2->Fill(ipt, ih, xpt, ken2, 0.);

            delete h_mul;
          }
        }

  qt_ken2->Save();
}
