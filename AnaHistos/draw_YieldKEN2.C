#include "GlobalVars.h"
#include "QueryTree.h"

void draw_YieldKEN2()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  QueryTree *qt_ken2_pion = new QueryTree("data/YieldKEN2-pion.root", "RECREATE");
  for(int ibg=0; ibg<2; ibg++)
    for(int beam=0; beam<3; beam++)
      for(int icr=0; icr<2; icr++)
      {
        char *bgname = ibg ? "bg" : "sig";
        int ih = beam + 3*icr;
        TH2 *h2_mul = (TH2*)f->Get(Form("h2_mul_pion_%s_%d",bgname,ih));

        for(int ipt=0; ipt<npT_pol; ipt++)
        {
          TH1 *h_mul = h2_mul->ProjectionY("h_mul", ipt+1,ipt+1);

          double stats[4];
          h_mul->GetStats(stats);
          double ken2 = stats[3]/stats[2];

          double xpt = (pTbin_pol[ipt] + pTbin_pol[ipt+1]) / 2.;
          int part = ibg + 2*beam + 2*3*icr;
          qt_ken2_pion->Fill(ipt, part, xpt, ken2, 0.);

          delete h_mul;
        }
      }
  qt_ken2_pion->Save();

  QueryTree *qt_ken2 = new QueryTree("data/YieldKEN2-isophoton.root", "RECREATE");
  for(int imul=0; imul<6; imul++)
    for(int beam=0; beam<3; beam++)
      for(int icr=0; icr<2; icr++)
        for(int checkmap=0; checkmap<2; checkmap++)
          for(int ical=0; ical<2; ical++)
          {
            int ih = imul + 6*beam + 6*3*icr + 6*3*2*checkmap + 6*3*2*2*ical;
            TH2 *h2_mul = (TH2*)f->Get(Form("h2_mul_photon_%d",ih));
            if(imul == 0)
              h2_mul->Add((TH2*)f->Get(Form("h2_mul_photon_%d",ih+1)));

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
