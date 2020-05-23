#include "QueryTree.h"

void draw_ERTEff_SM()
{
  const int npT = 21;
  const double pTbin[npT+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    30.0 };

  QueryTree *qt_ertsm = new QueryTree("data/ERTEff-SM.root", "RECREATE");
  qt_ertsm->SetQuiet();

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-DC3sigma.root");

  for(int sector=0; sector<8; sector++)
    for(int evtype=0; evtype<3; evtype++)
    {
      int ert_trig = 1;
      int ih = sector + 8*ert_trig + 8*2*evtype;
      TH2 *h2_passed = (TH2*)f->Get(Form("h2_ertsm_%d",ih));

      ert_trig = 0;
      ih = sector + 8*ert_trig + 8*2*evtype;
      TH2 *h2_total = (TH2*)f->Get(Form("h2_ertsm_%d",ih));
      h2_total->Add(h2_passed);

      const int nsm = sector < 6 ? 18 : 32;
      for(int sm=0; sm<nsm; sm++)
      {
        TH1 *h_passed_0 = h2_passed->ProjectionX("h_passed_0", sm+1,sm+1);
        TH1 *h_total_0 = h2_total->ProjectionX("h_total_0", sm+1,sm+1);
        TH1 *h_passed = h_passed_0->Rebin(npT, "h_passed", pTbin);
        TH1 *h_total = h_total_0->Rebin(npT, "h_total", pTbin);

        int part = sector + 8*sm + 8*32*evtype;
        qt_ertsm->Fill(h_passed, h_total, part);

        delete h_passed_0;
        delete h_total_0;
        delete h_passed;
        delete h_total;
      }

      delete h2_passed;
      delete h2_total;
    }

  for(int evtype=0; evtype<3; evtype++)
    for(int ipt=0; ipt<npT; ipt++)
      for(int sector=0; sector<8; sector++)
      {
        TH2 *h2_eff = new TH2F(Form("h2_eff_ev%d_ipt%d_sec%d",evtype,ipt,sector),
            "EMCal SM efficiency;zpos;ypos;", 8,-0.5,7.5, 4,-0.5,3.5);
        h2_eff->Sumw2();

        const int nsm = sector < 6 ? 18 : 32;
        for(int sm=0; sm<nsm; sm++)
        {
          double xpt, value, error;
          int part = sector + 8*sm + 8*32*evtype;
          qt_ertsm->Query(ipt, part, xpt, value, error);

          int ypos = sector < 6 ? sm/6 : sm/8;
          int zpos = sector < 6 ? sm%6 : sm%8;
          h2_eff->Fill((double)zpos, (double)ypos, value);
        }

        qt_ertsm->cd();
        h2_eff->Write();
        delete h2_eff;
      }

  qt_ertsm->Save();
}
