#include "GlobalVars.h"
#include "QueryTree.h"

void draw_PtShift()
{
  QueryTree *qt_pt = new QueryTree("data/PtShift.root", "RECREATE");

  QueryTree *qt_rbg = new QueryTree("data/BgRatio-isophoton.root");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-DC3sigma.root");

  int beam = 2;
  int checkmap = 1;
  int ical = 0;

  TH1 *h_1photon[2];  // isolated
  TH1 *h_2photon[2][3];  // pttype, isotype[inclusive|isoboth|isopair]

  TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_pol_0");
  h_1photon_t = (TH1*)h_1photon_t->Clone();
  h_1photon_t->Reset();

  for(int isolated=0; isolated<2; isolated++)
  {
    h_1photon[isolated] = (TH1*)h_1photon_t->Clone(Form("h_1photon_isolated%d",isolated));
    for(int icr=0; icr<2; icr++)
      for(int ipol=0; ipol<2; ipol++)
      {
        int ih = beam + 3*icr + 3*2*ipol + 3*2*2*checkmap + 3*2*2*2*isolated + 3*2*2*2*2*ical;
        TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_pol_%d",ih));
        h_1photon[isolated]->Add(h_tmp);
      } // icr, ipol
  } // isolated
  h_1photon[0]->Add(h_1photon[1]);

  for(int pttype=0; pttype<2; pttype++)
    for(int isotype=0; isotype<3; isotype++)
    {
      const char *ptname = pttype ? "2pt" : "";
      h_2photon[pttype][isotype] = (TH1*)h_1photon_t->Clone(Form("h_2photon%s_isotype%d",ptname,isotype));
      for(int isoboth=(isotype==1?1:0); isoboth<2; isoboth++)
        for(int isopair=(isotype==2?1:0); isopair<2; isopair++)
          for(int icr=0; icr<2; icr++)
            for(int ipol=0; ipol<2; ipol++)
            {
              int ih = beam + 3*icr + 3*2*ipol + 3*2*2*checkmap + 3*2*2*2*isoboth + 3*2*2*2*2*isopair + 3*2*2*2*2*2*ical;
              TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon%s_pol_%d",ptname,ih));
              TH1 *h_tmp = h2_tmp->ProjectionX("h_tmp", 111,160);
              h_2photon[pttype][isotype]->Add(h_tmp);
            } // isoboth, isopair, icr, ipol
    } // pttype, isotype

  for(int isolated=0; isolated<2; isolated++)
  {
    for(int ipt=0; ipt<npT; ipt++)
    {
      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      int ipt_pol = Get_ipt_pol(xpt);
      double dummy, tmp[2], rbg[3];
      for(int ibg=0; ibg<3; ibg++)
      {
        for(int icr=0; icr<2; icr++)
          qt_rbg->Query(ipt_pol, 3*icr+ibg, dummy, tmp[icr], dummy);
        rbg[ibg] = (tmp[0] + tmp[1]) / 2.;
      }

      double stats[4], ptbar[4];
      h_1photon[isolated]->GetXaxis()->SetRange(pTbin[ipt]*10+1, pTbin[ipt+1]*10);
      h_1photon[isolated]->GetStats(stats);
      ptbar[3] = stats[2]/stats[0];
      h_2photon[0][isolated]->GetXaxis()->SetRange(pTbin[ipt]*10+1, pTbin[ipt+1]*10);
      h_2photon[0][isolated]->GetStats(stats);
      ptbar[0] = stats[2]/stats[0];
      h_2photon[1][isolated*2]->GetXaxis()->SetRange(pTbin[ipt]*10+1, pTbin[ipt+1]*10);
      h_2photon[1][isolated*2]->GetStats(stats);
      ptbar[1] = stats[2]/stats[0];
      h_2photon[0][isolated*2]->GetXaxis()->SetRange(pTbin[ipt]*10+1, pTbin[ipt+1]*10);
      h_2photon[0][isolated*2]->GetStats(stats);
      ptbar[2] = stats[2]/stats[0];
      double ptshift = (ptbar[3] - rbg[0]*ptbar[0] - rbg[1]*ptbar[1] - rbg[2]*ptbar[2]) / (1 - rbg[0] - rbg[1] - rbg[2]);
      if(!TMath::Finite(ptshift))
        ptshift = 0.25;
      qt_pt->Fill(ipt, isolated, xpt, ptshift, 0.1);
    } // ipt
  } // isolated

  int isolated = 1;
  for(int ipt=0; ipt<npT_pol; ipt++)
  {
    double xpt = (pTbin_pol[ipt] + pTbin_pol[ipt+1]) / 2.;
    double dummy, tmp[2], rbg[3];
    for(int ibg=0; ibg<3; ibg++)
    {
      for(int icr=0; icr<2; icr++)
        qt_rbg->Query(ipt, 3*icr+ibg, dummy, tmp[icr], dummy);
      rbg[ibg] = (tmp[0] + tmp[1]) / 2.;
    }

    double stats[4], ptbar[4];
    h_1photon[isolated]->GetXaxis()->SetRange(pTbin_pol[ipt]*10+1, pTbin_pol[ipt+1]*10);
    h_1photon[isolated]->GetStats(stats);
    ptbar[3] = stats[2]/stats[0];
    h_2photon[0][isolated]->GetXaxis()->SetRange(pTbin_pol[ipt]*10+1, pTbin_pol[ipt+1]*10);
    h_2photon[0][isolated]->GetStats(stats);
    ptbar[0] = stats[2]/stats[0];
    h_2photon[1][isolated*2]->GetXaxis()->SetRange(pTbin_pol[ipt]*10+1, pTbin_pol[ipt+1]*10);
    h_2photon[1][isolated*2]->GetStats(stats);
    ptbar[1] = stats[2]/stats[0];
    h_2photon[0][isolated*2]->GetXaxis()->SetRange(pTbin_pol[ipt]*10+1, pTbin_pol[ipt+1]*10);
    h_2photon[0][isolated*2]->GetStats(stats);
    ptbar[2] = stats[2]/stats[0];
    double ptshift = (ptbar[3] - rbg[0]*ptbar[0] - rbg[1]*ptbar[1] - rbg[2]*ptbar[2]) / (1 - rbg[0] - rbg[1] - rbg[2]);
    if(!TMath::Finite(ptshift))
      ptshift = 0.25;
    qt_pt->Fill(ipt, 2, xpt, ptshift, 0.1);
  } // ipt

  qt_pt->Save();
}
