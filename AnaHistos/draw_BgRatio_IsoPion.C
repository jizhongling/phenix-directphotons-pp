#include "GlobalVars.h"
#include "QueryTree.h"
#include "BgGPR.h"

void draw_BgRatio_IsoPion()
{
  gSystem->Load("libGausProc.so");

  QueryTree *qt_rbg = new QueryTree("data/BgRatio-isopion.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  TH2 *h2_pion = (TH2*)f->Get("h2_2photon_pol_0");
  h2_pion = (TH2*)h2_pion->Clone();

  int beam = 2;
  int checkmap = 1;

  const unsigned nData = 45;
  vector<double> x(nData), y(nData), sigma_y(nData);

  for(int pttype=0; pttype<2; pttype++)
    for(int icr=0; icr<2; icr++)
    {
      char *ptname = pttype ? "2pt" : "";
      h2_pion->Reset();
      for(int isoboth=0; isoboth<2; isoboth++)
        for(int isopair=0; isopair<2; isopair++)
          for(int ipol=0; ipol<2; ipol++)
          {
            int ih = beam + 3*icr + 3*2*ipol + 3*2*2*checkmap + 3*2*2*2*isoboth + 3*2*2*2*2*isopair;
            TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon%s_pol_%d",ptname,ih));
            h2_pion->Add(h2_tmp);
          }

      for(int ipt=0; ipt<npT_pol; ipt++)
      {
        int ptbin_first = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt]);
        int ptbin_last = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt+1]) - 1;
        TH1 *h_minv = h2_pion->ProjectionY("h_minv", ptbin_first,ptbin_last);

        x.clear();
        y.clear();
        sigma_y.clear();

        for(int biny=68; biny<=87; biny++)
        {
          double xx = h_minv->GetXaxis()->GetBinCenter(biny);
          double yy = h_minv->GetBinContent(biny);
          double sigma_yy = h_minv->GetBinError(biny);
          x.push_back(xx);
          y.push_back(yy);
          sigma_y.push_back(sigma_yy);
        }

        for(int biny=188; biny<=212; biny++)
        {
          double xx = h_minv->GetXaxis()->GetBinCenter(biny);
          double yy = h_minv->GetBinContent(biny);
          double sigma_yy = h_minv->GetBinError(biny);
          x.push_back(xx);
          y.push_back(yy);
          sigma_y.push_back(sigma_yy);
        }

        double nbg, dnbg;
        BgGPR(x, y, sigma_y, nbg, dnbg);
        nbg /= 0.001;
        dnbg /= 0.001;

        double npeak = h_minv->Integral(113, 162);
        double rbg = nbg/npeak;
        double erbg = rbg*sqrt(dnbg*dnbg/nbg/nbg + 1./npeak);
        if( !TMath::Finite(rbg+erbg) || rbg < 0. || erbg > 1. )
        {
          rbg = 0.;
          erbg = 1.;
        }

        double xpt = (pTbin_pol[ipt] + pTbin_pol[ipt+1]) / 2.;
        int part = pttype + 2*icr;
        qt_rbg->Fill(ipt, part, xpt, rbg, erbg);
      } // ipt
    } // pttype, icr

  qt_rbg->Save();
}
