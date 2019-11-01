#include "GlobalVars.h"
#include "QueryTree.h"

void BgGPR(vector<double> &x, vector<double> &y, vector<double> &sigma_y,
    double &Integral, double &dIntegral)
{
  const int verbosity = 0;
  const char *outfile = "BgGPR.root";

  const double xmin = 0.047;
  const double xmax = 0.227;
  const int nPredictions = 180;

  GausProc a(x, y, sigma_y, xmin, xmax, nPredictions, outfile);
  //gSystem->Exec(Form("rm -f %s",outfile));

  a.SetVerbosity(verbosity);
  a.SetKernel(GausProc::RBF);
  //a.warp(0);  //should be used only if you have data spanning several orders of magnitude

  GPOptimizer c(&a, 2., 10.);
  c.GPoptimize(2,0);

  GausProc d(x, y, sigma_y, xmin, xmax, nPredictions, outfile);
  //d.warp(0);  //see comment above
  d.SetPar(0, c.getPar(0));
  d.SetPar(1, c.getPar(1));
  d.process();
  //d.Write(-1);
  //d.unwarp(0);  //see comment above
  //d.Write(-1, "_unwarp");

  d.Integral(0.112, 0.162, Integral, dIntegral);

  return;
}

void draw_BgRatio_Pion()
{
  gSystem->Load("libGausProc.so");

  QueryTree *qt_rbg = new QueryTree("data/BgRatio-pion.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  const unsigned nData = 45;
  vector<double> x(nData), y(nData), sigma_y(nData);

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

      double xpt = (pTbin_pol[ipt] + pTbin_pol[ipt+1]) / 2.;
      qt_rbg->Fill(ipt, icr, xpt, rbg, erbg);
    } // ipt
  } // icr

  qt_rbg->Save();
}
