#include "BgGPR.h"

void BgGPRMinv(TH1 *h_minv, double &npeak, double &enpeak,
    double &nbg, double &enbg,
    const char *outfile = "default.root", const int index = 0)
{
  const int binL1 = h_minv->GetXaxis()->FindBin(0.067);
  const int binL2 = h_minv->GetXaxis()->FindBin(0.087);
  const int binC1 = h_minv->GetXaxis()->FindBin(0.112);
  const int binC2 = h_minv->GetXaxis()->FindBin(0.162);
  const int binR1 = h_minv->GetXaxis()->FindBin(0.187);
  const int binR2 = h_minv->GetXaxis()->FindBin(0.212);
  const double binW = h_minv->GetXaxis()->GetBinWidth(binC1);

  npeak = h_minv->IntegralAndError(binC1, binC2-1, enpeak);

  const unsigned nData = 45;
  vector<double> x(nData), y(nData), sigma_y(nData);
  x.clear();
  y.clear();
  sigma_y.clear();

  for(int biny=binL1; biny<=binR2; biny++)
    if(biny<=binL2 || biny>=binR1)
    {
      double xx = h_minv->GetXaxis()->GetBinCenter(biny);
      double yy = h_minv->GetBinContent(biny);
      double sigma_yy = h_minv->GetBinError(biny);
      x.push_back(xx);
      y.push_back(yy);
      sigma_y.push_back(sigma_yy);
    }

  BgGPR(x,y,sigma_y, nbg,enbg, 0.047,0.227,180, 0.112,0.162, outfile,index);
  nbg /= binW;
  enbg /= binW;

  return;
}

void BgGPRMinv(TH1 *h_minv, double &npion, double &enpion,
    const char *outfile = "default.root", const int index = 0)
{
  double npeak, enpeak, nbg, enbg;
  BgGPRMinv(h_minv, npeak, enpeak, nbg, enbg, outfile,index);

  npion = npeak - nbg;
  enpion = sqrt(enpeak*enpeak + enbg*enbg);

  return;
}
