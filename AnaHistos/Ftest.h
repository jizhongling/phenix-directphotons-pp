#include "Chi2Fit.h"

void Ftest(const int ngroup, const vector<double> *data, const vector<double> *edata, double &F, double &p)
{
  int *n = new int[ngroup];
  double *bar = new double[ngroup];
  double *ebar = new double[ngroup];
  for(int i=0; i<ngroup; i++)
  {
    n[i] = data[i].size();
    if( n[i] != edata[i].size() )
    {
      cout << "Input has wrong dimension!" << endl;
      return;
    }
    Chi2Fit(n[i], (double*)&data[i][0], (double*)&edata[i][0], bar[i], ebar[i]);
  }

  double bar_all, ebar_all;
  Chi2Fit(ngroup, bar, ebar, bar_all, ebar_all);

  int d1 = ngroup - 1;
  int d2 = 0;
  double s1squre = 0.;
  double s2squre = 0.;
  for(int i=0; i<ngroup; i++)
  {
    d2 += n[i] - 1;
    s1squre += pow( (bar[i] - bar_all) / ebar[i], 2);
    for(int j=0; j<n[i]; j++)
      s2squre += pow( (data[i][j] - bar[i]) / edata[i][j], 2);
  }
  s1squre /= d1;
  s2squre /= d2;

  F = s1squre / s2squre;
  p = 1. - TMath::FDistI(F, d1, d2);
}
