double Chi2Fit(int n, double *x, double *ex, double &xbar, double &exbar)
{
  int sumN = 0;
  double sumw = 0.;
  double sumx = 0.;
  double sumx2 = 0.;

  for(int i=0; i<n; i++)
    if( TMath::Finite(ex[i]) && ex[i] > 0. )
    {
      sumN++;
      sumw += 1./ex[i]/ex[i];
      sumx += x[i]/ex[i]/ex[i];
      sumx2 += x[i]*x[i]/ex[i]/ex[i];
    }

  xbar = sumx/sumw;
  exbar = 1./sqrt(sumw);

  double chi2 = sumN > 1 ? ( sumx2 - sumw*xbar*xbar ) / ( sumN - 1 ) : 0.;
  return chi2;
}
