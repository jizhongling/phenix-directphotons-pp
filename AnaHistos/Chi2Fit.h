double Chi2Fit(int n, double *x, double *ex, double &xbar, double &exbar)
{
  int sumN = 0;
  double sumw = 0.;
  double sumx = 0.;
  double sumx2 = 0.;

  for(int i=0; i<n; i++)
    if( TMath::Finite(x[i]+ex[i]) && ex[i] > 0. )
    {
      sumN++;
      sumw += 1./ex[i]/ex[i];
      sumx += x[i]/ex[i]/ex[i];
      sumx2 += x[i]*x[i]/ex[i]/ex[i];
    }

  xbar = sumx/sumw;
  exbar = 1./sqrt(sumw);
  int ndf = sumN - 1;
  int chi2 = sumx2 - sumw*xbar*xbar;

  return (ndf>0 ? chi2/ndf : chi2);
}

double Pol0Fit(int n, double *x, double *ex, double &xbar, double &exbar)
{
  double *dummy = new double[n];
  for(int i=0; i<n; i++)
    dummy[i] = i;

  TGraphErrors *gr_fit = new TGraphErrors(n, dummy, x, 0, ex);
  TFitResultPtr r_fit = gr_fit->Fit("pol0", "QS");
  xbar = r_fit->Value(0);
  exbar = r_fit->ParError(0);
  unsigned ndf = r_fit->Ndf();
  double chi2 = r_fit->Chi2();

  delete[] dummy;
}
