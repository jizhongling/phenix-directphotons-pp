void Chi2Fit(Int_t n, Double_t *x, Double_t *ex, Double_t &xbar, Double_t &exbar)
{
  Double_t A = 0.;
  Double_t B = 0.;

  for(Int_t i=0; i<n; i++)
    if(ex[i] > 0.)
    {
      A += 1. / ( ex[i] * ex[i] );
      B += x[i] / ( ex[i] * ex[i] );
    }

  xbar = B / A;
  exbar = 1. / sqrt(A);

  return;
}
