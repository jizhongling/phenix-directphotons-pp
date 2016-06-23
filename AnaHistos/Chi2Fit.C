void Chi2Fit(Int_t n, Double_t *x, Double_t *ex, Double_t &xbar, Double_t &exbar)
{
  Double_t A = 0.;
  Double_t B = 0.;
  Double_t C = 0.;
  for(Int_t i=0; i<n; i++)
    if(ex[i] > 0.)
    {
      A += 1. / ( ex[i] * ex[i] );
      B += x[i] / ( ex[i] * ex[i] );
      C += ( x[i] * x[i] ) / ( ex[i] * ex[i] );
    }

  xbar = B / A;
  Double_t Chi2 = 0.;
  for(Int_t i=0; i<n; i++)
    Chi2 += pow( ( x[i] - xbar ) / ex[i] , 2. );

  Double_t D = ( C - Chi2 - 1. ) / A;
  exbar = sqrt( xbar * xbar - D );

  return;
}
