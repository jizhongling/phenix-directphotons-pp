Double_t Chi2Fit(Int_t n, Double_t *x, Double_t *ex, Double_t &xbar, Double_t &exbar)
{
  Int_t sumN = 0;
  Double_t sumw = 0.;
  Double_t sumx = 0.;
  Double_t sumx2 = 0.;

  for(Int_t i=0; i<n; i++)
    if( ex[i] > 0. && ex[i] < TMath::Infinity() )
    {
      sumN++;
      sumw += 1./ex[i]/ex[i];
      sumx += x[i]/ex[i]/ex[i];
      sumx2 += x[i]*x[i]/ex[i]/ex[i];
    }

  xbar = sumx/sumw;
  exbar = 1./sqrt(sumw);

  Double_t chi2 = sumN > 0 ? ( sumx2 - sumw*xbar*xbar ) / ( sumN - 1 ) : 1e9;
  return chi2;
}
