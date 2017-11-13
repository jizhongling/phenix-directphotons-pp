template <class GraphType>
Double_t GetMaximum(TMultiGraph *mg)
{
  Double_t max = -1e9;

  GraphType *gr;
  TIter iter( mg->GetListOfGraphs() );
  while( gr = (GraphType*)iter.Next() )
  {
    Int_t N = gr->GetN();
    if( N <= 0 ) continue;
    for(Int_t i=0; i<N; i++)
    {
      Double_t xx, yy;
      gr->GetPoint(i, xx, yy);
      if( yy > max )
        max = yy;
    }
  }

  return max;
}

template <class GraphType>
Double_t GetMinimum(TMultiGraph *mg)
{
  Double_t min = 1e9;

  GraphType *gr;
  TIter iter( mg->GetListOfGraphs() );
  while( gr = (GraphType*)iter.Next() )
  {
    Int_t N = gr->GetN();
    if( N <= 0 ) continue;
    for(Int_t i=0; i<N; i++)
    {
      Double_t xx, yy;
      gr->GetPoint(i, xx, yy);
      if( yy < min )
        min = yy;
    }
  }

  return min;
}

Double_t GetMeanError(TMultiGraph *mg, Double_t &mean, Double_t &emean)
{
  Int_t sumN = 0;
  Double_t sumw = 0.;
  Double_t sumy = 0.;
  Double_t sumy2 = 0.;

  TGraphErrors *gr;
  TIter iter( mg->GetListOfGraphs() );
  while( gr = (TGraphErrors*)iter.Next() )
  {
    Int_t N = gr->GetN();
    if( N <= 0 ) continue;
    for(Int_t i=0; i<N; i++)
    {
      Double_t xx, yy;
      gr->GetPoint(i, xx, yy);
      Double_t eyy = gr->GetErrorY(i);
      if( eyy > 0. && eyy < TMath::Infinity() )
      {
        sumN++;
        sumw += 1./eyy/eyy;
        sumy += yy/eyy/eyy;
        sumy2 += yy*yy/eyy/eyy;
      }
    }
  }

  mean = sumy/sumw;
  emean = 1./sqrt(sumw);

  Double_t chi2 = sumN > 0 ? ( sumy2 - sumw*mean*mean ) / ( sumN - 1 ) : 1e9;
  return chi2;
}
