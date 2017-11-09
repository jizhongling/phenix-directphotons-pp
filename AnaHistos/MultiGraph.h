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

template <class GraphType>
void GetMeanError(TMultiGraph *mg, Double_t &mean, Double_t &emean)
{
  Int_t sumN = 0;
  Double_t sumy = 0.;
  Double_t sumy2 = 0.;

  GraphType *gr;
  TIter iter( mg->GetListOfGraphs() );
  while( gr = (GraphType*)iter.Next() )
  {
    Int_t N = gr->GetN();
    if( N <= 0 ) continue;
    sumN += N;
    for(Int_t i=0; i<N; i++)
    {
      Double_t xx, yy;
      gr->GetPoint(i, xx, yy);
      sumy += yy;
      sumy2 += yy*yy;
    }
  }

  mean = sumy/sumN;
  Double_t rms2 = TMath::Abs( sumy2/sumN - mean*mean );
  emean = TMath::Sqrt( rms2/sumN );
  return;
}
