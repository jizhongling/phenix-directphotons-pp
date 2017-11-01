template <class GraphType>

void GetMeanRMS(TMultiGraph *mg, Double_t &mean, Double_t &rms)
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
  rms = TMath::Sqrt( rms2 );
  return;
}
