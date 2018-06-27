template <class GraphType>
double GetMaximum(TMultiGraph *mg)
{
  double max = -1e9;

  GraphType *gr;
  TIter iter( mg->GetListOfGraphs() );
  while( gr = (GraphType*)iter.Next() )
  {
    int N = gr->GetN();
    if( N <= 0 ) continue;
    for(int i=0; i<N; i++)
    {
      double xx, yy;
      gr->GetPoint(i, xx, yy);
      if( yy > max )
        max = yy;
    }
  }

  return max;
}

template <class GraphType>
double GetMinimum(TMultiGraph *mg)
{
  double min = 1e9;

  GraphType *gr;
  TIter iter( mg->GetListOfGraphs() );
  while( gr = (GraphType*)iter.Next() )
  {
    int N = gr->GetN();
    if( N <= 0 ) continue;
    for(int i=0; i<N; i++)
    {
      double xx, yy;
      gr->GetPoint(i, xx, yy);
      if( yy < min )
        min = yy;
    }
  }

  return min;
}

double GetMeanError(TMultiGraph *mg, double &mean, double &emean)
{
  int sumN = 0;
  double sumw = 0.;
  double sumy = 0.;
  double sumy2 = 0.;

  TGraphErrors *gr;
  TIter iter( mg->GetListOfGraphs() );
  while( gr = (TGraphErrors*)iter.Next() )
  {
    int N = gr->GetN();
    if( N <= 0 ) continue;
    for(int i=0; i<N; i++)
    {
      double xx, yy;
      gr->GetPoint(i, xx, yy);
      double eyy = gr->GetErrorY(i);
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

  double chi2 = sumN > 1 ? ( sumy2 - sumw*mean*mean ) / ( sumN - 1 ) : 0.;
  return chi2;
}
