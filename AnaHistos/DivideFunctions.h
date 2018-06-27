TGraphErrors *DivideHisto(TH1 *h1, TH1 *h2, double kh1 = 1., double kh2 = 1.)
{
  int gn = h1->GetXaxis()->GetNbins();

  vector<double> gx;
  vector<double> gy;
  vector<double> egx;
  vector<double> egy;

  for(int i=0; i<gn; i++)
  {
    double h1x = h1->GetXaxis()->GetBinCenter(i+1);
    double h1y = h1->GetBinContent(i+1) * kh1; 
    double h2y = h2->GetBinContent(i+1) * kh2; 
    double eh1x = h1->GetXaxis()->GetBinWidth(i+1) / 2.;
    double eh1y = h1->GetBinError(i+1) * kh1;
    double eh2y = h2->GetBinError(i+1) * kh2;

    if( h1y > 0. && h2y > 0. )
    {
      gx.push_back( h1x );
      gy.push_back( h1y / h2y );
      egx.push_back( eh1x );
      egy.push_back( gy[i] * sqrt( pow(eh1y/h1y,2.) + pow(eh2y/h2y,2.) ) );
    }
  }

  TGraphErrors *graph = new TGraphErrors( gx.size(), &(gx[0]), &(gy[0]), &(egx[0]), &(egy[0]) );
  return graph;
}

TGraphErrors *DivideGraph(TGraphErrors *gr1, TGraphErrors *gr2)
{
  int gn = gr1->GetN();
  double *gx = new double[gn];
  double *gy = new double[gn];
  double *egx = new double[gn];
  double *egy = new double[gn];

  for(int i=0; i<gn; i++)
  {
    gx[i] = gy[i] = egx[i] = egy[i] = 0.;

    double g1x, g1y, g2x, g2y;
    gr1->GetPoint(i, g1x, g1y);
    gr2->GetPoint(i, g2x, g2y);
    double eg1x = gr1->GetErrorX(i);
    double eg1y = gr1->GetErrorY(i);
    double eg2y = gr2->GetErrorY(i);

    if( g1y > 0. && g2y > 0. )
    {
      gx[i] = g1x;
      gy[i] = g1y / g2y;
      egx[i] = eg1x;
      egy[i] = gy[i] * sqrt( pow(eg1y/g1y,2.) + pow(eg2y/g2y,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(gn, gx, gy, egx, egy);
  return graph;
}

template <class GraphType>
void InvertGraph(GraphType *gr)
{
  int gn = gr->GetN();
  for(int i=0; i<gn; i++)
  {
    double gx, gy;
    gr->GetPoint(i, gx, gy);
    gr->SetPoint(i, gx, 1.-gy);
  }

  return;
}
