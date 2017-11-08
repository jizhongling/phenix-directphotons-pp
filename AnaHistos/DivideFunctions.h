TGraphErrors *DivideHisto(TH1 *h1, TH1 *h2)
{
  Int_t gn = h1->GetXaxis()->GetNbins();
  Double_t *gx = new Double_t[gn];
  Double_t *gy = new Double_t[gn];
  Double_t *egx = new Double_t[gn];
  Double_t *egy = new Double_t[gn];

  for(Int_t i=0; i<gn; i++)
  {
    gx[i] = gy[i] = egx[i] = egy[i] = 0.;

    Double_t h1x = h1->GetXaxis()->GetBinCenter(i+1);
    Double_t h1y = h1->GetBinContent(i+1); 
    Double_t h2y = h2->GetBinContent(i+1); 
    Double_t eh1x = h1->GetXaxis()->GetBinWidth(i+1) / 2.;
    Double_t eh1y = h1->GetBinError(i+1);
    Double_t eh2y = h2->GetBinError(i+1);

    if( h1y > 0. && h2y > 0. )
    {
      gx[i] = h1x;
      gy[i] = h1y / h2y;
      egx[i] = eh1x;
      egy[i] = gy[i] * sqrt( pow(eh1y/h1y,2.) + pow(eh2y/h2y,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(gn, gx, gy, egx, egy);
  return graph;
}

TGraphErrors *DivideHisto(TH1 *h1, TH1 *h2, Double_t kh1, Double_t kh2)
{
  Int_t gn = h1->GetXaxis()->GetNbins();
  Double_t *gx = new Double_t[gn];
  Double_t *gy = new Double_t[gn];
  Double_t *egx = new Double_t[gn];
  Double_t *egy = new Double_t[gn];

  for(Int_t i=0; i<gn; i++)
  {
    gx[i] = gy[i] = egx[i] = egy[i] = 0.;

    Double_t h1x = h1->GetXaxis()->GetBinCenter(i+1);
    Double_t h1y = h1->GetBinContent(i+1) * kh1; 
    Double_t h2y = h2->GetBinContent(i+1) * kh2; 
    Double_t eh1x = h1->GetXaxis()->GetBinWidth(i+1) / 2.;
    Double_t eh1y = h1->GetBinError(i+1) * kh1;
    Double_t eh2y = h2->GetBinError(i+1) * kh2;

    if( h1y > 0. && h2y > 0. )
    {
      gx[i] = h1x;
      gy[i] = h1y / h2y;
      egx[i] = eh1x;
      egy[i] = gy[i] * sqrt( pow(eh1y/h1y,2.) + pow(eh2y/h2y,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(gn, gx, gy, egx, egy);
  return graph;
}

TGraphErrors *DivideGraph(TGraphErrors *gr1, TGraphErrors *gr2)
{
  Int_t gn = gr1->GetN();
  Double_t *gx = new Double_t[gn];
  Double_t *gy = new Double_t[gn];
  Double_t *egx = new Double_t[gn];
  Double_t *egy = new Double_t[gn];

  for(Int_t i=0; i<gn; i++)
  {
    gx[i] = gy[i] = egx[i] = egy[i] = 0.;

    Double_t g1x, g1y, g2x, g2y;
    gr1->GetPoint(i, g1x, g1y);
    gr2->GetPoint(i, g2x, g2y);
    Double_t eg1x = gr1->GetErrorX(i);
    Double_t eg1y = gr1->GetErrorY(i);
    Double_t eg2y = gr2->GetErrorY(i);

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

void InvertGraph(TGraphErrors *gr)
{
  Int_t gn = gr->GetN();
  for(Int_t i=0; i<gn; i++)
  {
    Double_t gx, gy;
    gr->GetPoint(i, gx, gy);
    gr->SetPoint(i, gx, 1.-gy);
  }

  return;
}

void InvertGraph(TGraphAsymmErrors *gr)
{
  Int_t gn = gr->GetN();
  for(Int_t i=0; i<gn; i++)
  {
    Double_t gx, gy;
    gr->GetPoint(i, gx, gy);
    gr->SetPoint(i, gx, 1.-gy);
  }

  return;
}
