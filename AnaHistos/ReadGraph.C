void ReadGraph(TGraph *gr, Double_t *gx, Double_t *gy, Double_t *egy)
{
  const Int_t gn = 30;
  const Double_t vpt[gn+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

  for(Int_t ipt=0; ipt<gn; ipt++)
  {
    gx[ipt] = ( vpt[ipt] + vpt[ipt+1] ) / 2.;

    Double_t xx, yy;
    Int_t status = gr->GetPoint(ipt, xx, yy);
    if(status>0)
      gy[ipt] = yy;
    else
      gy[ipt] = gy[ipt-1];

    Double_t eyy = gr->GetErrorY(ipt);
    if(eyy>0)
      egy[ipt] = eyy;
    else
      egy[ipt] = egy[ipt-1];
  }

  return;
}

void ReadGraphErrors(const char *name, Int_t igr, Double_t *gx, Double_t *gy, Double_t *egy)
{
  TFile *f = new TFile(name);
  TGraphErrors *gr = (TGraphErrors*)f->Get(Form("gr_%d",igr));
  ReadGraph(gr, gx, gy, egy);
  f->Close();

  return;
}

void ReadGraphAsymmErrors(const char *name, Int_t igr, Double_t *gx, Double_t *gy, Double_t *egy)
{
  TFile *f = new TFile(name);
  TGraphAsymmErrors *gr = (TGraphAsymmErrors*)f->Get(Form("gr_%d",igr));
  ReadGraph(gr, gx, gy, egy);
  f->Close();

  return;
}
