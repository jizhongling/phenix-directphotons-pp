bool FitMinv(TH1 *h_minv, double &npion, double &enpion,
    const bool bsub = kTRUE,
    const double c1 = 0.11, const double c2 = 0.16,
    const double l1 = 0.06, const double r2 = 0.25)
{
  const int binC1 = h_minv->GetXaxis()->FindBin(c1);
  const int binC2 = h_minv->GetXaxis()->FindBin(c2);

  const double max = h_minv->GetMaximum();
  if( max <= 0. )
    return kFALSE;

  aset(h_minv, "m_{inv} [GeV]","", 0.,0.3);

  TF1 *fn_fit = new TF1("fn_fit", "gaus(0) + pol2(3)", l1, r2);
  TF1 *fn_bg = new TF1("fn_bg", "pol2", l1, r2);

  double par[10] = {max,0.137,0.010, 0.,0.,0.};
  fn_fit->SetParameters(par);
  h_minv->Fit(fn_fit, "RQ0");
  fn_fit->GetParameters(par);
  fn_bg->SetParameters(par+3);

  fn_fit->SetLineColor(kRed);
  fn_bg->SetLineColor(kGreen);
  fn_fit->SetLineWidth(1);
  fn_bg->SetLineWidth(1);
  h_minv->DrawCopy();
  fn_fit->DrawCopy("SAME");
  fn_bg->DrawCopy("SAME");

  double nsig = 0.;
  double nbg = 0.;
  for(int ib=binC1; ib<binC2; ib++)
  {
    nsig += h_minv->GetBinContent(ib);
    double bincenter = h_minv->GetXaxis()->GetBinCenter(ib);
    nbg += fn_bg->Eval(bincenter);
  }

  double ensig = sqrt(nsig);
  double rbg = nbg / nsig;
  double erbg = sqrt(nbg) / nsig;

  if(bsub)
  {
    npion = nsig * (1-rbg);
    enpion = sqrt( pow(ensig*(1-rbg),2.) + pow(nsig*erbg,2.) );
  }
  else
  {
    npion = nsig;
    enpion = ensig;
  }

  delete fn_fit;
  delete fn_bg;
  return kTRUE;
}
