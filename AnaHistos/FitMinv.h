const double bck[2][30] = {
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.10, 1.15, 1.20, 1.30, 1, 1, 1 },
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.08, 1.08, 1.11, 1.11, 1.11, 1.11, 1.11 }
};

const double meff[2][30] = {
  { 1, 1, 0.96, 0.97, 0.98, 0.985, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.985, 0.995, 0.995, 0.99, 0.98, 0.95, 1, 1, 1, 1, 1 },
  { 1, 1, 0.95, 0.97, 0.975, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.995, 0.995, 0.99, 0.99, 0.98, 0.98, 0.98, 0.98, 1, 1 }
};

bool FitMinv(TH1 *h_minv, double &npion, double &enpion,
    const bool bsub = true,
    const double c1 = 0.11, const double c2 = 0.16,
    const double l1 = 0.06, const double r2 = 0.25)
{
  const int binC1 = h_minv->GetXaxis()->FindBin(c1);
  const int binC2 = h_minv->GetXaxis()->FindBin(c2);

  const double max = h_minv->GetMaximum();
  if( max <= 0. )
    return false;

  double nsig = 0.;
  double e2nsig = 0.;
  for(int ib=binC1; ib<binC2; ib++)
  {
    nsig += h_minv->GetBinContent(ib);
    double ensig = h_minv->GetBinError(ib);
    e2nsig += ensig*ensig;
  }

  if(!bsub)
  {
    npion = nsig;
    enpion = sqrt(e2nsig);
    return true;
  }

  aset(h_minv, "m_{inv} [GeV]","", 0.,0.3, 0.,1.1*max);

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

  double nbg = 0.;
  double e2nbg = 0.;
  for(int ib=binC1; ib<binC2; ib++)
  {
    double bincenter = h_minv->GetXaxis()->GetBinCenter(ib);
    nbg += fn_bg->Eval(bincenter);
    if(nsig > 0.)
      e2nbg += nbg*e2nsig/nsig;
  }

  double rbg = nbg/nsig;
  npion = nsig*(1-rbg);
  enpion = sqrt(e2nsig*(1-rbg)*(1-rbg) + e2nbg);

  delete fn_fit;
  delete fn_bg;
  return true;
}
