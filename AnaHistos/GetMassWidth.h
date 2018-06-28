bool GetMassWidth(TH1 *h_minv, double &mass, double &emass, double &width, double &ewidth)
{
  double max = h_minv->GetMaximum();
  if( max <= 0. )
    return false;

  aset(h_minv, "m_{inv} [GeV]","", 0.,0.3);

  TF1 *fn_fit = new TF1("fn_fit", "gaus(0) + pol2(3)", 0.06, 0.25);

  double par[10] = {max,0.135,0.010, 0.,0.,0.};
  fn_fit->SetParameters(par);
  h_minv->Fit(fn_fit, "RQ0");

  fn_fit->SetLineColor(kRed);
  h_minv->DrawCopy("HIST");
  fn_fit->DrawCopy("SAME");

  double scale = 1.;
  if( fn_fit->GetNDF() > 0 )
    scale = sqrt( fn_fit->GetChisquare() / fn_fit->GetNDF() );

  mass = fn_fit->GetParameter(1);
  emass = fn_fit->GetParError(1) * scale;
  width = fn_fit->GetParameter(2);
  ewidth = fn_fit->GetParError(2) * scale;

  delete fn_fit;
  return true;
}
