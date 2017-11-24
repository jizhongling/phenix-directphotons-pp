Bool_t FitMinv(TH1 *h_minv, Double_t &npion, Double_t &enpion)
{
  Double_t max = h_minv->GetMaximum();
  if( max <= 0. )
    return kFALSE;

  aset(h_minv, "m_{inv} [GeV]","", 0.,0.3);

  TF1 *fn_fit = new TF1("fn_fit", "gaus(0) + pol2(3)", 0.06, 0.25);
  TF1 *fn_bg = new TF1("fn_bg", "pol2", 0.06, 0.25);

  Double_t par[10] = {max,0.135,0.010, 0.,0.,0.};
  fn_fit->SetParameters(par);
  h_minv->Fit(fn_fit, "RQ0");
  fn_fit->GetParameters(par);
  fn_bg->SetParameters(par+3);

  fn_fit->SetLineColor(kRed);
  fn_bg->SetLineColor(kGreen);
  h_minv->DrawCopy();
  fn_fit->DrawCopy("SAME");
  fn_bg->DrawCopy("SAME");

  Double_t nsig = 0.;
  Double_t nbg = 0.;
  for(Int_t ib=12; ib<=16; ib++)
  {
    nsig += h_minv->GetBinContent(ib);
    Double_t bincenter = h_minv->GetXaxis()->GetBinCenter(ib);
    nbg += fn_bg->Eval(bincenter);
  }
  if( fn_fit->GetNDF() < 10 ) 
    nbg = ( h_minv->Integral(6,10) + h_minv->Integral(19,23) ) / 2.;

  Double_t ensig = sqrt(nsig);
  Double_t rbg = nbg / nsig;
  Double_t erbg = sqrt(nbg) / nsig;

  npion = nsig * (1-rbg);
  enpion = sqrt( pow(ensig*(1-rbg),2.) + pow(nsig*erbg,2.) );

  delete fn_fit;
  delete fn_bg;
  return kTRUE;
}
