Bool_t FitMinv(TH1 *h_minv, Double_t &npion, Double_t &enpion,
    const Bool_t bsub = kTRUE,
    const Double_t c1 = 0.112, const Double_t c2 = 0.162,
    const Double_t l1 = 0.047, const Double_t l2 = 0.097,
    const Double_t r1 = 0.177, const Double_t r2 = 0.227)
{
  const Int_t binL1 = h_minv->GetXaxis()->FindBin(l1);
  const Int_t binL2 = h_minv->GetXaxis()->FindBin(l2);
  const Int_t binC1 = h_minv->GetXaxis()->FindBin(c1);
  const Int_t binC2 = h_minv->GetXaxis()->FindBin(c2);
  const Int_t binR1 = h_minv->GetXaxis()->FindBin(r1);
  const Int_t binR2 = h_minv->GetXaxis()->FindBin(r2);

  Double_t max = h_minv->GetMaximum();
  if( max <= 0. )
    return kFALSE;

  aset(h_minv, "m_{inv} [GeV]","", 0.,0.3);

  TF1 *fn_fit = new TF1("fn_fit", "gaus(0) + pol2(3)", l1, r2);
  TF1 *fn_bg = new TF1("fn_bg", "pol2", l1, r2);

  Double_t par[10] = {max,0.137,0.010, 0.,0.,0.};
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

  Double_t nsig = 0.;
  Double_t nbg = 0.;
  for(Int_t ib=binC1; ib<binC2; ib++)
  {
    nsig += h_minv->GetBinContent(ib);
    Double_t bincenter = h_minv->GetXaxis()->GetBinCenter(ib);
    nbg += fn_bg->Eval(bincenter);
  }
  if( fn_fit->GetNDF() < binR2-binL1-10 ) 
    nbg = ( h_minv->Integral(binL1,binL2-1) + h_minv->Integral(binR1,binR2-1) ) / 2.;

  Double_t ensig = sqrt(nsig);
  Double_t rbg = nbg / nsig;
  Double_t erbg = sqrt(nbg) / nsig;

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
