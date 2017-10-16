void draw_RawYieldCmp()
{
  TFile *f = new TFile("RawYieldCmp.root");
  TGraph *gr[3];
  for(Int_t part=0; part<3; part++)
    gr[part] = (TGraph*)f->Get(Form("gr_%d",part));

  mc();
  mcd();
  for(Int_t part=0; part<3; part++)
  {
    gr[part]->SetTitle("Raw yield ratio");
    aset(gr[part], "Runnumber","#frac{Mine}{Sasha}", 387000,398200, 0.,2.);
    style(gr[part], 20+part, 1+part);
    if(part==0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");
  }
  c0->Print("RawYieldCmp.pdf");
}
