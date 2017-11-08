void draw_RawYieldCmp()
{
  TFile *f = new TFile("RawYieldCmp.root");
  TGraph *gr[3];
  for(Int_t part=0; part<3; part++)
    gr[part] = (TGraph*)f->Get(Form("gr_%d",part));

  const char *name[3] = {"PbScW", "PbScE", "PbGlE"};
  mc();
  mcd();
  legi();
  for(Int_t part=0; part<3; part++)
  {
    gr[part]->SetTitle("Raw yield ratio");
    aset(gr[part], "Runnumber","#frac{Mine}{Sasha}", 387000,398200, 0.4999,0.5001);
    style(gr[part], 20+part, 1+part);
    if(part==0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");
    leg0->AddEntry(gr[part], name[part], "P");
  }
  leg0->Draw();
  c0->Print("RawYieldCmp.pdf");
}
