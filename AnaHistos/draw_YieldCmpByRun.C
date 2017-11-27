void draw_YieldCmpByRun()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};

  TFile *f = new TFile("data/YieldCmpByRun.root");
  TGraph *gr[3];
  for(Int_t part=0; part<3; part++)
    gr[part] = (TGraph*)f->Get(Form("gr_%d",part));

  mc();
  mcd();
  legi();

  for(Int_t part=0; part<3; part++)
  {
    gr[part]->SetTitle("Raw yield ratio");
    aset(gr[part], "Runnumber","#frac{Mine}{Sasha}", 387000,398200, 0.4999,0.5001);
    style(gr[part], part+24, part+1);
    if(part==0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");
    leg0->AddEntry(gr[part], pname[part], "P");
  }
  leg0->Draw();

  c0->Print("plots/YieldCmpByRun.pdf");
}
