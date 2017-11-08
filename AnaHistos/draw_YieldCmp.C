void draw_YieldCmp()
{
  gROOT->ProcessLine(".L ReadGraph.C");

  const Double_t gx[30];
  Double_t gy[3][30];
  Double_t egy[3][30];
  Double_t gy2[3][30];
  Double_t egy2[3][30];
  for(Int_t part=0; part<3; part++)
  {
    ReadGraphErrors("Yield.root", part, gx, (Double_t*)gy[part], (Double_t*)egy[part]);
    ReadGraphErrors("YieldSasha.root", part, gx, (Double_t*)gy2[part], (Double_t*)egy2[part]);
  }

  TGraphErrors *gr_ratio[3];
  Double_t rgy[3][30] = {};
  Double_t ergy[3][30] = {};

  for(Int_t part=0; part<3; part++)
  {
    for(Int_t ipt=2; ipt<21; ipt++)
    {
      rgy[part][ipt] = ( gy[part][ipt] - gy2[part][ipt] ) / ( gy2[part][ipt] + 1e-9 );
      ergy[part][ipt] = gy[part][ipt] / ( gy2[part][ipt] + 1e-9 ) * sqrt( pow(egy[part][ipt]/gy[part][ipt],2.) + pow(egy2[part][ipt]/gy2[part][ipt],2.) );
    }
    gr_ratio[part] =  new TGraphErrors(30, gx, rgy[part], 0, ergy[part]);
  }

  mc();
  mcd();
  legi(0, 0.4,0.7,0.6,0.9);

  const char *pname[3] = {"PbScW", "PbScE", "PbGlE"};
  for(Int_t part=0; part<3; part++)
  {
    gr_ratio[part]->SetTitle("Difference");
    aset(gr_ratio[part], "p_{T} [GeV]","Difference", 0.,20., -0.1,0.1);
    style(gr_ratio[part], part+20, part+1);
    if(part==0)
      gr_ratio[part]->Draw("APE");
    else
      gr_ratio[part]->Draw("PE");
    leg0->AddEntry(gr_ratio[part], Form("%s",pname[part]), "P");
  }
  leg0->Draw();

  c0->Print("YieldCmp.pdf");
}
