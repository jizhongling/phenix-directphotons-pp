#include "DivideFunctions.h"

void draw_ProdRatio()
{
  TFile *f = new TFile("data/ProdRatio-histo-200gev.root");
  TH2 *h2_prod = (TH2*)f->Get("h2_prod");

  TGraphErrors *gr_ratio[3];
  TH1 *h_pion = h2_prod->ProjectionX("h_pion", 1,1);

  /*   x: pT
   * [0]: dsigma/dy
   * [1]: m
   * [2]: n
   * [3]: T
   */
  TF1 *fn_prod[4];
  const Double_t par[4][4] = { {41.4,0.135,9.57,0.1142}, {4.47,0.548,9.84,0.1190}, {3.65,0.782,10.00,0.1218}, {0.62,0.958,10.11,0.1238} };
  for(Int_t id=0; id<4; id++)
  {
    fn_prod[id] = new TF1(Form("fn_prod%d",id), "[0]*([2]-1)*([2]-2)/([2]*[3]+[1]*([2]-1))/([2]*[3]+[1])*(([2]*[3]+sqrt(x**2+[1]**2))/([2]*[3]+[1]))**(-[2])", 0.,30.);
    fn_prod[id]->SetParameters(par[id]);
  }

  TF1 *fn_ratio[3];
  for(Int_t id=0; id<3; id++)
    fn_ratio[id] = new TF1(Form("fn_ratio%d",id), Form("fn_prod%d/fn_prod0",id+1), 0.,30.);

  mc();
  mcd();

  for(Int_t id=0; id<3; id++)
  {
    TH1 *h_other = h2_prod->ProjectionX("h_other", id+2,id+2);
    gr_ratio[id] = DivideHisto(h_other, h_pion);

    aset(gr_ratio[id], "","", 0.,30., 0.,2.);
    style(gr_ratio[id], id+20, id+1);
    if(id==0)
      gr_ratio[id]->Draw("AP");
    else
      gr_ratio[id]->Draw("P");

    fn_ratio[id]->SetLineColor(id+1);
    fn_ratio[id]->Draw("SAME");

    delete h_other;
  }

  c0->Print("plots/ProdRatio-200gev.pdf");
}
