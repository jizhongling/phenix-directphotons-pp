#include "DivideFunctions.h"

void draw_ProdRatio()
{
  // PYTHIA
  TFile *f[2];
  f[0] = new TFile("data/ProdRatio-histo-200gev.root");
  f[1] = new TFile("data/ProdRatio-histo-510gev.root");

  TH1 *h_prod[2][4];
  for(int ie=0; ie<2; ie++)
  {
    TH2 *h2_prod = (TH2*)f[ie]->Get("h2_prod");
    for(int id=0; id<4; id++)
      h_prod[ie][id] = h2_prod->ProjectionX(Form("h_meson%d_%d",ie,id), id+1,id+1);
  }

  // meson ratio
  TGraphErrors *gr_ratio1[2][3];
  for(int ie=0; ie<2; ie++)
    for(int id=0; id<3; id++)
      gr_ratio1[ie][id] = DivideHisto(h_prod[ie][id+1], h_prod[ie][0]);

  // r510/r200 ratio
  TGraphErrors *gr_ratio2[3];
  for(int id=0; id<3; id++)
    gr_ratio2[id] = DivideGraph(gr_ratio1[1][id], gr_ratio1[0][id]);

  // data for eta/meson [PHYSICAL REVIEW D 83, 032001 (2011)]
  TGraph *gr_eta = new TGraph("data/EtaPi0Ratio.txt");

  /* Tsallis parameters [PHYSICAL REVIEW D 83, 052004 (2011)]  
   *   x: pT
   * [0]: dsigma/dy
   * [1]: m
   * [2]: n
   * [3]: T
   */
  // par[bound][meson][ipar]
  const double par[3][4][4] = {
    { {41.4,0.135,9.57,0.1142}, {4.47,0.548,9.84,0.1190}, {3.65,0.782,10.00,0.1218}, {0.62,0.958,10.11,0.1238} },  // central
    { {35.6,0.135,9.67,0.1182}, {5.43,0.548,9.70,0.1136}, {4.42,0.782, 9.78,0.1151}, {0.79,0.958, 9.83,0.1161} },  // upper bound
    { {47.2,0.135,9.47,0.1102}, {3.51,0.548,9.98,0.1244}, {2.88,0.782,10.22,0.1285}, {0.45,0.958,10.39,0.1315} }   // lower bound
  };

  // Tsallis function
  TF1 *fn_prod[4][3];  // fn_prod[meson][bound]
  // meson id: pi0, eta, omega, eta'
  for(int id=0; id<4; id++)
    for(int iv=0; iv<3; iv++)
    {
      fn_prod[id][iv] = new TF1(Form("fn_prod%d_%d",id,iv), "[0]*([2]-1)*([2]-2)/([2]*[3]+[1]*([2]-1))/([2]*[3]+[1])*(([2]*[3]+sqrt(x**2+[1]**2))/([2]*[3]+[1]))**(-[2])", 2.,20.);
      fn_prod[id][iv]->SetParameters(par[iv][id]);
    }

  // Tsallis ratio
  TF1 *fn_ratio[3][3];  // fn_ratio[meson][bound]
  // meson id: eta, omega, eta' / pi0
  for(int id=0; id<3; id++)
    for(int iv=0; iv<3; iv++)
      fn_ratio[id][iv] = new TF1(Form("fn_ratio%d_%d",id,iv), Form("fn_prod%d_%d/fn_prod0_%d",id+1,iv,iv), 2.,20.);

  // convert Tsallis ratio functios to TGraphErrors to draw areas for uncertainties
  TGraphAsymmErrors *g_fn_ratio[3];  // g_fn_ratio[meson]
  for(int id=0; id<3; id++)
  {
    float xmin = 2;
    float xmax = 20;
    int xsteps = 100;

    g_fn_ratio[id] = new TGraphAsymmErrors( xsteps );

    for ( int step = 0; step < xsteps; step++ )
    {
      float x = xmin + step * (xmax-xmin)/((float)xsteps);
      float y = fn_ratio[id][0]->Eval(x);
      float exh = 0;
      float exl = 0;
      float eyh = fn_ratio[id][1]->Eval(x) - y;
      float eyl = y - fn_ratio[id][2]->Eval(x);

      g_fn_ratio[id]->SetPoint( step+1 , x , y );
      g_fn_ratio[id]->SetPointError( step+1 , exl, exh, eyl, eyh );
    }
  }

  // gSystem->ProcessLine(".L /phenix/u/zji/.rootrc.d/rootalias.C") to use this
  // create canvas
  mc(0, 2,1);

  for(int id=0; id<3; id++)
  {
    // cd to pad 1
    mcd(0, 1);
    // set axis and line stype
    gr_ratio1[0][id]->SetTitle("Prod. ratio for 200 GeV");
    aset(gr_ratio1[0][id], "p_{T} [GeV/c]","r_{200gev}", 2.,20., 0.,2.);
    style(gr_ratio1[0][id], id+20, id+1);
    gr_ratio1[0][id]->SetFillColor( gr_ratio1[0][id]->GetMarkerColor() );
    gr_ratio1[0][id]->SetFillStyle(3002);
    gr_ratio1[0][id]->SetLineColor( gr_ratio1[0][id]->GetMarkerColor() );
    gr_ratio1[0][id]->SetLineStyle(2);
    gr_ratio1[0][id]->SetLineWidth(2);
    //gr_ratio1[0][id]->SetMarkerStyle(1);

    if(id==0)
    {
      //gr_ratio1[0][id]->Draw("A3");
      gr_ratio1[0][id]->Draw("AXL");
    }
    else
    {
      //gr_ratio1[0][id]->Draw("3");
      gr_ratio1[0][id]->Draw("XL");
    }

    g_fn_ratio[id]->SetLineColor(id+1);
    g_fn_ratio[id]->SetLineWidth(3);
    g_fn_ratio[id]->SetFillColor(id+1);
    g_fn_ratio[id]->SetFillStyle(3002);
    g_fn_ratio[id]->Draw("3");
    g_fn_ratio[id]->Draw("XL");

    mcd(0, 2);
    gr_ratio2[id]->SetTitle("Ratio of prod. ratio for 510 to 200 GeV");
    aset(gr_ratio2[id], "p_{T} [GeV/c]","#frac{r_{510gev}}{r_{200gev}}", 2.,20., 0.5,2.);
    style(gr_ratio2[id], id+20, id+1);
    if(id==0)
      gr_ratio2[id]->Draw("AP");
    else
      gr_ratio2[id]->Draw("P");
  }

  mcd(0, 1);
  style(gr_eta, 20, 4);
  gr_eta->Draw("P");

  legi(0, 0.2,0.7,0.8,0.9);
  leg0->SetTextSize(0.03);
  leg0->SetNColumns(3);
  leg0->AddEntry(gr_ratio1[0][0], "Pythia #frac{#eta}{#pi^{0}}", "L");
  leg0->AddEntry(gr_ratio1[0][1], "Pythia #frac{#omega}{#pi^{0}}", "L");
  leg0->AddEntry(gr_ratio1[0][2], "Pythia #frac{#eta'}{#pi^{0}}", "L");
  leg0->AddEntry(g_fn_ratio[0], "Tsallis #frac{#eta}{#pi^{0}}", "L");
  leg0->AddEntry(g_fn_ratio[1], "Tsallis #frac{#omega}{#pi^{0}}", "L");
  leg0->AddEntry(g_fn_ratio[2], "Tsallis #frac{#eta'}{#pi^{0}}", "L");
  leg0->AddEntry(gr_eta, "Data #frac{#eta}{#pi^{0}}", "P");
  leg0->Draw();

  mcd(0, 2);
  legi(1, 0.2,0.8,0.9,0.9);
  leg0->SetTextSize(0.03);
  leg1->SetNColumns(3);
  leg1->AddEntry(gr_ratio2[0], "Pythia #frac{#eta}{#pi^{0}}", "P");
  leg1->AddEntry(gr_ratio2[1], "Pythia #frac{#omega}{#pi^{0}}", "P");
  leg1->AddEntry(gr_ratio2[2], "Pythia #frac{#eta'}{#pi^{0}}", "P");
  leg1->Draw();

  c0->Print("plots/ProdRatio.pdf");
}
