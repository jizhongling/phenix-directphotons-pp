void draw_pTDistro()
{
  gROOT->ProcessLine(".L ReadGraph.C");
  const Int_t trig = 1;

  const Double_t gx[30];
  Double_t gy[3][30];
  Double_t egy[3][30];
  for(Int_t part=0; part<3; part++)
    ReadGraphErrors("CrossSection-ertb.root", 12*trig+4+part, gx, (Double_t*)gy[part], (Double_t*)egy[part]);

  TGraphErrors *gr = new TGraphErrors(30, gx, gy[0], 0, egy[0]);

  TF1 *fn1 = new TF1("fn1", "expo", 5., 20.);
  //TF1 *fn2 = new TF1("fn2", "[0]*TMath::Exp([3]+[1]*x+[2]*x**2)", 5., 20.);
  TF1 *fn2 = new TF1("fn2", "expo(0)+expo(2)", 5., 20.);

  TCanvas *c = new TCanvas("c", "Canvas", 1200, 600);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  c->Divide(2,1);

  c->cd(1);
  gPad->SetLogy();
  gr->SetTitle("PbScW;p_{T} [GeV];d#sigma;");
  gr->GetXaxis()->SetRangeUser(5., 20.);
  gr->GetYaxis()->SetRangeUser(1., 2e4);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(1);
  //gr->SetMarkerSize(0.01);
  gr->Fit(fn1, "QR0");
  fn2->SetParameters( fn1->GetParameters() );
  gr->Fit(fn2, "QR0");
  gr->Fit(fn2, "QRE");
  gr->Draw("APE");

  TGraphErrors *gr_ratio[3];
  Double_t rgy[3][30] = {};
  Double_t ergy[3][30] = {};

  for(Int_t part=0; part<3; part++)
  {
    for(Int_t ipt=11; ipt<26; ipt++)
    {
      rgy[part][ipt] = ( gy[part][ipt] - fn2->Eval(gx[ipt]) ) / fn2->Eval(gx[ipt]);
      ergy[part][ipt] = rgy[part][ipt] * egy[part][ipt] / gy[part][ipt];
    }
    gr_ratio[part] =  new TGraphErrors(30, gx, rgy[part], 0, ergy[part]);
  }

  c->cd(2);
  for(Int_t part=0; part<3; part++)
  {
    gr_ratio[part]->SetTitle("Ratio;p_{T} [GeV];ratio;");
    gr_ratio[part]->GetXaxis()->SetRangeUser(5., 20.);
    gr_ratio[part]->GetYaxis()->SetRangeUser(-1., 1.);
    gr_ratio[part]->SetMarkerStyle(20+part);
    gr_ratio[part]->SetMarkerColor(1+part);
    //gr_ratio[part]->SetMarkerSize(0.01);
    if(part==0)
      gr_ratio[part]->Draw("APE");
    else
      gr_ratio[part]->Draw("PE");
  }

  c->Print("pTDistro.pdf");
}
