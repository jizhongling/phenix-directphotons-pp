void draw_pTDistro()
{
  gROOT->ProcessLine(".L ReadGraph.C");
  const Int_t trig = 2;
  const Int_t fig = 2;

  Double_t gx[30];
  Double_t gy[8][30];
  Double_t egy[8][30];
  for(Int_t im=0; im<1; im++)
    for(Int_t part=0; part<4; part++)
      ReadGraphErrors(Form("CrossSection-ert%c.root",97+trig), 36*im+12*trig+4*fig+part, gx, (Double_t*)gy[3*im+part], (Double_t*)egy[3*im+part]);

  Double_t Sasha_pT[30] = {0.25, 0.75, 1.1907e+0, 1.7345e+0, 2.1877e+0, 2.5953e+0, 3.2734e+0, 3.7256e+0, 4.1778e+0, 4.6748e+0,
    5.2618e+0, 5.7137e+0, 6.2555e+0, 6.7071e+0, 7.2488e+0, 7.7003e+0, 8.2418e+0, 8.7833e+0, 9.2798e+0, 9.8210e+0,
    1.0949e+1, 1.2844e+1, 1.4874e+1, 1.6903e+1, 1.8887e+1, 2.0961e+1, 2.2990e+1, 2.4929e+1, 2.7003e+1, 2.8987e+1};
  Double_t Sasha[30] = {1., 1., 7.2450e-1, 1.2425e-1, 2.9580e-2, 9.0062e-3, 2.7409e-3, 1.1117e-3, 4.5093e-4, 2.1547e-4,
    1.0725e-4, 5.1250e-5, 2.8849e-5, 1.6242e-5, 9.5250e-6, 5.5868e-6, 3.7051e-6, 2.4572e-6, 1.5014e-6, 1.1260e-6,
    4.7531e-7, 1.3858e-7, 4.5680e-8, 2.0061e-8, 8.4569e-9, 3.7137e-9, 1.8443e-9, 1.0792e-9, 5.3592e-10, 2.8890e-10};
  for(Int_t ipt=0; ipt<30; ipt++)
    Sasha[ipt] *= 1e9;

  TGraphErrors *gr = new TGraphErrors(30, gx, gy[1], 0, egy[1]);

  TF1 *fn1 = new TF1("fn1", "expo", 5., 20.);
  TF1 *fn2 = new TF1("fn2", "expo(0)+expo(2)", 5., 20.);
  //TF1 *fn2 = new TF1("fn2", "[0]*TMath::Exp([3]+[1]*x+[2]*x**2)", 5., 20.);

  mc(0, 2,1);

  mcd(0, 1);
  gPad->SetLogy();
  gr->SetTitle("PbScE;p_{T} [GeV];d#sigma;");
  aset(gr, "","", 0.,30., 1.,2e4);
  style(gr, 21, 2);
  gr->Fit(fn1, "QR0");
  fn2->SetParameters( fn1->GetParameters() );
  fn2->SetLineColor(3);
  gr->Fit(fn2, "QR0");
  gr->Fit(fn2, "QRE");
  gr->Draw("APE");

  TGraphErrors *gr_ratio[6];
  Double_t rgy[6][30] = {};
  Double_t ergy[6][30] = {};

  for(Int_t im=0; im<1; im++)
    for(Int_t part=0; part<3; part++)
    {
      for(Int_t ipt=0; ipt<30; ipt++)
      {
        rgy[3*im+part][ipt] = ( gy[3*im+part][ipt] - fn2->Eval(gx[ipt]) ) / fn2->Eval(gx[ipt]);
        ergy[3*im+part][ipt] = egy[3*im+part][ipt] / fn2->Eval(gx[ipt]); 
        //rgy[3*im+part][ipt] = ( gy[3*im+part][ipt] - gy[3*im+3][ipt] ) / gy[3*im+3][ipt];
        //ergy[3*im+part][ipt] = sqrt( pow(egy[3*im+part][ipt]/gy[3*im+part][ipt],2.) + pow(egy[3*im+3][ipt]/gy[3*im+3][ipt],2.) ); 
        //rgy[3*im+part][ipt] = ( gy[3*im+part][ipt] - Sasha[ipt] ) / Sasha[ipt];
        //ergy[3*im+part][ipt] = egy[3*im+part][ipt] / Sasha[ipt];
      }
      gr_ratio[3*im+part] =  new TGraphErrors(30, gx, rgy[3*im+part], 0, ergy[3*im+part]);
    }

  mcd(0, 2);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  const char *pname[3] = {"PbScW", "PbScE", "PbGlE"};
  for(Int_t im=0; im<1; im++)
    for(Int_t part=0; part<3; part++)
    {
      gr_ratio[3*im+part]->SetTitle("Ratio;p_{T} [GeV];ratio;");
      aset(gr_ratio[3*im+part], "","", 5.,20., -0.3,0.3);
      style(gr_ratio[3*im+part], 20+part+4*im, 1+part);
      if(im==0 && part==0)
        gr_ratio[3*im+part]->Draw("APE");
      else
        gr_ratio[3*im+part]->Draw("PE");
      if(im==0)
        leg0->AddEntry(gr_ratio[part], Form("%s",pname[part]), "P");
    }
  leg0->Draw();

  c0->Print(Form("pTDistro-ert%c-fig%d.pdf",97+trig,fig));
}
