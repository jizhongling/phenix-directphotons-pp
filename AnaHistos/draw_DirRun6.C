void draw_DirRun6()
{
  const double PI = TMath::Pi();
  const double jetphox_scale = 1./200.;  // combined 200 histograms

  //TGraph *gr_run6 = new TGraph("data/run6-cross.txt");
  TGraph *gr_run6 = new TGraph("data/atlas-cross.txt");

  mc();
  mcd();

  TGraphErrors *gr_ratio = new TGraphErrors(20);
  int igp = 0;

  //TFile *f_th = new TFile("data/isoprompt-run6.root");
  TFile *f_th = new TFile("data/isoprompt-atlas.root");
  TH1 *h_th = (TH1*)f_th->Get("hp41");
  h_th->Scale(jetphox_scale);

  //for(int ipt=0; ipt<18; ipt++)
  for(int ipt=0; ipt<8; ipt++)
  {
    double pT, run6;
    gr_run6->GetPoint(ipt, pT, run6);

    //double factor = 1. / (2*PI*pT*0.7);
    double factor = 1.;
    double nth = factor * h_th->GetBinContent(ipt+1);
    double enth = factor * h_th->GetBinError(ipt+1);

    double ratio = run6 / nth;
    double eratio = ratio * enth / nth;

    gr_ratio->SetPoint(igp, pT, ratio);
    gr_ratio->SetPointError(igp, 0., eratio);
    igp++;
  }

  gr_ratio->Set(igp);
  //aset(gr_ratio, "p_{T} [GeV]","#frac{data}{theory}");
  aset(gr_ratio, "p_{T} [GeV]","#frac{atlas}{mine}");
  style(gr_ratio, 20, 1);
  gr_ratio->Draw("AP");
  gr_ratio->Fit("pol0", "Q");

  //c0->Print("plots/DirRun6.pdf");
  c0->Print("plots/DirATLAS.pdf");
}
