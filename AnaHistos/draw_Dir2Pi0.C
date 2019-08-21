void draw_Dir2Pi0()
{
  const double PI = TMath::Pi();
  const double jetphox_scale = 1./400.;  // combined 400 histograms
  const char *fname[3] = {"halfpt", "onept", "twopt"};
  const char *gtitle[2] = {"LO", "NLO"};

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  mc(0, 2,1);
  legi(0, 0.2,0.7,0.4,0.9);

  for(int io=0; io<2; io++)
    for(int imu=0; imu<3; imu++)
    {
      TGraphErrors *gr_ratio = new TGraphErrors(20);
      int igp = 0;

      TFile *f_ph = new TFile( Form("data/isoprompt-x400-ct14-%s.root",fname[imu]) );
      TH1 *h_ph = (TH1*)f_ph->Get(Form("hp4%d",io));
      h_ph->Scale(jetphox_scale);

      for(int ipt=12; ipt<30; ipt++)
      {
        double sasha_pT, sasha;
        gr_sasha->GetPoint(ipt-2, sasha_pT, sasha);

        double factor = 1. / (2*PI*sasha_pT*0.5);
        int bin_th = h_ph->GetXaxis()->FindBin(sasha_pT);
        double ratio = factor * h_ph->GetBinContent(bin_th) / sasha;
        double eratio = factor * h_ph->GetBinError(bin_th) / sasha;

        gr_ratio->SetPoint(igp, sasha_pT, ratio);
        gr_ratio->SetPointError(igp, 0., eratio);
        igp++;
      }

      mcd(0, io+1);
      gr_ratio->Set(igp);
      gr_ratio->SetTitle(gtitle[io]);
      aset(gr_ratio, "p_{T} [GeV]","#gamma_{dir}/#pi^{0}");
      style(gr_ratio, 20+imu, 1+imu);
      if(imu==0)
        gr_ratio->Draw("AP");
      else
        gr_ratio->Draw("P");
      if(io==0)
      {
        leg0->AddEntry(gr_ratio, Form("#mu=%s",fname[imu]), "P");
        leg0->Draw();
      }

      delete h_ph;
      delete f_ph;
    }

  c0->Print("plots/Dir2Pi0.pdf");
}
