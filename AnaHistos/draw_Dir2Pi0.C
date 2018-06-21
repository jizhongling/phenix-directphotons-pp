void draw_Dir2Pi0()
{
  const Double_t PI = TMath::Pi();
  const Double_t jetphox_scale = 1./200.;  // combined 200 histograms
  const char *fname[3] = {"halfpt", "onept", "twopt"};
  const char *gtitle[2] = {"LO", "NLO"};

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  mc(0, 2,1);
  legi(0, 0.2,0.7,0.4,0.9);

  for(Int_t io=0; io<2; io++)
    for(Int_t imu=0; imu<3; imu++)
    {
      TGraphErrors *gr_ratio = new TGraphErrors(20);
      Int_t igp = 0;

      TFile *f_ph = new TFile( Form("data/isoprompt-%s.root",fname[imu]) );
      TH1 *h_ph = (TH1*)f_ph->Get(Form("hp4%d",io));
      h_ph->Scale(jetphox_scale);

      for(Int_t ipt=12; ipt<30; ipt++)
      {
        Double_t sasha_pT, sasha;
        gr_sasha->GetPoint(ipt-2, sasha_pT, sasha);
        Double_t npi0 = 2*PI*sasha_pT*0.7*sasha;

        Double_t ratio = h_ph->GetBinContent(igp+1) / npi0;
        Double_t eratio = h_ph->GetBinError(igp+1) / npi0;

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
