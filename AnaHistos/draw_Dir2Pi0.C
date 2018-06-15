void draw_Dir2Pi0()
{
  const Double_t PI = TMath::Pi();
  const Double_t jetphox_scale = 1./200.;  // combined 200 histograms
  const char *fname[3] = {"halfpt", "onept", "twopt"};

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  mc();
  mcd();

  for(Int_t imu=0; imu<3; imu++)
  {
    TGraphErrors *gr_ratio = new TGraphErrors(20);
    Int_t igp = 0;

    TFile *f_ph = new TFile( Form("data/isoprompt-%s.root",fname[imu]) );
    TH1 *h_ph = (TH1*)f_ph->Get("hp41");
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

    gr_ratio->Set(igp);
    aset(gr_ratio, "p_{T} [GeV]","#gamma_{dir}/#pi^{0}");
    style(gr_ratio, 20+imu, 1+imu);
    if(imu==0)
      gr_ratio->Draw("AP");
    else
      gr_ratio->Draw("P");

    delete h_ph;
    delete f_ph;
  }

  c0->Print("plots/Dir2Pi0-NLO.pdf");
}
