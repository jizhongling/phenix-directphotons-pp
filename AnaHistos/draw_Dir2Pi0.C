void draw_Dir2Pi0()
{
  const Double_t PI = TMath::Pi();
  const char *fname[3] = {"halfpt", "onept", "twopt"};
  const Double_t unit_dir[3] = {0.62364E+06, 0.47995E+06, 0.37847E+06};
  const Double_t unit_onef[3] = {0.76073E+05, 0.55060E+05, 0.42531E+05};

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  mc();
  mcd();

  for(Int_t imu=0; imu<3; imu++)
  {
    TGraphErrors *gr_ratio = new TGraphErrors(18);

    TFile *f_dir = new TFile( Form("data/dir-%s.root",fname[imu]) );
    TFile *f_onef = new TFile( Form("data/onef-%s.root",fname[imu]) );
    TH1 *h_dir = (TH1*)f_dir->Get("hp41");
    TH1 *h_onef = (TH1*)f_onef->Get("hp41");

    Double_t norm_dir = unit_dir[imu] / (200*100000);
    Double_t norm_onef = unit_onef[imu] / (200*100000);

    Int_t igp = 0;
    for(Int_t ipt=12; ipt<30; ipt++)
    {
      Double_t ndir = norm_dir * h_dir->GetBinContent(igp+1);
      Double_t nonef = norm_onef * h_onef->GetBinContent(igp+1);
      Double_t endir = norm_dir * h_dir->GetBinError(igp+1);
      Double_t enonef = norm_onef * h_onef->GetBinError(igp+1);
      Double_t ngamma = ndir + nonef;
      Double_t engamma = sqrt( endir*endir + enonef*enonef );

      Double_t sasha_pT, sasha;
      gr_sasha->GetPoint(ipt-2, sasha_pT, sasha);
      Double_t npi0 = 2*PI*sasha_pT*0.7*sasha;

      gr_ratio->SetPoint(igp, sasha_pT, ngamma/npi0);
      gr_ratio->SetPointError(igp, 0., engamma/npi0);
      igp++;
    }

    aset(gr_ratio, "p_{T} [GeV]","#gamma_{dir}/#pi^{0}");
    style(gr_ratio, 20+imu, 1+imu);
    if(imu==0)
      gr_ratio->Draw("AP");
    else
      gr_ratio->Draw("P");

    delete h_dir;
    delete h_onef;
    delete f_dir;
    delete f_onef;
  }

  c0->Print("plots/Dir2Pi0-NLO.pdf");
}
