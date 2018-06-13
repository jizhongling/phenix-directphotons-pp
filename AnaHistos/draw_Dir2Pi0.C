void draw_Dir2Pi0()
{
  const Double_t norm_dir = 0.47995E+06 / (200*100000);
  const Double_t norm_onef = 0.55060E+05 / (200*100000);

  TGraphErrors *gr_ratio = new TGraphErrors(20);
  Int_t igp = 0;

  TFile *f_dir = new TFile("data/dir.root");
  TFile *f_onef = new TFile("data/onef.root");
  TH1 *h_dir = (TH1*)f_dir->Get("hp41");
  TH1 *h_onef = (TH1*)f_onef->Get("hp41");

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  for(Int_t ipt=12; ipt<30; ipt++)
  {
    Double_t ndir = norm_dir * h_dir->GetBinContent(igp+1);
    Double_t nonef = norm_onef * h_onef->GetBinContent(igp+1);
    Double_t endir = norm_dir * h_dir->GetBinError(igp+1);
    Double_t enonef = norm_onef * h_onef->GetBinError(igp+1);
    Double_t ngamma = ndir + nonef;
    Double_t engamma = sqrt(endir*endir + enonef*enonef) / sqrt(norm_dir*norm_dir + norm_onef*norm_onef);

    Double_t sasha_pT, sasha;
    gr_sasha->GetPoint(ipt-2, sasha_pT, sasha);
    Double_t npi0 = 2*TMath::Pi()*sasha_pT*0.7*sasha;

    gr_ratio->SetPoint(igp, sasha_pT, ngamma/npi0);
    gr_ratio->SetPointError(igp, 0., engamma/npi0);
    igp++;
  }
  gr_ratio->Set(igp);

  mc();
  mcd();
  aset(gr_ratio, "p_{T} [GeV]","#gamma_{dir}/#pi^{0}", 4.,32., 0.,1.5);
  style(gr_ratio, 20, 1);
  gr_ratio->Draw("AP");
  c0->Print("plots/Dir2Pi0-NLO.pdf");
}
