
int
compare_cones()
{
  gStyle->SetOptStat(0);

  TFile *f_directg = new TFile("keep/anaphpythia_directphoton.root", "OPEN");
  TFile *f_minbias = new TFile("keep/anaphpythia.root", "OPEN");

  THnSparse *hn_cone_directg = (THnSparse*)f_directg->Get("hn_EConeDirectPhoton");
  THnSparse *hn_cone_minbias = (THnSparse*)f_minbias->Get("hn_EConeOtherPhoton");

  TH1F* h1_econe_directg = (TH1F*)hn_cone_directg->Projection(2);
  TH1F* h1_econe_minbias = (TH1F*)hn_cone_minbias->Projection(2);

  h1_econe_directg->SetLineColor(kBlue);
  //h1_econe_directg->SetFillColor(kBlue);

  h1_econe_minbias->SetLineColor(kGray+1);
  h1_econe_minbias->SetFillColor(kGray+1);

  h1_econe_directg->Scale( 1. / h1_econe_directg->Integral(0,-1) );
  h1_econe_minbias->Scale( 1. / h1_econe_minbias->Integral(0,-1) );

  int bin_thr = h1_econe_directg->FindBin(0.1);

  TH1F* h1_threshold = (TH1F*)h1_econe_minbias->Clone("threshold");
  h1_threshold->Reset();
  h1_threshold->SetBinContent(bin_thr,1);
  h1_threshold->SetFillColor(kRed);
  h1_threshold->SetLineColor(kRed);

  double n_directg_total = h1_econe_directg->Integral(0,-1);
  double n_directg_below = h1_econe_directg->Integral(0,bin_thr);
  double n_directg_above = h1_econe_directg->Integral(bin_thr+1,-1);

  double n_minbias_total = h1_econe_minbias->Integral(0,-1);
  double n_minbias_below = h1_econe_minbias->Integral(0,bin_thr);
  double n_minbias_above = h1_econe_minbias->Integral(bin_thr+1,-1);

  cout << "Threshold bin: " << bin_thr << " at center " << h1_econe_minbias->GetBinCenter(bin_thr) << endl;

  cout << "Direct Photons: " << n_directg_total << " " << n_directg_below << " " << n_directg_above << endl;

  cout << "Other Photons:  " << n_minbias_total << " " << n_minbias_below << " " << n_minbias_above << endl;

  TCanvas *c1 = new TCanvas();
  c1->SetLogy(1);
  h1_econe_minbias->GetYaxis()->SetRangeUser(1e-4,1);
  h1_econe_minbias->Draw("");
  h1_econe_directg->Draw("same");
  h1_threshold->Draw("same");
  gPad->RedrawAxis();

  return 0;
}
