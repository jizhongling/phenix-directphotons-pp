void draw_Theta_CV_Graph()
{
  const int n = 15;
  double x[n];
  for(int i=0; i<n; i++)
    x[i] = 0.32 + 0.08*i;
  double y1[n] = {0.010, 0.010, 0.010, 0.010, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014};
  double y2[n] = {0.050, 0.046, 0.042, 0.034, 0.026, 0.022, 0.026, 0.026, 0.026, 0.022, 0.022, 0.018, 0.018, 0.018, 0.018};

  TGraph *gr1 = new TGraph(n, x, y1);
  TGraph *gr2 = new TGraph(n, x, y2);
  TGraph *gr3 = new TGraph(2*n);
  for(int i=0; i<n; i++)
  {
    gr3->SetPoint(2*i, x[i], y1[i]);
    gr3->SetPoint(2*i+1, x[i], y2[i]);
  }

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);

  gr1->GetXaxis()->SetLabelSize(0.02);
  gr1->GetXaxis()->SetTitleSize(0.02);
  gr1->GetXaxis()->SetTitle("E_{cluster} (GeV)");
  gr1->GetYaxis()->SetLabelSize(0.02);
  gr1->GetYaxis()->SetTitleSize(0.02);
  gr1->GetYaxis()->SetTitle("#theta_{CV} [rad]");
  gr1->GetYaxis()->SetRangeUser(0., 0.06);
  gr1->SetLineColor(1);
  gr1->SetMarkerStyle(22);
  gr1->SetMarkerSize(1);
  gr1->Draw("ACP");

  gr2->SetLineColor(2);
  gr2->SetMarkerStyle(23);
  gr2->SetMarkerColor(2);
  gr2->SetMarkerSize(1);
  gr2->Draw("CP");

  gr3->SetFillStyle(3004);
  gr3->SetFillColor(kBlue);
  gr3->Draw("F");

  c->Print("plots/Theta_CV_PbSc_Graph.pdf");
}
