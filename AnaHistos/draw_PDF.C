TGraph *RelErr(TGraphErrors *gr1, TGraphErrors *gr2)
{
  int gn = gr1->GetN();
  TGraph *graph = new TGraph(gn);
  int igp = 0;

  for(int i=0; i<gn; i++)
  {
    double g1x, g1y, g2x, g2y;
    gr1->GetPoint(i, g1x, g1y);
    gr2->GetPoint(i, g2x, g2y);
    double eg1y = gr1->GetErrorY(i);
    double eg2y = gr2->GetErrorY(i);

    if( eg1y > 0. && eg2y > 0. )
    {
      double gx = g1x;
      double gy = eg1y / eg2y;

      graph->SetPoint(igp, gx, gy);
      igp++;
    }
  }

  graph->Set(igp);
  return graph;
}

void draw_PDF()
{
  TGraphErrors *gr_old = new TGraphErrors("data/reweighting-old.txt", "%lg %lg %lg");
  TGraphErrors *gr_new = new TGraphErrors("data/reweighting-new.txt", "%lg %lg %lg");
  TGraph *gr_ratio = RelErr(gr_new, gr_old);

  mc(0, 2,1);

  mcd(0, 1);
  gPad->SetLogx();
  legi(0, 0.60,0.25,0.95,0.40);
  leg0->SetTextSize(0.030);

  gr_old->SetTitle("x#Deltag(x)");
  aset(gr_old, "x","x#Deltag(x)", 1e-3,1., -0.3,0.3);
  gr_old->GetXaxis()->SetNdivisions(505);
  gr_old->GetYaxis()->SetNdivisions(505);
  style(gr_old, 1, kBlack);
  style(gr_new, 1, kRed);
  gr_old->SetFillColor(kCyan-7);
  gr_new->SetFillColor(kRed);
  //gr_old->SetFillStyle(3001);
  gr_new->SetFillStyle(3004);
  gr_old->Draw("A3");
  gr_old->Draw("CX");
  gr_new->Draw("3 SAME");
  gr_new->Draw("CX SAME");
  leg0->AddEntry(gr_old, "DSSV14", "LF");
  leg0->AddEntry(gr_new, "Reweighting", "LF");
  leg0->Draw();

  mcd(0, 2);
  gPad->SetLogx();

  gr_ratio->SetTitle("Rel err");
  aset(gr_ratio, "x","Rel err", 1e-3,1.);
  style(gr_ratio, 1, kBlack);
  gr_ratio->Draw("AC");

  c0->Print("plots/PDFReweighting.pdf");
}
