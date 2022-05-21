void draw_PDF()
{
  TGraphErrors *gr_old = new TGraphErrors("data/reweighting-old.txt", "%lg %lg %lg");
  TGraphErrors *gr_new = new TGraphErrors("data/reweighting-new.txt", "%lg %lg %lg");

  mc();
  mcd();
  gPad->SetLogx();
  legi(0, 0.60,0.25,0.95,0.40);
  leg0->SetTextSize(0.030);

  gr_old->SetTitle("");
  aset(gr_old, "x","x#Deltag(x)", 1e-5,1., -0.3,0.3);
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

  c0->Print("plots/PDFReweighting.pdf");
}
