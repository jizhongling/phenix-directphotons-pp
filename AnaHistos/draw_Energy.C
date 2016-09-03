void draw_Energy()
{
  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ertb/total.root");
  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");

  hn_1photon->GetAxis(0)->SetRange(1,4);
  TAxis *axis_e = hn_1photon->GetAxis(6);

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  TLegend *leg = new TLegend(0.7, 0.5, 0.9, 0.9);

  Int_t icon = 1;
  for(Int_t i=1; i<10; i++)
  {
    axis_e->SetRange(i+1, 30);
    Double_t E = axis_e->GetBinLowEdge(i+1);
    TH1 *h_1photon = hn_1photon->Projection(1)->Clone(Form("h_%d",i));
    h_1photon->GetXaxis()->SetRangeUser(0., 10.);
    h_1photon->GetYaxis()->SetRangeUser(1e4, 1e7);
    h_1photon->SetMarkerColor(icon);
    h_1photon->SetMarkerStyle(icon+20);
    if(icon == 1)
      h_1photon->Draw("P");
    else
      h_1photon->Draw("PSAME");
    leg->AddEntry(h_1photon, Form("E>%4.2f",E), "P");
    icon++;
  }
  leg->Draw();

  c->Print("Energy.pdf");
}
