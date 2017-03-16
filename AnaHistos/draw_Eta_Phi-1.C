void draw_Eta_Phi_1()
{
  TFile *f_sim = new TFile("Acceptance-signal.root");
  THnSparse *hn_photon = (THnSparse*)f_sim->Get("hn_photon");
  //hn_photon->GetAxis(1)->SetRange(5,6);
  TAxis *axis0_hn_photon = hn_photon->GetAxis(0);

  TFile *f_1sim = new TFile("AnaPHPythia-histo.root");
  THnSparse *hn_1photon = (THnSparse*)f_1sim->Get("hn_photon");
  //hn_photon->GetAxis(1)->SetRange(5,6);
  TAxis *axis1_hn_1photon = hn_1photon->GetAxis(0);

  TCanvas *c1 = new TCanvas("c1", "#eta distribution", 2400, 2400);
  gStyle->SetOptStat(0);
  c1->Divide(4,4);

  TCanvas *c2 = new TCanvas("c2", "#phi distribution", 2400, 2400);
  gStyle->SetOptStat(0);
  c2->Divide(4,4);

  TCanvas *c3 = new TCanvas("c3", "#eta and #phi distribution", 2400, 2400);
  gStyle->SetOptStat(0);
  c3->Divide(4,4);

  TCanvas *c4 = new TCanvas("c4", "#eta and #phi distribution", 2400, 2400);
  gStyle->SetOptStat(0);
  c4->Divide(4,4);

  Int_t ipad = 1;
  for(Int_t ipt=11; ipt<27; ipt++)
  {
    Double_t low = axis0_hn_photon->GetBinLowEdge(ipt);
    Double_t high = axis0_hn_photon->GetBinUpEdge(ipt);

    axis0_hn_photon->SetRange(ipt,ipt);
    TH1 *h_eta = (TH1*)hn_photon->Projection(2)->Clone("h_eta");
    TH1 *h_phi = (TH1*)hn_photon->Projection(3)->Clone("h_phi");
    TH2 *h2_eta_phi = (TH2*)hn_photon->Projection(2,3)->Clone("h2_eta_phi");

    axis1_hn_1photon->SetRange(ipt,ipt);
    TH1 *h_1eta = (TH1*)hn_1photon->Projection(2)->Clone("h_1eta");
    TH1 *h_1phi = (TH1*)hn_1photon->Projection(3)->Clone("h_1phi");
    TH2 *h2_1eta_phi = (TH2*)hn_1photon->Projection(2,3)->Clone("h2_1eta_phi");

    Double_t scale = h_1eta->GetEntries() / h_eta->GetEntries();
    h_eta->Scale(scale);
    h_phi->Scale(scale);
    h2_eta_phi->Scale(scale);

    c1->cd(ipad);
    h_1eta->SetTitle(Form("p_{T}: %4.2f-%4.2f",low,high));
    h_1eta->SetMarkerSize(2.);
    h_1eta->SetMarkerStyle(2);
    h_1eta->SetMarkerColor(2);
    h_1eta->Draw("P");
    h_eta->Draw("SAME");

    c2->cd(ipad);
    h_1phi->SetTitle(Form("p_{T}: %4.2f-%4.2f",low,high));
    h_1phi->SetMarkerSize(2.);
    h_1phi->SetMarkerStyle(2);
    h_1phi->SetMarkerColor(2);
    h_1phi->Draw("P");
    h_phi->Draw("SAME");

    c3->cd(ipad);
    h2_eta_phi->SetTitle(Form("p_{T}: %4.2f-%4.2f",low,high));
    h2_eta_phi->Draw("COLZ");

    c4->cd(ipad);
    h2_1eta_phi->SetTitle(Form("p_{T}: %4.2f-%4.2f",low,high));
    h2_1eta_phi->Draw("COLZ");

    ipad++;
  }

  c1->Print("Eta.pdf");
  c2->Print("Phi.pdf");
  c3->Print("Eta-Phi-sim.pdf");
  c4->Print("Eta-Phi-data.pdf");
}
