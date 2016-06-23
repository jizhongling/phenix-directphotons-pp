void draw_Eta_Phi()
{
  TFile *f_sim = new TFile("Acceptance-signal.root");
  THnSparse *hn_photon = (THnSparse*)f_sim->Get("hn_photon");
  hn_photon->GetAxis(1)->SetRange(5,6);
  TAxis *axis0_hn_photon = hn_photon->GetAxis(0);

  TFile *f_data = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/total.root");
  THnSparse *hn_1photon = (THnSparse*)f_data->Get("hn_1photon");
  hn_1photon->GetAxis(0)->SetRange(5,6);
  TAxis *axis1_hn_1photon = hn_1photon->GetAxis(1);

  TCanvas *c1 = new TCanvas("c1", "#eta distribution", 2400, 2400);
  gStyle->SetOptStat(0);
  c1->Divide(4,4);

  TCanvas *c2 = new TCanvas("c2", "#phi distribution", 2400, 2400);
  gStyle->SetOptStat(0);
  c2->Divide(4,4);

  Int_t ipad = 1;
  for(Int_t ipt=11; ipt<27; ipt++)
  {
    Double_t low = axis0_hn_photon->GetBinLowEdge(ipt);
    Double_t high = axis0_hn_photon->GetBinUpEdge(ipt);

    axis0_hn_photon->SetRange(ipt,ipt);
    TH1 *h_eta = (TH1*)hn_photon->Projection(2)->Clone("h_eta");
    TH1 *h_phi = (TH1*)hn_photon->Projection(3)->Clone("h_phi");

    axis1_hn_1photon->SetRange(ipt,ipt);
    TH1 *h_1eta = (TH1*)hn_1photon->Projection(3)->Clone("h_1eta");
    TH1 *h_1phi = (TH1*)hn_1photon->Projection(4)->Clone("h_1phi");

    Double_t scale_eta = h_1eta->GetEntries() / h_eta->GetEntries();
    h_eta->Scale(scale_eta);
    Double_t scale_phi = h_1phi->GetEntries() / h_phi->GetEntries();
    h_phi->Scale(scale_phi);

    c1->cd(ipad);
    h_1eta->SetTitle(Form("p_{T}: %4.2f-%4.2f",low,high));
    h_1eta->SetMarkerSize(2.);
    h_1eta->SetMarkerStyle(2);
    h_1eta->SetMarkerColor(2);
    h_1eta->Draw("P");
    h_eta->Draw("SAME");

    c2->cd(ipad++);
    h_1phi->SetTitle(Form("p_{T}: %4.2f-%4.2f",low,high));
    h_1phi->SetMarkerSize(2.);
    h_1phi->SetMarkerStyle(2);
    h_1phi->SetMarkerColor(2);
    h_1phi->Draw("P");
    h_phi->Draw("SAME");
  }

  c1->Print("Eta-PbScE.pdf");
  c2->Print("Phi-PbScE.pdf");
}
