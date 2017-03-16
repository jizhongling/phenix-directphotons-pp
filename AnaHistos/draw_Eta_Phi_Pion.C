void draw_Eta_Phi_Pion()
{
  gROOT->ProcessLine(".L ReadGraph.C");
  TH1::SetDefaultSumw2();

  TFile *f_sim = new TFile("MissingRatio-histo.root");
  THnSparse *hn_photon = (THnSparse*)f_sim->Get("hn_pion");
  TAxis *axis_pt_hn_photon = hn_photon->GetAxis(1);

  //TFile *f_data = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ertc-cv/total.root");
  TFile *f_data = new TFile("/phenix/plhf/zji/taxi/Run13pp510MinBias/10677/data/total.root");
  THnSparse *hn_1photon = (THnSparse*)f_data->Get("hn_pion");
  TAxis *axis_pt_hn_1photon = hn_1photon->GetAxis(1);

  Double_t gx[30], TrigE[2][30], eTrigE[2][30];
  Int_t ispion = 1;
  Int_t trig = 2;
  for(Int_t part=0; part<2; part++)
    ReadGraphAsymmErrors("TriggerEfficiency.root", 9*ispion+3*trig+part, gx, (Double_t*)TrigE[part], (Double_t*)eTrigE[part]);

  TCanvas *c1 = new TCanvas("c1", "#eta distribution", 2400, 2400);
  gStyle->SetOptStat(0);
  c1->Divide(4,4);

  TCanvas *c2_1 = new TCanvas("c2_1", "#phi distribution", 2400, 2400);
  gStyle->SetOptStat(0);
  c2_1->Divide(4,4);

  TCanvas *c2_2 = new TCanvas("c2_2", "#phi distribution", 2400, 2400);
  gStyle->SetOptStat(0);
  c2_2->Divide(4,4);

  TCanvas *c3 = new TCanvas("c3", "#eta and #phi distribution", 2400, 2400);
  gStyle->SetOptStat(0);
  c3->Divide(4,4);

  TCanvas *c4 = new TCanvas("c4", "#eta and #phi distribution", 2400, 2400);
  gStyle->SetOptStat(0);
  c4->Divide(4,4);

  Int_t secl[2] = {1, 7};
  Int_t sech[2] = {6, 8};

  //Double_t n_1phi = 0.;
  //for(Int_t ipt=11; ipt<30; ipt++)
  //{
  //  axis_pt_hn_1photon->SetRange(ipt,ipt);
  //  for(Int_t part=0; part<2; part++)
  //  {
  //    hn_1photon->GetAxis(0)->SetRange(secl[part],sech[part]);
  //    TH1 *h_1phi = (TH1*)hn_1photon->Projection(3)->Clone("h_1phi");
  //    if( TrigE[part][ipt] > 0. )
  //      n_1phi += h_1phi->GetEntries() / TrigE[part][ipt];
  //  }
  //}
  //
  //axis_pt_hn_photon->SetRange(11,29);
  //TH1 *h_phi = (TH1*)hn_photon->Projection(5)->Clone("h_phi");
  //Double_t scale0 = n_1phi / h_phi->GetEntries();

  Int_t ipad = 1;
  //for(Int_t ipt=11; ipt<27; ipt++)
  for(Int_t ipt=5; ipt<30; ipt+=3)
  {
    Double_t low = axis_pt_hn_photon->GetBinLowEdge(ipt);
    Double_t high = axis_pt_hn_photon->GetBinUpEdge(ipt+2);

    axis_pt_hn_photon->SetRange(ipt,ipt+2);
    hn_photon->GetAxis(2)->SetRange(112,162);
    TH1 *h_eta = (TH1*)hn_photon->Projection(4)->Clone("h_eta");
    TH1 *h_phi = (TH1*)hn_photon->Projection(5)->Clone("h_phi");
    TH2 *h2_eta_phi = (TH2*)hn_photon->Projection(4,5)->Clone("h2_eta_phi");

    axis_pt_hn_1photon->SetRange(ipt,ipt+2);
    hn_1photon->GetAxis(2)->SetRange(112,162);
    TH1 *h_1eta = (TH1*)hn_1photon->Projection(3)->Clone("h_1eta");
    TH1 *h_1phi = (TH1*)hn_1photon->Projection(4)->Clone("h_1phi");
    TH2 *h2_1eta_phi = (TH2*)hn_1photon->Projection(3,4)->Clone("h2_1eta_phi");

    TH1 *hv_1phi[2];
    TH1 *hv_phi[2];
    Double_t n_1phi = 0.;
    Double_t n_phi = 0.;
    for(Int_t part=0; part<2; part++)
    {
      hn_1photon->GetAxis(0)->SetRange(secl[part],sech[part]);
      hv_1phi[part] = (TH1*)hn_1photon->Projection(4)->Clone("h_1phi");
      n_1phi += hv_1phi[part]->GetEntries(); // / TrigE[part][ipt];
      hn_photon->GetAxis(3)->SetRange(secl[part],sech[part]);
      hv_phi[part] = (TH1*)hn_photon->Projection(5)->Clone("h_phi");
      n_phi += hv_phi[part]->GetEntries();
    }
    for(Int_t part=0; part<2; part++)
    {
      //Double_t scale = hv_1phi[part]->GetEntries() / hv_phi[part]->GetEntries();
      hv_phi[part]->Scale( n_1phi / n_phi );
      //if( TrigE[part][ipt] > 0. )
      //  hv_1phi[part]->Scale(1./TrigE[part][ipt]);
      hv_1phi[part]->SetMarkerSize(2.);
      hv_1phi[part]->SetMarkerStyle(2);
      hv_1phi[part]->SetMarkerColor(2);
      hv_1phi[part]->SetTitle(Form("p_{T}: %4.2f-%4.2f",low,high));
      hv_phi[part]->SetTitle(Form("p_{T}: %4.2f-%4.2f",low,high));
    }

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

    c2_1->cd(ipad);
    hv_1phi[0]->GetXaxis()->SetRangeUser(-0.6,0.95);
    hv_phi[0]->GetXaxis()->SetRangeUser(-0.6,0.95);
    hv_phi[0]->DrawCopy("");
    hv_1phi[0]->DrawCopy("PSAME");

    c2_2->cd(ipad);
    hv_1phi[1]->GetXaxis()->SetRangeUser(2.2,4.);
    hv_phi[1]->GetXaxis()->SetRangeUser(2.2,4.);
    hv_1phi[0]->GetXaxis()->SetRangeUser(2.2,4.);
    hv_phi[0]->GetXaxis()->SetRangeUser(2.2,4.);
    hv_1phi[1]->DrawCopy("P");
    hv_1phi[0]->DrawCopy("PSAME");
    for(Int_t part=0; part<2; part++)
      hv_phi[part]->DrawCopy("SAME");

    c3->cd(ipad);
    h2_eta_phi->SetTitle(Form("p_{T}: %4.2f-%4.2f",low,high));
    h2_eta_phi->Draw("COLZ");

    c4->cd(ipad);
    h2_1eta_phi->SetTitle(Form("p_{T}: %4.2f-%4.2f",low,high));
    h2_1eta_phi->Draw("COLZ");

    ipad++;
  }

  c1->Print("Eta.pdf");
  c2_1->Print("Phi-pion-1.pdf");
  c2_2->Print("Phi-pion-2.pdf");
  c3->Print("Eta-Phi-sim.pdf");
  c4->Print("Eta-Phi-data.pdf");
}
