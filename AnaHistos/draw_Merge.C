TGraphErrors *DivideHisto(TH1 *h1, TH1 *h2)
{
  Int_t gn = h1->GetXaxis()->GetNbins();
  Double_t *gx = new Double_t[gn];
  Double_t *gy = new Double_t[gn];
  Double_t *egx = new Double_t[gn];
  Double_t *egy = new Double_t[gn];

  for(Int_t i=0; i<gn; i++)
  {
    gx[i] = gy[i] = egx[i] = egy[i] = 0.;

    Double_t h1x = h1->GetXaxis()->GetBinCenter(i+1);
    Double_t h1y = h1->GetBinContent(i+1); 
    Double_t h2y = h2->GetBinContent(i+1); 
    Double_t eh1x = h1->GetXaxis()->GetBinWidth(i+1) / 2.;
    Double_t eh1y = h1->GetBinError(i+1);
    Double_t eh2y = h2->GetBinError(i+1);

    if( h1y > 0. && h2y > 0. )
    {
      gx[i] = h1x;
      gy[i] = h1y / h2y;
      egx[i] = eh1x;
      egy[i] = gy[i] * sqrt( pow(eh1y/h1y,2.) + pow(eh2y/h2y,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(gn, gx, gy, egx, egy);
  return graph;
}

void InvertGraph(TGraphAsymmErrors *gr)
{
  Int_t gn = gr->GetN();
  for(Int_t i=0; i<gn; i++)
  {
    Double_t gx, gy;
    gr->GetPoint(i, gx, gy);
    gr->SetPoint(i, gx, 1.-gy);
  }

  return;
}

void draw_Merge()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/pi0cross_run13pp510gev/fastMC/histos/total.root");
  TFile *f1 = new TFile("MissingRatio-histo.root");
  //TFile *f1 = new TFile("Pi0Sim-histo.root");
  //TFile *f1 = new TFile("MergeRate-histo.root");
  //TFile *f1 = new TFile("AnaPHPythia-histo.root");
  TObjArray *Glist = new TObjArray();

  TH1::SetDefaultSumw2();

  TH2 *h2_ph_asym = (TH2*)f->Get("h2_ph_asym");
  TH2 *h2_ph_emc_asym_noprob = (TH2*)f->Get("h2_ph_emc_asym_noprob");
  TH2 *h2_pi0_asym = (TH2*)f->Get("h2_pi0_asym");
  TH2 *h2_pi0_emc_asym_noprob = (TH2*)f->Get("h2_pi0_emc_asym_noprob");

  TH2 *h2_merge_photon = (TH2*)f1->Get("h2_merge_photon");
  TH2 *h2_merge_pion = (TH2*)f1->Get("h2_merge_pion");
  TH2 *h2_total_photon = (TH2*)f1->Get("h2_total_photon");
  TH2 *h2_total_pion = (TH2*)f1->Get("h2_total_pion");

  //TH2 *h2_merge_photon = (TH2*)f1->Get("h2_merge");
  //TH2 *h2_merge_pion = (TH2*)f1->Get("h2_merge");
  //TH2 *h2_total_photon = (TH2*)f1->Get("h2_total");
  //TH2 *h2_total_pion = (TH2*)f1->Get("h2_total");

  TCanvas *c = new TCanvas("c", "Canvas", 1200, 600);
  gStyle->SetOptStat(0);
  c->Divide(2,1);
  TLegend *leg = new TLegend(0.1,0.5,0.4,0.9);

  const Int_t icr = 1;  // icr = 0,1,2,3
  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};
  const char *pname[3] = {"PbScW", "PbScE", "PbGlE"};

  for(Int_t part=0; part<3; part++)
  {
    TH1 *h_ph_asym = h2_ph_asym->ProjectionX("h_ph_asym", secl[part], sech[part]);
    TH1 *h_ph_emc_asym_noprob = h2_ph_emc_asym_noprob->ProjectionX("h_ph_emc_asym_noprob", secl[part], sech[part]);
    TH1 *h_pi0_asym = h2_pi0_asym->ProjectionX("h_pi0_asym", secl[part], sech[part]);
    TH1 *h_pi0_emc_asym_noprob = h2_pi0_emc_asym_noprob->ProjectionX("h_pi0_emc_asym_noprob", secl[part], sech[part]);

    TH1 *h_merge_photon = h2_merge_photon->ProjectionX("h_merge_photon", 13+icr+4*part, 13+icr+4*part);
    TH1 *h_total_photon = h2_merge_photon->ProjectionX("h_total_photon", 37+icr+4*part, 37+icr+4*part);
    TH1 *h_merge_pion = h2_merge_photon->ProjectionX("h_merge_pion", 13+icr+4*part, 13+icr+4*part);
    TH1 *h_total_pion = h2_merge_photon->ProjectionX("h_total_pion", 37+icr+4*part, 37+icr+4*part);
    //TH1 *h_merge_photon = h2_merge_photon->ProjectionX("h_merge_photon", secl[part], sech[part]);
    //TH1 *h_total_photon = h2_total_photon->ProjectionX("h_total_photon", secl[part], sech[part]);
    //TH1 *h_merge_pion = h2_merge_photon->ProjectionX("h_merge_pion", secl[part], sech[part]);
    //TH1 *h_total_pion = h2_total_photon->ProjectionX("h_total_pion", secl[part], sech[part]);

    c->cd(1);

    TGraphAsymmErrors *gr_ph = new TGraphAsymmErrors(h_ph_emc_asym_noprob, h_ph_asym, "n");
    Glist->Add(gr_ph);
    gr_ph->SetName(Form("gr_%d",part));
    gr_ph->SetTitle("Merging rate for photons");
    gr_ph->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    gr_ph->GetYaxis()->SetTitle("rate");
    gr_ph->GetXaxis()->SetRangeUser(0., 30.);
    gr_ph->GetYaxis()->SetRangeUser(0., 1.);
    gr_ph->SetMarkerStyle(20+part);
    gr_ph->SetMarkerColor(1+part);
    if(part==0)
      gr_ph->Draw("AP");
    else
      gr_ph->Draw("P");

    TGraphAsymmErrors *gr_ph_1 = new TGraphAsymmErrors(h_merge_photon, h_total_photon, "n");
    gr_ph_1->SetMarkerStyle(24+part);
    gr_ph_1->SetMarkerColor(1+part);
    gr_ph_1->Draw("P");

    leg->AddEntry(gr_ph, Form("reconstructed p_{T}", "%s", pname[part]), "P");
    leg->AddEntry(gr_ph_1, Form("truth p_{T}", "%s", pname[part]), "P");
    leg->Draw();

    c->cd(2);

    TGraphAsymmErrors *gr_pi0 = new TGraphAsymmErrors(h_pi0_emc_asym_noprob, h_pi0_asym, "n");
    Glist->Add(gr_pi0);
    gr_pi0->SetName(Form("gr_%d",3+part));
    gr_pi0->SetTitle("Merging rate for #pi^{0}");
    gr_pi0->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    gr_pi0->GetYaxis()->SetTitle("rate");
    gr_pi0->GetXaxis()->SetRangeUser(0., 30.);
    gr_pi0->GetYaxis()->SetRangeUser(0., 1.);
    gr_pi0->SetMarkerStyle(20+part);
    gr_pi0->SetMarkerColor(1+part);
    if(part==0)
      gr_pi0->Draw("AP");
    else
      gr_pi0->Draw("P");

    TGraphAsymmErrors *gr_pion_1 = new TGraphAsymmErrors(h_merge_pion, h_total_pion, "n");
    gr_pion_1->SetMarkerStyle(24+part);
    gr_pion_1->SetMarkerColor(1+part);
    gr_pion_1->Draw("P");

    leg->Draw();
  }

  c->Print("Merge.pdf");

  TFile *fout = new TFile("Merge.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
