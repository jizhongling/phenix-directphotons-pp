void draw_ConversionRate()
{
  const char *cname[3] = {"vtxconv", "eeinconv", "eeoutconv"};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/MissingRatio-histo.root");

  TH2 *h2_total = (TH2*)f->Get("h2_photon");
  TH2 *h2_passed[3]; for(int ic=0; ic<3; ic++)
    h2_passed[ic] = (TH2*)f->Get( Form("h2_%s",cname[ic]) );

  mc(0, 3,1);

  TGraphAsymmErrors *gr[3][2];  // gr[ic][arm]
  for(int ic=0; ic<3; ic++)
    for(int arm=0; arm<2; arm++)
    {
      TH1 *h_total = h2_total->ProjectionX("h_total", arm+1, arm+1);
      TH1 *h_passed = h2_passed[ic]->ProjectionX("h_passed", arm+1, arm+1);
      gr[ic][arm] = new TGraphAsymmErrors(h_passed, h_total);
      delete h_total;
      delete h_passed;

      mcd(0, ic+1);
      aset(gr[ic][arm], "p_{T} [GeV]","rate", 0.,30., 0.,0.2);
      style(gr[ic][arm], arm+20, arm+1);
      if(arm == 0)
        gr[ic][arm]->Draw("AP");
      else
        gr[ic][arm]->Draw("P");

      if(ic == 0)
      {
        gr[ic][arm]->Fit("pol0", "Q","", 2.,28.);
        if(arm == 0)
        {
          gPad->Update();
          TPaveStats *st = (TPaveStats*)gr[ic][arm]->FindObject("stats");
          st->SetX1NDC(0.2);
          st->SetX2NDC(0.6);
        } // arm == 0
      } // ic == 0
    }

  mcd(0, 1);
  TLine *line = new TLine();
  line->DrawLine(0.,0.099,30.,0.099);
  line->DrawLine(0.,0.147,30.,0.147);

  gr[0][0]->SetTitle("All conversion inside magnetic field");
  gr[1][0]->SetTitle("e^{+}e^{-} conversion inside magnetic field");
  gr[2][0]->SetTitle("e^{+}e^{-} conversion outside magnetic field");

  c0->Print("plots/ConversionRate.pdf");
}
