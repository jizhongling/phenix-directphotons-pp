void draw_SelfVeto()
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  TFile *f_out = new TFile("data/eta-selfveto.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-warn-histo.root");
  TH3 *h3_isopair = (TH3*)f->Get("h3_isoeta");

  mc();
  mcd();

  for(int part=0; part<3; part++)
  {
    TH1 *h_total = h3_isopair->ProjectionX("h_total", secl[part],sech[part], 1,2);
    TH1 *h_passed = h3_isopair->ProjectionX("h_passed", secl[part],sech[part], 2,2);

    TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_passed, h_total);
    gr->SetName(Form("gr_%d",part));
    f_out->cd();
    gr->Write();

    aset(gr, "p_{T} [GeV]","#frac{isoboth}{isopair}", 0.,30., 0.,1.);
    style(gr, 20+part, 1+part);
    if(part==0)
      gr->Draw("AP");
    else
      gr->Draw("P");

    delete h_total;
    delete h_passed;
  }

  c0->Print("plots/SelfVeto.pdf");
  f_out->Close();
}
