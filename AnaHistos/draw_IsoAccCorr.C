void draw_IsoAccCorr()
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  TGraphAsymmErrors *gr[3];

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo.root");
  THnSparse *hn_isoprompt = (THnSparse*)f->Get("hn_isoprompt");

  mc();
  mcd();

  for(int part=0; part<3; part++)
  {
    hn_isoprompt->GetAxis(1)->SetRange(secl[part],sech[part]);
    hn_isoprompt->GetAxis(2)->SetRange(1,1);
    TH1 *h_isoall = hn_isoprompt->Projection(0);
    hn_isoprompt->GetAxis(2)->SetRange(2,2);
    TH1 *h_isoacc = hn_isoprompt->Projection(0);

    gr[part] = new TGraphAsymmErrors(h_isoall, h_isoacc);
    gr[part]->SetNameTitle(Form("gr_%d",part), "IsoAll/IsoAcc");
    aset(gr[part], "p_{T} [GeV/c]","IsoAll/IsoAcc", 5.,30., 0.8,1.);
    style(gr[part], 20+part, 1+part);
    if(part==0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");

    delete h_isoall;
    delete h_isoacc;
  }
  legi(0, 0.2,0.6,0.5,0.9);
  leg0->AddEntry(gr[0], "PbSc west", "P");
  leg0->AddEntry(gr[1], "PbSc east", "P");
  leg0->AddEntry(gr[2], "PbGl", "P");
  leg0->Draw();

  c0->Print("plots/IsoAll2IsoAcc.pdf");
  TFile *f_out = new TFile("data/IsoAll2IsoAcc.root", "RECREATE");
  for(int part=0; part<3; part++)
    gr[part]->Write();
  f_out->Close();
}
