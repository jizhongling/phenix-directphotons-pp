void draw_MergeAsym()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/MissingRatio-histo.root");
  //TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonNode-histo.root");
  THnSparse *hn_asym = (THnSparse*)f->Get("hn_asym");
  hn_asym->GetAxis(4)->SetRange(113,162);

  mc(0, 2,1);
  mc(1, 2,1);

  const char *name[2] = {"PbSc", "PbGl"};
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  TLine *line[2];
  line[0] = new TLine(0., 5.5, 30., 5.5);
  line[1] = new TLine(0., 4., 30., 4.);

  for(Int_t part=0; part<2; part++)
  {
    line[part]->SetLineColor(kRed);
    hn_asym->GetAxis(0)->SetRange(secl[part], sech[part]);

    mcd(0, part+1);
    gPad->SetLogz();
    TH2 *h2_dR_pT = hn_asym->Projection(3,1);
    h2_dR_pT->SetTitle(name[part]);
    aset(h2_dR_pT, "p_{T}^{#pi^{0}} [GeV]");
    h2_dR_pT->GetZaxis()->SetRangeUser(1e-7,1e1);
    h2_dR_pT->DrawCopy("colz");
    line[part]->Draw();
    delete h2_dR_pT;

    mcd(1, part+1);
    gPad->SetLogz();
    TH2 *h2_asym_pT = hn_asym->Projection(2,1);
    h2_asym_pT->SetTitle(name[part]);
    aset(h2_asym_pT, "p_{T}^{#pi^{0}} [GeV]");
    h2_asym_pT->GetZaxis()->SetRangeUser(1e-7,1e1);
    h2_asym_pT->DrawCopy("colz");
    line[part]->Draw();
    delete h2_asym_pT;
  }

  c0->Print("MergeAsym-dR-sim.pdf");
  c1->Print("MergeAsym-asym-sim.pdf");
}
