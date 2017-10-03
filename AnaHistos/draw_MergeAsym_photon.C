void draw_MergeAsym_photon()
{
  char name[100];

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/MissingRatio-histo.root");
  THnSparse *hn_photon = (THnSparse*)f->Get("hn_photon");

  mc(0, 2,3);

  const char *sec_name[2] = {"PbSc", "PbGl"};
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  const Int_t pTCut[3] = {11, 23, 27};

  for(Int_t part=0; part<2; part++)
  {
    hn_photon->GetAxis(2)->SetRange(secl[part], sech[part]);
    for(Int_t ipt=0; ipt<3; ipt++)
    {
      mcd(0, 2*ipt+part+1);

      hn_photon->GetAxis(0)->SetRange(pTCut[ipt], pTCut[ipt]);
      Double_t pTLow = hn_photon->GetAxis(0)->GetBinLowEdge(pTCut[ipt]);
      Double_t pTUp = hn_photon->GetAxis(0)->GetBinUpEdge(pTCut[ipt]);
      sprintf(name, "%s p_{T}^{#pi^{0}}: %.1f-%.1f [GeV]", sec_name[part], pTLow, pTUp);

      TH1 *h_pt = hn_photon->Projection(1);
      for(Int_t i=20; i<30; i++)
      {
        h_pt->SetBinContent(i+1, h_pt->GetBinContent(i+1)/4.);
        h_pt->SetBinError(i+1, h_pt->GetBinError(i+1)/4.);
      }

      h_pt->SetTitle(name);
      aset(h_pt, "Photon p_{T} [GeV]");
      h_pt->DrawCopy("HIST E");
      
      delete h_pt;
    }
  }

  c0->Print("MergeAsym-photon.pdf");
}
