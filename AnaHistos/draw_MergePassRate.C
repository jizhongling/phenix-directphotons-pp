void draw_MergePassRate()
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/MissingRatio-histo.root");
  THnSparse *hn_merge = (THnSparse*)f->Get("hn_merge");

  mc(0, 2,4);

  for(Int_t ieta=0; ieta<8; ieta++)
    for(Int_t part=0; part<2; part++)
    {
      mcd(0, ieta+1);
      hn_merge->GetAxis(3)->SetRange(ieta+1-ieta/7*7,ieta+1-ieta/7*3);  // |eta|
      hn_merge->GetAxis(1)->SetRange(secl[part],sech[part]);  // sector

      hn_merge->GetAxis(2)->SetRange(1,2);  // passed = 0 or 1
      TH1 *h_total = hn_merge->Projection(0);

      hn_merge->GetAxis(2)->SetRange(2,2);  // passed = 1
      TH1 *h_passed = hn_merge->Projection(0);

      TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_passed, h_total);
      gr->SetTitle( Form("|#eta|: %.2f - %.2f",(ieta-ieta/7*7)*0.05,(ieta+1-ieta/7*3)*0.05) );
      aset(gr, "p_{T} [GeV]","Bad Pass", 10.,30., 0.,1.);
      style(gr, part+20, part+1);
      if(part==0)
        gr->Draw("AP");
      else
        gr->Draw("P");

      delete h_total;
      delete h_passed;
    }

  c0->Print("plots/MergePassRate.pdf");
}
