void draw_TowerEnergy()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/MissingRatio-histo.root");
  THnSparse *hn_towers = (THnSparse*)f->Get("hn_towers");

  mc(0, 6,5);

  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};

  for(Int_t part=0; part<3; part++)
  {
    hn_towers->GetAxis(0)->SetRange(secl[part],sech[part]);
    Int_t ipad = 1;
    for(Int_t ipt=0; ipt<30; ipt++)
    {
      mcd(0, ipad++);
      hn_towers->GetAxis(1)->SetRange(ipt+1,ipt+1);
      Double_t pTlow = hn_towers->GetAxis(1)->GetBinLowEdge(ipt+1);
      Double_t pThigh = hn_towers->GetAxis(1)->GetBinUpEdge(ipt+1);
      TH1 *h_frac = hn_towers->Projection(2);
      h_frac->SetTitle(Form("pT: %3.1f-%3.1f GeV",pTlow,pThigh));
      aset(h_frac);
      h_frac->DrawCopy();
      delete h_frac;
    }
    c0->Print(Form("plots/TowerEnergy-part%d.pdf",part));
    c0->Clear("D");
  }
}
