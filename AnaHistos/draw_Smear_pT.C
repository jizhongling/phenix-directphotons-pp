void draw_Smear_pT()
{
  TFile *f[2];
  THnSparse *hn_pion[2];

  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};
  const char *pname[3] = {"PbScW", "PbScE", "PbGlE"};  // must be non-const if used directly in TLegend

  f[0] = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-warn-histo.root");
  f[1] = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/pi0cross_run13pp510gev/fastMC/eff-histo.root");

  mc();
  legi();
  for(Int_t i=0; i<2; i++)
  {
    hn_pion[i] = (THnSparse*)f[i]->Get("hn_pion");
    //mc(i, 6,5);
    for(Int_t part=0; part<1; part++)
    {
      hn_pion[i]->GetAxis(3)->SetRange(secl[part],sech[part]);
      Int_t ipad = 1;
      for(Int_t ipt=24; ipt<26; ipt+=2)
      {
        hn_pion[i]->GetAxis(0)->SetRange(ipt+1,ipt+2);
        TH1 *h_pt1 = hn_pion[i]->Projection(1);
        h_pt1->Scale( 1. / h_pt1->GetEntries() );

        //mcd(i, ipad++);
        Double_t pTlow = h_pt1->GetXaxis()->GetBinLowEdge(ipt+1);
        Double_t pTup = h_pt1->GetXaxis()->GetBinUpEdge(ipt+2);
        h_pt1->SetTitle(Form("pT: %3.1f-%3.1f GeV",pTlow,pTup));
        //aset(h_pt1);
        h_pt1->SetLineColor(i+1);
        //if(part==0)
        if(i==0)
        {
          h_pt1->Draw("HIST");
          leg0->AddEntry(h_pt1, "Mine");
        }
        else
        {
          h_pt1->Draw("HIST SAME");
          leg0->AddEntry(h_pt1, "Sasha's");
        }

        //delete h_pt1;
      }
    }
  }

  leg0->Draw();
  c0->Print("plots/Smear-pT.pdf");
  //c0->Print("plots/Smear-pT-Mine.pdf");
  //c1->Print("plots/Smear-pT-Sasha.pdf");
}
