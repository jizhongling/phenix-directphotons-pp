void draw_MergeAngle()
{
  TFile *f = new TFile("/phenix/plhf/zji/data/pisaRun13/AnaPHPythia-histo.root");
  TFile *f1 = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/MissingRatio-histo.root");

  const int reco = 1; // truth pT:0; reconstructed pT: 1

  TH1::SetDefaultSumw2();

  THnSparse *hn_pi0_total = (THnSparse*)f->Get("hn_total");
  THnSparse *hn_pi0_separate = (THnSparse*)f->Get("hn_separate");
  THnSparse *hn_total = (THnSparse*)f1->Get("hn_total");
  THnSparse *hn_separate = (THnSparse*)f1->Get("hn_separate");

  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  for(int part=0; part<2; part++)
  {
    hn_pi0_total->GetAxis(2)->SetRange(secl[part], sech[part]);
    hn_pi0_separate->GetAxis(2)->SetRange(secl[part], sech[part]);
    hn_total->GetAxis(2)->SetRange(secl[part], sech[part]);
    hn_separate->GetAxis(2)->SetRange(secl[part], sech[part]);

    mc(part, 5,5);
    int ipad = 1;
    for(int ipt=10; ipt<30; ipt++)
    {
      hn_pi0_total->GetAxis(reco)->SetRange(ipt+1,ipt+1);
      hn_pi0_separate->GetAxis(reco)->SetRange(ipt+1,ipt+1);
      hn_total->GetAxis(reco)->SetRange(ipt+1,ipt+1);
      hn_separate->GetAxis(reco)->SetRange(ipt+1,ipt+1);

      char name[100];
      double low =  hn_total->GetAxis(reco)->GetBinLowEdge(ipt+1);
      double high = hn_total->GetAxis(reco)->GetBinLowEdge(ipt+2);
      sprintf(name, "pT: %3.1f-%3.1f GeV", low, high);

      TH1 *h_total = hn_total->Projection(3);
      mcd(part,ipad);
      aset(h_total);
      h_total->SetTitle(name);
      h_total->GetXaxis()->SetRangeUser(0.,0.1);
      h_total->SetLineColor(kBlue);
      h_total->DrawCopy("E");
      delete h_total;

      TH1 *h_pi0_total = hn_pi0_total->Projection(3);
      mcd(part,ipad);
      aset(h_pi0_total);
      h_pi0_total->SetTitle(name);
      h_pi0_total->SetLineColor(kRed);
      h_pi0_total->DrawCopy("ESAME");
      delete h_pi0_total;

      TH1 *h_pi0_separate = hn_pi0_separate->Projection(3);
      delete h_pi0_separate;

      TH1 *h_separate = hn_separate->Projection(3);
      delete h_separate;
      
      ipad++;
    }
    gROOT->ProcessLine(Form("c%d->Print(\"Angle-part%d.pdf\");",part,part));
  }
}
