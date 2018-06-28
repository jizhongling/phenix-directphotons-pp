void draw_ClusterEShare()
{
  char name[100];
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/MissingRatio-histo.root");
  THnSparse *hn_sharing = (THnSparse*)f->Get("hn_sharing");

  mc(0, 6,5);
  mc(1, 6,5);
  
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};
  const double tsize[2] = {5.5, 4.};

  TLine *line = new TLine();
  line->SetLineColor(kRed);

  for(int part=0; part<2; part++)
  {
    hn_sharing->GetAxis(0)->SetRange(secl[part],sech[part]);
    for(int ipt=0; ipt<30; ipt++)
    {
      hn_sharing->GetAxis(1)->SetRange(ipt+1,ipt+1);
      double pTlow = hn_sharing->GetAxis(1)->GetBinLowEdge(ipt+1);
      double pTup = hn_sharing->GetAxis(1)->GetBinUpEdge(ipt+1);
      sprintf(name, "p_{T}: %.1f-%.1f", pTlow, pTup);

      mcd(0, ipt+1);
      TH1 *h_ratio = hn_sharing->Projection(2);
      h_ratio->SetTitle(name);
      aset(h_ratio);
      h_ratio->DrawCopy("HIST E");
      delete h_ratio;

      mcd(1, ipt+1);
      TH2 *h_dR = hn_sharing->Projection(3,4);
      h_dR->SetTitle(name);
      aset(h_dR);
      h_dR->DrawCopy("COLZ");
      line->DrawLine(-tsize[part], -tsize[part], -tsize[part], +tsize[part]);
      line->DrawLine(-tsize[part], +tsize[part], +tsize[part], +tsize[part]);
      line->DrawLine(+tsize[part], +tsize[part], +tsize[part], -tsize[part]);
      line->DrawLine(+tsize[part], -tsize[part], -tsize[part], -tsize[part]);
      delete h_dR;
    }

    c0->Print(Form("plots/ClusterESharing-part%d.pdf",part));
    c1->Print(Form("plots/ClusterDistance-part%d.pdf",part));
    c0->Clear("D");
    c1->Clear("D");
  }

}
