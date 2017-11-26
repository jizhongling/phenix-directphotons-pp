void draw_MergeAsym()
{
  char name[100];

  TFile *f_sim = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/MissingRatio-histo.root");
  TFile *f_data = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  THnSparse *hn_asym_sim = (THnSparse*)f_sim->Get("hn_asym");
  THnSparse *hn_asym_data = (THnSparse*)f_data->Get("hn_asym");

  hn_asym_sim->GetAxis(4)->SetRange(113,162);
  hn_asym_data->GetAxis(4)->SetRange(113,162);

  mc(0, 2,5);
  legi(0, 0.4,0.7,0.7,0.9);

  const char *sec_name[2] = {"PbSc", "PbGl"};
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  TGraphAsymmErrors *gr_sim[2][5];
  TGraphAsymmErrors *gr_data[2][5];

  TLine *line = new TLine(5., 0., 5., 1.);
  line->SetLineColor(kBlue);

  for(Int_t part=0; part<2; part++)
  {
    hn_asym_sim->GetAxis(0)->SetRange(secl[part], sech[part]);
    hn_asym_data->GetAxis(0)->SetRange(secl[part], sech[part]);

    hn_asym_sim->GetAxis(2)->SetRange(1,50);
    hn_asym_data->GetAxis(2)->SetRange(1,50);

    TH1 *h_pt_sim = (TH1*)hn_asym_sim->Projection(1)->Clone("h_pt_sim");
    TH1 *h_pt_data = (TH1*)hn_asym_data->Projection(1)->Clone("h_pt_data");

    for(Int_t ias=0; ias<5; ias++)
    {
      mcd(0, 2*ias+part+1);

      hn_asym_sim->GetAxis(2)->SetRange(10*ias+1,10*ias+10);
      hn_asym_data->GetAxis(2)->SetRange(10*ias+1,10*ias+10);

      Double_t AsymLow = hn_asym_sim->GetAxis(2)->GetBinLowEdge(10*ias+1);
      Double_t AsymUp = hn_asym_sim->GetAxis(2)->GetBinUpEdge(10*ias+10);
      sprintf(name, "%s Asym: %.1f-%.1f", sec_name[part], AsymLow, AsymUp);

      TH1 *h_1pt_sim = (TH1*)hn_asym_sim->Projection(1)->Clone("h_1pt_sim");
      TH1 *h_1pt_data = (TH1*)hn_asym_data->Projection(1)->Clone("h_1pt_data");

      gr_sim[part][ias] = new TGraphAsymmErrors(h_1pt_sim, h_pt_sim, "n");
      gr_data[part][ias] = new TGraphAsymmErrors(h_1pt_data, h_pt_data, "n");

      gr_sim[part][ias]->SetTitle(name);
      gr_data[part][ias]->SetTitle(name);

      aset(gr_sim[part][ias], "p_{T} [GeV]", "ratio", 0.,30., 0.,1.);
      aset(gr_data[part][ias], "p_{T} [GeV]", "ratio", 0.,30., 0.,1.);

      style(gr_sim[part][ias], 20, 1);
      style(gr_data[part][ias], 21, 2);

      gr_sim[part][ias]->Draw("AP");
      gr_data[part][ias]->Draw("P");
      line->Draw();

      if(ias==0 && part==0)
      {
        leg0->AddEntry(gr_sim[part][ias], "PISA", "LEP");
        leg0->AddEntry(gr_data[part][ias], "Data", "LEP");
      }
      leg0->Draw();

      delete h_1pt_sim;
      delete h_1pt_data;
    }

    delete h_pt_sim;
    delete h_pt_data;
  }

  c0->Print("plots/MergeAsym.pdf");
}
