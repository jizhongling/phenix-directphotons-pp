void draw_YieldPhoton()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  TH1 *h_1photon[2];  // 0: inclusive; 1: isolated

  int evtype = 2;
  int bbc10cm = 1;
  int checkmap = 1;
  int ival = 1;

  TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
  h_1photon_t = (TH1*)h_1photon_t->Clone();
  h_1photon_t->Reset();

  for(int isolated=0; isolated<2; isolated++)
  {
    h_1photon[isolated] = (TH1*)h_1photon_t->Clone(Form("h_1photon_isolated%d",isolated));
    for(int part=0; part<3; part++)
      for(int evenodd=0; evenodd<2; evenodd++)
        for(int pattern=0; pattern<3; pattern++)
        {
          int ih = part + 3*evenodd + 3*2*pattern + 3*2*3*evtype + 3*2*3*4*bbc10cm + 3*2*3*4*2*checkmap + 3*2*3*4*2*2*isolated + 3*2*3*4*2*2*2*ival;
          TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
          h_1photon[isolated]->Add(h_tmp);
          delete h_tmp;
        }
    if(isolated == 1)
      h_1photon[0]->Add(h_1photon[1]);
  }

  mc();
  mcd();
  legi(0, 0.4,0.7,0.9,0.9);
  gPad->SetLogy();
  h_1photon[0]->SetTitle("All photon yield");
  aset(h_1photon[0], "p_{T} [GeV]","Yield", 6.1,30.);
  style(h_1photon[0], 20, 1);
  style(h_1photon[1], 20, 2);
  h_1photon[0]->Draw("HISTO");
  h_1photon[1]->Draw("HISTO SAME");
  leg0->AddEntry(h_1photon[0], "Inclusive all photon", "L");
  leg0->AddEntry(h_1photon[1], "Isolated all photon", "L");
  leg0->Draw();
  c0->Print("plots/YieldPhoton.pdf");
}
