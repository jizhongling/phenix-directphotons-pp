void GenerateERTEff(TFile *f, TObjArray *Glist, Int_t ispion, Int_t trig)
{
  TH3 *h3_trig;
  if(ispion == 0)
    h3_trig = (TH3*)f->Get("h3_ert");
  else if(ispion == 1)
    h3_trig = (TH3*)f->Get("h3_ert_pion");

  mc();
  mcd();

  TGraphAsymmErrors *gr[2], *grc[2];
  Int_t sec_low[2] = {1, 7};
  Int_t sec_high[2] = {6, 8};
  for(Int_t part=0; part<2; part++)
  {
    TH1 *h_trig = h3_trig->ProjectionX("h_trig", sec_low[part],sec_high[part], 2+trig,2+trig);
    TH1 *h_total = h3_trig->ProjectionX("h_total", sec_low[part],sec_high[part], 1,1);
    gr[part] = new TGraphAsymmErrors(h_trig, h_total);
    if(ispion == 0)
      gr[part]->SetTitle(Form("ERT_4x4%c trigger efficiency for photon",97+trig));
    else if(ispion == 1)
      gr[part]->SetTitle(Form("ERT_4x4%c trigger efficeincy for #pi^{0}",97+trig));
    aset(gr[part], "p_{T} [GeV]", "efficiency", 0.,30., 0.,1.1);
    style(gr[part], 20+part, 1+part);
    if(part==0)
      gr[part]->Draw("APE");
    else
      gr[part]->Draw("PE");

    //gr[part]->SetName(Form("gr_%d",9*ispion+3*trig+part));
    //Glist->Add(gr[part]);
    grc[part] = (TGraphAsymmErrors*)gr[part]->Clone(Form("gr_%d",9*ispion+3*trig+part));
    Glist->Add(grc[part]);

    gr[part]->Fit("pol0", "QR", "", 12., 100.);
    TF1 *f_pol0 = (TF1*)gr[part]->GetListOfFunctions()->FindObject("pol0");
    f_pol0->SetLineColor(3+part);
    Double_t EffBar = f_pol0->GetParameter(0);
    Double_t eEffBar = f_pol0->GetParError(0);
    for(Int_t ipt=21; ipt<31; ipt++)
    {
      Double_t xx, yy;
      gr[part]->GetPoint(ipt, xx, yy);
      grc[part]->SetPoint(ipt, xx, EffBar);
      Double_t exxl = gr[part]->GetErrorXlow(ipt);
      Double_t exxh = gr[part]->GetErrorXhigh(ipt);
      grc[part]->SetPointError(ipt, exxl, exxh, eEffBar/2., eEffBar/2.);
    }
    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr[part]->FindObject("stats");
    st->SetX1NDC(0.5+0.2*part);
    st->SetY1NDC(0.4);
    st->SetX2NDC(0.7+0.2*part);
    st->SetY2NDC(0.6);
  }

  if( ispion == 1 && trig == 2 )
  {
    TGraph *gr_Sasha_PbSc = new TGraph("Sasha_PbSc_trig.txt");
    gr_Sasha_PbSc->Draw("C");
    TGraph *gr_Sasha_PbGl = new TGraph("Sasha_PbGl_trig.txt");
    gr_Sasha_PbGl->Draw("C");
  }

  legi(0, 0.7, 0.2, 0.9, 0.4);
  leg0->AddEntry(gr[0], "PbSc", "LPE");
  leg0->AddEntry(gr[1], "PbGl", "LPE");
  leg0->Draw();

  c0->Print(Form("TriggerEfficiency-pion%d-ert%c.pdf",ispion,97+trig));
  delete c0;

  return;
}

void draw_ERTEff()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");
  TObjArray *Glist = new TObjArray();

  for(Int_t ispion=0; ispion<2; ispion++)
    for(Int_t trig=0; trig<3; trig++)
    {
      cout << "\nispion " << ispion << endl;
      GenerateERTEff(f, Glist, ispion, trig);
    }

  TFile *fout = new TFile("ERTEff.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
