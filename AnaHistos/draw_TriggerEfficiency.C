void GenerateTriggerEfficiency(TFile *f, TObjArray *Glist, Int_t ispion, Int_t trig)
{
  TH3 *h3_trig;
  if(ispion == 0)
    h3_trig = (TH3*)f->Get("h3_trig");
  else if(ispion == 1)
    h3_trig = (TH3*)f->Get("h3_trig_pion");

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  TGraphAsymmErrors *gr[3];
  Int_t sec_low[3] = {1, 5, 7};
  Int_t sec_high[3] = {4, 6, 8};
  for(Int_t part=0; part<3; part++)
  {
    TH1 *h_trig = h3_trig->ProjectionX("h_trig", sec_low[part],sec_high[part], 5+4*part,5+4*part);
    TH1 *h_total = h3_trig->ProjectionX("h_total", sec_low[part],sec_high[part], 1, 1);
    gr[part] = new TGraphAsymmErrors(h_trig, h_total);
    Glist->Add(gr[part]);
    gr[part]->SetName(Form("gr_%d",9*ispion+3*trig+part));
    if(ispion == 0)
      gr[part]->SetTitle(Form("ERT_4x4%c trigger efficiency for photon",97+trig));
    else if(ispion == 1)
      gr[part]->SetTitle(Form("ERT_4x4%c trigger efficeincy for #pi^{0}",97+trig));
    gr[part]->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr[part]->GetYaxis()->SetTitle("efficiency");
    gr[part]->GetYaxis()->SetTitleOffset(1.2);
    gr[part]->GetXaxis()->SetRangeUser(0., 30.);
    gr[part]->GetYaxis()->SetRangeUser(0., 1.1);
    gr[part]->SetMarkerColor(1+part);
    gr[part]->SetMarkerStyle(20+part);
    gr[part]->SetMarkerSize(1.);
    if(part==0)
      gr[part]->Draw("APE");
    else
      gr[part]->Draw("PE");
  }

  TLegend *leg = new TLegend(0.6, 0.1, 0.9, 0.3);
  leg->AddEntry(gr[0], "PbScW", "LPE");
  leg->AddEntry(gr[1], "PbScE", "LPE");
  leg->AddEntry(gr[2], "PbGlE", "LPE");
  leg->Draw();

  c->Print(Form("TriggerEfficiency-pion%d-ert%c.pdf",ispion,97+trig));
  delete c;

  return;
}

void draw_TriggerEfficiency()
{
  //TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510MinBias/9161/data/total.root");
  TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510ERT/9473/data/total.root");
  TObjArray *Glist = new TObjArray();

  for(Int_t ispion=0; ispion<2; ispion++)
    for(Int_t trig=0; trig<3; trig++)
    {
      cout << "\nispion " << ispion << endl;
      GenerateTriggerEfficiency(f, Glist, ispion, trig);
    }

  TFile *fout = new TFile("TriggerEfficiency.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
