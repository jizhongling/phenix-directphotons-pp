void GenerateGraph(TFile *f, TObjArray *Glist, Int_t part)
{
  Double_t gx[30], gy[30] = {}, egy[30] = {};
  Double_t Acceptance[30], eAcceptance[30];
  ReadGraphErrors("Acceptance.root", 3+part, gx, Acceptance, eAcceptance);

  const Double_t PI = TMath::Pi();
  const Int_t sec_low[3] = {1, 5, 7};
  const Int_t sec_high[3] = {4, 6, 8};

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  hn_pion->GetAxis(3)->SetRange(2,2);

  TAxis *axis_sec = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);

  Int_t bin047 = axis_minv->FindBin(0.047);
  Int_t bin097 = axis_minv->FindBin(0.097);
  Int_t bin112 = axis_minv->FindBin(0.112);
  Int_t bin162 = axis_minv->FindBin(0.162);
  Int_t bin177 = axis_minv->FindBin(0.177);
  Int_t bin227 = axis_minv->FindBin(0.227);

  for(Int_t ipt=2; ipt<25; ipt++)
  {
    char title[100];
    Double_t low = axis_pt->GetBinLowEdge(ipt+1);
    Double_t high = axis_pt->GetBinLowEdge(ipt+2);
    sprintf(title, "p_{T} %3.1f-%3.1f", low, high);

    axis_sec->SetRange(sec_low[part], sec_high[part]);
    axis_pt->SetRange(ipt+1,ipt+1);
    TH1 *h_minv = hn_pion->Projection(2);
    h_minv->SetTitle(title);

    Double_t npi0 = h_minv->Integral(bin112,bin162) - ( h_minv->Integral(bin047,bin097) + h_minv->Integral(bin177,bin227) ) / 2.;
    gy[ipt] = npi0;
    egy[ipt] = sqrt(npi0);

    delete h_minv;
  }

  TGraphErrors *gr = new TGraphErrors(30, gx, gy, 0, egy);
  Glist->AddAtAndExpand(gr,part);
  gr->SetName(Form("gr_%d",part));
  gr->SetTitle("#pi^{0} yield");

  return;
}

void draw_Yield()
{
  gROOT->ProcessLine(".L ReadGraph.C");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");
  TObjArray *Glist = new TObjArray();

  for(Int_t part=0; part<3; part++)
    GenerateGraph(f, Glist, part);

  mc();
  mcd();
  legi(0, 0.4,0.7,0.7,0.9);
  const char *legname[3] = {"PbScW", "PbScE", "PbGlE"};
  for(Int_t part=0; part<3; part++)
  {
    TGraphErrors *gr = (TGraphErrors*)Glist->At(part);
    aset(gr, "p_{T} [GeV]","Yield", 0.,20., 1.,1e5);
    gPad->SetLogy();
    style(gr, part+20, part+1);
    leg0->AddEntry(gr, legname[part], "P");
    if(part == 0)
      gr->Draw("AP");
    else
      gr->Draw("P");
  }
  leg0->Draw();
  c0->Print("Yield.pdf");

  TFile *fout = new TFile("Yield.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
