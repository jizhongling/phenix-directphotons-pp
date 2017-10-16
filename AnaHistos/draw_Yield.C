void GenerateGraph(TFile *f, TObjArray *Glist, Int_t trig, Int_t part)
{
  cout << "\npart " << part << endl;
  Double_t gx[30], gy[30], egy[30];

  Double_t BR[30], MissR[30], Acceptance[30], Smear[30], TrigE[30], Merge[30];
  Double_t eBR[30], eMissR[30], eAcceptance[30], eSmear[30], eTrigE[30], eMerge[30];

  ReadGraphErrors("Acceptance.root", 3+part, gx, Acceptance, eAcceptance);

  const Double_t PI = TMath::Pi();
  const Int_t sec_low[3] = {1, 5, 7};
  const Int_t sec_high[3] = {4, 6, 8};

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");

  TAxis *axis_sec = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);

  Int_t bin110 = axis_minv->FindBin(0.110);
  Int_t bin160 = axis_minv->FindBin(0.160);

  for(Int_t ipt=2; ipt<21; ipt++)
  {
    char title[100];
    Double_t low = axis_pt->GetBinLowEdge(ipt+1);
    Double_t high = axis_pt->GetBinLowEdge(ipt+2);
    sprintf(title, "p_{T} %3.1f-%3.1f", low, high);

    axis_sec->SetRange(sec_low[part], sec_high[part]);
    axis_pt->SetRange(ipt+1,ipt+1);
    TH1 *h_minv = hn_pion->Projection(2);
    h_minv->SetTitle(title);

    Double_t npi0 = 0.;
    for(Int_t ib=bin110; ib<bin160; ib++)
      npi0 += h_minv->GetBinContent(ib);
    gy[ipt] = npi0;
    egy[ipt] = sqrt(npi0);

    delete h_minv;
  }

  TGraphErrors *gr = new TGraphErrors(30, gx, gy, 0, egy);
  Glist->AddAtAndExpand(gr,12*trig+8+part);
  gr->SetName(Form("gr_%d",12*trig+8+part));
  gr->SetTitle("#pi^{0} yield");

  return;
}

void draw_Yield()
{
  gROOT->ProcessLine(".L ReadGraph.C");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/PhotonNode-histo0.root");
  //TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/AnalysisTrain/pat/macro/DirectPhotonPP-num.root");
  TObjArray *Glist = new TObjArray();
  const Int_t trig = 2;

  for(Int_t part=0; part<3; part++)
    GenerateGraph(f, Glist, trig, part);

  mc();
  mcd();
  legi(0, 0.4,0.7,0.7,0.9);
  const char *legname[3] = {"PbScW", "PbScE", "PbGlE"};
  TGraphErrors *gr[3];
  for(Int_t part=0; part<3; part++)
  {
    gr[part] = (TGraphErrors*)Glist->At(12*trig+8+part);
    aset(gr[part], "p_{T} [GeV]","Yield", 0.,20., 1.,1e5);
    gPad->SetLogy();
    style(gr[part], part+20, part+1);
    leg0->AddEntry(gr[part], legname[part], "P");
    if(part == 0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");
  }
  leg0->Draw();
  c0->Print(Form("Yield-ert%c.pdf",97+trig));

  TFile *fout = new TFile(Form("Yield-ert%c.root",97+trig), "RECREATE");
  Glist->Write();
  fout->Close();
}
