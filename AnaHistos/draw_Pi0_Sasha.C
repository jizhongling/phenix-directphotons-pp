void draw_Pi0_Sasha()
{
  gROOT->ProcessLine(".L ReadGraph.C");

  TGraphErrors *gr;
  TObjArray *Glist = new TObjArray();
  const Int_t trig = 2;

  char hname[100];
  Double_t gx[30], gy[3][30] = {}, egy[3][30] = {};
  Double_t Acceptance[30], eAcceptance[30];
  ReadGraphErrors("Acceptance.root", 3, gx, Acceptance, eAcceptance);


  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/AnalysisTrain/pat/macro/Pi0PP-ERT.root");
  TH1 *mchist[3][40];  // mchist[is][ip]
  for(Int_t is=0; is<3; is++)
    for(Int_t ip=0; ip<40; ip++)
    {
      sprintf(hname,"mc_s%d_bcc0_pt_%03d_tp",is,5*ip);
      mchist[is][ip] = (TH1*)f->Get(hname);
    }

  Int_t bin110 = mchist[0][0]->GetXaxis()->FindBin(0.110);
  Int_t bin160 = mchist[0][0]->GetXaxis()->FindBin(0.160);

  TH1 *h_minv[3][25];  // h_minv[is][ipt]
  for(Int_t is=0; is<3; is++)
    for(Int_t ipt=0; ipt<25; ipt++)
    {
      sprintf(hname,"h_minv_sec%d_pt%d",is,ipt);
      h_minv[is][ipt] = new TH1F(hname,hname,1000,0.,1.);
    }

  Int_t ptl[25];
  Int_t ptr[25];
  Int_t ipt = 0;
  for(Int_t ip=0; ip<20; ip++)
  {
    ptl[ipt] = ip;
    ptr[ipt] = ip;
    ipt++;
  }
  for(Int_t ip=20; ip<40; ip+=4)
  {
    ptl[ipt] = ip;
    ptr[ipt] = ip+3;
    ipt++;
  }

  for(Int_t is=0; is<3; is++)
    for(Int_t ipt=2; ipt<25; ipt++)
      for(Int_t ip=ptl[ipt]; ip<=ptr[ipt]; ip++)
        h_minv[is][ipt]->Add(mchist[is][ip]);

  for(Int_t is=0; is<3; is++)
    for(Int_t ipt=0; ipt<25; ipt++)
    {
      Double_t npi0 = 0.;
      for(Int_t ib=bin110; ib<bin160; ib++)
        npi0 += h_minv[is][ipt]->GetBinContent(ib);
      gy[is][ipt] = npi0;
      egy[is][ipt] = sqrt(npi0);
    }

  for(Int_t part=0; part<3; part++)
  {
    gr = new TGraphErrors(30, gx, gy[part], 0, egy[part]);
    Glist->AddAtAndExpand(gr,12*trig+8+part);
    gr->SetName(Form("gr_%d",12*trig+8+part));
    gr->SetTitle("#pi^{0} yield");
  }

  TCanvas *c0 = new TCanvas("c0", "Canvas", 600, 600);
  gStyle->SetOptStat(0);
  TLegend *leg;

  const char *legname[3] = {"PbScW", "PbScE", "PbGlE"};
  for(Int_t part=0; part<3; part++)
  {
    gr = (TGraphErrors*)Glist->At(12*trig+8+part);
    gr->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr->GetXaxis()->SetRangeUser(0., 30.);
    gPad->SetLogy();
    gr->GetYaxis()->SetTitle("Yield");
    gr->GetYaxis()->SetRangeUser(1., 1e5);
    gr->GetYaxis()->SetTitleOffset(1.2);
    gr->SetMarkerColor(part+1);
    gr->SetMarkerStyle(part+20);
    gr->SetMarkerSize(1.);
    if(part == 0)
    {
      gr->Draw("AP");
      leg = new TLegend(0.4, 0.7, 0.7, 0.9);
      leg->Draw();
    }
    else
    {
      gr->Draw("P");
    }
    leg->AddEntry(gr, legname[part], "P");
  }
  c0->Print(Form("Yield-ert%c-sasha.pdf",97+trig));

  TFile *fout = new TFile(Form("Yield-ert%c-sasha.root",97+trig), "RECREATE");
  Glist->Write();
  fout->Close();
}
