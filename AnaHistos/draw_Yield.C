void GenerateGraph(TFile *f, TObjArray *Glist, Int_t trig, Int_t part)
{
  cout << "\npart " << part << endl;
  Double_t gx[30], gy[30], egy[30];
  TGraphErrors *gr;

  Double_t BR[30], MissR[30], Acceptance[30], Smear[30], TrigE[30], Merge[30];
  Double_t eBR[30], eMissR[30], eAcceptance[30], eSmear[30], eTrigE[30], eMerge[30];

  ReadGraphErrors("Acceptance.root", 3+part, gx, Acceptance, eAcceptance);
  ReadGraphAsymmErrors("TriggerEfficiency.root", 9+3*trig+part/2, gx, TrigE, eTrigE);

  const Double_t PI = 3.1415927;

  const Double_t sec_low[3] = {1, 5, 7};
  const Double_t sec_high[3] = {4, 6, 8};

  //THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");

  //hn_1photon->GetAxis(0)->SetRange(sec_low[part], sec_high[part]);
  //TH1 *h_pt = (TH1*)hn_1photon->Projection(1)->Clone(Form("h_pt_s%d",part));

  //gr = new TGraphErrors(h_pt);
  //Glist->AddAtAndExpand(gr,12*trig+8+part);
  //gr->SetName(Form("gr_%d",12*trig+8+part));
  //gr->SetTitle("Photon yield");

  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_pion");
  //hn_2photon->GetAxis(5)->SetRange(3,3);

  TAxis *axis_sec = hn_2photon->GetAxis(0);
  TAxis *axis_pt = hn_2photon->GetAxis(1);
  TAxis *axis_minv = hn_2photon->GetAxis(2);

  axis_sec->SetRange(sec_low[part], sec_high[part]);

  Int_t bin047 = axis_minv->FindBin(0.047);
  Int_t bin067 = axis_minv->FindBin(0.067);
  Int_t bin087 = axis_minv->FindBin(0.087);
  Int_t bin097 = axis_minv->FindBin(0.097);
  Int_t bin110 = axis_minv->FindBin(0.110);
  Int_t bin112 = axis_minv->FindBin(0.112);
  Int_t bin135 = axis_minv->FindBin(0.135);
  Int_t bin160 = axis_minv->FindBin(0.160);
  Int_t bin162 = axis_minv->FindBin(0.162);
  Int_t bin177 = axis_minv->FindBin(0.177);
  Int_t bin187 = axis_minv->FindBin(0.187);
  Int_t bin212 = axis_minv->FindBin(0.212);
  Int_t bin227 = axis_minv->FindBin(0.227);

  TCanvas *c = new TCanvas("c", "Canvas", 2400, 2000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  c->Divide(6,5);

  Int_t ipad = 1;
  for(Int_t ipt=2; ipt<21; ipt++)
  {
    c->cd(ipad++);

    char buf[100];
    Double_t low = axis_pt->GetBinLowEdge(ipt+1);
    Double_t high = axis_pt->GetBinLowEdge(ipt+2);
    sprintf(buf, "p_{T} %3.1f-%3.1f", low, high);

    axis_pt->SetRange(ipt+1,ipt+1);
    TH1 *h_inv_mass = hn_2photon->Projection(2);
    h_inv_mass->SetTitle(buf);

    Double_t npi0 = 0.;
    for(Int_t ib=bin110; ib<bin160; ib++)
      npi0 += h_inv_mass->GetBinContent(ib);
    gy[ipt] = npi0;
    egy[ipt] = sqrt(npi0);

    //TF1 *fn1 = new TF1("fn1", "gaus", 0., 0.5);
    //TF1 *fn2 = new TF1("fn2", "gaus(0)+pol2(3)", 0., 0.5);
    //TF1 *fn3 = new TF1("fn3", "pol2", 0., 0.5);

    //Double_t par[10];
    //h_inv_mass->Fit(fn1, "Q0", "", 0.112, 0.162);
    //fn2->SetParameters( fn1->GetParameters() );
    //h_inv_mass->Fit(fn2, "Q0", "", 0.047, 0.227);
    //TFitResultPtr fit = h_inv_mass->Fit(fn2, "QES0", "", 0.047, 0.227);
    //bool fitok = !fit->IsEmpty();
    //Double_t *covmat = 0;
    //if(fitok)
    //{
    //  TMatrixDSym submat(3);
    //  fit->GetCovarianceMatrix()->GetSub(3, 5, submat);
    //  covmat = submat.GetMatrixArray(); 
    //}
    //fn2->GetParameters(par);
    //fn3->SetParameters(par[3], par[4], par[5]);
    //fn2->SetLineColor(kRed);
    //fn3->SetLineColor(kGreen);
    //h_inv_mass->DrawCopy();
    //fn2->Draw("SAME");
    //fn3->Draw("SAME");

    //Double_t nbgside = 0.;
    //for(Int_t ib=bin047; ib<bin097; ib++)
    //  nbgside += h_inv_mass->GetBinContent(ib);
    //for(Int_t ib=bin177; ib<bin227; ib++)
    //  nbgside += h_inv_mass->GetBinContent(ib);
    //nbgside /= 2.;

    //Double_t npair = 0.;
    //for(Int_t ib=bin112; ib<bin162; ib++)
    //  npair += h_inv_mass->GetBinContent(ib);

    //Double_t npion;
    //Double_t enpion;

    //if(fitok)
    //{
    //  Double_t nbgfit = fn3->Integral(0.112, 0.162) / 0.001;
    //  Double_t enbgfit = fn3->IntegralError(0.112, 0.162, fn3->GetParameters(), covmat) / 0.001;
    //  Double_t nfit = npair - nbgfit;
    //  Double_t enfit = sqrt( npair + enbgfit*enbgfit );
    //}

    //Double_t nsub = npair - nbgside;
    //Double_t ensub = sqrt(npair + nbgside);

    //if(fitok && nfit > 0. && nfit < npair)
    //{
    //  npion = nfit;
    //  enpion = enfit;
    //}
    //else
    //{
    //  npion = nsub;
    //  enpion = ensub;
    //}

    //gy[ipt] = npion;
    //egy[ipt] = enpion;

    delete h_inv_mass;
  }

  gr = new TGraphErrors(30, gx, gy, 0, egy);
  Glist->AddAtAndExpand(gr,12*trig+8+part);
  gr->SetName(Form("gr_%d",12*trig+8+part));
  gr->SetTitle("#pi^{0} yield");

  c->Print(Form("Yield-ert%c-part%d.pdf",97+trig,part));
  delete c;

  return;
}

void draw_Yield()
{
  gROOT->ProcessLine(".L ReadGraph.C");
  gROOT->ProcessLine(".L Chi2Fit.C");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/PhotonNode-histo.root");
  //TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/AnalysisTrain/pat/macro/DirectPhotonPP-num.root");
  TObjArray *Glist = new TObjArray();
  const Int_t trig = 2;
  TGraphErrors *gr;
  TLegend *leg;

  for(Int_t part=0; part<3; part++)
    GenerateGraph(f, Glist, trig, part);

  TCanvas *c0 = new TCanvas("c0", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  const char *legname[3] = {"PbScW", "PbScE", "PbGlE"};
  for(Int_t part=0; part<3; part++)
  {
    gr = (TGraphErrors*)Glist->At(12*trig+8+part);
    gr->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr->GetXaxis()->SetRangeUser(0., 20.);
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
  c0->Print(Form("Yield-ert%c.pdf",97+trig));

  TFile *fout = new TFile(Form("Yield-ert%c.root",97+trig), "RECREATE");
  Glist->Write();
  fout->Close();
}
