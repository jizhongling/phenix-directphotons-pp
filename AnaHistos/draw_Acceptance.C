TGraphErrors *DivideHisto(TH1 *h1, TH1 *h2)
{
  Int_t gn = h1->GetXaxis()->GetNbins();
  Double_t *gx = new Double_t[gn];
  Double_t *gy = new Double_t[gn];
  Double_t *egx = new Double_t[gn];
  Double_t *egy = new Double_t[gn];

  for(Int_t i=0; i<gn; i++)
  {
    gx[i] = gy[i] = egx[i] = egy[i] = 0.;

    Double_t h1x = h1->GetXaxis()->GetBinCenter(i+1);
    Double_t h1y = h1->GetBinContent(i+1); 
    Double_t h2y = h2->GetBinContent(i+1); 
    Double_t eh1x = h1->GetXaxis()->GetBinWidth(i+1) / 2.;
    Double_t eh1y = h1->GetBinError(i+1);
    Double_t eh2y = h2->GetBinError(i+1);

    if( h1y > 0. && h2y > 0. )
    {
      gx[i] = h1x;
      gy[i] = h1y / h2y;
      egx[i] = eh1x;
      egy[i] = gy[i] * sqrt( pow(eh1y/h1y,2.) + pow(eh2y/h2y,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(gn, gx, gy, egx, egy);
  return graph;
}

TGraphErrors *DivideHisto(TH1 *h1, TH1 *h2, Double_t kh1, Double_t kh2)
{
  Int_t gn = h1->GetXaxis()->GetNbins();
  Double_t *gx = new Double_t[gn];
  Double_t *gy = new Double_t[gn];
  Double_t *egx = new Double_t[gn];
  Double_t *egy = new Double_t[gn];

  for(Int_t i=0; i<gn; i++)
  {
    gx[i] = gy[i] = egx[i] = egy[i] = 0.;

    Double_t h1x = h1->GetXaxis()->GetBinCenter(i+1);
    Double_t h1y = h1->GetBinContent(i+1) * kh1; 
    Double_t h2y = h2->GetBinContent(i+1) * kh2; 
    Double_t eh1x = h1->GetXaxis()->GetBinWidth(i+1) / 2.;
    Double_t eh1y = h1->GetBinError(i+1) * kh1;
    Double_t eh2y = h2->GetBinError(i+1) * kh2;

    if( h1y > 0. && h2y > 0. )
    {
      gx[i] = h1x;
      gy[i] = h1y / h2y;
      egx[i] = eh1x;
      egy[i] = gy[i] * sqrt( pow(eh1y/h1y,2.) + pow(eh2y/h2y,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(gn, gx, gy, egx, egy);
  return graph;
}

void GenerateAcceptance(TFile *fsig, TFile *ftot, TObjArray *Glist, Int_t ispion)
{
  const Int_t npT = 31;
  const Double_t vpT[npT] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

  TF1 *cross = new TF1("cross", "x*(1/(1+exp((x-[5])/[6]))*[0]/pow(1+x/[1],[2])+(1-1/(1+exp((x-[5])/[6])))*[3]/pow(x,[4]))", 0, 30);
  cross->SetParameters(2.02819e+04, 4.59173e-01, 7.51170e+00, 1.52867e+01, 7.22708e+00, 2.15396e+01, 3.65471e+00);

  TH1 *h_sig[3];
  for(Int_t part=0; part<3; part++)
    h_sig[part] = new TH1D(Form("h_sig%d",part), "#pi^{0} signal count; p_{T} [GeV/c];", npT-1,vpT);
  TH1 *h_tot;

  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};
  char *pname[3] = {"PbScW", "PbScE", "PbGlE"};  // must be non-const if used directly in TLegend

  if(ispion == 0)
  {
    THnSparse *hn_sig = (THnSparse*)fsig->Get("hn_photon");
    for(Int_t part=0; part<3; part++)
    {
      hn_sig->GetAxis(2)->SetRange(secl[part],sech[part]);
      h_sig[part] = (TH1*)hn_sig->Projection(0)->Clone();
    }
    h_tot = (TH1*)ftot->Get("h_photon");
  }
  else if(ispion == 1)
  {
    THnSparse *hn_pion = (THnSparse*)fsig->Get("hn_separate");
    for(Int_t part=0; part<3; part++)
    {
      hn_pion->GetAxis(2)->SetRange(secl[part],sech[part]);
      TH1 *h_pt = hn_pion->Projection(1);
      for(Int_t ipt=0; ipt<npT-1; ipt++)
      {
        //hn_pion->GetAxis(1)->SetRange(ipt+1, ipt+1);
        //TH1 *h_minv = hn_pion->Projection(2);
        //Double_t npion = h_minv->Integral(112,162) - ( h_minv->Integral(47,97) + h_minv->Integral(177,227) ) / 2.;
        Double_t npion = h_pt->Integral(ipt+1,ipt+1);
        h_sig[part]->SetBinContent(ipt+1, npion);
        h_sig[part]->SetBinError( ipt+1, sqrt(npion) * cross->Eval(vpT[ipt]) );
        //delete h_minv;
      }
      delete h_pt;
    }
    h_tot = (TH1*)ftot->Get("h_pion");
  }

  mc();
  mcd();
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);

  for(Int_t part=0; part<3; part++)
  {
    TGraphErrors *gr = DivideHisto(h_sig[part], h_tot);
    //TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_sig[part], h_tot, "n");
    Glist->Add(gr);
    gr->SetName(Form("gr_%d",3*ispion+part));
    if(ispion == 0)
      gr->SetTitle("Photon acceptance");
    else if(ispion == 1)
      gr->SetTitle("#pi^{0} acceptance");
    aset(gr, "p_{T} [GeV]", "acceptance", 0.,30., 0.,0.12);
    style(gr, 20+part, 1+part);
    if(part==0)
      gr->Draw("AP");
    else
      gr->Draw("P");
    leg0->AddEntry(gr, pname[part], "P");
    leg0->Draw();
  }

  if(ispion == 0)
  {
    c0->Print("Acceptance-photon.pdf");
    delete c0;
  }
  else if(ispion == 1)
  {
    c0->Print("Acceptance-pion.pdf");
    delete c0;
  }

  for(Int_t part=0; part<3; part++)
    delete h_sig[part];
  delete h_tot;

  return;
}

void draw_Acceptance()
{
  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-warn-histo.root");
  TObjArray *Glist = new TObjArray();

  for(Int_t ispion=1; ispion<2; ispion++)
  {
    cout << "\nispion " << ispion << endl;
    GenerateAcceptance(f, f, Glist, ispion);
  }

  TFile *fout = new TFile("Acceptance.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
