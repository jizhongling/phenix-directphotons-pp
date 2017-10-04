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

TGraphErrors *DivideGraph(TGraphErrors *gr1, TGraphErrors *gr2)
{
  Int_t gn = gr1->GetN();
  Double_t *gx = new Double_t[gn];
  Double_t *gy = new Double_t[gn];
  Double_t *egx = new Double_t[gn];
  Double_t *egy = new Double_t[gn];

  for(Int_t i=0; i<gn; i++)
  {
    gx[i] = gy[i] = egx[i] = egy[i] = 0.;

    Double_t g1x, g1y, g2x, g2y;
    gr1->GetPoint(i, g1x, g1y);
    gr2->GetPoint(i, g2x, g2y);
    Double_t eg1x = gr1->GetErrorX(i);
    Double_t eg1y = gr1->GetErrorY(i);
    Double_t eg2y = gr2->GetErrorY(i);

    if( g1y > 0. && g2y > 0. )
    {
      gx[i] = g1x;
      gy[i] = g1y / g2y;
      egx[i] = eg1x;
      egy[i] = gy[i] * sqrt( pow(eg1y/g1y,2.) + pow(eg2y/g2y,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(gn, gx, gy, egx, egy);
  return graph;
}

void draw_Smear_Sasha()
{
  TFile *f1 = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-warn-histo.root");
  THnSparse *hn_pion = (THnSparse*)f1->Get("hn_pion");

  TFile *f2 = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/pi0cross_run13pp510gev/fastMC/eff-histo.root");
  THnSparse *hn_pion2 = (THnSparse*)f2->Get("hn_pion");

  TFile *f[3];
  TH1 *hpt_acc1[3];
  TH1 *hpt_acc2[3];

  mc(0);
  mc(1);

  legi(0, 0.2,0.8,0.9,0.9);
  legi(1, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  leg1->SetNColumns(3);

  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};
  const char *pname[3] = {"PbScW", "PbScE", "PbGlE"};  // must be non-const if used directly in TLegend

  for(Int_t part=0; part<3; part++)
  {
    hn_pion->GetAxis(3)->SetRange(secl[part],sech[part]);
    TH1 *h_pt0 = hn_pion->Projection(0);
    TH1 *h_pt1 = hn_pion->Projection(1);
    TGraphErrors *gr1 = DivideHisto(h_pt1, h_pt0);
    delete h_pt0;
    delete h_pt1;

    //hn_pion2->GetAxis(3)->SetRange(secl[part],sech[part]);
    //TH1 *h_pt0 = hn_pion2->Projection(0);
    //TH1 *h_pt1 = hn_pion2->Projection(1);
    //TGraphErrors *gr2 = DivideHisto(h_pt1, h_pt0);
    //delete h_pt0;
    //delete h_pt1;

    f[part] = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/pi0cross_run13pp510gev/fastMC/eff-%d.root",part));
    hpt_acc1[part] = (TH1*)f[part]->Get("hpt_acc1");
    hpt_acc2[part] = (TH1*)f[part]->Get("hpt_acc2");
    TGraphErrors *gr2 = DivideHisto(hpt_acc2[part], hpt_acc1[part]);
    //TGraphAsymmErrors *gr2 = new TGraphAsymmErrors(hpt_acc2[part], hpt_acc1[part], "n");

    mcd(0);
    TGraphErrors *gr = DivideGraph(gr1, gr2);
    gr->SetTitle("Smear Ratio");
    aset(gr, "p_{T} [GeV]", "#frac{My smear}{Sasha's smear}", 0.,30., 0.7,1.3);
    style(gr, 20+part, 1+part);
    if(part==0)
      gr->Draw("AP");
    else
      gr->Draw("P");
    leg0->AddEntry(gr, pname[part], "P");
    leg0->Draw();

    mcd(1);
    gr1->SetTitle("Smear");
    aset(gr1, "p_{T} [GeV]", "Smear", 0.,30., 0.7,1.5);
    style(gr1, 20+part, 1+part);
    style(gr2, 20+part, 1+part);
    gr2->SetLineWidth(6.);
    if(part==0)
    {
      gr1->Draw("AP");
      gr2->Draw("L");
    }
    else
    {
      gr1->Draw("P");
      gr2->Draw("L");
    }
    leg1->AddEntry(gr1, pname[part], "P");
    leg1->Draw();

    gr2->RemovePoint(0);
    for(Int_t ip=1; ip<3; ip++)
    {
      gr1->RemovePoint(ip);
      gr2->RemovePoint(ip);
      gr->RemovePoint(ip);
    }
  }

  c0->Print("SmearRatio.pdf");
  c1->Print("Smear.pdf");
}
