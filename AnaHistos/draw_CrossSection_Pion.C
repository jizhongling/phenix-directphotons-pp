void GenerateGraph(TFile *f, TObjArray *Glist, Int_t ispion, Int_t trig, Int_t part)
{
  cout << "\npart " << part << endl;
  Double_t gx[30], gy[6][30], egy[6][30];
  TGraphErrors *gr;

  const Double_t PI = 3.1415927;

  if(trig == 0)
    const Double_t BBCCount = 2.107514 * pow(10,12);  // 2107514662585
  else if(trig == 1)
    const Double_t BBCCount = 2.107514 * pow(10,12);  // 2107514662585
  else if(trig == 2)
    const Double_t BBCCount =  6.817931 * pow(10,11);  // 681793156166

  const Double_t BBCRatio = 0.67;
  const Double_t eBBCRatio = 0.16;
  const Double_t BBCCross = 32.5 * pow(10,9);
  const Double_t eBBCCross = 3.25 * pow(10,9);
  const Double_t TrigBBC = 0.91;
  const Double_t eTrigBBC = 0.01;
  Double_t BR[30], MissR[30], Acceptance[30], Smear[30], TrigE[30], Merge[30];
  Double_t eBR[30], eMissR[30], eAcceptance[30], eSmear[30], eTrigE[30], eMerge[30];

  ReadGraphErrors("Acceptance.root", 3*ispion+part, gx, Acceptance, eAcceptance);
  ReadGraphAsymmErrors("TriggerEfficiency.root", 9*ispion+3*trig+part/2, gx, TrigE, eTrigE);

  const Double_t sec_low[3] = {1, 5, 7};
  const Double_t sec_high[3] = {4, 6, 8};

  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  hn_1photon->GetAxis(0)->SetRange(sec_low[part], sec_high[part]);
  TH1 *h_1photon = hn_1photon->Projection(1);

  TH2 *h2_sig_extra = (TH2*)f->Get("h2_sig_extra");
  TH1 *h_sig_extra = h2_sig_extra->ProjectionY("h_sig_extra", sec_low[part], sec_high[part]);

  TH2 *h2_bg_extra = (TH2*)f->Get("h2_bg_extra");
  TH1 *h_bg_extra = h2_bg_extra->ProjectionY("h_bg_extra", sec_low[part], sec_high[part]);

  //THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");
  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_pion");

  TAxis *axis_sec = hn_2photon->GetAxis(0);
  if(ispion == 0)
    TAxis *axis_pt = hn_2photon->GetAxis(1);
  else if(ispion == 1)
    TAxis *axis_pt = hn_2photon->GetAxis(1);
  TAxis *axis_minv = hn_2photon->GetAxis(2);

  axis_sec->SetRange(sec_low[part], sec_high[part]);

  Int_t bin047 = axis_minv->FindBin(0.047);
  Int_t bin067 = axis_minv->FindBin(0.067);
  Int_t bin087 = axis_minv->FindBin(0.087);
  Int_t bin097 = axis_minv->FindBin(0.097);
  Int_t bin112 = axis_minv->FindBin(0.112);
  Int_t bin162 = axis_minv->FindBin(0.162);
  Int_t bin177 = axis_minv->FindBin(0.177);
  Int_t bin187 = axis_minv->FindBin(0.187);
  Int_t bin212 = axis_minv->FindBin(0.212);
  Int_t bin227 = axis_minv->FindBin(0.227);

  const Int_t nData = 256;
  vector<Double_t> x(nData), y(nData), sigma_y(nData);

  TCanvas *c = new TCanvas("c", "Canvas", 2400, 2000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  c->Divide(6,5);

  Int_t ipad = 1;
  for(Int_t ipt=0; ipt<30; ipt++)
  {
    c->cd(ipad++);

    char buf[100];
    Double_t low = axis_pt->GetBinLowEdge(ipt+1);
    Double_t high = axis_pt->GetBinLowEdge(ipt+2);
    sprintf(buf, "p_{T} %3.1f-%3.1f", low, high);

    axis_pt->SetRange(ipt+1,ipt+1);
    TH1 *h_inv_mass = hn_2photon->Projection(2);
    h_inv_mass->SetTitle(buf);

    Double_t nphoton = h_1photon->GetBinContent(ipt+1);
    Double_t nsig_extra = h_sig_extra->GetBinContent(ipt+1);
    Double_t nbg_extra = h_bg_extra->GetBinContent(ipt+1);
    Double_t nextra = nsig_extra - nbg_extra/2.;

    x.clear();
    y.clear();
    sigma_y.clear();

    for(Int_t ib=bin067; ib<bin087; ib++)
    {
      Double_t xx = axis_minv->GetBinCenter(ib);
      Double_t yy = h_inv_mass->GetBinContent(ib);
      Double_t sigma_yy = h_inv_mass->GetBinError(ib);
      x.push_back(xx);
      y.push_back(yy);
      sigma_y.push_back(sigma_yy);
    }

    for(Int_t ib=bin187; ib<bin212; ib++)
    {
      Double_t xx = axis_minv->GetBinCenter(ib);
      Double_t yy = h_inv_mass->GetBinContent(ib);
      Double_t sigma_yy = h_inv_mass->GetBinError(ib);
      x.push_back(xx);
      y.push_back(yy);
      sigma_y.push_back(sigma_yy);
    }

    Double_t nbggpr, enbggpr;
    BgGPR(x, y, sigma_y, nbggpr, enbggpr);
    nbggpr /= 0.001;
    enbggpr /= 0.001;

    TF1 *fn1 = new TF1("fn1", "gaus", 0., 0.5);
    TF1 *fn2 = new TF1("fn2", "gaus(0)+pol2(3)", 0., 0.5);
    TF1 *fn3 = new TF1("fn3", "pol2", 0., 0.5);

    Double_t par[10];
    h_inv_mass->Fit(fn1, "Q0", "", 0.112, 0.162);
    fn2->SetParameters( fn1->GetParameters() );
    h_inv_mass->Fit(fn2, "Q0", "", 0.047, 0.227);
    TFitResultPtr fit = h_inv_mass->Fit(fn2, "QES", "", 0.047, 0.227);
    bool fitok = !fit->IsEmpty();
    Double_t *covmat = 0;
    if(fitok)
    {
      TMatrixDSym submat(3);
      fit->GetCovarianceMatrix()->GetSub(3, 5, submat);
      covmat = submat.GetMatrixArray(); 
    }
    fn2->GetParameters(par);
    fn3->SetParameters(par[3], par[4], par[5]);
    fn3->SetLineColor(kGreen);
    fn3->Draw("SAME");

    Double_t nbgside = 0.;
    for(Int_t ib=bin047; ib<bin097; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    for(Int_t ib=bin177; ib<bin227; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    nbgside /= 2.;

    Double_t npair = 0.;
    for(Int_t ib=bin112; ib<bin162; ib++)
      npair += h_inv_mass->GetBinContent(ib);

    Double_t npion[2];
    Double_t enpion[2];

    if(fitok)
    {
      Double_t nbgfit = fn3->Integral(0.112, 0.162) / 0.001;
      Double_t enbgfit = fn3->IntegralError(0.112, 0.162, fn3->GetParameters(), covmat) / 0.001;
      Double_t nfit = npair - nbgfit;
      Double_t enfit = sqrt( npair + enbgfit*enbgfit );
      Double_t ndiff = fabs(nbgfit - nbgside - nextra);
    }

    Double_t nsub = npair - nbgside;
    Double_t ensub = sqrt(npair + nbgside);
    Double_t ngpr = npair - nbggpr;
    Double_t engpr = sqrt( npair + enbggpr*enbggpr );

    if(fitok && nfit > 0. && nfit < npair)
    {
      npion[0] = nfit;
      enpion[0] = enfit;
    }
    else
    {
      npion[0] = nsub;
      enpion[0] = ensub;
    }

    npion[1] = ngpr;
    enpion[1] = engpr;

    for(Int_t im=0; im<2; im++)
    {
      gy[3*im+2][ipt] = (BBCCross/BBCCount/BBCRatio) * npion[im] / (2*PI*gx[ipt]) / ((high-low)*0.8)
        / Acceptance[ipt] / TrigBBC /  TrigE[ipt];
      egy[3*im+2][ipt] = gy[3*im+2][ipt] * sqrt(
          pow(enpion[im]/npion[im],2.) + pow(eAcceptance[ipt]/Acceptance[ipt],2.) + pow(eTrigE[ipt]/TrigE[ipt],2.)
          // + pow(eBBCCross/BBCCross,2.) + pow(eBBCRatio/BBCRatio,2.) + pow(eTrigBBC/TrigBBC,2.)
          );
    }
  }

  for(Int_t im=0; im<2; im++)
  {
    gr = new TGraphErrors(30, gx, gy[3*im+2], 0, egy[3*im+2]);
    Glist->AddAtAndExpand(gr,36*im+12*trig+8+part);
    gr->SetName(Form("gr_%d",36*im+12*trig+8+part));
    gr->SetTitle("#pi^{0} Cross Section");
  }

  c->Print(Form("CrossSection-pion%d-ert%c-part%d.pdf",ispion,97+trig,part));
  delete c;

  return;
}

void draw_CrossSection_Pion()
{
  gSystem->Load("libGausProc.so");
  gROOT->ProcessLine(".L ReadGraph.C");
  gROOT->ProcessLine(".L BgGPR.C");
  gROOT->ProcessLine(".L Chi2Fit.C");

  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ertc-cv/total.root");
  TObjArray *Glist = new TObjArray();
  const Int_t trig = 2;
  TGraphErrors *gr;
  TLegend *leg;

  for(Int_t part=0; part<3; part++)
    for(Int_t ispion=1; ispion<2; ispion++)
      GenerateGraph(f, Glist, ispion, trig, part);

  for(Int_t im=0; im<2; im++)
  {
    for(Int_t i=2; i<3; i++)
    {
      Double_t gx[30], gy[30], egy[30];
      Double_t gyy[3][30], egyy[3][30];

      for(Int_t part=0; part<3; part++)
        ReadGraph((TGraphErrors*)Glist->At(36*im+12*trig+4*i+part), gx, gyy[part], egyy[part]);

      for(Int_t ipt=0; ipt<30; ipt++)
      {
        Double_t yy[3];
        Double_t eyy[3];
        for(Int_t part=0; part<3; part++)
        {
          yy[part] = gyy[part][ipt];
          eyy[part] = egyy[part][ipt];
        }
        Chi2Fit(3, yy, eyy, gy[ipt], egy[ipt]);
      }

      gr = new TGraphErrors(30, gx, gy, 0, egy);
      Glist->AddAtAndExpand(gr, 36*im+12*trig+4*i+3);
      gr->SetName(Form("gr_%d",36*im+12*trig+4*i+3));
      gr->SetTitle(((TGraphErrors*)Glist->At(36*im+12*trig+4*i))->GetTitle());
    }
  }

  TCanvas *c0 = new TCanvas("c0", "Canvas", 1800, 1200);
  gStyle->SetOptStat(0);
  c0->Divide(3,2);

  const char *legname[4] = {"PbScW", "PbScE", "PbGlE", "Combined"};
  for(Int_t im=0; im<2; im++)
  {
    for(Int_t i=2; i<3; i++)
      for(Int_t part=0; part<4; part++)
      {
        c0->cd(i+3*(part/3)+1);
        gr = (TGraphErrors*)Glist->At(36*im+12*trig+4*i+part);
        gr->GetXaxis()->SetTitle("p_{T} [GeV]");
        gr->GetXaxis()->SetRangeUser(0., 30.);
        if(i==0)
        {
          gr->GetYaxis()->SetTitle("Ratio");
          gr->GetYaxis()->SetRangeUser(0., 1.2);
        }
        else
        {
          gPad->SetLogy();
          gr->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3} [pb*GeV^{-2}*c^{-3}]");
          gr->GetYaxis()->SetRangeUser(0.5, 1e5);
        }
        gr->GetYaxis()->SetTitleOffset(1.2);
        gr->SetMarkerColor(part+1);
        gr->SetMarkerStyle(part+20);
        gr->SetMarkerSize(1.);
        if(part%3 == 0)
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
    c0->Print(Form("CrossSection-method-%d-ert%c.pdf",im,97+trig));
    c0->Clear("D");
  }

  TFile *fout = new TFile(Form("CrossSection-ert%c.root",97+trig), "RECREATE");
  Glist->Write();
  fout->Close();
}
