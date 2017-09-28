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
  else if(trig == 3)
    const Double_t BBCCount =  5.118194 * pow(10,8);  // 511819399

  const Double_t BBCRatio = 0.67;
  const Double_t eBBCRatio = 0.16;
  const Double_t BBCCross = 32.5 * pow(10,9);
  const Double_t eBBCCross = 3.25 * pow(10,9);
  const Double_t TrigBBC = 0.91;
  const Double_t eTrigBBC = 0.01;
  Double_t BR[30], MissR[30], Acceptance[30], Conversion[30], Smear[30], TrigE[30], Merge[30];
  Double_t eBR[30], eMissR[30], eAcceptance[30], eConversion[30], eSmear[30], eTrigE[30], eMerge[30];

  ReadGraphErrors("Acceptance.root", 3*ispion+part, gx, Acceptance, eAcceptance);
  ReadGraphErrors("ConversionRate.root", (part+1)/2, gx, Conversion, eConversion);
  ReadGraphAsymmErrors("TriggerEfficiency.root", 9*ispion+3*trig+part/2, gx, TrigE, eTrigE);

  for(Int_t ipt=0; ipt<30; ipt++)
  {
    if(part==0)
    {
      Conversion[ipt] = 0.128;
      eConversion[ipt] = 0.0019;
    }
    else
    {
      Conversion[ipt] = 0.;
      eConversion[ipt] = 0.;
    }
  }

  const Double_t sec_low[3] = {1, 5, 7};
  const Double_t sec_high[3] = {4, 6, 8};

  //THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");
  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_pion");

  TAxis *axis_sec = hn_2photon->GetAxis(0);
  TAxis *axis_pt = hn_2photon->GetAxis(1);
  TAxis *axis_minv = hn_2photon->GetAxis(2);

  axis_sec->SetRange(sec_low[part], sec_high[part]);

  Int_t bin047 = axis_minv->FindBin(0.047);
  Int_t bin097 = axis_minv->FindBin(0.097);
  Int_t bin112 = axis_minv->FindBin(0.112);
  Int_t bin162 = axis_minv->FindBin(0.162);
  Int_t bin177 = axis_minv->FindBin(0.177);
  Int_t bin227 = axis_minv->FindBin(0.227);

  mc(1, 6,5);

  Int_t ipad = 1;
  for(Int_t ipt=1; ipt<30; ipt++)
  {
    mcd(1, ipad++);

    char buf[100];
    Double_t low = axis_pt->GetBinLowEdge(ipt+1);
    Double_t high = axis_pt->GetBinLowEdge(ipt+2);
    sprintf(buf, "p_{T} %3.1f-%3.1f", low, high);

    axis_pt->SetRange(ipt+1,ipt+1);
    TH1 *h_inv_mass = hn_2photon->Projection(2);
    h_inv_mass->SetTitle(buf);

    TF1 *fn1 = new TF1("fn1", "gaus", 0., 0.5);
    TF1 *fn2 = new TF1("fn2", "gaus(0)+pol2(3)", 0., 0.5);
    TF1 *fn3 = new TF1("fn3", "pol2", 0., 0.5);

    Double_t par[10];
    h_inv_mass->Fit(fn1, "Q0", "", 0.112, 0.162);
    fn2->SetParameters( fn1->GetParameters() );
    h_inv_mass->Fit(fn2, "Q0", "", 0.047, 0.227);
    TFitResultPtr fit = h_inv_mass->Fit(fn2, "QES0", "", 0.047, 0.227);
    bool fitok = !fit->IsEmpty();
    //fitok = false;  // !!!
    Double_t *covmat = 0;
    if(fitok)
    {
      TMatrixDSym submat(3);
      fit->GetCovarianceMatrix()->GetSub(3, 5, submat);
      covmat = submat.GetMatrixArray(); 
    }
    fn2->GetParameters(par);
    fn3->SetParameters(par[3], par[4], par[5]);
    fn2->SetLineColor(kRed);
    fn3->SetLineColor(kGreen);
    h_inv_mass->DrawCopy();
    fn2->Draw("SAME");
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
    }

    Double_t nsub = npair - nbgside;
    Double_t ensub = sqrt(npair + nbgside);

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

    for(Int_t im=0; im<1; im++)
    {
      gy[3*im+2][ipt] = (BBCCross/BBCCount/BBCRatio) * npion[im] / (2*PI*gx[ipt]) / ((high-low)*1.0)
        / Acceptance[ipt] / pow(1.-Conversion[ipt],2.) / TrigBBC / TrigE[ipt];
      egy[3*im+2][ipt] = gy[3*im+2][ipt] * sqrt(
          pow(enpion[im]/npion[im],2.) + pow(eAcceptance[ipt]/Acceptance[ipt],2.) + pow(eTrigE[ipt]/TrigE[ipt],2.)
          + pow(2.*eConversion[ipt]/(1.-Conversion[ipt]),2.)
          // + pow(eBBCCross/BBCCross,2.) + pow(eBBCRatio/BBCRatio,2.) + pow(eTrigBBC/TrigBBC,2.)
          );
    }

    delete h_inv_mass;
  }

  for(Int_t im=0; im<1; im++)
  {
    gr = new TGraphErrors(30, gx, gy[3*im+2], 0, egy[3*im+2]);
    Glist->AddAtAndExpand(gr,36*im+12*trig+8+part);
    gr->SetName(Form("gr_%d",36*im+12*trig+8+part));
    gr->SetTitle("#pi^{0} Cross Section");
  }

  c1->Print(Form("CrossSection-pion%d-ert%c-part%d.pdf",ispion,97+trig,part));
  delete c1;

  return;
}

void draw_CrossSection_Pion()
{
  gROOT->ProcessLine(".L ReadGraph.C");
  gROOT->ProcessLine(".L Chi2Fit.C");

  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonNode-histo-ertc.root");
  TObjArray *Glist = new TObjArray();
  const Int_t trig = 2;
  TGraphErrors *gr;
  TLegend *leg;

  for(Int_t part=0; part<3; part++)
    for(Int_t ispion=1; ispion<2; ispion++)
      GenerateGraph(f, Glist, ispion, trig, part);

  for(Int_t im=0; im<1; im++)
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
        Chi2Fit(2, yy, eyy, gy[ipt], egy[ipt]);
      }

      gr = new TGraphErrors(30, gx, gy, 0, egy);
      Glist->AddAtAndExpand(gr, 36*im+12*trig+4*i+3);
      gr->SetName(Form("gr_%d",36*im+12*trig+4*i+3));
      gr->SetTitle(((TGraphErrors*)Glist->At(36*im+12*trig+4*i))->GetTitle());
    }
  }

  mc(0, 3,2);
  legi(0, 0.4, 0.7, 0.7, 0.9);

  const char *legname[4] = {"PbScW", "PbScE", "PbGlE", "Combined"};
  for(Int_t im=0; im<1; im++)
  {
    for(Int_t i=2; i<3; i++)
      for(Int_t part=0; part<3; part++)
      {
        mcd(0, i+3*(part/3)+1);
        gr = (TGraphErrors*)Glist->At(36*im+12*trig+4*i+part);
        gr->GetXaxis()->SetTitle("p_{T} [GeV]");
        gr->GetXaxis()->SetRangeUser(0., 30.);
        if(i==0)
        {
          aset(gr, "p_{T} [GeV]", "Ratio", 0.,30., 0.,1.2);
        }
        else
        {
          gPad->SetLogy();
          aset(gr, "p_{T} [GeV]", "Ed^{3}#sigma/dp^{3} [pb*GeV^{-2}*c^{-3}]", 5.,20., 1.,1e6);
        }
        style(gr, part+20, part+1);
        if(part%3 == 0)
        {
          gr->Draw("AP");
          leg0->Draw();
        }
        else
        {
          gr->Draw("P");
        }
        leg0->AddEntry(gr, legname[part], "P");
      }
    c0->Print(Form("CrossSection-method-%d-ert%c.pdf",im,97+trig));
    c0->Clear("D");
  }

  TFile *fout = new TFile(Form("CrossSection-ert%c.root",97+trig), "RECREATE");
  Glist->Write();
  fout->Close();
}
