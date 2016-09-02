void GenerateGraph(TFile *f, Int_t ispion, Int_t trig, Int_t part, Int_t gn, Double_t **gx, Double_t **gy, Double_t **egy)
{
  const Double_t PI = 3.1415927;
  const Double_t BBCCount = 4.193255 * pow(10,12);  // 4193255288313
  const Double_t BBCRatio = 0.67;
  const Double_t eBBCRatio = 0.16;
  const Double_t BBCCross = 32.5 * pow(10,9);
  const Double_t eBBCCross = 3.25 * pow(10,9);

  Double_t gxx[30];
  Double_t BR[30], MissR[30], Acceptance[30], Smear[30], TrigE[30];
  Double_t eBR[30], eMissR[30], eAcceptance[30], eSmear[30], eTrigE[30];

  ReadGraphErrors("BR.root", 0, gxx, BR, eBR);
  ReadGraphErrors("MissingRatio.root", 24*ispion+15+4*part, gxx, MissR, eMissR);
  ReadGraphErrors("Acceptance.root", 3*ispion+part, gxx, Acceptance, eAcceptance);
  ReadGraphErrors("Smear.root", 24*ispion+15+4*part, gxx, Smear, eSmear);
  ReadGraphAsymmErrors("TriggerEfficiency.root", 9*ispion+3*trig+part, gxx, TrigE, eTrigE);

  Double_t sec_low[3] = {1, 5, 7};
  Double_t sec_high[3] = {4, 6, 8};

  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  hn_1photon->GetAxis(0)->SetRange(sec_low[part], sec_high[part]);
  TH1 *h_1photon = hn_1photon->Projection(1);

  TH2 *h2_sig_extra = (TH2*)f->Get("h2_sig_extra");
  TH1 *h_sig_extra = h2_sig_extra->ProjectionY("h_sig_extra", sec_low[part], sec_high[part]);

  TH2 *h2_bg_extra = (TH2*)f->Get("h2_bg_extra");
  TH1 *h_bg_extra = h2_bg_extra->ProjectionY("h_bg_extra", sec_low[part], sec_high[part]);

  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");

  TAxis *axis0 = hn_2photon->GetAxis(0);
  if(ispion == 0)
    TAxis *axis1_2 = hn_2photon->GetAxis(1);
  else if(ispion == 1)
    TAxis *axis1_2 = hn_2photon->GetAxis(2);
  TAxis *axis3 = hn_2photon->GetAxis(3);

  axis0->SetRange(sec_low[part], sec_high[part]);

  Int_t bin047 = axis3->FindBin(0.047);
  Int_t bin067 = axis3->FindBin(0.067);
  Int_t bin087 = axis3->FindBin(0.087);
  Int_t bin097 = axis3->FindBin(0.097);
  Int_t bin112 = axis3->FindBin(0.112);
  Int_t bin162 = axis3->FindBin(0.162);
  Int_t bin177 = axis3->FindBin(0.177);
  Int_t bin187 = axis3->FindBin(0.187);
  Int_t bin212 = axis3->FindBin(0.212);
  Int_t bin227 = axis3->FindBin(0.227);

  const Int_t nData = 256;
  vector<Double_t> x(nData), y(nData), sigma_y(nData);

  TCanvas *c = new TCanvas("c", "Canvas", 2400, 2000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  c->Divide(6,5);

  Int_t ipad = 1;
  for(Int_t ipt=0; ipt<gn; ipt++)
  {
    c->cd(ipad++);

    char buf[100];
    Double_t low = axis1_2->GetBinLowEdge(ipt+1);
    Double_t high = axis1_2->GetBinLowEdge(ipt+2);
    sprintf(buf, "p_{T} %3.1f-%3.1f", low, high);

    axis1_2->SetRange(ipt+1,ipt+1);
    TH1 *h_inv_mass = hn_2photon->Projection(3);
    h_inv_mass->SetTitle(buf);

    Double_t nphoton = h_1photon->GetBinContent(ipt+1);
    Double_t nsig_extra = h_sig_extra->GetBinContent(ipt+1);
    Double_t nbg_extra = h_bg_extra->GetBinContent(ipt+1);
    Double_t nextra = nsig_extra - nbg_extra/2.;
    Double_t npair = 0.;
    Double_t nbgfit = 0.;
    Double_t nbgside = 0.;
    Double_t nbggpr = 0.;
    Double_t dnbggpr = 0.;

    x.clear();
    y.clear();
    sigma_y.clear();

    for(Int_t ib=bin067; ib<bin087; ib++)
    {
      Double_t xx = axis3->GetBinCenter(ib);
      Double_t yy = h_inv_mass->GetBinContent(ib);
      Double_t sigma_yy = h_inv_mass->GetBinError(ib);
      x.push_back(xx);
      y.push_back(yy);
      sigma_y.push_back(sigma_yy);
    }

    for(Int_t ib=bin187; ib<bin212; ib++)
    {
      Double_t xx = axis3->GetBinCenter(ib);
      Double_t yy = h_inv_mass->GetBinContent(ib);
      Double_t sigma_yy = h_inv_mass->GetBinError(ib);
      x.push_back(xx);
      y.push_back(yy);
      sigma_y.push_back(sigma_yy);
    }

    BgGPR(x, y, sigma_y, nbggpr, dnbggpr);
    nbggpr /= 0.001;
    dnbggpr /= 0.001;

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

    for(Int_t ib=bin047; ib<bin097; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    for(Int_t ib=bin177; ib<bin227; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    nbgside /= 2.;

    for(Int_t ib=bin112; ib<bin162; ib++)
    {
      npair += h_inv_mass->GetBinContent(ib);
      Double_t bincenter = axis3->GetBinCenter(ib);
      nbgfit += fn3->Eval(bincenter);
    }

    if(fitok)
    {
      nbgfit = fn3->Integral(0.112, 0.162) * 1000.;
      Double_t enbgfit = fn3->IntegralError(0.112, 0.162, 0, covmat) * 1000.;
    }

    Double_t nfit = npair - nbgfit;
    Double_t nsub = npair - nbgside;
    Double_t ensub = sqrt( 1./npair + 1./nbgside );
    Double_t ndiff = fabs(nbgfit - nbgside - nextra);
    Double_t ngpr = npair - nbggpr;
    Double_t nerror = fabs(nfit - ngpr);

    if(fitok && nfit > 0. && nfit < npair)
    {
      Double_t enfit = sqrt( 1./npair + pow(enbgfit/nbgfit,2.) );
      Double_t npion = nfit * ( 1. + MissR[ipt] );
      Double_t enpion = npion * sqrt( pow(enfit/nfit,2.) + pow(eMissR[ipt]/(1.+MissR[ipt]),2.) );
    }
    else
    {
      Double_t npion = nsub * ( 1. + MissR[ipt] );
      Double_t enpion = npion * sqrt( pow(ensub/nsub,2.) + pow(eMissR[ipt]/(1.+MissR[ipt]),2.) );
    }
    Double_t ndecay = npion * ( 1. + BR[ipt] );
    Double_t endecay = ndecay * sqrt( pow(enpion/npion,2.) + pow(eBR[ipt]/(1.+BR[ipt]),2.) );
    Double_t ndirect = nphoton - ndecay;
    Double_t endirect = sqrt( 1./nphoton + pow(endecay/ndecay,2.) );

    if(ispion == 0)
    {
      gx[0][ipt] = gxx[ipt];
      gy[0][ipt] = ndirect / nphoton;
      egy[0][ipt] = gy[0][ipt] * sqrt( pow(endirect/ndirect,2.) + 1./nphoton );

      gx[1][ipt] = gxx[ipt];
      gy[1][ipt] = BBCCross * (ndirect/(BBCCount*BBCRatio)) / (2*PI*gx[1][ipt]) / ((high-low)*0.8) / Acceptance[ipt] / Smear[ipt] / TrigE[ipt];
      //egy[1][ipt] = gy[1][ipt] * sqrt( pow(eBBCCross/BBCCross,2.) + pow(endirect/ndirect,2.) + 1./BBCCount + pow(eBBCRatio/BBCRatio,2.)
      //    + pow(eSmear[ipt]/Smear[ipt],2.) + pow(eAcceptance[ipt]/Acceptance[ipt],2.) + pow(eTrigE[ipt]/TrigE[ipt],2.) );
      egy[1][ipt] = gy[1][ipt] * sqrt( pow(endirect/ndirect,2.)
          + pow(eSmear[ipt]/Smear[ipt],2.) + pow(eAcceptance[ipt]/Acceptance[ipt],2.) + pow(eTrigE[ipt]/TrigE[ipt],2.) );

      cout << "pT=" << low << "-" << high << "\tnphoton=" << nphoton
        << "\tnpair=" << npair << "\tnfit=" << nfit << "\tnsub=" << nsub
        << "\tndiff=" << ndiff/nfit*100. << "%\tngpr=" << ngpr
        << "\tnerror=" << nerror/nfit*100. << "%" << endl;
    }
    else if(ispion == 1)
    {
      gx[2][ipt] = gxx[ipt];
      gy[2][ipt] = BBCCross * (npion/2./(BBCCount*BBCRatio)) / (2*PI*gx[2][ipt]) / ((high-low)*0.8) / Acceptance[ipt] / Smear[ipt] / TrigE[ipt];
      //egy[2][ipt] = gy[2][ipt] * sqrt( pow(eBBCCross/BBCCross,2.) + pow(enpion/npion,2.) + 1./BBCCount + pow(eBBCRatio/BBCRatio,2.)
      //    + pow(eSmear[ipt]/Smear[ipt],2.) + pow(eAcceptance[ipt]/Acceptance[ipt],2.) + pow(eTrigE[ipt]/TrigE[ipt],2.) );
      egy[2][ipt] = gy[2][ipt] * sqrt( pow(enpion/npion,2.)
          + pow(eSmear[ipt]/Smear[ipt],2.) + pow(eAcceptance[ipt]/Acceptance[ipt],2.) + pow(eTrigE[ipt]/TrigE[ipt],2.) );
    }
  }

  c->Print(Form("CrossSection-pion%d-ert%c-part%d.pdf",ispion,97+trig,part));
  delete c;

  return;
}

TGraphErrors **CreateGraph(TFile *f, Int_t part)
{
  const Int_t gn = 30;
  Double_t *gx[3];
  Double_t *gy[3];
  Double_t *egy[3];
  for(Int_t i=0; i<3; i++)
  {
    gx[i] = new Double_t[gn];
    gy[i] = new Double_t[gn];
    egy[i] = new Double_t[gn];
    for(Int_t j=0; j<gn; j++)
    {
      gx[i][j] = 0.;
      gy[i][j] = 0.;
      egy[i][j] = 0.;
    }
  }

  for(Int_t ispion=0; ispion<2; ispion++)
    GenerateGraph(f, ispion, 1, part, gn, gx, gy, egy);

  TGraphErrors **graph = new TGraphErrors*[3];
  for(Int_t i=0; i<3; i++)
    graph[i] = new TGraphErrors(gn, gx[i], gy[i], 0, egy[i]);
  return graph;
}

void draw_CrossSection()
{
  gSystem->Load("libGausProc.so");
  gROOT->ProcessLine(".L ReadGraph.C");
  gROOT->ProcessLine(".L BgGPR.C");
  gROOT->ProcessLine(".L Chi2Fit.C");

  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ertb-cv/total.root");
  TObjArray *Glist = new TObjArray();

  TCanvas *c0 = new TCanvas("c0", "Canvas", 1800, 1200);
  gStyle->SetOptStat(0);
  c0->Divide(3,2);

  TGraphErrors **gr[3];
  for(Int_t part=0; part<3; part++)
  {
    cout << "part " << part << endl;
    gr[part] = CreateGraph(f, part);
    for(Int_t i=0; i<3; i++)
    {
      Glist->Add(gr[part][i]);
      gr[part][i]->SetName(Form("gr_%d",3*i+part));
    }
    gr[part][0]->SetTitle("DirectPhoton fraction");
    gr[part][1]->SetTitle("DirectPhoton Cross Section");
    gr[part][2]->SetTitle("#pi^{0} Cross Section");
    for(Int_t i=0; i<3; i++)
    {
      c0->cd(i+1);
      if(i==0)
      {
        gr[part][i]->GetYaxis()->SetTitle("Ratio");
        gr[part][i]->GetYaxis()->SetRangeUser(0., 1.2);
      }
      else
      {
        gPad->SetLogy();
        gr[part][i]->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3} [pb*GeV^{-2}*c^{-3}]");
        gr[part][i]->GetYaxis()->SetRangeUser(0.5, 1e5);
      }
      gr[part][i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      gr[part][i]->GetYaxis()->SetTitleOffset(1.2);
      gr[part][i]->GetXaxis()->SetRangeUser(0., 30.);
      gr[part][i]->SetMarkerColor(part+1);
      gr[part][i]->SetMarkerStyle(part+20);
      gr[part][i]->SetMarkerSize(1.);
      if(part == 0)
      {
        gr[part][i]->Draw("AP");
      }
      else
      {
        gr[part][i]->Draw("P");
      }
    }
  }

  for(Int_t i=0; i<3; i++)
  {
    c0->cd(i+1);
    TLegend *leg = new TLegend(0.4, 0.7, 0.7, 0.9);
    leg->AddEntry(gr[0][i], "PbScW", "P");
    leg->AddEntry(gr[1][i], "PbScE", "P");
    leg->AddEntry(gr[2][i], "PbGlE", "P");
    leg->Draw();
  }

  const Int_t gn = 30;
  Double_t *gx[3];
  Double_t *gy[3];
  Double_t *egy[3];
  for(Int_t i=0; i<3; i++)
  {
    gx[i] = new Double_t[gn];
    gy[i] = new Double_t[gn];
    egy[i] = new Double_t[gn];
    for(Int_t j=0; j<gn; j++)
    {
      gx[i][j] = 0.;
      gy[i][j] = 0.;
      egy[i][j] = 0.;
    }
  }

  for(Int_t i=0; i<3; i++)
  {
    Double_t xx[gn][3];
    Double_t yy[gn][3];
    Double_t eyy[gn][3];
    for(Int_t ipt=0; ipt<gn; ipt++)
      for(Int_t part=0; part<3; part++)
      {
        gr[part][i]->GetPoint(ipt, xx[ipt][part], yy[ipt][part]);
        eyy[ipt][part] = gr[part][i]->GetErrorY(ipt);
        gx[i][ipt] = xx[ipt][part];
      }
    for(Int_t ipt=0; ipt<gn; ipt++)
      Chi2Fit(3, (Double_t*)yy[ipt], (Double_t*)eyy[ipt], gy[i][ipt], egy[i][ipt]);
  }

  for(Int_t i=0; i<3; i++)
  {
    c0->cd(i+4);
    TGraphErrors *grt = new TGraphErrors(gn, gx[i], gy[i], 0, egy[i]);
    if(i==0)
      grt->SetTitle("DirectPhoton fraction");
    else if(i==1)
      grt->SetTitle("DirectPhoton Cross Section");
    else
      grt->SetTitle("#pi^{0} Cross Section");
    if(i==0)
    {
      grt->GetYaxis()->SetTitle("Ratio");
      grt->GetYaxis()->SetRangeUser(0., 1.2);
    }
    else
    {
      gPad->SetLogy();
      grt->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3} [pb*GeV^{-2}*c^{-3}]");
      grt->GetYaxis()->SetRangeUser(0.1, pow(10,5));
    }
    grt->GetXaxis()->SetTitle("p_{T} [GeV]");
    grt->GetYaxis()->SetTitleOffset(1.2);
    grt->GetXaxis()->SetRangeUser(0., 30.);
    grt->SetMarkerColor(1);
    grt->SetMarkerStyle(20);
    grt->SetMarkerSize(1.);
    grt->Draw("AP");
  }

  c0->Print("CrossSection-ertb.pdf");

  TFile *fout = new TFile("CrossSection-ertb.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
