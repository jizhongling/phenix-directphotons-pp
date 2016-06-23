TGraphErrors **CreateGraph(TFile *f, Int_t part)
{
  const Double_t PI = 3.1415927;
  const Double_t BBCCount = 4.193255 * pow(10,12);  // 4193255288313
  const Double_t BBCRatio = 0.54;
  const Double_t eBBCRatio = 0.09;
  const Double_t BBCCross = 32.5 * pow(10,9);
  const Double_t eBBCCross = 3.25 * pow(10,9);
  const Double_t BR[30] = {0.0444823,0.155126,0.189877,0.209775,0.227416,0.235346,0.242619,0.258304,0.254693,0.219828,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313};
  const Double_t eBR[30] = {1.60171e-05,0.00010065,0.000293243,0.000651539,0.00129034,0.00232753,0.00398618,0.00663705,0.0101993,0.0144502,0.0228892,0.0291354,0.0348798,0.0512812,0.0699279,0.078914,0.111523,0.142477,0.247947,0.213559,0.156142,0.154185,0.30816,0.8415,0.8415};

  const Int_t gn = 30;
  Double_t *gx = new Double_t[gn];
  Double_t *gy[3];
  Double_t *egy[3];
  for(Int_t i=0; i<3; i++)
  {
    gy[i] = new Double_t[gn];
    egy[i] = new Double_t[gn];
    for(Int_t j=0; j<gn; j++)
    {
      gx[j] = 0.;
      gy[i][j] = 0.;
      egy[i][j] = 0.;
    }
  }

  if(part == 0)
  {
    const Int_t sector_low = 1;
    const Int_t sector_high = 4;
    const Double_t MissR[30] = {0.576394, 1.19587, 1.61874, 1.24187, 1.43121, 1.12439, 1.12484, 1.16195, 0.977662, 1.06683, 0.992115, 0.991668, 0.986668, 0.904868, 0.844696, 0.881821, 0.823503, 0.673003, 0.713694, 0.601754, 0.734083, 0.920672, 0.973042, 0.735017, 1.0551, 1.23027, 1.24595, 2.16717, 2.67301, 3.33668};
    const Double_t eMissR[30] = {0.427523, 0.0848828, 0.119692, 0.101558, 0.125594, 0.108244, 0.115694, 0.131099, 0.114497, 0.126486, 0.117178, 0.113272, 0.119361, 0.099389, 0.0946024, 0.101291, 0.0989192, 0.0789351, 0.0930587, 0.0814176, 0.0580832, 0.0843133, 0.0910573, 0.073874, 0.116943, 0.14423, 0.161901, 0.317037, 0.420699, 0.622627};
    const Double_t Acceptance[30] = {0.00270575,0.218638,0.197805,0.178892,0.16642,0.155844,0.145005,0.136193,0.13182,0.122922,0.117749,0.110509,0.107532,0.104263,0.0995625,0.0961354,0.0890214,0.089119,0.0846746,0.0838922,0.0751721,0.0641918,0.0576509,0.0537159,0.0548628,0.0617032,0.0758218,0.109818,0.191812,0.403237};
    const Double_t eAcceptance[30] = {0.000178645,0.00146455,0.00153235,0.00155725,0.00158723,0.00161168,0.00163235,0.00163417,0.00166831,0.00166341,0.00168533,0.00168287,0.00171492,0.00172876,0.00175484,0.00175923,0.00173022,0.00178385,0.00178154,0.00181729,0.000905022,0.000922492,0.000968387,0.00104042,0.00117003,0.00140013,0.00178229,0.0025309,0.00418803,0.00921339};
    const Double_t TrigE[30] = {8.09532e-06,1.88351e-05,9.80708e-05,0.000389597,0.00105885,0.00353447,0.0152754,0.073618,0.227824,0.470076,0.679126,0.884265,0.909871,0.953846,0.952096,0.954128,0.897059,0.958333,0.965517,1,0.953488,0.916667,1,1,1,1,1,1,1,1};
    const Double_t eTrigE[30] = {2.47946e-06,7.6374e-07,5.0765e-06,2.25507e-05,7.22226e-05,0.000230904,0.000770765,0.00250532,0.00587997,0.00991908,0.012883,0.0127204,0.0153001,0.016987,0.0227901,0.0298569,0.0510686,0.0523567,0.074913,0.0615421,0.0581073,0.166603,0.264349,0.369032,0.417333,0.601879,0.8415,0.8415,0.601879,0.8415};
  }
  else if(part == 1)
  {
    const Int_t sector_low = 5;
    const Int_t sector_high = 6;
    const Double_t MissR[30] = {2.58414, 1.15151, 1.2534, 1.11775, 1.10035, 1.03849, 0.852603, 1.23974, 0.843985, 0.957385, 0.881202, 1.09943, 0.565636, 0.6154, 0.63001, 0.518561, 0.35348, 0.409721, 0.307164, 0.35066, 0.578968, 0.562609, 0.608813, 0.489362, 0.791451, 0.600303, 0.891448, 1.84624, 1.83777, 3.5327};
    const Double_t eMissR[30] = {1.91341, 0.112286, 0.12483, 0.123574, 0.1316, 0.139969, 0.119626, 0.185305, 0.142144, 0.157711, 0.147912, 0.193396, 0.0949982, 0.097322, 0.109026, 0.0866015, 0.0663085, 0.0762913, 0.0613931, 0.0704716, 0.0644004, 0.0741542, 0.0888177, 0.0775301, 0.130682, 0.106328, 0.166021, 0.359809, 0.341299, 0.880618};
    const Double_t Acceptance[30] = {0.000856647,0.0757418,0.0709951,0.068296,0.0648773,0.0604169,0.0583629,0.0563723,0.0534483,0.052139,0.0478959,0.0478119,0.042894,0.0433451,0.0391754,0.0371019,0.0375537,0.0386295,0.0366648,0.0342728,0.0287382,0.0254594,0.0244883,0.0223833,0.0220564,0.0253672,0.0329152,0.05056,0.102242,0.232299};
    const Double_t eAcceptance[30] = {0.000105447,0.000942802,0.00099396,0.00103134,0.00105654,0.0010662,0.0010947,0.00110698,0.00111808,0.00113515,0.0011265,0.00115509,0.00113288,0.00116323,0.00114978,0.0011414,0.00116787,0.00121978,0.00121647,0.00120787,0.000577929,0.000598059,0.000648083,0.000690069,0.000764125,0.000926882,0.00121644,0.0017943,0.00325304,0.00801897};
    const Double_t TrigE[30] = {8.45417e-06,2.15626e-05,0.00012377,0.00043417,0.00106439,0.0057713,0.0278967,0.0999358,0.240773,0.414296,0.604027,0.805699,0.828,0.888889,0.95122,0.941176,0.934783,0.96875,0.944444,0.909091,0.9375,0.888889,1,1,1,1,1,1,1,1};
    const Double_t eTrigE[30] = {3.86319e-06,1.17435e-06,7.86683e-06,3.27608e-05,0.000100533,0.000399679,0.00139129,0.00394678,0.00820872,0.0137363,0.0187326,0.0222699,0.0272703,0.0315277,0.0369117,0.0440953,0.0594246,0.068259,0.116475,0.179384,0.076593,0.21166,0.601879,0.601879,0.8415,0.8415,0.8415,0.8415,0.8415,0.8415};
  }
  else if(part == 2)
  {
    const Int_t sector_low = 7;
    const Int_t sector_high = 8;
    const Double_t MissR[30] = {0.579534, 0.957785, 1.72653, 1.22189, 0.697323, 0.850314, 0.994086, 1.06729, 1.01035, 0.748298, 0.672339, 0.598552, 0.717492, 0.56897, 0.582423, 0.468526, 0.35638, 0.391049, 0.45989, 0.516432, 0.504037, 0.526388, 0.374648, 0.47018, 0.380933, 0.404627, 0.427724, 0.583086, 0.42485, 0.799791};
    const Double_t eMissR[30] = {0.626164, 0.0852248, 0.158208, 0.125625, 0.0846392, 0.106273, 0.132056, 0.152662, 0.148727, 0.119227, 0.111418, 0.0945248, 0.110301, 0.0920305, 0.0834909, 0.0725711, 0.0584358, 0.0597026, 0.0784037, 0.0916945, 0.052552, 0.0587681, 0.0471476, 0.0601228, 0.053128, 0.0601148, 0.0705125, 0.0999384, 0.0913654, 0.309123};
    const Double_t Acceptance[30] = {0.00113871,0.0870023,0.0862187,0.0871609,0.0832814,0.0831192,0.0813722,0.0785139,0.0765003,0.0722308,0.0711005,0.0695953,0.0694432,0.0663325,0.0671815,0.0648339,0.0668015,0.0625377,0.0644579,0.0611808,0.0582463,0.0534358,0.0494103,0.0467121,0.0452257,0.0463664,0.0513841,0.0599462,0.078631,0.121376};
    const Double_t eAcceptance[30] = {0.000119849,0.00100339,0.00108485,0.00115114,0.00118284,0.00123202,0.0012729,0.00128699,0.00131656,0.00131739,0.00134984,0.00137175,0.00141318,0.00141419,0.00147323,0.00147579,0.00152182,0.00152221,0.00157656,0.00157718,0.000805297,0.000847709,0.000901676,0.000975209,0.00107009,0.00122798,0.00149432,0.00193841,0.00290379,0.00629987};
    const Double_t TrigE[30] = {3.86368e-05,6.81334e-05,0.000125886,0.000260003,0.000887477,0.0019329,0.0187945,0.0556734,0.128918,0.239865,0.351272,0.488696,0.662824,0.638393,0.690141,0.681818,0.760563,0.790698,0.709677,0.757576,0.808511,0.866667,0.625,0.5,0.5,1,1,1,1,1};
    const Double_t eTrigE[30] = {7.23467e-06,2.04833e-06,8.44487e-06,2.56335e-05,8.06052e-05,0.000140612,0.00106058,0.00276853,0.00577574,0.0105686,0.0155742,0.0217299,0.027293,0.0349367,0.0436609,0.0573316,0.06149,0.0808377,0.103954,0.0976443,0.0751301,0.149614,0.235052,0.417333,0.417333,0.601879,0.8415,0.8415,0.8415,0.8415};
  }
  else
  {
    cerr << "Wrong part " << part << endl;
    exit(1);
  }

  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  hn_1photon->GetAxis(0)->SetRange(sector_low, sector_high);
  TH1 *h_1photon = hn_1photon->Projection(1);

  TH2 *h2_sig_extra = (TH2*)f->Get("h2_sig_extra");
  TH1 *h_sig_extra = h2_sig_extra->ProjectionY("h_sig_extra", sector_low, sector_high);

  TH2 *h2_bg_extra = (TH2*)f->Get("h2_bg_extra");
  TH1 *h_bg_extra = h2_bg_extra->ProjectionY("h_bg_extra", sector_low, sector_high);

  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");

  TAxis *axis0 = hn_2photon->GetAxis(0);
  TAxis *axis1 = hn_2photon->GetAxis(1);
  TAxis *axis2 = hn_2photon->GetAxis(2);

  axis0->SetRange(sector_low, sector_high);

  Int_t bin047 = axis2->FindBin(0.047);
  Int_t bin067 = axis2->FindBin(0.067);
  Int_t bin087 = axis2->FindBin(0.087);
  Int_t bin097 = axis2->FindBin(0.097);
  Int_t bin112 = axis2->FindBin(0.112);
  Int_t bin162 = axis2->FindBin(0.162);
  Int_t bin177 = axis2->FindBin(0.177);
  Int_t bin187 = axis2->FindBin(0.187);
  Int_t bin212 = axis2->FindBin(0.212);
  Int_t bin227 = axis2->FindBin(0.227);

  const Int_t nData = 256;
  vector<Double_t> x(nData), y(nData), sigma_y(nData);

  TCanvas *c = new TCanvas("c", "Canvas", 2400, 2000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c->Divide(6,5);

  Int_t ipad = 1;
  for(Int_t ipt=0; ipt<30; ipt++)
  {
    c->cd(ipad++);

    char buf[100];
    Double_t low = axis1->GetBinLowEdge(ipt+1);
    Double_t high = axis1->GetBinLowEdge(ipt+2);
    sprintf(buf, "p_{T} %3.1f-%3.1f", low, high);

    axis1->SetRange(ipt+1,ipt+1);
    TH1 *h_inv_mass = hn_2photon->Projection(2);
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
      Double_t xx = axis2->GetBinCenter(ib);
      Double_t yy = h_inv_mass->GetBinContent(ib);
      Double_t sigma_yy = h_inv_mass->GetBinError(ib);
      x.push_back(xx);
      y.push_back(yy);
      sigma_y.push_back(sigma_yy);
    }

    for(Int_t ib=bin187; ib<bin212; ib++)
    {
      Double_t xx = axis2->GetBinCenter(ib);
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
    TF1 *fn2 = new TF1("fn2", "gaus(0)+pol3(3)", 0., 0.5);
    TF1 *fn3 = new TF1("fn3", "pol3", 0., 0.5);

    Double_t par[10];
    h_inv_mass->Fit(fn1, "Q0", "", 0.112, 0.162);
    fn2->SetParameters( fn1->GetParameters() );
    h_inv_mass->Fit(fn2, "Q", "", 0.047, 0.227);
    fn2->GetParameters(par);
    fn3->SetParameters(par[3], par[4], par[5], par[6]);

    for(Int_t ib=bin047; ib<bin097; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    for(Int_t ib=bin177; ib<bin227; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    nbgside /= 2.;

    for(Int_t ib=bin112; ib<bin162; ib++)
    {
      npair += h_inv_mass->GetBinContent(ib);
      Double_t bincenter = axis2->GetBinCenter(ib);
      nbgfit += fn3->Eval(bincenter);
    }

    Double_t nfit = npair - nbgfit;
    Double_t enfit = sqrt(1./npair + 1./nbgfit);
    Double_t nsub = npair - nbgside;
    Double_t ndiff = fabs(nbgfit - nbgside - nextra);
    Double_t ngpr = npair - nbggpr;
    Double_t nerror = fabs(nfit - ngpr);

    Double_t npion = nfit * ( 1. + MissR[ipt] );
    Double_t enpion = npion * sqrt( pow(enfit/nfit,2.) + pow(eMissR[ipt]/(1.+MissR[ipt]),2.) );
    Double_t ndecay = npion * ( 1. + BR[ipt] );
    Double_t endecay = ndecay * sqrt( pow(enpion/npion,2.) + pow(eBR[ipt]/(1.+BR[ipt]),2.) );
    Double_t ndirect = nphoton - ndecay;
    Double_t endirect = sqrt( 1./nphoton + pow(endecay/ndecay,2.) );

    gx[ipt] = (low + high) / 2.;
    gy[0][ipt] = ndirect / nphoton;
    egy[0][ipt] = gy[0][ipt] * sqrt( pow(endirect/ndirect,2.) + 1./nphoton );
    gy[1][ipt] = BBCCross * (ndirect/(BBCCount/BBCRatio)) / (2*PI*gx[ipt]) / ((high-low)*0.5) / Acceptance[ipt] / TrigE[ipt];
    egy[1][ipt] = gy[1][ipt] * sqrt( pow(eBBCCross/BBCCross,2.) + pow(endirect/ndirect,2.) + 1./BBCCount + pow(eBBCRatio/BBCRatio,2.)
        + pow(eAcceptance[ipt]/Acceptance[ipt],2.) + pow(eTrigE[ipt]/TrigE[ipt],2.) );
    gy[2][ipt] = BBCCross * (npion/2./(BBCCount/BBCRatio)) / (2*PI*gx[ipt]) / ((high-low)*0.5) / Acceptance[ipt] / TrigE[ipt];
    egy[2][ipt] = gy[2][ipt] * sqrt( pow(eBBCCross/BBCCross,2.) + pow(enpion/npion,2.) + 1./BBCCount + pow(eBBCRatio/BBCRatio,2.)
        + pow(eAcceptance[ipt]/Acceptance[ipt],2.) + pow(eTrigE[ipt]/TrigE[ipt],2.) );

    cout << "pT=" << low << "-" << high << "\tnphoton=" << nphoton
      << "\tnpair=" << npair << "\tnfit=" << nfit << "\tnsub=" << nsub
      << "\tndiff=" << ndiff/nfit*100. << "%\tngpr=" << ngpr
      << "\tnerror=" << nerror/nfit*100. << "%" << endl;
  }

  if(part == 0)
  {
    c->Print("CrossSection-PbScW.pdf");
    delete c;
  }
  if(part == 1)
  {
    c->Print("CrossSection-PbScE.pdf");
    delete c;
  }
  else if(part == 2)
  {
    c->Print("CrossSection-PbGlE.pdf");
    delete c;
  }

  TGraphErrors **graph = new TGraphErrors*[3];
  for(Int_t i=0; i<3; i++)
    graph[i] = new TGraphErrors(gn, gx, gy[i], 0, egy[i]);
  return graph;
}

void draw_CrossSection()
{
  gSystem->Load("libGausProc.so");
  gROOT->ProcessLine(".L BgGPR.C");
  gROOT->ProcessLine(".L Chi2Fit.C");

  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/total.root");

  TCanvas *c0 = new TCanvas("c0", "Canvas", 1800, 1200);
  gStyle->SetOptStat(0);
  c0->Divide(3,2);

  TGraphErrors **gr[3];
  for(Int_t part=0; part<3; part++)
  {
    cout << "part " << part << endl;
    gr[part] = CreateGraph(f, part);
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
        gr[part][i]->GetYaxis()->SetRangeUser(0.1, pow(10,5));
      }
      gr[part][i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      gr[part][i]->GetYaxis()->SetTitleOffset(1.2);
      gr[part][i]->GetXaxis()->SetRangeUser(0., 20.);
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
  Double_t *gx = new Double_t[gn];
  Double_t *gy[3];
  Double_t *egy[3];
  for(Int_t i=0; i<3; i++)
  {
    gy[i] = new Double_t[gn];
    egy[i] = new Double_t[gn];
    for(Int_t j=0; j<gn; j++)
    {
      gx[j] = 0.;
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
        gx[ipt] = xx[ipt][part];
      }
    for(Int_t ipt=0; ipt<gn; ipt++)
      Chi2Fit(3, (Double_t*)yy[ipt], (Double_t*)eyy[ipt], gy[i][ipt], egy[i][ipt]);
  }

  for(Int_t i=0; i<3; i++)
  {
    c0->cd(i+4);
    TGraphErrors *grt = new TGraphErrors(gn, gx, gy[i], 0, egy[i]);
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
    grt->GetXaxis()->SetRangeUser(0., 20.);
    grt->SetMarkerColor(1);
    grt->SetMarkerStyle(20);
    grt->SetMarkerSize(1.);
    grt->Draw("AP");
  }

  c0->Print("CrossSection.pdf");
}
