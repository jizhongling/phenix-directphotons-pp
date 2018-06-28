TGraphErrors **CreateGraph(TFile *f, int part, int calib)
{
  const int gn = 30;
  double *gx[2];
  double *gy[2];
  double *egy[2];
  for(int i=0; i<2; i++)
  {
    gx[i] = new double[gn];
    gy[i] = new double[gn];
    egy[i] = new double[gn];
    for(int j=0; j<gn; j++)
    {
      gx[i][j] = 0.;
      gy[i][j] = 0.;
      egy[i][j] = 0.;
    }
  }

  if(calib == 0)
    TH3 *h3_tof = (TH3*)f->Get("h3_tof_raw");
  else if(calib == 1)
    TH3 *h3_tof = (TH3*)f->Get("h3_tof");
  TAxis *axis_pt = h3_tof->GetYaxis();

  mc(1, 6,5);
  int ipad = 1;
  for(int ipt=1; ipt<gn; ipt++)
  {
    mcd(1, ipad++);

    if(part == 0)
      TH1 *h_tof = (TH1*)h3_tof->ProjectionZ("h_tof", 1,6, ipt+1,ipt+1)->Clone();
    else if(part == 1)
      TH1 *h_tof = (TH1*)h3_tof->ProjectionZ("h_tof", 7,8, ipt+1,ipt+1)->Clone();
    double low = axis_pt->GetBinLowEdge(ipt+1);
    double high = axis_pt->GetBinUpEdge(ipt+1);
    h_tof->SetTitle(Form("pT: %3.1f-%3.1f GeV",low, high));

    TF1 *fn2 = new TF1("fn2", "gaus", -30., 30.);

    const double range_left[2] = {-30., -10.};
    const double range_right[2] = {10., 20.};
    fn2->SetParameters(1e3, 0., 1.);
    h_tof->Fit(fn2, "RQ0", "", range_left[part], range_right[part]);
    fn2->SetParameters( fn2->GetParameters() );
    h_tof->Fit(fn2, "RQE", "", range_left[part], range_right[part]);
    h_tof->DrawCopy();

    double scale = 1.;
    if( fn2->GetNDF() > 0 )
      scale = sqrt( fn2->GetChisquare() / fn2->GetNDF() );

    gx[0][ipt] = axis_pt->GetBinCenter(ipt+1);
    gy[0][ipt] = fn2->GetParameter(1);
    egy[0][ipt] = fn2->GetParError(1) * scale;

    gx[1][ipt] = axis_pt->GetBinCenter(ipt+1);
    gy[1][ipt] = fn2->GetParameter(2);
    egy[1][ipt] = fn2->GetParError(2) * scale;

    delete h_tof;
  }

  c1->Print(Form("plots/ToF_Calib-part%d-calib%d.pdf",part,calib));
  delete c1;

  TGraphErrors **graph = new TGraphErrors*[2];
  for(int i=0; i<2; i++)
    graph[i] = new TGraphErrors(gn, gx[i], gy[i], 0, egy[i]);
  return graph;
}

void draw_ToF_Calib_pT()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  mc(0, 2,2);
  gStyle->SetStatY(1.);
  gStyle->SetStatH(0.13);

  TGraphErrors **gr_raw[2];
  TGraphErrors **gr_calib[2];
  for(int part=0; part<2; part++)
  {
    gr_raw[part] = CreateGraph(f, part, 0);
    gr_calib[part] = CreateGraph(f, part, 1);
    for(int i=0; i<2; i++)
    {
      mcd(0, 2*i+part+1);
      aset(gr_raw[part][0], "p_{T} [GeV]","ToF [ns]", 0.,30., -20.,20.);
      aset(gr_raw[part][1], "p_{T} [GeV]","#sigma_{ToF} [ns]", 0.,30., 0.,10.);
      style(gr_raw[part][i], 4, kBlue);
      style(gr_calib[part][i], 4, kRed);
      gr_raw[part][i]->Draw("AP");
      gr_calib[part][i]->Draw("P");
    }
  }
  gr_raw[0][0]->SetTitle("PbSc ToF");
  gr_raw[1][0]->SetTitle("PbGl ToF");
  gr_raw[0][1]->SetTitle("PbSc #sigma_{ToF}");
  gr_raw[1][1]->SetTitle("PbGl #sigma_{ToF}");

  //TF1 *f1 = new TF1("f1", "[0]*sqrt(x-[1])+pol2(2)", 0.5,30.);
  //gr_calib[0][0]->Fit(f1, "R");

  c0->Print("plots/ToF_Calib.pdf");
}
