void draw_InvMass_Check()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/total.root");
  THnSparse *hn_minv = (THnSparse*)f->Get("hn_minv");

  TAxis *axis_sec = hn_minv->GetAxis(0);
  TAxis *axis_pt = hn_minv->GetAxis(1);

  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  double pTbins[30];
  double E1[2][30], E2[2][30], Theta[2][30], M[2][30];
  double eE1[2][30], eE2[2][30], eTheta[2][30], eM[2][30];

  for(int part=0; part<2; part++)
  {
    axis_sec->SetRange(secl[part],sech[part]);

    for(int ipt=0; ipt<30; ipt++)
    {
      axis_pt->SetRange(ipt+1,ipt+1);
      pTbins[ipt] = axis_pt->GetBinCenter(ipt+1);

      TH1 *h_tmp = hn_minv->Projection(2);
      E1[part][ipt] = h_tmp->GetMean();
      eE1[part][ipt] = h_tmp->GetMeanError();
      delete h_tmp;

      TH1 *h_tmp = hn_minv->Projection(3);
      E2[part][ipt] = h_tmp->GetMean();
      eE2[part][ipt] = h_tmp->GetMeanError();
      delete h_tmp;

      TH1 *h_tmp = hn_minv->Projection(4);
      Theta[part][ipt] = h_tmp->GetMean();
      eTheta[part][ipt] = h_tmp->GetMeanError();
      delete h_tmp;

      M[part][ipt] = sqrt( 2. * E1[part][ipt] * E2[part][ipt] ) * Theta[part][ipt];
      eM[part][ipt] = sqrt( pow(eE1[part][ipt]/E1[part][ipt],2) + pow(eE2[part][ipt]/E2[part][ipt],2) + pow(eTheta[part][ipt]/Theta[part][ipt],2) ) / 2.;
    }
  }

  TGraphErrors *gr[2][4];
  for(int part=0; part<2; part++)
  {
    TCanvas *c = new TCanvas("c", "", 2400, 2400);
    gStyle->SetOptStat(0);
    c->Divide(2,2);

    gr[part][0] = new TGraphErrors(30, pTbins, (double*)E1[part], 0, (double*)eE1[part]);
    gr[part][1] = new TGraphErrors(30, pTbins, (double*)E2[part], 0, (double*)eE2[part]);
    gr[part][2] = new TGraphErrors(30, pTbins, (double*)Theta[part], 0, (double*)eTheta[part]);
    gr[part][3] = new TGraphErrors(30, pTbins, (double*)M[part], 0, (double*)eM[part]);

    gr[part][0]->SetTitle("#bar{E_{1}} (GeV)");
    gr[part][1]->SetTitle("#bar{E_{2}} (GeV)");
    gr[part][2]->SetTitle("#bar{#sqrt{1-cos(#theta)}}");
    gr[part][3]->SetTitle("#sqrt{2#bar{E_{1}}#bar{E_{2}}}#bar{#sqrt{1-cos#theta}} (GeV)");

    int ipad = 1;
    for(int igr=0; igr<3; igr++)
    {
      c->cd(ipad++);
      gr[part][igr]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      gr[part][igr]->SetMarkerStyle(20);
      gr[part][igr]->SetMarkerColor(1);
      gr[part][igr]->SetMarkerSize(2);
      gr[part][igr]->Draw("AP");
    }

    c->cd(ipad++);
    TLatex Tl;
    Tl.SetTextSize(0.066);
    Tl.DrawLatex(0.2, 0.5, "m_{inv} = #sqrt{2E_{1}E_{2}(1-cos#theta)}");

    c->Print(Form("plots/InvMass-Check-part%d.pdf",part));
    delete c;
  }
}
