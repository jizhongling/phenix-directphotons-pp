TGraphErrors **CreateGraph(TFile *f, Int_t part, Int_t data)
{
  const Int_t gn = 30;
  Double_t *gx[2];
  Double_t *gy[2];
  Double_t *egy[2];
  for(Int_t i=0; i<2; i++)
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

  Int_t olddata = data;
  data = 0;

  if(data == 0)
  {
    THnSparse *hn_minv = (THnSparse*)f->Get("hn_pion");
    TAxis *axis_pt = hn_minv->GetAxis(1);
    if(part == 0)
      hn_minv->GetAxis(3)->SetRange(1,6);
    else if(part == 1)
      hn_minv->GetAxis(3)->SetRange(7,8);
  }
  else if(data == 1)
  {
    TH3 *h3_minv = (TH3*)f->Get("h3_inv_mass_pi0calib");
    TAxis *axis_pt = h3_minv->GetYaxis();
  }

  mc(1, 6,5);
  mc(2, 6,5);

  Int_t ipad = 1;
  for(Int_t ipt=1; ipt<gn; ipt++)
  {
    mcd(1, ipad);

    axis_pt->SetRange(ipt+1, ipt+1);
    if(data == 0)
      TH1 *h_minv = hn_minv->Projection(2);
    else if(data == 1)
    {
      if(part ==  0)
        TH1 *h_minv = (TH1*)h3_minv->ProjectionZ("h_minv", 1,6, ipt+1,ipt+1)->Clone();
      else if(part ==  1)
        TH1 *h_minv = (TH1*)h3_minv->ProjectionZ("h_minv", 7,8, ipt+1,ipt+1)->Clone();
    }
    aset(h_minv);
    Double_t low = axis_pt->GetBinLowEdge(ipt+1);
    Double_t high = axis_pt->GetBinUpEdge(ipt+1);
    h_minv->SetTitle(Form("pT: %4.2f-%4.2f",low, high));

    TF1 *fn1 = new TF1("fn1", "gaus", 0., 0.5);
    TF1 *fn2 = new TF1("fn2", "gaus(0)+pol1(3)", 0., 0.5);

    if(olddata==0)
    {
      h_minv->Fit(fn1, "WLQ0", "", 0.112, 0.162);
      fn2->SetParameters( fn1->GetParameters() );
      h_minv->Fit(fn2, "WLQ0", "", 0.047, 0.227);
      h_minv->Fit(fn2, "WLQE", "", 0.047, 0.227);
    }
    else if(olddata==1)
    {
      h_minv->Fit(fn1, "Q0", "", 0.112, 0.162);
      fn2->SetParameters( fn1->GetParameters() );
      h_minv->Fit(fn2, "Q0", "", 0.047, 0.227);
      h_minv->Fit(fn2, "QE", "", 0.047, 0.227);
    }
    h_minv->DrawCopy("E");

    TH1 *h_pull = new TH1F("h_pull", "Pull statistics for fitting", 180,0.047,0.227);
    for(Int_t i=1; i<=180; i++)
    {
      Double_t meas = h_minv->GetBinContent(47+i);
      Double_t emeas = h_minv->GetBinError(47+i);
      Double_t est = fn2->Eval(0.047+0.001*i);
      if(emeas>0.)
        Double_t pull = (meas - est) / emeas;
      h_pull->SetBinContent(i, pull);
    }

    mcd(2, ipad);
    aset(h_pull);
    h_pull->SetTitle(Form("pT: %4.2f-%4.2f",low, high));
    h_pull->GetYaxis()->SetRangeUser(-10.,10.);
    h_pull->DrawCopy();
    delete h_pull;

    Double_t scale = 1.;
    if( fn2->GetNDF() > 0 )
      scale = sqrt( fn2->GetChisquare() / fn2->GetNDF() );

    gx[0][ipt] = axis_pt->GetBinCenter(ipt+1);
    gy[0][ipt] = fn2->GetParameter(1);
    egy[0][ipt] = fn2->GetParError(1) * scale;

    gx[1][ipt] = axis_pt->GetBinCenter(ipt+1);
    gy[1][ipt] = fn2->GetParameter(2);
    egy[1][ipt] = fn2->GetParError(2) * scale;

    ipad++;
    delete h_minv;
  }

  c1->Print(Form("Pi0Peak-part%d-data%d.pdf",part,olddata));
  delete c1;

  c2->Print(Form("Pi0Peak-part%d-data%d-pull.pdf",part,olddata));
  delete c2;

  if(data == 0)
    delete hn_minv;
  else if(data == 1)
    delete h3_minv;

  TGraphErrors **graph = new TGraphErrors*[2];
  for(Int_t i=0; i<2; i++)
    graph[i] = new TGraphErrors(gn, gx[i], gy[i], 0, egy[i]);
  return graph;
}

void draw_Pi0Peak()
{
  TFile *f_sim = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/MissingRatio-histo.root");
  TFile *f_data = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  mc(0, 2,2);

  TGraphErrors **gr_sim[2];
  TGraphErrors **gr_data[2];
  for(Int_t part=0; part<2; part++)
  {
    gr_sim[part] = CreateGraph(f_sim, part, 0);
    gr_data[part] = CreateGraph(f_data, part, 1);
    for(Int_t i=0; i<2; i++)
    {
      mcd(0, 2*i+part+1);
      aset(gr_sim[part][i]);
      gr_sim[part][i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      gr_sim[part][i]->GetXaxis()->SetRangeUser(0., 20.);
      if(i == 0)
      {
        gr_sim[part][i]->GetYaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
        gr_sim[part][i]->GetYaxis()->SetRangeUser(0.13, 0.145);
      }
      else if(i == 1)
      {
        gr_sim[part][i]->GetYaxis()->SetTitle("#sigma_{#gamma#gamma} [GeV]");
        gr_sim[part][i]->GetYaxis()->SetRangeUser(0., 0.015);
      }
      style(gr_sim[part][i], 4, kBlue);
      style(gr_data[part][i], 4, kRed);
      gr_sim[part][i]->Draw("AP");
      gr_data[part][i]->Draw("P");
    }
  }
  gr_sim[0][0]->SetTitle("PbSc m_{#gamma#gamma}");
  gr_sim[1][0]->SetTitle("PbGl m_{#gamma#gamma}");
  gr_sim[0][1]->SetTitle("PbSc #sigma_{#gamma#gamma}");
  gr_sim[1][1]->SetTitle("PbGl #sigma_{#gamma#gamma}");

  c0->Print("Pi0Peak.pdf");
}
