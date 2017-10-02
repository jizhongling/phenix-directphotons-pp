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

  if(data == 0)
    TH3 *h3_tof = (TH3*)f->Get("h3_tof_raw");
  else if(data == 1)
    TH3 *h3_tof = (TH3*)f->Get("h3_tof");
  TAxis *axis_pt = h3_tof->GetYaxis();

  TCanvas *c = new TCanvas("c", "Canvas", 3600, 3000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  c->Divide(6,5);

  Int_t ipad = 1;
  for(Int_t ipt=1; ipt<gn; ipt++)
  {
    c->cd(ipad++);

    if(part == 0)
      TH1 *h_tof = (TH1*)h3_tof->ProjectionZ("h_tof", 1,6, ipt+1,ipt+1)->Clone();
    else if(part == 1)
      TH1 *h_tof = (TH1*)h3_tof->ProjectionZ("h_tof", 7,8, ipt+1,ipt+1)->Clone();
    Double_t low = axis_pt->GetBinLowEdge(ipt+1);
    Double_t high = axis_pt->GetBinUpEdge(ipt+1);
    h_tof->SetTitle(Form("pT: %4.2f-%4.2f",low, high));

    TF1 *fn2 = new TF1("fn2", "gaus", -100., 100.);

    h_tof->Fit(fn2, "Q0");
    fn2->SetParameters( fn2->GetParameters() );
    h_tof->Fit(fn2, "QE");
    h_tof->DrawCopy();

    Double_t scale = 1.;
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

  c->Print(Form("ToF_Calib-part%d-data%d.pdf",part,data));
  delete c;

  TGraphErrors **graph = new TGraphErrors*[2];
  for(Int_t i=0; i<2; i++)
    graph[i] = new TGraphErrors(gn, gx[i], gy[i], 0, egy[i]);
  return graph;
}

void draw_ToF_Calib()
{
  TFile *f_sim = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonNode-histo.root");
  TFile *f_data = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonNode-histo.root");

  mc(0, 2,2);
  gStyle->SetStatY(1.);
  gStyle->SetStatH(0.13);

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
      gr_sim[part][i]->GetXaxis()->SetRangeUser(0., 30.);
      if(i == 0)
      {
        gr_sim[part][i]->GetYaxis()->SetTitle("ToF [ns]");
        gr_sim[part][i]->GetYaxis()->SetRangeUser(-20., 20.);
      }
      else if(i == 1)
      {
        gr_sim[part][i]->GetYaxis()->SetTitle("#sigma_{ToF} [ns]");
        gr_sim[part][i]->GetYaxis()->SetRangeUser(0., 10.);
      }
      gr_sim[part][i]->GetYaxis()->SetTitleOffset(1.5);
      gr_sim[part][i]->SetMarkerStyle(4);
      gr_sim[part][i]->SetMarkerColor(kBlue);
      gr_data[part][i]->SetMarkerStyle(4);
      gr_data[part][i]->SetMarkerColor(kRed);
      gr_sim[part][i]->Draw("AP");
      gr_data[part][i]->Draw("P");
    }
  }
  gr_sim[0][0]->SetTitle("PbSc ToF");
  gr_sim[1][0]->SetTitle("PbGl ToF");
  gr_sim[0][1]->SetTitle("PbSc #sigma_{ToF}");
  gr_sim[1][1]->SetTitle("PbGl #sigma_{ToF}");

  //TF1 *f1 = new TF1("f1", "[0]*sqrt(x-[1])+pol2(2)", 0.5,30.);
  //gr_data[0][0]->Fit(f1, "R");

  c0->Print("ToF_Calib.pdf");
}
