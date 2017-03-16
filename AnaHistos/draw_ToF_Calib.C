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

  //TH1::SetDefaultSumw2();

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
      TH1 *h_tof = h3_tof->ProjectionZ("h_tof", 1,6, ipt+1,ipt+1);
    else if(part == 1)
      TH1 *h_tof = h3_tof->ProjectionZ("h_tof", 7,8, ipt+1,ipt+1);
    Double_t low = axis_pt->GetBinLowEdge(ipt+1);
    Double_t high = axis_pt->GetBinUpEdge(ipt+1);
    h_tof->SetTitle(Form("pT: %4.2f-%4.2f",low, high));

    TF1 *fn2 = new TF1("fn2", "gaus", -100., 100.);

    h_tof->Fit(fn2, "Q0");
    fn2->SetParameters( fn2->GetParameters() );
    h_tof->Fit(fn2, "QE");

    Double_t scale = 1.;
    if( fn2->GetNDF() > 0 )
      scale = sqrt( fn2->GetChisquare() / fn2->GetNDF() );

    gx[0][ipt] = axis_pt->GetBinCenter(ipt+1);
    gy[0][ipt] = fn2->GetParameter(1);
    egy[0][ipt] = fn2->GetParError(1) * scale;

    gx[1][ipt] = axis_pt->GetBinCenter(ipt+1);
    gy[1][ipt] = fn2->GetParameter(2);
    egy[1][ipt] = fn2->GetParError(2) * scale;
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
  TFile *f_sim = new TFile("/phenix/plhf/zji/taxi/Run13pp510ERT/10853/data/total.root");
  TFile *f_data = new TFile("/phenix/plhf/zji/taxi/Run13pp510ERT/10853/data/total.root");

  TCanvas *c0 = new TCanvas("c0", "Canvas", 1200, 1200);
  gStyle->SetOptStat(0);
  c0->Divide(2,2);

  TGraphErrors **gr_sim[2];
  TGraphErrors **gr_data[2];
  for(Int_t part=0; part<2; part++)
  {
    gr_sim[part] = CreateGraph(f_sim, part, 0);
    gr_data[part] = CreateGraph(f_data, part, 1);
    for(Int_t i=0; i<2; i++)
    {
      c0->cd(2*i+part+1);
      gr_sim[part][i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      gr_sim[part][i]->GetXaxis()->SetRangeUser(0., 30.);
      if(i == 0)
      {
        gr_sim[part][i]->GetYaxis()->SetTitle("Tof [ns]");
        gr_sim[part][i]->GetYaxis()->SetRangeUser(-20., 20.);
      }
      else if(i == 1)
      {
        gr_sim[part][i]->GetYaxis()->SetTitle("#sigma_{Tof} [ns]");
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
  gr_sim[0][0]->SetTitle("PbSc m_{#gamma#gamma}");
  gr_sim[1][0]->SetTitle("PbGl m_{#gamma#gamma}");
  gr_sim[0][1]->SetTitle("PbSc #sigma_{#gamma#gamma}");
  gr_sim[1][1]->SetTitle("PbGl #sigma_{#gamma#gamma}");

  c0->Print("ToF_Calib.pdf");
}
