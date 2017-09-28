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
  {
    THnSparse *hn_minv = (THnSparse*)f->Get("hn_pion");
    TAxis *axis_pt = hn_minv->GetAxis(0);
    if(part == 0)
      hn_minv->GetAxis(2)->SetRange(1,6);
    else if(part == 1)
      hn_minv->GetAxis(2)->SetRange(7,8);
  }
  else if(data == 1)
  {
    THnSparse *hn_minv = (THnSparse*)f->Get("hn_2photon");
    TAxis *axis_pt = hn_minv->GetAxis(1);
    if(part == 0)
      hn_minv->GetAxis(0)->SetRange(1,6);
    else if(part == 1)
      hn_minv->GetAxis(0)->SetRange(7,8);
  }

  TCanvas *c = new TCanvas("c", "Canvas", 3600, 3000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  c->Divide(6,5);

  Int_t ipad = 1;
  for(Int_t ipt=0; ipt<gn; ipt++)
  {
    c->cd(ipad++);

    axis_pt->SetRange(ipt+1, ipt+1);
    if(data == 0)
      TH1 *h_minv = hn_minv->Projection(1);
    else if(data == 1)
      TH1 *h_minv = hn_minv->Projection(3);
    Double_t low = axis_pt->GetBinLowEdge(ipt+1);
    Double_t high = axis_pt->GetBinUpEdge(ipt+1);
    h_minv->SetTitle(Form("pT: %4.2f-%4.2f",low, high));

    TF1 *fn1 = new TF1("fn1", "gaus", 0., 0.5);
    TF1 *fn2 = new TF1("fn2", "gaus(0)+pol3(3)", 0., 0.5);

    h_minv->Fit(fn1, "Q0", "", 0.112, 0.162);
    fn2->SetParameters( fn1->GetParameters() );
    h_minv->Fit(fn2, "Q0", "", 0.047, 0.227);
    h_minv->Fit(fn2, "QE", "", 0.047, 0.227);

    gx[0][ipt] = axis_pt->GetBinCenter(ipt+1);
    gy[0][ipt] = fn2->GetParameter(1);
    egy[0][ipt] = fn2->GetParError(1);

    gx[1][ipt] = axis_pt->GetBinCenter(ipt+1);
    gy[1][ipt] = fn2->GetParameter(2);
    egy[1][ipt] = fn2->GetParError(2);
  }

  c->Print(Form("Pi0Peak-part%d-data%d.pdf",part,data));
  delete c;

  TGraphErrors **graph = new TGraphErrors*[2];
  for(Int_t i=0; i<2; i++)
    graph[i] = new TGraphErrors(gn, gx[i], gy[i], 0, egy[i]);
  return graph;
}

void draw_Pi0Peak()
{
  TFile *f_sim = new TFile("MissingRatio-histo.root");
  TFile *f_data = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ertb-cv/total.root");

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
        gr_sim[part][i]->GetYaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
        gr_sim[part][i]->GetYaxis()->SetRangeUser(0.13, 0.145);
      }
      else if(i == 1)
      {
        gr_sim[part][i]->GetYaxis()->SetTitle("#sigma_{#gamma#gamma} [GeV]");
        gr_sim[part][i]->GetYaxis()->SetRangeUser(0., 0.015);
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

  c0->Print("Pi0Peak.pdf");
}
