void draw_DCCheck()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h2/h3[ns][we]
  TH3 *h3_dphiz[2][2];
  TH2 *h2_board[2][2];

  TH3 *h3_dphiz_t = (TH3*)f->Get("h3_dcdphiz_0");
  TH2 *h2_board_t = (TH2*)f->Get("h2_alphaboard_0");
  TH3 *h3_live = (TH3*)f->Get("h3_dclive_0");
  h3_dphiz_t = (TH3*)h3_dphiz_t->Clone();
  h2_board_t = (TH2*)h2_board_t->Clone();
  h3_live = (TH3*)h3_live->Clone();
  h3_dphiz_t->Reset();
  h2_board_t->Reset();
  h3_live->Reset();

  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      h3_dphiz[ns][we] = (TH3*)h3_dphiz_t->Clone( Form("h3_dphiz_ns%d_we_%d",ns,we) );
      h2_board[ns][we] = (TH2*)h2_board_t->Clone( Form("h2_board_ns%d_we_%d",ns,we) );
      for(int qual=4; qual<64; qual++)
      {
        int ih = ns + 2*we + 2*2*qual;
        TH3 *h3_dphiz_tmp = (TH3*)f->Get( Form("h3_dcdphiz_%d",ih) );
        TH2 *h2_board_tmp = (TH2*)f->Get( Form("h2_alphaboard_%d",ih) );
        h3_dphiz[ns][we]->Add(h3_dphiz_tmp);
        h2_board[ns][we]->Add(h2_board_tmp);
        delete h3_dphiz_tmp;
        delete h2_board_tmp;
      }
    }

  for(int qual=4; qual<64; qual++)
  {
    int ih = qual;
    TH3 *h3_tmp = (TH3*)f->Get( Form("h3_dclive_%d",ih) );
    h3_live->Add(h3_tmp);
    delete h3_tmp;
  }

  mc(0, 2,2);
  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      mcd(0, 2*ns+we+1);
      TString NS = ns ? "S" : "N";
      TString WE = we ? "E" : "W";
      h2_board[ns][we]->SetTitle(NS+WE);
      h2_board[ns][we]->DrawCopy("COLZ");
    }

  mc(1);
  mcd(1);
  h3_live->GetZaxis()->SetRange(3,-1);
  h3_live->Project3D("yx")->Draw("COLZ");

  TF1 *f_fit = new TF1("f_fit", "gaus(0)+pol2(3)", -100., 100.);
  for(int id=0; id<2; id++)
  {
    mc(id+2, 2,2);
    for(int ns=0; ns<2; ns++)
      for(int we=0; we<2; we++)
      {
        TString ax = id ? "y" : "x";
        TString NS = ns ? "S" : "N";
        TString WE = we ? "E" : "W";

        mcd(id+2, 2*ns+we+1);
        h3_dphiz[ns][we]->GetZaxis()->SetRange(3,14);
        TH1 *h_diff = h3_dphiz[ns][we]->Project3D(ax);
        if(id==0)
        {
          double par[10] = {6e6,0.,0.005, 0.,0.,0.};
          f_fit->SetParameters(par);
          h_diff->Fit(f_fit, "Q0","", -0.025,0.025);
          h_diff->GetXaxis()->SetRangeUser(-0.015,0.015);
        }
        else
        {
          double par[10] = {1.5e6,0.,5., 0.,0.,0.};
          f_fit->SetParameters(par);
          h_diff->Fit(f_fit, "Q0","", -15.,15.);
          h_diff->GetXaxis()->SetRangeUser(-10,10.);
        }

        cout << NS+WE+" RMS: " << h_diff->GetRMS() << endl;
        h_diff->GetXaxis()->SetRange(0,-1);
        h_diff->SetTitle(NS+WE);
        aset(h_diff);
        f_fit->SetLineColor(kRed);
        h_diff->DrawCopy();
        f_fit->DrawCopy("SAME");
      }
  }

  TMultiGraph *mg[2][2];  // mg[ns][we]
  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
      mg[ns][we] = new TMultiGraph();

  for(int i=0; i<44; i++)
  {
    TFile *f = new TFile(Form("histos/DCCheck-%d.root",i));
    if( f->IsZombie() ) continue;
    for(int ns=0; ns<2; ns++)
      for(int we=0; we<2; we++)
      {
        TGraph *gr = (TGraph*)f->Get( Form("gr_yield_ns%d_we%d",ns,we) );
        gr->SetMarkerStyle(20);
        if( gr->GetN() > 0 )
          mg[ns][we]->Add(gr);
      }
  }

  mc(4, 2,2);
  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      mcd(4, 2*ns+we+1);
      TString NS = ns ? "S" : "N";
      TString WE = we ? "E" : "W";
      mg[ns][we]->Draw("AP");  // must before GetXaxis()
      mg[ns][we]->SetTitle(NS+WE);
      mg[ns][we]->GetXaxis()->SetTitle("runnumber");
      mg[ns][we]->GetYaxis()->SetTitle("Ntrack/Nevent");
      //mg[ns][we]->GetXaxis()->SetLimits(387000., 398200.);  // Do not use SetRangeUser()
      mg[ns][we]->GetYaxis()->SetRangeUser(0., 1.5);  // Do not use SetLimits()
    }

  mc(5, 2,1);
  TFile *f_combine = new TFile("histos/DCCheck.root");
  TH2 *h2_phi[2];  // h2_phi[ns]
  for(int ns=0; ns<2; ns++)
  {
    mcd(5, ns+1);
    TString name = ns ? "h2_phi_zm" : "h2_phi_zp";
    TString title = ns ? "DC phi vs run in South" : "DC phi vs run in North";
    h2_phi[ns] = (TH2*)f_combine->Get(name);
    h2_phi[ns]->SetTitle(title);
    h2_phi[ns]->Draw("COLZ");
  }
}
