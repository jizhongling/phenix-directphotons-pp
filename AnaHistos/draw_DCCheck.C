void draw_DCCheck()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h2[ns][we]
  TH2 *h2_sd[2][2];
  TH2 *h2_board[2][2];

  TH2 *h2_sd_t = (TH2*)f->Get("h2_dcsd_0");
  TH2 *h2_board_t = (TH2*)f->Get("h2_alphaboard_0");
  TH3 *h3_live = (TH3*)f->Get("h3_dclive_0");
  h2_sd_t = (TH2*)h2_sd_t->Clone();
  h2_board_t = (TH2*)h2_board_t->Clone();
  h3_live = (TH3*)h3_live->Clone();
  h2_sd_t->Reset();
  h2_board_t->Reset();
  h3_live->Reset();

  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      h2_sd[ns][we] = (TH2*)h2_sd_t->Clone( Form("h2_sd_ns%d_we_%d",ns,we) );
      h2_board[ns][we] = (TH2*)h2_board_t->Clone( Form("h2_board_ns%d_we_%d",ns,we) );
      for(int qual=4; qual<64; qual++)
      {
        int ih = ns + 2*we + 2*2*qual;
        TH2 *h2_sd_tmp = (TH2*)f->Get( Form("h2_dcsd_%d",ih) );
        TH2 *h2_board_tmp = (TH2*)f->Get( Form("h2_alphaboard_%d",ih) );
        h2_sd[ns][we]->Add(h2_sd_tmp);
        h2_board[ns][we]->Add(h2_board_tmp);
        delete h2_sd_tmp;
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
      mcd(0, ns+2*we+1);
      TString NS = ns ? "S" : "N";
      TString WE = we ? "E" : "W";
      h2_board[ns][we]->SetTitle(NS+WE);
      h2_board[ns][we]->DrawCopy("COLZ");
    }

  mc(1);
  mcd(1);

  h3_live->GetZaxis()->SetRange(3,-1);
  h3_live->Project3D("yx")->Draw();
}
