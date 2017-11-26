void TowerLocation(UInt_t towerid, Int_t &sector, Int_t &ytower, Int_t &ztower)
{
  Int_t itower = 0;

  if( towerid < 15552 )
  {// PbSc
    sector = towerid / 2592;
    itower = towerid % 2592;
    ztower = itower % 72;
    ytower = itower / 72;
  }
  else
  {// PbGl
    sector = 6 + ( towerid - 15552 ) / 4608;
    itower = ( towerid - 15552 ) % 4608;
    ztower = itower % 96;
    ytower = itower / 96;
  }

  return;
}

void draw_TowerHits()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/total.root");
  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  TH1 *h_tower = hn_1photon->Projection(5);

  TH2 *h2_map[8];
  for(Int_t i=0; i<=5; i++)
    h2_map[i] =  new TH2F(Form("h2_map%d",i), Form("Sector %d;iypos;izpos;",i), 36, -0.5, 35.5, 72, -0.5, 71.5);
  for(Int_t i=6; i<=7; i++)
    h2_map[i] =  new TH2F(Form("h2_map%d",i), Form("Sector %d;iypos;izpos;",i), 48, -0.5, 47.5, 96, -0.5, 95.5);

  TH1 *h_hits = new TH1F("h_hits", "Hits distribution;nhits;ntowers;", 1000, 0., 100000.);

  for(Int_t it=0; it<24768; it++)
  {
    Int_t sector, iypos, izpos;
    TowerLocation(it, sector, iypos, izpos);
    Double_t nhits = h_tower->GetBinContent(it+1);
    h2_map[sector]->Fill(iypos, izpos, nhits);
    if(nhits > 100000.)
    {
      cout << sector << " " << iypos-1 << " " << izpos << " 40" << endl;
      cout << sector << " " << iypos << " " << izpos-1 << " 40" << endl;
      cout << sector << " " << iypos << " " << izpos << " 50" << endl;
      cout << sector << " " << iypos << " " << izpos+1 << " 40" << endl;
      cout << sector << " " << iypos+1 << " " << izpos << " 40" << endl;
    }
    if(nhits > 0.)
      h_hits->Fill(nhits);
  }

  TCanvas *c = new TCanvas("c", "Canvas", 2400, 1200);
  gStyle->SetOptStat(0);
  c->Divide(4,2);

  for(Int_t i=0; i<8; i++)
  {
    c->cd(i+1);
    h2_map[i]->Draw("colz");
  }

  c->Print("plots/TowerHits.pdf");

  TCanvas *c0 = new TCanvas("c0", "Canvas", 600, 600);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  c0->cd();
  h_hits->Fit("gaus", "Q", "", 0., 10000.);

  c0->Print("plots/HitsDistribution.pdf");
}
