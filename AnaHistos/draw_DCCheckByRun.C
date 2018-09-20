#include "BBCCounts.h"

void draw_DCCheckByRun()
{
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");
  ReadClockCounts();

  TGraph *gr_yield = new TGraph(1000);
  //TH2 *h_phi[2];
  //h2_phi[0] = new TH2F("h2_phi_zp", "DC phi vs run;runnumber;#phi [rad];", 1000,0.,1000., 50,-1.,4.);
  //h2_phi[1] = (TH2*)h2_phi[0]->Clone("h2_phi_zm");

  while( fin >> runnumber )
  {
    thread++;
    if(thread >2) break;
    if( thread%10 == 0 ) cout << "NFiles = " << thread << endl;

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/13888/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH3 *h3_live = (TH3*)f->Get("h3_dclive_0");
    h3_live = (TH3*)h3_live->Clone();
    h3_live->Reset();

    for(int qual=4; qual<64; qual++)
    {
      int ih = qual;
      TH3 *h3_tmp = (TH3*)f->Get( Form("h3_dclive_%d",ih) );
      h3_live->Add(h3_tmp);
      delete h3_tmp;
    }

    ULong64_t nmb = GetBBCNarrowLive(runnumber);
    ULong_t scaledown = GetERT4x4cScaledown(runnumber) + 1;
    double nev = nmb / scaledown;

    double nyield = h3_live->GetEntries();
    gr_yield->SetPoint(thread, (double)runnumber, nyield/nev);
  }

  mc();
  mcd();
  gr_yield->Set(thread+1);
  aset(gr_yield, "runnumber","yield/nmb", 11200,38700.,398200.);
  style(gr_yield, 20, 1);
  gr_yield->Draw("AP");
}
