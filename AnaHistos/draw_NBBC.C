#include "DataBase.h"

const char *dname[2] = {"NarrowVTX", "NoVTX"};

void draw_BBCCounts()
{
  const char *fname[2] = {"runlist-Sasha", "runlist-DC3sigma"};
  DataBase *db = new DataBase();

  for(int id=0; id<2; id++)
  {
    ULong64_t sum_db = 0;

    int runnumber;
    ifstream fin( Form("/phenix/plhf/zji/taxi/Run13pp510MinBias/%s.txt",fname[id]) );

    while( fin >> runnumber )
    {
      ULong64_t N_bbc_live = db->GetBBCNarrowLive(runnumber) / ( db->GetERT4x4cScaledown(runnumber) + 1 );
      sum_db += N_bbc_live;
    }

    cout << "BBCCounts for " << fname[id] << ": " << sum_db << endl;
  }
}

void draw_Ratio()
{
  // NoVTX, NarrowVTX, MB
  const double mean_r_bbc[2] = {575056., 1.00209e+06};
  const double mean_rej_bbc[2] = {309495., 164739.};

  double sum_r_bbc[2] = {}, sum_r_10cm[2] = {};
  double sum_rej_bbc[2] = {}, sum_rej_ert_c[2] = {};

  TH1::SetDefaultSumw2();

  TH1 *h_r10cm[2];
  h_r10cm[0] = new TH1F("h_r10cm_0", "BBC10cm/BBCNarrow", 100, 0.45, 0.65);
  h_r10cm[1] = new TH1F("h_r10cm_1", "BBC10cm/BBCNoVTX", 100, 0.07, 0.27);

  TH1 *h_rej[2];
  for(int id=0; id<2; id++)
    h_rej[id] = new TH1F(Form("h_rej_%d",id), Form("Rejection power from %s",dname[id]), 1600, 0.5, 1600.5);

  int thread = -1;
  int irun = 0;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  TTree *t1 = new TTree("t1", "Run group information");
  double r10cm[2], rej[2];
  t1->Branch("runnumber", &runnumber, "runnumber/I");
  t1->Branch("r10cm", r10cm, "r10cm[2]/D");
  t1->Branch("rej", rej, "rej[2]/D");

  while( fin >> runnumber )
  {
    thread++;
    if( thread%10 == 0 ) cout << "Nfiles = " << thread << endl;

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-MB/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events= (TH1*)f->Get("h_events");

    double N_r_bbc[2], N_r_10cm[2];
    N_r_bbc[0] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_narrow") );
    N_r_10cm[0] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_narrow_10cm") );
    N_r_bbc[1] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_novtx") );
    N_r_10cm[1] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_novtx_narrow_10cm") );

    for(int id=0; id<2; id++)
    {
      sum_r_bbc[id] += N_r_bbc[id];
      sum_r_10cm[id] += N_r_10cm[id];
      r10cm[id] = N_r_10cm[id] / N_r_bbc[id];
      h_r10cm[id]->Fill(r10cm[id], N_r_bbc[id]/mean_r_bbc[id]);
    }

    double N_rej_bbc[2], N_rej_ert_c[2];
    N_rej_bbc[0] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_narrow_10cm") );
    N_rej_ert_c[0] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_narrow_10cm_ert_c") );
    N_rej_bbc[1] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_novtx_narrow_10cm") );
    N_rej_ert_c[1] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_novtx_narrow_10cm_ert_c") );

    for(int id=0; id<2; id++)
    {
      sum_rej_bbc[id] += N_rej_bbc[id];
      sum_rej_ert_c[id] += N_rej_ert_c[id];
      rej[id] = N_rej_bbc[id] / N_rej_ert_c[id];
      h_rej[id]->Fill(rej[id], N_rej_bbc[id]/mean_rej_bbc[id]);
    }

    t1->Fill();

    delete f;
    irun++;
  }

  for(int id=0; id<2; id++)
    cout << "r10cm in " << dname[id] << endl
      << "mean_r_bbc = " << sum_r_bbc[id]/irun << endl
      << "sum_10cm/sum_bbc = " << sum_r_10cm[id]/sum_r_bbc[id] << endl
      << "r10cm mean = " << h_r10cm[id]->GetMean() << endl
      << "r10cm error = " << h_r10cm[id]->GetMeanError() << endl << endl;

  for(int id=0; id<2; id++)
    cout << "Rejection power in " << dname[id] << endl
      << "mean_rej_bbc = " << sum_rej_bbc[id]/irun << endl
      << "sum_bbc/sum_ert_c = " << sum_rej_bbc[id]/sum_rej_ert_c[id] << endl
      << "rej mean = " << h_rej[id]->GetMean() << endl
      << "rej error = " << h_rej[id]->GetMeanError() << endl << endl;

  TF1 *fn_gaus = new TF1("fn_gaus", "gaus");
  mc(0, 2,2);

  mcd(0, 1);
  fn_gaus->SetParameters(h_r10cm[0]->GetMaximum(), 0.55, 0.020);
  h_r10cm[0]->Fit(fn_gaus, "Q", "", 0.45, 0.65);

  mcd(0, 2);
  fn_gaus->SetParameters(h_r10cm[1]->GetMaximum(), 0.17, 0.015);
  h_r10cm[1]->Fit(fn_gaus, "Q", "", 0.07, 0.27);

  for(int id=0; id<2; id++)
  {
    mcd(0, id+3);
    h_rej[id]->Rebin(10);
    fn_gaus->SetParameters(h_rej[id]->GetMaximum(), 800., 100.);
    h_rej[id]->Fit(fn_gaus, "Q", "", 0.5, 1600.5);
  }

  c0->Print("plots/NBBC-r10cm-rej.pdf");

  TFile *f_out = new TFile("data/NBBC.root", "RECREATE");
  t1->Write();
  f_out->Close();
}

void draw_NBBC()
{
  double sum_db = 0.;
  double sum_rej = 0.;

  TGraphErrors *gr_rlum = new TGraphErrors(1000);
  int igp = 0;

  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist-DC3sigma.txt");

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    thread++;
    if( thread%10 == 0 ) cout << "Nfile = " << thread << endl;

    TFile *f_ert = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    TFile *f_mb = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-MB/PhotonNode-%d.root",runnumber));
    if( f_ert->IsZombie() || f_mb->IsZombie() ) continue;

    TH1 *h_events_ert = (TH1*)f_ert->Get("h_events");
    TH1 *h_events_mb = (TH1*)f_mb->Get("h_events");

    ULong64_t N_bbc_live = db->GetBBCNarrowLive(runnumber) / ( db->GetERT4x4cScaledown(runnumber) + 1 );
    double N_r_bbc = h_events_mb->GetBinContent( h_events_mb->GetXaxis()->FindBin("bbc_narrow") );
    double N_r_10cm = h_events_mb->GetBinContent( h_events_mb->GetXaxis()->FindBin("bbc_narrow_10cm") );
    double r10cm = N_r_10cm / N_r_bbc;
    //double N_db = (double)N_bbc_live * r10cm;
    //double re2N_db = 1./N_bbc_live + 1./N_r_bbc + 1./N_r_10cm;
    double N_db = N_r_10cm * ( db->GetBBCNarrowScaledown(runnumber) + 1 ) / ( db->GetERT4x4cScaledown(runnumber) + 1 );
    double re2N_db = 1./N_r_10cm;
    sum_db += N_db;

    double N_ert_scaled = h_events_ert->GetBinContent( h_events_ert->GetXaxis()->FindBin("ert_c") );
    double N_rej_bbc = h_events_mb->GetBinContent( h_events_mb->GetXaxis()->FindBin("bbc_narrow_10cm") );
    double N_rej_ert_c = h_events_mb->GetBinContent( h_events_mb->GetXaxis()->FindBin("bbc_narrow_10cm_ert_c") );
    double rej = N_rej_bbc / N_rej_ert_c;
    double N_rej = N_ert_scaled * rej;
    double re2N_rej = 1./N_ert_scaled + 1./N_rej_bbc + 1./N_rej_ert_c;
    sum_rej += N_rej;

    double yy = N_rej / N_db;
    double eyy = yy * sqrt( re2N_db + re2N_rej );
    if( TMath::Finite(yy+eyy) )
    {
      gr_rlum->SetPoint(igp, runnumber, yy);
      gr_rlum->SetPointError(igp, 0., eyy);
      igp++;
    }
    else
      cout << "Problem in Run " << runnumber << endl;

    delete f_ert;
    delete f_mb;
  }

  gr_rlum->Set(igp);

  cout << "number of runs = " << gr_rlum->GetN() << endl
    << "sum_db = " <<  sum_db << endl
    << "sum_rej = " <<  sum_rej << endl
    << "sum_rej/sum_db = " << sum_rej/sum_db << endl;

  mc();
  mcd();

  gr_rlum->SetTitle("Luminosity ratio");
  aset(gr_rlum, "Runnumber","#frac{L_{Rej}}{L_{DB}}", 386700.,398200., 0.6,1.8);
  style(gr_rlum, 20, 1);
  gr_rlum->Draw("AP");

  c0->Print("plots/NBBC-DBRej.pdf");
}
