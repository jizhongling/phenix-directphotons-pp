#include "BBCCounts.h"

void draw_Ratio()
{
  Double_t Nnovtx = 0.;
  Double_t N10cm = 0.;
  Double_t Nlive = 0.;
  Double_t Nert_c = 0.;
  TH1 *h_r10cm = new TH1F("h_r10cm", "BBC10cm/BBCNoVTX", 100, 0.07, 0.27);
  TH1 *h_rej = new TH1F("h_rej", "Rejection power", 1600, 0.5, 1600.5);

  Int_t thread = -1;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  while( fin >> runnumber )
  {
    thread++;
    if( thread%10 == 0 ) cout << "Nfile = " << thread << endl;

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-MB/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events= (TH1*)f->Get("h_events");

    Double_t Nbbc_novtx = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_novtx") );
    Double_t Nbbc_10cm = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_10cm_novtx") );
    Nnovtx += Nbbc_novtx;
    N10cm += Nbbc_10cm;
    Double_t r10cm = Nbbc_10cm / Nbbc_novtx;
    h_r10cm->Fill(r10cm);

    Double_t Nbbc_live = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_live") );
    Double_t Nbbc_ert_c = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_ert_c") );
    Nlive += Nbbc_live;
    Nert_c += Nbbc_ert_c;
    Double_t rej = Nbbc_live / Nbbc_ert_c;
    h_rej->Fill(rej);

    delete f;
  }

  cout << "N10cm/Nnovtx = " << N10cm/Nnovtx << endl
    << "r10cm mean = " << h_r10cm->GetMean() << endl
    << "r10cm error = " << h_r10cm->GetMeanError() << endl
    << "Nlive/Nert_c = " << Nlive/Nert_c << endl
    << "rej mean = " << h_rej->GetMean() << endl
    << "rej error = " << h_rej->GetMeanError() << endl;

  mc(0, 2,1);

  TF1 *fn_gaus = new TF1("fn_gaus", "gaus");

  mcd(0, 1);
  fn_gaus->SetParameters(h_r10cm->GetMaximum(), 0.16, 0.015);
  h_r10cm->Fit(fn_gaus, "Q", "", 0.07, 0.27);

  mcd(0, 2);
  h_rej->Rebin(10);
  fn_gaus->SetParameters(h_rej->GetMaximum(), 800., 100.);
  h_rej->Fit(fn_gaus, "Q", "", 0.5, 1600.5);

  c0->Print("NBBC-r10cm-rej.pdf");
}

void draw_NBBC()
{
  Double_t Ndb = 0.;
  Double_t Nrej = 0.;
  TGraph *gr_ert = new TGraph(1000);

  Int_t thread = -1;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  ReadClockCounts();

  while( fin >> runnumber )
  {
    thread++;
    if( thread%10 == 0 ) cout << "Nfile = " << thread << endl;

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events = (TH1*)f->Get("h_events");

    const Double_t r10cm = 0.164395;
    const Double_t rej = 824.2;
    Double_t NBBC_ert_db = r10cm * (Double_t)( GetBBCNovtxLive(runnumber) / (GetERT4x4cScaledown(runnumber)+1) );
    Double_t NBBC_ert_rej = rej * h_events->GetBinContent( h_events->GetXaxis()->FindBin("ert_c") );

    if( NBBC_ert_db > 0. && NBBC_ert_rej > 0. )
    {
      Ndb += NBBC_ert_db;
      Nrej += NBBC_ert_rej;
      gr_ert->SetPoint(thread, runnumber, NBBC_ert_db/NBBC_ert_rej);
    }
    else
    {
      cout << "Problem in Run " << runnumber << endl;
    }

    delete f;
  }

  cout << "Ndb = " <<  Ndb << endl
    << "Nrej = " <<  Nrej << endl
    << "Ndb/Nrej = " <<  Ndb/Nrej << endl;

  mc();
  mcd();

  aset(gr_ert, "Runnumber","#frac{DB}{Rej}", 386700.,398200., 0.,3.);
  style(gr_ert, 20, 1);
  gr_ert->Draw("AP");

  c0->Print("NBBC-DBRej.pdf");
}
