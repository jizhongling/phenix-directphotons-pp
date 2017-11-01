#include "BBCCounts.h"

const char *dname[3] = {"NoVTX", "NarrowVTX", "MB"};

void draw_Ratio()
{
  // NoVTX, NarrowVTX, MB
  const Double_t mean_r_bbc[2] = {1.00209e+06, 575056.};
  const Double_t mean_rej_bbc[3] = {164739., 309495., 589550.};

  Double_t sum_r_bbc[2] = {}, sum_r_10cm[2] = {};
  Double_t sum_rej_bbc[3] = {}, sum_rej_ert_c[3];

  TH1::SetDefaultSumw2();

  TH1 *h_r10cm[2];
  h_r10cm[0] = new TH1F("h_r10cm_0", "BBC10cm/BBCNoVTX", 100, 0.07, 0.27);
  h_r10cm[1] = new TH1F("h_r10cm_1", "BBC10cm/BBCNarrow", 100, 0.45, 0.65);
  
  TH1 *h_rej[3];
  for(Int_t id=0; id<3; id++)
    h_rej[id] = new TH1F(Form("h_rej_%d",id), Form("Rejection power from %s",dname[id]), 1600, 0.5, 1600.5);

  Int_t thread = -1;
  Int_t irun = 0;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  while( fin >> runnumber )
  {
    thread++;
    if( thread%10 == 0 ) cout << "Nfile = " << thread << endl;

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-MB/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events= (TH1*)f->Get("h_events");

    Double_t N_r_bbc[2], N_r_10cm[2];
    N_r_bbc[0] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_novtx") );
    N_r_10cm[0] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_novtx_narrow_10cm") );
    N_r_bbc[1] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_narrow") );
    N_r_10cm[1] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_narrow_10cm") );

    for(Int_t id=0; id<2; id++)
    {
      sum_r_bbc[id] += N_r_bbc[id];
      sum_r_10cm[id] += N_r_10cm[id];
      Double_t r10cm = N_r_10cm[id] / N_r_bbc[id];
      h_r10cm[id]->Fill(r10cm, N_r_bbc[id]/mean_r_bbc[id]);
    }

    Double_t N_rej_bbc[3], N_rej_ert_c[3];
    N_rej_bbc[0] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_novtx_narrow_10cm") );
    N_rej_ert_c[0] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_novtx_narrow_10cm_ert_c") );
    N_rej_bbc[1] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_narrow_10cm") );
    N_rej_ert_c[1] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_narrow_10cm_ert_c") );
    N_rej_bbc[2] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_mb_narrow_10cm") );
    N_rej_ert_c[2] = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_mb_narrow_10cm_ert_c") );

    for(Int_t id=0; id<3; id++)
    {
      sum_rej_bbc[id] += N_rej_bbc[id];
      sum_rej_ert_c[id] += N_rej_ert_c[id];
      Double_t rej = N_rej_bbc[id] / N_rej_ert_c[id];
      h_rej[id]->Fill(rej, N_rej_bbc[id]/mean_rej_bbc[id]);
    }

    delete f;
    irun++;
  }

  for(Int_t id=0; id<2; id++)
    cout << "r10cm in " << dname[id] << endl
      << "mean_r_bbc = " << sum_r_bbc[id]/irun << endl
      << "sum_10cm/sum_bbc = " << sum_r_10cm[id]/sum_r_bbc[id] << endl
      << "r10cm mean = " << h_r10cm[id]->GetMean() << endl
      << "r10cm error = " << h_r10cm[id]->GetMeanError() << endl << endl;

  for(Int_t id=0; id<3; id++)
    cout << "Rejection power in " << dname[id] << endl
      << "mean_rej_bbc = " << sum_rej_bbc[id]/irun << endl
      << "sum_bbc/sum_ert_c = " << sum_rej_bbc[id]/sum_rej_ert_c[id] << endl
      << "rej mean = " << h_rej[id]->GetMean() << endl
      << "rej error = " << h_rej[id]->GetMeanError() << endl << endl;

  TF1 *fn_gaus = new TF1("fn_gaus", "gaus");
  mc(0, 2,3);

  mcd(0, 1);
  fn_gaus->SetParameters(h_r10cm[0]->GetMaximum(), 0.17, 0.015);
  h_r10cm[0]->Fit(fn_gaus, "Q", "", 0.07, 0.27);

  mcd(0, 2);
  fn_gaus->SetParameters(h_r10cm[0]->GetMaximum(), 0.55, 0.020);
  h_r10cm[1]->Fit(fn_gaus, "Q", "", 0.45, 0.65);

  for(Int_t id=0; id<3; id++)
  {
    mcd(0, id+3);
    h_rej[id]->Rebin(10);
    fn_gaus->SetParameters(h_rej[id]->GetMaximum(), 800., 100.);
    h_rej[id]->Fit(fn_gaus, "Q", "", 0.5, 1600.5);
  }

  c0->Print("NBBC-r10cm-rej.pdf");
}

void draw_NBBC()
{
  // NoVTX, NarrowVTX
  const Double_t r10cm[2] = {0.1644, 0.5389};
  const Double_t er10cm[2] = {0.00065, 0.00097};
  const Double_t rej = 826.7;
  const Double_t erej = 4.1;

  Double_t sum_db[2] = {};
  Double_t sum_rej[2] = {};

  TGraphErrors *gr_rlum[2];
  Int_t igp[2] = {};
  for(Int_t id=0; id<2; id++)
    gr_rlum[id] = new TGraphErrors(1000);

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

    Double_t N_scaled = h_events->GetBinContent( h_events->GetXaxis()->FindBin("ert_c") );
    Double_t N_rej = rej * N_scaled;
    Double_t re2N_rej = pow(erej/rej,2.) + 1./N_scaled;

    const ULong64_t N_live[2] = {GetBBCNovtxLive(runnumber), GetBBCNarrowLive(runnumber)};
    const ULong_t scaledown = GetERT4x4cScaledown(runnumber) + 1;
    Double_t N_db[2], re2N_db[2];

    for(Int_t id=0; id<2; id++)
    {
      N_db[id] = r10cm[id] * (Double_t)(N_live[id]/scaledown);
      re2N_db[id] = pow(er10cm[id]/r10cm[id],2.) + 1./N_live[id];

      sum_db[id] += N_db[id];
      sum_rej[id] += N_rej[id];

      Double_t yy = N_db[id] / N_rej;
      Double_t eyy = yy * sqrt( re2N_db[id] + re2N_rej );
      if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
      {
        gr_rlum[id]->SetPoint(igp[id], runnumber, yy);
        gr_rlum[id]->SetPointError(igp[id], 0., eyy);
        igp[id]++;
      }
      else
      {
        cout << "Problem in Run " << runnumber << endl;
      }
    }

    delete f;
  }

  mc(0, 2,1);

  for(Int_t id=0; id<2; id++)
  {
    gr_rlum[id]->Set(igp[id]);

    cout << dname[id] << endl
      << "number of runs = " << gr_rlum[id]->GetN() << endl
      << "mean_rlum = " << gr_rlum[id]->GetMean(2) << endl
      << "emean_rlum = " << gr_rlum[id]->GetRMS(2) << endl
      << "sum_db = " <<  sum_db[id] << endl
      << "sum_rej = " <<  sum_rej[id] << endl
      << "sum_db/sum_rej = " <<  sum_db[id]/sum_rej[id] << endl << endl;

    mcd(0, id+1);
    gr_rlum[id]->SetTitle( Form("r_{10cm} from %s scaled events",dname[id]) );
    aset(gr_rlum[id], "Runnumber","#frac{L_{db}}{L_{rej}}", 386700.,398200., 0.,3.);
    style(gr_rlum[id], 20, 1);
    gr_rlum[id]->Draw("AP");
  }

  c0->Print("NBBC-DBRej.pdf");
}
