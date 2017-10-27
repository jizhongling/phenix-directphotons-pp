#include "BBCCounts.h"

void draw_Ratio()
{
  TH1 *h_ratio = new TH1F("h_ratio", "BBC10cm/BBCNoVTX", 40,0.54,0.62);

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

    Double Nbbc_novtx = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_novtx") );
    Double Nbbc_10cm = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_10cm") );
    Double_t ratio = Nbbc_10cm / Nbbc_novtx;
    h_ratio->Fill(ratio);

    delete f;
  }

  mc();
  mcd();

  TF1 *fn_gaus = new TF1("fn_gaus", "gaus", 0.54, 0.62);
  fn_gaus->SetParameters(h_ratio->GetMaximum(), 0.58, 0.01);
  h_ratio->Fit(fn_gaus, "RQ");

  cout << "Mean = " << h_ratio->GetMean() << endl
    << "Error = " << h_ratio->GetMeanError() << endl;

  c0->Print("NBBC-r10cm.pdf");
}

void draw_NBBC()
{
  Double_t Ndb[3];
  Double_t Nrej[3];
  TGraph *gr_ert[3];
  for(Int_t igr=0; igr<3; igr++)
  {
    Ndb[igr] = 0.;
    Nrej[igr] = 0.;
    gr_ert[igr] = new TGraph(1000);
  }


  Int_t thread = -1;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  ReadClockCounts();

  while( fin >> runnumber )
  {
    thread++;
    if( thread%10 == 0 ) cout << "Nfile = " << thread << endl;

    TFile *f_ert = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    TFile *f_mb = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-MB/PhotonNode-%d.root",runnumber));
    if( f_ert->IsZombie() || f_mb->IsZombie() ) continue;

    TH1 *h_events_ert = (TH1*)f_ert->Get("h_events");
    TH1 *h_events_mb = (TH1*)f_mb->Get("h_events");

    const Double_t ratio = 0.577109;
    const Int_t bin_ert = h_events_ert->GetXaxis()->FindBin("ert_a");
    const Int_t bin_mb = h_events_mb->GetXaxis()->FindBin("bbc");
    Double_t NBBC_ert_db[3];
    Double_t NBBC_ert_rej[3];
    NBBC_ert_db[0] = ratio * (Double_t)( GetBBCNarrowLive(runnumber) / (GetERT4x4aScaledown(runnumber)+1) );
    NBBC_ert_db[1] = ratio * (Double_t)( GetBBCNarrowLive(runnumber) / (GetERT4x4bScaledown(runnumber)+1) );
    NBBC_ert_db[2] = ratio * (Double_t)( GetBBCNarrowLive(runnumber) / (GetERT4x4cScaledown(runnumber)+1) );
    NBBC_ert_rej[0] = h_events_ert->GetBinContent(bin_ert)   * h_events_mb->GetBinContent(bin_mb) / h_events_mb->GetBinContent(bin_mb+1);
    NBBC_ert_rej[1] = h_events_ert->GetBinContent(bin_ert+1) * h_events_mb->GetBinContent(bin_mb) / h_events_mb->GetBinContent(bin_mb+2);
    NBBC_ert_rej[2] = h_events_ert->GetBinContent(bin_ert+2) * h_events_mb->GetBinContent(bin_mb) / h_events_mb->GetBinContent(bin_mb+3);

    for(Int_t igr=0; igr<3; igr++)
    {
      if( NBBC_ert_db[igr] < TMath::Infinity() && NBBC_ert_rej[igr] < TMath::Infinity() )
      {
        Ndb[igr] += NBBC_ert_db[igr];
        Nrej[igr] += NBBC_ert_rej[igr];
        gr_ert[igr]->SetPoint(thread, runnumber, NBBC_ert_db[igr]/NBBC_ert_rej[igr]);
      }
      else
      {
        cout << "Problem in Run " << runnumber << Form(" for ERT_4x4%c",97+igr) << endl;
      }
    }

    delete f_ert;
    delete f_mb;
  }

  mc();
  mcd();
  legi();

  for(Int_t igr=0; igr<3; igr++)
  {
    cout << Form("Ndb/Nrej for ERT_4x4%c: ",97+igr) <<  Ndb[igr]/Nrej[igr] << endl;
    aset(gr_ert[igr], "Runnumber","#frac{DB}{Rej}", 387000,399000, 0.5,1.5);
    style(gr_ert[igr], 20+igr, 1+igr);
    if(igr==0)
      gr_ert[igr]->Draw("AP");
    else
      gr_ert[igr]->Draw("P");
    leg0->AddEntry(gr_ert[igr], Form("ERT_4x4%c",97+igr), "P");
  }
  leg0->Draw();

  c0->Print("NBBC-DBRej.pdf");
}
