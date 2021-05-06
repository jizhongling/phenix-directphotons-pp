#include "GlobalVars.h"

void draw_ToFEff_Photon()
{
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};
  const char *name[2] = {"PbSc", "PbGl"};

  TFile *f_twr = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-DC3sigma.root");
  THnSparse *hn_etwr = (THnSparse*)f_twr->Get("hn_etwr");
  hn_etwr->GetAxis(7)->SetRange(1,1);
  hn_etwr->GetAxis(2)->SetRange(4,4);
  hn_etwr->GetAxis(3)->SetRange(4,4);
  hn_etwr->GetAxis(1)->SetRange(26,30);

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");
  //TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/PhotonEff-histo.root");

  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  TAxis *axis_sec = hn_1photon->GetAxis(0);
  TAxis *axis_pt = hn_1photon->GetAxis(1);
  TAxis *axis_pattern = hn_1photon->GetAxis(2);
  TAxis *axis_cut = hn_1photon->GetAxis(3);
  TAxis *axis_type = hn_1photon->GetAxis(4);

  TGraphAsymmErrors *gr[2];
  for(int part=0; part<2; part++)
    gr[part] =  new TGraphAsymmErrors(npT);

  mc(1, 2,1);

  for(int part=0; part<2; part++)
  {
    hn_etwr->GetAxis(0)->SetRange(secl[part],sech[part]);
    TH2 *h2_ct = hn_etwr->Projection(6,5);
    TH2 *h2_tof = hn_etwr->Projection(6,4);
    TH1 *h_tof = h2_tof->ProjectionX("h_tof");
    h_tof->Reset();

    double ngap_bp = 0.;
    double ngap_bpt = 0.;
    for(int ic=2; ic<=8; ic+=2)
    {
      for(int it=1; it<=16; it++)
      {
        ngap_bp += h2_ct->GetBinContent(it,ic);
        if( (ic==4 || ic==8) )
          ngap_bpt += h2_ct->GetBinContent(it,ic);
      }

      TH1 *h_tmp = h2_tof->ProjectionX("h_tmp",ic,ic);
      h_tof->Add(h_tmp);
      delete h_tmp;
    }

    mcd(1, part+1);
    aset(h_tof);
    style(h_tof, 20, 2);
    h_tof->DrawCopy();

    axis_type->SetRange(2,2);
    axis_sec->SetRange(secl[part],sech[part]);
    axis_cut->SetRange(3,3);
    TH1 *h_total = hn_1photon->Projection(1);
    axis_cut->SetRange(4,4);
    TH1 *h_passed = hn_1photon->Projection(1);
    gr[part]->Divide(h_passed, h_total);

    double npassed = h_passed->Integral(26,30);
    double ntotal = h_total->Integral(26,30);
    cout << ( npassed - ngap_bpt*120./9. ) / ( ntotal - ngap_bp*120./9. ) << endl;
    cout << ngap_bpt*120./9. << "\t" << ngap_bp*120./9. << "\t" << ntotal - npassed << endl; 
    cout << npassed << "\t" << ntotal << endl; 
    cout << h2_ct->GetEntries() << "\t" << h2_tof->GetEntries() << endl; 
    delete h2_ct;
    delete h2_tof;
    delete h_passed;
    delete h_total;
  }

  mc(0, 2,1);

  for(int part=0; part<2; part++)
  {
    mcd(0, part+1);
    gr[part]->SetTitle( Form("ToF efficeincy for %s", name[part]) );
    aset(gr[part], "p_{T} (GeV/c)", "Eff", 2.,30., 0.,1.1);
    style(gr[part], 24, kRed);
    gr[part]->Draw("APE");
    //gr[part]->Fit("pol0", "Q","", 9.,30.);

    //gPad->Update();
    //TPaveStats *st = (TPaveStats*)gr[part]->FindObject("stats");
    //st->SetY1NDC(0.4);
    //st->SetY2NDC(0.6);
  }

  c0->Print("plots/ToFEff-photon.pdf");
  c1->Print("plots/ToF-gap.pdf");
}
