#include "Pileup.h"
#include "BBCCounts.h"
#include "ReadGraph.h"

void anaPileup_Photon(const int process = 0)
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  const int nThread = 100;
  int thread = -1;
  int irun = 0;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  int id = 0;

  TGraphErrors *gr[npT*2*2*3];
  TGraphErrors *gr_run[npT*2*2*3];
  int igp[npT*2*2*3] = {};
  for(int ipt=0; ipt<npT; ipt++)
    for(int ic=0; ic<2; ic++)
      for(int part=0; part<3; part++)
      {
        int ig = part + 3*ic + 2*3*id + 2*2*3*ipt;
        gr[ig] = new TGraphErrors(nThread);
        gr_run[ig] = new TGraphErrors(nThread);
        gr[ig]->SetName(Form("gr_%d",ig));
        gr_run[ig]->SetName(Form("gr_run_%d",ig));
      }

  const double A = 0.22;
  const double eA = 0.04;
  double xMiss[3][30] = {}, Miss[3][30] = {}, eMiss[3][30] = {};
  double xMerge[3][30] = {}, Merge[3][30] = {}, eMerge[3][30] = {};
  double xBadPass[3][30] = {}, BadPass[3][30] = {}, eBadPass[3][30] = {};
  for(int part=0; part<3; part++)
  {
    ReadGraph<TGraphErrors>("data/MissingRatio.root", part, xMiss[part], Miss[part], eMiss[part]);
    ReadGraph<TGraphAsymmErrors>("data/Merge-photon.root", part, xMerge[part], Merge[part], eMerge[part]);
    ReadGraph<TGraphErrors>("data/MergePassRate.root", part/2, xBadPass[part], BadPass[part], eBadPass[part]);
  }

  ReadClockCounts();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/13664/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events = (TH1*)f->Get("h_events");

    // h[ic][part]
    TH1 *h_1photon[2][3];
    TH2 *h2_2photon[2][3];

    int bbc10cm = 1;
    int evtype = 2;

    TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
    h_1photon_t->Reset();
    for(int ic=0; ic<2; ic++)
      for(int part=0; part<3; part++)
      {
        h_1photon[ic][part] = (TH1*)h_1photon_t->Clone(Form("h_1photon_ic%d_part%d",ic,part));
        for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
          for(int pattern=0; pattern<3; pattern++)
            for(int isolated=0; isolated<2; isolated++)
            {
              int cut = ic + 2;
              int ih = sector + 8*pattern + 3*8*isolated + 2*3*8*cut + 4*2*3*8*evtype + 3*4*2*3*8*bbc10cm;
              TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
              h_1photon[ic][part]->Add(h_tmp);
              delete h_tmp;
            }
        if(ic==1)
          h_1photon[0][part]->Add(h_1photon[1][part]);
      }

    TH2 *h2_2photon_t = (TH2*)f->Get("h2_2photon_0");
    h2_2photon_t->Reset();
    for(int ic=0; ic<2; ic++)
      for(int part=0; part<3; part++)
      {
        h2_2photon[ic][part] = (TH2*)h2_2photon_t->Clone(Form("h2_2photon_ic%d_part%d",ic,part));
        for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
          for(int pattern=0; pattern<3; pattern++)
            for(int isoboth=0; isoboth<2; isoboth++)
              for(int isopair=0; isopair<2; isopair++)
              {
                int cut = ic + 2;
                int ih = sector + 8*pattern + 3*8*isoboth + 2*3*8*isopair + 2*2*3*8*cut + 4*2*2*3*8*evtype + 3*4*2*2*3*8*bbc10cm;
                TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
                h2_2photon[ic][part]->Add(h2_tmp);
                delete h2_tmp;
              }
        if(ic==1)
          h2_2photon[0][part]->Add(h2_2photon[1][part]);
      }

    ULong64_t nclock = GetClockLive(runnumber);
    ULong64_t nmb = GetBBCNarrowLive(runnumber);
    ULong_t scaledown = GetERT4x4cScaledown(runnumber) + 1;

    //double nev = h_events->GetBinContent( h_events->GetXaxis()->FindBin("ert_c_10cm") );
    double nev = nmb / scaledown;

    for(int ipt=0; ipt<npT; ipt++)
      for(int ic=0; ic<2; ic++)
        for(int part=0; part<3; part++)
        {
          int ig = part + 3*ic + 2*3*id + 2*2*3*ipt;

          double nphoton = h_1photon[ic][part]->Integral(pTlow[id][ipt], pThigh[id][ipt]);
          double npion = h2_2photon[ic][part]->ProjectionY("h_minv", pTlow[id][ipt], pThigh[id][ipt])->Integral(111,160);
          npion /= 1.1;

          double xpT = pTbinC[id][ipt];
          int ipMiss = Get_ipt(xMiss[part], xpT);
          double aMiss = Miss[part][ipMiss];
          int ipMerge = Get_ipt(xMerge[part], xpT);
          double aMerge = Merge[part][ipMerge];
          int ipBadPass = Get_ipt(xBadPass[part], xpT);
          double aBadPass = BadPass[part][ipBadPass];
          if(ipt<23)
            aBadPass = 0.;
          double aMissPass = aMiss + aMerge * aBadPass;
          double aMissAll = aMiss + aMerge * 2.;
          double ndir = nphoton - aMissPass * npion - A * ( 1. + aMissAll ) * npion;
          double endir = sqrt( nphoton + pow(aMissPass + A * ( 1. + aMissAll ), 2.) * npion );

          double xx = (double)nmb / (double)nclock;
          double yy = ndir / nev;
          double eyy = endir / nev;
          if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
          {
            gr[ig]->SetPoint(igp[ig], xx, yy);
            gr[ig]->SetPointError(igp[ig], 0., eyy);
            gr_run[ig]->SetPoint(igp[ig], runnumber, yy);
            gr_run[ig]->SetPointError(igp[ig], 0., eyy);
            igp[ig]++;
          }
        } // ipt, id, ic, part

    delete f;
    irun++;
  }

  TFile *f_out = new TFile(Form("pileup/Pileup-photon-%d.root",process), "RECREATE");
  for(int ipt=0; ipt<npT; ipt++)
    for(int ic=0; ic<2; ic++)
      for(int part=0; part<3; part++)
      {
        int ig = part + 3*ic + 2*3*id + 2*2*3*ipt;
        gr[ig]->Set(igp[ig]);
        gr_run[ig]->Set(igp[ig]);
        gr[ig]->Write();
        gr_run[ig]->Write();
      }
  f_out->Close();
}
