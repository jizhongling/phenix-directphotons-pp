#include "GlobalVars.h"
#include "BBCCounts.h"

void anaBBCEff(const Int_t process = 0)
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  TGraphErrors *gr[npT*2];
  Int_t igp[npT*2] = {};
  for(Int_t ipt=0; ipt<npT; ipt++)
    for(Int_t part=0; part<2; part++)
    {
      Int_t ig = ipt*2+part;
      gr[ig] = new TGraphErrors(1000);
      gr[ig]->SetName(Form("gr_%d",ig));
    }

  const Int_t nThread = 1000;
  Int_t thread = -1;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  ReadClockCounts();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;
    if( thread%10 == 0 ) cout << "Nfiles = " << thread << endl;

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH3 *h3_trig = (TH3*)f->Get("h3_bbc_pion");

    ULong64_t nclock = GetClockLive(runnumber);
    ULong64_t nmb = GetBBCNovtxLive(runnumber);

    for(Int_t part=0; part<2; part++)
    {
      TH1 *h_trig = h3_trig->ProjectionX("h_trig", secl[part],sech[part], 2,2);
      TH1 *h_total = h3_trig->ProjectionX("h_total", secl[part],sech[part], 1,1);
      TGraphAsymmErrors *gr_trig = new TGraphAsymmErrors(h_trig, h_total);
      Double_t *gry = gr_trig->GetY();
      Double_t *egryl = gr_trig->GetEYlow();
      Double_t *egryh = gr_trig->GetEYhigh();

      for(Int_t ipt=0; ipt<npT; ipt++)
      {
        Int_t ig = ipt*2+part;
        Double_t xx = (Double_t)nmb / (Double_t)nclock;
        Double_t yy = gry[ipt];
        Double_t eyy = egryl[ipt] > egryh[ipt] ? egryl[ipt] : egryh[ipt];
        if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
        {
          gr[ig]->SetPoint(igp[ig], xx, yy);
          gr[ig]->SetPointError(igp[ig], 0., eyy);
          igp[ig]++;
        }
      }

      delete gr_trig;
      delete h_trig;
      delete h_total;
    }

    delete f;
  }

  TFile *f_out = new TFile("BBCEff.root", "RECREATE");
  for(Int_t ipt=0; ipt<npT; ipt++)
    for(Int_t part=0; part<2; part++)
    {
      Int_t ig = ipt*2+part;
      gr[ig]->Set(igp[ig]);
      gr[ig]->Write();
    }
  f_out->Close();
}

void draw_BBCEff_nocorr()
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};
  const char *name[2] = {"PbSc", "PbGl"};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  TH3 *h3_trig;
  h3_trig = (TH3*)f->Get("h3_bbc_pion");

  mc(0, 2,1);

  for(Int_t part=0; part<2; part++)
  {
    TH1 *h_trig = h3_trig->ProjectionX("h_trig", secl[part],sech[part], 2,2);
    TH1 *h_total = h3_trig->ProjectionX("h_total", secl[part],sech[part], 1,1);
    TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_trig, h_total);

    mcd(0, part+1);
    gr->SetTitle( Form("BBC trigger efficeincy for %s", name[part]) );
    aset(gr, "p_{T} [GeV]", "Eff", 0.,30., 0.,1.);
    style(gr, 24, kRed);
    gr->Draw("APE");
    gr->Fit("pol0", "Q","", 0.,20.);

    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr->FindObject("stats");
    st->SetX1NDC(0.7);
    st->SetY1NDC(0.6);
    st->SetX2NDC(1.0);
    st->SetY2NDC(0.8);
  }

  c0->Print("BBCEff.pdf");
}

void draw_BBCEff()
{
  TFile *f = new TFile("BBCEff.root");

  TF1 *fn_pol1 = new TF1("fn_pol1", "pol1");
  Int_t igp[2] = {};
  TGraphErrors *gr_trig[2];
  for(Int_t part=0; part<2; part++)
    gr_trig[part] = new TGraphErrors(30);

  mc(0, 5,6);
  mc(1, 5,6);
  mc(2, 2,1);

  for(Int_t ipt=0; ipt<npT; ipt++)
    for(Int_t part=0; part<2; part++)
    {
      Int_t ig = ipt*2+part;
      mcd(part, ipt+1);
      TGraphErrors *gr = (TGraphErrors*)f->Get(Form("gr_%d",ig));
      gr->SetTitle( Form("pT: %4.2f-%4.2f", pTbin[ipt], pTbin[ipt+1]) );
      aset(gr, "Nmb/Nclock","Eff", 0.,0.6, 0.85,1.05);
      gr->Draw("AP");
      gr->Fit(fn_pol1, "Q");

      Double_t scale = sqrt( fn_pol1->GetChisquare() / fn_pol1->GetNDF() );
      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t yy = fn_pol1->GetParameter(0);
      Double_t eyy = fn_pol1->GetParError(0) * scale;
      if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
      {
        gr_trig[part]->SetPoint(igp[part], xx, yy);
        gr_trig[part]->SetPointError(igp[part], 0., eyy);
        igp[part]++;
      }
    }

  for(Int_t part=0; part<2; part++)
  {
    mcd(2, part+1);
    gr_trig[part]->Set(igp[part]);
    aset(gr_trig[part], "pT [GeV]","Eff", 0.,30., 0.,1.);
    style(gr_trig[part], 24, kRed);
    gr_trig[part]->Draw("AP");
    gr_trig[part]->Fit("pol0", "Q","", 0.,20.);

    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr_trig[part]->FindObject("stats");
    st->SetX1NDC(0.7);
    st->SetY1NDC(0.6);
    st->SetX2NDC(1.0);
    st->SetY2NDC(0.8);
  }

  TFile *f_out = new TFile("BBCEff-pileup.root", "RECREATE");
  mcw(0, "PbSc");
  mcw(1, "PbGl");
  f_out->Close();

  c2->Print("BBCEff.pdf");
}
