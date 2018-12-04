#include "GlobalVars.h"
#include "Ftest.h"

void draw_RawALL()
{
  const char *pattern_list[5] = {"SOOSSOO", "OSSOOSS", "SSOO", "OOSS", "NONE"};
  const char *crossing_list[2] = {"even", "odd"};
  const int ngroup = 2;

  ofstream fout_txt("data/raw-asym.txt", ofstream::trunc);
  vector<double> *vp_ALL = new vector<double>[8];
  vector<double> *vp_eALL = new vector<double>[8];

  TFile *f_tree = new TFile("data/raw-asym.root");
  TTree *t1 = (TTree*)f_tree->Get("t1");
  for(int ipt=0; ipt<11; ipt++)
  {
    for(int pattern=0; pattern<4; pattern++)
      for(int icr=0; icr<2; icr++)
      {
        int id = icr + 2*pattern;
        vp_ALL[id].clear();
        vp_eALL[id].clear();

        TSQLResult *res = t1->Query( "ALL:eALL",
            Form("ipt==%d&&spin_pattern==%d&&crossing==%d",ipt,pattern,icr) );
        TSQLRow *row;
        while( row = res->Next() )
        {
          TString field0 = row->GetField(0);
          TString field1 = row->GetField(1);
          vp_ALL[id].push_back(field0.Atof());
          vp_eALL[id].push_back(field1.Atof());
          fout_txt << field0 << ",";
          delete row;
        }
        fout_txt << endl;
        delete res;
      }

    double F, p;
    Ftest(8, vp_ALL, vp_eALL, F, p);
    cout << ipt << ": " << F << ", " << p << endl;
  }
  fout_txt.close();

  TFile *f_out = new TFile("data/RawALL.root", "RECREATE");

  TMultiGraph *mg[4*2*npT/ngroup];
  for(int pattern=0; pattern<4; pattern++)
    for(int icr=0; icr<2; icr++)
      for(int ipt=0; ipt<npT/ngroup; ipt++)
      {
        int ig = ipt + npT/ngroup*icr + 2*npT/ngroup*pattern;
        mg[ig] = new TMultiGraph();
      }

  for(int i=0; i<20; i++)
  {
    TFile *f = new TFile(Form("histos/raw-asym-%d.root",i));
    if( f->IsZombie() ) continue;

    for(int pattern=0; pattern<4; pattern++)
      for(int icr=0; icr<2; icr++)
        for(int ipt=0; ipt<npT/ngroup; ipt++)
        {
          int ig = ipt + npT/ngroup*icr + 2*npT/ngroup*pattern;
          TGraphErrors *gr = (TGraphErrors*)f->Get(Form("gr_%d",ig));
          if( gr->GetN() > 0)
            mg[ig]->Add(gr);
        }
  }

  TGraphErrors *gr_all[4*2];
  int igp[4*2] = {};
  for(int pattern=0; pattern<4; pattern++)
    for(int icr=0; icr<2; icr++)
    {
      int igr = icr + 2*pattern;
      gr_all[igr] = new TGraphErrors(npT/ngroup);
      mc(igr, 6,5);
    }

  TF1 *fn_mean = new TF1("fn_mean", "pol0");

  for(int pattern=0; pattern<4; pattern++)
    for(int icr=0; icr<2; icr++)
    {
      int igr = icr + 2*pattern;

      for(int ipt=0; ipt<npT/ngroup; ipt++)
      {
        int ig = ipt + npT/ngroup*icr + 2*npT/ngroup*pattern;
        if( !mg[ig]->GetListOfGraphs() ) continue;

        mcd(igr, ipt+1);
        mg[ig]->Draw("AP");  // must before GetXaxis()
        mg[ig]->SetTitle(Form("p_{T}: %.1f-%.1f GeV",pTbin[ipt*ngroup],pTbin[(ipt+1)*ngroup]));
        mg[ig]->GetXaxis()->SetTitle("runnumber");
        mg[ig]->GetYaxis()->SetTitle("A_{LL}");
        //mg[ig]->GetXaxis()->SetLimits(386700., 398200.);  // Do not use SetRangeUser()
        //mg[ig]->GetYaxis()->SetRangeUser(-1., 1.);  // Do not use SetLimits()

        mg[ig]->Fit(fn_mean, "Q");
        double mean = fn_mean->GetParameter(0);
        double emean = fn_mean->GetParError(0);

        if( TMath::Finite(mean+emean) && emean > 0. )
        {
          double xpt = ( pTbin[ipt*ngroup] + pTbin[(ipt+1)*ngroup] ) / 2.;
          gr_all[igr]->SetPoint(igp[igr], xpt, mean);
          gr_all[igr]->SetPointError(igp[igr], 0., emean);
          igp[igr]++;
        }
      }

      f_out->cd();
      mcw( igr, Form("pattern%d-cross%d", pattern, icr) );
    }

  mc(8);
  mcd(8);
  legi(0, 0.2,0.8,0.7,0.9);
  leg0->SetNColumns(2);
  leg0->SetTextSize(0.02);
  for(int pattern=0; pattern<4; pattern++)
    for(int icr=0; icr<2; icr++)
    {
      int igr = icr + 2*pattern;
      if(igp[igr] <= 0) continue;
      gr_all[igr]->Set(igp[igr]);

      gr_all[igr]->SetTitle("Raw A_{LL}");
      aset(gr_all[igr], "p_{T} [GeV]","A_{LL}", 0.,30., -1.,1.);
      style(gr_all[igr], 1, 1+igr);
      gr_all[igr]->SetMarkerSize(0.);
      if(igr==0)
        gr_all[igr]->Draw("AP");
      else
        gr_all[igr]->Draw("P");
      leg0->AddEntry(gr_all[igr], Form("%s %s",pattern_list[pattern],crossing_list[icr]), "L");
    }
  leg0->Draw();
  c8->Print("plots/RawALL.pdf");

  f_out->Close();
}
