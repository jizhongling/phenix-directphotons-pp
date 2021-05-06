#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "GetEfficiency.h"

void draw_ProbEff_Pion()
{
  const int secl[2] = {0, 2};
  const int sech[2] = {1, 2};

  QueryTree *qt_prob = new QueryTree("data/ProbEff-pion.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-Sasha.root");

  // h[icut][part]
  TH2 *h2_pion[2][2];

  int evtype = 2;
  int bbc10cm = 1;
  int tof = 1;
  int checkmap = 1;
  int ival = 1;

  TH2 *h2_pion_t = (TH2*)f->Get("h2_pion_0");
  h2_pion_t = (TH2*)h2_pion_t->Clone();
  h2_pion_t->Reset();

  for(int part=0; part<2; part++)
  {
    h2_pion[0][part] = (TH2*)h2_pion_t->Clone(Form("h2_pion_icut0_part%d",part));
    h2_pion[1][part] = (TH2*)h2_pion_t->Clone(Form("h2_pion_icut1_part%d",part));
    for(int sector=secl[part]; sector<=sech[part]; sector++)
      for(int evenodd=0; evenodd<2; evenodd++)
        for(int pattern=0; pattern<3; pattern++)
          for(int isolated=0; isolated<2; isolated++)
            for(int prob=0; prob<2; prob++)
            {
              int ih = sector + 3*evenodd + 3*2*pattern + 3*2*3*evtype + 3*2*3*4*tof + 3*2*3*4*2*prob + 3*2*3*4*2*2*bbc10cm + 3*2*3*4*2*2*2*checkmap + 3*2*3*4*2*2*2*2*isolated + 3*2*3*4*2*2*2*2*2*ival;
              TH2 *h2_tmp = (TH2*)f->Get(Form("h2_pion_%d",ih));
              h2_pion[0][part]->Add(h2_tmp);
              if(prob == 1)
                h2_pion[1][part]->Add(h2_tmp);
              delete h2_tmp;
            } // isolated
  }

  for(int part=0; part<2; part++)
    for(int ic=0; ic<2; ic++)
      mc(part*2+ic, 6,5);

  for(int part=0; part<2; part++)
    for(int ipt=0; ipt<npT; ipt++)
    {
      TH1 *h_minv;

      double nt, ent;
      mcd(part*2, ipt+1);
      h_minv = (TH1*)h2_pion[0][part]->ProjectionY("_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, nt, ent);
      if( ipt > 22 )
        nt = h_minv->Integral(12,16);
      delete h_minv;

      double np, enp;
      mcd(part*2+1, ipt+1);
      h_minv = (TH1*)h2_pion[1][part]->ProjectionY("_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, np, enp);
      if( ipt > 22 )
        np = h_minv->Integral(12,16);
      delete h_minv;

      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double yy, eyyl, eyyh;
      if( !GetEfficiency(nt,np, yy,eyyl,eyyh) )
      {
        eyyl = yy * sqrt( pow(ent/nt,2) + pow(enp/np,2) );
        eyyh = 0.;
      }
      if( TMath::Finite(yy+eyyl+eyyh) )
        qt_prob->Fill(ipt, part, xx, yy, eyyl, eyyh);
    }

  mc(5, 2,1);
  TFile *f_sasha = new TFile("data/sasha-prob.root");
  for(int part=0; part<2; part++)
  {
    mcd(5, part+1);
    TGraphAsymmErrors *gr = qt_prob->GraphAsymm(part);
    if(part == 0)
      gr->SetTitle("Prob Eff for PbSc");
    else if(part == 1)
      gr->SetTitle("Prob Eff for PbGl");
    aset(gr, "p_{T} (GeV/c)","Prob Eff", 0.,30., 0.8,1.1);
    style(gr, part+20, part+1);
    gr->Draw("AP");
    TGraphErrors *gr_sasha = (TGraphErrors*)f_sasha->Get( Form("gr_%d",part) );
    gr_sasha->Draw("C");
  }
  c5->Print("plots/ProbEff-pion.pdf");

  qt_prob->Write();
  for(int part=0; part<2; part++)
  {
    mcw( part*2, Form("part%d-total",part) );
    mcw( part*2+1, Form("part%d-passed",part) );
  }
  qt_prob->Close();
}
