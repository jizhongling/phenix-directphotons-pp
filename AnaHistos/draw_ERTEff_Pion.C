#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "GetEfficiency.h"

void draw_ERTEff_Pion()
{
  const char *pname[2] = {"PbSc", "PbGl"};
  const int secl[2] = {0, 2};
  const int sech[2] = {1, 2};

  QueryTree *qt_ert = new QueryTree("data/ERTEff-pion.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-Sasha.root");

  // h[part][ert_trig]
  TH2 *h2_ert_pion[2][2];

  int evtype = 2;
  int checkmap = 1;
  int ival = 1;

  TH2 *h2_ert_pion_t = (TH2*)f->Get("h2_ert_pion_0");
  h2_ert_pion_t->Reset();
  for(int part=0; part<2; part++)
  {
    for(int ert_trig=0; ert_trig<2; ert_trig++)
    {
      h2_ert_pion[part][ert_trig] = (TH2*)h2_ert_pion_t->Clone(Form("h2_ert_pion_part%d_cond%d",part,ert_trig));
      for(int sector=secl[part]; sector<=sech[part]; sector++)
        for(int isolated=0; isolated<2; isolated++)
        {
          int ih = sector + 3*ert_trig + 3*2*evtype + 3*2*3*checkmap + 3*2*3*2*isolated + 3*2*3*2*2*ival;
          TH2 *h2_tmp = (TH2*)f->Get(Form("h2_ert_pion_%d",ih));
          h2_ert_pion[part][ert_trig]->Add(h2_tmp);
          delete h2_tmp;
        }
    }
    h2_ert_pion[part][0]->Add(h2_ert_pion[part][1]);
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
      h_minv = (TH1*)h2_ert_pion[part][0]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, nt, ent);
      delete h_minv;

      double np, enp;
      mcd(part*2+1, ipt+1);
      h_minv = (TH1*)h2_ert_pion[part][1]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, np, enp);
      delete h_minv;

      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double yy, eyyl, eyyh;
      if( !GetEfficiency(nt,np, yy,eyyl,eyyh) )
      {
        eyyl = yy * sqrt( pow(ent/nt,2) + pow(enp/np,2) );
        eyyh = 0.;
      }
      if( TMath::Finite(yy+eyyl+eyyh) )
        qt_ert->Fill(ipt, part, xx, yy, eyyl, eyyh);
    }

  mc(4);
  mcd(4);
  legi(0, 0.7,0.2,0.9,0.4);

  for(int part=0; part<2; part++)
  {
    TGraphAsymmErrors *gr = qt_ert->GraphAsymm(part);
    gr->SetTitle("ERT_4x4c trigger efficeincy for #pi^{0}");
    aset(gr, "p_{T} (GeV/c)","Eff", 1.,30., 0.,1.1);
    style(gr, part+20, part+1);
    if(part==0)
      gr->Draw("APE");
    else
      gr->Draw("PE");

    gr->Fit("pol0", "Q","", 10.,20.);
    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr->FindObject("stats");
    st->SetY1NDC(0.3-part*0.2);
    st->SetY2NDC(0.5-part*0.2);

    TGraph *gr_sasha =  new TGraph( Form("data/sasha-trig-part%d.txt",part) );
    gr_sasha->Draw("C");

    leg0->AddEntry(gr, pname[part], "LPE");
  }
  leg0->Draw();

  c4->Print("plots/ERTEff-pion.pdf");

  qt_ert->Write();
  for(int part=0; part<2; part++)
  {
    mcw( part*2, Form("part%d-total",part) );
    mcw( part*2+1, Form("part%d-passed",part) );
  }
  qt_ert->Close();
}
