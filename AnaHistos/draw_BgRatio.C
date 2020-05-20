#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "BgGPRMinv.h"

void draw_BgRatio()
{
  gSystem->Load("libGausProc.so");

  const char *GPR_outfile = "data/BgGPR.root";
  const char *pname[2] = {"PbSc", "PbGl"};

  gSystem->Exec(Form("rm -f %s",GPR_outfile));
  QueryTree *qt_rbg = new QueryTree("data/BgRatio.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[evtype][part]
  TH2 *h2_pion[3][3];

  int tof = 1;
  int prob = 1;
  int checkmap = 1;
  int ival = 1;

  TH2 *h2_pion_t = (TH2*)f->Get("h2_pion_0");
  h2_pion_t = (TH2*)h2_pion_t->Clone();
  h2_pion_t->Reset();

  for(int evtype=0; evtype<3; evtype++)
  {
    for(int part=0; part<3; part++)
    {
      h2_pion[evtype][part] = (TH2*)h2_pion_t->Clone(Form("h2_pion_type%d_part%d",evtype,part));
      for(int isolated=0; isolated<2; isolated++)
      {
        int ih = part + 3*evtype + 3*3*tof + 3*3*2*prob + 3*3*2*2*checkmap + 3*3*2*2*2*isolated + 3*3*2*2*2*2*ival;
        TH2 *h2_tmp = (TH2*)f->Get(Form("h2_pion_%d",ih));
        h2_pion[evtype][part]->Add(h2_tmp);
      } // isolated
    }
    h2_pion[evtype][0]->Add(h2_pion[evtype][1]);
    h2_pion[evtype][1] = h2_pion[evtype][2];
  }

  mc(3);
  legi(0, 0.2,0.8,0.7,0.9);
  leg0->SetNColumns(2);
  for(int part=0; part<2; part++)
  {
    mc(part, 6,5);

    for(int ipt=0; ipt<npT; ipt++)
    {
      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      double rbg[2], erbg[2];
      int evtype = ipt<22 ? 2 : 1;

      TH1 *h_minv;
      mcd(part, ipt+1);
      double npeak = 1., enpeak = 1., nbg = 1., enbg = 1.;
      h_minv = (TH1*)h2_pion[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");

      BgGPRMinv(h_minv, npeak, enpeak, nbg, enbg, GPR_outfile, part+2*ipt);
      rbg[1] = nbg/npeak;
      if(TMath::Finite(enbg))
        erbg[1] = rbg[1]*sqrt(enpeak*enpeak/npeak/npeak + enbg*enbg/nbg/nbg);
      else
        erbg[1] = fabs(rbg[1]);

      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      h_minv->Rebin(10);
      FitMinv(h_minv, npeak, enpeak, nbg, enbg, true, 0.11,0.16);
      rbg[0] = nbg/npeak;
      erbg[0] = rbg[0]*sqrt(enpeak*enpeak/npeak/npeak + enbg*enbg/nbg/nbg);

      TFile *f_gpr = new TFile(GPR_outfile);
      TH1 *ho = (TH1*)f_gpr->Get(Form("ho_%d",part+2*ipt));
      ho->Scale(10.);
      ho->DrawCopy("SAME");
      f_gpr->Close();

      double rsys = rbg[0] - rbg[1];
      double ersys = sqrt(erbg[0]*erbg[0] + erbg[1]*erbg[1]);
      qt_rbg->Fill(ipt, part, xpt, rsys, ersys);

      delete h_minv;
    }

    mcd(3);
    TGraphErrors *gr = qt_rbg->Graph(part);
    aset(gr, "p_{T} [GeV]","#Delta r_{sig}", 6.1,16., -0.02,0.1);
    style(gr, part+20, part+1);
    char *opt = part==0 ? "AP" : "P";
    gr->Draw(opt);
    leg0->AddEntry(gr, pname[part], "P");
  }
  leg0->Draw();
  c3->Print("plots/BgRatio.pdf");

  qt_rbg->Write();
  for(int part=0; part<2; part++)
    mcw( part, Form("Minv-part%d",part) );
  qt_rbg->Close();
}
