#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "BgGPRMinv.h"

void draw_BgRatio()
{
  gSystem->Load("libGausProc.so");

  const char *GPR_outfile = "data/BgGPR.root";
  const char *method[2] = {"gaus+pol3", "GPR"};
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
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(2);
  for(int part=0; part<2; part++)
  {
    mc(part, 4,4);

    for(int ipt=0; ipt<npT; ipt+=2)
    {
      int evtype = ipt<22 ? 2 : 1;

      TH1 *h_minv;
      mcd(part, ipt/2+1);
      double npeak = 1., enpeak = 1., nbg = 1., enbg = 1.;
      h_minv = (TH1*)h2_pion[evtype][part]->ProjectionY("h_py", ipt+1,ipt+2)->Clone("h_minv");

      BgGPRMinv(h_minv, npeak, enpeak, nbg, enbg, GPR_outfile, part+ipt);
      double xpt = (pTbin[ipt] + pTbin[ipt+2]) / 2.;
      double rbg = nbg/npeak;
      double erbg = rbg*sqrt(enpeak*enpeak/npeak/npeak + enbg*enbg/nbg/nbg);
      qt_rbg->Fill(ipt, part+2, xpt, 1-rbg, erbg);

      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+2]) );
      h_minv->Rebin(10);
      double minv_shift = ipt<22 ? 0. : 0.01;
      FitMinv(h_minv, npeak, enpeak, nbg, enbg, true, 0.11-minv_shift,0.16+minv_shift);
      xpt = (pTbin[ipt] + pTbin[ipt+2]) / 2.;
      rbg = nbg/npeak;
      erbg = rbg*sqrt(enpeak*enpeak/npeak/npeak + enbg*enbg/nbg/nbg);
      qt_rbg->Fill(ipt, part, xpt, 1-rbg, erbg);

      TFile *f_gpr = new TFile(GPR_outfile);
      TH1 *ho = (TH1*)f_gpr->Get(Form("ho_%d",part+ipt));
      ho->Scale(10.);
      ho->DrawCopy("SAME");
      f_gpr->Close();

      delete h_minv;
    }

    mcd(3);
    for(int im=0; im<2; im++)
    {
      TGraphErrors *gr = qt_rbg->Graph(part+2*im);
      aset(gr, "p_{T} [GeV]","r_{sig}", 6.1,30., 0.8,1.1);
      style(gr, part+4*im+20, part+1);
      char *opt = im==0&&part==0 ? "AP" : "P";
      gr->Draw(opt);
      cout << "im " << im << ", part " << part << ", low pT >>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
      gr->Fit("pol0","N","", 6.,14.);
      cout << "im " << im << ", part " << part << ", high pT >>>>>>>>>>>>>>>>>>>>>>>>" << endl;
      gr->Fit("pol0","N","", 14.,30.);
      leg0->AddEntry(gr, Form("%s %s",method[im],pname[part]), "P");
    }
  }
  leg0->Draw();
  c3->Print("plots/BgRatio.pdf");

  qt_rbg->Write();
  for(int part=0; part<2; part++)
    mcw( part, Form("Minv-part%d",part) );
  qt_rbg->Close();
}
