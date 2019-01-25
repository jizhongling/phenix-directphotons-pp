#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"

void draw_PtRatio()
{
  const char *pname[2] = {"PbSc", "PbGl"};
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  QueryTree *qt_ptratio = new QueryTree("data/PtRatio.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[evtype][part]
  TH2 *h2_2photon[3][2];
  TH2 *h2_2photon2pt[3][2];

  int bbc10cm = 1;
  int ival = 1;

  TH2 *h2_2photon_t = (TH2*)f->Get("h2_2photon_0");
  h2_2photon_t = (TH2*)h2_2photon_t->Clone();
  h2_2photon_t->Reset();
  for(int evtype=1; evtype<3; evtype++)
    for(int part=0; part<2; part++)
    {
      h2_2photon[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_2photon_type%d_part%d",evtype,part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int pattern=0; pattern<3; pattern++)
          for(int evenodd=0; evenodd<2; evenodd++)
            for(int isoboth=0; isoboth<2; isoboth++)
              for(int isopair=0; isopair<2; isopair++)
              {
                int ih = sector + 8*evenodd + 2*8*pattern + 3*2*8*isoboth + 2*3*2*8*isopair + 2*2*3*2*8*evtype + 3*2*2*3*2*8*bbc10cm + 2*3*2*2*3*2*8*ival;
                TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
                h2_2photon[evtype][part]->Add(h2_tmp);
                delete h2_tmp;
              }
    }

  for(int evtype=1; evtype<3; evtype++)
    for(int part=0; part<2; part++)
    {
      h2_2photon2pt[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_2photon2pt_type%d_part%d",evtype,part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int pattern=0; pattern<3; pattern++)
          for(int evenodd=0; evenodd<2; evenodd++)
            for(int isoboth=0; isoboth<2; isoboth++)
              for(int isopair=0; isopair<2; isopair++)
              {
                int ih = sector + 8*evenodd + 2*8*pattern + 3*2*8*isoboth + 2*3*2*8*isopair + 2*2*3*2*8*evtype + 3*2*2*3*2*8*bbc10cm + 2*3*2*2*3*2*8*ival;
                TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon2pt_%d",ih));
                h2_2photon2pt[evtype][part]->Add(h2_tmp);
                delete h2_tmp;
              }
    }

  for(int part=0; part<2; part++)
  {
    mc(part, 6,5);
    mc(part+2, 6,5);
    for(int ipt=0; ipt<npT; ipt++)
    {
      int evtype = 2;
      if(ipt < 22)  // <14GeV use ERT_4x4c
        evtype = 2;
      else  // >14GeV use ERT_4x4b
        evtype = 1;

      mcd(part, ipt+1);
      double n2photon = 1., en2photon = 1.;
      TH1 *h_minv = (TH1*)h2_2photon[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, n2photon, en2photon, true, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, n2photon, en2photon, true, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, n2photon, en2photon, false, 0.10,0.17);
      n2photon /= bck[part/2][ipt] * meff[part/2][ipt];
      delete h_minv;

      mcd(part+2, ipt+1);
      double n2photon2pt = 1., en2photon2pt = 1.;
      h_minv = (TH1*)h2_2photon2pt[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, n2photon2pt, en2photon2pt, true, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, n2photon2pt, en2photon2pt, true, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, n2photon2pt, en2photon2pt, false, 0.10,0.17);
      n2photon2pt /= bck[part/2][ipt] * meff[part/2][ipt];
      delete h_minv;

      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      double ratio = n2photon2pt / n2photon;
      double eratio = ratio * sqrt( pow(en2photon2pt/n2photon2pt,2) + pow(en2photon/n2photon,2) );
      if( TMath::Finite(ratio+eratio) )
        qt_ptratio->Fill(ipt, part, xpt, ratio, eratio);
    }
  }

  mc(4);
  legi(0, 0.4,0.2,0.7,0.4);

  for(int part=0; part<2; part++)
  {
    TGraphErrors *gr = qt_ptratio->Graph(part);
    mcd(4);
    gr->SetTitle("N_{#gamma}(p_{T}^{2#gamma})/N_{#gamma}(p_{T}^{1#gamma})");
    aset(gr, "p_{T} [GeV]","Ratio", 0.,30., 0.,10.);
    style(gr, part+20, part+1);
    if(part==0)
      gr->Draw("AP");
    else
      gr->Draw("P");
    leg0->AddEntry(gr, pname[part], "P");
  }

  mcd(4);
  leg0->Draw();
  c4->Print("plots/PtRatio.pdf");

  qt_ptratio->Write();
  for(int part=0; part<2; part++)
  {
    mcw( part, Form("Minv-2photon-part%d",part) );
    mcw( part+2, Form("Minv-2photon2pt-part%d",part) );
  }
  qt_ptratio->Close();
}
