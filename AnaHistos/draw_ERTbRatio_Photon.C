#include "GlobalVars.h"
#include "QueryTree.h"

void draw_ERTbRatio_Photon()
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  QueryTree *qt_ratio = new QueryTree("data/ERTbRatio-photon.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[evtype][part]
  TH1 *h_1photon[3][3];

  int checkmap = 1;
  int ival = 1;

  TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
  h_1photon_t = (TH1*)h_1photon_t->Clone();
  h_1photon_t->Reset();

  for(int evtype=0; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h_1photon[evtype][part] = (TH1*)h_1photon_t->Clone(Form("h_1photon_type%d_part%d",evtype,part));
      for(int isolated=0; isolated<2; isolated++)
      {
        int ih = part + 3*evtype + 3*3*checkmap + 3*3*2*isolated[ival] + 3*3*2*2*ival;
        TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
        h_1photon[evtype][part]->Add(h_tmp);
        delete h_tmp;
      }
    }

  for(int part=0; part<3; part++)
    for(int ipt=22; ipt<30; ipt++)
    {
      double nphoton_ertc = h_1photon[2][part]->GetBinContent(ipt+1);
      double nphoton_ertb = h_1photon[1][part]->GetBinContent(ipt+1);

      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      double ratio = nphoton_ertc / nphoton_ertb;
      double eratio = ratio * sqrt( 1./nphoton_ertc + 1./nphoton_ertb );
      if( TMath::Finite(ratio+eratio) )
        qt_ratio->Fill(ipt, part, xpt, ratio, eratio);
    }

  mc(0, 3,1);
  for(int part=0; part<3; part++)
  {
    mcd(0, part+1);
    TGraphErrors *gr = qt_ratio->Graph(part);
    aset(gr, "p_{T} [GeV]");
    style(gr, part+20, part+1);
    gr->Draw("AP");
    gr->Fit("pol0", "Q","", 10.,30.);
  }
  c0->Print("plots/ERTbRatio-photon.pdf");

  qt_ratio->Save();
}
