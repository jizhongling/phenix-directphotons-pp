#include "GlobalVars.h"

void draw_IsoPhotonShuffle()
{
  const int ntype = 5*2*npT_pol;
  const char *type_names[5] = {"photon", "pion-pt1-peak", "pion-pt1-side", "pion-pt2-peak", "pion-pt2-side"};

  TH1 *h_ndf[ntype];
  TH1 *h_chi2[ntype];
  for(int imul=0; imul<5; imul++)
    for(int icr=0; icr<2; icr++)
      for(int ipt=0; ipt<npT_pol; ipt++)
      {
        int ih = imul + 5*icr + 5*2*ipt;
        h_ndf[ih] = new TH1F(Form("h_ndf_%d",ih), "NDF of #chi^{2}", 1000, -0.5, 999.5);
        h_chi2[ih] = new TH1F(Form("h_chi2_%d",ih), "#chi_{reduced}^{2}", 200, 0., 2.);
      }

  TFile *f = new TFile("data/isophoton-shuffle.root");
  TTree *t1 = (TTree*)f->Get("T");
  int ALLtype, ALLndf;
  double ALLchi2;
  t1->SetBranchAddress("Type", &ALLtype);
  t1->SetBranchAddress("NDF", &ALLndf);
  t1->SetBranchAddress("Chi2", &ALLchi2);
  for(int ien=0; ien<t1->GetEntries(); ien++)
  {
    t1->GetEntry(ien);
    int ih = ALLtype;
    h_ndf[ih]->Fill((double)ALLndf);
    h_chi2[ih]->Fill(ALLchi2);
  }

  TF1 *fn_chi2 = new TF1("fn_chi2", "ROOT::Math::chisquared_pdf(x*[0],[0])*[1]", 0., 2.);

  mc(0, 5,3);
  gStyle->SetOptStat(1110);
  c0->Print("plots/IsoPhoton-shuffle.pdf[");
  for(int imul=0; imul<5; imul++)
    for(int icr=0; icr<2; icr++)
    {
      for(int ipt=0; ipt<npT_pol; ipt++)
      {
        int ih = imul + 5*icr + 5*2*ipt;
        int ndf = TMath::Nint( h_ndf[ih]->GetMean() );
        fn_chi2->SetParameters((double)ndf, 1.);
        double norm = h_chi2[ih]->GetMaximum() / fn_chi2->Eval(1.);
        fn_chi2->SetParameters((double)ndf, norm);

        mcd(0, ipt+1);
        h_chi2[ih]->SetTitle( Form("%s-cr%d p_{T}: %3.1f-%3.1f",
              type_names[imul], icr, pTbin_pol[ipt], pTbin_pol[ipt+1]) );
        h_chi2[ih]->SetLineColor(kBlack);
        aset(h_chi2[ih], "#chi_{reduced}^{2}","", 0.,2.);
        fn_chi2->SetLineColor(kRed);
        fn_chi2->SetLineWidth(1);
        h_chi2[ih]->Draw();
        fn_chi2->DrawCopy("SAME");
      }
      c0->Print("plots/IsoPhoton-shuffle.pdf");
      c0->Clear("D");
    }
  c0->Print("plots/IsoPhoton-shuffle.pdf]");
}
