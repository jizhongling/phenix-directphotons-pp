#include "GlobalVars.h"
#include "QueryTree.h"

void draw_MissingRatio()
{
  const int secl[5] = {1, 5, 7, 1, 1};
  const int sech[5] = {4, 6, 8, 8, 8};
  const char *pname[5] = {"PbSc west", "PbSc east", "PbGl", "Combined", "Combined"};

  QueryTree *qt_miss = new QueryTree("data/MissingRatio.root", "RECREATE");
  QueryTree *qt_miss_eta = new QueryTree("data/MissingRatio-eta.root", "RECREATE");
  QueryTree *qt_merge1 = new QueryTree("data/Merge-1photon.root", "RECREATE");
  QueryTree *qt_merge2 = new QueryTree("data/Merge-2photon.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo.root");
  THnSparse *hn_missing = (THnSparse*)f->Get("hn_missing");
  THnSparse *hn_missing_eta = (THnSparse*)f->Get("hn_missing_eta");

  mc(0, 2,1);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  mc(1, 2,1);
  legi(1, 0.2,0.7,0.4,0.9);

  for(int part=0; part<5; part++)
  {
    int npT_rebin = part<4 ? npT : npT_pol;
    double *pTrebin = part<4 ? pTbin : pTbin_pol;

    mcd(0, 1);

    hn_missing->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_missing->GetAxis(3)->SetRange(3,3);
    hn_missing->GetAxis(4)->SetRange(3,3);
    TH1 *h_2photon = hn_missing->Projection(1);
    hn_missing->GetAxis(3)->SetRange(2,2);
    hn_missing->GetAxis(4)->SetRange(2,2);
    TH1 *h_1photon = hn_missing->Projection(1);
    h_2photon = h_2photon->Rebin(npT_rebin, "h_2photon", pTrebin);
    h_1photon = h_1photon->Rebin(npT_rebin, "h_1photon", pTrebin);

    qt_miss->Fill(h_1photon, h_2photon, part);
    if(part < 3)
    {
      TGraphErrors *gr_miss = qt_miss->Graph(part);
      gr_miss->SetNameTitle(Form("gr_%d",part), "Missing Ratio for #pi^{0}");
      aset(gr_miss, "p_{T}^{1#gamma} (GeV/c)","Missing ratio R", 5.,30., 0.,1.5);
      style(gr_miss, 20+part, 1+part);
      if(part == 0)
        gr_miss->Draw("AP");
      else
        gr_miss->Draw("P");
      leg0->AddEntry(gr_miss, pname[part], "P");
      leg0->Draw();
    }

    delete h_2photon;
    delete h_1photon;

    mcd(0, 2);

    hn_missing_eta->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_missing_eta->GetAxis(3)->SetRange(3,3);
    hn_missing_eta->GetAxis(4)->SetRange(3,3);
    h_2photon = hn_missing_eta->Projection(1);
    hn_missing_eta->GetAxis(3)->SetRange(2,2);
    hn_missing_eta->GetAxis(4)->SetRange(2,2);
    h_1photon = hn_missing_eta->Projection(1);
    h_2photon = h_2photon->Rebin(npT_rebin, "h_2photon", pTrebin);
    h_1photon = h_1photon->Rebin(npT_rebin, "h_1photon", pTrebin);

    qt_miss_eta->Fill(h_1photon, h_2photon, part);
    if(part < 3)
    {
      TGraphErrors *gr_miss_eta = qt_miss_eta->Graph(part);
      gr_miss_eta->SetNameTitle(Form("gr_%d",part), "Missing Ratio for #eta");
      aset(gr_miss_eta, "p_{T}^{1#gamma} (GeV/c)","Missing ratio R", 5.,30., 0.,1.5);
      style(gr_miss_eta, 20+part, 1+part);
      if(part == 0)
        gr_miss_eta->Draw("AP");
      else
        gr_miss_eta->Draw("P");
    }

    mcd(1, 1);
    gPad->SetLogy();

    hn_missing->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_missing->GetAxis(3)->SetRange(3,3);
    hn_missing->GetAxis(4)->SetRange(3,3);
    TH1 *h_separated = hn_missing->Projection(1);
    hn_missing->GetAxis(4)->SetRange(2,2);
    TH1 *h_merged = hn_missing->Projection(1);
    h_separated = h_separated->Rebin(npT_rebin, "h_separated", pTrebin);
    h_merged = h_merged->Rebin(npT_rebin, "h_merged", pTrebin);

    qt_merge1->Fill(h_merged, h_separated, part);
    if(part < 3)
    {
      TGraphErrors *gr_merge1 = qt_merge1->Graph(part);
      gr_merge1->SetNameTitle(Form("gr_%d",part+3), "Merging Ratio for one photon p_{T}");
      aset(gr_merge1, "p_{T}^{1#gamma} (GeV/c)","Merging ratio", 5.,30., 1e-4,10.);
      style(gr_merge1, 20+part, 1+part);
      if(part == 0)
        gr_merge1->Draw("AP");
      else
        gr_merge1->Draw("P");
      leg1->AddEntry(gr_merge1, pname[part], "P");
      leg1->Draw();
    }

    delete h_separated;
    delete h_merged;

    mcd(1, 2);
    gPad->SetLogy();

    hn_missing->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_missing->GetAxis(3)->SetRange(3,3);
    hn_missing->GetAxis(4)->SetRange(3,3);
    h_separated = hn_missing->Projection(0);
    hn_missing->GetAxis(4)->SetRange(2,2);
    h_merged = hn_missing->Projection(0);
    h_separated = h_separated->Rebin(npT_rebin, "h_separated", pTrebin);
    h_merged = h_merged->Rebin(npT_rebin, "h_merged", pTrebin);

    qt_merge2->Fill(h_merged, h_separated, part);
    if(part < 3)
    {
      TGraphErrors *gr_merge2 = qt_merge2->Graph(part);
      gr_merge2->SetNameTitle(Form("gr_%d",part), "Merging Ratio for two photon p_{T}");
      aset(gr_merge2, "p_{T}^{2#gamma} (GeV/c)","Merging ratio", 5.,30., 1e-4,20.);
      style(gr_merge2, 20+part, 1+part);
      if(part == 0)
        gr_merge2->Draw("AP");
      else
        gr_merge2->Draw("P");
    }

    delete h_separated;
    delete h_merged;
  }

  c0->Print("plots/MissingRatio.pdf");
  c1->Print("plots/Merge-photon.pdf");
  qt_miss->Save();
  qt_miss_eta->Save();
  qt_merge1->Save();
  qt_merge2->Save();
}
