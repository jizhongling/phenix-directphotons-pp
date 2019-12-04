#include "GlobalVars.h"
#include "QueryTree.h"

void draw_MergePassRate()
{
  const int secl[3] = {1, 7, 1};
  const int sech[3] = {6, 8, 8};
  const char *pname[3] = {"PbSc", "PbGl", "Combined"};

  QueryTree *qt_badpass = new QueryTree("data/MergePassRate.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/MissingRatio-histo.root");
  THnSparse *hn_merge = (THnSparse*)f->Get("hn_merge");

  //mc(0, 2,4);
  mc();
  mcd();
  legi(0, 0.5,0.7,0.8,0.9);

  for(int ieta=7; ieta<8; ieta++)
    for(int part=0; part<3; part++)
    {
      //mcd(0, ieta+1);
      hn_merge->GetAxis(3)->SetRange(ieta+1-ieta/7*7,ieta+1-ieta/7*3);  // |eta|
      hn_merge->GetAxis(1)->SetRange(secl[part],sech[part]);  // sector

      hn_merge->GetAxis(2)->SetRange(1,2);  // passed = 0 or 1
      TH1 *h_total = hn_merge->Projection(0);

      hn_merge->GetAxis(2)->SetRange(2,2);  // passed = 1
      TH1 *h_passed = hn_merge->Projection(0);

      if(part == 2)
      {
        h_total = h_total->Rebin(npT_pol, "h_total_pol", pTbin_pol);
        h_passed = h_passed->Rebin(npT_pol, "h_passed_pol", pTbin_pol);
      }

      if(ieta==7)
        qt_badpass->Fill(h_passed, h_total, part);
      TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_passed, h_total);
      gr->SetTitle( Form("|#eta|: %.2f - %.2f",(ieta-ieta/7*7)*0.05,(ieta+1-ieta/7*3)*0.05) );
      aset(gr, "p_{T} [GeV]","Bad Pass", 16.,30., 0.,0.3);
      style(gr, part+20, part+1);
      if(part==0)
        gr->Draw("AP");
      else
        gr->Draw("P");
      leg0->AddEntry(gr, pname[part], "P");

      delete h_total;
      delete h_passed;
    }
  leg0->Draw();

  c0->Print("plots/MergePassRate.pdf");
  qt_badpass->Save();
}
