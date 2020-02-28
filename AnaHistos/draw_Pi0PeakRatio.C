#include "GlobalVars.h"
#include "QueryTree.h"

void draw_Pi0PeakRatio()
{
  const int secl[3] = {1, 7};
  const int sech[3] = {6, 8};

  QueryTree *qt_ratio = new QueryTree("data/Pi0PeakRatio.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo.root");
  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");

  mc();
  mcd();
  for(int part=0; part<2; part++)
  {
    hn_pion->GetAxis(3)->SetRange(secl[part],sech[part]);
    for(int ipt=0; ipt<npT; ipt++)
    {
      hn_pion->GetAxis(0)->SetRange(ipt+1,ipt+1);
      TH1 *h_minv = hn_pion->Projection(2);

      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      int minv_shift = ipt<22 ? 0 : 10;
      double ent, nt = h_minv->IntegralAndError(0,-1,ent);
      double enp, np = h_minv->IntegralAndError(111-minv_shift,160+minv_shift,enp);
      double ratio = np/nt;
      double eratio = sqrt(enp*enp/np/np + ent*ent/nt/nt);
      if( TMath::Finite(ratio+eratio) )
        qt_ratio->Fill(ipt, part, xpt, ratio, eratio);

      delete h_minv;
    }

    TGraphErrors *gr = qt_ratio->Graph(part);
    aset(gr, "p_{T} [GeV]","Peak ratio", 5.,30.);
    style(gr, 20+part, 1+part);
    char *opt = part==0 ? "AP" : "P";
    gr->Draw(opt);
  }
  c0->Print("plots/Pi0PeakRatio.pdf");

  qt_ratio->Save();
}
