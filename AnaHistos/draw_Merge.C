#include "GlobalVars.h"
#include "QueryTree.h"
#include "GetEfficiency.h"

void draw_Merge()
{
  const char *pname[2] = {"PbSc", "PbGl"};
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  SetWeight();

  QueryTree *qt_merge = new QueryTree("data/Merge.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo.root");

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  TAxis *axis_pt0 = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);
  TAxis *axis_sec = hn_pion->GetAxis(3);
  TAxis *axis_peak = hn_pion->GetAxis(4);

  for(int part=0; part<2; part++)
    for(int ipt=0; ipt<npT; ipt++)
    {
      double xpt = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double ww = cross_pi0->Eval(xpt);

      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);
      TH1 *h_minv;

      double nt, ent;
      axis_peak->SetRange(2,3);  // Require one or two peaks
      h_minv = hn_pion->Projection(2);
      h_minv->Scale(1./ww);
      nt = h_minv->Integral(0,-1);
      ent = sqrt(nt);
      delete h_minv;

      double np, enp;
      axis_peak->SetRange(3,3);  // Require two peaks
      h_minv = hn_pion->Projection(2);
      h_minv->Scale(1./ww);
      np = h_minv->Integral(0,-1);
      enp = sqrt(np);
      delete h_minv;

      double yy, eyyl, eyyh;
      if( !GetEfficiency(nt,np, yy,eyyl,eyyh) )
      {
        eyyl = yy * sqrt( pow(ent/nt,2) + pow(enp/np,2) );
        eyyh = 0.;
      }
      if( TMath::Finite(yy+eyyl+eyyh) )
        qt_merge->Fill(ipt, part, xpt, yy, eyyl, eyyh);
    }


  mc();
  mcd();
  legi(0, 0.2,0.2,0.4,0.4);
  for(int part=0; part<2; part++)
  {
    TGraphAsymmErrors *gr = qt_merge->GraphAsymm(part);
    gr->SetTitle("Separating rate");
    aset(gr, "p_{T} [GeV/c]","rate", 0.,30., 0.,1.1);
    style(gr, part+20, part+1);
    if(part==0)
      gr->Draw("AP");
    else
      gr->Draw("P");
    leg0->AddEntry(gr, pname[part], "P");
    TGraph *gr_sasha =  new TGraph( Form("data/sasha-merge-part%d.txt",part) );
    gr_sasha->Draw("C");
  }
  leg0->Draw();
  c0->Print("plots/Merge.pdf");

  qt_merge->Save();
}
