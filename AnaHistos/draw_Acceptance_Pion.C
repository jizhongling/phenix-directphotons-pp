#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "GetEfficiency.h"

void draw_Acceptance_Pion()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  SetWeight();

  QueryTree *qt_acc = new QueryTree("data/Acceptance-pion.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo.root");

  TH1 *h_total = (TH1*)f->Get("h_pion");

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  TAxis *axis_pt0 = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);
  TAxis *axis_sec = hn_pion->GetAxis(3);
  TAxis *axis_peak = hn_pion->GetAxis(4);

  for(int part=0; part<3; part++)
    for(int ipt=0; ipt<npT; ipt++)
    {
      double xpt = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double ww = cross_pi0->Eval(xpt);

      TH1 *h_tt = (TH1*)h_total->Clone();
      h_tt->Scale(1./ww);
      double nt = h_tt->GetBinContent(ipt+1);
      double ent = sqrt(nt);
      delete h_tt;

      double np, enp;
      mcd(part, ipt+1);
      //axis_peak->SetRange(3,3);  // Require two peaks
      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);
      TH1 *h_minv = hn_pion->Projection(2);
      h_minv->Rebin(10);
      h_minv->Scale(1./ww);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, np, enp);
      // Do not use fit result, use all invariant mass instead
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
        qt_acc->Fill(ipt, part, xpt, yy, eyyl, eyyh);
    }


  mc(3);
  mcd(3);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);

  for(int part=0; part<3; part++)
  {
    TGraphAsymmErrors *gr = qt_acc->GraphAsymm(part);
    gr->SetTitle("#pi^{0} acceptance");
    aset(gr, "p_{T} [GeV]","acceptance", 0.,30., 0.,0.12);
    style(gr, part+20, part+1);
    if(part==0)
      gr->Draw("AP");
    else
      gr->Draw("P");
    leg0->AddEntry(gr, pname[part], "P");
    TGraph *gr_sasha =  new TGraph( Form("data/sasha-acc-part%d.txt",part) );
    gr_sasha->Draw("C");
  }

  leg0->Draw();
  c3->Print("plots/Acceptance-pion.pdf");

  qt_acc->Write();
  for(int part=0; part<3; part++)
    mcw( part, Form("part%d",part) );
  qt_acc->Close();
}
