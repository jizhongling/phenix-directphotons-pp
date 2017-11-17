#include "GlobalVars.h"
#include "FitMinv.h"

void draw_ToFEff()
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};
  const char *name[2] = {"PbSc", "PbGl"};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonNode-histo.root");

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");

  TGraphErrors *gr_tof[2];
  Int_t igr[2] = {};
  for(Int_t part=0; part<2; part++)
    gr_tof[part] =  new TGraphErrors(npT);

  mc(0, 5,6);
  mc(1, 5,6);
  mc(2, 2,1);

  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      hn_pion->GetAxis(0)->SetRange(secl[part],sech[part]);
      TH1 *h_minv;

      mcd(0, ipt+1);
      Double_t nt, ent;
      hn_pion->GetAxis(3)->SetRange(1,1);
      hn_pion->GetAxis(1)->SetRange(ipt+1,ipt+1);
      h_minv = hn_pion->Projection(2); 
      h_minv->Rebin(10);
      h_minv->SetTitle(Form("pT: %4.2f-%4.2f",pTbin[ipt],pTbin[ipt+1]));
      FitMinv(h_minv, nt, ent);
      delete h_minv;

      mcd(1, ipt+1);
      Double_t np, enp;
      hn_pion->GetAxis(3)->SetRange(2,2);
      hn_pion->GetAxis(1)->SetRange(ipt+1,ipt+1);
      h_minv = hn_pion->Projection(2); 
      h_minv->Rebin(10);
      h_minv->SetTitle(Form("pT: %4.2f-%4.2f",pTbin[ipt],pTbin[ipt+1]));
      FitMinv(h_minv, np, enp);
      delete h_minv;

      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t yy = np / nt;
      Double_t eyy = yy * sqrt( pow(ent/nt,2.) + pow(enp/np,2.) );
      if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
      {
        gr_tof[part]->SetPoint(igr[part], xx, yy);
        gr_tof[part]->SetPointError(igr[part], 0., eyy);
        igr[part]++;
      }
    }

  for(Int_t part=0; part<2; part++)
  {
    mcd(2, part+1);
    gr_tof[part]->SetTitle( Form("ToF efficeincy for %s", name[part]) );
    aset(gr_tof[part], "p_{T} [GeV]", "Eff", 0.,30., 0.,1.1);
    style(gr_tof[part], 24, kRed);
    gr_tof[part]->Draw("APE");
    gr_tof[part]->Fit("pol0", "Q","", 9.,30.);

    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr_tof[part]->FindObject("stats");
    st->SetX1NDC(0.7);
    st->SetY1NDC(0.6);
    st->SetX2NDC(1.0);
    st->SetY2NDC(0.8);
  }

  c0->Print("ToFEff-PbSc.pdf");
  c1->Print("ToFEff-PbGl.pdf");
  c2->Print("ToFEff.pdf");
}
