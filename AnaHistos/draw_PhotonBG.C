#include "GlobalVars.h"

void draw_PhotonBG()
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");
  THnSparse *hn_photonbg = (THnSparse*)f->Get("hn_photonbg");

  for(int part=0; part<3; part++)
  {
    for(int ic=0; ic<2; ic++)
      mc(ic*3+part, 5,6);
    mc(part+6, 5,6);
  }

  double ncounts[2] = {};
  for(int part=0; part<3; part++)
    for(int ic=1; ic>=0; ic--)
      for(int ipt=0; ipt<npT; ipt++)
      {
        hn_photonbg->GetAxis(4)->SetRange(ic+1,ic+1);
        hn_photonbg->GetAxis(0)->SetRange(secl[part],sech[part]);
        hn_photonbg->GetAxis(1)->SetRange(ipt+1,ipt+1);

        mcd(ic*3+part, ipt+1);
        TH2 *h2_corr = hn_photonbg->Projection(2,3);
        h2_corr->SetTitle( Form("p_{T}: %.1f-%.1f",pTbin[ipt],pTbin[ipt+1]) );
        h2_corr->GetXaxis()->SetRangeUser(0.93,1.);
        h2_corr->GetYaxis()->SetRangeUser(pTbin[ipt]*0.9,pTbin[ipt+1]*1.1);
        h2_corr->DrawCopy("COLZ");
        delete h2_corr;

        mcd(part+6, ipt+1);
        gPad->SetLogy();
        TH1 *h_angle = hn_photonbg->Projection(3);
        if( ic==0 && ipt>=20 )
          ncounts[part/2] += h_angle->Integral(0,-1);
        //TH1 *h_angle = new TH1F("h_angle", "#eta distributions;#eta;", 10, 0., 0.5);
        //for(int ib=1; ib<=h_rt->GetNbinsX(); ib++)
        //{
        //  double xx = h_rt->GetXaxis()->GetBinCenter(ib);
        //  double yy = h_rt->GetBinContent(ib);
        //  double theta = TMath::ASin(xx);
        //  double eta = -TMath::Log( TMath::Tan(theta/2) );
        //  h_angle->Fill(eta, yy);
        //}
        h_angle->SetTitle( Form("p_{T}: %.1f-%.1f",pTbin[ipt],pTbin[ipt+1]) );
        aset(h_angle);
        style(h_angle, ic+20, ic+1);
        if(ic==1)
          h_angle->DrawCopy();
        else
          h_angle->DrawCopy("SAME");
        delete h_angle;
      }
  cout << ncounts[0] << endl;
  cout << ncounts[1] << endl;

  TFile *f_out = new TFile("data/PhotonBG.root", "RECREATE");
  for(int part=0; part<3; part++)
  {
    for(int ic=0; ic<2; ic++)
      mcw( ic*3+part, Form("part%d-bbc%d",part,ic) );
    mcw( part+6, Form("angle-part%d",part) );
  }
  f_out->Close();
}
