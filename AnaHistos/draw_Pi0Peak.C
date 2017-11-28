#include "GlobalVars.h"
#include "GetMassWidth.h"

void draw_Pi0Peak()
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  TGraphErrors *gr_data_mass[2], *gr_data_width[2];
  TGraphErrors *gr_sim_mass[2], *gr_sim_width[2];
  Int_t igp_data[2] = {};
  Int_t igp_sim[2] = {};
  for(Int_t part=0; part<2; part++)
  {
    gr_data_mass[part] =  new TGraphErrors(npT);
    gr_data_width[part] =  new TGraphErrors(npT);
    gr_sim_mass[part] =  new TGraphErrors(npT);
    gr_sim_width[part] =  new TGraphErrors(npT);
    for(Int_t id=0; id<2; id++)
      mc(part*2+id, 6,5);
  }

  SetWeight();

  TFile *f_data = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");
  TFile *f_sim = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-warn-histo.root");
  //TFile *f_sim = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/MissingRatio-histo.root");

  THnSparse *hn_data = (THnSparse*)f_data->Get("hn_pion");
  THnSparse *hn_sim = (THnSparse*)f_sim->Get("hn_pion");

  mc(4, 2,2);

  for(Int_t part=0; part<2; part++)
  {
    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t ww = cross->Eval(xx);
      Double_t mass, emass, width, ewidth;
      TH1 *h_minv;

      mcd(part*2, ipt+1);
      hn_data->GetAxis(4)->SetRange(3,3);
      hn_data->GetAxis(3)->SetRange(4,4);
      hn_data->GetAxis(0)->SetRange(secl[part],sech[part]);
      hn_data->GetAxis(1)->SetRange(ipt+1,ipt+1);
      h_minv = hn_data->Projection(2);
      h_minv->Rebin(10);
      if( GetMassWidth(h_minv, mass, emass, width, ewidth) )
      {
        gr_data_mass[part]->SetPoint(igp_data[part], xx, mass);
        gr_data_mass[part]->SetPointError(igp_data[part], 0., emass);
        gr_data_width[part]->SetPoint(igp_data[part], xx, width);
        gr_data_width[part]->SetPointError(igp_data[part], 0., ewidth);
        igp_data[part]++;
      }
      delete h_minv;

      mcd(part*2+1, ipt+1);
      hn_sim->GetAxis(3)->SetRange(secl[part],sech[part]);
      hn_sim->GetAxis(1)->SetRange(ipt+1,ipt+1);
      h_minv = hn_sim->Projection(2);
      h_minv->Scale(1./ww);
      h_minv->Rebin(10);
      if( GetMassWidth(h_minv, mass, emass, width, ewidth) )
      {
        gr_sim_mass[part]->SetPoint(igp_sim[part], xx, mass);
        gr_sim_mass[part]->SetPointError(igp_sim[part], 0., emass);
        gr_sim_width[part]->SetPoint(igp_sim[part], xx, width);
        gr_sim_width[part]->SetPointError(igp_sim[part], 0., ewidth);
        igp_sim[part]++;
      }
      delete h_minv;
    }

    gr_data_mass[part]->Set(igp_data[part]);
    gr_data_width[part]->Set(igp_data[part]);
    gr_sim_mass[part]->Set(igp_sim[part]);
    gr_sim_width[part]->Set(igp_sim[part]);

    aset(gr_data_mass[part], "p_{T} [GeV]","mass [GeV]", 0.,20., 0.13,0.145);
    aset(gr_data_width[part], "p_{T} [GeV]","width [GeV]", 0.,20., 0.007,0.015);
    style(gr_data_mass[part], 24, kRed);
    style(gr_data_width[part], 24, kRed);
    style(gr_sim_mass[part], 24, kBlue);
    style(gr_sim_width[part], 24, kBlue);

    mcd(4, part+1);
    gr_data_mass[part]->Draw("AP");
    gr_sim_mass[part]->Draw("P");

    mcd(4, part+3);
    gr_data_width[part]->Draw("AP");
    gr_sim_width[part]->Draw("P");
  }

  gr_data_mass[0]->SetTitle("PbSc m_{#gamma#gamma}");
  gr_data_mass[1]->SetTitle("PbGl m_{#gamma#gamma}");
  gr_data_width[0]->SetTitle("PbSc #sigma_{#gamma#gamma}");
  gr_data_width[1]->SetTitle("PbGl #sigma_{#gamma#gamma}");

  TFile *f_out = new TFile("data/Pi0Peak.root", "RECREATE");
  for(Int_t part=0; part<2; part++)
  {
    mcw( part*2, Form("part%d-data",part) );
    mcw( part*2+1, Form("part%d-sim",part) );
  }
  f_out->Write();
  f_out->Close();

  c4->Print("plots/Pi0Peak.pdf");
}
