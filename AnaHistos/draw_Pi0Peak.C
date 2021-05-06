#include "GlobalVars.h"
#include "GetMassWidth.h"

void draw_Pi0Peak()
{
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  TGraphErrors *gr_data_mass[2], *gr_data_width[2];
  TGraphErrors *gr_sim_mass[2], *gr_sim_width[2];
  int igp_data[2] = {};
  int igp_sim[2] = {};
  for(int part=0; part<2; part++)
  {
    gr_data_mass[part] =  new TGraphErrors(npT);
    gr_data_width[part] =  new TGraphErrors(npT);
    gr_sim_mass[part] =  new TGraphErrors(npT);
    gr_sim_width[part] =  new TGraphErrors(npT);
    for(int id=0; id<2; id++)
      mc(part*2+id, 6,5);
  }

  SetWeight();

  TFile *f_data = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-Sasha.root");
  TFile *f_sim = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-syserr.root");

  TH2 *h2_pion_data[2][3];

  int tof = 1;
  int prob = 1;
  int evtype = 2;
  int checkmap = 1;
  int ival = 1;

  TH2 *h2_2photon_t = (TH2*)f_data->Get("h2_2photon_0");
  h2_2photon_t = (TH2*)h2_2photon_t->Clone();
  h2_2photon_t->Reset();

  for(int pttype=0; pttype<2; pttype++)
  {
    for(int part=0; part<3; part++)
    {
      const char *ptname = pttype ? "2pt" : "";
      h2_pion_data[pttype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_pion_data_type%d_part%d",pttype,part));

      for(int isoboth=0; isoboth<2; isoboth++)
        for(int isopair=0; isopair<2; isopair++)
        {
          int ih = part + 3*evtype + 3*3*checkmap + 3*3*2*isoboth + 3*3*2*2*isopair + 3*3*2*2*2*ival;
          TH2 *h2_tmp = (TH2*)f_data->Get(Form("h2_2photon%s_%d",ptname,ih));
          h2_pion_data[pttype][part]->Add(h2_tmp);
        } // isoboth, isopair
    }
    h2_pion_data[pttype][0]->Add(h2_pion_data[pttype][1]);
    h2_pion_data[pttype][1] = h2_pion_data[pttype][2];
  }

  //TH2 *h2_pion_t = (TH2*)f_data->Get("h2_pion_0");
  //h2_pion_t = (TH2*)h2_pion_t->Clone();
  //h2_pion_t->Reset();

  //for(int part=0; part<3; part++)
  //{
  //  h2_pion_data[1][part] = (TH2*)h2_pion_t->Clone(Form("h2_pion_data_part%d",part));

  //  for(int isolated=0; isolated<2; isolated++)
  //  {
  //    int ih = part + 3*evtype + 3*3*tof + 3*3*2*prob + 3*3*2*2*checkmap + 3*3*2*2*2*isolated + 3*3*2*2*2*2*ival;
  //    TH2 *h2_tmp = (TH2*)f_data->Get(Form("h2_pion_%d",ih));
  //    h2_pion_data[1][part]->Add(h2_tmp);
  //  } // isolated
  //}
  //h2_pion_data[1][0]->Add(h2_pion_data[1][1]);
  //h2_pion_data[1][1] = h2_pion_data[1][2];

  THnSparse *hn_sim = (THnSparse*)f_sim->Get("hn_pion");
  hn_sim->GetAxis(7)->SetRange(1,1); // isys
  hn_sim->GetAxis(6)->SetRange(3,3); // econe_trk[ival]: EMCal, nomap, withmap

  mc(4, 2,2);

  for(int part=0; part<2; part++)
  {
    for(int ipt=0; ipt<npT; ipt++)
    {
      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double ww = cross_pi0->Eval(xx);
      double mass, emass, width, ewidth;
      TH1 *h_minv;

      mcd(part*2, ipt+1);
      h2_pion_data[1][part]->GetXaxis()->SetRange(ipt+1,ipt+1);
      h_minv = (TH1*)h2_pion_data[1][part]->ProjectionY()->Clone();
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

    aset(gr_data_mass[part], "p_{T} (GeV/c)","mass (GeV)", 0.,20., 0.13,0.145);
    aset(gr_data_width[part], "p_{T} (GeV/c)","width (GeV)", 0.,20., 0.007,0.015);
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
  for(int part=0; part<2; part++)
  {
    mcw( part*2, Form("part%d-data",part) );
    mcw( part*2+1, Form("part%d-sim",part) );
  }
  f_out->Write();
  f_out->Close();

  c4->Print("plots/Pi0Peak.pdf");
}
