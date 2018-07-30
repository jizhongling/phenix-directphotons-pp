#include "Pileup.h"
#include "MultiGraph.h"

void draw_Pileup_Photon()
{
  const char *dname[2] = {"ERT", "MB"};
  const char *mname[2] = {"pol1", "log"};
  const char *cname[6] = {"PbScW without ToF", "PbScE without ToF", "PbGl without ToF", "PbScW with ToF", "PbScE with ToF", "PbGl with ToF"};

  TFile *f_out = new TFile("data/Pileup-isophoton-fit.root", "RECREATE");

  int id = 0;

  TMultiGraph *mg[npT*2*2*3];
  for(int ipt=0; ipt<npT; ipt++)
    for(int ic=0; ic<2; ic++)
      for(int part=0; part<3; part++)
      {
        int ig = part + 3*ic + 2*3*id + 2*2*3*ipt;
        mg[ig] = new TMultiGraph();
      }

  for(int i=0; i<10; i++)
  {
    TFile *f = new TFile(Form("pileup/Pileup-isophoton-%d.root",i));
    if( f->IsZombie() ) continue;

    for(int ipt=0; ipt<npT; ipt++)
      for(int ic=0; ic<2; ic++)
        for(int part=0; part<3; part++)
        {
          int ig = part + 3*ic + 2*3*id + 2*2*3*ipt;
          TGraphErrors *gr = (TGraphErrors*)f->Get(Form("gr_%d",ig));
          if( gr->GetN() > 0)
            mg[ig]->Add(gr);
        }
  }

  TGraphErrors *gr_ratio[2*2*3];
  TGraphErrors *gr_tof[2*2*3];
  int igp[2*2*3] = {};
  int igp2[2*2*3] = {};
  for(int im=0; im<2; im++)
    for(int part=0; part<3; part++)
    {
      int igr = part + 3*im + 2*3*id;
      gr_ratio[igr] = new TGraphErrors(npT);
      gr_tof[igr] = new TGraphErrors(npT);
    }

  TF1 *fn_mean = new TF1("fn_mean", "pol0", 0., 0.2);
  TF1 *fn_fit[2];
  fn_fit[0] = new TF1("fn_pol1", "pol1", 0., 0.2);
  fn_fit[1] = new TF1("fn_log", "-[0]*log(1-[1]*x)/([1]*x)", 1e-3, 0.2);

  for(int i=0; i<5; i++)
    mc(i, 3,2);

  for(int ipt=0; ipt<npT; ipt++)
    for(int im=0; im<2; im++)
    {
      for(int part=0; part<3; part++)
      {
        int igr = part + 3*im + 2*3*id;
        double p0[2], ep0[2]; // p0[ic]
        double mean[2], emean[2]; // mean[ic]

        for(int ic=0; ic<2; ic++)
        {
          int cond = part + 3*ic;
          int ig = part + 3*ic + 2*3*id + 2*2*3*ipt;

          mcd(0, cond+1);
          mg[ig]->Draw("AP");  // must before GetXaxis()
          mg[ig]->SetTitle(cname[cond]);
          mg[ig]->GetXaxis()->SetTitle("Nmb/Nclock");
          mg[ig]->GetYaxis()->SetTitle("Npi0/Nevent");
          mg[ig]->GetXaxis()->SetLimits(0., 0.2);  // Do not use SetRangeUser()
          //mg[ig]->GetYaxis()->SetRangeUser(0., 1e-3);  // Do not use SetLimits()

          fn_fit[im]->SetParameters(GetMaximum<TGraphErrors>(mg[ig]), 1.);
          mg[ig]->Fit(fn_fit[im], "RQ");

          double scale = sqrt( fn_fit[im]->GetChisquare() / fn_fit[im]->GetNDF() );
          if(scale < 1.) scale = 1.;
          p0[ic] = fn_fit[im]->GetParameter(0);
          ep0[ic] = fn_fit[im]->GetParError(0) * scale;

          mg[ig]->Fit(fn_mean, "RQN");
          scale = sqrt( fn_mean->GetChisquare() / fn_mean->GetNDF() );
          if(scale < 1.) scale = 1.;
          mean[ic] = fn_mean->GetParameter(0);
          emean[ic] = fn_mean->GetParError(0) * scale;
        } // ic

        if( ipt > 0 )
        {
          double xx = ( pTbin[id][ipt-1] + pTbin[id][ipt] ) / 2.;
          double yy = p0[1] / mean[1];
          double eyy = yy * sqrt( pow(emean[1]/mean[1],2.) + pow(ep0[1]/p0[1],2.) );
          if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
          {
            gr_ratio[igr]->SetPoint(igp[igr], xx, yy);
            gr_ratio[igr]->SetPointError(igp[igr], 0., eyy);
            igp[igr]++;
          }

          yy = p0[1] / p0[0];
          eyy = yy * sqrt( pow(ep0[0]/p0[0],2.) + pow(ep0[1]/p0[1],2.) );
          if( yy > 0. && eyy > 0. && eyy < 1. )
          {
            gr_tof[igr]->SetPoint(igp2[igr], xx, yy);
            gr_tof[igr]->SetPointError(igp2[igr], 0., eyy);
            igp2[igr]++;
          }
        } // ipt > 0
      } // part

      f_out->cd();
      mcw( 0, Form("data%d-pt%d-%d-%s", id, pTlow[id][ipt], pThigh[id][ipt], mname[im]) );
      c0->Clear("D");
    } // ipt, id, im

  for(int im=0; im<2; im++)
    for(int part=0; part<3; part++)
    {
      int igr = part + 3*im + 2*3*id;
      gr_ratio[igr]->Set(igp[igr]);
      gr_tof[igr]->Set(igp2[igr]);

      mcd(id+1, part+3*im+1);
      gr_ratio[igr]->SetTitle( Form("%s %s fit by %s", dname[id], cname[part+3], mname[im]) );
      aset(gr_ratio[igr], "pT [GeV]", "#frac{p0}{mean}");
      style(gr_ratio[igr], 20, kRed);
      gr_ratio[igr]->Draw("AP");
      gr_ratio[igr]->Fit("pol0", "Q");

      mcd(id+3, part+3*im+1);
      gr_tof[igr]->SetTitle( Form("%s eff fit by %s in %s", cname[part+3], mname[im], dname[id]) );
      aset(gr_tof[igr], "pT [GeV]", "#ToF Eff");
      style(gr_tof[igr], 20, kRed);
      gr_tof[igr]->Draw("AP");
      gr_tof[igr]->Fit("pol0", "Q");
    }

  f_out->Close();
  c1->Print("plots/Pileup-isophoton-ratio-ERT.pdf");
  c2->Print("plots/Pileup-isophoton-ratio-MB.pdf");
  c3->Print("plots/ToFEff-isophoton-ERT.pdf");
  c4->Print("plots/ToFEff-isophoton-MB.pdf");
}
