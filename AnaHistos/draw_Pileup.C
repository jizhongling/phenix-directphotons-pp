#include "Pileup.h"
#include "MultiGraph.h"

void draw_Pileup()
{
  const char *dname[2] = {"ERT", "MB"};
  const char *mname[2] = {"pol1", "log"};
  const char *cname[4] = {"PbSc without ToF", "PbGl without ToF", "PbSc with ToF", "PbGl with ToF"};

  TMultiGraph *mg[npT*8];
  for(Int_t ipt=0; ipt<npT; ipt++)
    for(Int_t id=0; id<2; id++)
      for(Int_t ic=0; ic<2; ic++)
        for(Int_t is=0; is<2; is++)
        {
          Int_t ig = ipt*8+id*4+ic*2+is;
          mg[ig] = new TMultiGraph();
        }

  for(Int_t i=0; i<44; i++)
  {
    TFile *f = new TFile(Form("pileup/Pileup-%d.root",i));
    if( f->IsZombie() ) continue;

    for(Int_t ipt=0; ipt<npT; ipt++)
      for(Int_t id=0; id<2; id++)
        for(Int_t ic=0; ic<2; ic++)
          for(Int_t is=0; is<2; is++)
          {
            Int_t ig = ipt*8+id*4+ic*2+is;
            TGraphErrors *gr = (TGraphErrors*)f->Get(Form("gr_%d",ig));
            if( gr->GetN() > 0)
              mg[ig]->Add(gr);
          }
  }

  TGraphErrors *gr_ratio[8];
  TGraphErrors *gr_tof[8];
  Int_t igp[8] = {};
  Int_t igp2[8] = {};
  for(Int_t id=0; id<2; id++)
    for(Int_t im=0; im<2; im++)
      for(Int_t is=0; is<2; is++)
      {
        Int_t igr = id*4+im*2+is;
        gr_ratio[igr] = new TGraphErrors(npT);
        gr_tof[igr] = new TGraphErrors(npT);
      }

  TF1 *fn_mean = new TF1("fn_mean", "pol0", 0., 0.2);
  TF1 *fn_fit[2];
  fn_fit[0] = new TF1("fn_pol1", "pol1", 0., 0.2);
  fn_fit[1] = new TF1("fn_log", "-[0]*log(1-[1]*x)/([1]*x)", 1e-3, 0.2);

  for(Int_t i=0; i<6; i++)
    mc(i, 2,2);

  for(Int_t ipt=0; ipt<npT; ipt++)
    for(Int_t id=0; id<2; id++)
      for(Int_t im=0; im<2; im++)
      {
        for(Int_t is=0; is<2; is++)
        {
          Int_t igr = id*4+im*2+is;
          Double_t p0[2], ep0[2]; // p0[ic]
          Double_t mean[2], emean[2]; // mean[ic]

          for(Int_t ic=0; ic<2; ic++)
          {
            Int_t cond = ic*2+is;
            Int_t ig = ipt*8+id*4+ic*2+is;

            mcd(im, cond+1);
            mg[ig]->Draw("AP");  // must before GetXaxis()
            mg[ig]->SetTitle(cname[cond]);
            mg[ig]->GetXaxis()->SetTitle("Nmb/Nclock");
            mg[ig]->GetYaxis()->SetTitle("Npi0/Nevent");
            mg[ig]->GetXaxis()->SetLimits(0., 0.2);  // Do not use SetRangeUser()
            //mg[ig]->GetYaxis()->SetRangeUser(0., 1e-3);  // Do not use SetLimits()

            fn_fit[im]->SetParameters(GetMaximum<TGraphErrors>(mg[ig]), 1.);
            for(Int_t ifit=0; ifit<5; ifit++)
            {
              mg[ig]->Fit(fn_fit[im], "RQN");
              fn_fit[im]->SetParameters( fn_fit[im]->GetParameters() );
            }
            mg[ig]->Fit(fn_fit[im], "RQ");

            Double_t scale = sqrt( fn_fit[im]->GetChisquare() / fn_fit[im]->GetNDF() );
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
            Double_t xx = ( pTbin[id][ipt-1] + pTbin[id][ipt] ) / 2.;
            Double_t yy = p0[1] / mean[1];
            Double_t eyy = yy * sqrt( pow(emean[1]/mean[1],2.) + pow(ep0[1]/p0[1],2.) );
            if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
            {
              gr_ratio[igr]->SetPoint(igp[igr], xx, yy);
              gr_ratio[igr]->SetPointError(igp[igr], 0., eyy);
              igp[igr]++;
            }

            yy = p0[1] / p0[0];
            eyy = yy * sqrt( pow(ep0[0]/p0[0],2.) + pow(ep0[1]/p0[1],2.) );
            if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
            {
              gr_tof[igr]->SetPoint(igp2[igr], xx, yy);
              gr_tof[igr]->SetPointError(igp2[igr], 0., eyy);
              igp2[igr]++;
            }
          } // ipt > 0
          else if( ipt == 0 )
          {
            Double_t yy = p0[1] / p0[0];
            Double_t eyy = yy * sqrt( pow(ep0[0]/p0[0],2.) + pow(ep0[1]/p0[1],2.) );
            cout << Form("%s eff fit by %s in %s: ", cname[2+is], mname[im], dname[id]) << yy << " +- " << eyy << endl;
          } // ipt == 0
        } // is

        if( im == 0 )
        {
          c0->Print( Form("pileup/Pileup-data%d-pt%d-%d-pol1.pdf", id, pTlow[id][ipt], pThigh[id][ipt]) );
          c0->Clear("D");
        }
        else if( im == 1 )
        {
          c1->Print( Form("pileup/Pileup-data%d-pt%d-%d-log.pdf", id, pTlow[id][ipt], pThigh[id][ipt]) );
          c1->Clear("D");
        }
      } // ipt, id, im

  for(Int_t id=0; id<2; id++)
    for(Int_t im=0; im<2; im++)
      for(Int_t is=0; is<2; is++)
      {
        Int_t igr = id*4+im*2+is;
        gr_ratio[igr]->Set(igp[igr]);
        gr_tof[igr]->Set(igp2[igr]);

        mcd(id+2, im*2+is+1);
        gr_ratio[igr]->SetTitle( Form("%s %s fit by %s", dname[id], cname[2+is], mname[im]) );
        aset(gr_ratio[igr], "pT [GeV]", "#frac{p0}{mean}");
        style(gr_ratio[igr], 20, kRed);
        gr_ratio[igr]->Draw("AP");
        gr_ratio[igr]->Fit("pol0", "Q");

        mcd(id+4, im*2+is+1);
        gr_tof[igr]->SetTitle( Form("%s eff fit by %s in %s", cname[2+is], mname[im], dname[id]) );
        aset(gr_tof[igr], "pT [GeV]", "#ToF Eff");
        style(gr_tof[igr], 20, kRed);
        gr_tof[igr]->Draw("AP");
        gr_tof[igr]->Fit("pol0", "Q");
      }

  c2->Print("pileup/Pileup-ratio-ERT.pdf");
  c3->Print("pileup/Pileup-ratio-MB.pdf");
  c4->Print("pileup/ToF-ERT.pdf");
  c5->Print("pileup/ToF-MB.pdf");
}
