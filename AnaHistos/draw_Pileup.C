#include "Pileup.h"
#include "QueryTree.h"

void draw_Pileup()
{
  const char *dname[2] = {"ERT", "MB"};
  const char *mname[2] = {"pol1", "log"};
  const char *cname[4] = {"PbSc without ToF", "PbGl without ToF", "PbSc with ToF", "PbGl with ToF"};

  QueryTree *qt_fit = new QueryTree("data/Pileup-fit.root", "RECREATE");

  QueryTree *qt_pile = new QueryTree("data/Pileup.root");

  TF1 *fn_mean = new TF1("fn_mean", "pol0", 0., 0.2);
  TF1 *fn_fit[2];
  fn_fit[0] = new TF1("fn_pol1", "pol1", 0., 0.2);
  fn_fit[1] = new TF1("fn_log", "-[0]*log(1-[1]*x)/([1]*x)", 1e-3, 0.2);

  for(int i=0; i<5; i++)
    mc(i, 2,2);

  for(int ipt=0; ipt<npT; ipt++)
    for(int id=0; id<2; id++)
      for(int im=0; im<2; im++)
      {
        for(int is=0; is<2; is++)
        {
          int igr = id*4+im*2+is;
          double p0[2], ep0[2]; // p0[ic]
          double mean[2], emean[2]; // mean[ic]

          for(int ic=0; ic<2; ic++)
          {
            int cond = ic*2+is;
            int ig = ipt*8+id*4+ic*2+is;

            mcd(0, cond+1);
            TGraphErrors *gr = qt_pile->Graph(ig);
            gr->Draw("AP");  // must before GetXaxis()
            gr->SetTitle(cname[cond]);
            aset(gr, "Nmb/Nclock","Npi0/Nevent", 0.,0.2);
            style(gr, 20, 1);

            fn_fit[im]->SetParameters(gr->GetMaximum(), 1.);
            gr->Fit(fn_fit[im], "RQ");

            double scale = sqrt( fn_fit[im]->GetChisquare() / fn_fit[im]->GetNDF() );
            if(scale < 1.) scale = 1.;
            p0[ic] = fn_fit[im]->GetParameter(0);
            ep0[ic] = fn_fit[im]->GetParError(0) * scale;

            gr->Fit(fn_mean, "RQN");
            scale = sqrt( fn_mean->GetChisquare() / fn_mean->GetNDF() );
            if(scale < 1.) scale = 1.;
            mean[ic] = fn_mean->GetParameter(0);
            emean[ic] = fn_mean->GetParError(0) * scale;
          } // ic

          if( ipt > 0 )
          {
            double xpt = ( pTbin[id][ipt-1] + pTbin[id][ipt] ) / 2.;
            double yy = p0[1] / mean[1];
            double eyy = yy * sqrt( pow(emean[1]/mean[1],2) + pow(ep0[1]/p0[1],2) );
            if( TMath::Finite(yy + eyy) )
              qt_fit->Fill(ipt, igr, xpt, yy, eyy);

            yy = p0[1] / p0[0];
            eyy = yy * sqrt( pow(ep0[0]/p0[0],2) + pow(ep0[1]/p0[1],2) );
            if( TMath::Finite(yy + eyy) )
              qt_fit->Fill(ipt, 8+igr, xpt, yy, eyy);
          } // ipt > 0
        } // is

        f_out->cd();
        mcw( 0, Form("data%d-pt%d-%d-%s", id, pTlow[id][ipt], pThigh[id][ipt], mname[im]) );
        c0->Clear("D");
      } // ipt, id, im

  for(int id=0; id<2; id++)
    for(int im=0; im<2; im++)
      for(int is=0; is<2; is++)
      {
        int igr = id*4+im*2+is;
        TGraphErrors *gr_ratio = qt_fit->Graph(igr);
        TGraphErrors *gr_tof = qt_fit->Graph(8+igr);

        mcd(id+1, im*2+is+1);
        gr_ratio->SetTitle( Form("%s %s fit by %s", dname[id], cname[2+is], mname[im]) );
        aset(gr_ratio, "p_{T} [GeV/c]", "#frac{p0}{mean}");
        style(gr_ratio, 20, kRed);
        gr_ratio->Draw("AP");
        gr_ratio->Fit("pol0", "Q");

        mcd(id+3, im*2+is+1);
        gr_tof->SetTitle( Form("%s eff fit by %s in %s", cname[2+is], mname[im], dname[id]) );
        aset(gr_tof, "p_{T} [GeV/c]", "#ToF Eff");
        style(gr_tof, 20, kRed);
        gr_tof->Draw("AP");
        gr_tof->Fit("pol0", "Q");
      }

  c1->Print("plots/Pileup-ratio-ERT.pdf");
  c2->Print("plots/Pileup-ratio-MB.pdf");
  c3->Print("plots/ToFEff-ERT.pdf");
  c4->Print("plots/ToFEff-MB.pdf");
  qt_fit->Save();
}
