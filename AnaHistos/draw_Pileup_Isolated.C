#include "Pileup.h"
#include "QueryTree.h"

void draw_Pileup_Isolated()
{
  const char *type[3] = {"photon", "isoboth", "isopair"};
  const char *mname[2] = {"pol1", "log"};
  const char *pname[3] = {"PbSc west", "PbSc east", "PbGl"};

  QueryTree *qt_fit = new QueryTree("data/Pileup-isolated-fit.root", "RECREATE");

  QueryTree *qt_pile = new QueryTree("data/Pileup-isolated.root");

  TF1 *fn_mean = new TF1("fn_mean", "pol0", 0., 0.2);
  TF1 *fn_fit[2];
  fn_fit[0] = new TF1("fn_pol1", "pol1", 0., 0.2);
  fn_fit[1] = new TF1("fn_log", "-[0]*log(1-[1]*x)/([1]*x)", 1e-3, 0.2);

  mc(0, 3,1);
  c0->Print("plots/Pileup-isolated.pdf(", "pdf");
  for(int iph=0; iph<3; iph++)
    for(int ipt=0; ipt<npT; ipt++)
      for(int im=0; im<2; im++)
      {
        for(int part=0; part<3; part++)
        {
          mcd(0, part+1);
          int ig = iph + 3*part + 3*3*ipt;
          TGraphErrors *gr = qt_pile->Graph(ig);
          gr->Draw("AP");  // must before GetXaxis()
          gr->SetTitle( Form("%s-%s-%s: %d-%d", type[iph], mname[im], pname[part], pTlow[0][ipt], pThigh[0][ipt]) );
          aset(gr, "Nmb/Nclock",Form("%s/Nevent",type[iph]), 0.,0.2);
          style(gr, 20, 1);

          fn_fit[im]->SetParameters(gr->GetMaximum(), 1.);
          gr->Fit(fn_fit[im], "RQ");

          double scale = sqrt( fn_fit[im]->GetChisquare() / fn_fit[im]->GetNDF() );
          if(scale < 1.) scale = 1.;
          double p0 = fn_fit[im]->GetParameter(0);
          double ep0 = fn_fit[im]->GetParError(0) * scale;

          gr->Fit(fn_mean, "RQN");
          scale = sqrt( fn_mean->GetChisquare() / fn_mean->GetNDF() );
          if(scale < 1.) scale = 1.;
          double mean = fn_mean->GetParameter(0);
          double emean = fn_mean->GetParError(0) * scale;

          if( ipt > 0 )
          {
            double xx = ( pTbin[0][ipt-1] + pTbin[0][ipt] ) / 2.;
            double yy = p0 / mean;
            double eyy = yy * sqrt( pow(emean/mean,2) + pow(ep0/p0,2) );
            if( TMath::Finite(yy + eyy) )
            {
              int igr = iph + 3*part + 3*3*im;
              qt_fit->Fill(ipt, igr, xx, yy, eyy);
            }
          } // ipt > 0
        } // part

        c0->Print("plots/Pileup-isolated.pdf", "pdf");
        c0->Clear("D");
      } // iph, ipt, im

  mc(1, 3,3);
  for(int iph=0; iph<3; iph++)
  {
    for(int im=0; im<2; im++)
      for(int part=0; part<3; part++)
      {
        int igr = iph + 3*part + 3*3*im;
        TGraphErrors *gr_ratio = qt_fit->Graph(igr);

        mcd(1, part+3*im+1);
        gr_ratio->SetTitle( Form("%s: %s fit by %s", type[iph], pname[part], mname[im]) );
        aset(gr_ratio, "pT [GeV]", "#frac{p0}{mean}");
        style(gr_ratio, 20, kRed);
        gr_ratio->Draw("AP");
        gr_ratio->Fit("pol0", "Q");
      }
    c1->Print("plots/Pileup-isolated.pdf", "pdf");
    c1->Clear("D");
  }
  c1->Print("plots/Pileup-isolated.pdf)", "pdf");

  qt_fit->Save();
}
