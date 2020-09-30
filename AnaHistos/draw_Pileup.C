#include "Pileup.h"
#include "QueryTree.h"

void draw_Pileup()
{
  const char *type[5] = {"photon", "isophoton", "2photon", "isoboth", "isopair"};
  const char *mname[2] = {"pol1", "log"};
  const char *pname[3] = {"PbSc west", "PbSc east", "PbGl"};

  QueryTree *qt_fit = new QueryTree("data/Pileup-fit.root", "RECREATE");

  QueryTree *qt_pile = new QueryTree("data/Pileup.root");

  TF1 *fn_mean = new TF1("fn_mean", "pol0", 0., 0.2);
  TF1 *fn_fit[2];
  fn_fit[0] = new TF1("fn_pol1", "pol1", 0., 0.2);
  fn_fit[1] = new TF1("fn_log", "-[0]*log(1-[1]*x)/([1]*x)", 1e-3, 0.2);

  mc(0, 3,2);
  c0->Print("plots/Pileup.pdf[");
  for(int iph=0; iph<5; iph++)
    for(int ipt=0; ipt<npT; ipt++)
    {
      for(int part=0; part<3; part++)
        for(int im=0; im<2; im++)
        {
          mcd(0, part+3*im+1);
          int ig = iph + 5*part + 5*3*ipt;
          TGraphErrors *gr = qt_pile->Graph(ig);
          gr->Draw("AP");  // must before GetXaxis()
          gr->SetTitle( Form("%s-%s-%s: %d-%d GeV", type[iph], mname[im], pname[part], pTlow[0][ipt], pThigh[0][ipt]) );
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
              int igr = iph + 5*part + 5*3*im;
              qt_fit->Fill(ipt, igr, xx, yy, eyy);
            }
          } // ipt > 0
        } // part, im

      c0->Print("plots/Pileup.pdf");
      if(iph<2 && ipt==0)
        c0->Print(Form("plots/Pileup-%s.pdf",type[iph]));
      c0->Clear("D");
    } // iph, ipt

  mc(1, 3,2);
  for(int iph=0; iph<5; iph++)
  {
    for(int im=0; im<2; im++)
      for(int part=0; part<3; part++)
      {
        int igr = iph + 5*part + 5*3*im;
        TGraphErrors *gr_ratio = qt_fit->Graph(igr);

        mcd(1, part+3*im+1);
        gr_ratio->SetTitle( Form("%s: %s fit by %s", type[iph], pname[part], mname[im]) );
        aset(gr_ratio, "p_{T} [GeV/c]", "#frac{p0}{mean}");
        style(gr_ratio, 20, kRed);
        gr_ratio->Draw("AP");
        gr_ratio->Fit("pol0", "Q");
      }
    c1->Print("plots/Pileup.pdf");
    if(iph<2)
      c1->Print(Form("plots/Pileup-%s-ratio.pdf",type[iph]));
    c1->Clear("D");
  }
  c1->Print("plots/Pileup.pdf]");

  qt_fit->Save();
}
