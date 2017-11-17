
int
compare_cones()
{
  gStyle->SetOptStat(0);

  TFile *f_directg = new TFile("keep/anaphpythia_directphoton.root", "OPEN");
  TFile *f_minbias = new TFile("keep/anaphpythia.root", "OPEN");

  THnSparse *hn_cone_directg = (THnSparse*)f_directg->Get("hn_EConePhoton");
  THnSparse *hn_cone_minbias = (THnSparse*)f_minbias->Get("hn_EConePhoton");

  /* select direct photon / other photons */
  hn_cone_directg->GetAxis(4)->SetRange(2,2);
  hn_cone_minbias->GetAxis(4)->SetRange(1,1);

  /* vector to store cut efficiencies */
  vector< float > v_cone_radius;
  vector< float > v_directg_below_thr;
  vector< float > v_minbias_below_thr;

  /* Name histograms and output files based on nbin_rad_X */
  TString name_directg_base = "h_directg";
  TString name_minbias_base = "h_minbias";

  /* select cone radius - loop over cone radii */
  for ( unsigned i = 2; i < 12; i++ )
    {
      /* Update names */
      TString name_directg = name_directg_base;
      TString name_minbias = name_minbias_base;

      name_directg.Append("_conebin_");
      name_directg+=i;

      name_minbias.Append("_conebin_");
      name_minbias+=i;

      /* Pick correct cone bin */
      hn_cone_directg->GetAxis(3)->SetRange(i,i);
      hn_cone_minbias->GetAxis(3)->SetRange(i,i);

      TH1F* temp_proj = (TH1F*)hn_cone_directg->Projection(3);
      double r_cone = temp_proj->GetMean();
      delete temp_proj;

      cout << "** Select cone radius: " << r_cone << " rad" << endl;
      v_cone_radius.push_back( r_cone );

      /* Project histograms */
      TH1F* h1_econe_directg = (TH1F*)hn_cone_directg->Projection(2);
      h1_econe_directg->SetName( name_directg );

      TH1F* h1_econe_minbias = (TH1F*)hn_cone_minbias->Projection(2);
      h1_econe_minbias->SetName( name_minbias );

      h1_econe_directg->SetLineColor(kBlue);
      //h1_econe_directg->SetFillColor(kBlue);

      h1_econe_minbias->SetLineColor(kGray+1);
      h1_econe_minbias->SetFillColor(kGray+1);

      h1_econe_directg->Scale( 1. / h1_econe_directg->Integral(0,-1) );
      h1_econe_minbias->Scale( 1. / h1_econe_minbias->Integral(0,-1) );

      /* set threshold */
      double thr = 0.1;
      int bin_thr = h1_econe_directg->FindBin( thr );

      /* Calculate integrals above and below threshold-
       * Include underflow bin (0) and overflow bin (-1) */
      double n_directg_total = h1_econe_directg->Integral(0,-1);
      double n_directg_below = h1_econe_directg->Integral(0,bin_thr);
      double n_directg_above = h1_econe_directg->Integral(bin_thr+1,-1);

      double n_minbias_total = h1_econe_minbias->Integral(0,-1);
      double n_minbias_below = h1_econe_minbias->Integral(0,bin_thr);
      double n_minbias_above = h1_econe_minbias->Integral(bin_thr+1,-1);

      cout << "Threshold bin: " << bin_thr << " at center " << h1_econe_minbias->GetBinCenter(bin_thr) << endl;

      cout << "Direct Photons: " << n_directg_total << " " << n_directg_below << " " << n_directg_above << endl;

      cout << "Other Photons:  " << n_minbias_total << " " << n_minbias_below << " " << n_minbias_above << endl;

      v_directg_below_thr.push_back( n_directg_below );
      v_minbias_below_thr.push_back( n_minbias_below );

      TCanvas *c1 = new TCanvas();
      c1->SetLogy(1);
      h1_econe_minbias->GetXaxis()->SetRangeUser(0,1);
      h1_econe_minbias->GetYaxis()->SetRangeUser(1e-4,1);
      h1_econe_minbias->DrawClone("");
      h1_econe_directg->DrawClone("same");

      TLine *l_threshold = new TLine(thr, 1e-4, thr, 1);
      l_threshold->SetLineStyle(2);
      l_threshold->SetLineWidth(2);
      l_threshold->SetLineColor(kRed);
      l_threshold->Draw();

      gPad->RedrawAxis();
      TString c1_plotname = "plots/";
      c1_plotname.Append( name_directg );
      c1_plotname.Append("_zoom.eps");
      cout << c1_plotname << endl;
      c1->Print( c1_plotname );

      //      delete c1;

      TCanvas *c2 = new TCanvas();
      c2->SetLogy(1);
      h1_econe_minbias->Rebin(100);
      h1_econe_directg->Rebin(100);
      h1_econe_minbias->GetXaxis()->SetRangeUser(1,30);
      h1_econe_minbias->GetYaxis()->SetRangeUser(1e-4,1);
      h1_econe_minbias->Draw("");
      h1_econe_directg->Draw("same");

      gPad->RedrawAxis();
      TString c2_plotname = "plots/";
      c2_plotname.Append( name_directg );
      c2_plotname.Append(".eps");
      c2->Print( c2_plotname );

      //      delete c2;
    }

  /* Plot cut efficiencies */
  TGraph* g_directg_below_thr = new TGraph( v_cone_radius.size(), &v_cone_radius[0], &v_directg_below_thr[0] );
  TGraph* g_minbias_below_thr = new TGraph( v_cone_radius.size(), &v_cone_radius[0], &v_minbias_below_thr[0] );

  g_directg_below_thr->SetLineColor(kBlue);
  g_directg_below_thr->SetMarkerColor(kBlue);
  g_directg_below_thr->SetMarkerStyle(20);

  g_minbias_below_thr->SetMarkerStyle(20);


  TH1F* hframe = new TH1F("hframe",";isolation cone angle [rad];efficiency;",10,0,1.1);
  hframe->GetYaxis()->SetRangeUser(1e-4,1);

  TCanvas *c3 = new TCanvas();
  hframe->Draw();
  g_directg_below_thr->Draw("LP");
  g_minbias_below_thr->Draw("LPsame");

  c3->SetLogy(1);
  gPad->RedrawAxis();

  c3->Print("plots/efficinecy_vs_conesize.eps");

  return 0;
}
