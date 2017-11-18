int extract_raw_yield_pi0( unsigned run_index, int run_number, string histfile, string histname, bool visualize )
{

  gStyle->SetOptStat(1);

  /* Get number of events in run */
  int nevents = 0;

  /* Open THnSparse histogram file */
  TFile *f_in = new TFile( histfile.c_str(), "OPEN" );
  THnSparse* hn_pi0_data = (THnSparse*) f_in->Get( histname.c_str() );
  cout << "Entries (all sectors): " << hn_pi0_data->GetEntries() << endl;

  /* Print histogram information
   */
  hn_pi0_data->Print("a");

  /* Apply slices and selection on THnSparse
   */
  // @TODO
  double ptmin = 2.0;
  double ptmax = 29.0;
  Int_t ptmin_bin = hn_pi0_data->GetAxis(1)->FindBin(ptmin);
  Int_t ptmax_bin = hn_pi0_data->GetAxis(1)->FindBin(ptmax);

  cout << "ptmin_bin: " << ptmin_bin << " center at " << hn_pi0_data->GetAxis(1)->GetBinCenter(ptmin_bin) << endl;
  cout << "ptmax_bin: " << ptmax_bin << " center at " << hn_pi0_data->GetAxis(1)->GetBinCenter(ptmax_bin) << endl;

  /* Apply cut */
  hn_pi0_data->GetAxis(1)->SetRange(ptmin_bin,ptmax_bin); //pT
  hn_pi0_data->GetAxis(5)->SetRange(4,4); //ERT trigger

  /* Project sector bins into single 1-D histogram:
   */
  TH1F* h_inv_mass_select = (TH1F*)hn_pi0_data->Projection(2);
  cout << "X axis of 1-D histogram: " << h_inv_mass_select->GetXaxis()->GetTitle() << endl;
  cout << "Total entries: " << h_inv_mass_select->GetEntries() << endl;;
  cout << "Underflow entries: " << h_inv_mass_select->GetBinContent(0) << endl;;
  cout << "Overflow entries: " << h_inv_mass_select->GetBinContent(-1) << endl;;

  /* DEBUG step 1: Count entries, do side band subtraction */
  Int_t bin047 =  47+1; //(h_inv_mass_select)->FindBin(0.047);
  Int_t bin097 =  97; //(h_inv_mass_select)->FindBin(0.097);
  Int_t bin112 = 112+1; //(h_inv_mass_select)->FindBin(0.112);
  Int_t bin162 = 162; //(h_inv_mass_select)->FindBin(0.162);
  Int_t bin187 = 187+1; //(h_inv_mass_select)->FindBin(0.187);
  Int_t bin227 = 227; //(h_inv_mass_select)->FindBin(0.227);

  cout << "centermin_bin: " << bin112 << endl;
  cout << "centermax_bin: " << bin162 << endl;

  float raw_yield_count = (h_inv_mass_select)->Integral(bin112,bin162);

  //  float raw_yield_count = (h_inv_mass_select)->Integral(bin112,bin162) - ( (h_inv_mass_select)->Integral(bin047,bin097) + (h_inv_mass_select)->Integral(bin187,bin227) ) / 2.;
  /* END side band subtraction */

  /* print results */
  cout << "Raw pi0 yield count (run, yield): " << run_number << " " << raw_yield_count << endl;

  TCanvas *c1 = new TCanvas();
  h_inv_mass_select->Draw();

  return 0;



  ///////////////////////////////////////////////////////////
  /* Below: separate by sector */

  /* parameters for Gaussian fit to pi0 mass peak */
  //float fitrange_min = 0.105;
  //float fitrange_max = 0.165;
  //
  ///* collect sector information in vectors */
  //vector< int > v_sector;
  //vector< float > v_mean;
  //vector< float > v_dmean;
  //vector< float > v_sigma;
  //vector< float > v_dsigma;
  //vector< float > v_chisquare;
  //
  //vector< float > v_raw_pi0_count;

  //
  //  TH1F* h_inv_mass_sector[8];
  //  for ( unsigned int s = 0; s < 8; s++ )
  //    {
  //      TString hname_s("h_inv_mass_sector");
  //      hname_s+=s;
  //
  //      h_inv_mass_sector[s] = (TH1F*)hn_pi0_data->ProjectionZ( hname_s, s+1, s+1, ptmin_bin, ptmax_bin, "e" ); // sector/iterator count from 0, ROOT bins count from 1
  //    }
  //
  //  /* combine sectors */
  //  TH1F* h_inv_mass_combined[3];
  //  h_inv_mass_combined[0] = (TH1F*)h_inv_mass_sector[0]->Clone("h_inv_mass_PbScW");
  //  h_inv_mass_combined[0]->Add( h_inv_mass_sector[1] );
  //  h_inv_mass_combined[0]->Add( h_inv_mass_sector[2] );
  //  h_inv_mass_combined[0]->Add( h_inv_mass_sector[3] );
  //
  //  h_inv_mass_combined[1] = (TH1F*)h_inv_mass_sector[4]->Clone("h_inv_mass_PbScE");
  //  h_inv_mass_combined[1]->Add( h_inv_mass_sector[5] );
  //
  //  h_inv_mass_combined[2] = (TH1F*)h_inv_mass_sector[6]->Clone("h_inv_mass_PbGlE");
  //  h_inv_mass_combined[2]->Add( h_inv_mass_sector[7] );
  //
  //  for ( unsigned int s = 0; s < 3; s++ )
  //    {
  //      /* DEBUG step 1: Count entries, do side band subtraction */
  //      Int_t bin047 = (h_inv_mass_combined[s])->FindBin(0.047);
  //      Int_t bin097 = (h_inv_mass_combined[s])->FindBin(0.097);
  //      Int_t bin112 = (h_inv_mass_combined[s])->FindBin(0.112);
  //      Int_t bin162 = (h_inv_mass_combined[s])->FindBin(0.162);
  //      Int_t bin187 = (h_inv_mass_combined[s])->FindBin(0.187);
  //      Int_t bin227 = (h_inv_mass_combined[s])->FindBin(0.227);
  //
  //      float raw_yield_count = (h_inv_mass_combined[s])->Integral(bin112,bin162) - ( (h_inv_mass_combined[s])->Integral(bin047,bin097) + (h_inv_mass_combined[s])->Integral(bin187,bin227) ) / 2.;
  //      v_raw_pi0_count.push_back( raw_yield_count );
  //
  //
  //      /* END side band subtraction */
  //
  //      /* RE-bin histogram for higher counts per pin for better fit stability */
  //      (h_inv_mass_combined[s])->Rebin(5);
  //
  //      /* Fit: Gaussian around peak region to get initial parameters for later fit */
  //      TF1 *f_fit_pre = new TF1("f_fit_pre","gaus");
  //
  //      (h_inv_mass_combined[s])->Fit(f_fit_pre, "Q", "", fitrange_min, fitrange_max);
  //
  //      /* Fit: Gaussian + background (polynomial) */
  //      TF1 *f_fit = new TF1("f_fit","[0]*exp(-0.5*((x-[1])/[2])**2)+[3]+[4]*x+[5]*x*x+[6]*x*x*x");
  //      f_fit->SetParameter(0, f_fit_pre->GetParameter(0) );
  //      f_fit->SetParameter(1, f_fit_pre->GetParameter(1) );
  //      f_fit->SetParameter(2, f_fit_pre->GetParameter(2) );
  //      f_fit->SetParameter(3,0.3);
  //      f_fit->SetParameter(4,0.3);
  //      f_fit->SetParameter(5,0.3);
  //      f_fit->SetParameter(6,0.3);
  //
  //      /* Fit twice because sometimes this helops ROOT improve fit result */
  //      (h_inv_mass_combined[s])->Fit(f_fit, "Q", "");
  //      (h_inv_mass_combined[s])->Fit(f_fit, "Q", "");
  //
  //      /* separate out background contribution */
  //      TF1 *f_fit_bg = new TF1("f_fit_bg","[0]+[1]*x+[2]*x*x+[3]*x*x*x");
  //      f_fit->SetParameter(0, f_fit->GetParameter(3) );
  //      f_fit->SetParameter(1, f_fit->GetParameter(4) );
  //      f_fit->SetParameter(2, f_fit->GetParameter(5) );
  //      f_fit->SetParameter(3, f_fit->GetParameter(6) );
  //
  //      /* Get pi0 yeild from fit minus background */
  //      double pi0_yield = f_fit->Integral(0,0.3) - f_fit_bg->Integral(0,0.3);
  //
  //      //      cout << "Pi0 yield: " << pi0_yield << endl;
  //
  //      v_sector.push_back( s );
  //      v_mean.push_back( f_fit->GetParameter(1) );
  //      v_dmean.push_back( f_fit->GetParError(1) );
  //      v_sigma.push_back( f_fit->GetParameter(2) );
  //      v_dsigma.push_back( f_fit->GetParError(2) );
  //      v_chisquare.push_back( f_fit->GetChisquare() );
  //    }
  //
  //  /* print results */
  //  cout << "Raw pi0 yield count (run, PbScW, PbScE, PbGl): " << run_number << " " << v_raw_pi0_count.at(0) << " " << v_raw_pi0_count.at(1) << " " << v_raw_pi0_count.at(2) << endl;
  //
  //  for ( unsigned int s = 0; s < 3; s++ )
  //    {
  //      cout << "pi0peak "
  //       << run_index << " "
  //       << run_number << " "
  //       << nevents << " "
  //       << v_sector.at( s ) << " "
  //       << v_mean.at( s ) << " "
  //       << v_dmean.at( s ) << " "
  //       << v_sigma.at( s ) << " "
  //       << v_dsigma.at( s ) << " "
  //       << v_chisquare.at( s )
  //       << endl;
  //    }
  //
  //  /* Visual output */
  //  if ( visualize )
  //    {
  //      //      gStyle->SetOptStat(0);
  //
  //      /* Plot energy spectrum */
  //      TCanvas *c1 = new TCanvas();
  //      h_inv_mass_combined[0]->Draw();
  //
  //      TCanvas *c2 = new TCanvas();
  //      h_inv_mass_combined[1]->Draw();
  //
  //      TCanvas *c3 = new TCanvas();
  //      h_inv_mass_combined[2]->Draw();
  //
  //      c1->Print("plots-multicollision/pi0_example_1run_PbScW.png");
  //      c2->Print("plots-multicollision/pi0_example_1run_PbScE.png");
  //      c3->Print("plots-multicollision/pi0_example_1run_PbGlE.png");
  //    }
  //
  //  return 0;
}
