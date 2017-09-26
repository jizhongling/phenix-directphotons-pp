const int _nhistos = 5;

merge_energyranges( string infile, string outfile )
{
  TFile *fin = new TFile( infile.c_str() , "OPEN" );

  TH1I* histarray[_nhistos];

  vector< string > histnames;
  histnames.push_back("hitmap_erange_0");
  histnames.push_back("hitmap_erange_1");
  histnames.push_back("hitmap_erange_2");
  histnames.push_back("hitmap_erange_3");
  histnames.push_back("hitmap_erange_4");

  for ( unsigned i = 0; i < histnames.size(); i++ )
    {
      histarray[i] = (TH1I*) fin->Get( histnames.at(i).c_str() );
      cout << "Read entries: " << histarray[i]->GetEntries() << endl;
    }

  cout << "\n";

  TFile *fout = new TFile( outfile.c_str() , "RECREATE" );
  TH1I* hist_sum_a = histarray[0]->Clone("hitmap_erange_0_to_4");
  TH1I* hist_sum_b = histarray[0]->Clone("hitmap_erange_0_to_2");
  TH1I* hist_sum_c = histarray[3]->Clone("hitmap_erange_3_to_4");

  for ( unsigned i = 1; i < histnames.size(); i++ )
    {
      hist_sum_a->Add(histarray[i]);
      cout << "Summed entries: " << hist_sum_a->GetEntries() << endl;

      if ( i < 3 )
	hist_sum_b->Add(histarray[i]);
      else
	hist_sum_c->Add(histarray[i]);
    }

  /* Write histograms to output file */
  fout->cd();
  hist_sum_a->Write();
  hist_sum_b->Write();
  hist_sum_c->Write();
  fout->Close();
  fin->Close();
}
