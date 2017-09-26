print_TH2_ybins()
{

  TFile *fin = new TFile("warnmap-data/WarnmapData_Run13pp510MinBias.root","OPEN");

  TH2F* hin = (TH2F*)fin->Get("hitmap_energy");

  for ( unsigned j = 1; j < hin->GetNbinsY(); j++ )
    {

      float jcenter = hin->GetYaxis()->GetBinCenter( j );
      float jwidth = hin->GetYaxis()->GetBinWidth( j );

      cout << "Bin " << j << " from " << jcenter - jwidth / 2. << " to " << jcenter + jwidth / 2. << endl;
    }

}
