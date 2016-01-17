
extract_cluster_per_run( TString infilename )
{

  TFile *fin = new TFile(infilename, "OPEN");

  //  TH1I* evt_counter = (TH1I*)fin->Get("evt_counter");
  TH1I* cluster_counter = (TH1I*)fin->Get("cluster_counter");

  cout << "Cluster / event in run: " << cluster_counter->GetMean() << endl;

}
