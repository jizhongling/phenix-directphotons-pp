
extract_events_per_run( TString infilename )
{

  TFile *fin = new TFile(infilename, "OPEN");

  TH1I* evt_counter = (TH1I*)fin->Get("evt_counter");

  cout << "Events in run: " << evt_counter->GetBinContent(1) << endl;

}
