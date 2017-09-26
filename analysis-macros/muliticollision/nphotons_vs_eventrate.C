nphotons_vs_eventrate()
{
  // see Run13 pi0 cross section AN-1206

  gStyle->SetOptStat(0);

  TTree *trates = new TTree();
  trates->ReadFile("rates_MB.txt","run_index/F:run_number:nevents:sector:");

  // plot: # direct photons / # MB events vs Live BBC_LL1_narrow / clock trigger rate



}
