
int print_warnmap_summary_paul( )
{
  gStyle->SetOptStat(0);

  string warnfile1="warnmap-paul/iter10_rms10/warn_Run9pp200MinBias.txt";
  string warnfile2="warnmap-paul/iter10_rms10/warn_Run9pp500MinBias.txt";
  string warnfile3="warnmap-paul/iter10_rms10/warn_Run9pp200ERT.txt";
  string warnfile4="warnmap-paul/iter10_rms10/warn_Run9pp500ERT.txt";

  /* Read warnmap as tree */
  TTree *twarn1 = new TTree();
  twarn1->ReadFile( warnfile1.c_str(), "towerid/I:status" );

  TTree *twarn2 = new TTree();
  twarn2->ReadFile( warnfile2.c_str(), "towerid/I:status" );

  TTree *twarn3 = new TTree();
  twarn3->ReadFile( warnfile3.c_str(), "towerid/I:status" );

  TTree *twarn4 = new TTree();
  twarn4->ReadFile( warnfile4.c_str(), "towerid/I:status" );

  /* Count total towers */
  int count_total1 = twarn1->GetEntries();
  int count_total2 = twarn2->GetEntries();
  int count_total3 = twarn3->GetEntries();
  int count_total4 = twarn4->GetEntries();

  /* Count live towers */
  int count_live1 = twarn1->GetEntries("status <= 10");
  int count_live2 = twarn2->GetEntries("status <= 10");
  int count_live3 = twarn3->GetEntries("status <= 10");
  int count_live4 = twarn4->GetEntries("status <= 10");


  /* Print summary */
  cout << warnfile1 << ": " << count_live1 << " out of " << count_total1 << " towers alive = " << (float)count_live1 / (float)count_total1 << endl;
  cout << warnfile2 << ": " << count_live2 << " out of " << count_total2 << " towers alive = " << (float)count_live2 / (float)count_total2 << endl;
  cout << warnfile3 << ": " << count_live3 << " out of " << count_total3 << " towers alive = " << (float)count_live3 / (float)count_total3 << endl;
  cout << warnfile4 << ": " << count_live4 << " out of " << count_total4 << " towers alive = " << (float)count_live4 / (float)count_total4 << endl;

  return 0;
}
