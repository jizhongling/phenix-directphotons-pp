// 0: Runnumber; 1: CLOCK live trigger count; 2: BBC novtx live count; 3: BBC narrow live count;
// 4: BBC novtx scaledown; 5: BBC narrow scaledown
unsigned long long n_clock_bbc[6][1020];

void ReadClockCounts()
{
  // initialize the CLOCK Live Counts and BBC narrow Live counts
  for(int i=0; i<6; i++)
    for(int j=0; j<1020; j++)
      n_clock_bbc[i][j] = 0;

  TFile *fin = new TFile("/phenix/plhf/zji/install/share/PhotonNode/clock-counts.root");
  TTree *t1 = (TTree*)fin->Get("t1");
  long long runno, clock, bbcnovtx_live, bbcnarrow_live;
  int bbcnovtx_scaledown, bbcnarrow_scaledown;
  t1->SetBranchAddress("runnumber", &runno);
  t1->SetBranchAddress("clock_live", &clock);
  t1->SetBranchAddress("bbcnovtx_live", &bbcnovtx_live);
  t1->SetBranchAddress("bbcnarrow_live", &bbcnarrow_live);
  t1->SetBranchAddress("bbcnovtx_scaledown", &bbcnovtx_scaledown);
  t1->SetBranchAddress("bbcnarrow_scaledown", &bbcnarrow_scaledown);
  int nentries = t1->GetEntries();
  for(int i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    n_clock_bbc[0][i] = runno;
    n_clock_bbc[1][i] = clock;
    n_clock_bbc[2][i] = bbcnovtx_live;
    n_clock_bbc[3][i] = bbcnarrow_live;
    n_clock_bbc[4][i] = bbcnovtx_scaledown;
    n_clock_bbc[5][i] = bbcnarrow_scaledown;
  }
  delete fin;

  return;
}

unsigned long long GetClockLive(unsigned runnumber)
{
  for(int i=0; i<1020; i++)
    if(runnumber == n_clock_bbc[0][i])
      return n_clock_bbc[1][i];

  return 0;
}

unsigned long long GetBBCNovtxLive(unsigned runnumber)
{
  for(int i=0; i<1020; i++)
    if(runnumber == n_clock_bbc[0][i])
      return n_clock_bbc[2][i];

  return 0;
}

unsigned long long GetBBCNarrowLive(unsigned runnumber)
{
  for(int i=0; i<1020; i++)
    if(runnumber == n_clock_bbc[0][i])
      return n_clock_bbc[3][i];

  return 0;
}

unsigned long GetBBCNovtxScaledown(unsigned runnumber)
{
  for(int i=0; i<1020; i++)
    if(runnumber == n_clock_bbc[0][i])
      return n_clock_bbc[4][i];

  return 0;
}

unsigned long GetBBCNarrowScaledown(unsigned runnumber)
{
  for(int i=0; i<1020; i++)
    if(runnumber == n_clock_bbc[0][i])
      return n_clock_bbc[5][i];

  return 0;
}
