// 0: Runnumber; 1: CLOCK live count;
// 2: BBCLL1 novertex live count; 3: BBCLL1 narrowvtx live count;
// 4: BBCLL1 novertex scaledown; 5: BBC narrowvtx scaledown;
// 6-8: ERT_4x4a/b/c scaledown;

#define NVAL 9
#define NRUN 1020
ULong64_t n_db[NVAL][NRUN];

void ReadClockCounts()
{
  // initialization
  for(int i=0; i<NVAL; i++)
    for(int j=0; j<NRUN; j++)
      n_db[i][j] = 0;

  TFile *fin = new TFile("data/clock-counts.root");
  TTree *t1 = (TTree*)fin->Get("t1");

  ULong64_t runno, clock, bbcnovtx_live, bbcnarrow_live;
  ULong64_t bbcnovtx_scaledown, bbcnarrow_scaledown;
  ULong64_t erta_scaledown, ertb_scaledown, ertc_scaledown;

  t1->SetBranchAddress("runnumber", &runno);
  t1->SetBranchAddress("clock_live", &clock);
  t1->SetBranchAddress("bbcnovtx_live", &bbcnovtx_live);
  t1->SetBranchAddress("bbcnarrow_live", &bbcnarrow_live);
  t1->SetBranchAddress("bbcnovtx_scaledown", &bbcnovtx_scaledown);
  t1->SetBranchAddress("bbcnarrow_scaledown", &bbcnarrow_scaledown);
  t1->SetBranchAddress("erta_scaledown", &erta_scaledown);
  t1->SetBranchAddress("ertb_scaledown", &ertb_scaledown);
  t1->SetBranchAddress("ertc_scaledown", &ertc_scaledown);

  int nentries = t1->GetEntries();
  for(int i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    n_db[0][i] = runno;
    n_db[1][i] = clock;
    n_db[2][i] = bbcnovtx_live;
    n_db[3][i] = bbcnarrow_live;
    n_db[4][i] = bbcnovtx_scaledown;
    n_db[5][i] = bbcnarrow_scaledown;
    n_db[6][i] = erta_scaledown;
    n_db[7][i] = ertb_scaledown;
    n_db[8][i] = ertc_scaledown;
  }

  delete fin;
  return;
}

ULong64_t GetClockLive(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[0][i])
      return n_db[1][i];

  return 0;
}

ULong64_t GetBBCNovtxLive(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[0][i])
      return n_db[2][i];

  return 0;
}

ULong64_t GetBBCNarrowLive(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[0][i])
      return n_db[3][i];

  return 0;
}

ULong64_t GetBBCNovtxScaledown(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[0][i])
      return n_db[4][i];

  return 0;
}

ULong64_t GetBBCNarrowScaledown(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[0][i])
      return n_db[5][i];

  return 0;
}

ULong64_t GetERT4x4aScaledown(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[0][i])
      return n_db[6][i];

  return 0;
}

ULong64_t GetERT4x4bScaledown(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[0][i])
      return n_db[7][i];

  return 0;
}

ULong64_t GetERT4x4cScaledown(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[0][i])
      return n_db[8][i];

  return 0;
}
