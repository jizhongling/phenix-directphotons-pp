class DataBase
{
  public:
    DataBase();

    ULong64_t GetClockLive(ULong64_t runnumber);
    ULong64_t GetBBCNovtxLive(ULong64_t runnumber);
    ULong64_t GetBBCNarrowLive(ULong64_t runnumber);
    ULong64_t GetBBCNarrowScaledown(ULong64_t runnumber);
    ULong64_t GetBBCNarrowScaledown(ULong64_t runnumber);
    ULong64_t GetERT4x4aScaledown(ULong64_t runnumber);
    ULong64_t GetERT4x4bScaledown(ULong64_t runnumber);
    ULong64_t GetERT4x4cScaledown(ULong64_t runnumber);
    void GetSpinPattern(int runnumber, char pattern[]);
    void GetPol(int runnumber, double pb[], double py[]);
    void GetRelLum(int runnumber, double rlum[]);
    void GetPolPattern(int runnumber, int pol[]);

  protected:
    static const int NRUN = 1020;
    static const int NVAL = 9;
    static const char *pattern_list[5];

    /* 0: Runnumber; 1: CLOCK live count;
     * 2: BBCLL1 novertex live count; 3: BBCLL1 narrowvtx live count;
     * 4: BBCLL1 novertex scaledown; 5: BBC narrowvtx scaledown;
     * 6-8: ERT_4x4a/b/c scaledown; */
    ULong64_t n_db[NRUN][NVAL];
    int spin_run[NRUN];  // runnumber from spin database
    char spin_pattern[NRUN][10];  // e.g. SSOO
    double pol_blue[NRUN][2];  // pb, pbstat
    double pol_yellow[NRUN][2];  // py, pystat
    double rel_lum[NRUN][2];  // N++/N+- for even and odd crossings
    int pol_pattern[NRUN][120];  // same: 1; opposite: -1; other: 0
};

char *DataBase::pattern_list[5] = {"SOOSSOO", "OSSOOSS", "SSOO", "OOSS", "NONE"};

void DataBase::DataBase()
{
  /* Initialization */
  for(int i=0; i<NRUN; i++)
  {
    spin_run[i] = 0;
    strcpy(spin_pattern[i], "NONE");
    for(int j=0; j<NVAL; j++)
      n_db[i][j] = 0;
  }

  TFile *fin = new TFile("data/clock-counts.root");
  TTree *t1 = (TTree*)fin->Get("t1");
  t1->AddFriend("t1f = t1", "data/SpinPattern.root");

  ULong64_t runno, clock, bbcnovtx_live, bbcnarrow_live;
  ULong64_t bbcnovtx_scaledown, bbcnarrow_scaledown;
  ULong64_t erta_scaledown, ertb_scaledown, ertc_scaledown;
  int runnumber;
  char pattern[10];
  double pb[3];  // pb, pbstat, pbsyst;
  double py[3];  // py, pystat, pysyst;
  double rlum[2];  // N++/N+- for even and odd crossings
  int pol[120];  // same: 1; opposite: -1; other: 0

  t1->SetBranchAddress("runnumber", &runno);
  t1->SetBranchAddress("clock_live", &clock);
  t1->SetBranchAddress("bbcnovtx_live", &bbcnovtx_live);
  t1->SetBranchAddress("bbcnarrow_live", &bbcnarrow_live);
  t1->SetBranchAddress("bbcnovtx_scaledown", &bbcnovtx_scaledown);
  t1->SetBranchAddress("bbcnarrow_scaledown", &bbcnarrow_scaledown);
  t1->SetBranchAddress("erta_scaledown", &erta_scaledown);
  t1->SetBranchAddress("ertb_scaledown", &ertb_scaledown);
  t1->SetBranchAddress("ertc_scaledown", &ertc_scaledown);
  t1->SetBranchAddress("t1f.runnumber", &runnumber);
  t1->SetBranchAddress("t1f.pattern", &pattern);
  t1->SetBranchAddress("t1f.blue_pol", &pb);
  t1->SetBranchAddress("t1f.yellow_pol", &py);
  t1->SetBranchAddress("t1f.rel_lum", &rlum);
  t1->SetBranchAddress("t1f.polarization", &pol);

  int nentries = t1->GetEntries();
  for(int i=0; i<nentries; i++)
  {
    t1->GetEntry(i);
    n_db[i][0] = runno;
    n_db[i][1] = clock;
    n_db[i][2] = bbcnovtx_live;
    n_db[i][3] = bbcnarrow_live;
    n_db[i][4] = bbcnovtx_scaledown;
    n_db[i][5] = bbcnarrow_scaledown;
    n_db[i][6] = erta_scaledown;
    n_db[i][7] = ertb_scaledown;
    n_db[i][8] = ertc_scaledown;
    spin_run[i] = runnumber;
    strcpy(spin_pattern[i], pattern);
    for(int j=0; j<2; j++)
    {
      pol_blue[i][j] = pb[j];
      pol_yellow[i][j] = py[j];
      rel_lum[i][j] = rlum[j];
    }
    for(int j=0; j<120; j++)
      pol_pattern[i][j] = pol[j];
  }

  delete fin;
}

ULong64_t DataBase::GetClockLive(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[i][0])
      return n_db[i][1];

  return 0;
}

ULong64_t DataBase::GetBBCNovtxLive(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[i][0])
      return n_db[i][2];

  return 0;
}

ULong64_t DataBase::GetBBCNarrowLive(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[i][0])
      return n_db[i][3];

  return 0;
}

ULong64_t DataBase::GetBBCNovtxScaledown(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[i][0])
      return n_db[i][4];

  return 0;
}

ULong64_t DataBase::GetBBCNarrowScaledown(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[i][0])
      return n_db[i][5];

  return 0;
}

ULong64_t DataBase::GetERT4x4aScaledown(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[i][0])
      return n_db[i][6];

  return 0;
}

ULong64_t DataBase::GetERT4x4bScaledown(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[i][0])
      return n_db[i][7];

  return 0;
}

ULong64_t DataBase::GetERT4x4cScaledown(ULong64_t runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == n_db[i][0])
      return n_db[i][8];

  return 0;
}

void DataBase::GetSpinPattern(int runnumber, char pattern[])
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == spin_run[i])
      strcpy(pattern, spin_pattern[i]);

  return;
}

void DataBase::GetPol(int runnumber, double pb[], double py[])
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == spin_run[i])
      for(int j=0; j<2; j++)
      {
        pb[j] = pol_blue[i][j];
        py[j] = pol_yellow[i][j];
      }

  return;
}

void DataBase::GetRelLum(int runnumber, double rlum[])
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == spin_run[i])
      for(int j=0; j<2; j++)
        rlum[j] = rel_lum[i][j];

  return;
}

void DataBase::GetPolPattern(int runnumber, int pol[])
{
  for(int i=0; i<NRUN; i++)
    if(runnumber == spin_run[i])
      for(int j=0; j<120; j++)
        pol[j] = pol_pattern[i][j];

  return;
}
