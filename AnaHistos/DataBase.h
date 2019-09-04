class DataBase
{
  public:
    DataBase();

    unsigned long long GetClockLive(unsigned long long runnumber);
    unsigned long long GetBBCNovtxLive(unsigned long long runnumber);
    unsigned long long GetBBCNarrowLive(unsigned long long runnumber);
    unsigned long long GetBBCNarrowScaledown(unsigned long long runnumber);
    unsigned long long GetBBCNarrowScaledown(unsigned long long runnumber);
    unsigned long long GetERT4x4aScaledown(unsigned long long runnumber);
    unsigned long long GetERT4x4bScaledown(unsigned long long runnumber);
    unsigned long long GetERT4x4cScaledown(unsigned long long runnumber);

    int GetSpinPattern(int runnumber);
    void GetPol(int runnumber, double pb[], double py[]);
    void GetRelLum(int runnumber, double rlum[], double erlum[]);

  protected:
    static const int NRUN = 1020;
    static const int NVAL = 9;

    /* 0: Runnumber; 1: CLOCK live count;
     * 2: BBCLL1 novertex live count; 3: BBCLL1 narrowvtx live count;
     * 4: BBCLL1 novertex scaledown; 5: BBC narrowvtx scaledown;
     * 6-8: ERT_4x4a/b/c scaledown; */
    unsigned long long n_db[NRUN][NVAL];

    int spin_run[NRUN];  // runnumber from spin database
    char spin_pattern[NRUN][10];  // e.g. SSOO
    double pol_blue[NRUN][2];  // pb, pbstat
    double pol_yellow[NRUN][2];  // py, pystat
    double rel_lum[NRUN][2];  // N++/N+- for even and odd crossings
    double erel_lum[NRUN][2];  // uncertainty for rel_lum
};

void DataBase::DataBase()
{
  /* Initialization */
  for(int i=0; i<NRUN; i++)
  {
    for(int j=0; j<NVAL; j++)
      n_db[i][j] = 0;

    spin_run[i] = 0;
    strcpy(spin_pattern[i], "NONE");
    for(int j=0; j<2; j++)
    {
      pol_blue[i][j] = 0.;
      pol_yellow[i][j] = 0.;
      rel_lum[i][j] = 0.;
    }
  }

  TFile *fin = new TFile("data/clock-counts.root");
  TTree *t1 = (TTree*)fin->Get("t1");
  t1->AddFriend("t1f = t1", "data/SpinPattern.root");

  unsigned long long runno, clock, bbcnovtx_live, bbcnarrow_live;
  unsigned long long bbcnovtx_scaledown, bbcnarrow_scaledown;
  unsigned long long erta_scaledown, ertb_scaledown, ertc_scaledown;

  int runnumber;
  char pattern[10];
  double pb[3];  // pb, pbstat, pbsyst
  double py[3];  // py, pystat, pysyst
  double rlum[2];  // N++/N+- for even and odd crossings

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
  }

  delete fin;
}

unsigned long long DataBase::GetClockLive(unsigned long long runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(n_db[i][0] == runnumber)
      return n_db[i][1];

  return 0;
}

unsigned long long DataBase::GetBBCNovtxLive(unsigned long long runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(n_db[i][0] == runnumber)
      return n_db[i][2];

  return 0;
}

unsigned long long DataBase::GetBBCNarrowLive(unsigned long long runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(n_db[i][0] == runnumber)
      return n_db[i][3];

  return 0;
}

unsigned long long DataBase::GetBBCNovtxScaledown(unsigned long long runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(n_db[i][0] == runnumber)
      return n_db[i][4];

  return 0;
}

unsigned long long DataBase::GetBBCNarrowScaledown(unsigned long long runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(n_db[i][0] == runnumber)
      return n_db[i][5];

  return 0;
}

unsigned long long DataBase::GetERT4x4aScaledown(unsigned long long runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(n_db[i][0] == runnumber)
      return n_db[i][6];

  return 0;
}

unsigned long long DataBase::GetERT4x4bScaledown(unsigned long long runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(n_db[i][0] == runnumber)
      return n_db[i][7];

  return 0;
}

unsigned long long DataBase::GetERT4x4cScaledown(unsigned long long runnumber)
{
  for(int i=0; i<NRUN; i++)
    if(n_db[i][0] == runnumber)
      return n_db[i][8];

  return 0;
}

int DataBase::GetSpinPattern(int runnumber)
{
  const char *pattern_list[5] = {"SOOSSOO", "OSSOOSS", "SSOO", "OOSS", "NONE"};

  for(int i=0; i<NRUN; i++)
    if(spin_run[i] == runnumber)
    {
      for(int j=0; j<4; j++)
        if(strcmp(pattern_list[j], spin_pattern[i]) == 0)
          return j;

      return 4;
    }

  return 4;
}

void DataBase::GetPol(int runnumber, double pb[], double py[])
{
  for(int i=0; i<NRUN; i++)
    if(spin_run[i] == runnumber)
    {
      for(int j=0; j<2; j++)
      {
        pb[j] = pol_blue[i][j];
        py[j] = pol_yellow[i][j];
      }

      return;
    }

  return;
}

void DataBase::GetRelLum(int runnumber, double rlum[], double erlum[])
{
  for(int i=0; i<NRUN; i++)
    if(spin_run[i] == runnumber)
    {
      for(int j=0; j<2; j++)
      {
        rlum[j] = rel_lum[i][j];
        erlum[j] = erel_lum[i][j];
      }

      return;
    }

  return;
}
