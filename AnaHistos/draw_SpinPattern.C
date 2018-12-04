void draw_SpinPattern(const int process = 0)
{
  gSystem->Load("libuspin.so");

  SpinDBOutput spin_out;
  SpinDBContent spin_cont;

  /* Initialize object to access spin DB */
  spin_out.Initialize();
  spin_out.SetUserName("phnxrc");
  spin_out.SetTableName("spin");

  const int nThread = 1020;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runnumber.txt");

  const char *pattern_list[5] = {"SOOSSOO", "OSSOOSS", "SSOO", "OOSS", "NONE"};
  vector<int> assoc_run[5];
  for(int ipat=0; ipat<5; ipat++)
    assoc_run[ipat].clear();

  TFile *f_out = new TFile("data/SpinPattern.root", "RECREATE");
  TTree *t1 = new TTree("t1", "Spin pattern");
  char pattern[10];
  double pb[3];  // pb, pbstat, pbsyst
  double py[3];  // py, pystat, pysyst
  double rlum[2];  // N++/N+- for even and odd crossings
  double erlum[2];  // uncertainty for rlum
  t1->Branch("runnumber", &runnumber, "runnumber/I");
  t1->Branch("pattern", &pattern, "pattern/C");
  t1->Branch("blue_pol", &pb, "pb[3]/D");
  t1->Branch("yellow_pol", &py, "py[3]/D");
  t1->Branch("rel_lum", &rlum, "rlum[2]/D");
  t1->Branch("error_rel_lum", &erlum, "erlum[2]/D");

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    /* Retrieve entry from Spin DB */
    int qa_level = spin_out.GetDefaultQA(runnumber);
    spin_out.StoreDBContent(runnumber, runnumber, qa_level);
    spin_out.GetDBContentStore(spin_cont, runnumber);

    spin_cont.GetPolarizationBlue(1, pb[0], pb[1], pb[2]);
    spin_cont.GetPolarizationYellow(1, py[0], py[1], py[2]);

    string pat[2];  // for even and odd crossings
    ULong64_t lum_same[2] = {};  // for even and odd crossings
    ULong64_t lum_opp[2] = {};  // for even and odd crossings

    if( spin_out.CheckRunRow(runnumber,qa_level) == 1 &&
        spin_cont.GetRunNumber() == runnumber )
      for(int ib=0; ib<120; ib++)
      {
        int blue = spin_cont.GetSpinPatternBlue(ib);
        int yellow = spin_cont.GetSpinPatternYellow(ib);
        int value = blue * yellow;
        ULong64_t bbc_narrow = spin_cont.GetScalerBbcNoCut(ib);
        char cond = 'N';
        if(value == 1)
        {
          cond = 'S';
          lum_same[ib%2] += bbc_narrow;
        }
        else if(value == -1)
        {
          cond = 'O';
          lum_opp[ib%2] += bbc_narrow;
        }
        else
          value = 0;
        pat[ib%2] += cond;
      }

    for(int itype=0; itype<2; itype++)
    {
      rlum[itype] = (double)lum_same[itype] / (double)lum_opp[itype];
      erlum[itype] = rlum[itype] * sqrt(1./(double)lum_same[itype] + 1./(double)lum_opp[itype]);
    }

    size_t pattern_start[4] = {};
    for(int ipat=0; ipat<4; ipat++)
      pattern_start[ipat] = pat[0].find(pattern_list[ipat]);

    int assoc_pattern = 4;
    size_t min_start = string::npos;
    for(int ipat=0; ipat<4; ipat++)
      if(pattern_start[ipat] < min_start)
      {
        assoc_pattern = ipat;
        min_start = pattern_start[ipat];
      }

    if(min_start != string::npos)
      assoc_run[assoc_pattern].push_back(runnumber);
    strcpy(pattern, pattern_list[assoc_pattern]);

    cout << runnumber << ": " << pat[0] << "\t";
    for(int ipat=0; ipat<4; ipat++)
      cout << pattern_start[ipat] << " ";
    cout << "\n" << runnumber << ": " << pat[1] << endl;
    t1->Fill();
  }

  for(int ipat=0; ipat<4; ipat++)
  {
    size_t ngroup = assoc_run[ipat].size();
    cout << pattern_list[ipat] << "(" << ngroup <<"): ";
    for(size_t i=0; i<ngroup; i++)
      cout << assoc_run[ipat].at(i) << " ";
    cout << endl;
  }
  t1->Write();
  f_out->Close();
}
