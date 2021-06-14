void SpinDB(const int runnumber)
{
  gSystem->Load("libuspin.so");

  SpinDBOutput spin_out;
  SpinDBContent spin_cont;

  /* Initialize object to access spin DB */
  spin_out.Initialize();
  spin_out.SetUserName("phnxrc");
  spin_out.SetTableName("spin");

  /* Retrieve entry from Spin DB and get fill number */
  int qa_level = spin_out.GetDefaultQA(runnumber);
  spin_out.StoreDBContent(runnumber, runnumber, qa_level);
  spin_out.GetDBContentStore(spin_cont, runnumber);
  int fillnumber = spin_cont.GetFillNumber();

  /* Get spin info */
  if( spin_out.CheckRunRow(runnumber,qa_level) == 1 )
  {
    int runnumber = spin_cont.GetRunNumber();
    int qa_level = spin_cont.GetQALevel();
    int fillnumber = spin_cont.GetFillNumber();
    int badrunqa = spin_cont.GetBadRunFlag();
    int crossing_shift = spin_cont.GetCrossingShift();

    double pb, pbstat, pbsyst;
    double py, pystat, pysyst;

    spin_cont.GetPolarizationBlue(1, pb, pbstat, pbsyst);
    spin_cont.GetPolarizationYellow(1, py, pystat, pysyst);

    int badbunch[120];
    int spinpattern_blue[120];
    int spinpattern_yellow[120];

    long long bbc_narrow[120];
    long long bbc_wide[120];
    long long zdc_narrow[120];
    long long zdc_wide[120];

    for(int i=0; i<120; i++)
    {
      badbunch[i] = spin_cont.GetBadBunchFlag(i);
      spinpattern_blue[i] = spin_cont.GetSpinPatternBlue(i);
      spinpattern_yellow[i] = spin_cont.GetSpinPatternYellow(i);
      cout << spinpattern_blue[i]*spinpattern_yellow[i] << " ";

      bbc_narrow[i] = spin_cont.GetScalerBbcVertexCut(i);
      bbc_wide[i] = spin_cont.GetScalerBbcNoCut(i);
      zdc_narrow[i] = spin_cont.GetScalerZdcNarrow(i);
      zdc_wide[i] = spin_cont.GetScalerZdcWide(i);
    }
    cout << endl;
  }
}
