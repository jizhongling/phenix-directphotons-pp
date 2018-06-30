void draw_SpinPattern()
{
  gSystem->Load("libuspin.so");

  const int runnumber = 391870;

  SpinDBOutput spin_out;
  SpinDBContent spin_cont;

  /* Initialize object to access spin DB */
  spin_out.Initialize();
  spin_out.SetUserName("phnxrc");
  spin_out.SetTableName("spin");

  /* Retrieve entry from Spin DB */
  int qa_level = spin_out.GetDefaultQA(runnumber);
  spin_out.StoreDBContent(runnumber, runnumber, qa_level);
  spin_out.GetDBContentStore(spin_cont, runnumber);

  string pattern;
  if( spin_out.CheckRunRow(runnumber,qa_level) == 1 &&
      spin_cont.GetRunNumber() == runnumber )
    for(int ib=0; ib<120; ib++)
    {
      int blue = spin_cont.GetSpinPatternBlue(ib);
      int yellow = spin_cont.GetSpinPatternYellow(ib);
      int value = blue * yellow;
      if(value == 1)  pattern += "S";
      else if(value==-1) pattern += "O";
      else pattern += "N";
    }
  cout << pattern << endl;

  vector<int> SOOSSOO;
  vector<int> OSSOOSS;
  vector<int> SSOO;
  vector<int> OOSS;

  size_t found = -1;
  while(true)
  {
    found = pattern.find("SOOSSOO", found+1);
    if(found != string::npos)
    {
      SOOSSOO.push_back(found);
      continue;
    }
    break;
  }

  found = -1;
  while(true)
  {
    found = pattern.find("OSSOOSS", found+1);
    if(found != string::npos)
    {
      OSSOOSS.push_back(found);
      continue;
    }
    break;
  }

  found = -1;
  while(true)
  {
    found = pattern.find("SSOO", found+1);
    if(found != string::npos)
    {
      SSOO.push_back(found);
      continue;
    }
    break;
  }

  found = -1;
  while(true)
  {
    found = pattern.find("OOSS", found+1);
    if(found != string::npos)
    {
      OOSS.push_back(found);
      continue;
    }
    break;
  }

  vector<int>::iterator ip;
  for( ip = SSOO.begin(); ip != SSOO.end(); ip++ )
    cout << *ip << " ";
  cout << endl;
  for( ip = OOSS.begin(); ip != OOSS.end(); ip++ )
    cout << *ip << " ";
  cout << endl;
}
