void draw_SpinPattern()
{
  gSystem->Load("libDirectPhotonPP.so");

  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist.txt");
  Int_t nrun = 0;
  Int_t runnumber[1024];
  while(fin >> runnumber[nrun]) nrun++;
  fin.close();

  TFile *f = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510ERT/8847/data/DirectPhotonPP_PhotonNode-%d.root", runnumber[600]));
  TTree *T1 = (TTree*)f->Get("T1");

  SpinPattern *spinpat = new SpinPattern;
  T1->SetBranchAddress("RUN/SpinPattern", &spinpat);
  T1->GetEntry(0);

  string pattern;
  for(Int_t ib=0; ib<120; ib++)
  {
    int blue = spinpat->get_spinpattern_blue(ib);
    int yellow = spinpat->get_spinpattern_yellow(ib);
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
