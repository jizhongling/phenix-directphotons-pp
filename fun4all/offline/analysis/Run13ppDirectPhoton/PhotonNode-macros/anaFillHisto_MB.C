void anaFillHisto_MB(const int process = 0)
{
  gSystem->Load("libfun4all.so");	// framework + reco modules
  gSystem->Load("libPhotonNode.so");

  const int nThread = 10;
  int thread = -1;
  int runNumber;

  ifstream inFiles("/phenix/plhf/zji/taxi/Run13pp510MinBias/runnumber.txt");
  if(!inFiles)
  {
    cerr << "\nUnable to open input file list!" << endl;
    return;
  }
  cout << "\nUsing input files list..." << endl << endl;

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",0);

  // Server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // My Reconstruction Module
  FillHisto *my1 = new FillHisto("FillHisto_TAXI");
  my1->SelectMB();
  se->registerSubsystem(my1);

  // Real input from DST files
  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // Loop over input DST files
  while( inFiles >> runNumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    char dstFileName[1000];
    sprintf(dstFileName, "/phenix/spin/phnxsp01/zji/taxi/Run13pp510MinBias/12358/data/PhotonNode-%d.root", runNumber);

    cout << "\nfileopen for " << dstFileName << endl; 
    int openReturn = se->fileopen("DSTin1", dstFileName);
    if(openReturn)
    {
      cout << "\nAbnormal return: openReturn from Fun4All fileopen method = " << openReturn << endl;
      continue;
    }

    // Do the analysis for this DST file
    se->run(0);

    cout << "\nClosing input file, and a No Input file open message from Fun4All should appear" << endl;
    int closeReturn = se->fileclose("DSTin1");
    if(closeReturn)
      cout << "\nAbnormal return: closeReturn from Fun4All fileclose = " << closeReturn << endl;
  }

  // Write histograms
  se->End();

  // Delete Fun4All server
  delete se;
}
