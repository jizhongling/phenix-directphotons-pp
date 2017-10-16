void anaFillHisto(const int process=64)
{
  // Set up Fun4All libraries
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("libPhotonNode.so");

  const int nThread = 10;
  int thread = -1;
  int runNumber;
  char dstFileName[1000];

  //ifstream inFiles("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist.txt");
  ifstream inFiles("/phenix/plhf/zji/taxi/Run13pp510MinBias/runnumber.txt");
  if(!inFiles)
  {
    cerr << "\nUnable to open input file list!" << endl;
    return;
  }
  cout << "\nUsing input files list..." << endl << endl;

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Reconstruction Module
  //SubsysReco *my1 = new FillHisto("FILLHISTO", Form("histo%d.root",process));
  SubsysReco *my1 = new FillHistoMB("FILLHISTOMB", Form("histo%d.root",process));
  se->registerSubsystem(my1);

  // Input Manager
  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // Loop over input DST files
  while( inFiles >> runNumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    //sprintf(dstFileName, "/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/11465/data/DirectPhotonPP_PhotonNode--%d.root", runNumber);
    sprintf(dstFileName, "/phenix/plhf/zji/taxi/Run13pp510MinBias/11343/data/DirectPhotonPP_PhotonNode-%d.root", runNumber);

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

  // Write out the histogram file
  se->End();

  delete my1;
  delete se;
  inFiles.close();
}
