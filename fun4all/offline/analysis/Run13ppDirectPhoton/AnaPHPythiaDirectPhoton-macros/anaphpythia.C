void anaphpythia(const int process = 0)
{
  gSystem->Load("libfun4allfuncs.so");	// framework + reco modules
  gSystem->Load("libAnaPHPythiaDirectPhoton.so");

  const int nThread = 10;
  char dstFileName[1000];

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",0);

  // Server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Reconstruction Modules
  SubsysReco *my1 = new AnaPHPythiaHistos("AnaPHPythiaHistos", Form("histo%d.root",process));
  se->registerSubsystem(my1);

  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // Loop over input DST files
  for(int thread = process*nThread; thread < (process+1)*nThread; thread++)
  {
    sprintf(dstFileName, "/phenix/spin/phnxsp01/zji/data/pisaRun13/phpythia-isodecay/phpythia%d.root", thread);

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

  delete se;
}
