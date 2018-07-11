void anaFastMC_PH(const int process = 0)
{
  // Set up Fun4All libraries
  gSystem->Load("libfun4allfuncs.so");	// framework + reco modules
  gSystem->Load("libAnaFastMC.so");

  // Switches for MCMethod
  enum MCMethod {PHParticleGen, FastMC};

  // Used for input DST files
  const int nThread = 10;
  char dstFileName[1000];

  // Setup recoConsts
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",0);

  // Fun4All server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Reconstruction Modules
  AnaFastMC *my1 = new AnaFastMC("AnaFastMC");
  my1->set_outfile( Form("histos/AnaFastMC-PH-histo%d.root",process) );
  my1->set_mcmethod(PHParticleGen);

  // Register reconstruction modules
  se->registerSubsystem(my1);

  // Real input from DST files
  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // Loop over input DST files
  for(int thread=process*nThread; thread<(process+1)*nThread; thread++)
  {
    sprintf(dstFileName, "/phenix/spin/phnxsp01/zji/data/pisaRun13/phpythia-ckin3/phpythia%d.root", thread);

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
  se->unregisterSubsystem(my1);

  // Delete Fun4All server
  delete se;
}
