void anaFastMC(const int process=0)
{
  // Set up Fun4All libraries
  gSystem->Load("libfun4allfuncs.so");	// framework + reco modules
  gSystem->Load("libAnaFastMC.so");

  // Switches for MCMethod and WarnMap
  enum MCMethod {PHParticleGen, FastMC};
  enum WarnMap {None, Nils, Sasha};

  // Used for input DST files and output files
  const int nThread = 20;
  char dstFileName[1000];
  char outFileName[1000];

  // Setup recoConsts
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", 390039);

  // Fun4All server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Reconstruction Modules
  sprintf(outFileName, "histos/AnaFastMC-PH-nowarn-histo%d.root", process);
  AnaFastMC *reco_ph_nowarn = new AnaFastMC();
  reco_ph_nowarn->set_outfile(outFileName);
  reco_ph_nowarn->set_mcmethod(PHParticleGen);
  reco_ph_nowarn->set_warnmap(None);

  sprintf(outFileName, "histos/AnaFastMC-PH-warn-histo%d.root", process);
  AnaFastMC *reco_ph_warn = new AnaFastMC();
  reco_ph_warn->set_outfile(outFileName);
  reco_ph_warn->set_mcmethod(PHParticleGen);
  reco_ph_warn->set_warnmap(Sasha);

  sprintf(outFileName, "histos/AnaFastMC-Fast-nowarn-histo%d.root", process);
  AnaFastMC *reco_fast_nowarn = new AnaFastMC();
  reco_fast_nowarn->set_outfile(outFileName);
  reco_fast_nowarn->set_mcmethod(FastMC);
  reco_fast_nowarn->set_warnmap(None);

  sprintf(outFileName, "histos/AnaFastMC-Fast-warn-histo%d.root", process);
  AnaFastMC *reco_fast_warn = new AnaFastMC();
  reco_fast_warn->set_outfile(outFileName);
  reco_fast_warn->set_mcmethod(FastMC);
  reco_fast_warn->set_warnmap(Sasha);

  // Register reconstruction modules for FastMC generater
  se->registerSubsystem(reco_fast_nowarn);
  se->registerSubsystem(reco_fast_warn);

  // A dummy (null) input is needed for the Fun4All framework
  Fun4AllDummyInputManager *dummy_in = new Fun4AllDummyInputManager("DST_dummy_in", "DST");
  se->registerInputManager(dummy_in);

  // Run over all events
  se->run(2000000);

  // Write histograms and unregister reconstruction modules for FastMC generator
  se->End();
  se->unregisterSubsystem(reco_fast_nowarn);
  se->unregisterSubsystem(reco_fast_warn);

  // Register reconstruction modules for PHParticleGen generater
  se->registerSubsystem(reco_ph_nowarn);
  se->registerSubsystem(reco_ph_warn);

  // Real input from DST files
  Fun4AllInputManager *real_in = new Fun4AllDstInputManager("DST_real_in", "DST");
  se->registerInputManager(real_in);

  // Loop over input DST files
  for(int thread=process*nThread; thread<process*nThread+nThread; thread++)
  {
    sprintf(dstFileName, "/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/phparticlegen/phparticlegen%d.root", thread);

    cout << "\nfileopen for " << dstFileName << endl; 
    int openReturn = se->fileopen("DST_real_in", dstFileName);
    if(openReturn)
    {
      cout << "\nAbnormal return: openReturn from Fun4All fileopen method = " << openReturn << endl;
      continue;
    }

    // Do the analysis for this DST file
    se->run(0);

    cout << "\nClosing input file, and a No Input file open message from Fun4All should appear" << endl;
    int closeReturn = se->fileclose("DST_real_in");
    if(closeReturn)
      cout << "\nAbnormal return: closeReturn from Fun4All fileclose = " << closeReturn << endl;
  }

  // Write histograms and unregister reconstruction modules for PHParticleGen generator
  se->End();
  se->unregisterSubsystem(reco_fast_nowarn);
  se->unregisterSubsystem(reco_fast_warn);

  // Delete Fun4All server
  delete se;
}
