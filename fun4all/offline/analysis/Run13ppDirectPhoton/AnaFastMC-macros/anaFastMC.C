void anaFastMC(const char *filelist="phparticlegen.txt", const char *outfilename = "histo.root")
{
  // Set up input file location
  gSystem->Load("libFROG.so");
  FROG fr;

  // Set up Fun4All libraries
  gSystem->Load("libfun4allfuncs.so");	// framework + reco modules
  gSystem->Load("libAnaFastMC.so");

  // Switches for MCMethod and WarnMap
  enum MCMethod {PHParticleGen, FastMC};
  enum WarnMap {None, Nils, Sasha};

  // Used for input DST files
  char dstString[5000];
  char *dstFileName;

  // Output filename prefiex
  string outfile_prefix = "histos/AnaFastMC-";

  // Setup recoConsts
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", 390039);

  // Fun4All server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Reconstruction Modules
  AnaFastMC *reco_ph_nowarn = new AnaFastMC();
  reco_ph_nowarn->set_outfile(outfile_prefix + "PH-nowarn-" + outfilename);
  reco_ph_nowarn->set_mcmethod(PHParticleGen);
  reco_ph_nowarn->set_warnmap(None);

  AnaFastMC *reco_ph_warn = new AnaFastMC();
  reco_ph_warn->set_outfile(outfile_prefix + "PH-warn-" + outfilename);
  reco_ph_warn->set_mcmethod(PHParticleGen);
  reco_ph_warn->set_warnmap(Sasha);

  AnaFastMC *reco_fast_nowarn = new AnaFastMC();
  reco_fast_nowarn->set_outfile(outfile_prefix + "Fast-nowarn-" + outfilename);
  reco_fast_nowarn->set_mcmethod(FastMC);
  reco_fast_nowarn->set_warnmap(None);

  AnaFastMC *reco_fast_warn = new AnaFastMC();
  reco_fast_warn->set_outfile(outfile_prefix + "Fast-warn-" + outfilename);
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

  ifstream inFiles(filelist);
  if(!inFiles)
  {
    cerr << "\nUnable to open input file list " << filelist << endl;
    return;
  }
  cout << "\nUsing input files list " << filelist << endl << endl;

  // Loop over input DST files
  while( inFiles >> dstString )
  {
    dstFileName = fr.location(dstString);  // get FROG location
    if(!dstFileName)
    {
      cerr << "\nBad return from FROG, null dstFileName" << endl;
      return;
    }

    cout << "\nfileopen for " << dstFileName << endl; 
    int openReturn = se->fileopen("DST_real_in", dstFileName);
    if(openReturn)
      cout << "\nAbnormal return: openReturn from Fun4All fileopen method = " << openReturn << endl; 

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

  // Close input filelist and delete Fun4All server
  inFiles.close();
  delete se;
}
