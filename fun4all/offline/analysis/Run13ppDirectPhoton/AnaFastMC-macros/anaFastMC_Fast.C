void anaFastMC_Fast(const int process = 0)
{
  // Set up Fun4All libraries
  gSystem->Load("libfun4allfuncs.so");	// framework + reco modules
  gSystem->Load("libAnaFastMC.so");

  // Switches for MCMethod
  enum MCMethod {PHParticleGen, FastMC};

  // Setup recoConsts
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",0);

  // Fun4All server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Reconstruction Modules
  AnaFastMC *my1 = new AnaFastMC("AnaFastMC");
  my1->set_outfile( Form("histos/AnaFastMC-Fast-histo%d.root",process) );
  my1->set_mcmethod(FastMC);

  // Register reconstruction modules
  se->registerSubsystem(my1);

  // A dummy (null) input is needed for the Fun4All framework
  Fun4AllDummyInputManager *in1 = new Fun4AllDummyInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // Run over events
  se->run(5000000);

  // Write histograms
  se->End();
  se->unregisterSubsystem(my1);

  // Delete Fun4All server
  delete se;
}
