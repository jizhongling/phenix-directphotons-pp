void anaFastMC_Fast(const int process = 0)
{
  gSystem->Load("libfun4all.so");	// framework + reco modules
  gSystem->Load("libAnaFastMC.so");

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",0);

  // Server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // My Reconstruction Module
  enum MCMethod {FastMC, PHParticleGen};
  AnaFastMC *my1 = new AnaFastMC("AnaFastMC");
  my1->set_outfile( Form("histos/AnaFastMC-Fast-histo%d.root",process) );
  my1->set_mcmethod(FastMC);
  my1->enable_calcsys();
  se->registerSubsystem(my1);

  // A dummy (null) input is needed for the Fun4All framework
  Fun4AllDummyInputManager *in1 = new Fun4AllDummyInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // Run over all events
  se->run(5000000);

  // Write histograms
  se->End();

  // Delete Fun4All server
  delete se;
}
