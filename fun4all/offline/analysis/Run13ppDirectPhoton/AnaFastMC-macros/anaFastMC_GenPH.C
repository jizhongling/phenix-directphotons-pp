void anaFastMC_GenPH(const int process = 0,
    const int scale = 40,
    const char *outputname = "phpythia.root",
    const char *histoname = "histo.root")
{
  //gSystem->Load("libfun4allfuncs.so");	// framework only
  gSystem->Load("libfun4all.so");	// framework + reco modules
  gSystem->Load("libPHPythiaEventGen.so");
  gSystem->Load("libPHPythia.so");
  gSystem->Load("libPHPyTrigger.so");		// For PHPyTrigger derived classes
  gSystem->Load("libPHPyParticleSelect.so");	// For PHPyParticleSelect derived classes
  gSystem->Load("libsimreco.so");	// framework + reco modules
  gSystem->Load("libAnaFastMC.so");

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",0);

  /////////////////////////////////////////////////////////////////
  //  Server...
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  /////////////////////////////////////////////////////////////////
  //  Reconstruction Modules...

  SubsysReco *sync = new SyncSimreco();
  se->registerSubsystem(sync);

  PHPythia *phpythia = new PHPythia();

  // Set your own seed, otherwise, seeds from /dev/random
  //phpythia->SetSeed(1999);			

  se->registerSubsystem(phpythia);

  //** You can force the generated particles to use a vertex read from a file,
  //** in place of the default (z=0) value
  //** this is needed for instance when you want to have matching vertices between 
  //** different types of simulated files, prior to sending that to PISA
  // se->registerSubsystem( new PHPyVertexShift( "PHPyVertexShift", "./events.txt") );

  //** You can use dedicated triggers, derived from the PHPyTrigger base class
  // se->registerSubsystem( new PHPyJPsiMuonTrigger() );

  // My Reconstruction Module
  enum MCMethod {FastMC, PHParticleGen};
  AnaFastMC *my1 = new AnaFastMC("AnaFastMC");
  my1->set_outfile(histoname);
  my1->set_mcmethod(PHParticleGen);
  my1->use_xsec_weight();
  se->registerSubsystem(my1);

  //** You can select only particular particles to write out
  PHPyParticleSelect *pselect = new PHPySelectStable();
  se->registerSubsystem( pselect );

  //** A dummy (null) input is needed for the Fun4All framework
  Fun4AllDummyInputManager *in1 = new Fun4AllDummyInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // DST output manager
  TString OUTPUT = outputname;
  Fun4AllDstOutputManager *dst_output_mgr  = new Fun4AllDstOutputManager("PHPYTHIA",OUTPUT.Data());
  dst_output_mgr->AddNode("Sync");
  dst_output_mgr->AddNode("PHPythiaHeader");
  dst_output_mgr->AddNode("PHPythia");
  se->registerOutputManager(dst_output_mgr);

  // OSCAR output manager
  // with following output manager, one can write the PHPythia output in an oscar formated output text file
  // PHPyOscarOutputManager *oscar_manager  = new PHPyOscarOutputManager( "OSCAR", oscar_outputname );
  // se->registerOutputManager(oscar_manager);

  // Run over all events
  //double pt_start = 3. + process/scale * 0.1;
  //PtWeights *ptweights = new PtWeights();
  //double weight_pythia = ptweights->Integral(pt_start, pt_start+1., "MinBias");
  //my1->SetWeightPythia(weight_pythia);
  se->run(10000);

  // Write histograms
  se->End();
}
