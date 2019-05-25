void anaFastMC_PH(const int process = 0)
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

  //** You can select only particular particles to write out
  //PHPyParticleSelect *pselect = new PHPyParticleSelect();
  //se->registerSubsystem( pselect );

  // My Reconstruction Module
  enum MCMethod {FastMC, PHParticleGen};
  AnaFastMC *my1 = new AnaFastMC( Form("AnaFastMC%d",process) );
  my1->set_outfile( Form("../AnaFastMC-PH-histo%d.root",process) );
  my1->set_mcmethod(PHParticleGen);
  se->registerSubsystem(my1);

  // Real input from DST files
  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  // DST output manager
  //TString OUTPUT = outputname;
  //Fun4AllDstOutputManager *dst_output_mgr  = new Fun4AllDstOutputManager("PHPYTHIA",OUTPUT.Data());
  //dst_output_mgr->AddNode("Sync");
  //dst_output_mgr->AddNode("PHPythiaHeader");
  //dst_output_mgr->AddNode("PHPythia");
  //se->registerOutputManager(dst_output_mgr);

  // OSCAR output manager
  // with following output manager, one can write the PHPythia output in an oscar formated output text file
  // PHPyOscarOutputManager *oscar_manager  = new PHPyOscarOutputManager( "OSCAR", oscar_outputname );
  // se->registerOutputManager(oscar_manager);

  // Input DST file
  char dstFileName[1000];
  sprintf(dstFileName, "/phenix/spin/phnxsp01/zji/data/pisaRun13/phpythia-BG/phpythia%d.root", process);

  cout << "\nfileopen for " << dstFileName << endl; 
  int openReturn = se->fileopen("DSTin1", dstFileName);
  if(openReturn)
    cout << "\nAbnormal return: openReturn from Fun4All fileopen method = " << openReturn << endl;

  // Do the analysis for this DST file
  se->run(0);

  cout << "\nClosing input file, and a No Input file open message from Fun4All should appear" << endl;
  int closeReturn = se->fileclose("DSTin1");
  if(closeReturn)
    cout << "\nAbnormal return: closeReturn from Fun4All fileclose = " << closeReturn << endl;

  // Write histograms
  se->End();
}
