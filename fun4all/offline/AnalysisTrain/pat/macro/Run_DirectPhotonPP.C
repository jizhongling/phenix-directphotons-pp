void Run_DirectPhotonPP(const char *outFile = "HISTOS.root")
{
  // Load all your necessary libraries here
  gSystem->Load("libDirectPhotonPP.so");

  // To get access to external files
  gSystem->Load("libTOAD");

  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  toad_loader->SetVerbosity(1);
  string file_tofmap = toad_loader->location("Run13pp510_EMC_TOF_Correction.root");
  string file_ecal_run = toad_loader->location("Run13pp_RunbyRun_Calib.dat");

  // EMCal (re-)calibration class
  EmcLocalRecalibrator *emclocal = new EmcLocalRecalibrator();
  emclocal->SetEnergyCorrectionFile( file_ecal_run );
  emclocal->SetTofCorrectionFile( file_tofmap );

  // Analysis module
  // Put your SubSysReco derived analysis class here
  // outFile can be passed in and used to open your output file
  // Note, however, you can do whatever you like with what is
  // passed in
  DirectPhotonPP *dp = new DirectPhotonPP(outFile);
  dp->SetEmcLocalRecalibrator( emclocal );
  //dp->anaSetRunList(file_runlist.c_str());
  //dp->anaSelectGAMMA(); // Select GAMMA or MB data
  //dp->anaSetPtmin(10.);

  // create server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  se->registerSubsystem( dp );

  // clean up
  delete toad_loader;
}

void
InputData(vector<string> &indata)
{
  // Put all the data objects that you need here
  indata.push_back("CNT");
  // indata.push_back("DST_SVX"); //this should be on for run12 data should be left commented out for run11
  return;
}
