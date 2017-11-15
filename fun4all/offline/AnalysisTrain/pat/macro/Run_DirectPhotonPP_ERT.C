void Run_DirectPhotonPP_ERT(const char *outFile = "HISTOS.root")
{
  /* Load all your necessary libraries here
   */
  gSystem->Load("libDirectPhotonPP.so");

  /* To get access to external files
   */
  gSystem->Load("libTOAD");

  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  toad_loader->SetVerbosity(1);

  /* Analysis module
   * Put your SubSysReco derived analysis class here
   * outFile can be passed in and used to open your output file
   * Note, however, you can do whatever you like with what is
   * passed in
   */
  DirectPhotonPP *dp = new DirectPhotonPP(outFile);

  /* Set which type of DST is being used- MinBias or ERT
   */
  dp->SetDstDataType("ERT");

  /* Set minimum energy threshold for direct photons
   */
  dp->SetDirectPhotonEnergyMin( 5.0 );

  /* You can only register ONE EmcLocalRecalibrator or the module will refuse to run.
   */

  /* EMCal (re-)calibration class
   */
  //EmcLocalRecalibrator *emclocal = new EmcLocalRecalibrator();
  //string file_tofmap = toad_loader->location("Run13pp510_EMC_TOF_Correction.root");
  //string file_ecal_run = toad_loader->location("Run13pp_RunbyRun_Calib.dat");
  //emclocal->SetEnergyCorrectionFile( file_ecal_run );
  //emclocal->SetTofCorrectionFile( file_tofmap );
  //dp->SetEmcLocalRecalibrator( emclocal );

  /* EMCal (re-)calibration class (Sasha calibration)
   */
  EmcLocalRecalibratorSasha *emcrecalib_sasha = new EmcLocalRecalibratorSasha();
  string file_ecal = toad_loader->location("ecorr_run13pp500gev.txt");
  string file_ecal_run = toad_loader->location("ecorr_run_run13pp500gev.txt");
  string file_tcal = toad_loader->location("tcorr_run13pp500gev.txt");
  emcrecalib_sasha->anaGetCorrCal( file_ecal.c_str() );
  emcrecalib_sasha->anaGetCorrCal_run( file_ecal_run.c_str() );
  emcrecalib_sasha->anaGetCorrTof( file_tcal.c_str() );
  dp->SetEmcLocalRecalibrator( emcrecalib_sasha );

  /* Select which warn map file to read (filename, number of columns)
   */
  //dp->ReadTowerStatus( "Warnmap_Run13pp510.txt", 4 );
  dp->ReadTowerStatus( "warn_all_run13pp500gev.dat", 2 );

  //dp->SetClusterDebugMode(true);
  //dp->anaSetRunList(file_runlist.c_str());
  //dp->anaSelectGAMMA(); // Select GAMMA or MB data
  //dp->anaSetPtmin(10.);

  /* create server
   */
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  se->registerSubsystem( dp );

  /* clean up
   */
  delete toad_loader;

  cout << "All done." << endl;
}

void
InputData(vector<string> &indata)
{
  /* Put all the data objects that you need here
   */
  indata.push_back("CNT");
  // indata.push_back("DST_SVX"); //this should be on for run12 data should be left commented out for run11
  return;
}
