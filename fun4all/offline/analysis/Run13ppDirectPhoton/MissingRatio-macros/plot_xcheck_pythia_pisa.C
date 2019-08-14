int
plot_xcheck_pythia_pisa()
{
  // Settings
  gStyle->SetOptStat(0);

  // Open files
  TFile *fpisa = new TFile("histos/pisa-anaphpyhtia-isolation-0-5kevents.root", "OPEN");
  TFile *fpythia = new TFile("histos/anaphpythia-isolation-0-5kevents.root", "OPEN");

  // Get trees
  TTree *tpisa2 = (TTree*)fpisa->Get("mcparticles");
  TTree *tpythia = (TTree*)fpythia->Get("event_truth");

  tpisa2->SetScanField(0);
  tpythia->SetScanField(0);

  vector<float> v_evtcount;
  vector<float> v_eta;
  vector<float> v_phi;
  vector<float> v_ptot;

  TString cut_1072("t_pid == 1072 && t_parentpid == 0");

  tpisa2->Draw("eventcounter:t_eta:t_phi:t_ptot",cut_1072,"Q");

  TH1F* h_pid_pythia = new TH1F("h_pid_pythia","",2001, -1000.5, 1000.5);
  h_pid_pythia->GetXaxis()->SetTitle("Pythia6 PID");
  h_pid_pythia->GetYaxis()->SetTitle("Counts");

  for ( unsigned i = 0; i < tpisa2->GetEntries(cut_1072); i++ )
    {
      cout << "PRocessing event " << i << endl;
      float event = tpisa2->GetV1()[i];
      float eta   = tpisa2->GetV2()[i];
      float phi   = tpisa2->GetV3()[i];
      float ptot  = tpisa2->GetV4()[i];

      TString cut_pythia = TString::Format("eventcounter == %f && fabs(t_eta - %f) < 0.00001 && fabs(t_phi - %f) < 0.00001 && fabs(t_ptot - %f) < 0.00001 ",
					   event, eta, phi, ptot);
      unsigned matches = tpythia->GetEntries( cut_pythia );

      if ( matches == 1 )
	{
	  tpythia->Draw("t_pid",cut_pythia,"Q");
	  h_pid_pythia->Fill( tpythia->GetV1()[0] );
	}
      else
	{
	  cout << "Event " << i << " -> " << matches << " matches!" << endl;
	  h_pid_pythia->Fill( -999 );
	}

      //      tpythia->Scan("eventcounter:t_eta:t_phi:t_ptot:t_parentpid:t_pid",cut_pythia);
    }

  h_pid_pythia->Draw();
//////  cout << "Set size: " << evts_check.size() << endl;
//////
//////  set<int>::iterator it;
//////  for (it = evts_check.begin(); it != evts_check.end(); ++it)
//////    {
//////      cout << "######## Event " << *it << " #########" << endl;
//////      TString evt_cut(" eventcounter == " );
//////      evt_cut+=*it;
//////
//////      cout << "PISA particles:" << endl;
//////      //      tpisa2->Scan("eventcounter:t_eta:t_phi:t_ptot:t_parentpid:t_pid",evt_cut && "t_pid==1072");
//////      tpisa2->Scan("eventcounter:t_eta:t_phi:t_ptot:t_parentpid:t_pid",evt_cut);
//////
//////      cout << "PHPythia particles:" << endl;
//////      tpythia->Scan("eventcounter:t_eta:t_phi:t_ptot:t_parentpid:t_pid",evt_cut);
//////    }
//////
//  /* use set to print event records */
//
//  /* get list of events that have a 'pid--1072' particle */
//  tpisa2->Draw("eventcounter>>h1(5001,-0.5,5000.5","t_pid==1072 && t_anclvl == 0");
//
//  unsigned ncheck = h1->GetEntries();
//
//  set<int> evts_check;
//
//  for ( unsigned i = 0; i < 100; i++ )
//    {
//      evts_check.insert( tpisa2->GetV1()[i] );
//    }
//
//  cout << "Set size: " << evts_check.size() << endl;
//
//  set<int>::iterator it;
//  for (it = evts_check.begin(); it != evts_check.end(); ++it)
//    {
//      cout << "######## Event " << *it << " #########" << endl;
//      TString evt_cut(" eventcounter == " );
//      evt_cut+=*it;
//
//      cout << "PISA particles:" << endl;
//      //      tpisa2->Scan("eventcounter:t_eta:t_phi:t_ptot:t_parentpid:t_pid",evt_cut && "t_pid==1072");
//      tpisa2->Scan("eventcounter:t_eta:t_phi:t_ptot:t_parentpid:t_pid",evt_cut);
//
//      cout << "PHPythia particles:" << endl;
//      tpythia->Scan("eventcounter:t_eta:t_phi:t_ptot:t_parentpid:t_pid",evt_cut);
//
//      cout << endl;
//      cout << endl;
//      cout << endl;
//    }

  return 0;
}
