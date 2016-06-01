void GetRawAsymm(TFile *f, Int_t arm, Int_t pattern, Long64_t *yield)
{
  const Double_t BR[30] = {0.0444823,0.155126,0.189877,0.209775,0.227416,0.235346,0.242619,0.258304,0.254693,0.219828,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313};

  if(arm == 0)
  {
    const Int_t sector_low = 1;
    const Int_t sector_high = 4;
    //const Double_t MissR[30] = {1.17052, 1.02665, 1.06256, 0.869176, 0.860046, 0.771741, 0.742194, 0.687999, 0.632733, 0.792072, 0.813304, 0.699638, 0.774388, 0.657568, 0.558923, 0.732799, 0.649979, 0.54994, 0.572371, 0.476033, 0.524398, 0.658759, 0.725874, 0.477244, 0.837308, 1.13577, 1.1215, 1.69883, 1.80851, 1.98848};
    const Double_t MissR[30] = {0.576394, 1.5215, 2.09614, 1.55921, 1.98752, 1.47444, 1.54181, 1.73786, 1.36634, 1.40074, 1.27356, 1.14209, 1.11697, 1.03645, 0.886776, 1.00954, 0.928018, 0.729948, 0.753175, 0.634566, 0.928788, 1.14879, 1.18122, 0.905212, 1.30759, 1.91196, 2.2468, 3.98329, 5.55577, 7.05814};
    const Double_t TrigE[30] = {0.000164287,0.000226765,0.000554624,0.00132205,0.00273391,0.00653341,0.0205823,0.0836921,0.247432,0.506061,0.711566,0.897269,0.920601,0.926923,0.952096,0.963303,0.897059,0.916667,0.931034,1,0.930233,0.916667,1,1,1,1,1,1,1,1};
  }
  else if(arm == 1)
  {
    const Int_t sector_low = 5;
    const Int_t sector_high = 8;
    //const Double_t MissR[30] = {0.968893, 0.814013, 0.887797, 0.794183, 0.61442, 0.658933, 0.691201, 0.627664, 0.553573, 0.585781, 0.526425, 0.565531, 0.375501, 0.469886, 0.398121, 0.398228, 0.281594, 0.346736, 0.348726, 0.358025, 0.375395, 0.394456, 0.325099, 0.34044, 0.398916, 0.383546, 0.586689, 0.94373, 0.874969, 1.76884};
    const Double_t MissR[30] = {0.761121, 1.21671, 1.77948, 1.53485, 1.23296, 1.18563, 1.27982, 1.4427, 1.188, 1.22016, 0.984438, 0.974272, 0.872964, 0.72481, 0.726905, 0.57041, 0.421232, 0.448614, 0.453693, 0.491227, 0.692569, 0.688636, 0.597951, 0.62367, 0.800111, 0.857026, 1.14432, 1.86413, 2.59804, 6.26555};
    const Double_t TrigE[30] = {0.000175027,0.000228947,0.000544174,0.00116377,0.002314,0.00449421,0.0283686,0.0850631,0.195189,0.340642,0.491794,0.645161,0.753769,0.774536,0.794643,0.833333,0.837607,0.866667,0.897959,0.840909,0.860759,0.916667,1,1,1,1,1,1,1,1};
  }
  else
  {
    cerr << "Wrong arm " << arm << endl;
    exit(1);
  }

  if(pattern == 1)
  {
    const Int_t pattern_bin = 3;
  }
  else if(pattern == -1)
  {
    const Int_t pattern_bin = 1;
  }
  else
  {
    cout << "Irrelevant pattern " << pattern << endl;
    return;
  }

  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  hn_1photon->GetAxis(0)->SetRange(sector_low, sector_high);
  hn_1photon->GetAxis(2)->SetRange(pattern_bin, pattern_bin);
  TH1 *h_1photon = hn_1photon->Projection(1);

  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");

  TAxis *axis0 = hn_2photon->GetAxis(0);
  TAxis *axis1 = hn_2photon->GetAxis(1);
  TAxis *axis2 = hn_2photon->GetAxis(2);
  TAxis *axis3 = hn_2photon->GetAxis(3);

  axis0->SetRange(sector_low, sector_high);
  axis3->SetRange(pattern_bin, pattern_bin);

  Int_t bin047 = axis2->FindBin(0.047);
  Int_t bin097 = axis2->FindBin(0.097);
  Int_t bin112 = axis2->FindBin(0.112);
  Int_t bin162 = axis2->FindBin(0.162);
  Int_t bin177 = axis2->FindBin(0.177);
  Int_t bin227 = axis2->FindBin(0.227);

  for(Int_t ipt=0; ipt<30; ipt++)
  {
    axis1->SetRange(ipt+1,ipt+1);
    TH1 *h_inv_mass = hn_2photon->Projection(2);

    Double_t nphoton = h_1photon->GetBinContent(ipt+1);
    Double_t npair = 0.;
    Double_t nbgside = 0.;

    for(Int_t ib=bin047; ib<bin097; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    for(Int_t ib=bin177; ib<bin227; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    nbgside /= 2.;

    for(Int_t ib=bin112; ib<bin162; ib++)
      npair += h_inv_mass->GetBinContent(ib);

    Double_t nsub = npair - nbgside;
    Double_t npion = nsub * ( 1. + MissR[ipt] );
    Double_t ndecay = npion * ( 1. + BR[ipt] );
    Double_t ndirect = nphoton - ndecay;
    yield[ipt] = ndirect;

    delete h_inv_mass;
  }

  delete h_1photon;
  return;
}

void generate_ALL(Int_t Process = 0)
{
  gSystem->Load("libDirectPhotonPP.so");

  const Int_t npT = 30;
  const Float_t pTbins[npT+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0};

  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist.txt");
  Int_t nrun = 0;
  Int_t runnumber[1024];
  while(fin >> runnumber[nrun]) nrun++;
  fin.close();

  const Int_t nP = 50;
  const Int_t firstRun = Process * nP;
  const Int_t lastRun = (Process + 1) * nP;

  Float_t ALL[npT][1024] = {};
  Float_t eALL[npT][1024] = {};

  for(Int_t ir=firstRun; ir<lastRun; ir++)
  {
    TFile *f1 = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510ERT/8847/data/DirectPhotonPP_PhotonNode-%d.root", runnumber[ir]));
    TFile *f2 = new TFile(Form("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/PhotonNode-histo-%d.root", runnumber[ir]));
    TTree *T1 = (TTree*)f1->Get("T1");

    SpinPattern *spinpat = new SpinPattern;
    T1->SetBranchAddress("RUN/SpinPattern", &spinpat);
    T1->GetEntry(0);

    Long64_t bbc_same = 0;
    Long64_t bbc_opposite = 0;
    for(Int_t ib=0; ib<120; ib++)
    {
      Int_t pattern = spinpat->get_spinpattern_blue(ib) * spinpat->get_spinpattern_yellow(ib);
      if(pattern == 1)
        bbc_same += spinpat->get_bbc_wide(ib);
      else if(pattern == -1)
        bbc_opposite += spinpat->get_bbc_wide(ib);
    }

    Float_t pb = spinpat->get_pb();
    Float_t py = spinpat->get_py();

    if( bbc_same <= 0 || bbc_opposite <= 0 ||
        pb <= 0. || py <= 0. )
      continue;

    Long64_t yield_same_west[npT] = {};
    Long64_t yield_same_east[npT] = {};
    Long64_t yield_opposite_west[npT] = {};
    Long64_t yield_opposite_east[npT] = {};

    GetRawAsymm(f2, 0, 1, yield_same_west);
    GetRawAsymm(f2, 1, 1, yield_same_east);
    GetRawAsymm(f2, 0, -1, yield_opposite_west);
    GetRawAsymm(f2, 1, -1, yield_opposite_east);

    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      Float_t yield_same = ( yield_same_west[ipt] + yield_same_east[ipt] ) / 2.;
      Float_t yield_opposite = ( yield_opposite_west[ipt] + yield_opposite_east[ipt] ) / 2.;
      if( yield_same <= 0. || yield_opposite <= 0. ) continue;
      Float_t R = (Float_t)bbc_same / bbc_opposite;
      Float_t eR = R * sqrt( 1./yield_same + 1./yield_opposite );
      ALL[ipt][ir] = ( yield_same - R * yield_opposite ) / ( yield_same + R * yield_opposite ) / ( pb * py );
      eALL[ipt][ir] = sqrt( pow( 2. * R * yield_same * yield_opposite / (pb * py) / pow(yield_same + R * yield_opposite, 2.) , 2. )
          * pow( 1./yield_same + 1./yield_opposite + (eR*eR)/(R*R), 2.)
          + (1./pb + 1./py) * ALL[ipt][ir] * ALL[ipt][ir] );
    }

    delete spinpat;
    delete f1;
    delete f2;
  }

  for(Int_t ipt=0; ipt<npT; ipt++)
  {
    ofstream fout(Form("ALL/Process%d-pT%0d.txt", Process, ipt), ofstream::trunc);
    for(Int_t ir=firstRun; ir<lastRun; ir++)
      if( eALL[ipt][ir] > 0. )
        fout << runnumber[ir] << "\t" << ALL[ipt][ir]
          << "\t" << eALL[ipt][ir] << endl;
    fout.close();
  }
}
