#include <iostream>
#include <cmath>
#include <map>
#include <cstring>
#include <vector>
#include <TStyle.h>
#include <Pythia8/Pythia.h>
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>

#ifdef _CPPPWHGHOOKS
#include <Pythia8Plugins/PowhegHooks.h>
#define _CPPHOOKS PowhegHooks
#else
#include <Pythia8Plugins/QEDQCDPowhegHooks.h>
#define _CPPHOOKS QEDQCDPowhegHooks
#endif

using namespace Pythia8;

int main(int, char **);
int main(int argc, char **argv) {

  double CorrectPhiDelta(double a, double b);
  void Fill_For_Each_Weight(vector<TH1D> &vec_h, double val, vector<double> &vec_weights);

  int nFiles;
  string fileName;
  Pythia pythia;
  _CPPHOOKS *powhegHooks = 0; // POWHEG UserHooks
  bool loadhooks;
  pythia.readFile("anaLHEF.cmnd");

  //---read commandline args----------------------------------------
  if (argc < 3) {
    cout << endl << "Usage: " << argv[0]
      << " outputfile.root eventfile1.lhe eventfile2.lhe ..." << endl;
    exit(EXIT_FAILURE);
  } 
  const char *rootFileName = argv[1]; // output file
  nFiles = argc - 2; // number of event files to process

  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);

  // fill histograms with measured data points to compare with
  const int nPtBins = 30;
  const double etaAbsMin[3] = {0.0, 0.0, 0.0};
  const double etaAbsMax[3] = {0.25, 0.5, 1.0};
  const double ptBins[nPtBins+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

  TH1D *h_photon = new TH1D("h_photon", "Cross section of direct photon;p_{T} [GeV];#frac{d#sigma}{dp_{T}} [pb]", nPtBins, ptBins);
  vector< vector<TH1D> > vec_sim[2][2]; // two vector dimensions for different rap. bins & for different weights (e.g. for scale/pdf variation)


  // prepare bookkeeping of weights
  //----------------------------------------------------------------------
  const string sudaWeightID = "sudakovwgt";
  bool isSudaWeight = false; // was photon radiation enhanced?
  double sudaWeight = 1.;    // reweighting factor associated with radiation enhancement
  vector<double> vec_weights;   // shall later contain: sudaWeight * primary event weight (using vector to store multiple weights, e.g for scale/pdf variation)
  vector<string> vec_weightsID;// vector storing descriptive id of weights

  // pythia settings required for usage with powheg
  //----------------------------------------------------------------------
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Next:numberCount = 50000");

  // read in from conf file
  int vetoMode    = pythia.settings.mode("POWHEG:veto");
  int MPIvetoMode = pythia.settings.mode("POWHEG:MPIveto");
  loadhooks = (vetoMode > 0 || MPIvetoMode > 0);

  if (loadhooks) { // if NOT use SCALUP as starting scale
    if (vetoMode > 0) { // use kinematical limit as starting scale and veto
      pythia.readString("SpaceShower:pTmaxMatch = 2");
      pythia.readString("TimeShower:pTmaxMatch = 2");
    }
    if (MPIvetoMode > 0) {
      pythia.readString("MultipartonInteractions:pTmaxMatch = 2");
    }
    // activate POWHEG compliance
    powhegHooks = new _CPPHOOKS();
    pythia.setUserHooksPtr((UserHooks *) powhegHooks);
  }


  // variables to keep track of
  //----------------------------------------------------------------------
  int nEvents = 0,
      iPhoton = -1;     // index of photon in pythia event

  double ptMax = 0., // pT of hardest photon
         ptTemp = 0.,
         etaAbsPhoton = 999.;

  const double isoConeRadius = 0.5;
  double isoCone_mom;
  double isoCone_dR;


  // loop over lhef files showering each event
  //----------------------------------------------------------------------
  for (int iFile = 0; iFile < nFiles; iFile++) {

    fileName = argv[2 + iFile];
    cout << "Showering events in " << fileName << endl;

    // tell Pythia to use several lhe files while initializing once
    if (iFile == 1) pythia.readString("Beams:newLHEFsameInit = on");
    pythia.readString("Beams:LHEF = " + fileName);
    pythia.init();

    // skip pythia errors and break, when showering has reached the end of the LHE file
    //----------------------------------------------------------------------
    while (true) {
      if (!pythia.next()) {
        if (pythia.info.atEndOfFile()) break;
        continue;
      }

      nEvents++;

      // at very first event read in weight IDs and book histograms for each weight
      //----------------------------------------------------------------------
      if (nEvents == 1 && iFile == 0) {

        // check if the sudakov weight from enhanced radiation is present
        for (map<string,double>::iterator it = pythia.info.weights_detailed->begin();
            it != pythia.info.weights_detailed->end(); ++it) {
          if (it->first == sudaWeightID.c_str()) {
            isSudaWeight = true;
            printf("Sudakov reweighting of hard process is taken into account.\n");
            continue;
          }
        }

        // // if more weights at the same time are used,
        // // e.g. for scale or pdf variation, you can  access them like this
        // for (map<string,double>::iterator it = pythia.info.weights_detailed->begin();
        //     it != pythia.info.weights_detailed->end(); ++it) {
        //   if (it->first.find("scales") != std::string::npos){
        //     vec_weightsID.push_back(it->first);
        //   }
        // }

        // insert central value always at first position for convenience 
        // for (map<string,double>::iterator it = pythia.info.weights_detailed->begin();
        //     it != pythia.info.weights_detailed->end(); ++it) {
        //   if (it->first == "central"){ // NB: these strings follow 'lhrwgt_id' in powheg-input.save
        //     vec_weightsID.insert(vec_weightsID.begin(), it->first);
        //   }
        // }
        vec_weightsID.insert(vec_weightsID.begin(), "central");

        printf("Number of weights = %lu\n", vec_weightsID.size());
        for(long unsigned int i = 0; i < vec_weightsID.size(); i++)
          printf("weight description at position %lu: %s\n", i, vec_weightsID.at(i).c_str());

        // book histograms for each weight
        for(int iH=0; iH<2; iH++)
          for(int iso=0; iso<2; iso++)
            for(int i = 0; i < 3; i++){ // consider each rapidity bin
              vector <TH1D> vec_histo_temp;
              for(long unsigned int j = 0; j < vec_weightsID.size(); j++){
                vec_histo_temp.push_back( *(TH1D*)h_photon->Clone(Form("hard%d_iso%d_rap%d_%s",iH,iso,i,vec_weightsID.at(j).c_str())) );
              }
              vec_sim[iH][iso].push_back(vec_histo_temp);
            }	 
      } // back to the general event loop...


      // if Sudakov reweighting is activated, get corresponding weight for this event
      if (isSudaWeight) sudaWeight = pythia.info.getWeightsDetailedValue(sudaWeightID.c_str());

      // reload vector with regular weights * sudaWeight for this event 
      if(vec_weights.size() != 0) vec_weights.clear();
      //for(long unsigned int i = 0; i < vec_weightsID.size(); i++){
      //  vec_weights.push_back(pythia.info.getWeightsDetailedValue(vec_weightsID.at(i)) * sudaWeight);
      //}
      vec_weights.push_back(pythia.info.weight() * sudaWeight);

      // The actual event analysis starts here.
      ptMax  = 0.;
      ptTemp = 0.;
      iPhoton = -1;
      etaAbsPhoton = 999.;

      // search for hardest photon in this event
      //----------------------------------------------------------------------
      for (int i = 5; i < pythia.event.size(); i++) {
        if (pythia.event[i].id() == 22 && pythia.event[i].isFinal() && // final photon
            pythia.event[i].status() < 90 &&                      // no decay photons allowed, only direct photons
            TMath::Abs(pythia.event[i].eta()) < etaAbsMax[2] &&   // in maximal acceptance
            pythia.event[i].pT() > ptBins[4]){                    // in the pt reach of interest

          // find ptMax
          ptTemp = pythia.event[i].pT();
          if (ptTemp > ptMax) {
            ptMax = ptTemp;
            iPhoton = i; // remember index of hardest photon
          }
        }
      }

      // skip to next event, if no photon was found
      if(iPhoton < 0) continue;

      etaAbsPhoton = TMath::Abs(pythia.event[iPhoton].eta());

      // use following line to ignore events with extreme weights that can cause ugly fluctuations
      // but make sure the cross section does not decrease significantly
      if(ptMax > pythia.info.getScalesAttribute("uborns")*2.5){
        nEvents--;
        continue;
      }

      // loop over all direct photons in this event
      //----------------------------------------------------------------------
      for (int iDir = 5; iDir < pythia.event.size(); iDir++)
        if (pythia.event[iDir].id() == 22 && pythia.event[iDir].isFinal() && // final photon
            pythia.event[iDir].status() < 90 &&                      // no decay photons allowed, only direct photons
            TMath::Abs(pythia.event[iDir].eta()) < etaAbsMax[2] &&   // in maximal acceptance
            pythia.event[iDir].pT() > ptBins[4]){                    // in the pt reach of interest

          // check whether it is the hardest photon and get its pt and eta
          int iH = (iDir == iPhoton ? 1 : 0);
          double ptDir = pythia.event[iDir].pT();
          double eDir = pythia.event[iDir].e();
          double etaAbsDir = TMath::Abs(pythia.event[iDir].eta());

          // isolation cut: sum energy around photon and check whether threshold is reached
          //----------------------------------------------------------------------
          isoCone_mom = 0.; // reset sum of energy in cone
          for (int i = 5; i < pythia.event.size(); i++) {
            if ( !pythia.event[i].isFinal() ) continue;
            if ( !pythia.event[i].isVisible() ) continue;
            if ( TMath::Abs(pythia.event[i].eta()) > etaAbsMax[2]+isoConeRadius ) continue;
            if ( i == iDir ) continue;

            // distance between photon and particle at index i
            isoCone_dR = sqrt( pow2(CorrectPhiDelta(pythia.event[i].phi(), pythia.event[iDir].phi()))
                + pow2(pythia.event[i].eta() - pythia.event[iDir].eta()) );

            // sum energy in isolation cone
            if(isoCone_dR < isoConeRadius) isoCone_mom += pythia.event[i].pAbs();
          }

          // check whether threshold is reached
          int iso = (isoCone_mom < 0.1*eDir ? 1 : 0);

          // Fill histograms
          //----------------------------------------------------------------------
          for( int i = 0; i < 3; i++)
            if( etaAbsMin[i] < etaAbsDir &&
                etaAbsDir < etaAbsMax[i] )
              Fill_For_Each_Weight(vec_sim[iH][iso].at(i), ptDir, vec_weights);

        } // end of direct photon loop
    } // end of while loop; break if next file
  } // end of file loop

  // statistics on event generation
  pythia.stat();

  // write histograms to file ----------------------------------------
  TFile outFile(rootFileName, "RECREATE");

  // normalize simulated spectra for nEvents and pt bin width, then write
  for(int i = 0; i < 3; i++)
    for(unsigned long int j = 0; j < vec_weights.size(); j++){

      for(int iH=0; iH<2; iH++)
        vec_sim[iH][0].at(i).at(j).Add(&vec_sim[iH][1].at(i).at(j));
      for(int iso=0; iso<2; iso++)
        vec_sim[0][iso].at(i).at(j).Add(&vec_sim[1][iso].at(i).at(j));

      for(int iH=0; iH<2; iH++)
        for(int iso=0; iso<2; iso++)
        {
          vec_sim[iH][iso].at(i).at(j).Scale( 1./nEvents, "width");
          vec_sim[iH][iso].at(i).at(j).Write();
        }
    }

  outFile.Close();

  if (powhegHooks) delete powhegHooks;
  return 0;

}

// PYTHIA8's phi goes from -pi to pi; compute correct angle difference
//----------------------------------------------------------------------
double CorrectPhiDelta(double angle1, double angle2){
  double pi = TMath::Pi();
  double phi = TMath::Abs(angle1 - angle2);
  if(phi >= pi) return 2*pi-phi;
  else return phi;
}

//----------------------------------------------------------------------
void Fill_For_Each_Weight(vector<TH1D> &vec_h, double val, vector<double> &vec_weights){

  for(unsigned long int i = 0; i < vec_weights.size(); i++)
    vec_h.at(i).Fill(val,vec_weights.at(i));

  return;
}

//----------------------------------------------------------------------
