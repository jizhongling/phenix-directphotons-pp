#include "HistogramBooker.h"

#include "Fun4AllHistoManager.h"

/* ROOT header */
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>

using namespace std;

HistogramBooker::HistogramBooker()
{
}

/* ----------------------------------------------- */

HistogramBooker::~HistogramBooker()
{
}

/* ----------------------------------------------- */

Fun4AllHistoManager* HistogramBooker::GetHistoManager( std::string managername )
{

  /* Create new HistogramManager */
  Fun4AllHistoManager* hm = new Fun4AllHistoManager( managername );

  /* Parameters common among multiple histograms */
  const int n_pTbins = 31;
  const double pTbins[n_pTbins+1] = { 0.0,
                                      0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
                                      5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
                                      12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0,
                                      100.0};

  /* phi bins */
  const double PI = 3.1415927;
  const double phi_sec[8] = {
    -PI/8, 0, PI/8, 2*PI/8,
    PI-2*PI/8, PI-PI/8, PI, PI+PI/8
  };

  const int n_phibins = 163;
  double phi_twr[n_phibins+1];
  for(int is=0; is<6; is++)
    for(int it=0; it<19; it++)
      phi_twr[19*is+it] = phi_sec[is] + 0.02 * ( it - 9 );
  for(int is=6; is<8; is++)
    for(int it=0; it<25; it++)
      phi_twr[114+25*(is-6)+it] = phi_sec[is] + 0.016 * ( it - 12 );


  /* ---------------------------------------------------
   * BEGIN: Event count and trigger checks ====>
   * --------------------------------------------------- */
  {
    /*
     * Histogram to count number of events
     */
    TH1* h1_events = new TH1I("h1_events", "number of events;trigger;# events", 8,0.5,8.5);
    h1_events->GetXaxis()->SetBinLabel(1, "all");
    h1_events->GetXaxis()->SetBinLabel(2, "BBCNarrow");
    h1_events->GetXaxis()->SetBinLabel(3, "BBCNarrow_zcut");
    h1_events->GetXaxis()->SetBinLabel(4, "BBCNarrow_zcut_ERT4x4a");
    h1_events->GetXaxis()->SetBinLabel(5, "BBCNarrow_zcut_ERT4x4b");
    h1_events->GetXaxis()->SetBinLabel(6, "BBCNarrow_zcut_ERT4x4c");
    h1_events->GetXaxis()->SetBinLabel(7, "BBCNarrow_zcut_ERT4x4or");
    h1_events->GetXaxis()->SetBinLabel(8, "ERT4x4or");
    hm->registerHisto( h1_events , 1 );

    /*
     * 3D histogram to count for trigger efficiency
     */
    //    TH3 *h3_trig = new TH3F("h3_trig", "number of clusters;p_{T} [GeV];sector;trigger;", n_pTbins,0.,0., 8,-0.5,7.5, 4,0.5,4.5);
    //    h3_trig->GetXaxis()->Set(n_pTbins, pTbins);
    //    h3_trig->GetZaxis()->SetBinLabel(1, "all");
    //    h3_trig->GetZaxis()->SetBinLabel(2, "ERT4x4a");
    //    h3_trig->GetZaxis()->SetBinLabel(3, "ERT4x4b");
    //    h3_trig->GetZaxis()->SetBinLabel(4, "ERT4x4c");
    //    hm->registerHisto(h3_trig, 1);

    /*
     * Trigger efficiency for pion
     */
    //    TH3 *h3_trig_pion = static_cast<TH3*>( h3_trig->Clone("h3_trig_pion") );
    //    hm->registerHisto(h3_trig_pion, 1);
  }
  /* ---------------------------------------------------
   * <==== END: Event count and trigger checks
   * --------------------------------------------------- */


  /* ---------------------------------------------------
   * BEGIN: Warnmap checks ====>
   * --------------------------------------------------- */
  {
    /*
     * 2D histogram storing transverse momentum of each
     * cluster in each sector with warnmap applied
     */
    TH2* h2_pT_1cluster = new TH2F("h2_pT_1cluster", "Single cluster transverse momentum (warnmap applied);p_{T} [GeV];EMCal sector;", n_pTbins,pTbins, 8,-0.5,7.5);
    hm->registerHisto( h2_pT_1cluster , 1 );

    /*
     * 2D histogram storing transverse momentum of each
     * cluster in each sector without warnmap applied
     */
    TH2* h2_pT_1cluster_nowarn = new TH2F("h2_pT_1cluster_nowarn", "Single cluster transverse momentum (warnmap not applied);p_{T} [GeV];EMCal sector;", n_pTbins,pTbins, 8,-0.5,7.5);
    hm->registerHisto( h2_pT_1cluster_nowarn , 1 );
  }
  /* ---------------------------------------------------
   * <==== END: Warnmap checks
   * --------------------------------------------------- */


  /* ---------------------------------------------------
   * BEGIN: TOF crosschecks ====>
   * --------------------------------------------------- */
  {
    /*
     * 3D histogram of sector, pT and TOF to check TOF calibration
     */
    int ndim_tof = 3;
    int nbins_tof[] = { 8, n_pTbins, 1001 };
    double xmin_tof[] = { -0.5, 0, -100.5 };
    double xmax_tof[] = {  7.5, 0, 100.5 };
    THnSparse* hn_tof = new THnSparseF("hn_tof",
                                       "TOF;EMCal sector;p_{T} [GeV];TOF [ns];",
                                       ndim_tof,
                                       nbins_tof,
                                       xmin_tof,
                                       xmax_tof );

    hn_tof->SetBinEdges(1,pTbins);
    hn_tof->GetAxis(0)->SetName("EMCalSector");
    hn_tof->GetAxis(1)->SetName("pT");
    hn_tof->GetAxis(2)->SetName("tof");
    hm->registerHisto( hn_tof , 1 );

    /*
     * Using TOF from DST file without local recalibration.
     */
    THnSparse* hn_tof_raw = static_cast<THnSparse*>(hn_tof->Clone("hn_tof_raw"));
    hn_tof_raw->SetTitle("TOF (uncalibrated)");
    hm->registerHisto( hn_tof_raw , 1 );
  }
  /* ---------------------------------------------------
   * <==== END: TOF crosschecks
   * --------------------------------------------------- */


  /* ---------------------------------------------------
   * BEGIN: Pi0 crosschecks (calibration etc) ====>
   * --------------------------------------------------- */
  {
    /*
     * Histogram storing invariant mass of photon
     * candidate pairs in different sectors and pT bins.
     * Used to check sector-by-sector EMCal energy calibration.
     *
     * Storing number of identified direct photon candidates paired with other photon in event in bins of
     *
     * - sector
     * - transverse momentum of photon pair
     * - invariant mass (paired with any other photon in event) (default: 300, 0, 0.3)
     * - eta
     * - phi
     * - trigger
     *
     */
    int ndim_hn_pi0calib = 6;
    int nbins_hn_pi0calib[] = {8, n_pTbins, 1200, 70, n_phibins, 4};
    double xmin_hn_pi0calib[] = {-0.5, 0., 0., -0.35, 0., -0.5};
    double xmax_hn_pi0calib[] = {7.5, 0., 1.2, 0.35, 0., 3.5};

    THnSparse* hn_pi0calib = new THnSparseF("hn_pi0",
                                            "Photon pair invariant mass;EMCal sector;p_{T} [GeV];m_{inv} [GeV];#eta;#phi [rad];ERT trigger;",
                                            ndim_hn_pi0calib,
                                            nbins_hn_pi0calib,
                                            xmin_hn_pi0calib,
                                            xmax_hn_pi0calib );
    hn_pi0calib->SetBinEdges(1,pTbins);
    hn_pi0calib->GetAxis(0)->SetName("EMCalSector");
    hn_pi0calib->GetAxis(1)->SetName("pT");
    hn_pi0calib->GetAxis(2)->SetName("minv");
    hn_pi0calib->SetBinEdges(4, phi_twr);
    hm->registerHisto( hn_pi0calib , 1 );

    /*
     * Using energies from DST file without local recalibration.
     */
    THnSparse* hn_pi0calib_raw = static_cast<THnSparse*>(hn_pi0calib->Clone("hn_pi0_raw"));
    hn_pi0calib_raw->SetTitle("Photon pair invariant mass (uncalibrated)");
    hm->registerHisto( hn_pi0calib_raw , 1 );

    /*
     * Using energies from DST file without TOF cut.
     */
    THnSparse* hn_pi0calib_notof = static_cast<THnSparse*>(hn_pi0calib->Clone("hn_pi0_notof"));
    hn_pi0calib_notof->SetTitle("Photon pair invariant mass (calibrated, no TOF cut)");
    hm->registerHisto( hn_pi0calib_notof , 1 );

  }
  /* ---------------------------------------------------
   * <==== END: Pi0 crosschecks (calibration etc) ====>
   * --------------------------------------------------- */


  /* ---------------------------------------------------
   * BEGIN: Direct Photon Candidates ====>
   * --------------------------------------------------- */
  {
    /*
     * Histogram to count number of direct photon candidates per events
     */
    TH1* h1_nphotons = new TH1I("h1_nphotons","frequency;# photons / event", 11,-0.5,10.5);
    hm->registerHisto( h1_nphotons , 1 );

    /*
     * storing number of identified direct photon candidates in bins of
     *
     * - sector
     * - transverse momentum of photon
     * - photon energy
     * - eta
     * - phi
     * - trigger? MB, ERT-a, ERT-b, ERT-c
     * - isolated photon? 0 = no, 1 = yes
     */
    int ndim_hn_1photon = 7;
    int nbins_hn_1photon[] = {8, n_pTbins, n_pTbins, 70, n_phibins, 4, 2};
    double xmin_hn_1photon[] = {-0.5, 0., 0., -0.35, 0., -0.5, -0.5};
    double xmax_hn_1photon[] = {7.5, 0., 0., 0.35, 0., 3.5, 1.5};

    THnSparse* hn_1photon = new THnSparseF("hn_1photon",
                                           "Single Photon Spectrum;EMCal sector;p_{T} [GeV];E_{#gamma} [GeV];#eta;#phi [rad];ERT trigger;IsolationStatus;",
                                           ndim_hn_1photon,
                                           nbins_hn_1photon,
                                           xmin_hn_1photon,
                                           xmax_hn_1photon );

    hn_1photon->SetBinEdges(1,pTbins);
    hn_1photon->SetBinEdges(2,pTbins);
    hn_1photon->GetAxis(0)->SetName("EMCalSector");
    hn_1photon->GetAxis(1)->SetName("pT");
    hn_1photon->GetAxis(2)->SetName("Egamma");
    hn_1photon->SetBinEdges(4, phi_twr);
    hm->registerHisto( hn_1photon , 1 );

    /*
     * storing number of identified direct photon candidates paired with other photon in event in bins of
     *
     * - sector
     * - transverse momentum of photon
     * - invariant mass of photon pair
     * - eta
     * - phi
     * - trigger? MB, ERT-a, ERT-b, ERT-c
     * - isolated photon? 0 = no, 1 = yes
     *
     */
    int ndim_hn_2photon = 7;
    int nbins_hn_2photon[] = {8, n_pTbins, 300, 70, n_phibins, 4, 2};
    double xmin_hn_2photon[] = {-0.5, 0., 0., -0.35, 0., -0.5, -0.5};
    double xmax_hn_2photon[] = {7.5, 0., 0.3, 0.35, 0., 3.5, 1.5};

    THnSparse* hn_2photon = new THnSparseF("hn_2photon",
                                           "Photon Pair Spectrum;EMCal sector;p_{T} [GeV];m_{inv} [GeV];#eta;#phi [rad];ERT trigger;IsolationStatus;",
                                           ndim_hn_2photon,
                                           nbins_hn_2photon,
                                           xmin_hn_2photon,
                                           xmax_hn_2photon );

    hn_2photon->SetBinEdges(1,pTbins);
    hn_2photon->GetAxis(0)->SetName("EMCalSector");
    hn_2photon->GetAxis(1)->SetName("pT");
    hn_2photon->GetAxis(2)->SetName("minv");
    hn_2photon->SetBinEdges(4, phi_twr);
    hm->registerHisto( hn_2photon , 1 );

  }
  /* ---------------------------------------------------
   * <==== END: Direct Photon Candidates ====>
   * --------------------------------------------------- */

  return hm;
}
