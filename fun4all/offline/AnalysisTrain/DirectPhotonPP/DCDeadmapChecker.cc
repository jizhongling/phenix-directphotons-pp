#include "DCDeadmapChecker.h"

#include <PHCentralTrack.h>
#include <emcClusterContent.h>

#include <TMath.h>
#include <TVector3.h>

#include <boost/foreach.hpp>

#define INSERT_ARRAY(vec,array) \
  do { \
    vec.insert( vec.end(), array, array+sizeof(array)/sizeof(array[0]) ); \
  } while(0)

using namespace std;

DCDeadmapChecker::DCDeadmapChecker(int eventsmod):
  imap(-1),
  nevents(0),
  nmod(eventsmod)
{
  const double cum_rRuns[nmap] = {0.06429, 0.1486, 0.1767, 0.214, 0.3447, 0.3501, 0.4817, 0.7094, 0.7141, 0.8178, 0.8253, 0.8963, 0.9003, 0.9978, 1};

  m_kbb.clear();
  for(int i=0; i<nmap; i++)
  {
    cum_nRuns[i] = TMath::Nint( cum_rRuns[i] * nmod );
    v_deadmap[i].clear();
  }

  m_kbb.insert( make_pair( "W0_0",   KBB( 7.075, -99.99, -1.113, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "W0_1",   KBB( 7.075,  79.67,  99.99, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "E0_0",   KBB(-7.049, -99.99, -0.699, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "E0_1",   KBB(-7.049,  80.19,  99.99, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "W1",     KBB( 0.000,  37.94,  43.18, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "W2",     KBB( 1.815,  25.18,  26.18, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "NE1",    KBB( 0.000,  67.42,  71.53, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "NE2",    KBB( 0.000,  29.72,  31.90, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "W3",     KBB(-2.363,  38.23,  39.23, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "SW1",    KBB( 6.400,  19.40,  29.09, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "W4",     KBB(-1.803,  62.98,  67.41, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "E1",     KBB(-2.348,  8.515,  11.68, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "E2",     KBB(-2.350,  38.62,  40.17, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "SE1",    KBB(-2.213,  54.55,  55.39, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "NE3",    KBB(-6.527,  19.74,  29.52, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "W5",     KBB( 1.386,  47.38,  51.05, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "NWh1_0", KBB( 27.52,  8.486,  13.15,  0.00,  0.15) ) );
  m_kbb.insert( make_pair( "NWh1_1", KBB( 27.52,  15.45,  19.53,  0.00,  0.15) ) );
  m_kbb.insert( make_pair( "NWh1_2", KBB( 27.52,  76.57,  79.36, -0.20,  0.00) ) );
  m_kbb.insert( make_pair( "SWh1_0", KBB( 5.774,  68.80,  72.06, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "SWh1_1", KBB( 24.33,  44.09,  46.89, -0.10,  0.10) ) );
  m_kbb.insert( make_pair( "NEh1_0", KBB(-33.23,  59.77,  67.44, -0.20,  0.20) ) );
  m_kbb.insert( make_pair( "NEh1_1", KBB(-32.56,  55.19,  58.07, -0.30,  0.30) ) );
  m_kbb.insert( make_pair( "SEh1_0", KBB(-40.68,  75.27,  80.18,  0.00,  0.20) ) );
  m_kbb.insert( make_pair( "SEh1_1", KBB(-29.85,  51.26,  54.74, -0.40,  0.40) ) );
  m_kbb.insert( make_pair( "Wh1",    KBB( 2.097,  49.84,  56.80, -10.0,  10.0) ) );
  m_kbb.insert( make_pair( "SEh2",   KBB( 0.000,  59.38,  60.69, -10.0,  10.0) ) );

  string common_edge[] = {"W0_0", "W0_1", "E0_0", "E0_1"};
  string common_hot[] = {"NWh1_0", "NWh1_1", "NWh1_2", "SWh1_0", "SWh1_1", "NEh1_0", "NEh1_1", "SEh1_0", "SEh1_1"};
  string w12[] = {"W1", "W2"};
  string w23[] = {"W2", "W3"};
  string w234[] = {"W2", "W3", "W4"};
  string e12se1[] = {"E1", "E2", "SE1"};

  for(int i=0; i<nmap; i++)
  {
    INSERT_ARRAY(v_deadmap[i], common_hot);
    INSERT_ARRAY(v_deadmap[i], common_edge);
  }
  for(int i=0; i<4; i++)
    INSERT_ARRAY(v_deadmap[i], w12);
  for(int i=4; i<7; i++)
    INSERT_ARRAY(v_deadmap[i], w23);
  for(int i=7; i<nmap; i++)
    INSERT_ARRAY(v_deadmap[i], w234);
  v_deadmap[0].push_back("NE1");
  v_deadmap[2].push_back("NE2");
  v_deadmap[5].push_back("SW1");
  INSERT_ARRAY(v_deadmap[8], e12se1);
  v_deadmap[10].push_back("NE3");
  v_deadmap[11].push_back("W5");
  v_deadmap[12].push_back("W5");
  v_deadmap[12].push_back("Wh1");
  v_deadmap[14].push_back("SEh2");
}

void DCDeadmapChecker::SetMapByRunnumber(int runnumber)
{
  const int run1[nmap] = {387027, 388261, 389558, 389588, 389904, 391442, 391465, 393066, 396067, 396268, 396889, 397049, 397531, 397577, 397737};
  const int run2[nmap] = {388052, 389557, 389587, 389768, 391377, 391450, 393064, 396054, 396075, 397000, 396910, 397534, 397534, 398149, 397738};

  imap = -1;
  for(int i=0; i<nmap; i++)
    if(runnumber >= run1[i] && runnumber <= run2[i])
      imap = i;

  return;
}

void DCDeadmapChecker::SetMapByEvent()
{
  if(imap < 0 || imap >= nmap || nevents >= nmod)
  {
    imap = 0;
    nevents = 0;
  }

  if(nevents >= cum_nRuns[imap])
    imap++;
  nevents++;

  return;
}

bool DCDeadmapChecker::IsDead(string nswe, double board, double alpha)
{
  if(imap < 0 || imap >= nmap)
    return false;

  BOOST_FOREACH(const string &smap, v_deadmap[imap])
  {
    if( smap[0] == nswe[1] ||
        (smap[0] == nswe[0] && smap[1] == nswe[1]) )
    {
      KBB kbb = m_kbb[smap];
      double b1 = kbb.k * alpha + kbb.b1;
      double b2 = kbb.k * alpha + kbb.b2;
      if( board < 0. || board > 80. ||
          (board > b1 && board < b2 &&
           alpha > kbb.alpha1 && alpha < kbb.alpha2) )
        return true;
    }
  }

  return false;
}

bool DCDeadmapChecker::IsDead(const PHCentralTrack *tracks, int itrk)
{
  /* phi and zed distributions */
  double phi = tracks->get_phi(itrk);
  double zed = tracks->get_zed(itrk);
  double alpha = tracks->get_alpha(itrk);
  int arm = tracks->get_dcarm(itrk);
  string ns = zed > 0. ? "N" : "S";
  string we = arm == 1 ? "W" : "E";
  string nswe = ns + we;
  double board = 0.;
  if( we == "W" )
    board = ( 0.573231 + phi - 0.0046 * cos( phi + 0.05721 ) ) / 0.01963496;
  else
    board = ( 3.72402 - phi + 0.008047 * cos( phi + 0.87851 ) ) / 0.01963496;

  return IsDead(nswe, board, alpha);
}

bool DCDeadmapChecker::ChargeVeto(const emcClusterContent *cluster, const PHCentralTrack *tracks)
{
  /* 3 sigma charge veto */
  int itrk_match = GetEmcMatchTrack(cluster, tracks);
  if( itrk_match >= 0 )
    return true;
  else 
    return false;
}

int DCDeadmapChecker::GetEmcMatchTrack(const emcClusterContent *cluster, const PHCentralTrack *tracks)
{
  int itrk_match = -1;
  double dzmin = 9999.;
  double mommax = 0.;

  TVector3 v3_cluster(cluster->x(), cluster->y(), cluster->z());

  int npart = tracks->get_npart();
  for(int itrk=0; itrk<npart; itrk++)
  {
    TVector3 v3_track(tracks->get_pemcx(itrk), tracks->get_pemcy(itrk), tracks->get_pemcz(itrk));
    double dphi = fabs((v3_track-v3_cluster).Phi());
    double dz = fabs((v3_track-v3_cluster).Z());
    double mom = tracks->get_mom(itrk);
    if( dphi > 0.015 ||
        !TMath::Finite(mom) ||
        IsDead(tracks, itrk) )
      continue;

    if( itrk_match != -1 )
    {
      if( dz < 8. && dz < dzmin )
      {
        itrk_match = itrk;
        dzmin = dz;
        mommax = mom;
      }
      else if( dzmin >= 8. && mom > mommax )
      {
        itrk_match = itrk;
        dzmin = dz;
        mommax = mom;
      }
    }
    else
    {
      itrk_match = itrk;
      dzmin = dz;
      mommax = mom;
    }
  }

  return itrk_match;
}
