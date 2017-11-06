#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdlib>
#include <ctime>
//#include <cmath>

#include <TPythia6.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include <PHIODataNode.h>
#include <PHObject.h>
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHTimeStamp.h>
#include <Fun4AllReturnCodes.h>
#include <getClass.h>

#include <PHPythiaHeader.h>
#include <PHPythiaContainer.h>
#include <AnaPHPythiaDirectPhoton.h>

using namespace std;

AnaPHPythiaDirectPhoton::AnaPHPythiaDirectPhoton(const std::string &name): SubsysReco(name)
{
  phpythia = 0;         // array of pythia particles
  phpythiaheader = 0;   // pythia header
}

AnaPHPythiaDirectPhoton::~AnaPHPythiaDirectPhoton()
{
}

// Done at the beginning of processing
int AnaPHPythiaDirectPhoton::InitRun(PHCompositeNode *topNode)
{
  return EVENT_OK;
}

// Done at the end of processing
int AnaPHPythiaDirectPhoton::End(PHCompositeNode *topNode)
{
  return EVENT_OK;
}

// Process each event
int AnaPHPythiaDirectPhoton::process_event(PHCompositeNode *topNode)
{
  // Get PYTHIA Header
  phpythiaheader = findNode::getClass<PHPythiaHeader>(topNode,"PHPythiaHeader");
  if (!phpythiaheader)
    {
      cout << PHWHERE << "Unable to get PHPythiaHeader, is Node missing?" << endl;
      return ABORTEVENT;
    }

  // Get PYTHIA Particles
  phpythia = findNode::getClass<PHPythiaContainer>(topNode,"PHPythia");
  if (!phpythia)
    {
      cout << PHWHERE << "Unable to get PHPythia, is Node missing?" << endl;
      return ABORTEVENT;
    }

  // Print some header information
  cout << "Event: " << phpythiaheader->GetEvt() << "\t" << phpythiaheader->GetNpart() << endl;

  cout << "KS\tKF\tpid\tname\thistory" << endl;
  // Print information from each
  int npart = phpythia->size();
  for (int ipart=0; ipart<npart; ipart++)
    {
      TMCParticle *part = phpythia->getParticle(ipart);
      //TMCParticle *parent = phpythia->getParent(part);
      cout << setw(4) << phpythia->getLineNumber(part)
           << "\t" << part->GetKS()
           << setw(8)  <<  part->GetKF() // << "\t" << part->GetEnergy()
           << setw(12) << part->GetName()
        //<< "\t" << ipart
        //<< "\t" << part->GetParent()
        //<< "\t" << part->GetFirstChild()
        //<< "\t" << part->GetLastChild()
        //<< "\t" << phpythia->getChildNumber(part)
           << "\t" << phpythia->getHistoryString(part)
           << endl;
    }

  cout << endl << "Example 1 : Get ancestor Ks0(pid=310)'s energy" << endl;
  for (int ipart=0; ipart<npart; ipart++)  {
    TMCParticle *part = phpythia->getParticle(ipart);
    if (part->GetKS() != 1) continue;

    TMCParticle *anc = phpythia->hasAncestor(part, 310);
    float energy = 0.0;
    if (anc) {
      energy = anc->GetEnergy();
      cout << setw(12) << phpythia->getLineNumber(part)
           << setw(12) << part->GetName()
           << " 's anc "
           << setw(12) << anc->GetName() << setw(12) << energy << endl;
    }
  }

  cout << endl << "Example 2 : Get stable particle pi-, which has ancestor K0(311) and Ks0(pid=310),  get K0's energy" << endl;
  for (int ipart=0; ipart<npart; ipart++)  {
    TMCParticle *part = phpythia->getParticle(ipart);
    if (part->GetKF() != -211) continue;

    if (phpythia->hasAncestor(part, 311, 310)) {
      TMCParticle *anc = phpythia->hasAncestor(part, 311);
      float energy = 0.0;
      if (anc) {
        energy = anc->GetEnergy();
        cout << setw(12) << phpythia->getLineNumber(part)
             << setw(12) << part->GetName()
             << " 's anc "
             << setw(12) << anc->GetName() << setw(12) << energy << endl;
      }
    }
  }
  cout << endl;

  // some useful functions
  // phpythia->getHistoryString(part)               // return a TString with part's all ancestor
  // phpythia->getParent(part)                      // return pointer to part's parent TMCParticle
  // phpythia->getChild(part, i)                    // return pointer to part's ith child TMCParticle
  // phpythia->getChildNumber(part)                 // return part's child number
  // phpythia->hasAncestor(part, pid1, pid2, pid3)  // whether part has all ancestor pid1, pid2, pid3
  // phpythia->hasAncestorAbs(part, pid)            // return the nearest ancestor pointer with abs(KF)==pid
  // phpythia->hasAncestorRange(part, pid1, pid2)   // return the nearest ancestor whose pid1<=KF<=pid2


  return EVENT_OK;
}

