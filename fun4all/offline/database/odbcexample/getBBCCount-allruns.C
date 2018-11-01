#include <sql.h>
#include <odbc++/setup.h>
#include <odbc++/types.h>
#include <odbc++/drivermanager.h>
#include <odbc++/errorhandler.h>
#include <odbc++/connection.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>
#include <odbc++/databasemetadata.h>
#include <odbc++/statement.h>
#include <odbc++/preparedstatement.h>

#include <phool.h>
#include <TFile.h>
#include <TTree.h>

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

using namespace std;
using namespace odbc;

int main()
{
  Connection *con = 0;
  Statement *stmt = 0;
  ResultSet *rs = 0;
  ostringstream cmd;

  // Connect to database daq as user phnxrc. Database server
  // name comes from /opt/phenix/etc/odbc.ini config file. The third
  // argument is the password and left blank because it's automatically
  // picked up from odbc.ini file 
  try
  {
    con = DriverManager::getConnection("daq", "phnxrc", "");
  }
  catch(SQLException &e)
  {
    cerr << PHWHERE << "Exception caught during DriverManager::getConnection" << endl;
    cerr << e.getMessage() << endl;
    return 1;
  }

  if(con)
  {
    cout << "connected" << endl;
  }

  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runnumber.txt");
  int nrun = 0;
  int runnumber[1024];
  while(fin >> runnumber[nrun]) nrun++;
  fin.close();

  // runnumber and clock counts to fill the TTree
  int entry_count[10] = {};
  ULong64_t runno;
  ULong64_t clock;
  ULong64_t bbcnovtx_count;
  ULong64_t bbcnarrow_count;
  ULong64_t bbcnovtx_scaledown;
  ULong64_t bbcnarrow_scaledown;
  ULong64_t erta_scaledown;
  ULong64_t ertb_scaledown;
  ULong64_t ertc_scaledown;
  ULong64_t ertc_sum = 0;

  TFile *fout = new TFile("clock-counts.root", "RECREATE");
  TTree *t1 = new TTree("t1", "Clock and BBC counts");
  t1->Branch("runnumber", &runno, "runnumber/l");
  t1->Branch("clock_live", &clock, "clock_live/l");
  t1->Branch("bbcnovtx_live", &bbcnovtx_count, "bbcnovtx_live/l");
  t1->Branch("bbcnarrow_live", &bbcnarrow_count, "bbcnarrow_live/l");
  t1->Branch("bbcnovtx_scaledown", &bbcnovtx_scaledown, "bbcnovtx_scaledown/l");
  t1->Branch("bbcnarrow_scaledown", &bbcnarrow_scaledown, "bbcnarrow_scaledown/l");
  t1->Branch("erta_scaledown", &erta_scaledown, "erta_scaledown/l");
  t1->Branch("ertb_scaledown", &ertb_scaledown, "ertb_scaledown/l");
  t1->Branch("ertc_scaledown", &ertc_scaledown, "ertc_scaledown/l");

  for(int i=0; i<nrun; i++)
  {
    // get scalerglive from daq table
    cmd.str("");
    cmd << "select * from trigger where runnumber = " << runnumber[i] << " order by bitnb";
    stmt = con->createStatement();
    try
    {
      rs = stmt->executeQuery(cmd.str().c_str());
    }
    catch(SQLException &e)
    {
      cerr << e.getMessage() << endl; 
      return 1;
    }

    // show attributes: runnumber, name, scalerglive
    //cout << rs->getMetaData()->getColumnName(1).c_str() << "   " ;
    //cout << rs->getInt(1) << endl;
    //cout << rs->getMetaData()->getColumnName(2).c_str() << "   " ;
    //cout << rs->getString(2) << endl;
    //cout << rs->getMetaData()->getColumnName(26).c_str() << "   " ;
    //cout << rs->getLong(26) << endl;

    ostringstream name;
    do
    {
      name.str("");
      rs->next();
      try
      {
        name << rs->getString(2);
      }
      catch(SQLException &e)
      {
        cerr << "Exception caught during ResultSet::getString: " << e.getMessage() << endl;
        break;
      }

      if( name.str().compare("BBCLL1(>0 tubes) novertex") == 0 )
      {
        entry_count[0]++;
        bbcnovtx_count = rs->getLong(26);
        bbcnovtx_scaledown = rs->getLong(5);
      }
      else if( name.str().compare("BBCLL1(>0 tubes) narrowvtx") == 0 )
      {
        entry_count[1]++;
        bbcnarrow_count = rs->getLong(26);
        bbcnarrow_scaledown = rs->getLong(5);
      }
      else if( name.str().compare(0, 8, "ERT_4x4b") == 0 )
      {
        entry_count[2]++;
        ertb_scaledown = rs->getLong(5);
      }
      else if( name.str().compare(0, 8, "ERT_4x4a") == 0 ||
          name.str().compare(0, 11, "ERTLL1_4x4a") == 0 )
      {
        entry_count[3]++;
        erta_scaledown = rs->getLong(5);
      }
      else if( name.str().compare(0, 8, "ERT_4x4c") == 0 )
      {
        entry_count[4]++;
        ertc_scaledown = rs->getLong(5);
        ertc_sum += rs->getLong(26);
      }
      else if( name.str().compare(0, 5, "CLOCK") == 0 )
      {
        entry_count[5]++;
        clock = rs->getLong(26);
        break;
      }
    } while( !name.str().empty() );

    runno = runnumber[i];
    t1->Fill();
  }

  t1->Write();
  for(int ien=0; ien<10; ien++)
    if( entry_count[ien] > 0 )
      cout << entry_count[ien] << " ";
  cout << endl;
  cout << "ERT4x4c Sum: " << ertc_sum << endl;

  delete t1;
  delete fout;
  delete rs;
  delete con;
  return 0;
}
