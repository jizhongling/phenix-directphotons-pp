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

#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;
using namespace odbc;

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    cerr << "Usage: " << argv[0] << " <runnumber>" << endl;
    return 1;
  }

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
    //return 1;
  }

  if(con)
  {
    cout << "connected" << endl;
  }

  // get scalerglive from daq table
  cmd.str("");
  cmd << "select * from trigger where runnumber = " << argv[1] << " order by bitnb";
  stmt = con->createStatement();
  try
  {
    rs = stmt->executeQuery(cmd.str().c_str());
  }
  catch(SQLException &e)
  {
    cerr << e.getMessage() << endl; 
    //return 1;
  }

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
    // show name and scalerglive attributes for BBCLL1
    if( name.str().compare(0, 6, "BBCLL1") == 0 )
    {
      cout << rs->getMetaData()->getColumnName(2).c_str() << "   "
        << name.str() << endl;
      cout << rs->getMetaData()->getColumnName(26).c_str() << "   "
        << rs->getLong(26) << endl;
      cout << endl;
    }
    // show all attributes for BBCLL1(>0 tubes) narrowvtx
    if( name.str().compare("BBCLL1(>0 tubes) narrowvtx") == 0 )
    {
      for(int ien=1; ien<=33; ien++)
        cout << "Entry " << ien << ": " << rs->getMetaData()->getColumnName(ien).c_str()
          << " = " << (ien==2?name.str():rs->getString(ien)) << endl;
    }
  } while( !name.str().empty() );

  delete rs;
  delete con;
  return 0;
}
