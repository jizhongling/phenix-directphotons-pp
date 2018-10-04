#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <odbc++/connection.h>
#include <odbc++/setup.h>
#include <odbc++/types.h>
#include <odbc++/errorhandler.h>
#include <sql.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>
#include <odbc++/preparedstatement.h>
#include <odbc++/databasemetadata.h>
#include "phool.h"

using namespace odbc;
using namespace std;

#define LEN 7

int main () {
  
  Connection* con = 0;
  Statement* stmt = 0;
  ResultSet* rs = 0;
  ostringstream cmd;
  //array will hold values of postgres arrays as strings, then we use
  //atof() and atoi() to get numeric values
  std::vector <std::string> array;
  size_t position;
  string ss;
  
  //Connect to database sandbox as user phnxrc. Database server
  //name comes from /opt/phenix/etc/odbc.ini config file. The third
  //argument is the password and left blank because it's automatically
  //picked up from odbc.ini file 
  try
    {
      con = DriverManager::getConnection("sandbox", "phnxrc", "");
    }
  catch (SQLException& e)
    {
      cout << PHWHERE
	   << " Exception caught during DriverManager::getConnection" << endl;
      cout << e.getMessage() << endl;
      return 1;
    }
  
  if(con){
    cout << "connected" << endl;
  }
  //  exit(0);
  //set command string to an empty string
  cmd.str("");
  
  //set simple (not array) attributes
  int runnum = 123;
  int fillnum = 10;
  int badrun = 0;
  int crossingshift = 1;
  float transvers = 0.5;
  float transverserr = 0.01;
  
  //declare array attributes
  float pb[120];
  float pbe[120];
  float py[120];
  float pye[120];
  int patternb[120];
  int patterny[120];
  int bbcvcut[120];
  int bbcwocut[120];
  int zdcnarrow[120];
  int zdcwide[120];
  int badbunchqa[120];
  
  //first insert simple attributes and only the first element for arrays
  cmd << "insert into spin values("  << runnum << "," << fillnum << "," << badrun << "," << crossingshift << ",'{0.1}','{0.1}','{0.1}','{0.1}','{1}','{1}','{1}','{1}','{1}','{1}','{1}'," << transvers <<  "," << transverserr << ")";
  stmt = con->createStatement();
  try{
    rs = stmt->executeQuery(cmd.str().c_str());
  }
  catch (SQLException& e)
    {
      cout << e.getMessage() << endl; 
      //  return 1;
    }
  
  
  //now fill array attributes, eventually LEN would be 120 
  for (int i = 1; i<LEN; i++){
    pb[i] = 0.5;
    pbe[i] = 0.1;
    py[i] = 1.23e-13;
    pye[i] = 0.01;
    patternb [i] = 1;
    patterny [i] = 2;
    bbcvcut[i] = 3; 
    bbcwocut[i] = 4; 
    zdcnarrow [i] = 5;
    zdcwide [i] = 6;
    badbunchqa[i] = 7;
    
    cmd.str("");
    cmd << "update spin set polarblue[" << i << "]=" << pb[i] <<
      ",polarblueerror[" << i << "]=" << pbe[i] <<  
      ",polaryellow[" << i << "]=" << py[i] <<  
      ",polaryellowerror[" << i << "]=" << pye[i] <<  
      ",spinpatternblue [" << i << "]=" << patternb[i] <<  
      ",spinpatternyellow[" << i << "]=" << patterny[i] <<  
      ",bbcvertexcut[" << i << "]=" << bbcvcut[i] <<  
      ",bbcwithoutcut[" << i << "]=" << bbcwocut[i] <<  
      ",zdcnarrow[" << i << "]=" << zdcnarrow[i] <<  
      ",zdcwide[" << i << "]=" << zdcwide[i] << 
      ",badbunchqa[" << i << "]=" << badbunchqa[i] << "where runnumber=123;" ;
    
    cout << cmd.str() << endl;
    stmt = con->createStatement();
    try{
      rs = stmt->executeQuery(cmd.str().c_str());
    }
    catch (SQLException& e)
      {
	cout << e.getMessage() << endl; 
	//  return 1;
      }
  }
  //*/
  
  
  /// see what's in the first row of spin table
  cmd.str("");
  cmd << "select * from spin where runnumber=123";
  stmt = con->createStatement();
  try{
    rs = stmt->executeQuery(cmd.str().c_str());
  }
  catch (SQLException& e)
    {
      cout << e.getMessage() << endl; 
      //  return 1;
    }
  rs->next();
  cout << endl;
  
  //show first 4 attributes: runnumber, fillnumber, badrunqa, crossingshift
  for ( unsigned i = 1; i < 5; i++){
    cout << ((rs->getMetaData())->getColumnName(i)).c_str() << "   " ;
    cout << rs->getInt(i) << endl;
  }
  //show array attributes
  for ( unsigned i = 5; i < 16; i++){
    cout << ((rs->getMetaData())->getColumnName(i)).c_str() << "    " ;
    ss = "";
    array.clear();
    position = 0;
    ss = rs->getString(i);
    //remove opening {
    ss = ss.substr(1,ss.size());
    //remove closing }
    ss = ss.substr(0,ss.size()-1);
    std::string _separator = ",";
    
    position = ss.find_first_of(_separator); // returns npos when fails
    while ( position != ss.npos )
      {
	
	// This thing here checks that we dont push empty strings
	// to the array
	if ( position != 0 )
	  array.push_back( ss.substr( 0, position ) );
	
	// When the cutted part is pushed into the array we
	// remove it and the separator from the ss
	ss.erase( 0, position + _separator.length() );
	
	// And the we look for a new place for the _separator
	position = ss.find_first_of( _separator );
      }
    // We will push the rest of the stuff in to the array
    if ( ss.empty() == false )
      array.push_back( ss );
    for ( vector<string>::iterator it = array.begin(); it != array.end(); ++it ){
      const char * cc = (*it).c_str();
      if(i < 9){ // attributes from 5 to 8 are arrays of reals
	cout << atof(cc) << "  ";
      }
      else { // attr from 9 to 15 are arrays of ints
	cout << atoi(cc) << "  ";
      }
    }
    cout << endl;
  }
  //finally print the  last 2 attributes: 
  for ( unsigned i = 16; i < 18; i++){
    cout << ((rs->getMetaData())->getColumnName(i)).c_str() << "   " ;
    cout << rs->getFloat(i) << endl;
  }

  // delete the table entry we just inserted
  cmd.str("");
  cmd << "delete from spin where runnumber = 123";
  stmt = con->createStatement();
  try{
    int n = stmt->executeUpdate(cmd.str().c_str());
    cout << "deleted " << n << " row" << endl;
  }
  catch (SQLException& e)
    {
      cout << e.getMessage() << endl; 
      //  return 1;
    }
  
  delete rs;
  delete con;
  return 0;
}
