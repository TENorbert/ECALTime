#include "AnaInput.h"

AnaInput::AnaInput( string datacardInput ) {

  datacardfile = datacardInput ;

}

AnaInput::~AnaInput(){

   cout<<" close input "<<endl ;

}

vector<TTree*> forestData ;

// GetForest , run first
void AnaInput::GetForest( string DataSet, TString treeName ) {

    cout<<"  =>>> getting a forest of "<< treeName <<endl ;
    vector<string> fileList;
    GetParameters( DataSet , &fileList );

    vector<TTree*> forest ;
    for ( size_t i =0 ; i< fileList.size(); i++ ) {
        TTree* tr = GetTree( fileList[i], treeName ) ;
        forestData.push_back( tr );
    }
}

TTree* AnaInput::TreeMap( string fileName ){

    vector<string> f0Names ;
    GetParameters( "TheData", &f0Names );

    TTree* theTr = 0;
    for ( size_t i=0; i< f0Names.size(); i++ ) {
        if ( f0Names[i] == fileName ) theTr = forestData[i] ;
    }

    return theTr ;
}

// get the tree from a specific file(s)
TTree* AnaInput::GetTree( string fName, TString treeName, TFile* file  ) {
  
  TTree* tr = 0;

  string filePath ;
  GetParameters( "RootFiles", &filePath );

  TString theFileName ;
  TChain* theChain = new TChain( treeName ) ;

  if ( fName[ fName.size()-1 ] == '+'  ) {
     string ChainName = fName.substr( 0, fName.size()-1 ) + "Chain"  ;
     vector<string> chainlist;
     GetParameters( ChainName, &chainlist );
     cout<<" * fileName+ = "<< ChainName <<endl;
     for ( size_t j=0; j< chainlist.size(); j++) {
         theFileName = filePath + chainlist[j]+".root" ;
         //cout<<" ** fileName = "<< theFileName <<endl;
         theChain->Add( theFileName );
     }
  } else {
    theFileName = filePath + fName+".root" ;
    cout<<" * fileName = "<< theFileName <<endl;
    theChain->Add( theFileName );
    //if ( file == NULL ) file = TFile::Open( theFileName );
    //tr = (TTree*) file->Get( treeName );

  }
  tr = theChain ;

  return tr ;
}


double AnaInput::NormalizeComponents(  string theChannel, string cfgFile ){

  if ( cfgFile == "" ) cfgFile = datacardfile ;

  double lumi ;
  double Scal = 1 ;
  GetParameters("Lumi", &lumi, cfgFile );

  vector<double> nEvents ;
  GetParameters( "nEvents" , &nEvents, cfgFile );
  vector<double> xsec;
  GetParameters("xsec", &xsec, cfgFile );
  vector<double> Eff;
  GetParameters("EffHLT", &Eff, cfgFile )  ;
  vector<string>  channel;
  GetParameters( "channel" , &channel, cfgFile );

  int idx = -1;
  for (size_t i =0; i< channel.size(); i++) {
      if ( channel[i] == theChannel )  idx = i ;
  }

  if ( idx >= 0 ) {
     double nBase = xsec[idx]*Eff[idx];
     Scal = (nBase*lumi) / nEvents[idx] ;
     cout<<" Scal of "<< channel[idx]<< " = " << Scal <<endl;
  
  } else {
     cout <<" No matched componenet !! " <<endl;
  }

  return Scal;
}

// Methods to read DataCard.txt
void AnaInput::GetParameters(string paraName, int* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;

     bool gotIt = false ;
     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 2;
           if ( pos < line.npos ) {
              string str_end = line.substr(vpos-1, 1) ;
              if ( str_end == ' ' || str_end == '=') {
                 getName  = line.substr( pos, paraName.size() );
                 getValue = line.substr( vpos );
                 *thePara = atoi( getValue.c_str() );
                 //cout<< paraName <<" = "<< *thePara << endl;
                 gotIt = true;
              }
           }
           if ( gotIt ) break ;
     }
     paraFile.close();
}

void AnaInput::GetParameters(string paraName, double* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;

     bool gotIt = false ;
     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 2;
           if ( pos < line.npos ) {
              string str_end = line.substr(vpos-1, 1) ;
              if ( str_end == ' ' || str_end == '=') {
                 getName  = line.substr( pos, paraName.size() );
                 getValue = line.substr( vpos );
                 *thePara = atof( getValue.c_str() );
                 //cout<< paraName <<" = "<< *thePara << endl;
                 gotIt = true ;
              }
           }
           if ( gotIt ) break ;
     }
     paraFile.close();
}

void AnaInput::GetParameters(string paraName, string* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     size_t  pos ;
     size_t  vpos ;

     bool gotIt = false ;
     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 2;

           if ( pos < line.npos ) {
              string str_end = line.substr(vpos-1, 1) ;
              if ( str_end == ' ' || str_end == '=') {
                 //cout<<" pos = "<< pos <<endl;
                 getName  = line.substr( pos, paraName.size() );
                 //*thePara = line.substr( vpos );
                 //cout<< paraName <<" = "<< *thePara << endl;
                 string strTmp = line.substr( vpos );
                 for (string::iterator it = strTmp.begin(); it< strTmp.end(); it++) {
                     if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') thePara->push_back( *it );
                 }
                 gotIt = true ;
              }
           }
           if ( gotIt ) break;
     }
     paraFile.close();
}

void AnaInput::GetParameters(string paraName, vector<double>* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;
     vector<double>  vvec;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 1;
           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              string arrVal = line.substr( vpos );
              if ( arrVal[0] != '=' && arrVal[0] != ' ' ) continue;
	      int vidx = 0;
	      string vtemp ;
	      //cout<< paraName <<" = ( " ;
              for (string::iterator it = arrVal.begin(); it< arrVal.end(); it++) {
                  if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') vtemp.push_back( *it );
                  if ( (*it) == ',' || (*it) == ')' ) { 
                     if ( vtemp.size() > 0 ) vvec.push_back( atof( vtemp.c_str() ) ) ;
		     //cout<< vtemp << *it;
		     vidx++ ;
		     vtemp.clear() ;
                  }
              }
              *thePara = vvec ;
           }
     }
     paraFile.close();

} 

void AnaInput::GetParameters(string paraName, vector<string>* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );

     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;
     vector<string>  vvec;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() ;
           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              string arrVal = line.substr( vpos );
              if ( arrVal[0] != '=' && arrVal[0] != ' ' ) continue;
	      int vidx = 0;
	      string vtemp ;
	      //cout<< paraName <<" = ( " ;
              for (string::iterator it = arrVal.begin(); it< arrVal.end(); it++) {
                  if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') vtemp.push_back( *it );
                  if ( (*it) == ',' || (*it) == ')' ) { 
                     if ( vtemp.size() > 0 ) vvec.push_back( vtemp ) ;
		     //cout<< vtemp << *it;
		     vidx++ ;
		     vtemp.clear() ;
                  }
              }
              *thePara = vvec ;
           }
     }
     paraFile.close();

}
 
void AnaInput::GetParameters(string paraName, vector<int>* thePara, string cfgFile ){

     if ( cfgFile == "" ) cfgFile = datacardfile ;

     fstream paraFile( cfgFile.c_str() );
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;
     vector<int>  vvec;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() ;
           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              string arrVal = line.substr( vpos );
              if ( arrVal[0] != '=' && arrVal[0] != ' ' ) continue;
	      //int vidx = 0;
	      string vtemp ;
	      //cout<< paraName <<" = ( " ;
              for (string::iterator it = arrVal.begin(); it< arrVal.end(); it++) {
                  if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') vtemp.push_back( *it );
                  if ( (*it) == ',' || (*it) == ')' ) { 
                     if ( vtemp.size() > 0 ) vvec.push_back( atoi( vtemp.c_str() ) ) ;
		     //cout<< vtemp << *it;
		     //vidx++ ;
		     vtemp.clear() ;
                  }
              }
              *thePara = vvec ;
           }
     }
     paraFile.close();

}
 

