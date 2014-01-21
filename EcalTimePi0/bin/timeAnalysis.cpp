#include <iostream>
#include <string>

#include <map>
#include <vector>
#include <algorithm> 
#include <functional>
#include <set>
#include <boost/tokenizer.hpp>

#include <iostream>
#include <math.h> 
#include <assert.h>

#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/PluginManager/interface/standard.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/LorentzVector.h"

//<<<<<<< timeAnalysis.cpp
#include "ECALTime/EcalTimePi0/interface/EcalTimePhyTreeContent.h"
// remove line below, when established
//#include "CalibCalorimetry/EcalTiming/interface/EcalTimeTreeContent.h"

//=======
//included via EcalObjectTime.h  ->  analysisHisograms.h
//#include "CalibCalorimetry/EcalTiming/interface/EcalTimeTreeContent.h"
//>>>>>>> 1.3
#include "ECALTime/EcalTimePi0/interface/timeVsAmpliCorrector.h"
	//<<<<<<< timeAnalysis.cpp
//#include "ECALTime/EcalTimePi0/interface/EcalObjectTime.h"
#include "ECALTime/EcalTimePi0/interface/DPSelection.h"
	//=======

#include "ECALTime/EcalTimePi0/interface/analysisHistograms.h"
// following is commented out because already included through analysisHisograms.h
//#include "ECALTime/EcalTimePi0/interface/EcalObjectTime.h"

//>>>>>>> 1.3

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"


typedef std::set<std::pair<int,int> > SetOfIntPairs;

// authors: S. Cooper and G. Franzoni (UMN)
// Second Author: Tambe E. Norbert (UMN)
#define BarrelLimit  1.479
#define EndcapLimit  3.0

#define ADCtoGeVEB   0.039
#define ADCtoGeVEE   0.063

#define numAeffBins     35
#define numAoSigmaBins  25

#define NSlices 54   // # of slices in log(A)
#define LogStep 0.05 // size of Slices log(A)

#define lightSpeed 299792458



// -------- Globals ----------------------------------------
EcalTimePhyTreeContent treeVars_; 
TFile* saving_;
std::vector<std::string> listOfFiles_;
bool speak_=false;
char buffer_ [75];
std::string bufferTitle_; 
AnaInput    * Input_;    // imports values from a .txt datacard
DPSelection * Selector_; // 
// default settings
std::string outputRootName_ = "outputHistos.root";
int   numEvents_      = -1;
unsigned int  minRun_ = 0;
unsigned int  maxRun_ = 9999999;
unsigned int  minLS_ = 0;
unsigned int  maxLS_ = 9999999;
float eTGammaMinEB_   = 0.2;
float s4s9GammaMinEB_ = 0.85;
float eTPi0MinEB_     = 0.65;
float eTGammaMinEE_   = 0.250;
float s4s9GammaMinEE_ = 0.85;
float eTPi0MinEE_     = 0.800;
float swissCrossMaxEB_ = 0.95; // 1-E4/E1
float swissCrossMaxEE_ = 0.95; // 1-E4/E1
// based on range and bins, the bin width is 50 ps
//float rangeTDistro_ = 3; // 1-E4/E1
//int   binsTDistro_  = 120; // 1-E4/E1
std::vector<std::vector<double> > trigIncludeVector;
std::vector<std::vector<double> > trigExcludeVector;
std::vector<std::vector<double> > ttrigIncludeVector;
std::vector<std::vector<double> > ttrigExcludeVector;
// expected standard path for datacard when running e.g. on pcminn03
//std::string   datacardfile_=std::string("/data/franzoni/cmssw/44x/CMSSW_4_4_2_forTimeHardware/src/ECALTime/EcalTimePi0/data/DataCard.txt");
std::string   datacardfile_=std::string("/local/cms/user/norbert/EcalTiming/DataSet/DataCard.txt");

int  minEntriesForFit_ = 7;
int  flagOneVertex_ = 0;
bool limitFit_(true); 
//std::string fitOption_(""); // default: use chi2 method
std::string fitOption_("L"); // use likelihood method

// Ao dependent timing corrections
// By the definition of these corrections, the timing should be zero for the hits in Module 1 
// or Low eta EE within the valid A/sigma ranges.  Earlier data will have positive time due 
// to the graduate timing shifts in the positive direction.
TF1* timeCorrectionEB_ = new TF1("timeCorrectionEB_","pol4(0)",0,1.2);
TF1* timeCorrectionEE_ = new TF1("timeCorrectionEE_","pol4(0)",0,1.2);


// -------- Histograms -------------------------------------
TH1 * nVertices_;
TH1F* mass_;
TH1F* dZvertices_;
TH1F* Zvertices_;


// ---------------------------------------------------------------------------------------
// - Function to decide to include/exclude event based on the vectors passed for triggers 
bool includeEvent(int* triggers,
    int numTriggers,
    std::vector<std::vector<double> > includeVector,
    std::vector<std::vector<double> > excludeVector)
{
 bool keepEvent = false;
  if(includeVector.size()==0) keepEvent = true;
  for (int ti = 0; ti < numTriggers; ++ti) {
    for(uint i=0; i!=includeVector.size();++i){
      if(includeVector[i].size()==1 && triggers[ti]==includeVector[i][0]) keepEvent=true;
      else if(includeVector[i].size()==2 && (triggers[ti]>=includeVector[i][0] && triggers[ti]<=includeVector[i][1])) keepEvent=true;
    }
  }
  if(!keepEvent)
    return false;

  keepEvent = true;
  for (int ti = 0; ti < numTriggers; ++ti) {
    for(uint i=0; i!=excludeVector.size();++i){
      if(excludeVector[i].size()==1 && triggers[ti]==excludeVector[i][0]) keepEvent=false;
      else if(excludeVector[i].size()==2 && (triggers[ti]>=excludeVector[i][0] && triggers[ti]<=excludeVector[i][1])) keepEvent=false;
    }
  }

  return keepEvent;
}

// ---------------------------------------------------------------------------------------
// ------- Function to decide to include/exclude event based on the vectors passed -------
bool includeEvent(double eventParameter,
    std::vector<std::vector<double> > includeVector,
    std::vector<std::vector<double> > excludeVector)
{
  bool keepEvent = false;
  if(includeVector.size()==0) keepEvent = true;
  for(uint i=0; i!=includeVector.size();++i){
    if(includeVector[i].size()==1 && eventParameter==includeVector[i][0])
      keepEvent=true;
    else if(includeVector[i].size()==2 && (eventParameter>=includeVector[i][0] && eventParameter<=includeVector[i][1]))
      keepEvent=true;
  }
  if(!keepEvent) // if it's not in our include list, skip it
    return false;

  keepEvent = true;
  for(uint i=0; i!=excludeVector.size();++i){
    if(excludeVector[i].size()==1 && eventParameter==excludeVector[i][0])
      keepEvent=false;
    else if(excludeVector[i].size()==2 && (eventParameter>=excludeVector[i][0] && eventParameter<=excludeVector[i][1]))
      keepEvent=false;
  }

  return keepEvent; // if someone includes and excludes, exseedion will overrule

}



// ---------------------------------------------------------------------------------------
// ------------------ Function to split arg input lists by comma -------------------------
std::vector<std::string> split(std::string msg, std::string separator)
{
  boost::char_separator<char> sep(separator.c_str());
  boost::tokenizer<boost::char_separator<char> > tok(msg, sep );
  std::vector<std::string> token;
  for ( boost::tokenizer<boost::char_separator<char> >::const_iterator i = tok.begin(); i != tok.end(); ++i ) {
    token.push_back(std::string(*i));
  }
  return token;
}


// ---------------------------------------------------------------------------------------
// ------------------ Function to generate include/exclude vectors -----------------------
void genIncludeExcludeVectors(std::string optionString,
			      std::vector<std::vector<double> >& includeVector,
			      std::vector<std::vector<double> >& excludeVector)
{
  std::vector<std::string> rangeStringVector;
  std::vector<double> rangeIntVector;

  if(optionString != "-1"){
    std::vector<std::string> stringVector = split(optionString,",") ;

    for (uint i=0 ; i<stringVector.size() ; i++) {
      bool exclude = false;

      if(stringVector[i].at(0)=='x'){
        exclude = true;
        stringVector[i].erase(0,1);
      }
      rangeStringVector = split(stringVector[i],"-") ;

      rangeIntVector.clear();
      for(uint j=0; j<rangeStringVector.size();j++) {
        rangeIntVector.push_back(atof(rangeStringVector[j].c_str()));
      }
      if(exclude) excludeVector.push_back(rangeIntVector);
      else includeVector.push_back(rangeIntVector);

    }
  }
}


// ---------------------------------------------------------------------------------------
// ------------------ Function to parse the command-line arguments------------------------
void parseArguments(int argc, char** argv)
{
  std::string stringGenericOption    = "--";
  std::string stringHelp             = "--help";
  std::string stringInputFileName    = "--i";
  std::string stringOutFileName      = "--o";
  std::string stringETGammaMinEB     = "--eTGammaMinEB";
  std::string strings4s9GammaMinEB   = "--s4s9GammaMinEB";
  std::string stringeTPi0MinEB       = "--eTPi0MinEB";
  std::string stringETGammaMinEE     = "--eTGammaMinEE";
  std::string strings4s9GammaMinEE   = "--s4s9GammaMinEE";
  std::string stringeTPi0MinEE       = "--eTPi0MinEE";
  std::string stringminAmpliOverSigma= "--minAOverSigma";
  std::string stringNumEvents        = "--n";
  std::string stringMinRun           = "--minRun";
  std::string stringMaxRun           = "--maxRun";
  std::string stringMinLS            = "--minLS";
  std::string stringMaxLS            = "--maxLS";
  std::string vertex                 = "--vertex";
  std::string stringTriggers         = "--trig";
  std::string stringTechTriggers     = "--techTrig";
  std::string stringDataCard         = "--datacard";

  // if no arguments are passed, suggest help
  if (argc < 2){
    std::cerr << "\n\tERROR: specify arguments, or use --help\n" << std::endl ;
    exit (1) ;  
  }

  // loop over input options
  for (int v=1; v<argc; v++ )
  {
    //std::cout << "argv number " << v << " is: " << argv[v] << std::endl;

    if (argv[v] == stringHelp) { // help message
      std::cout << " --help : display help" << std::endl ;
      std::cout << " --o : set name of output root file name (e.g. histograms.root)" << std::endl ;
      std::cout << " --n : number of events" << std::endl ;
      std::cout << " --i <list of strings> list of input files" << std::endl ;     
      std::cout << " --eTGammaMinEB: min eT for EB gammas" << std::endl;
      std::cout << " --s4s9GammaMinEB: min EB shower shape" << std::endl;
      std::cout << " --eTPi0MinEB min eT for EB pi0 candidate" << std::endl;
      std::cout << " --eTGammaMinEE: min eT for EE gammas" << std::endl;
      std::cout << " --s4s9GammaMinEE: min EE shower shape" << std::endl;
      std::cout << " --eTPi0MinEE: min eT for EE pi0 candidate" << std::endl;
      std::cout << " --minAOverSigma: min ampli considered for time" << std::endl;
      std::cout << " --minRun: lowest run number considered" << std::endl;
      std::cout << " --maxRun: highest run number considered" << std::endl;
      std::cout << " --minLS: lowest lumi section number considered" << std::endl;
      std::cout << " --maxLS: highest lumi section number considered" << std::endl;
      std::cout << " --vertex: require vertex@IP (1), veto it (2) or either (0, or unset)" << std::endl;
      std::cout << " --trig: L1 triggers to include (exclude with x)" << std::endl;
      std::cout << " --techTrig: L1 technical triggers to include (exclude with x)" << std::endl;
      std::cout << " --datacard: data card holding Delayed Photon analysis selections" << std::endl;
      exit(1);      }


    else if (argv[v] == stringNumEvents) { // set number of events
      std::cout << "events number" << std::endl;
      numEvents_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMaxLS) { // set last LS of interval to be considered 
      std::cout << "max LS number" << std::endl;
      maxLS_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMinLS) { // set first LS of interval to be considered 
      std::cout << "min LS number" << std::endl;
      minLS_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMaxRun) { // set last run of interval to be considered 
      std::cout << "max LS number" << std::endl;
      maxRun_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringMinRun) { // set first run of interval to be considered 
      std::cout << "min run number" << std::endl;
      minRun_=atoi(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringETGammaMinEB) { // choose et cut for EB single cluster
      eTGammaMinEB_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == strings4s9GammaMinEB) { // choose cut for EB shower shape
      s4s9GammaMinEB_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringeTPi0MinEB) { // choose et cut for EB pi0 candidate
      eTPi0MinEB_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringETGammaMinEE) { // choose et cut for EE single cluster
      eTGammaMinEE_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == strings4s9GammaMinEE) { // choose cut for EE shower shape
      s4s9GammaMinEE_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringeTPi0MinEE) { // choose et cut for EE pi0 candidate
      eTPi0MinEE_ = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == stringOutFileName) { // set output file
      outputRootName_ = argv[v+1];
      v++;
    }
    else if (argv[v] == stringminAmpliOverSigma) { // set min amplitude considered for time measurement
      // removed as a option since this variable is not globally defined in EcalObjectTime.h
      // minAmpliOverSigma_  = atof(argv[v+1]);
      v++;
    }
    else if (argv[v] == vertex) { // collect requirement for one vertex only or not
      flagOneVertex_  = atof(argv[v+1]);
       if (flagOneVertex_!=0 && flagOneVertex_!=1 && flagOneVertex_!=2){
         std::cout << "Not a valid value for flagOneVertex_ (0,1,2). Returning." << std::endl;
	 exit (1);}
       v++;
    } 
    else if (argv[v] == stringTriggers) { // set L1 triggers to include/exclude
      genIncludeExcludeVectors(std::string(argv[v+1]),trigIncludeVector,trigExcludeVector);
      v++;
    }
    else if (argv[v] == stringTechTriggers) { // set L1 technical triggers to include/exclude
      genIncludeExcludeVectors(std::string(argv[v+1]),ttrigIncludeVector,ttrigExcludeVector);
      v++;
    }
    else if (argv[v] == stringDataCard) { // set datacard holding selections for DP analysis
      datacardfile_ = std::string(argv[v+1]);
      v++;
    }
    // handle here the case of multiple arguments for input files
    else if (argv[v] == stringInputFileName){// && v<(argc-1) ) 

      for (int u=v+1; u<argc; u++) {

        if ( 0==std::string(argv[u]).find( stringGenericOption ) ){
          if ( 0==listOfFiles_.size())  {std::cout << "no input files listed" << std::cout;}
          //else  {std::cout << "no more files listed, found: " << argv[u] << std::cout;}
          break;
        }

        else {  listOfFiles_.push_back(argv[u]);
          v++;
        }

      }// loop on arguments following --i

      continue;

    }//end 'if input files'

    else
    {std::cout << "input format unrecognized" << std::endl; exit(1);}

    }// loop over arguments input to the program
}


// ---------------------------------------------------------------------------------------
// ------------------ Function to initialize the histograms ------------------------------
void initializeHists(TFileDirectory subDir){

  mass_         = subDir.make<TH1F>("mass global","mass (global); m(ele,ele) [GeV]",80,50,130);
  dZvertices_   = subDir.make<TH1F>("dZvertices global","dZvertices (global); #DeltaZ(ele_{1},ele_{2}) [cm]",250,0,25);
  Zvertices_    = subDir.make<TH1F>("Zvertices global","Zvertices (global); z vertex [cm]",250,-25,25);
  nVertices_=subDir.make<TH1F>("num vertices global","num vertices (global); num vertices",41,-0.5,40.5);

}//end initializeHists




// ---------------------------------------------------------------------------------------
//! main program
int main (int argc, char** argv)
{
  // First parse arguments
  parseArguments(argc, argv);

  if (listOfFiles_.size()==0){
    std::cout << "\tno input file found" << std::endl;
    return(1);
  }
  else{
    std::cout << "\tfound " << listOfFiles_.size() << " input files: " << std::endl;
    for(std::vector<std::string>::const_iterator  file_itr=listOfFiles_.begin(); file_itr!=listOfFiles_.end(); file_itr++){
      std::cout << "\t" << (*file_itr) << std::endl;
    }
  }

  // Tree construction
  // FIX should turn this string into a configurable 
  TChain * chain            = new TChain ("EcalTimeAnalysis") ;        // ntuple producer in CMSSW CVS
  //TChain * chainForSelector = new TChain ("EcalTimeAnalysisForSel") ;  // ntuple producer in CMSSW CVS
  std::vector<std::string>::const_iterator file_itr;
  for(file_itr=listOfFiles_.begin(); file_itr!=listOfFiles_.end(); file_itr++){
    chain           ->Add( (*file_itr).c_str() );
    //chainForSelector->Add( (*file_itr).c_str() );
  }
  // a separate instane for SC's DP selector
  TChain * chainForSelector = (TChain *) chain->Clone();
  int nEntries = chain->GetEntries () ;

  if (numEvents_==-1) numEvents_ = nEntries;
  std::cout << "\n\tFOUND "         <<  listOfFiles_.size() << " input files" << std::endl ;    
  std::cout << "\n\tFOUND "         <<  nEntries << " events" << std::endl ;    
  std::cout << "\tWILL run on: "    <<  numEvents_ << " events" << std::endl;
  std::cout << "\tOutput file: "    <<  outputRootName_ << std::endl;
  std::cout << "\tminAOverSigma: "  <<  minAmpliOverSigma_ << std::endl;
  std::cout << "\teTGammaMinEB: "   <<  eTGammaMinEB_ << std::endl;
  std::cout << "\ts4s9GammaMinEB: " <<  s4s9GammaMinEB_ << std::endl;
  std::cout << "\teTPi0MinEB: "     <<  eTPi0MinEB_ << std::endl;
  std::cout << "\teTGammaMinEE: "   <<  eTGammaMinEE_ << std::endl;
  std::cout << "\ts4s9GammaMinEE: " <<  s4s9GammaMinEE_ << std::endl;
  std::cout << "\teTPi0MinEE: "     <<  eTPi0MinEE_ << std::endl;
  std::cout << "\tminRun: "         <<  minRun_ << std::endl;
  std::cout << "\tmaxRun: "         <<  maxRun_ << std::endl;
  std::cout << "\tminLS: "          <<  minLS_ << std::endl;
  std::cout << "\tmaxLS: "          <<  maxLS_ << std::endl;
  std::cout << "\tdatacard: "       <<  datacardfile_ << std::endl;
  
  // establish link between ttree variables and local variables
  setBranchAddresses (chain, treeVars_);

  Input_    = new AnaInput( datacardfile_ );
  Selector_ = new DPSelection( datacardfile_ ) ;
  int nEntriesSelector = chainForSelector->GetEntries () ;
  std::cout << "nEntriesSelector: " << nEntriesSelector << " nEntries: " << nEntries << std::endl;

  // link the variables for selection module
  Selector_->Init( chainForSelector ) ;

  string hfolder ;
  string plotType ;
  int ProcessEvents ;
  string debugStr ;
  int selectBackground ;
  int split ;
  int counter[10] ;
  for (int i=0; i<10; i++)  counter[i] = 0 ;

  Input_->GetParameters("PlotType",      &plotType ) ; 
  Input_->GetParameters("Path",          &hfolder ) ; 
  Input_->GetParameters("ProcessEvents", &ProcessEvents ) ; 

  Input_->GetParameters( "SelectBackground",  &selectBackground );
  Input_->GetParameters( "SplitEvent",        &split );
  
  // setting up the TFileService in the ServiceRegistry;
  edmplugin::PluginManager::Config config;
  edmplugin::PluginManager::configure(edmplugin::standard::config());
  std::vector<edm::ParameterSet> psets;
  edm::ParameterSet pSet;
  pSet.addParameter("@service_type",std::string("TFileService"));
  pSet.addParameter("fileName",outputRootName_); // this is the file TFileService will write into    //pSet.addParameter("fileName",std::string("TimePerf-plots.root"));
  psets.push_back(pSet);
  static edm::ServiceToken services(edm::ServiceRegistry::createSet(psets));
  static edm::ServiceRegistry::Operate operate(services);
  edm::Service<TFileService> fs;

  TFileDirectory subDirECALECAL=fs->mkdir("ECALECAL");  
  HistSet plotsECALECAL;   plotsECALECAL.book(subDirECALECAL,std::string("ECALECAL"));
  plotsECALECAL.setTree(&treeVars_);

  TFileDirectory subDirEB=fs->mkdir("EB");  
  HistSet plotsEB;   plotsEB.book(subDirEB,std::string("EB"));
  plotsEB.setTree(&treeVars_);


  TFileDirectory subDirEBchi2=fs->mkdir("EBchi2");  
  HistSet plotsEBchi2;   plotsEBchi2.book(subDirEBchi2,std::string("EBchi2"));
  plotsEBchi2.setTree(&treeVars_);


  TFileDirectory subDirEBchi2loose=fs->mkdir("EBchi2loose");  
  HistSet plotsEBchi2loose;   plotsEBchi2loose.book(subDirEBchi2loose,std::string("EBchi2loose"));
  plotsEBchi2loose.setTree(&treeVars_);


  TFileDirectory subDirEBchi2tight=fs->mkdir("EBchi2tight");  
  HistSet plotsEBchi2tight;   plotsEBchi2tight.book(subDirEBchi2tight,std::string("EBchi2tight"));
  plotsEBchi2tight.setTree(&treeVars_);


  TFileDirectory subDirEBEE=fs->mkdir("EBEE");  
  HistSet plotsEBEE;   plotsEBEE.book(subDirEBEE,std::string("EBEE"));
  plotsEBEE.setTree(&treeVars_);

  timeCorrector theCorr;
  std::cout << "\ncreated object theCorr to be used for timeVsAmpliCorrections" << std::endl;
  std::cout << "\ninitializing theCorr" << std::endl;
  //theCorr.initEB( std::string("EBmod4") );
  //theCorr.initEE( std::string("EElow") );
  theCorr.initEB( "EB" );
  theCorr.initEE( "EElow" );


  //Initialize output root file
  //saving_ = new TFile(outputRootName_.c_str (),"recreate");

  // Initialize the histograms
  TFileDirectory subDirGeneral=fs->mkdir("General");  
  initializeHists(subDirGeneral);

  int eventCounter     = 0;
  int eventPassCounter = 0;

  /////////////////////////////////////////////////////  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////  /////////////////////////////////////////////////////
  // Main loop over entries
  for (int entry = 0 ; (entry < nEntries && eventCounter < numEvents_); ++entry)
  {
    chain            ->GetEntry (entry) ;
    chainForSelector ->GetEntry (entry) ;
    // Keep the event?
    bool keepEvent = includeEvent(treeVars_.l1ActiveTriggers,
        treeVars_.l1NActiveTriggers,trigIncludeVector,trigExcludeVector)
            && includeEvent(treeVars_.l1ActiveTechTriggers,
                treeVars_.l1NActiveTechTriggers,ttrigIncludeVector,ttrigExcludeVector);
    if(!keepEvent)
      continue;

    
    // do analysis if the run is in the desired range  
    if( treeVars_.runId<minRun_  || maxRun_<treeVars_.runId) continue;
    
    // do analysis if the LS is in the desired range  
    if( treeVars_.lumiSection<minLS_  || maxLS_<treeVars_.lumiSection) continue;
    
    // if evet being actually processed, increment counter of analyzed events
    eventCounter++;
/*
    /////////////////// DP analysis events selection ////////////////////////
    bool passEvent = true ;

    // FIXME
    if ( entry%2 == split ) continue ;    // FIXME
    // FIXME
    counter[0]++ ;   
    // put back the cuts to be the same as they are in the datacard
    Selector_->ResetCollection() ;
    Selector_->ResetCuts() ;

    bool passPhoton = Selector_->PhotonFilter( false );
    // DEBUG //    //std::cout << "passPhoton Loose No Iso: " << passPhoton << std::endl; // DEBUG //
    passEvent = ( passPhoton && passEvent ) ? true : false ;
    if ( passEvent ) counter[1]++ ;   

    bool passVertex = Selector_->VertexFilter();
    passEvent = ( passVertex && passEvent ) ? true : false ;
    if ( passEvent ) counter[2]++ ;   

    Selector_->ResetCuts( "PhotonCuts", 1, 1.4 ) ;
    bool passIso = Selector_->PhotonFilter( true );
    passEvent = ( passIso && passEvent ) ? true : false ;
    if ( passEvent ) counter[3]++ ;   

    bool isGJets = Selector_->GammaJetsBackground() ;
    passEvent = ( !isGJets && passEvent ) ? true : false ;
    if ( passEvent ) counter[4]++ ;   

    Selector_->ResetCuts( "PhotonCuts", 0, 100. ) ;
    bool passTight = Selector_->PhotonFilter( true );
    passEvent = ( passTight && passEvent ) ? true : false ;
    if ( passEvent ) counter[5]++ ;   

    bool passJet    = Selector_->JetMETFilter();
    passEvent = ( passJet && passEvent ) ? true : false ;
    if ( passEvent ) counter[6]++ ;   

    bool passHLT    = Selector_->HLTFilter();
    passEvent = ( passHLT && passEvent ) ? true : false ;
    if ( passEvent ) counter[7]++ ;   
    /////////////////// DP analysis events selection ////////////////////////

   
    // if selection according to datacard not passed, move on to next event 
    if ( ! passEvent ) continue;    
    
    // define samples here /////////////////////////////////// (to be updated/FIXED including SC's implementation )
    // control region
    // 1: gamma+jets background control sample - FIXME
    if( selectBackground==1 && (treeVars_.nJets<1 || treeVars_.nJets>2) ) continue;
    // 2: QCD background control sample - FIXME
    if( selectBackground==2 && (treeVars_.nJets<1 || treeVars_.nJets>2) ) continue;
 */
     // look only at half of the events 
    if ( (treeVars_.orbit % 2)==1 ) continue;

    bool verticesAreOnlyNextToNominalIP;
    int  count=0;
    
    for(int v=0; v<treeVars_.nVertices; v++  )
	//{        if (fabs(treeVars_.vtxZ[0])<15) count++; }
	{        if (fabs(treeVars_.vtxZ[0])<20) count++; }
    
    if ( treeVars_.nVertices >0 && count==treeVars_.nVertices ) verticesAreOnlyNextToNominalIP = true;
    else                                                        verticesAreOnlyNextToNominalIP = false;
    
    //    --vertex: require vertex@IP (1), veto it (2) or either (0, or unset)
    if (flagOneVertex_ ==1 && (!verticesAreOnlyNextToNominalIP) ) continue;
    if (flagOneVertex_ ==2 && (verticesAreOnlyNextToNominalIP) )  continue;
    
    eventPassCounter++;
    
    speak_=false;
    if (entry<10 || entry%10000==0) speak_=true;

    if (speak_)  std::cout << "\n\n------> reading entry " << entry << "\tLS: " << treeVars_.lumiSection << " <------\n" ; 
    if (speak_)  std::cout << "  found " << treeVars_.nSuperClusters << " superclusters" << std::endl ;
    if (speak_)  std::cout << "  found " << treeVars_.nClusters << " basic clusters" << std::endl ;

    nVertices_->Fill(treeVars_.nVertices);
    
    ///////////////////////////////////////////////////////////////////////
    // outer loop on supercluster
    for (int sc1=0; sc1<treeVars_.nSuperClusters; sc1++){

      // require object FIRST, since for jets there may not BE superClusterRawEnergy or superClusterEta 
  /*  if ( fabs ( treeVars_.SCPIdx[sc1] - 22) > 2) continue; */    

      float et1 = treeVars_.superClusterRawEnergy[sc1]/cosh( treeVars_.superClusterEta[sc1] );
      //std::cout << "gf1 num SC n SC" << treeVars_.nSuperClusters << "\t et: " << et1 << "\t" << treeVars_.superClusterRawEnergy[sc1] << "\t" << treeVars_.superClusterEta[sc1] << "\t" << treeVars_.SCPIdx[sc1] << std::endl; 
      if (et1<70) continue; // match choice of control region; this should be translated into cut on the photon, not SC
      // select only photons

      float tmpEne=-9999;
	// loop on BC and match to sc1  ===============
	int bc1=-1;
	for (int bc=0; bc<treeVars_.nClusters; bc++){
	  if ( (pow(treeVars_.superClusterEta[sc1]-treeVars_.clusterEta[bc],2)+ pow(treeVars_.superClusterPhi[sc1]-treeVars_.clusterPhi[bc],2) ) < 0.02 
	       && treeVars_.clusterEnergy[bc]>tmpEne) {
	    tmpEne=treeVars_.clusterEnergy[bc];
	    bc1=bc;
	  }// end - if good bc candidate
	}// end - loop over BC

	
	// protect in case of no matching
	if(bc1==-1) continue;
	
	ClusterTime bcTime1 = timeAndUncertSingleCluster(bc1,treeVars_);
	if(! (bcTime1.isvalid ) ) continue;

	if      ( fabs(treeVars_.clusterEta[bc1])<1.4   ){
	  plotsEB.fillSingle(sc1, bc1,  bcTime1);
	  
	  int type=3; float cut=6;
	  plotsEBchi2loose.fillSingle(sc1, bc1,  bcTime1, type, cut);

	  type=3; cut=4;
	  plotsEBchi2.fillSingle(sc1, bc1,  bcTime1, type, cut);

	  type=3; cut=3;
	  plotsEBchi2tight.fillSingle(sc1, bc1,  bcTime1, type, cut);

	} 
	else if( fabs(treeVars_.clusterEta[bc1])>1.5 && fabs(treeVars_.clusterEta[bc1])<2.5 ){
	plotsEBEE.fillSingle(sc1, bc1,  bcTime1);
	}

	break;	
	
    }// end loop sc1
    
  }   // end of loop over entries
  
  /////////////////// summary of DP analysis events selection ////////////////////////
  std::cout<<"\n final event summary -  all:"<< counter[0] <<" pho:"<<counter[1]<<" vtx:"<< counter[2]<<" Iso:"<<counter[3]<<" noGjets:"<<counter[4] ;
  std::cout<<" tight:"<< counter[5] <<" jet:"<<counter[6]<<" hlt:"<<counter[7]<< " - events that passed preselection;     inclusive plots: " << eventPassCounter << " (that have passed the Sample type you want: sig, control, control QCD )" <<std::endl; 
  
  delete chain ;
  
  return 0 ;
}
