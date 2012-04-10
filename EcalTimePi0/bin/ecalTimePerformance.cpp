#include <iostream>
#include <string>

#include <map>
#include <vector>
#include <algorithm> 
#include <functional>
#include <set>
#include <boost/tokenizer.hpp>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/PluginManager/interface/standard.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// following is commented out because already included through analysisHisograms.h
//#include "ECALTime/EcalTimePi0/interface/EcalObjectTime.h"
//included via EcalObjectTime.h  ->  analysisHisograms.h
//#include "CalibCalorimetry/EcalTiming/interface/EcalTimeTreeContent.h"

#include "ECALTime/EcalTimePi0/interface/analysisHistograms.h"
#include "ECALTime/EcalTimePi0/interface/timeVsAmpliCorrector.h"


#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"


typedef std::set<std::pair<int,int> > SetOfIntPairs;
//float travelDistance(int sc_num, EcalTimeTreeContent treeVars_);
//float extraTravelTime(int sc_num, EcalTimeTreeContent & treeVars_);
//float extraTravelTime(int sc_num, int vtx_num, EcalTimeTreeContent & treeVars_);

// authors: S. Cooper and G. Franzoni (UMN)

#define BarrelLimit  1.479
#define EndcapLimit  3.0

#define ADCtoGeVEB   0.039
#define ADCtoGeVEE   0.063

#define numAeffBins     35
#define numAoSigmaBins  25

#define NSlices 54   // # of slices in log(A)
#define LogStep 0.05 // size of Slices log(A)


// -------- Globals ----------------------------------------
EcalTimePhyTreeContent treeVars_; 
TFile* saving_;
std::vector<std::string> listOfFiles_;
bool speak_=false;
char buffer_ [75];
std::string bufferTitle_; 
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

std::vector<std::vector<double> > trigIncludeVector;
std::vector<std::vector<double> > trigExcludeVector;
std::vector<std::vector<double> > ttrigIncludeVector;
std::vector<std::vector<double> > ttrigExcludeVector;

int  minEntriesForFit_ = 7;
int  flagOneVertex_ = 0;
int  stringDoAveragePerRun_ = 1;
bool limitFit_(true); 
//std::string fitOption_(""); // default: use chi2 method
std::string fitOption_("L"); // use likelihood method

// Ao dependent timing corrections
// By the definition of these corrections, the timing should be zero for the hits in Module 1 
// or Low eta EE within the valid A/sigma ranges.  Earlier data will have positive time due 
// to the graduate timing shifts in the positive direction.
TF1* timeCorrectionEB_ = new TF1("timeCorrectionEB_","pol4(0)",0,1.2);
TF1* timeCorrectionEE_ = new TF1("timeCorrectionEE_","pol4(0)",0,1.2);

double AeffBins_[36] = {0,                          // set of fine bins for large stats
 			1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,
 			11,12,13,14,15,16,17,18,19,20,
 			22,24,26,28,30,
                        32,36,40,46,52,
                        58,74,86,110,130};
int    AeffNBins_    = 35;
int    AeffMax_      = 120;

// defining arrays large enough; number of actual bins to define TGraphErrors comes from AeffNBins_ 
float  AeffBinCentersAny_[256]; float  AeffBinCentersErrAny_[256]; float  sigmaAeffAny_[256]; float  sigmaAeffErrAny_[256];
float  AeffBinCentersEB_[256];  float  AeffBinCentersErrEB_[256];  float  sigmaAeffEB_[256];  float  sigmaAeffErrEB_[256];
float  AeffBinCentersEE_[256];  float  AeffBinCentersErrEE_[256];  float  sigmaAeffEE_[256];  float  sigmaAeffErrEE_[256];

// set of fine bins for AoSigmaBins, for t vs Ampli bias study 
double AoSigmaBins_[33] = {0,  10,  18,  
			   24, 28,  32,  36,  40,
			   44, 52,  60,  72,  80,
                           92, 102, 116, 144, 172,
                           240,320, 400, 480, 560, 
			   650, 800, 1000, 1200, 1500,
                           1900, 2400, 3000, 3700, 4000}; //25 bins in total
int    AoSigmaNBins_     = 32;    // (counting to be updated) 300 combinations between different bins; + 25 self-combinations =-> 325 in total
                                  // check above that numAeffBins is set to a value equal or larger than this (==325)
int    AoSigmaNPairs_    = (AoSigmaNBins_)*(AoSigmaNBins_-1)/2 + AoSigmaNBins_;
int    AoSigmaMax_       = 4000;   //up to about 8 GeV

float  AoSigmaBinCentersEB_[32][32];  float  AoSigmaBinCentersErrEB_[32][32];  float  sigmaAoSigmaEB_[32][32];  float  sigmaAoSigmaErrEB_[32][32];
float  AoSigmaBinCentersEE_[32][32];  float  AoSigmaBinCentersErrEE_[32][32];  float  sigmaAoSigmaEE_[32][32];  float  sigmaAoSigmaErrEE_[32][32];

int numDtBins_  = 75;
int DtMax_      = 15; // useful to catch tails also at low Aeff (<10)

// -------- Histograms -------------------------------------
TH1 * nVertices_;
TH1F* mass_;
TH1F* dZvertices_;
TH1F* Zvertices_;

TGraphErrors *ebTimeVsRun;
TGraphErrors *eeTimeVsRun;
TGraphErrors *eeTimeVSebTime;
float        ebTimes[1000];
float        ebTimeErrs[1000];
float        eeTimes[1000];
float        eeTimeErrs[1000];
float        runs[1000];
float        runErrs[1000];
int          runCounter;

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
  std::string stringDoAveragePerRun  = "--stringDoAveragePerRun";

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
      std::cout << " --stringDoAveragePerRun: do (1) or don't (0) subtrct average time by run" << std::endl;
      std::cout << " --trig: L1 triggers to include (exclude with x)" << std::endl;
      std::cout << " --techTrig: L1 technical triggers to include (exclude with x)" << std::endl;
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
    else if (argv[v] == stringDoAveragePerRun) { // collect requirement for whether averate run time should be subtracted or not
      stringDoAveragePerRun_  = atof(argv[v+1]);
       if (stringDoAveragePerRun_!=0 && stringDoAveragePerRun_!=1 ){
         std::cout << "Not a valid value for stringDoAveragePerRun_ (0,1). Returning." << std::endl;
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
  
  // initialize TGraphErrors-related stuff
  runCounter = 0;
  for (int v=0; v<1000; v++){
    ebTimes[v]=0; ebTimeErrs[v]=0;
    eeTimes[v]=0; eeTimeErrs[v]=0;
    runs[v]=0;    runErrs[v]=0;
  }

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
  TChain * chain = new TChain ("EcalTimeAnalysis") ;  // ntuple producer in CMSSW CVS
  //TChain * chain = new TChain ("EcalTimePi0Analysis") ;  // ntuple producer in UserCode/UMN space
  std::vector<std::string>::const_iterator file_itr;
  for(file_itr=listOfFiles_.begin(); file_itr!=listOfFiles_.end(); file_itr++){
    chain->Add( (*file_itr).c_str() );
  }
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
	
  setBranchAddresses (chain, treeVars_);
  
  // setting up the TFileService in the ServiceRegistry;
  edmplugin::PluginManager::Config config;
  edmplugin::PluginManager::configure(edmplugin::standard::config());
  std::vector<edm::ParameterSet> psets;
  edm::ParameterSet pSet;
  pSet.addParameter("@service_type",std::string("TFileService"));
  pSet.addParameter("fileName",outputRootName_); // this is the file TFileService will write into
  psets.push_back(pSet);
  static edm::ServiceToken services(edm::ServiceRegistry::createSet(psets));
  static edm::ServiceRegistry::Operate operate(services);
  edm::Service<TFileService> fs;
  // Initialize the histograms
  TFileDirectory subDirGeneral=fs->mkdir("General");  
  initializeHists(subDirGeneral);

 
  TFileDirectory subDirECALECAL=fs->mkdir("ECALECAL");  
  HistSet plotsECALECAL;   plotsECALECAL.book(subDirECALECAL,std::string("ECALECAL"));
  plotsECALECAL.setTree(&treeVars_);
  
  TFileDirectory subDirEBEB=fs->mkdir("EBEB");  
  HistSet plotsEBEB;   plotsEBEB.book(subDirEBEB,std::string("EBEB"));
  plotsEBEB.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBnear=fs->mkdir("EBEBnear");  
  HistSet plotsEBEBnear;   plotsEBEBnear.book(subDirEBEBnear,std::string("EBEBnear"));
  plotsEBEBnear.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBEtaNear=fs->mkdir("EBEBEtaNear");  
  HistSet plotsEBEBEtaNear;   plotsEBEBEtaNear.book(subDirEBEBEtaNear,std::string("EBEBEtaNear"));
  plotsEBEBEtaNear.setTree(&treeVars_);

  TFileDirectory subDirEBEBfar=fs->mkdir("EBEBfar");  
  HistSet plotsEBEBfar;   plotsEBEBfar.book(subDirEBEBfar,std::string("EBEBfar"));
  plotsEBEBfar.setTree(&treeVars_);

  TFileDirectory subDirEBEBVeryFar=fs->mkdir("EBEBVeryFar");  
  HistSet plotsEBEBVeryFar;   plotsEBEBVeryFar.book(subDirEBEBVeryFar,std::string("EBEBVeryFar"));
  plotsEBEBVeryFar.setTree(&treeVars_);

  TFileDirectory subDirEBEBforward=fs->mkdir("EBEBforward");  
  HistSet plotsEBEBforward;   plotsEBEBforward.book(subDirEBEBforward,std::string("EBEBforward"));
  plotsEBEBforward.setTree(&treeVars_);

  TFileDirectory subDirEBEBmod1=fs->mkdir("EBEBmod1");  
  HistSet plotsEBEBmod1;   plotsEBEBmod1.book(subDirEBEBmod1,std::string("EBEBmod1"));
  plotsEBEBmod1.setTree(&treeVars_);

  TFileDirectory subDirEBEBmod2=fs->mkdir("EBEBmod2");  
  HistSet plotsEBEBmod2;   plotsEBEBmod2.book(subDirEBEBmod2,std::string("EBEBmod2"));
  plotsEBEBmod2.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBmod3=fs->mkdir("EBEBmod3");  
  HistSet plotsEBEBmod3;   plotsEBEBmod3.book(subDirEBEBmod3,std::string("EBEBmod3"));
  plotsEBEBmod3.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBmod4=fs->mkdir("EBEBmod4");  
  HistSet plotsEBEBmod4;   plotsEBEBmod4.book(subDirEBEBmod4,std::string("EBEBmod4"));
  plotsEBEBmod4.setTree(&treeVars_);
  
  TFileDirectory subDirEEEEnear=fs->mkdir("EEEEnear");  
  HistSet plotsEEEEnear;   plotsEEEEnear.book(subDirEEEEnear,std::string("EEEEnear"));
  plotsEEEEnear.setTree(&treeVars_);
  
  TFileDirectory subDirEEEEfar=fs->mkdir("EEEEfar");  
  HistSet plotsEEEEfar;   plotsEEEEfar.book(subDirEEEEfar,std::string("EEEEfar"));
  plotsEEEEfar.setTree(&treeVars_);
  
  TFileDirectory subDirEEEEVeryFar=fs->mkdir("EEEEVeryFar");  
  HistSet plotsEEEEVeryFar;   plotsEEEEVeryFar.book(subDirEEEEVeryFar,std::string("EEEEVeryFar"));
  plotsEEEEVeryFar.setTree(&treeVars_);
  
  TFileDirectory subDirEEEEEtaNear=fs->mkdir("EEEEEtaNear");  
  HistSet plotsEEEEEtaNear;   plotsEEEEEtaNear.book(subDirEEEEEtaNear,std::string("EEEEEtaNear"));
  plotsEEEEEtaNear.setTree(&treeVars_);
  
  TFileDirectory subDirEBMEEM=fs->mkdir("EBMEEM");  
  HistSet plotsEBMEEM;   plotsEBMEEM.book(subDirEBMEEM,std::string("EBMEEM"));
  plotsEBMEEM.setTree(&treeVars_);
  
  TFileDirectory subDirEBPEEP=fs->mkdir("EBPEEP");  
  HistSet plotsEBPEEP;   plotsEBPEEP.book(subDirEBPEEP,std::string("EBPEEP"));
  plotsEBPEEP.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBzsmall=fs->mkdir("EBEBzsmall");  
  HistSet plotsEBEBzsmall;   plotsEBEBzsmall.book(subDirEBEBzsmall,std::string("EBEBzsmall"));
  plotsEBEBzsmall.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBzbigP=fs->mkdir("EBEBzbigP");  
  HistSet plotsEBEBzbigP;   plotsEBEBzbigP.book(subDirEBEBzbigP,std::string("EBEBzbigP"));
  plotsEBEBzbigP.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBzbigN=fs->mkdir("EBEBzbigN");  
  HistSet plotsEBEBzbigN;   plotsEBEBzbigN.book(subDirEBEBzbigN,std::string("EBEBzbigN"));
  plotsEBEBzbigN.setTree(&treeVars_);
  
  // separate folders for different number of reconstructed vertices
  TFileDirectory subDirEBEBlowPU=fs->mkdir("EBEBlowPU");  
  HistSet plotsEBEBlowPU;   plotsEBEBlowPU.book(subDirEBEBlowPU,std::string("EBEBlowPU"));
  plotsEBEBlowPU.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBmidPU=fs->mkdir("EBEBmidPU");  
  HistSet plotsEBEBmidPU;   plotsEBEBmidPU.book(subDirEBEBmidPU,std::string("EBEBmidPU"));
  plotsEBEBmidPU.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBhighPU=fs->mkdir("EBEBhighPU");  
  HistSet plotsEBEBhighPU;   plotsEBEBhighPU.book(subDirEBEBhighPU,std::string("EBEBhighPU"));
  plotsEBEBhighPU.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBsuperPU=fs->mkdir("EBEBsuperPU");  
  HistSet plotsEBEBsuperPU;   plotsEBEBsuperPU.book(subDirEBEBsuperPU,std::string("EBEBsuperPU"));
  plotsEBEBsuperPU.setTree(&treeVars_);
  
  TFileDirectory subDirEEEE=fs->mkdir("EEEE");  
  HistSet plotsEEEE;   plotsEEEE.book(subDirEEEE,std::string("EEEE"));
  plotsEEEE.setTree(&treeVars_);
  
  TFileDirectory subDirEEEEzsmall=fs->mkdir("EEEEzsmall");  
  HistSet plotsEEEEzsmall;   plotsEEEEzsmall.book(subDirEEEEzsmall,std::string("EEEEzsmall"));
  plotsEEEEzsmall.setTree(&treeVars_);

  TFileDirectory subDirEEEEzbigP=fs->mkdir("EEEEzbigP");  
  HistSet plotsEEEEzbigP;   plotsEEEEzbigP.book(subDirEEEEzbigP,std::string("EEEEzbigP"));
  plotsEEEEzbigP.setTree(&treeVars_);
  
  TFileDirectory subDirEEEEzbigN=fs->mkdir("EEEEzbigN");  
  HistSet plotsEEEEzbigN;   plotsEEEEzbigN.book(subDirEEEEzbigN,std::string("EEEEzbigN"));
  plotsEEEEzbigN.setTree(&treeVars_);
    
  TFileDirectory subDirEEPEEP=fs->mkdir("EEPEEP");  
  HistSet plotsEEPEEP;   plotsEEPEEP.book(subDirEEPEEP,std::string("EEPEEP"));
  plotsEEPEEP.setTree(&treeVars_);
  
  TFileDirectory subDirEEMEEM=fs->mkdir("EEMEEM");  
  HistSet plotsEEMEEM;   plotsEEMEEM.book(subDirEEMEEM,std::string("EEMEEM"));
  plotsEEMEEM.setTree(&treeVars_);
    
  TFileDirectory subDirEEPEEM=fs->mkdir("EEPEEM");  
  HistSet plotsEEPEEM;   plotsEEPEEM.book(subDirEEPEEM,std::string("EEPEEM"));
  plotsEEPEEM.setTree(&treeVars_);
    
  TFileDirectory subDirEBEE=fs->mkdir("EBEE");  
  HistSet plotsEBEE;   plotsEBEE.book(subDirEBEE,std::string("EBEE"));
  plotsEBEE.setTree(&treeVars_);
    
  TFileDirectory subDirEBEBequalShare=fs->mkdir("EBEBequalShare");  
  HistSet plotsEBEBequalShare;   plotsEBEBequalShare.book(subDirEBEBequalShare,std::string("EBEBequalShare"));
  plotsEBEBequalShare.setTree(&treeVars_);
    
  TFileDirectory subDirEBEBunevenShare=fs->mkdir("EBEBunevenShare");  
  HistSet plotsEBEBunevenShare;   plotsEBEBunevenShare.book(subDirEBEBunevenShare,std::string("EBEBunevenShare"));
  plotsEBEBunevenShare.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBchi2=fs->mkdir("EBEBchi2");  
  HistSet plotsEBEBchi2;   plotsEBEBchi2.book(subDirEBEBchi2,std::string("EBEBchi2"));
  plotsEBEBchi2.setTree(&treeVars_);

  TFileDirectory subDirEBEEchi2=fs->mkdir("EBEEchi2");  
  HistSet plotsEBEEchi2;   plotsEBEEchi2.book(subDirEBEEchi2,std::string("EBEEchi2"));
  plotsEBEEchi2.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBchi2loose=fs->mkdir("EBEBchi2loose");  
  HistSet plotsEBEBchi2loose;   plotsEBEBchi2loose.book(subDirEBEBchi2loose,std::string("EBEBchi2loose"));
  plotsEBEBchi2loose.setTree(&treeVars_);

  TFileDirectory subDirEBEEchi2loose=fs->mkdir("EBEEchi2loose");  
  HistSet plotsEBEEchi2loose;   plotsEBEEchi2loose.book(subDirEBEEchi2loose,std::string("EBEEchi2loose"));
  plotsEBEEchi2loose.setTree(&treeVars_);

  TFileDirectory subDirEBEBchi2tight=fs->mkdir("EBEBchi2tight");  
  HistSet plotsEBEBchi2tight;   plotsEBEBchi2tight.book(subDirEBEBchi2tight,std::string("EBEBchi2tight"));
  plotsEBEBchi2tight.setTree(&treeVars_);

  TFileDirectory subDirEBEEchi2tight=fs->mkdir("EBEEchi2tight");  
  HistSet plotsEBEEchi2tight;   plotsEBEEchi2tight.book(subDirEBEEchi2tight,std::string("EBEEchi2tight"));
  plotsEBEEchi2tight.setTree(&treeVars_);
  
  TFileDirectory subDirEBEBseed2sec=fs->mkdir("EBEBseed2sec");  
  HistSet plotsEBEBseed2sec;   plotsEBEBseed2sec.book(subDirEBEBseed2sec,std::string("EBEBseed2sec"));
  plotsEBEBseed2sec.setTree(&treeVars_);

  TFileDirectory subDirEBEEseed2sec=fs->mkdir("EBEEseed2sec");  
  HistSet plotsEBEEseed2sec;   plotsEBEEseed2sec.book(subDirEBEEseed2sec,std::string("EBEEseed2sec"));
  plotsEBEEseed2sec.setTree(&treeVars_);

  int eventCounter = 0;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // auxiliary loop over the events to find time average per run and avearge Z
  std::map <int, std::pair<float,float> >  averageTimePerRun;   // runNumber, (timeEB   , timeEE)
  std::map <int, std::pair<float,float> >  averageTimeSqPerRun; // runNumber, (time^2EB , time^2EE)
  std::map <int, std::pair<int,int> >      countersPerRun;      // runNumber, (count EB , count EE)
  
  std::map  <int, std::pair<float,float> >::iterator run_time;
  std::map <int, std::pair<float,float> >::iterator  run_timeSq;
  std::map <int, std::pair<int,int> >::iterator      run_count;

  for (int entry = 0 ; (entry < nEntries && eventCounter < numEvents_); ++entry)
    {
      chain->GetEntry (entry) ;
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
      
      bool verticesAreOnlyNextToNominalIP;
      int  count=0;
      
      for(int v=0; v<treeVars_.nVertices; v++  )
	{        if (fabs(treeVars_.vtxZ[0])<15) count++; }
    
      if ( treeVars_.nVertices >0 && count==treeVars_.nVertices ) verticesAreOnlyNextToNominalIP = true;
      else                                                        verticesAreOnlyNextToNominalIP = false;
      
      //    --vertex: require vertex@IP (1), veto it (2) or either (0, or unset)
      if (flagOneVertex_ ==1 && (!verticesAreOnlyNextToNominalIP) ) continue;
      if (flagOneVertex_ ==2 && (verticesAreOnlyNextToNominalIP) )  continue;
      
      // if evet being actually processed, increment counter of analyzed events
      eventCounter++;
      
      speak_=false;
      if (entry<10 || entry%10000==0) speak_=true;
      
      if (speak_)  std::cout << "\n\n------> PRELIMINARY loop reading entry " << entry << "\tLS: " << treeVars_.lumiSection << " <------\n" ; 
      if (speak_)  std::cout << "  found " << treeVars_.nSuperClusters << " superclusters" << std::endl ;
      if (speak_)  std::cout << "  found " << treeVars_.nClusters << " basic clusters" << std::endl ;
      
      
      ///////////////////////////////////////////////////////////////////////
      // outer loop on supercluster
      for (int sc1=0; sc1<treeVars_.nSuperClusters; sc1++){
	
	float et1 = treeVars_.superClusterRawEnergy[sc1]/cosh( treeVars_.superClusterEta[sc1] );
	if (et1<20) continue;
	
	math::PtEtaPhiELorentzVectorD  el1(et1  ,
					   treeVars_.superClusterEta[sc1],
					   treeVars_.superClusterPhi[sc1],
					   treeVars_.superClusterRawEnergy[sc1] );
	
	///////////////////////////////////////////////////////////////////////
	// inner loop on supercluster
	for (int sc2=(sc1+1); sc2<treeVars_.nSuperClusters; sc2++){
	  
	  float et2 = treeVars_.superClusterRawEnergy[sc2]/cosh( treeVars_.superClusterEta[sc2] );
	  if (et2<20) continue;
	  
	  math::PtEtaPhiELorentzVectorD  el2(et2 ,
					     treeVars_.superClusterEta[sc2],
					     treeVars_.superClusterPhi[sc2],
					     treeVars_.superClusterRawEnergy[sc2] );
	  
	  // there seems to be a problem with vertexing - since nearly none of the electrons have the same vertex... CHECK!
	  float dvertex = pow(treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterVertexZ[sc2],2);
	  dvertex       = sqrt(dvertex);
	  
	  math::PtEtaPhiELorentzVectorD diEle = el1;
	  diEle += el2;
	  
	  // require invariant mass
	  if( fabs( diEle.M() - 91 ) > 20 ) continue;
	  // require two electrons from the same vertex
	  if ( dvertex > 0.01 )             continue;

	  // at this stage I have a suitable di-electron system for time studies
	  
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
	  
	
	  tmpEne=-9999;
	  // loop on BC and match to sc2 ==============
	  int bc2=-1;
	  for (int bc=0; bc<treeVars_.nClusters; bc++){
	    if ( pow(treeVars_.superClusterEta[sc2]-treeVars_.clusterEta[bc],2)+ pow(treeVars_.superClusterPhi[sc2]-treeVars_.clusterPhi[bc],2) < 0.02 
		 && treeVars_.clusterEnergy[bc]>tmpEne) {
	      tmpEne=treeVars_.clusterEnergy[bc];
	      bc2=bc;
	    }// end - if good bc candidate
	  }// end - loop over BC
	  
	  // protect in case of no matching
	  if(bc1==-1 || bc2==-1) continue;
	  if(0) {
	    std::cout << "\n\nsc1 : " << treeVars_.superClusterEta[sc1] << " " << treeVars_.superClusterPhi[sc1] << " " << treeVars_.superClusterRawEnergy[sc1] << std::endl;
	    std::cout << "bc1 : " << treeVars_.clusterEta[bc1] << " " << treeVars_.clusterPhi[bc1] << " " << treeVars_.clusterEnergy[bc1] << "\n"<< std::endl;
	    std::cout << "sc2 : " << treeVars_.superClusterEta[sc2] << " " << treeVars_.superClusterPhi[sc2] << " " << treeVars_.superClusterRawEnergy[sc2] << std::endl;
	    std::cout << "bc2 : " << treeVars_.clusterEta[bc2] << " " << treeVars_.clusterPhi[bc2] << " " << treeVars_.clusterEnergy[bc2] << std::endl;
	  }
	  

	  ClusterTime bcTime1 = timeAndUncertSingleCluster(bc1,treeVars_);
	  ClusterTime bcTime2 = timeAndUncertSingleCluster(bc2,treeVars_);
	  if(! (bcTime1.isvalid && bcTime2.isvalid) ) continue;
	

	  int run = treeVars_.runId;
	  
	  run_time = averageTimePerRun.find(  run  ) ;
	  // if this is a new run, add it to both maps 
	  if (run_time == averageTimePerRun.end() )
	    {
	      averageTimePerRun[ run ]   = std::pair<float,float>(0.,0.);
	      averageTimeSqPerRun[ run ] = std::pair<float,float>(0.,0.);
	      countersPerRun   [ run ]   = std::pair<int,int>(0,0);
	    }
	  run_time   = averageTimePerRun  .find( run ) ;
	  run_timeSq = averageTimeSqPerRun.find( run ) ;
	  run_count  = countersPerRun     .find( run ) ;
	  
	  
	  if      ( fabs(treeVars_.clusterEta[bc1])<1.4    &&  fabs(treeVars_.clusterEta[bc2])<1.4 ){
	    // if EBEB, and subcase 
	    run_count ->second.first += 2;
	    run_time  ->second.first += bcTime1.time;
	    run_time  ->second.first += bcTime2.time;
	    run_timeSq->second.first += (bcTime1.time)*(bcTime1.time);
	    run_timeSq->second.first += (bcTime2.time)*(bcTime2.time);
	  }
	  
	  else if ( fabs(treeVars_.clusterEta[bc1])>1.5    &&  fabs(treeVars_.clusterEta[bc2])>1.5 ){
	    // if EEEE, and subcase 
	    run_count ->second.second += 2;
	    run_time  ->second.second += bcTime1.time;
	    run_time  ->second.second += bcTime2.time;
	    run_timeSq->second.second += (bcTime1.time)*(bcTime1.time);
	    run_timeSq->second.second += (bcTime2.time)*(bcTime2.time);
	  }
	  
	  else if ( (fabs(treeVars_.clusterEta[bc1])<1.4 && fabs(treeVars_.clusterEta[bc2])>1.5) ||
		    (fabs(treeVars_.clusterEta[bc1])>1.5 && fabs(treeVars_.clusterEta[bc2])<1.4)    ) {
	    // if EBEE, and subcase 
	    run_count ->second.first  ++;
	    run_count ->second.second ++; 
	    if ( fabs(treeVars_.clusterEta[bc1])<1.4 ) {
	      run_time  ->second.first  += bcTime1.time; 
	      run_time  ->second.second += bcTime2.time; 
	      run_timeSq->second.first  += (bcTime1.time)*(bcTime1.time); 
	      run_timeSq->second.second += (bcTime2.time)*(bcTime2.time);    }
	    else {
	      run_time  ->second.first  += bcTime2.time;
	      run_time  ->second.second += bcTime1.time;
	      run_timeSq->second.first  += (bcTime2.time)*(bcTime2.time); 
	      run_timeSq->second.second += (bcTime1.time)*(bcTime1.time);    }
	  }
	  // if I've found a pair of supercluster, bail out of loop to repeat using twice the same supercluster
	  break;	
	  
	}// end loop sc2
      }// end loop sc1
      
    } // end of  auxiliary  PRELIMINARY loop over entries
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // check that the size of the three maps is the same
  if ( averageTimePerRun.size()==averageTimeSqPerRun.size() && averageTimeSqPerRun.size()==countersPerRun.size() && countersPerRun.size()<=1000 ){
    std::cout << "++ averageTimePerRun, averageTimeSqPerRun and countersPerRun have all the same size: " << countersPerRun.size() << std::endl;   
    std::cout << "++ value of stringDoAveragePerRun_ : " << stringDoAveragePerRun_ << std::endl;
  }
  else if (averageTimePerRun.size()>=1000 || averageTimeSqPerRun.size()>=1000 || countersPerRun.size()>=1000) {
    std::cout << "++ one of : averageTimePerRun, averageTimeSqPerRun and countersPerRun has size > 1000 (" << countersPerRun.size() << ") - bailing out" << std::endl;  assert(0);  }
  else{
    std::cout << "++ averageTimePerRun, averageTimeSqPerRun and countersPerRun DO NOT have all the same size - bailing out " << std::endl; assert(0);  }
  
  runCounter = 0;
  run_time   = averageTimePerRun.begin();
  run_timeSq = averageTimeSqPerRun.begin();
  run_count  = countersPerRun.begin();
  for( ;
       run_time!=averageTimePerRun.end() && run_timeSq!=averageTimeSqPerRun.end() && run_count!=countersPerRun.end() ;
         
       )
    {
      if( run_count->second.first < 2) continue;
      
      float avEb   = run_time  ->second.first  / run_count->second.first;
      float avEe   = run_time  ->second.second / run_count->second.second;
      float varEb  = ( run_timeSq->second.first   - avEb * run_time->second.first  ) / (run_count->second.first  - 1 );
      float varEe  = ( run_timeSq->second.second  - avEe * run_time->second.second ) / (run_count->second.second - 1 );
     
      std::cout << "\t\t run: " << run_count->first // << " (" << run_time->first << " , " << run_timeSq->first << ") "
		<< "\t\t entries EB: " << run_count->second.first 
		<< "\t av time EB: " << avEb
		<< "\t err time EB: " << sqrt(varEb)
		<< "\t\t entries EE: " << run_count->second.second
		<< "\t av. time EE: " << avEe
		<< "\t err time EB: " << sqrt(varEe)
		<< std::endl;
      
      runs[runCounter]      = run_count->first;
      runErrs[runCounter]   = 0.;
      ebTimes[runCounter]   = avEb;
      ebTimeErrs[runCounter]= sqrt(varEb/run_count->second.first);
      eeTimes[runCounter]   = avEe;
      eeTimeErrs[runCounter]= sqrt(varEe/run_count->second.second);

      run_count++;
      run_timeSq++;
      run_time++;
      runCounter++;
    }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  // save the history of times VS run number into a TGraph
  ebTimeVsRun = new TGraphErrors(runCounter, runs, ebTimes, runErrs, ebTimeErrs );   
  ebTimeVsRun = subDirGeneral.make <TGraphErrors> (runCounter, runs, ebTimes, runErrs, ebTimeErrs );   
  ebTimeVsRun->SetTitle("EB time VS run"); ebTimeVsRun->GetXaxis()->SetTitle("run number"); ebTimeVsRun->GetYaxis()->SetTitle("EB per-run average time [ns]");
  ebTimeVsRun->SetName("EB time VS run");   ebTimeVsRun->Write();

  eeTimeVsRun = new TGraphErrors (runCounter, runs, eeTimes, runErrs, eeTimeErrs );
  eeTimeVsRun->SetTitle("EE time VS run");  eeTimeVsRun->GetXaxis()->SetTitle("run number"); eeTimeVsRun->GetYaxis()->SetTitle("EE per-run average time [ns]");
  eeTimeVsRun->SetName("EE time VS run");   eeTimeVsRun->Write();  

  eeTimeVSebTime = new TGraphErrors (runCounter, ebTimes, eeTimes, ebTimeErrs, eeTimeErrs );
  eeTimeVSebTime->SetTitle("EE time VS EB time");  eeTimeVSebTime->GetXaxis()->SetTitle("EB per-run avegate time [ns]"); eeTimeVSebTime->GetYaxis()->SetTitle("EE per-run average time [ns]");
  eeTimeVSebTime->SetName("EE time VS EB time");    eeTimeVSebTime->Write();  

  timeCorrector theCorr;
  std::cout << "\ncreated object theCorr to be used for timeVsAmpliCorrections" << std::endl;
  std::cout << "\ninitializing theCorr" << std::endl;
  theCorr.initEB( std::string("EB") );
  theCorr.initEE( std::string("EElow") );



  // reset counter
  eventCounter = 0;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Main loop over entries
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (int entry = 0 ; (entry < nEntries && eventCounter < numEvents_); ++entry)
  {
    chain->GetEntry (entry) ;
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
    
    bool verticesAreOnlyNextToNominalIP;
    int  count=0;
    
    for(int v=0; v<treeVars_.nVertices; v++  )
	{        if (fabs(treeVars_.vtxZ[0])<15) count++; }
    
    if ( treeVars_.nVertices >0 && count==treeVars_.nVertices ) verticesAreOnlyNextToNominalIP = true;
    else                                                        verticesAreOnlyNextToNominalIP = false;
    
    //    --vertex: require vertex@IP (1), veto it (2) or either (0, or unset)
    if (flagOneVertex_ ==1 && (!verticesAreOnlyNextToNominalIP) ) continue;
    if (flagOneVertex_ ==2 && (verticesAreOnlyNextToNominalIP) )  continue;
    
    // if evet being actually processed, increment counter of analyzed events
    eventCounter++;

    // find the average time for EE and EB based on the auxiliary loop
    run_time   = averageTimePerRun  .find( treeVars_.runId );
    run_timeSq = averageTimeSqPerRun.find( treeVars_.runId );
    run_count  = countersPerRun     .find( treeVars_.runId );
    float eBtimeForRun=0.;
    float eEtimeForRun=0.;
    if( stringDoAveragePerRun_==1             && 
	run_time!=averageTimePerRun.end()     && 
	run_timeSq!=averageTimeSqPerRun.end() && 
	run_count!=countersPerRun.end() )
      {
	eBtimeForRun = (run_time->second.first  / run_count->second.first);
	eEtimeForRun = (run_time->second.second / run_count->second.second);
      }
    std::pair <float,float> thePhases ( eBtimeForRun, eEtimeForRun);


    speak_=false;
    if (entry<10 || entry%10000==0) speak_=true;

    if (speak_)  std::cout << "\n\n------> reading entry " << entry << "\tLS: " << treeVars_.lumiSection << " <------\n" ; 
    if (speak_)  std::cout << "  found " << treeVars_.nSuperClusters << " superclusters" << std::endl ;
    if (speak_)  std::cout << "  found " << treeVars_.nClusters << " basic clusters" << std::endl ;


    ///////////////////////////////////////////////////////////////////////
    // outer loop on supercluster
    for (int sc1=0; sc1<treeVars_.nSuperClusters; sc1++){

      float et1 = treeVars_.superClusterRawEnergy[sc1]/cosh( treeVars_.superClusterEta[sc1] );
      if (et1<20) continue;

      math::PtEtaPhiELorentzVectorD  el1(et1  ,
					 treeVars_.superClusterEta[sc1],
					 treeVars_.superClusterPhi[sc1],
					 treeVars_.superClusterRawEnergy[sc1] );

      ///////////////////////////////////////////////////////////////////////
      // inner loop on supercluster
      for (int sc2=(sc1+1); sc2<treeVars_.nSuperClusters; sc2++){

	float et2 = treeVars_.superClusterRawEnergy[sc2]/cosh( treeVars_.superClusterEta[sc2] );
	if (et2<20) continue;
	
	math::PtEtaPhiELorentzVectorD  el2(et2 ,
					   treeVars_.superClusterEta[sc2],
					   treeVars_.superClusterPhi[sc2],
					   treeVars_.superClusterRawEnergy[sc2] );
      
	// there seems to be a problem with vertexing - since nearly none of the electrons have the same vertex... CHECK!
	float dvertex = pow(treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterVertexZ[sc2],2);
	//dvertex       += pow(treeVars_.superClusterVertexY[sc1]-treeVars_.superClusterVertexY[sc2],2);
	//dvertex       += pow(treeVars_.superClusterVertexX[sc1]-treeVars_.superClusterVertexX[sc2],2);
	dvertex       = sqrt(dvertex);
	
	math::PtEtaPhiELorentzVectorD diEle = el1;
	diEle += el2;

	// ////////////////////////
	mass_      ->Fill(diEle.M());
	dZvertices_->Fill(dvertex);
	Zvertices_->Fill( (treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterVertexZ[sc2])/2 );
	nVertices_->Fill(treeVars_.nVertices);

	// require invariant mass
	if( fabs( diEle.M() - 91 ) > 20 ) continue;
	// require two electrons from the same vertex
	if ( dvertex > 0.01 )             continue;

	if(0) std::cout << "di-electron system mass: " << diEle.M() << " vertex distance: " << dvertex << std::endl;

	// at this stage I have a suitable di-electron system for time studies

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

	
	tmpEne=-9999;
	// loop on BC and match to sc2 ==============
	int bc2=-1;
	for (int bc=0; bc<treeVars_.nClusters; bc++){
	  if ( pow(treeVars_.superClusterEta[sc2]-treeVars_.clusterEta[bc],2)+ pow(treeVars_.superClusterPhi[sc2]-treeVars_.clusterPhi[bc],2) < 0.02 
	       && treeVars_.clusterEnergy[bc]>tmpEne) {
	    tmpEne=treeVars_.clusterEnergy[bc];
	    bc2=bc;
	  }// end - if good bc candidate
	}// end - loop over BC
	
	// protect in case of no matching
	if(bc1==-1 || bc2==-1) continue;
	if(0) {
	std::cout << "\n\nsc1 : " << treeVars_.superClusterEta[sc1] << " " << treeVars_.superClusterPhi[sc1] << " " << treeVars_.superClusterRawEnergy[sc1] << std::endl;
	std::cout << "bc1 : " << treeVars_.clusterEta[bc1] << " " << treeVars_.clusterPhi[bc1] << " " << treeVars_.clusterEnergy[bc1] << "\n"<< std::endl;
	std::cout << "sc2 : " << treeVars_.superClusterEta[sc2] << " " << treeVars_.superClusterPhi[sc2] << " " << treeVars_.superClusterRawEnergy[sc2] << std::endl;
	std::cout << "bc2 : " << treeVars_.clusterEta[bc2] << " " << treeVars_.clusterPhi[bc2] << " " << treeVars_.clusterEnergy[bc2] << std::endl;
	}
	
	ClusterTime bcTime1 = timeAndUncertSingleCluster(bc1,treeVars_);
	ClusterTime bcTime2 = timeAndUncertSingleCluster(bc2,treeVars_);

	if(! (bcTime1.isvalid && bcTime2.isvalid) ) continue;

	// fill the structures which hold all the plots
	plotsECALECAL.fill(sc1,sc2, bc1,bc2,thePhases);
	if      ( fabs(treeVars_.clusterEta[bc1])<1.4    &&  fabs(treeVars_.clusterEta[bc2])<1.4 ){
 	  plotsEBEB    .fill(sc1,sc2, bc1,bc2,thePhases);
	  
	  if           ( fabs(treeVars_.clusterEta[bc1])<0.44 && fabs(treeVars_.clusterEta[bc2])<0.44)  plotsEBEBmod1.fill(sc1,sc2, bc1,bc2,thePhases);
	  else if      ( fabs(treeVars_.clusterEta[bc1])<0.79 && fabs(treeVars_.clusterEta[bc2])<0.79)  plotsEBEBmod2.fill(sc1,sc2, bc1,bc2,thePhases);
	  else if      ( fabs(treeVars_.clusterEta[bc1])<1.13 && fabs(treeVars_.clusterEta[bc2])<1.13)  plotsEBEBmod3.fill(sc1,sc2, bc1,bc2,thePhases);
	  else if      ( fabs(treeVars_.clusterEta[bc1])<1.4  && fabs(treeVars_.clusterEta[bc2])<1.4 )  plotsEBEBmod4.fill(sc1,sc2, bc1,bc2,thePhases);

	  if      ( (treeVars_.clusterEta[bc1]>1. && treeVars_.clusterEta[bc2]>1. ) ||
		    (treeVars_.clusterEta[bc1]<-1. && treeVars_.clusterEta[bc2]<-1. ) )          plotsEBEBforward.fill(sc1,sc2, bc1,bc2,thePhases);

	  if( sqrt( pow(treeVars_.clusterEta[bc1]-treeVars_.clusterEta[bc2],2) +  pow(treeVars_.clusterPhi[bc1]-treeVars_.clusterPhi[bc2],2) ) < 1.5 ) {
	    plotsEBEBnear.fill(sc1,sc2, bc1,bc2,thePhases);}
	  else                                                                        plotsEBEBfar.fill(sc1,sc2, bc1,bc2,thePhases);

	  if( sqrt( pow(treeVars_.clusterEta[bc1]-treeVars_.clusterEta[bc2],2) +  pow(treeVars_.clusterPhi[bc1]-treeVars_.clusterPhi[bc2],2) ) > 3. ) {
	    plotsEBEBVeryFar.fill(sc1,sc2, bc1,bc2,thePhases);}

	  if (fabs ( treeVars_.clusterEta[bc1]-treeVars_.clusterEta[bc2] ) < 0.2 ){
	    plotsEBEBEtaNear.fill(sc1,sc2, bc1,bc2,thePhases);}


	  if       ( fabs( treeVars_.superClusterVertexZ[sc1] ) <  2.) {plotsEBEBzsmall.fill(sc1,sc2, bc1,bc2, thePhases);}
	  else if  (  treeVars_.superClusterVertexZ[sc1]        < -5.) {plotsEBEBzbigN .fill(sc1,sc2, bc1,bc2, thePhases);}
	  else if  (  treeVars_.superClusterVertexZ[sc1]        >  5.) {plotsEBEBzbigP .fill(sc1,sc2, bc1,bc2, thePhases);}

	  int type=3; float cut=6;
	  plotsEBEBchi2loose.fill(sc1,sc2, bc1,bc2, type, cut, thePhases); // cutting on chi2/ndf of 2
	  type=3; cut=4.;
	  plotsEBEBchi2    .fill(sc1,sc2, bc1,bc2, type, cut, thePhases); // cutting on chi2/ndf of 2.5
	  type=3; cut=3.;
	  plotsEBEBchi2tight.fill(sc1,sc2, bc1,bc2, type, cut, thePhases); // cutting on chi2/ndf of 2
	  type=1;  cut=1.5;
	  plotsEBEBseed2sec.fill(sc1,sc2, bc1,bc2, type, cut, thePhases); // cutting on agreement (<1.5ns) between seed and second, within a cluster

	  // fill different according to different pile up (selected on number of reco vertices)
	  if      (treeVars_.nVertices<6)  plotsEBEBlowPU  .fill(sc1,sc2, bc1,bc2, thePhases);
	  else if (treeVars_.nVertices<11) plotsEBEBmidPU  .fill(sc1,sc2, bc1,bc2, thePhases);
	  else if (treeVars_.nVertices<16) plotsEBEBhighPU .fill(sc1,sc2, bc1,bc2, thePhases);
	  else                             plotsEBEBsuperPU.fill(sc1,sc2, bc1,bc2, thePhases);

	  float energyRatio1 = treeVars_.xtalInBCEnergy[bc1][bcTime1.seed];
	  if(bcTime1.second>-1) {energyRatio1 /= treeVars_.xtalInBCEnergy[bc1][bcTime1.second]; }
	  else { energyRatio1 /= 99999; }
	  float energyRatio2 = treeVars_.xtalInBCEnergy[bc2][bcTime2.seed];
	  if(bcTime2.second>-1) {energyRatio2 /= treeVars_.xtalInBCEnergy[bc2][bcTime2.second]; }
	  else { energyRatio2 /= 99999; }
	  
	  float minRatio = 0.7; float maxRatio = 1.3;
	  if(minRatio<energyRatio1 && minRatio<energyRatio2 && energyRatio1<maxRatio && energyRatio2<maxRatio) 	  plotsEBEBequalShare.fill(sc1,sc2, bc1,bc2, thePhases);  
	  
	  minRatio = 2; maxRatio = 10;
	  if(minRatio<energyRatio1 && minRatio<energyRatio2 && energyRatio1<maxRatio && energyRatio2<maxRatio) 	  plotsEBEBunevenShare.fill(sc1,sc2, bc1,bc2, thePhases);  
	  
	}// if EBEB, and subcases
	else if ( fabs(treeVars_.clusterEta[bc1])>1.5    &&  fabs(treeVars_.clusterEta[bc2])>1.5 ) 	  {
	  plotsEEEE.fill(sc1,sc2, bc1,bc2, thePhases);
	  
	  if( sqrt( pow(treeVars_.clusterEta[bc1]-treeVars_.clusterEta[bc2],2) + pow(treeVars_.clusterPhi[bc1]-treeVars_.clusterPhi[bc2],2) ) < 0.75 ){
	    plotsEEEEnear.fill(sc1,sc2, bc1,bc2, thePhases);}
	  else                                                                      plotsEEEEfar .fill(sc1,sc2, bc1,bc2, thePhases);

	  if( sqrt( pow(treeVars_.clusterEta[bc1]-treeVars_.clusterEta[bc2],2) + pow(treeVars_.clusterPhi[bc1]-treeVars_.clusterPhi[bc2],2) ) > 1.5 ){
	    plotsEEEEVeryFar.fill(sc1,sc2, bc1,bc2, thePhases);}

	  if (fabs ( treeVars_.clusterEta[bc1]-treeVars_.clusterEta[bc2] ) < 0.2 ){
	    plotsEEEEEtaNear.fill(sc1,sc2, bc1,bc2,thePhases);}
      


	  if       ( fabs( treeVars_.superClusterVertexZ[sc1] ) <  2. ){plotsEEEEzsmall.fill(sc1,sc2, bc1,bc2, thePhases);}
	  else if  (  treeVars_.superClusterVertexZ[sc1]        < -5.) {plotsEEEEzbigN .fill(sc1,sc2, bc1,bc2, thePhases);}
	  else if  (  treeVars_.superClusterVertexZ[sc1]        >  5.) {plotsEEEEzbigP .fill(sc1,sc2, bc1,bc2, thePhases);}


	  if(  treeVars_.clusterEta[bc1]>1.5 && treeVars_.clusterEta[bc2]>1.5        ) 	  {
	    plotsEEPEEP.fill(sc1,sc2, bc1,bc2, thePhases);}
	  else if(  (treeVars_.clusterEta[bc1]>1.5 && treeVars_.clusterEta[bc2]<-1.5 ) ||  (treeVars_.clusterEta[bc2]>1.5 && treeVars_.clusterEta[bc1]<-1.5 ) ){
	    plotsEEPEEM.fill(sc1,sc2, bc1,bc2, thePhases); }
	  else if(  treeVars_.clusterEta[bc1]<-1.5    &&  treeVars_.clusterEta[bc2]<-1.5   )  {   
	    plotsEEMEEM.fill(sc1,sc2, bc1,bc2, thePhases); }
	}

	else if ( (fabs(treeVars_.clusterEta[bc1])<1.4 && fabs(treeVars_.clusterEta[bc2])>1.5) ||
		  (fabs(treeVars_.clusterEta[bc1])>1.5 && fabs(treeVars_.clusterEta[bc2])<1.4)    ) {
	  plotsEBEE.fill(sc1,sc2, bc1,bc2, thePhases);
	  
	  int type=3; float cut=6.; 
	  plotsEBEEchi2loose    .fill(sc1,sc2, bc1,bc2, type, cut, thePhases); // cutting on chi2/ndf of 2.
	  type=3; cut=4.;
	  plotsEBEEchi2    .fill(sc1,sc2, bc1,bc2, type, cut, thePhases); // cutting on chi2/ndf of 2.5
	  type=3; cut=3.;
	  plotsEBEEchi2tight    .fill(sc1,sc2, bc1,bc2, type, cut, thePhases); // cutting on chi2/ndf of 3.
	  type=1; cut=1.5;
	  plotsEBEEseed2sec.fill(sc1,sc2, bc1,bc2, type, cut, thePhases); // cutting on agreement (<1.5ns) between seed and second, within a cluster

	  if ( 
	      (treeVars_.clusterEta[bc1]<1.4 && treeVars_.clusterEta[bc1]>0 && treeVars_.clusterEta[bc2]>1.5 ) ||
	      (treeVars_.clusterEta[bc2]<1.4 && treeVars_.clusterEta[bc2]>0 && treeVars_.clusterEta[bc1]>1.5 ) 
	      ) plotsEBPEEP.fill(sc1,sc2, bc1,bc2, type, cut, thePhases);
	  
	  if ( 
	      (treeVars_.clusterEta[bc1]>-1.4 && treeVars_.clusterEta[bc1]<0 && treeVars_.clusterEta[bc2]<-1.5 ) ||
	      (treeVars_.clusterEta[bc2]>-1.4 && treeVars_.clusterEta[bc2]<0 && treeVars_.clusterEta[bc1]<-1.5 ) 
	      ) plotsEBMEEM.fill(sc1,sc2, bc1,bc2, type, cut, thePhases);
	  
	}
	// if I've found a pair of supercluster, bail out of loop to repeat using twice the same supercluster
	break;	
	
      }// end loop sc2
    }// end loop sc1
    
  }   // end of loop over entries
  

  delete chain ;
  
  return 0 ;
}



/*
	// computing extra time of flight due to bending of electrons
	float r_curv1 = et1 / 0.3 / 3.8 ; // in meters
	float distShowerToVertex =        // IP - shower distance in transverse plane (~ ECAL Radius)
	pow( (treeVars_.superClusterVertexX[sc1]-treeVars_.superClusterX[sc1]), 2) +
	pow( (treeVars_.superClusterVertexY[sc1]-treeVars_.superClusterY[sc1]), 2);  
	// pow( (treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterZ[sc1]), 2);
	distShowerToVertex = sqrt(distShowerToVertex) / 100;  // in meters
	//std::cout << "distShowerToVertex: " << distShowerToVertex << std::endl;
	float alpha = asin( 0.5 * distShowerToVertex / r_curv1 ) ;
	float road  = r_curv1 * 2 * alpha;
	road        = sqrt( road*road + pow(  (treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterZ[sc1])  ,2)/100./100. );
	// actual time of flight of electron subtracted of straight line 
	float roadMoreThanPhoton = road - 
	sqrt( pow( (treeVars_.superClusterVertexX[sc1]-treeVars_.superClusterX[sc1]), 2) +
	pow( (treeVars_.superClusterVertexY[sc1]-treeVars_.superClusterY[sc1]), 2) +
	pow( (treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterZ[sc1]), 2)
	)/100.;
	std::cout << "eta: " << treeVars_.superClusterEta[sc1] << " pt: " << treeVars_.superClusterRawEnergy[sc1] 
	<< " r_curv1: " << r_curv1 << " transv-distShowerToVertex: " << distShowerToVertex 
	<< " path: " << road << " roadMoreThanPhoton: " << roadMoreThanPhoton << std::endl;
*/
