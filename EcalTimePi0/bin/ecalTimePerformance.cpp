#include <iostream>
#include <string>

#include <map>
#include <vector>
#include <algorithm> 
#include <functional>
#include <set>
#include <boost/tokenizer.hpp>

#include "CalibCalorimetry/EcalTiming/interface/EcalTimeTreeContent.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"


typedef std::set<std::pair<int,int> > SetOfIntPairs;
float travelDistance(int sc_num, EcalTimeTreeContent treeVars_);
float extraTravelTime(int sc_num, EcalTimeTreeContent & treeVars_);

// authors: S. Cooper and G. Franzoni (UMN)

#define BarrelLimit  1.479
#define EndcapLimit  3.0

#define ADCtoGeVEB   0.039
#define ADCtoGeVEE   0.063

#define numAeffBins     35
#define numAoSigmaBins  25

#define NSlices 54 // # of slices in log(A)
#define LogStep 0.05 // size of Slices log(A)

#define lightSpeed 299792458

struct ClusterTime {
  int   numCry;
  float time;
  float timeErr;
  float chi2;
} ;


// -------- Globals ----------------------------------------
EcalTimeTreeContent treeVars_; 
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


float minAmpliOverSigma_   = 30;    // dimensionless

float maxChi2NDF_ = 20;  //TODO: gf configurable

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

// Consts
//const float sigmaNoiseEB        = 0.75;  // ADC ; using high frequency noise
//const float sigmaNoiseEE        = 1.58;  // ADC ; using high frequency noise
const float sigmaNoiseEB          = 1.06;  // ADC ; using total single-sample noise
const float sigmaNoiseEE          = 2.10;  // ADC ; using total single-sample noise
// const float timingResParamN       = 35.1; // ns ; Fig. 2 from CFT-09-006
// const float timingResParamConst   = 0.020; //ns ;   "
const float timingResParamNEB     = 28.51;   // ns ; plots approved http://indico.cern.ch/conferenceDisplay.py?confId=92739
const float timingResParamConstEB = 0.35;     // ns ; rough, probably conservative estimate
const float timingResParamNEE     = 31.84;   // ns ; Fig. 2 from CFT-09-006
const float timingResParamConstEE = 0.5;     // ns ; rough, probably conservative estimate

// -------- Histograms -------------------------------------
TH1F* xtalEnergyHist_;
TH1F* mass_;
TH1F* dZvertices_;
TH1F* seedTimeDiffHist_, *seedTimeDiffHistEE_, *seedTimeDiffHistEB_;
TH1F* clusTimeDiffHist_, *clusTimeDiffHistEE_, *clusTimeDiffHistEB_;
TH1F* clusTimeDiffHistTOF_, *clusTimeDiffHistEETOF_, *clusTimeDiffHistEBTOF_;
TH1F* numCryBC1, *numCryBC2;

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
      minAmpliOverSigma_  = atof(argv[v+1]);
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
void initializeHists(){
  saving_->cd();
  // Initialize histograms -- xtals
  xtalEnergyHist_   = new TH1F("XtalEnergy","Crystal energy;GeV",110,-1,10);
  mass_         = new TH1F("mass","mass; GeV",80,50,130);
  dZvertices_   = new TH1F("dZvertices","dZvertices; #DeltaZ(ele_{1},ele_{2}) [cms]",250,0,25);
  seedTimeDiffHist_ = new TH1F("seed time difference","seed time difference; #Deltat [ns]",1000,-5,5);
  seedTimeDiffHistEE_ = new TH1F("EE seed time difference","EE seed time difference; #Deltat [ns]",1000,-5,5);
  seedTimeDiffHistEB_ = new TH1F("EB seed time difference","EB seed time difference; #Deltat [ns]",1000,-5,5);
  clusTimeDiffHist_ = new TH1F("cluster time difference","cluster time difference; #Deltat [ns]",1000,-5,5);
  clusTimeDiffHistEE_ = new TH1F("EE cluster time difference","EE cluster time difference; #Deltat [ns]",1000,-5,5);
  clusTimeDiffHistEB_ = new TH1F("EB cluster time difference","EB cluster time difference; #Deltat [ns]",1000,-5,5);
  clusTimeDiffHistTOF_ = new TH1F("TOF-corr cluster time difference","TOF-corr cluster time difference; #Deltat [ns]",1000,-5,5);
  clusTimeDiffHistEETOF_ = new TH1F("TOF-corr EE cluster time difference","TOF-corr EE cluster time difference; #Deltat [ns]",1000,-5,5);
  clusTimeDiffHistEBTOF_ = new TH1F("TOF-corr EB cluster time difference","TOF-corr EB cluster time difference; #Deltat [ns]",1000,-5,5);
  numCryBC1 = new TH1F("num cry in bc1","num cry in bc1; num cry",25,0,25);
  numCryBC2 = new TH1F("num cry in bc2","num cry in bc2; num cry",25,0,25);
}//end initializeHists

// ---------------------------------------------------------------------------------------
// ------------------ Function to compute time and error for a cluster -------------------
ClusterTime timeAndUncertSingleCluster(int bClusterIndex)
{
  ClusterTime theResult; //initialize
  theResult.numCry = -999999;   theResult.time   = -999999;
  theResult.timeErr= -999999;   theResult.chi2   = -999999;

  float weightTsum  = 0;
  float weightSum   = 0;
  int   numCrystals = 0;
  float timingResParamN    =0;
  float timingResParamConst=0;
  
  //std::cout << "\n++ BC statrs (eta: " << treeVars_.clusterEta[bClusterIndex] << ") : "  << std::endl;
  // loop on the cry components of a basic cluster; get timeBest and uncertainty 
  for(int thisCry=0; thisCry<treeVars_.nXtalsInCluster[bClusterIndex]; thisCry++)
  {
    bool  thisIsInEB=false;
    float sigmaNoiseOfThis=0;
    if(treeVars_.xtalInBCIEta[bClusterIndex][thisCry]!=-999999)       {
      sigmaNoiseOfThis   =sigmaNoiseEB;
      timingResParamN    =timingResParamNEB;
      timingResParamConst=timingResParamConstEB;
      thisIsInEB=true;    }
    else if(treeVars_.xtalInBCIy[bClusterIndex][thisCry]!=-999999)    {
      sigmaNoiseOfThis=sigmaNoiseEE;
      timingResParamN    =timingResParamNEE;
      timingResParamConst=timingResParamConstEE;
      thisIsInEB=false;    }
    else    {  std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}

    // remove hits beyond gain switch
    if ( treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] > 3950 )  continue;

    float ampliOverSigOfThis = treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] / sigmaNoiseOfThis; 
    // minimum amplitude and spike rejection (could be updated)
    if( ampliOverSigOfThis < minAmpliOverSigma_) continue;
    if( treeVars_.xtalInBCSwissCross[bClusterIndex][thisCry] > 0.95) continue;

    numCrystals++;
    float timeOfThis  = treeVars_.xtalInBCTime[bClusterIndex][thisCry];
    //    old estimated: fully parameterized
    //    float sigmaOfThis = sqrt(pow(timingResParamN/ampliOverSigOfThis,2)+pow(timingResParamConst,2));
    //    new estimate: error from ratio + constant term  
    // float sigmaOfThis = pow( treeVars_.xtalInBCTimeErr[bClusterIndex][thisCry], 2 ) + pow( timingResParamConst, 2);
    // supposedly a large time constant is already included in the cxtalInBCTimeErr of 600 ps; rough estimate is that it's actually 300 (EB)

    // there appears to be a sizable fraction of hits which have xtalInBCTimeErr=0! hm?
    // float sigmaOfThis = pow( treeVars_.xtalInBCTimeErr[bClusterIndex][thisCry], 2); // - 0.6*0.6 + 0.3*0.3 ;
    // temporarily use the ROUGH parameterization 
    float sigmaOfThis = pow( timingResParamN / ampliOverSigOfThis , 2) + timingResParamConst*timingResParamConst;
    sigmaOfThis       = sqrt(sigmaOfThis);

    //std::cout << "t: " << timeOfThis << " a/s: " << ampliOverSigOfThis << " sig: " << sigmaOfThis << "\t\t";
    weightTsum+=(timeOfThis/pow(sigmaOfThis,2));
    weightSum+=1/pow(sigmaOfThis,2);
  }


  float bestTime = weightTsum/weightSum;
  //  if (weightTsum!=0 ) std::cout << "bestTime = " << bestTime << std::endl;
  //  else               {
  //    std::cout << "bestTime n.a. " << std::endl;
  //    return theResult;
  //  }
  

  float chi2 = -999999;
  // loop on the cry components to get chi2
  // do this only if you have at least 2 crystals over threshold and not spiky
  if(numCrystals>1){
    chi2=0;
    for(int thisCry=0; thisCry<treeVars_.nXtalsInCluster[bClusterIndex]; thisCry++)
      {
  	//bool  thisIsInEB=false;
  	float sigmaNoiseOfThis=0;
  	if(treeVars_.xtalInBCIEta[bClusterIndex][thisCry]!=-999999)       {
  	  sigmaNoiseOfThis=sigmaNoiseEB;
  	  //thisIsInEB=true;
  	}
  	else if(treeVars_.xtalInBCIy[bClusterIndex][thisCry]!=-999999)    {
  	  sigmaNoiseOfThis=sigmaNoiseEE;
  	  //thisIsInEB=false;    
  	}
  	else    {  std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}
  	
  	float ampliOverSigOfThis = treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] / sigmaNoiseOfThis; 
  	if( ampliOverSigOfThis < minAmpliOverSigma_) continue;
  	
  	float timeOfThis  = treeVars_.xtalInBCTime[bClusterIndex][thisCry];
  	float sigmaOfThis = sqrt(pow(timingResParamN/ampliOverSigOfThis,2)+pow(timingResParamConst,2));
  	
  	chi2 += pow( (timeOfThis-bestTime)/sigmaOfThis, 2);
  	
      }// end loop on cry
  }//end if

  
  if(weightSum <= 0) {
    return theResult;}
  else{
    //std::cout << "-- GFdeb time: " << bestTime << " error: " << sqrt(1/weightSum) << std::endl;//gfdebug
    theResult.numCry = numCrystals;
    theResult.time   = bestTime;
    theResult.timeErr= sqrt(1/weightSum);
    theResult.chi2   = chi2;
    return theResult;
  }

}// end timeAndUncertSingleCluster


// ---------------------------------------------------------------------------------------
// ------------------ Function to write hists --------------------------------------------
void writeHists()
{
  saving_->cd();
  // write out control histograms
  TDirectory *controlPlots = saving_->mkdir("ecalTimePerformance-plots");
  controlPlots->cd();
  xtalEnergyHist_   ->Write(); 
  mass_   ->Write(); 
  dZvertices_->Write(); 
  seedTimeDiffHist_ ->Write();   seedTimeDiffHistEE_ ->Write();   seedTimeDiffHistEB_ ->Write(); 
  clusTimeDiffHist_ ->Write();   clusTimeDiffHistEE_ ->Write();   clusTimeDiffHistEB_ ->Write(); 
  clusTimeDiffHistTOF_ ->Write();   clusTimeDiffHistEETOF_ ->Write();   clusTimeDiffHistEBTOF_ ->Write(); 
  numCryBC1 ->Write();
  numCryBC2 ->Write();
}



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

  // Initialize output root file
  saving_ = new TFile(outputRootName_.c_str (),"recreate");

  // Initialize the histograms
  initializeHists();

  int eventCounter = 0;
  /////////////////////////////////////////////////////
  // Main loop over entries
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

	// require invariant mass
	if( fabs( diEle.M() - 91 ) > 20 ) continue;
	// require two electrons from the same vertex
	if ( dvertex > 0.01 )             continue;

	if(0) std::cout << "di-electron system mass: " << diEle.M() << " vertex distance: " << dvertex << std::endl;
	
	// at this stage I have a suitable di-electron system for time studies
	
	  
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
	
	int seed1(-1); tmpEne=-9999;
	for (int cry=0; cry<treeVars_.nXtalsInCluster[bc1]; cry++){
	  if(treeVars_.xtalInBCEnergy[bc1][cry]>tmpEne){
	    tmpEne=treeVars_.xtalInBCEnergy[bc1][cry];
	    seed1=cry;
	  } 	}
	float timeSeed1=treeVars_.xtalInBCTime[bc1][seed1];

	int seed2(-1); tmpEne=-9999;
	for (int cry=0; cry<treeVars_.nXtalsInCluster[bc2]; cry++){
	  if(treeVars_.xtalInBCEnergy[bc2][cry]>tmpEne){
	    tmpEne=treeVars_.xtalInBCEnergy[bc2][cry];
	    seed2=cry;
	  } 	}
	float timeSeed2=treeVars_.xtalInBCTime[bc2][seed2];
	
	ClusterTime bcTime1 = timeAndUncertSingleCluster(bc1);
	ClusterTime bcTime2 = timeAndUncertSingleCluster(bc2);
	
	seedTimeDiffHist_ -> Fill(timeSeed1 - timeSeed2 );
	
	clusTimeDiffHist_    -> Fill(bcTime1.time - bcTime2.time );
	clusTimeDiffHistTOF_ -> Fill( (bcTime1.time - extraTravelTime(sc1,treeVars_) ) - (bcTime2.time -extraTravelTime(sc2,treeVars_)) );
	
	if ( fabs(treeVars_.clusterEta[bc1]) < 1.4 &&  fabs(treeVars_.clusterEta[bc2]) < 1.4){
	  clusTimeDiffHistEB_ -> Fill(bcTime1.time - bcTime2.time );
	  seedTimeDiffHistEB_ -> Fill(timeSeed1 - timeSeed2 );
	  clusTimeDiffHistEBTOF_ -> Fill( (bcTime1.time - extraTravelTime(sc1,treeVars_) ) - (bcTime2.time -extraTravelTime(sc2,treeVars_)) );
	}
	else{
 	clusTimeDiffHistEE_ -> Fill(bcTime1.time - bcTime2.time );
	seedTimeDiffHistEE_ -> Fill(timeSeed1 - timeSeed2 );
	clusTimeDiffHistEETOF_ -> Fill( (bcTime1.time - extraTravelTime(sc1,treeVars_) ) - (bcTime2.time -extraTravelTime(sc2,treeVars_)) );
	}

	numCryBC1->Fill(bcTime1.numCry);
	numCryBC2->Fill(bcTime2.numCry);
	

	
	// if I've found a pair of supercluster, bail out of loop to repeat using twice the same supercluster
	break;	

      }// end loop sc2
    }// end loop sc1

  }   // end of loop over entries


  // now save the plots
  writeHists();
  saving_->Close();

  delete chain ;
  
  return 0 ;
}


float travelDistance(int sc_num, EcalTimeTreeContent treeVars_) {
  return
    sqrt (	  pow( (treeVars_.superClusterVertexX[sc_num]-treeVars_.superClusterX[sc_num]), 2) +
		  pow( (treeVars_.superClusterVertexY[sc_num]-treeVars_.superClusterY[sc_num]), 2) +   
		  pow( (treeVars_.superClusterVertexZ[sc_num]-treeVars_.superClusterZ[sc_num]), 2)
		  );
}

float extraTravelTime(int sc_num, EcalTimeTreeContent & treeVars_) { // extra travel time with respect to collision at IP, in ns
  
  float travelled = sqrt (	  pow( (treeVars_.superClusterX[sc_num]-treeVars_.superClusterVertexX[sc_num]), 2) +
				  pow( (treeVars_.superClusterY[sc_num]-treeVars_.superClusterVertexY[sc_num]), 2) +   
				  pow( (treeVars_.superClusterZ[sc_num]-treeVars_.superClusterVertexZ[sc_num]), 2)
				  );
  float nominal = sqrt (	  pow( (treeVars_.superClusterX[sc_num]), 2) +
				  pow( (treeVars_.superClusterY[sc_num]), 2) +   
				  pow( (treeVars_.superClusterZ[sc_num]), 2)
				  );

  //std::cout << "extraTravelTime [ns]: " <<  (travelled-nominal)/100./lightSpeed*1e9 << std::endl;
  return  (travelled-nominal)/100./lightSpeed*1e9;

}

