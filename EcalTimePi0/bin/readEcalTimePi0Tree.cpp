#include <iostream>
#include <string>

#include <map>
#include <vector>
#include <algorithm> 
#include <functional>
#include <set>

#include "ECALTime/EcalTimePi0/interface/EcalTimePi0TreeContent.h"
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

// initial authors P. Govoni et al
// authors: S. Cooper and G. Franzoni (UMN)

#define BarrelLimit  1.479
#define EndcapLimit  3.0

#define ADCtoGeVEB   0.039
#define ADCtoGeVEE   0.063

#define numAeffBins  40

struct ClusterTime {
  int   numCry;
  float time;
  float timeErr;
  float chi2;
} ;


// -------- Globals ----------------------------------------
EcalTimePi0TreeContent treeVars_; 
TFile* saving_;
std::vector<std::string> listOfFiles_;
bool speak_=false;
char buffer_ [75];
std::string bufferTitle_; 
// mass spectraVars
float pi0MassEB_=0;
float pi0WidthEB_=0;
float pi0MassEE_=0;
float pi0WidthEE_=0;
float pi0MassEEP_=0;
float pi0WidthEEP_=0;
float pi0MassEEM_=0;
float pi0WidthEEM_=0;
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

float minAmpliOverSigma_   = 7;    // dimensionless

float maxChi2NDF_ = 20;  //TODO: gf configurable

int  minEntriesForFit_ = 7;
bool limitFit_(true); 
//std::string fitOption_(""); // default: use chi2 method
std::string fitOption_("L"); // use likelihood method



//parameters for histograms binning and ranges ///////////////////////////////////////////

// double AeffBins_[41] = {0,                          // set of fine bins for large stats
//  			1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,
//  			11,12,13,14,15,16,17,18,19,20,
//  			22,24,26,28,30,32,34,36,38,40,
//  			44,48,52,56,60,68,76,84,110,150};
// int    AeffNBins_    = 40;
// int    AeffMax_      = 150;

double AeffBins_[35] = {0,                          // set of fine bins for large stats
 			1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,
 			11,12,13,14,15,16,17,18,19,20,
 			22,24,26,28,30,
                        32,36,40,46,52,
                        58,74,86,120};
int    AeffNBins_    = 34;
int    AeffMax_      = 120;

// defining arrays large enough; number of actual bins to define TGraphErrors comes from AeffNBins_ 
float  AeffBinCentersAny_[256]; float  AeffBinCentersErrAny_[256]; float  sigmaAeffAny_[256]; float  sigmaAeffErrAny_[256];
float  AeffBinCentersEB_[256];  float  AeffBinCentersErrEB_[256];  float  sigmaAeffEB_[256];  float  sigmaAeffErrEB_[256];
float  AeffBinCentersEE_[256];  float  AeffBinCentersErrEE_[256];  float  sigmaAeffEE_[256];  float  sigmaAeffErrEE_[256];

// double AeffBins_[21] = {0,                          // set of coarser bins for smaller stats
//  			2 ,4 ,6 ,8 ,10,
//  			12,14,16,18,20,
//  			24,28,32,36,40,
//  			48,56,68,84,150};
// int    AeffNBins_    = 20;
// int    AeffMax_      = 150;



int numDtBins_  = 75;
int DtMax_      = 15; // useful to catch tails also at low Aeff (<10)



// Consts
//const float sigmaNoiseEB        = 0.75;  // ADC ; using high frequency noise
//const float sigmaNoiseEE        = 1.58;  // ADC ; using high frequency noise
const float sigmaNoiseEB        = 1.06;  // ADC ; using total single-sample noise
const float sigmaNoiseEE        = 2.10;  // ADC ; using total single-sample noise
const float timingResParamN     = 35.1; // ns ; Fig. 2 from CFT-09-006
const float timingResParamConst = 0.020; //ns ;   "
//const float timingResParamN     = 31.5; // ns ; Fig. 5 from CFT-09-006
//const float timingResParamConst = 0.38; //ns ;   "
// EB noise from Jean: https://espace.cern.ch/cmsccecal/ECAL%20PFG%20and%20offline%20weekly/default.aspx?InstanceID=35&Paged=Next&p_StartTimeUTC=20090603T140000Z&View={38FE356C-17A7-4C7D-987B-8302CABFAD4F}
// EE noise from Jean: https://espace.cern.ch/cmsccecal/ECAL%20PFG%20and%20offline%20weekly/default.aspx?InstanceID=44&Paged=Next&p_StartTimeUTC=20090603T140000Z&View={38FE356C-17A7-4C7D-987B-8302CABFAD4F}
// -------- Histograms -------------------------------------
// xtals
TH1F* xtalEnergyHist_;
TH1F* xtalTimeHist_;
TH1F* xtalIEtaHist_;
TH1F* xtalIPhiHist_;
TH1F* xtalIXHist_;
TH1F* xtalIYHist_;
// TH1F* xtalStatusHist_;
// TH2F* xtalOccupancyHistEB_;

// BasicClusters
TH1F* BCNumPerEventHist_;
TH1F* BCNumCrysHist_;
TH1F* BCNumCrysOverThrHist_;
TH1F* BCNumCrysOverThrEBHist_;
TH1F* BCNumCrysOverThrEEHist_;
TH1F* BCEnergyHist_;
TH1F* BCEtHist_;
TH2F* BCOccupancyEBHist_;
TH2F* BCOccupancyEEHist_;
TH2F* BCOccupancyEEPHist_;
TH2F* BCOccupancyEEMHist_;
TH2F* BCOccupancyHistAny_;
TH1F* BCEtaHist_;
TH1F* BCPhiHist_;
TH1F* BCClusterShapeEEHist_;
TH1F* BCClusterShapeEEPHist_;
TH1F* BCClusterShapeEEMHist_;
TH1F* BCClusterShapeEBHist_;
// diphotons control plots
TH1F* massDiGammaHist_;
TH1F* massDiGammaEBHist_;
TH1F* massDiGammaEEHist_;
TH1F* massDiGammaEEPHist_;
TH1F* massDiGammaEEMHist_;
TH2F* diPhotonOccupancyAny_;
// single cluster resolution
TH1F*   dtUpToQuarterGeVEB_;
TH1F*   dtUpToHalfGeVEB_;
TH1F*   dtUpToOneGeVEB_;
TH1F*   dtUpToTwoGeVEB_;
TH1F*   dtUpOverTwoGeVEB_;

TH1F*   dtUpToThreeQuarterGeVEE_;
TH1F*   dtUpToOneAndHalfGeVEE_;
TH1F*   dtUpToThreeGeVEE_;
TH1F*   dtUpToSixGeVEE_;
TH1F*   dtUpOverSixGeVEE_;

TH2F*   dtVSAeffHistAny_;
TH1F*   dtSliceVSAeffAny_[numAeffBins];
TH1F*   AeffSliceAny_[numAeffBins];
TH2F*   dtVSAeffHistEB_;
TH1F*   dtSliceVSAeffEB_[numAeffBins];
TH1F*   AeffSliceEB_[numAeffBins];
TH2F*   dtVSAeffHistEE_;
TH1F*   dtSliceVSAeffEE_[numAeffBins];
TH1F*   AeffSliceEE_[numAeffBins];
TProfile* dtVSAeffProfAny_;
TProfile* dtVSAeffProfEB_;
TProfile* dtVSAeffProfEE_;

TH1F*  dtRMSVSAeffAny_; TH1F*  dtSigmaAeffAny_;
TH1F*  dtRMSVSAeffEB_;  TH1F*  dtSigmaAeffEB_;
TH1F*  dtRMSVSAeffEE_;  TH1F*  dtSigmaAeffEE_;
TH1F*  singleClusterChi2HistAny_;
TH1F*  singleClusterChi2NDFHistAny_;
TH1F*  singleClusterChi2HistEB_;
TH1F*  singleClusterChi2NDFHistEB_;
TH1F*  singleClusterChi2HistEE_;
TH1F*  singleClusterChi2NDFHistEE_;

TH2F* deltaTCrysVsAmplitudeEB_;
TH2F* deltaTCrysVsAmplitudeEE_;
TH2F* deltaTCrysVsAmplitudeAny_;

// selection on pi0 candidates 
TH2F* diPhotonPeakOccupancyAny_;
TH2F* diPhotonSidesOccupancyAny_;
// single cluster resolution in pi0 peak
TH1F*     dtUpToQuarterGeVEBPeak_;
TH1F*     dtUpToHalfGeVEBPeak_;
TH1F*     dtUpToOneGeVEBPeak_;
TH1F*     dtUpToTwoGeVEBPeak_;
TH1F*     dtUpOverTwoGeVEBPeak_;

TH1F*     dtUpToThreeQuarterGeVEEPeak_;
TH1F*     dtUpToOneAndHalfGeVEEPeak_;
TH1F*     dtUpToThreeGeVEEPeak_;
TH1F*     dtUpToSixGeVEEPeak_;
TH1F*     dtUpOverSixGeVEEPeak_;

TH2F*     dtVSAeffHistAnyPeak_;
TH1F*     dtSliceVSAeffAnyPeak_[numAeffBins];
TH2F*     dtVSAeffHistEBPeak_;
TH1F*     dtSliceVSAeffEBPeak_[numAeffBins];
TH2F*     dtVSAeffHistEEPeak_;
TH1F*     dtSliceVSAeffEEPeak_[numAeffBins];
TH1F*     dtRMSVSAeffAnyPeak_; TH1F*    dtSigmaAeffAnyPeak_;
TH1F*     dtRMSVSAeffEBPeak_;  TH1F*    dtSigmaAeffEBPeak_;
TH1F*     dtRMSVSAeffEEPeak_;  TH1F*    dtSigmaAeffEEPeak_;
TProfile* dtVSAeffProfAnyPeak_;
TProfile* dtVSAeffProfEBPeak_;
TProfile* dtVSAeffProfEEPeak_;
TH1F*  singleClusterChi2HistAnyPeak_;
TH1F*  singleClusterChi2NDFHistAnyPeak_;
TH1F*  singleClusterChi2HistEBPeak_;
TH1F*  singleClusterChi2NDFHistEBPeak_;
TH1F*  singleClusterChi2HistEEPeak_;
TH1F*  singleClusterChi2NDFHistEEPeak_;
TH1F*  dtPullSingleClusterHistPi0Peak_;
TH1F*  dtPullSingleClusterHistPi0PeakEE_;
TH1F*  dtPullSingleClusterHistPi0PeakEB_;

// double cluster resolution
TH1F* dtDoubleClusterHistAny_;
TH1F* dtPullDoubleClusterHistAny_;
TH1F* dtPullChi2CutDoubleClusterHistAny_;
TH2F* dtVsPtDoubleClusterHistAny_;

TH1F* dtDoubleClusterHistPi0Peak_;
TH1F* dtDoubleClusterHistPi0PeakEE_;
TH1F* dtDoubleClusterHistPi0PeakEB_;
TH1F* dtPullDoubleClusterHistPi0Peak_;
TH1F* dtPullDoubleClusterHistPi0PeakEE_;
TH1F* dtPullDoubleClusterHistPi0PeakEB_;
TH1F* dtPullChi2CutDoubleClusterHistPi0Peak_;
TH1F* dtPullChi2CutDoubleClusterHistPi0PeakEE_;
TH1F* dtPullChi2CutDoubleClusterHistPi0PeakEB_;
TH2F* dtVsPtDoubleClusterHistPi0Peak_;
TH2F* dtVsPtDoubleClusterHistPi0PeakEE_;
TH2F* dtVsPtDoubleClusterHistPi0PeakEB_;

// double cluster studies VS Aeff
TH2F* dtSeedsVsAeffDoubleClusterHist_;             TH1F*    dtSeedsSigmaAeffAny_;
TH2F* dtSeedsVsAeffDoubleClusterHistEE_;           TH1F*    dtSeedsSigmaAeffEE_;
TH2F* dtSeedsVsAeffDoubleClusterHistEB_;           TH1F*    dtSeedsSigmaAeffEB_;
TH2F* dtSeedsVsAeffDoubleClusterHistPi0Peak_;      TH1F*    dtSeedsSigmaAeffPi0Peak_;
TH2F* dtSeedsVsAeffDoubleClusterHistPi0PeakEE_;    TH1F*    dtSeedsSigmaAeffPi0PeakEE_;
TH2F* dtSeedsVsAeffDoubleClusterHistPi0PeakEB_;    TH1F*    dtSeedsSigmaAeffPi0PeakEB_;
TH2F* dtVsAeffDoubleClusterHist_;                  TH1F*    dtClustersSigmaAeffAny_;
TH2F* dtVsAeffDoubleClusterHistEE_;                TH1F*    dtClustersSigmaAeffEE_;
TH2F* dtVsAeffDoubleClusterHistEB_;                TH1F*    dtClustersSigmaAeffEB_;
TH2F* dtVsAeffDoubleClusterHistPi0Peak_;           TH1F*    dtClustersSigmaAeffPi0Peak_;
TH2F* dtVsAeffDoubleClusterHistPi0PeakEE_;         TH1F*    dtClustersSigmaAeffPi0PeakEE_;
TH2F* dtVsAeffDoubleClusterHistPi0PeakEB_;         TH1F*    dtClustersSigmaAeffPi0PeakEB_;

TH1F* dtSeedsSlicesVsAeffDoubleClusterHist_[numAeffBins];         
TH1F* dtSeedsSlicesVsAeffDoubleClusterHistEE_[numAeffBins];       
TH1F* dtSeedsSlicesVsAeffDoubleClusterHistEB_[numAeffBins];       
TH1F* dtSeedsSlicesVsAeffDoubleClusterHistPi0Peak_[numAeffBins];  
TH1F* dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEE_[numAeffBins];
TH1F* dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEB_[numAeffBins];
TH1F* dtSlicesVsAeffDoubleClusterHist_[numAeffBins];              
TH1F* dtSlicesVsAeffDoubleClusterHistEE_[numAeffBins];            
TH1F* dtSlicesVsAeffDoubleClusterHistEB_[numAeffBins];            
TH1F* dtSlicesVsAeffDoubleClusterHistPi0Peak_[numAeffBins];       
TH1F* dtSlicesVsAeffDoubleClusterHistPi0PeakEE_[numAeffBins];     
TH1F* dtSlicesVsAeffDoubleClusterHistPi0PeakEB_[numAeffBins];     



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


  //gf: support development
  //std::cout << "\nargc:       " << argc << std::endl;
  //for (int v=0; v<argc; v++ ){      std::cout << "argument: " << v << " argv: " << argv[v] << std::endl;    }

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
      std::cout << " --eTPi0MinEE min eT for EE pi0 candidate" << std::endl;
      std::cout << " --minAOverSigma min ampli considered for time" << std::endl;
      std::cout << " --minRun lowest run number considered" << std::endl;
      std::cout << " --maxRun highest run number considered" << std::endl;
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
//initializeHists
  int numChi2Bins = 400;
  int chi2Max = 200;
  int numChi2NDFBins = 200;
  int chi2NDFMax = 100;

  saving_->cd();
  // Initialize histograms -- xtals
  xtalEnergyHist_ = new TH1F("XtalEnergy","Crystal energy;GeV",110,-1,10);
  xtalTimeHist_ = new TH1F("XtalTime","Time of all crystals;ns",150,-75,75);
  xtalIEtaHist_ = new TH1F("xtalIEta","i#eta of crystal",171,-85,86);
  xtalIPhiHist_ = new TH1F("xtalIPhi","i#phi of crystal",361,1,361);
  xtalIXHist_ = new TH1F("xtalIX","ix of crystal",101,1,101);
  xtalIYHist_ = new TH1F("xtalIY","iy of crystal",101,1,101);
  //  xtalStatusHist_ = new TH1F("XtalStatus","Crystal status flag",16,0,15);
  // TH2F* xtalOccupancyHistEB_ = new TH2F("XtalOccupancyEB","Crystal occupancy;i#phi;i#eta",360,1.,361.,172,-86,86);

  // Initialize histograms -- BasicClusters
  BCNumPerEventHist_ = new TH1F("BCNumPerEvent","Number of BC per event",100,0,100);
  BCNumCrysHist_ = new TH1F("BCNumCrys","Number of crystals per BC",10,0,10);
  BCNumCrysOverThrHist_  = new TH1F("BCNumCrys over threshold","Number of crystals per BC over threshold",10,0,10);
  BCNumCrysOverThrEBHist_= new TH1F("EB BCNumCrys over threshold","EB Number of crystals per BC over threshold",10,0,10);
  BCNumCrysOverThrEEHist_= new TH1F("EE BCNumCrys over threshold","EE Number of crystals per BC over threshold",10,0,10);

  BCEnergyHist_ = new TH1F("BCEnergy","Energy of BCs;GeV",100,0,25);
  BCEtHist_ = new TH1F("BCEt","E_{T} of BCs;GeV",100,0,25);
  BCOccupancyEBHist_  = new TH2F("BCOccupancyEB","BC occupancy;i#eta;i#phi",171,-85,86,361,1.,361.);
  BCOccupancyEEHist_  = new TH2F("BCOccupancyEE","BC occupancy;ix;iy",101,1.,101.,101,1,101);
  BCOccupancyEEPHist_  = new TH2F("BCOccupancyEEP","BC occupancy EE+;ix;iy",101,1.,101.,101,1,101);
  BCOccupancyEEMHist_  = new TH2F("BCOccupancyEEM","BC occupancy EE-;ix;iy",101,1.,101.,101,1,101);
  BCOccupancyHistAny_ = new TH2F("BCOccupancyAny","BC occupancy;#eta;#phi",50,-3.5,3.5,50,-1*TMath::Pi(),TMath::Pi());
  BCEtaHist_ = new TH1F("Cluster #eta","#eta of cluster",171,-3.5,3.5);
  BCPhiHist_ = new TH1F("Cluster #phi","#phi of cluster",50,-1*TMath::Pi(),TMath::Pi());
  BCClusterShapeEEHist_ = new TH1F("EE cluster shape","e2x2 / e3x3",65,-0.1,1.2);
  BCClusterShapeEEPHist_ = new TH1F("EEP cluster shape","e2x2 / e3x3",65,-0.1,1.2);
  BCClusterShapeEEMHist_ = new TH1F("EEM cluster shape","e2x2 / e3x3",65,-0.1,1.2);
  BCClusterShapeEBHist_  = new TH1F("EB cluster shape","e2x2 / e3x3",65,-0.1,1.2);
  // Initialize histograms -- diphotons control plots
  massDiGammaHist_      = new TH1F("massDiGamma","m(#gamma#gamma)",50,0,0.500);
  massDiGammaEBHist_    = new TH1F("massDiGamma EB","m(#gamma#gamma) EB",50,0,0.500);
  massDiGammaEEHist_    = new TH1F("massDiGamma EE","m(#gamma#gamma) EE",50,0,0.500);
  massDiGammaEEPHist_   = new TH1F("massDiGamma EEP","m(#gamma#gamma) EEP",50,0,0.500);
  massDiGammaEEMHist_   = new TH1F("massDiGamma EEM","m(#gamma#gamma) EEM",50,0,0.500);
  diPhotonOccupancyAny_ = new TH2F("di-photon occupancy","di-photon occupancy;#eta;#phi",50,-3.5,3.5,50,-1*TMath::Pi(),TMath::Pi());
  // Initialize histograms -- single cluster resolution
  dtUpToQuarterGeVEB_   = new TH1F("EB #Delta(t),   A_{eff} up to 1/4 GeV", "EB #Delta(t),   ~5<A_{eff}/#sigma_{N}<6", 400, -20, 20); 
  dtUpToHalfGeVEB_      = new TH1F("EB #Delta(t),   A_{eff} up to half GeV", "EB #Delta(t),   6<A_{eff}/#sigma_{N}<12", 400, -20, 20); 
  dtUpToOneGeVEB_       = new TH1F("EB #Delta(t),   A_{eff} up to one GeV", "EB #Delta(t),   12<A_{eff}/#sigma_{N}<24", 200, -DtMax_, DtMax_); 
  dtUpToTwoGeVEB_       = new TH1F("EB #Delta(t),   A_{eff} up to two GeV", "EB #Delta(t),   24<A_{eff}/#sigma_{N}<48", 200, -DtMax_, DtMax_); 
  dtUpOverTwoGeVEB_     = new TH1F("EB #Delta(t),   A_{eff} over two GeV", "EB #Delta(t),   A_{eff}/#sigma_{N}>48", 200, -DtMax_, DtMax_); 

  dtUpToThreeQuarterGeVEE_ = new TH1F("EE #Delta(t),   A_{eff} up to 3/4 GeV", "EE #Delta(t),   A_{eff}/#sigma_{N}<6", 200, -DtMax_, DtMax_); 
  dtUpToOneAndHalfGeVEE_   = new TH1F("EE #Delta(t),   A_{eff} up to one&1/2 GeV", "EE #Delta(t),   6<A_{eff}/#sigma_{N}<12", 200, -DtMax_, DtMax_); 
  dtUpToThreeGeVEE_        = new TH1F("EE #Delta(t),   A_{eff} up to three GeV", "EE #Delta(t),   12<A_{eff}/#sigma_{N}<24", 200, -DtMax_, DtMax_); 
  dtUpToSixGeVEE_          = new TH1F("EE #Delta(t),   A_{eff} up to six GeV", "EE #Delta(t),   24<A_{eff}/#sigma_{N}<48", 200, -DtMax_, DtMax_); 
  dtUpOverSixGeVEE_        = new TH1F("EE #Delta(t),   A_{eff} over six GeV", "EE #Delta(t),   A_{eff}/#sigma_{N}>48", 200, -DtMax_, DtMax_); 

  //  dtVSAeffHistAny_ = new TH2F("#Delta(t) VS A_{eff}/#sigma_{N}","#Delta(t) VS A_{eff}/#sigma_{N}; A_{eff}/#sigma_{N}; #Delta(t) [ns]",AeffNBins_,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);//gfold
  dtVSAeffHistAny_ = new TH2F("#Delta(t) VS A_{eff}/#sigma_{N}","#Delta(t) VS A_{eff}/#sigma_{N}; A_{eff}/#sigma_{N}; #Delta(t) [ns]", AeffNBins_  ,AeffBins_,numDtBins_,-DtMax_,DtMax_);
  //  dtVSAeffHistEB_  = new TH2F("EB:  #Delta(t)  VS  A_{eff}/#sigma_{N}","EB:  #Delta(t)  VS  A_{eff}/#sigma_{N}; A_{eff}/#sigma_{N}; #Delta(t) [ns]",AeffNBins_,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);//gfold
  dtVSAeffHistEB_  = new TH2F("EB:  #Delta(t)  VS  A_{eff}/#sigma_{N}","EB:  #Delta(t)  VS  A_{eff}/#sigma_{N}; A_{eff}/#sigma_{N}; #Delta(t) [ns]", AeffNBins_  ,AeffBins_,numDtBins_,-DtMax_,DtMax_);
  //  dtVSAeffHistEE_  = new TH2F("EE:  #Delta(t)  VS  A_{eff}/#sigma_{N}","EE:  #Delta(t)  VS  A_{eff}/#sigma_{N}; A_{eff}/#sigma_{N}; #Delta(t) [ns]",AeffNBins_,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);//gfold
  dtVSAeffHistEE_  = new TH2F("EE:  #Delta(t)  VS  A_{eff}/#sigma_{N}","EE:  #Delta(t)  VS  A_{eff}/#sigma_{N}; A_{eff}/#sigma_{N}; #Delta(t) [ns]",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);
  dtVSAeffProfAny_ = new TProfile("#Delta(t)  VS  A_{eff}/#sigma_{N} prof","#Delta(t) VS A_{eff}/#sigma_{N} prof",AeffNBins_,0.,AeffMax_,-DtMax_,DtMax_);
  for (int v=0; v<AeffNBins_; v++){// build histograms for RMS and sigma of DeltaT for Any
    //float binLeft=(v*AeffMax_/AeffNBins_); float binRight=((v+1)*AeffMax_/AeffNBins_);//gfold
    float binLeft=AeffBins_[v]; float binRight=AeffBins_[v+1];
    sprintf (buffer_, "Aeff bin %d, [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSliceVSAeffAny_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_); 
    AeffBinCentersAny_[v]=0; AeffBinCentersErrAny_[v]=0; sigmaAeffAny_[v]=0; sigmaAeffErrAny_[v]=0;
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; A_{eff} bin center");
    AeffSliceAny_[v] = new TH1F(bufferTitle_.c_str(),bufferTitle_.c_str(),20,binLeft,binRight);
  }//end loop
  dtRMSVSAeffAny_  = new TH1F("RMS(#Delta(t)) VS   A_{eff}", "RMS(#Delta(t)) VS   A_{eff}; A_{eff}/#sigma_{N}; RMS(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtSigmaAeffAny_  = new TH1F("#sigma(#Delta(t)) VS   A_{eff}", "#sigma(#Delta(t)) VS   A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  for (int v=0; v<AeffNBins_; v++){// build histograms for RMS and sigma of DeltaT for EB
    //    float binLeft=(v*AeffMax_/AeffNBins_); float binRight=((v+1)*AeffMax_/AeffNBins_);//gfold
    float binLeft=AeffBins_[v]; float binRight=AeffBins_[v+1];
    sprintf (buffer_, "EB: Aeff bin %d, [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSliceVSAeffEB_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  
    AeffBinCentersEB_[v]=0; AeffBinCentersErrEB_[v]=0; sigmaAeffEB_[v]=0;  sigmaAeffErrEB_[v]=0;
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; A_{eff} bin center");
    AeffSliceEB_[v] = new TH1F(bufferTitle_.c_str(),bufferTitle_.c_str(),20,binLeft,binRight);
  }//end loop
  dtRMSVSAeffEB_  = new TH1F("EB: RMS(#Delta(t)) VS   A_{eff}", "EB: RMS(#Delta(t)) VS   A_{eff}; A_{eff}/#sigma_{N}; RMS(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtSigmaAeffEB_  = new TH1F("EB: #sigma(#Delta(t)) VS   A_{eff}", "EB: #sigma(#Delta(t)) VS   A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  for (int v=0; v<AeffNBins_; v++){// build histograms for RMS and sigma of DeltaT for EE
    //float binLeft=(v*AeffMax_/AeffNBins_); float binRight=((v+1)*AeffMax_/AeffNBins_);//gfold
    float binLeft=AeffBins_[v]; float binRight=AeffBins_[v+1];
    sprintf (buffer_, "EE: Aeff bin %d, [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSliceVSAeffEE_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  
    AeffBinCentersEE_[v]=0; AeffBinCentersErrEE_[v]=0; sigmaAeffEE_[v]=0;  sigmaAeffErrEE_[v]=0;
    AeffSliceEE_[v] = new TH1F(bufferTitle_.c_str(),bufferTitle_.c_str(),20,binLeft,binRight);
  }//end loop
  dtRMSVSAeffEE_  = new TH1F("EE: RMS(#Delta(t)) VS   A_{eff}", "EE: RMS(#Delta(t)) VS   A_{eff}; A_{eff}/#sigma_{N}; RMS(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtSigmaAeffEE_  = new TH1F("EE: #sigma(#Delta(t)) VS   A_{eff}", "EE: #sigma(#Delta(t)) VS   A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtVSAeffProfEB_  = new TProfile("EB:  #Delta(t)   VS  A_{eff}/#sigma_{N} prof","EB:  #Delta(t)  VS  A_{eff}/#sigma_{N} prof",AeffNBins_,AeffBins_,-DtMax_,DtMax_);
  dtVSAeffProfEE_  = new TProfile("EE:  #Delta(t)   VS  A_{eff}/#sigma_{N} prof","EE:  #Delta(t)  VS  A_{eff}/#sigma_{N} prof",AeffNBins_,AeffBins_,-DtMax_,DtMax_);
  singleClusterChi2HistAny_ = new TH1F("clusterChi2Any","#Chi^{2} of crystal times in a cluster (any)",numChi2Bins,0,chi2Max);
  singleClusterChi2NDFHistAny_ = new TH1F("clusterChi2NDFAny","#Chi^{2}/NDF of crystal times in a cluster (any)",numChi2NDFBins,0,chi2NDFMax);
  singleClusterChi2HistEB_ = new TH1F("clusterChi2EB","#Chi^{2} of crystal times in a cluster (EB)",numChi2Bins,0,chi2Max);
  singleClusterChi2NDFHistEB_ = new TH1F("clusterChi2NDFEB","#Chi^{2}/NDF of crystal times in a cluster (EB)",numChi2NDFBins,0,chi2NDFMax);
  singleClusterChi2HistEE_ = new TH1F("clusterChi2EE","#Chi^{2} of crystal times in a cluster (EE)",numChi2Bins,0,chi2Max);
  singleClusterChi2NDFHistEE_ = new TH1F("clusterChi2NDFEE","#Chi^{2}/NDF of crystal times in a cluster (EE)",numChi2NDFBins,0,chi2NDFMax);
  //variable bins stop here. FOR NOW.

  deltaTCrysVsAmplitudeAny_ = new TH2F("deltaTCrysVsAmplitudeAny","#Delta(t_{cry2}-t_{cry1}), A/#sigma_{cry1} > 30 vs. A/#sigma_{cry2} (EB/EE);A_{2}/#sigma_{2};ns",50,0,100,1000,-50,50);
  deltaTCrysVsAmplitudeEB_ = new TH2F("deltaTCrysVsAmplitudeEB","#Delta(t_{cry2}-t_{cry1}), A/#sigma_{cry1} > 30 vs. A/#sigma_{cry2} (EB);A_{2}/#sigma_{2};ns",50,0,100,1000,-50,50);
  deltaTCrysVsAmplitudeEE_ = new TH2F("deltaTCrysVsAmplitudeEE","#Delta(t_{cry2}-t_{cry1}), A/#sigma_{cry1} > 30 vs. A/#sigma_{cry2} (EE);A_{2}/#sigma_{2};ns",50,0,100,1000,-50,50);

  // Initialize histograms -- selection on pi0 candidates 
  diPhotonPeakOccupancyAny_     = new TH2F("#pi_{0} occupancy (di-photon peak)","di-photon peak;#eta;#phi",50,-3.5,3.5,50,-1*TMath::Pi(),TMath::Pi());
  diPhotonSidesOccupancyAny_    = new TH2F("di-photon side-bands","di-photon side-bands;#eta;#phi",50,-3.5,3.5,50,-1*TMath::Pi(),TMath::Pi());
  // Initialize histograms -- single cluster resolution in pi0 peak
  dtUpToQuarterGeVEBPeak_   = new TH1F("EBPeak: #Delta(t) A_{eff} up to 1/4 GeV", "EBPeak: #Delta(t) ~5<A_{eff}/#sigma_{N}<6", 400, -20, 20); 
  dtUpToHalfGeVEBPeak_      = new TH1F("EBPeak: #Delta(t) A_{eff} up to half GeV", "EBPeak: #Delta(t) 6<A_{eff}/#sigma_{N}<12", 400, -20, 20); 
  dtUpToOneGeVEBPeak_       = new TH1F("EBPeak: #Delta(t) A_{eff} up to one GeV", "EBPeak: #Delta(t) 12<A_{eff}/#sigma_{N}<24", 200, -DtMax_, DtMax_); 
  dtUpToTwoGeVEBPeak_       = new TH1F("EBPeak: #Delta(t) A_{eff} up to two GeV", "EBPeak: #Delta(t) 24<A_{eff}/#sigma_{N}<48", 200, -DtMax_, DtMax_); 
  dtUpOverTwoGeVEBPeak_     = new TH1F("EBPeak: #Delta(t) A_{eff} over two GeV", "EBPeak: #Delta(t) A_{eff}/#sigma_{N}>48", 200, -DtMax_, DtMax_); 

  dtUpToThreeQuarterGeVEEPeak_ = new TH1F("EEPeak: #Delta(t) A_{eff} up to 3/4 GeV", "EEPeak: #Delta(t) A_{eff}/#sigma_{N}<6", 200, -DtMax_, DtMax_); 
  dtUpToOneAndHalfGeVEEPeak_   = new TH1F("EEPeak: #Delta(t) A_{eff} up to one&1/2 GeV", "EEPeak: #Delta(t) 6<A_{eff}/#sigma_{N}<12", 200, -DtMax_, DtMax_); 
  dtUpToThreeGeVEEPeak_        = new TH1F("EEPeak: #Delta(t) A_{eff} up to three GeV", "EEPeak: #Delta(t) 12<A_{eff}/#sigma_{N}<24", 200, -DtMax_, DtMax_); 
  dtUpToSixGeVEEPeak_          = new TH1F("EEPeak: #Delta(t) A_{eff} up to six GeV", "EEPeak: #Delta(t) 24<A_{eff}/#sigma_{N}<48", 200, -DtMax_, DtMax_); 
  dtUpOverSixGeVEEPeak_        = new TH1F("EEPeak: #Delta(t) A_{eff} over six GeV", "EEPeak: #Delta(t) A_{eff}/#sigma_{N}>48", 200, -DtMax_, DtMax_); 

  dtVSAeffHistAnyPeak_ = new TH2F("Peak: #Delta(t) VS A_{eff}/#sigma_{N}","Peak: #Delta(t) VS A_{eff}/#sigma_{N}",AeffNBins_,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);
  dtVSAeffHistEBPeak_  = new TH2F("EBPeak: #Delta(t) VS A_{eff}/#sigma_{N}","EBPeak #Delta(t) VS A_{eff}/#sigma_{N}",AeffNBins_,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);
  dtVSAeffHistEEPeak_  = new TH2F("EEPeak: E#Delta(t) VS A_{eff}/#sigma_{N}","EEPeak: #Delta(t) VS A_{eff}/#sigma_{N}",AeffNBins_,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);
  for (int v=0; v<AeffNBins_; v++){// build histograms for RMS and sigma of DeltaT for Any
    //    float binLeft=(v*AeffMax_/AeffNBins_); float binRight=((v+1)*AeffMax_/AeffNBins_);//gfold
    float binLeft=AeffBins_[v]; float binRight=AeffBins_[v+1];
    sprintf (buffer_, "Peak: #Deltat bin %d, [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSliceVSAeffAnyPeak_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  }
  dtRMSVSAeffAnyPeak_  = new TH1F("RMS(#Delta(t)) VS   A_{eff}", "Peak: RMS(#Delta(t)) VS   A_{eff}; A_{eff}/#sigma_{N}; RMS(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtSigmaAeffAnyPeak_  = new TH1F("#sigma(#Delta(t)) VS   A_{eff}", "Peak: #sigma(#Delta(t)) VS   A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  for (int v=0; v<AeffNBins_; v++){// build histograms for RMS and sigma of DeltaT for EB
    //    float binLeft=(v*AeffMax_/AeffNBins_); float binRight=((v+1)*AeffMax_/AeffNBins_);//gfold
    float binLeft=AeffBins_[v]; float binRight=AeffBins_[v+1];
    sprintf (buffer_, "EBPeak: #Deltat bin %d, [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSliceVSAeffEBPeak_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  }
  dtRMSVSAeffEBPeak_  = new TH1F("EBpeak: RMS(#Delta(t)) VS   A_{eff}", "EBPeak: RMS(#Delta(t)) VS   A_{eff}; A_{eff}/#sigma_{N}; RMS(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtSigmaAeffEBPeak_  = new TH1F("EBpeak: #sigma(#Delta(t)) VS   A_{eff}", "EBPeak: #sigma(#Delta(t)) VS   A_{eff} peak; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);
  for (int v=0; v<AeffNBins_; v++){// build histograms for RMS and sigma of DeltaT for EE
    //float binLeft=(v*AeffMax_/AeffNBins_); float binRight=((v+1)*AeffMax_/AeffNBins_);//gfold
    float binLeft=AeffBins_[v]; float binRight=AeffBins_[v+1];
    sprintf (buffer_, "EEPeak: #Deltat bin %d, [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSliceVSAeffEEPeak_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  }
  dtRMSVSAeffEEPeak_  = new TH1F("EEPeak: RMS(#Delta(t)) VS   A_{eff}", "EE: RMS(#Delta(t)) VS   A_{eff}; A_{eff}/#sigma_{N}; RMS(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtSigmaAeffEEPeak_  = new TH1F("EEPeak: #sigma(#Delta(t)) VS   A_{eff}", "EE: #sigma(#Delta(t)) VS   A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtVSAeffProfAnyPeak_ = new TProfile("Peak: #Delta(t) VS A_{eff}/#sigma_{N} prof","Peak: #Delta(t) VS A_{eff}/#sigma_{N} prof",AeffNBins_,AeffBins_,-DtMax_,DtMax_);
  dtVSAeffProfEBPeak_  = new TProfile("EBPeak: #Delta(t) VS A_{eff}/#sigma_{N} prof","EBPeak #Delta(t) VS A_{eff}/#sigma_{N} prof",AeffNBins_,AeffBins_,-DtMax_,DtMax_);
  dtVSAeffProfEEPeak_  = new TProfile("EEPeak: #Delta(t) VS A_{eff}/#sigma_{N} prof","EEPeak: #Delta(t) VS A_{eff}/#sigma_{N} prof",AeffNBins_,AeffBins_,-DtMax_,DtMax_);
  singleClusterChi2HistAnyPeak_ = new TH1F("clusterChi2AnyPeak","#Chi^{2} of crystal times in a cluster (any peak)",numChi2Bins,0,chi2Max);
  singleClusterChi2NDFHistAnyPeak_ = new TH1F("clusterChi2NDFAnyPeak","#Chi^{2}/NDF of crystal times in a cluster (any peak)",numChi2NDFBins,0,chi2NDFMax);
  singleClusterChi2HistEBPeak_ = new TH1F("clusterChi2EBPeak","#Chi^{2} of crystal times in a cluster (EB peak)",numChi2Bins,0,chi2Max);
  singleClusterChi2NDFHistEBPeak_ = new TH1F("clusterChi2NDFEBPeak","#Chi^{2}/NDF of crystal times in a cluster (EB peak)",numChi2NDFBins,0,chi2NDFMax);
  singleClusterChi2HistEEPeak_ = new TH1F("clusterChi2EEPeak","#Chi^{2} of crystal times in a cluster (EE peak)",numChi2Bins,0,chi2Max);
  singleClusterChi2NDFHistEEPeak_ = new TH1F("clusterChi2NDFEEPeak","#Chi^{2}/NDF of crystal times in a cluster (EE peak)",numChi2NDFBins,0,chi2NDFMax);
  dtPullSingleClusterHistPi0Peak_  = new TH1F("DeltaTPullSingleClusterPeak","#Delta(t)/#sigma(t) between two crystals EB/EE  under #pi^{0} mass peak (pull)",100,-5,5);
  dtPullSingleClusterHistPi0PeakEB_  = new TH1F("DeltaTPullSingleClusterPeakEB","#Delta(t)/#sigma(t) between two crystals EB  under #pi^{0} mass peak (pull)",100,-5,5);
  dtPullSingleClusterHistPi0PeakEE_  = new TH1F("DeltaTPullSingleClusterPeakEE","#Delta(t)/#sigma(t) between two crystals EE  under #pi^{0} mass peak (pull)",100,-5,5);

  // should these DeltaT be vary in [-DtMax_, DtMax_] ? Once fixed/understood
  // Initialize histograms -- double cluster resolution
  dtDoubleClusterHistAny_     = new TH1F("DeltaTDoubleClusterAny","#Delta(t) between two clusters EB/EE",100,-25,25);
  dtPullDoubleClusterHistAny_ = new TH1F("DeltaTPullDoubleClusterAny","#Delta(t)/#sigma(t) between two clusters EB/EE (pull)",100,-5,5);
  dtPullChi2CutDoubleClusterHistAny_ = new TH1F("DeltaTPullChi2CutDoubleClusterAny","#Delta(t)/#sigma(t) between two clusters EB/EE (pull) with #Chi^{2} cut",100,-5,5);
  dtVsPtDoubleClusterHistAny_ = new TH2F("DeltaTVSPtDoubleClusterAny","#Delta(t)  between two clusters EB/EE VS P_{t}(di-photon) ",50,0,10,50,-25,25);

  dtDoubleClusterHistPi0Peak_  = new TH1F("DeltaTDoubleClusterPi0Peak","#Delta(t) between two clusters under #pi^{0} mass peak EB/EE; #Delta(#Delta(t) [ns]t) [ns]",100,-25,25);
  dtDoubleClusterHistPi0PeakEE_= new TH1F("DeltaTDoubleClusterPi0PeakEE","#Delta(t) between two EE clusters under #pi^{0} mass peak; #Delta(t) [ns]",100,-25,25);
  dtDoubleClusterHistPi0PeakEB_= new TH1F("DeltaTDoubleClusterPi0PeakEB","#Delta(t) between two EB clusters under #pi^{0} mass peak; #Delta(t) [ns]",100,-25,25);
  dtPullDoubleClusterHistPi0Peak_  = new TH1F("DeltaTPullDoubleClusterPeak","#Delta(t)/#sigma(t) between two clusters EB/EE  under #pi^{0} mass peak  (pull)",100,-5,5);
  dtPullDoubleClusterHistPi0PeakEE_= new TH1F("DeltaTPullDoubleClusterPeakEE","#Delta(t)/#sigma(t) between two clusters EE  under #pi^{0} mass peak (pull)",100,-5,5);
  dtPullDoubleClusterHistPi0PeakEB_= new TH1F("DeltaTPullDoubleClusterPeakEB","#Delta(t)/#sigma(t) between two clusters EB  under #pi^{0} mass peak (pull)",100,-5,5);
  dtPullChi2CutDoubleClusterHistPi0Peak_  = new TH1F("DeltaTPullCHi2CutDoubleClusterPeak","#Delta(t)/#sigma(t) between two clusters EB/EE (pull) with #Chi^{2} cut",100,-5,5);
  dtPullChi2CutDoubleClusterHistPi0PeakEE_= new TH1F("DeltaTPullCHi2CutDoubleClusterPeakEE","#Delta(t)/#sigma(t) between two clusters EE (pull) with #Chi^{2} cut",100,-5,5);
  dtPullChi2CutDoubleClusterHistPi0PeakEB_= new TH1F("DeltaTPullCHi2CutDoubleClusterPeakEB","#Delta(t)/#sigma(t) between two clusters EB (pull) with #Chi^{2} cut",100,-5,5);
  dtVsPtDoubleClusterHistPi0Peak_  = new TH2F("DeltaTVSPtDoubleClusterPeak","#Delta(t) between two clusters EB/EE VS P_{t}(#pi_{0}) ",50,0,10,50,-25,25);
  dtVsPtDoubleClusterHistPi0PeakEE_= new TH2F("DeltaTVSPtDoubleClusterPeakEE","#Delta(t) between two clusters EE VS P_{t}(#pi_{0}) ",50,0,10,50,-25,25);
  dtVsPtDoubleClusterHistPi0PeakEB_= new TH2F("DeltaTVSPtDoubleClusterPeakEB","#Delta(t) between two clusters EB VS P_{t}(#pi_{0}) ",50,0,10,50,-25,25);

  dtSeedsVsAeffDoubleClusterHist_  = new TH2F("DeltaTSeedsVSAeffDoubleClusterPeak","#Delta(t) between clusters seeds VS A_{eff} of seeds) ",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);
  dtSeedsVsAeffDoubleClusterHistEE_= new TH2F("DeltaTSeedsVSAeffDoubleClusterPeakEE","#Delta(t) between clusters seeds EE VS A_{eff} of seeds",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);
  dtSeedsVsAeffDoubleClusterHistEB_= new TH2F("DeltaTSeedsVSAeffDoubleClusterPeakEB","#Delta(t) between clusters seeds EB VS VS A_{eff} of seeds",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);

  dtSeedsVsAeffDoubleClusterHistPi0Peak_  = new TH2F("DeltaTSeedsVSAeffDoubleClusterPeak","#Delta(t) between clusters seeds VS A_{eff} of seeds",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);
  dtSeedsVsAeffDoubleClusterHistPi0PeakEE_= new TH2F("DeltaTSeedsVSAeffDoubleClusterPeakEE","#Delta(t) between clusters seeds VS A_{eff} of seeds",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);
  dtSeedsVsAeffDoubleClusterHistPi0PeakEB_= new TH2F("DeltaTSeedsVSAeffDoubleClusterPeakEB","#Delta(t) between clusters seeds VS A_{eff} of seeds",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);


  dtVsAeffDoubleClusterHist_  = new TH2F("DeltaTVSAeffDoubleClusterPeak","#Delta(t) between clusters VS A_{eff} of seeds) ",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);
  dtVsAeffDoubleClusterHistEE_= new TH2F("DeltaTVSAeffDoubleClusterPeakEE","#Delta(t) between clusters EE VS A_{eff} of seeds",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);
  dtVsAeffDoubleClusterHistEB_= new TH2F("DeltaTVSAeffDoubleClusterPeakEB","#Delta(t) between clusters EB VS A_{eff} of seeds",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);

  dtVsAeffDoubleClusterHistPi0Peak_  = new TH2F("DeltaTVSAeffDoubleClusterPeak","#Delta(t) between clusters VS A_{eff} of seeds",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);
  dtVsAeffDoubleClusterHistPi0PeakEE_= new TH2F("DeltaTVSAeffDoubleClusterPeakEE","#Delta(t) between clusters VS A_{eff} of seeds",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);
  dtVsAeffDoubleClusterHistPi0PeakEB_= new TH2F("DeltaTVSAeffDoubleClusterPeakEB","#Delta(t) between clusters VS A_{eff} of seeds",AeffNBins_,AeffBins_,numDtBins_,-DtMax_,DtMax_);

  
  // double cluster studies VS Aeff
  dtSeedsSigmaAeffAny_= new TH1F("#sigma(#Delta(t)) VS A_{eff} (seeds)", "#sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtSeedsSigmaAeffEE_= new TH1F("EE #sigma(#Delta(t)) VS A_{eff} (seeds)", "EE #sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtSeedsSigmaAeffEB_= new TH1F("EB #sigma(#Delta(t)) VS A_{eff} (seeds)", "EB #sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtSeedsSigmaAeffPi0Peak_= new TH1F("Peak #sigma(#Delta(t)) VS A_{eff} (seeds)", "Peak #sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtSeedsSigmaAeffPi0PeakEE_= new TH1F("EE Peak #sigma(#Delta(t)) VS A_{eff} (seeds)", "EE Peak #sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtSeedsSigmaAeffPi0PeakEB_= new TH1F("EB Peak #sigma(#Delta(t)) VS A_{eff} (seeds)", "EB Peak #sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtClustersSigmaAeffAny_= new TH1F("#sigma(#Delta(t)) VS A_{eff} (clusters)", "#sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtClustersSigmaAeffEE_= new TH1F("EE #sigma(#Delta(t)) VS A_{eff} (clusters)", "EE #sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtClustersSigmaAeffEB_= new TH1F("EB #sigma(#Delta(t)) VS A_{eff} (clusters)", "EB #sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtClustersSigmaAeffPi0Peak_= new TH1F("EE Peak #sigma(#Delta(t)) VS A_{eff} (clusters)", "EE Peak #sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtClustersSigmaAeffPi0PeakEE_= new TH1F("EE Peak #sigma(#Delta(t)) VS A_{eff} (clusters)", "EE Peak #sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  
  dtClustersSigmaAeffPi0PeakEB_= new TH1F("EB Peak #sigma(#Delta(t)) VS A_{eff} (clusters)", "EB Peak #sigma(#Delta(t)) VS A_{eff}; A_{eff}/#sigma_{N}; #sigma(#Delta(t)) [ns]",AeffNBins_,AeffBins_);  

  
  for (int v=0; v<AeffNBins_; v++){// build slice histograms DeltaT Seeds/Cluster VS Aeff
    //    float binLeft=(v*AeffMax_/AeffNBins_); float binRight=((v+1)*AeffMax_/AeffNBins_);//gfold
    float binLeft=AeffBins_[v]; float binRight=AeffBins_[v+1];

    // slices for seeds
    sprintf (buffer_, "#Deltat seeds bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSeedsSlicesVsAeffDoubleClusterHist_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  

    sprintf (buffer_, "EE #Deltat seeds bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSeedsSlicesVsAeffDoubleClusterHistEE_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  

    sprintf (buffer_, "EB #Deltat seeds bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSeedsSlicesVsAeffDoubleClusterHistEB_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  

    sprintf (buffer_, "Peak #Deltat seeds bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSeedsSlicesVsAeffDoubleClusterHistPi0Peak_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  

    sprintf (buffer_, "Peak EE #Deltat seeds bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEE_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  

    sprintf (buffer_, "Peak EB #Deltat seeds bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEB_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  


    // slices for clusters
    sprintf (buffer_, "#Deltat cluster bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSlicesVsAeffDoubleClusterHist_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  

    sprintf (buffer_, "EE #Deltat cluster bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSlicesVsAeffDoubleClusterHistEE_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  

    sprintf (buffer_, "EB #Deltat cluster bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSlicesVsAeffDoubleClusterHistEB_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  

    sprintf (buffer_, "Peak #Deltat cluster bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSlicesVsAeffDoubleClusterHistPi0Peak_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  

    sprintf (buffer_, "Peak EE #Deltat cluster bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSlicesVsAeffDoubleClusterHistPi0PeakEE_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  

    sprintf (buffer_, "Peak EB #Deltat cluster bin %d, A_{eff} [%4.1f,%4.1f)", v+1, binLeft, binRight);
    bufferTitle_.erase(); bufferTitle_=std::string(buffer_)+std::string("; #Deltat [ns]");
    dtSlicesVsAeffDoubleClusterHistPi0PeakEB_[v] = new TH1F(buffer_,bufferTitle_.c_str(),numDtBins_,-DtMax_,DtMax_);  

  }//end loop to inizialize cluster dt slices histograms 

}//end initializeHists



// ---------------------------------------------------------------------------------------
// ------------------ Function to plot the control hists ---------------------------------
void doControlHists()
{
  BCNumPerEventHist_->Fill(treeVars_.nClusters);
  // loop on basic clusters
  for (int bCluster=0; bCluster < treeVars_.nClusters; bCluster++)
  {

    float eBC=0; // calculate energy of BC for validation
    for (int cryInBC=0; cryInBC < treeVars_.nXtalsInCluster[bCluster]; cryInBC++){
      eBC+= treeVars_.xtalInBCEnergy[bCluster][cryInBC];}

    BCEnergyHist_->Fill(treeVars_.clusterEnergy[bCluster]);
    BCEtHist_->Fill(treeVars_.clusterTransverseEnergy[bCluster]);
    BCNumCrysHist_->Fill(treeVars_.nXtalsInCluster[bCluster]);

    // basic cluster occupancy in physics coordinates
    BCOccupancyHistAny_ -> Fill(treeVars_.clusterEta[bCluster],treeVars_.clusterPhi[bCluster]);
    BCEtaHist_ -> Fill(treeVars_.clusterEta[bCluster]);
    BCPhiHist_ -> Fill(treeVars_.clusterPhi[bCluster]);
    
    bool thisIsEB = false;
    //  basic cluster occupancy in detector coordinates, using first cry of BC as a representative
    if(treeVars_.xtalInBCIEta[bCluster][0] != -999999){                                        // this is EB; ieta=-999999 tags EE
      BCOccupancyEBHist_->Fill(treeVars_.xtalInBCIEta[bCluster][0],treeVars_.xtalInBCIPhi[bCluster][0]); thisIsEB = true;   }
    else if (treeVars_.xtalInBCIx[bCluster][0] != -999999 && treeVars_.clusterEta[bCluster]>0){ // this is EEP; ix=-999999 tags EB
      BCOccupancyEEHist_ ->Fill(treeVars_.xtalInBCIx[bCluster][0],treeVars_.xtalInBCIy[bCluster][0]);
      BCOccupancyEEPHist_->Fill(treeVars_.xtalInBCIx[bCluster][0],treeVars_.xtalInBCIy[bCluster][0]);
      thisIsEB = false;   }
    else if (treeVars_.xtalInBCIx[bCluster][0] != -999999 && treeVars_.clusterEta[bCluster]<0){ // this is EEM; ix=-999999 tags EB
      BCOccupancyEEHist_ ->Fill(treeVars_.xtalInBCIx[bCluster][0],treeVars_.xtalInBCIy[bCluster][0]);
      BCOccupancyEEMHist_->Fill(treeVars_.xtalInBCIx[bCluster][0],treeVars_.xtalInBCIy[bCluster][0]);
      thisIsEB = false;   }

      if (speak_)  std::cout << "\tbCluster: num "               << bCluster 
        << "\t eBC: "                      << treeVars_.clusterEnergy[bCluster]
          << "\t eBC_predicted: "            << eBC
          //<< "\n\t et: "                     << treeVars_.clusterTransverseEnergy[bCluster]
          //<< "\t predicted et: "             << treeVars_.clusterEnergy[bCluster]*sin(2*atan(exp(-1* treeVars_.clusterEta[bCluster] )) )
          << " eta: "                        << treeVars_.clusterEta[bCluster]
          << "\n\t num crystals: "           << treeVars_.nXtalsInCluster[bCluster]
          << "\n\t\tfirst crystal:  \tieta " << treeVars_.xtalInBCIEta[bCluster][0] 
          << "\teta "                        << treeVars_.xtalInBCEta[bCluster][0] 
          << " \t energy "                   << treeVars_.xtalInBCEnergy[bCluster][0] 
          << " \t ADC "                      << treeVars_.xtalInBCAmplitudeADC[bCluster][0] 
          << " \t time "                     << treeVars_.xtalInBCTime[bCluster][0] 
          << std::endl;

      // count number of crystals in a BC over threshold
      int numCryOverThreshold=0; 
      for(int thisCry=0; thisCry<treeVars_.nXtalsInCluster[bCluster]; thisCry++)
      {
        if (treeVars_.xtalInBCIEta[bCluster][thisCry]!=-999999)  xtalIEtaHist_ -> Fill (treeVars_.xtalInBCIEta[bCluster][thisCry]);
        if (treeVars_.xtalInBCIPhi[bCluster][thisCry]!=-999999)  xtalIPhiHist_ -> Fill (treeVars_.xtalInBCIPhi[bCluster][thisCry]);
        if (treeVars_.xtalInBCIx[bCluster][thisCry]  !=-999999)  xtalIXHist_   -> Fill (treeVars_.xtalInBCIx[bCluster][thisCry]);
        if (treeVars_.xtalInBCIy[bCluster][thisCry]  !=-999999)  xtalIYHist_   -> Fill (treeVars_.xtalInBCIy[bCluster][thisCry]);
	xtalEnergyHist_                                                        -> Fill (treeVars_.xtalInBCEnergy[bCluster][thisCry]);
	
	if (thisIsEB &&                                                           // this is barrel
	    (treeVars_.xtalInBCAmplitudeADC[bCluster][thisCry]/sigmaNoiseEB) > minAmpliOverSigma_ ) {  
	  xtalIEtaHist_ -> Fill (treeVars_.xtalInBCIEta[bCluster][thisCry]);
	  numCryOverThreshold++;	 } 
	else if ( (!thisIsEB) &&                                                   // this is endcap
		  (treeVars_.xtalInBCAmplitudeADC[bCluster][thisCry]/sigmaNoiseEE) > minAmpliOverSigma_ ) {  
	  xtalIEtaHist_ -> Fill (treeVars_.xtalInBCIEta[bCluster][thisCry]);
	  numCryOverThreshold++;	  }
	
      }//end loop on crystals

  
      if(thisIsEB){
	BCNumCrysOverThrHist_  ->Fill(numCryOverThreshold);
	BCNumCrysOverThrEBHist_->Fill(numCryOverThreshold);}
      else{
	BCNumCrysOverThrHist_  ->Fill(numCryOverThreshold);
	BCNumCrysOverThrEEHist_->Fill(numCryOverThreshold);}
     
      
  }//end loop on basic clusters
}// end doControlHistograms

// ---------------------------------------------------------------------------------------
// ------------------ Function to write hists --------------------------------------------
void writeHists()
{
  saving_->cd();
  // write out control histograms
  TDirectory *controlPlots = saving_->mkdir("control");
  controlPlots->cd();
  BCNumPerEventHist_->Write();
  BCEnergyHist_->Write();
  BCEtHist_->Write();
  BCNumCrysHist_->Write();
  BCNumCrysOverThrHist_  ->Write();
  BCNumCrysOverThrEBHist_->Write();
  BCNumCrysOverThrEEHist_->Write();
  BCOccupancyEBHist_->Write();
  BCOccupancyEEHist_->Write();
  BCOccupancyEEPHist_->Write();
  BCOccupancyEEMHist_->Write();
  BCOccupancyHistAny_->Write();
  BCEtaHist_->Write();
  BCPhiHist_->Write();
  BCClusterShapeEEHist_->Write();
  BCClusterShapeEEPHist_->Write();
  BCClusterShapeEEMHist_->Write();
  BCClusterShapeEBHist_->Write();

  xtalEnergyHist_->Write(); 
  xtalTimeHist_->Write();
  xtalIEtaHist_->Write();
  xtalIPhiHist_->Write();
  xtalIXHist_ ->Write();
  xtalIYHist_ ->Write();
  
  // xtalStatusHist_->Write();
  // xtalOccupancyHistEB_->Write();
  massDiGammaHist_-> Write(); 
  massDiGammaEBHist_-> Write(); 
  massDiGammaEEHist_-> Write(); 
  massDiGammaEEPHist_-> Write(); 
  massDiGammaEEMHist_-> Write(); 

  diPhotonOccupancyAny_      -> Write();
  diPhotonPeakOccupancyAny_  -> Write();
  diPhotonSidesOccupancyAny_ -> Write();

  // write out single cluster resolution plots 
  TDirectory *singleClusResolution = saving_->mkdir("single-resolution");
  singleClusResolution->cd();
  
  dtVSAeffHistAny_-> Write(); 
  dtVSAeffHistEB_ -> Write(); 
  dtVSAeffHistEE_ -> Write(); 
  dtVSAeffProfAny_-> Write(); 
  dtVSAeffProfEB_ -> Write(); 
  dtVSAeffProfEE_ -> Write(); 
  
  dtUpToQuarterGeVEB_-> Write(); 
  dtUpToHalfGeVEB_   -> Write(); 
  dtUpToOneGeVEB_    -> Write(); 
  dtUpToTwoGeVEB_    -> Write(); 
  dtUpOverTwoGeVEB_  -> Write(); 
				    
  dtUpToThreeQuarterGeVEE_-> Write(); 
  dtUpToOneAndHalfGeVEE_  -> Write(); 
  dtUpToThreeGeVEE_       -> Write(); 
  dtUpToSixGeVEE_         -> Write(); 
  dtUpOverSixGeVEE_       -> Write(); 

  dtRMSVSAeffAny_    -> Write();
  dtSigmaAeffAny_    -> Write();
  TGraphErrors *     dtSigmaAeffGraphAny_ = new TGraphErrors(AeffNBins_ , AeffBinCentersAny_, sigmaAeffAny_, AeffBinCentersErrAny_, sigmaAeffErrAny_);
  dtSigmaAeffGraphAny_->SetName("sigmaVsAeff Any");
  dtSigmaAeffGraphAny_->SetTitle(dtSigmaAeffAny_->GetTitle());
  dtSigmaAeffGraphAny_->GetXaxis()->SetTitle(dtSigmaAeffAny_->GetXaxis()->GetTitle());
  dtSigmaAeffGraphAny_->GetYaxis()->SetTitle(dtSigmaAeffAny_->GetYaxis()->GetTitle());
  dtSigmaAeffGraphAny_->Write();

  dtRMSVSAeffEB_     -> Write();
  dtSigmaAeffEB_     -> Write();
  TGraphErrors *     dtSigmaAeffGraphEB_ = new TGraphErrors(AeffNBins_ , AeffBinCentersEB_, sigmaAeffEB_, AeffBinCentersErrEB_, sigmaAeffErrEB_);
  dtSigmaAeffGraphEB_->SetName("sigmaVsAeff EB");
  dtSigmaAeffGraphEB_->SetTitle(dtSigmaAeffEB_->GetTitle());
  dtSigmaAeffGraphEB_->GetXaxis()->SetTitle(dtSigmaAeffEB_->GetXaxis()->GetTitle());
  dtSigmaAeffGraphEB_->GetYaxis()->SetTitle(dtSigmaAeffEB_->GetYaxis()->GetTitle());
  dtSigmaAeffGraphEB_->Write();

  dtRMSVSAeffEE_     -> Write();
  dtSigmaAeffEE_     -> Write();
  TGraphErrors *     dtSigmaAeffGraphEE_ = new TGraphErrors(AeffNBins_ , AeffBinCentersEE_, sigmaAeffEE_, AeffBinCentersErrEE_, sigmaAeffErrEE_);
  dtSigmaAeffGraphEE_->SetName("sigmaVsAeff EE");
  dtSigmaAeffGraphEE_->SetTitle(dtSigmaAeffEE_->GetTitle());
  dtSigmaAeffGraphEE_->GetXaxis()->SetTitle(dtSigmaAeffEE_->GetXaxis()->GetTitle());
  dtSigmaAeffGraphEE_->GetYaxis()->SetTitle(dtSigmaAeffEE_->GetYaxis()->GetTitle());
  dtSigmaAeffGraphEE_->Write();


  singleClusterChi2HistAny_->Write();
  singleClusterChi2NDFHistAny_->Write();
  singleClusterChi2HistEB_->Write();
  singleClusterChi2NDFHistEB_->Write();
  singleClusterChi2HistEE_->Write();
  singleClusterChi2NDFHistEE_->Write();

  deltaTCrysVsAmplitudeAny_->Write();
  deltaTCrysVsAmplitudeEB_->Write();
  deltaTCrysVsAmplitudeEE_->Write();

  // write out 1-d control plots for DeltaT RMS and sigma for Any
  TDirectory *singleClusResolutionSlices = singleClusResolution->mkdir("dtslices-any");
  singleClusResolutionSlices->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSliceVSAeffAny_[v] -> Write();  }

  // write out 1-d control plots for DeltaT RMS and sigma for EB
  TDirectory *singleClusResolutionEBSlices = singleClusResolution->mkdir("dtslices-EB");
  singleClusResolutionEBSlices->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSliceVSAeffEB_[v] -> Write();  }

  // write out 1-d control plots for DeltaT RMS and sigma for EE
  TDirectory *singleClusResolutionEESlices = singleClusResolution->mkdir("dtslices-EE");
  singleClusResolutionEESlices->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSliceVSAeffEE_[v] -> Write();  }


  // write out single cluster resolution plots after pi0 selection
  TDirectory *singleClusResolutionPi0Clusters = saving_->mkdir("single-resolution-pi0cluster");
  singleClusResolutionPi0Clusters->cd();

  dtUpToQuarterGeVEBPeak_ -> Write(); 
  dtUpToHalfGeVEBPeak_    -> Write(); 
  dtUpToOneGeVEBPeak_     -> Write(); 
  dtUpToTwoGeVEBPeak_     -> Write(); 
  dtUpOverTwoGeVEBPeak_  -> Write(); 
				    
  dtUpToThreeQuarterGeVEEPeak_ -> Write(); 
  dtUpToOneAndHalfGeVEEPeak_   -> Write(); 
  dtUpToThreeGeVEEPeak_        -> Write(); 
  dtUpToSixGeVEEPeak_          -> Write(); 
  dtUpOverSixGeVEEPeak_        -> Write(); 
  
  dtVSAeffHistAnyPeak_ -> Write(); 
  dtVSAeffHistEBPeak_  -> Write(); 
  dtVSAeffHistEEPeak_  -> Write(); 
  dtVSAeffProfAnyPeak_ -> Write(); 
  dtVSAeffProfEBPeak_  -> Write(); 
  dtVSAeffProfEEPeak_  -> Write(); 
  
  dtRMSVSAeffAnyPeak_-> Write();
  dtSigmaAeffAnyPeak_-> Write();
  singleClusterChi2HistAnyPeak_->Write();
  singleClusterChi2NDFHistAnyPeak_->Write();
  singleClusterChi2HistEBPeak_->Write();
  singleClusterChi2NDFHistEBPeak_->Write();
  singleClusterChi2HistEEPeak_->Write();
  singleClusterChi2NDFHistEEPeak_->Write();
  dtPullSingleClusterHistPi0Peak_->Write();
  dtPullSingleClusterHistPi0PeakEB_->Write();
  dtPullSingleClusterHistPi0PeakEE_->Write();

  dtSigmaAeffEBPeak_ -> Write();
  dtSigmaAeffEEPeak_ -> Write();

  // write out 1-d control plots for DeltaT RMS and sigma for any peak
  TDirectory *singleClusResolutionSlicesPeak = singleClusResolutionPi0Clusters->mkdir("dtslices-anyPeak");
  singleClusResolutionSlicesPeak->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSliceVSAeffAnyPeak_[v] -> Write();  }

  // write out 1-d control plots for DeltaT RMS and sigma for EB peak
  TDirectory *singleClusResolutionEBSlicesPeak = singleClusResolutionPi0Clusters->mkdir("dtslices-EBPeak");
  singleClusResolutionEBSlicesPeak->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSliceVSAeffEBPeak_[v] -> Write();  }

  // write out 1-d control plots for DeltaT RMS and sigma for EE peak
  TDirectory *singleClusResolutionEESlicesPeak = singleClusResolutionPi0Clusters->mkdir("dtslices-EEPeak");
  singleClusResolutionEESlicesPeak->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSliceVSAeffEEPeak_[v] -> Write();  }
  
  // write out double cluster resolution plots
  TDirectory *doubleClusResolution = saving_->mkdir("double-resolution");
  doubleClusResolution->cd();

  dtDoubleClusterHistAny_    ->Write();
  dtPullDoubleClusterHistAny_->Write();
  dtPullChi2CutDoubleClusterHistAny_->Write();
  dtVsPtDoubleClusterHistAny_->Write();

  dtDoubleClusterHistPi0Peak_->Write();
  dtDoubleClusterHistPi0PeakEE_->Write();
  dtDoubleClusterHistPi0PeakEB_->Write();
  dtPullDoubleClusterHistPi0Peak_  ->Write();
  dtPullDoubleClusterHistPi0PeakEE_->Write();
  dtPullDoubleClusterHistPi0PeakEB_->Write();
  dtPullChi2CutDoubleClusterHistPi0Peak_  ->Write();
  dtPullChi2CutDoubleClusterHistPi0PeakEE_->Write();
  dtPullChi2CutDoubleClusterHistPi0PeakEB_->Write();
  dtVsPtDoubleClusterHistPi0Peak_  ->Write();
  dtVsPtDoubleClusterHistPi0PeakEE_->Write();
  dtVsPtDoubleClusterHistPi0PeakEB_->Write();

  dtSeedsVsAeffDoubleClusterHist_  ->Write();
  dtSeedsVsAeffDoubleClusterHistEE_->Write();
  dtSeedsVsAeffDoubleClusterHistEB_->Write();
  dtVsAeffDoubleClusterHist_  ->Write();
  dtVsAeffDoubleClusterHistEE_->Write();
  dtVsAeffDoubleClusterHistEB_->Write();
  dtSeedsSigmaAeffAny_->Write();
  dtSeedsSigmaAeffEE_->Write();
  dtSeedsSigmaAeffEB_->Write();
  dtClustersSigmaAeffAny_->Write();
  dtClustersSigmaAeffEE_->Write();
  dtClustersSigmaAeffEB_->Write();
 
  // write out 1-d control plots for DeltaT cluster-SEEDS RMS and sigma
  TDirectory *doubleClusSeedsResolutionSlices = doubleClusResolution->mkdir("dtSeedsSlices-Any");
  doubleClusSeedsResolutionSlices->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSeedsSlicesVsAeffDoubleClusterHist_[v] -> Write();  }

  TDirectory *doubleClusSeedsResolutionSlicesEE = doubleClusResolution->mkdir("dtSeedsSlices-EE");
  doubleClusSeedsResolutionSlicesEE->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSeedsSlicesVsAeffDoubleClusterHistEE_[v] -> Write();  }

  TDirectory *doubleClusSeedsResolutionSlicesEB = doubleClusResolution->mkdir("dtSeedsSlices-EB");
  doubleClusSeedsResolutionSlicesEB->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSeedsSlicesVsAeffDoubleClusterHistEB_[v] -> Write();  }


  // write out 1-d control plots for DeltaT cluster (combining cry info) RMS and sigma
  TDirectory *doubleClusSeedsResolution = doubleClusResolution->mkdir("dtClustersSlices-Any");
  doubleClusSeedsResolution->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSlicesVsAeffDoubleClusterHist_[v] -> Write();  }

  TDirectory *doubleClusSeedsResolutionEE = doubleClusResolution->mkdir("dtClustersSlices-EE");
  doubleClusSeedsResolutionEE->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSlicesVsAeffDoubleClusterHistEE_[v] -> Write();  }

  TDirectory *doubleClusSeedsResolutionEB = doubleClusResolution->mkdir("dtClustersSlices-EB");
  doubleClusSeedsResolutionEB->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSlicesVsAeffDoubleClusterHistEB_[v] -> Write();  }




  // write out double cluster pi0peak resolution plots
  TDirectory *doubleClusResolutionPi0Clusters = saving_->mkdir("double-resolution-pi0");
  doubleClusResolutionPi0Clusters->cd();
  dtSeedsVsAeffDoubleClusterHistPi0Peak_  ->Write();
  dtSeedsVsAeffDoubleClusterHistPi0PeakEE_->Write();
  dtSeedsVsAeffDoubleClusterHistPi0PeakEB_->Write();
  dtVsAeffDoubleClusterHistPi0Peak_  ->Write();
  dtVsAeffDoubleClusterHistPi0PeakEE_->Write();
  dtVsAeffDoubleClusterHistPi0PeakEB_->Write();

  dtSeedsSigmaAeffPi0Peak_->Write();
  dtSeedsSigmaAeffPi0PeakEE_->Write();
  dtSeedsSigmaAeffPi0PeakEB_->Write();
  dtClustersSigmaAeffPi0Peak_->Write();
  dtClustersSigmaAeffPi0PeakEE_->Write();
  dtClustersSigmaAeffPi0PeakEB_->Write();

  // write out 1-d control plots for DeltaT cluster-SEEDS RMS and sigma under the peak 
  TDirectory *doubleClusSeedsResolutionSlicesPeak = doubleClusResolutionPi0Clusters->mkdir("dtSeedsSlices-Any");
  doubleClusSeedsResolutionSlicesPeak->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSeedsSlicesVsAeffDoubleClusterHistPi0Peak_[v] -> Write();  }

  TDirectory *doubleClusSeedsResolutionSlicesEEPeak = doubleClusResolutionPi0Clusters->mkdir("dtSeedsSlices-EE");
  doubleClusSeedsResolutionSlicesEEPeak->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEE_[v] -> Write();  }

  TDirectory *doubleClusSeedsResolutionSlicesEBPeak = doubleClusResolutionPi0Clusters->mkdir("dtSeedsSlices-EB");
  doubleClusSeedsResolutionSlicesEBPeak->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEB_[v] -> Write();  }


  // write out 1-d control plots for DeltaT cluster (combining cry info) RMS and sigma under the peak 
  TDirectory *doubleClusSeedsResolutionPeak = doubleClusResolutionPi0Clusters->mkdir("dtClustersSlices-Any");
  doubleClusSeedsResolutionPeak->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSlicesVsAeffDoubleClusterHistPi0Peak_[v] -> Write();  }

  TDirectory *doubleClusSeedsResolutionEEPeak = doubleClusResolutionPi0Clusters->mkdir("dtClustersSlices-EE");
  doubleClusSeedsResolutionEEPeak->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSlicesVsAeffDoubleClusterHistPi0PeakEE_[v] -> Write();  }

  TDirectory *doubleClusSeedsResolutionEBPeak = doubleClusResolutionPi0Clusters->mkdir("dtClustersSlices-EB");
  doubleClusSeedsResolutionEBPeak->cd();
  for (int v=0; v<AeffNBins_; v++){    dtSlicesVsAeffDoubleClusterHistPi0PeakEB_[v] -> Write();  }




}

// ---------------------------------------------------------------------------------------
// ------------------ Function to compute time and error for a cluster -------------------
//std::pair<float,float> timeAndUncertSingleCluster(int bClusterIndex)
ClusterTime timeAndUncertSingleCluster(int bClusterIndex)
{
  float weightTsum  = 0;
  float weightSum   = 0;
  int   numCrystals = 0;

  // loop on the cry components of a basic cluster; get timeBest and uncertainty 
  for(int thisCry=0; thisCry<treeVars_.nXtalsInCluster[bClusterIndex]; thisCry++)
  {
    bool  thisIsInEB=false;
    float sigmaNoiseOfThis=0;
    if(treeVars_.xtalInBCIEta[bClusterIndex][thisCry]!=-999999)       {
      sigmaNoiseOfThis=sigmaNoiseEB;
      thisIsInEB=true;    }
    else if(treeVars_.xtalInBCIy[bClusterIndex][thisCry]!=-999999)    {
      sigmaNoiseOfThis=sigmaNoiseEE;
      thisIsInEB=false;    }
    else    {  std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}
    float ampliOverSigOfThis = treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] / sigmaNoiseOfThis; 
    if( ampliOverSigOfThis < minAmpliOverSigma_) continue;

    numCrystals++;
    float timeOfThis  = treeVars_.xtalInBCTime[bClusterIndex][thisCry];
    float sigmaOfThis = sqrt(pow(timingResParamN/ampliOverSigOfThis,2)+pow(timingResParamConst,2));

    //std::cout << "GFdeb eampli: " << treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] //gfdebug
    //          << " ampliOverSigOfThis: " << ampliOverSigOfThis
    //          << " timeOfThis: " << timeOfThis
    //          << " sigmaOfThis: " << sigmaOfThis
    //          << std::endl;//gfdebug

    weightTsum+=(timeOfThis/pow(sigmaOfThis,2));
    weightSum+=1/pow(sigmaOfThis,2);
  }
  float bestTime = weightTsum/weightSum;

  float chi2 = -999999;
  // loop on the cry components to get chi2
  // do this only if you have at least 2 crystals over threshold
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


  ClusterTime theResult; //initialize
  theResult.numCry = -999999;   theResult.time   = -999999;
  theResult.timeErr= -999999;   theResult.chi2   = -999999;
  
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
// ------------------ Function to do single BasicCluster resolution studies  -------------
void doSingleClusterResolutionPlots(std::set<int> bcIndicies, bool isAfterPi0Selection)
{
  // single-cluster resolution studies: DeltaT VS Aeff
  // filling different plots for all clusters compared to clusters from pi0 peak selection

  for (std::set<int>::const_iterator bcItr = bcIndicies.begin();
      bcItr!=bcIndicies.end(); ++bcItr)
  {// loop on bc

    int bCluster = *bcItr;
    ClusterTime myClusterTime = timeAndUncertSingleCluster(bCluster);
    bool isEB=false;
    if( fabs(treeVars_.clusterEta[bCluster]) < BarrelLimit) isEB=true;
    if(!isAfterPi0Selection && myClusterTime.numCry > 1)
    {
      singleClusterChi2HistAny_->Fill(myClusterTime.chi2);
      singleClusterChi2NDFHistAny_->Fill(myClusterTime.chi2/(myClusterTime.numCry-1));
      if(isEB)
      {
        singleClusterChi2HistEB_->Fill(myClusterTime.chi2);
        singleClusterChi2NDFHistEB_->Fill(myClusterTime.chi2/(myClusterTime.numCry-1));
      }
      else
      {
        singleClusterChi2HistEE_->Fill(myClusterTime.chi2);
        singleClusterChi2NDFHistEE_->Fill(myClusterTime.chi2/(myClusterTime.numCry-1));
      }
    }
    else if(myClusterTime.numCry > 1) // these are the pi-zero clusters
    {
      singleClusterChi2HistAnyPeak_->Fill(myClusterTime.chi2);
      singleClusterChi2NDFHistAnyPeak_->Fill(myClusterTime.chi2/(myClusterTime.numCry-1));
      if(isEB)
      {
        singleClusterChi2HistEBPeak_->Fill(myClusterTime.chi2);
        singleClusterChi2NDFHistEBPeak_->Fill(myClusterTime.chi2/(myClusterTime.numCry-1));
      }
      else
      {
        singleClusterChi2HistEEPeak_->Fill(myClusterTime.chi2);
        singleClusterChi2NDFHistEEPeak_->Fill(myClusterTime.chi2/(myClusterTime.numCry-1));
      }
    }

    
    // loop on the cry components of a basic cluster
    for(int thisCry=0; thisCry<treeVars_.nXtalsInCluster[bCluster]; thisCry++)
    {
      bool  thisIsInEB=false;
      float sigmaNoiseOfThis=0;
      if(treeVars_.xtalInBCIEta[bCluster][thisCry]      !=-999999)   {sigmaNoiseOfThis=sigmaNoiseEB; thisIsInEB=true;}
      else if (treeVars_.xtalInBCIy[bCluster][thisCry]  !=-999999)   {sigmaNoiseOfThis=sigmaNoiseEE; thisIsInEB=false;}
      else {std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}

      float ampliOverSigOfThis = treeVars_.xtalInBCAmplitudeADC[bCluster][thisCry] / sigmaNoiseOfThis; 
      float ampliOfThis = treeVars_.xtalInBCAmplitudeADC[bCluster][thisCry]; 
      float sigmaOfThis = sqrt(pow(timingResParamN/ampliOverSigOfThis,2)+pow(timingResParamConst,2));
      if( ampliOverSigOfThis < minAmpliOverSigma_) continue;

      // loop on the other cry among the components of a basic cluster
      for(int thatCry=thisCry+1; thatCry<treeVars_.nXtalsInCluster[bCluster]; thatCry++)
      {
        float sigmaNoiseOfThat=0;
        if(treeVars_.xtalInBCIEta[bCluster][thatCry]      !=-999999)   sigmaNoiseOfThat=sigmaNoiseEB;
        else if (treeVars_.xtalInBCIy[bCluster][thatCry]  !=-999999)   sigmaNoiseOfThat=sigmaNoiseEE;
        else {std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}

        float ampliOverSigOfThat = treeVars_.xtalInBCAmplitudeADC[bCluster][thatCry] / sigmaNoiseOfThat; 
        float ampliOfThat = treeVars_.xtalInBCAmplitudeADC[bCluster][thatCry];
        float sigmaOfThat = sqrt(pow(timingResParamN/ampliOverSigOfThis,2)+pow(timingResParamConst,2));

        float Aeff = ampliOfThis * ampliOfThat / sqrt( pow(ampliOfThis,2) + pow(ampliOfThat,2) );
        float dt  = treeVars_.xtalInBCTime[bCluster][thisCry] - treeVars_.xtalInBCTime[bCluster][thatCry]; 
        float errorDt = sqrt( pow(sigmaOfThis,2) + pow(sigmaOfThat,2));

        // Insert deltaT vs. ampli plots here
        if(ampliOverSigOfThis > 30 && ampliOverSigOfThat > 1) // at least require the 2nd cry to have 1 ampli/sigma
        {
          deltaTCrysVsAmplitudeAny_->Fill(ampliOverSigOfThat,-1*dt);
          if(thisIsInEB)
            deltaTCrysVsAmplitudeEB_->Fill(ampliOverSigOfThat,-1*dt);
          else
            deltaTCrysVsAmplitudeEE_->Fill(ampliOverSigOfThat,-1*dt);
        }
        
        // If cry below amp. threshold, skip it
        if( ampliOverSigOfThat < minAmpliOverSigma_) continue;


        // for debug
        //std::cout << "ampliOverSigOfThis: " << ampliOverSigOfThis << "\tampliOverSigOfThat: " << ampliOverSigOfThat
        //          << "\n timeOfThis: " << treeVars_.xtalInBCTime[bCluster][thisCry] << "\ttimeOfThat: " << treeVars_.xtalInBCTime[bCluster][thatCry]
        //          << "\n Aeff: " << Aeff << "\tdt: " << dt 
        //          << std::endl;

        if(!isAfterPi0Selection)
        {
          if (thisIsInEB) {
            if      (Aeff < 6)   dtUpToQuarterGeVEB_->Fill(dt);
            else if (Aeff < 12)  dtUpToHalfGeVEB_   ->Fill(dt);
            else if (Aeff < 24)  dtUpToOneGeVEB_    ->Fill(dt);
            else if (Aeff < 48)  dtUpToTwoGeVEB_    ->Fill(dt);
            else                 dtUpOverTwoGeVEB_  ->Fill(dt);

            dtVSAeffHistEB_ -> Fill(Aeff/sigmaNoiseEB, dt); 
            dtVSAeffProfEB_ -> Fill(Aeff/sigmaNoiseEB, dt); 
	    
	    dtVSAeffHistAny_  -> Fill(Aeff/sigmaNoiseEB, dt);
	    dtVSAeffProfAny_  -> Fill(Aeff/sigmaNoiseEB, dt);

	    int bin = dtVSAeffHistAny_ -> FindBin(Aeff/sigmaNoiseEB,-DtMax_); // finding bin at the minumum Y
	    //int binTMP=bin; 
	    bin-=(AeffNBins_+2);// removing bins of Y underflows
	    bin-=1; // removing bin of X underflow for minumum valid Y 
	    bin = ((bin-1)% AeffNBins_)+1;
	    //std::cout << "Aeff/sigmaNoiseEB: " << Aeff/sigmaNoiseEB << " " << binTMP << " " << bin << std::endl;
	    if(0<=bin && bin<AeffNBins_) AeffSliceAny_[bin]->Fill(Aeff/sigmaNoiseEB);

	    bin = dtVSAeffHistEB_ -> FindBin(Aeff/sigmaNoiseEB,-DtMax_); // finding bin at the minumum Y
	    //binTMP=bin; 
	    bin-=(AeffNBins_+2);// removing bins of Y underflows
	    bin-=1; // removing bin of X underflow for minumum valid Y 
	    bin = ((bin-1)% AeffNBins_)+1;
	    //std::cout << "Aeff/sigmaNoiseEB: " << Aeff/sigmaNoiseEB << " " << binTMP << " " << bin << std::endl;
	    if(0<=bin && bin<AeffNBins_) AeffSliceEB_[bin]->Fill(Aeff/sigmaNoiseEB);
          }
          else      {
            if      (Aeff < 6)   dtUpToThreeQuarterGeVEE_->Fill(dt);
            else if (Aeff < 12)  dtUpToOneAndHalfGeVEE_  ->Fill(dt);
            else if (Aeff < 24)  dtUpToThreeGeVEE_       ->Fill(dt);
            else if (Aeff < 48)  dtUpToSixGeVEE_         ->Fill(dt);
            else                 dtUpOverSixGeVEE_       ->Fill(dt);

            dtVSAeffHistEE_ -> Fill(Aeff/sigmaNoiseEE, dt); 
            dtVSAeffProfEE_ -> Fill(Aeff/sigmaNoiseEE, dt);

	    dtVSAeffHistAny_  -> Fill(Aeff/sigmaNoiseEE, dt);
	    dtVSAeffProfAny_  -> Fill(Aeff/sigmaNoiseEE, dt);

	    int bin = dtVSAeffHistAny_ -> FindBin(Aeff/sigmaNoiseEE,-DtMax_); // finding bin at the minumum Y
	    //int binTMP=bin; 
	    bin-=(AeffNBins_+2);// removing bins of Y underflows
	    bin-=1; // removing bin of X underflow for minumum valid Y 
	    bin = ((bin-1)% AeffNBins_)+1;
	    //std::cout << "Aeff/sigmaNoiseEE: " << Aeff/sigmaNoiseEE << " " << binTMP << " " << bin << std::endl;
	    if(0<=bin && bin<AeffNBins_) AeffSliceAny_[bin]->Fill(Aeff/sigmaNoiseEE);

	    bin = dtVSAeffHistEE_ -> FindBin(Aeff/sigmaNoiseEE,-DtMax_);  // finding bin at the minumum Y
	    //binTMP=bin; 
	    bin-=(AeffNBins_+2);// removing bins of Y underflows
	    bin-=1; // removing bin of X underflow for minumum valid Y 
	    bin = ((bin-1)% AeffNBins_)+1;
	    //std::cout << "Aeff/sigmaNoiseEE: " << Aeff/sigmaNoiseEE << " " << binTMP << " " << bin <<  " " << AeffSliceEE_[bin]->FindBin(Aeff/sigmaNoiseEE) << std::endl;
	    if(0<=bin && bin<AeffNBins_) AeffSliceEE_[bin]->Fill(Aeff/sigmaNoiseEE);
          }
        }
        else // clusters matching the pi0 mass
        {
          dtVSAeffHistAnyPeak_  -> Fill(Aeff, dt); //FIXME: average of sigmaEB/EE?
          dtVSAeffProfAnyPeak_  -> Fill(Aeff, dt); //FIXME: average of sigmaEB/EE?
          dtPullSingleClusterHistPi0Peak_->Fill(dt/errorDt);
          if (thisIsInEB) {
            if      (Aeff < 6)   dtUpToQuarterGeVEBPeak_->Fill(dt);
            else if (Aeff < 12)  dtUpToHalfGeVEBPeak_   ->Fill(dt);
            else if (Aeff < 24)  dtUpToOneGeVEBPeak_    ->Fill(dt);
            else if (Aeff < 48)  dtUpToTwoGeVEBPeak_    ->Fill(dt);
            else                 dtUpOverTwoGeVEBPeak_  ->Fill(dt);

            dtVSAeffHistEBPeak_ -> Fill(Aeff/sigmaNoiseEB, dt); 
            dtVSAeffProfEBPeak_ -> Fill(Aeff/sigmaNoiseEB, dt); 
            dtPullSingleClusterHistPi0PeakEB_->Fill(dt/errorDt);
          }
          else      {
            if      (Aeff < 6)   dtUpToThreeQuarterGeVEEPeak_->Fill(dt);
            else if (Aeff < 12)  dtUpToOneAndHalfGeVEEPeak_  ->Fill(dt);
            else if (Aeff < 24)  dtUpToThreeGeVEEPeak_       ->Fill(dt);
            else if (Aeff < 48)  dtUpToSixGeVEEPeak_         ->Fill(dt);
            else                 dtUpOverSixGeVEEPeak_       ->Fill(dt);

            dtVSAeffHistEEPeak_ -> Fill(Aeff/sigmaNoiseEE, dt); 
            dtVSAeffProfEEPeak_ -> Fill(Aeff/sigmaNoiseEE, dt);
            dtPullSingleClusterHistPi0PeakEE_->Fill(dt/errorDt);
          }
        } // else-if pi0 selection

      }// loop on thatCry

    }// loop on thisCry
  }// end loop on bc
}// end doSingleClusterResolutionPlots


// ---------------------------------------------------------------------------------------
// ------------------ Function to do select pi-zero candidates ---------------------------
SetOfIntPairs selectPi0Candidates()
{
  SetOfIntPairs returnPairs;
  float eTA, eTB ;
  float e22A, e33A,    e22B, e33B;
  float eTGammaMinA,   eTGammaMinB;
  float s4s9GammaMinA, s4s9GammaMinB;
  bool  AisEB,         BisEB;
  float eTPi0Min;

  // (FIRST) loop on basic cluster - to build pi0 candidates and get the mass
  for (int bClusterA=0; bClusterA < treeVars_.nClusters; bClusterA++)
    {
      eTA = treeVars_.clusterTransverseEnergy[bClusterA];
      
      e22A = treeVars_.clusterE2x2[bClusterA];
      e33A = treeVars_.clusterE3x3[bClusterA];
      
      // discriminate between EE and EB and set thresholds accordingly
      if ( fabs(treeVars_.clusterEta[bClusterA]) < BarrelLimit) {
        AisEB         = true;
        eTGammaMinA   = eTGammaMinEB_;
        s4s9GammaMinA = s4s9GammaMinEB_;
      }
      else{
        AisEB         = false;
        eTGammaMinA   = eTGammaMinEE_;
        s4s9GammaMinA = s4s9GammaMinEE_;
      }

      if(treeVars_.clusterEta[bClusterA]<-BarrelLimit)     {
        BCClusterShapeEEHist_  -> Fill(e22A/e33A);
        BCClusterShapeEEMHist_ -> Fill(e22A/e33A);}
      else if(treeVars_.clusterEta[bClusterA]>BarrelLimit) {
        BCClusterShapeEEHist_  -> Fill(e22A/e33A);
        BCClusterShapeEEPHist_ -> Fill(e22A/e33A);}
      else	                                      
        BCClusterShapeEBHist_  -> Fill(e22A/e33A);

      // first selecton cut: photon candidate Et
      if( eTA < eTGammaMinA ) continue;

      // second selection cut: cluster shape
      if ( e22A/e33A < s4s9GammaMinA ) continue;

      for (int bClusterB=(bClusterA+1); bClusterB < treeVars_.nClusters; bClusterB++)
      {

        eTB = treeVars_.clusterTransverseEnergy[bClusterB];

        e22B = treeVars_.clusterE2x2[bClusterB];
        e33B = treeVars_.clusterE3x3[bClusterB];

        // discriminate between EE and EB and set thresholds accordingly
        if ( fabs(treeVars_.clusterEta[bClusterB]) < BarrelLimit) {
          BisEB         = true;
          eTGammaMinB   = eTGammaMinEB_;
          s4s9GammaMinB = s4s9GammaMinEB_;
        }
        else{
          BisEB         = false;
          eTGammaMinB   = eTGammaMinEE_;
          s4s9GammaMinB = s4s9GammaMinEE_;
        }

        // first selecton cut: photon candidate Et
        if( eTB < eTGammaMinB ) continue;

        // second selection cut: cluster shape
        if ( e22B/e33B < s4s9GammaMinB ) continue;

        math::PtEtaPhiMLorentzVectorD gammaA (eTA, treeVars_.clusterEta[bClusterA], treeVars_.clusterPhi[bClusterA], 0);
        math::PtEtaPhiMLorentzVectorD gammaB (eTB, treeVars_.clusterEta[bClusterB], treeVars_.clusterPhi[bClusterB], 0);

        math::PtEtaPhiMLorentzVectorD pi0Candidate = gammaA + gammaB;

        // std::cout << "gammaA: " << gammaA << " " << gammaA.M() << "\t\t gammaB: " << gammaB << " " << gammaB.M() << std::endl;
        // std::cout << "pi0Candidate: " << pi0Candidate << " " << pi0Candidate.M() << std::endl;

        if ( fabs(pi0Candidate.Eta()) < BarrelLimit) {
          eTPi0Min      = eTPi0MinEB_;}
        else{
          eTPi0Min      = eTPi0MinEE_;}

        // third selection cut: pi0 candidate Et
        if(pi0Candidate.Et() < eTPi0Min ) continue;

        massDiGammaHist_ -> Fill(pi0Candidate.M());

        if(treeVars_.clusterEta[bClusterA]<-BarrelLimit){
          massDiGammaEEHist_  -> Fill(pi0Candidate.M());
          massDiGammaEEMHist_ -> Fill(pi0Candidate.M());}
        else if(treeVars_.clusterEta[bClusterA]>BarrelLimit) {
          massDiGammaEEHist_  -> Fill(pi0Candidate.M());
          massDiGammaEEPHist_ -> Fill(pi0Candidate.M());}
        else	     
          massDiGammaEBHist_  -> Fill(pi0Candidate.M());

        // occupancy of all candidates (this is FIRST loop)
        diPhotonOccupancyAny_ -> Fill(pi0Candidate.Eta(), pi0Candidate.Phi());
	
	// TODO: don't you want to insert the mass cut here? //gf

	
	
	
        /////////////////////////////////////////////////////////////
	  // here I have di-gamma pairs that pass cuts
	  // now select pi0's based on the mass
	  /////////////////////////////////////////////////////////////
	  
	  float nSigma =1;
	  // reject sidebands in EB
	  //std::cout << "select pi0 BisEB: " << BisEB << " pi0MassEB_ " << pi0MassEB_ << " pi0WidthEB_ " << pi0WidthEB_ << " mass: " << pi0Candidate.M() << std::endl; //gfdebu
	  if( BisEB &&
	      (pi0Candidate.M() < pi0MassEB_-nSigma*pi0WidthEB_ ||
	       pi0MassEB_+nSigma*pi0WidthEB_ < pi0Candidate.M())
	      ) {
	    diPhotonSidesOccupancyAny_ -> Fill(pi0Candidate.Eta(), pi0Candidate.Phi());
	    //std::cout << "sideband found in EB" << std::endl; //gf debug
	    continue;}
	  
          // reject sidebands in EE
          if( (!BisEB) &&
              (pi0Candidate.M() < pi0MassEE_-nSigma*pi0WidthEE_ ||
               pi0MassEE_+nSigma*pi0WidthEE_ < pi0Candidate.M())
	      ) {
            diPhotonSidesOccupancyAny_ -> Fill(pi0Candidate.Eta(), pi0Candidate.Phi());
	    //std::cout << "sideband found in EE" << std::endl; //gf debug
            continue;}
	  
	  
	  

        returnPairs.insert(std::make_pair<int,int>(bClusterA,bClusterB));

      }//loop on candidateB

    }//loop on candidateA - (FIRST) to build pi0 candidates and get the mass
  
  return returnPairs;
}

// ---------------------------------------------------------------------------------------
// ------------------ Function to fit mass spectra ---------------------------------------
void fitMassSpectra()
{
  // get mass and width for EB 
  TF1 *massEB = new TF1("massEB","gaus(0)+pol3(3)",0,1); 
  //massEB->SetParameter(0,500);
  massEB->SetParameter(1,0.12);
  massEB->SetParameter(2,0.015);  
  massDiGammaEBHist_ ->Fit("massEB",fitOption_.c_str(),"",0,1);
  pi0MassEB_  = massEB->GetParameter(1);
  pi0WidthEB_ = massEB->GetParameter(2);

  // get mass and width for EE 
  TF1 *massEE = new TF1("massEE","gaus(0)+pol3(3)",0,1); 
  massEE->SetParameter(0,500);
  massEE->SetParameter(1,0.12);
  massEE->SetParameter(2,0.02);  
  massDiGammaEEHist_ ->Fit("massEE",fitOption_.c_str(),"",0,1);
  pi0MassEE_  = massEE->GetParameter(1);
  pi0WidthEE_ = massEE->GetParameter(2);

  // get mass and width for EP 
  TF1 *massEEP = new TF1("massEEP","gaus(0)+pol3(3)",0,1); 
  massEEP->SetParameter(0,500);
  massEEP->SetParameter(1,0.12);
  massEEP->SetParameter(2,0.02);  
  massDiGammaEEPHist_ ->Fit("massEEP",fitOption_.c_str(),"",0,1);
  pi0MassEEP_  = massEEP->GetParameter(1);
  pi0WidthEEP_ = massEEP->GetParameter(2);
  
  // get mass and width for EM 
  TF1 *massEEM = new TF1("massEEM","gaus(0)+pol3(3)",0,1); 
  massEEM->SetParameter(0,500);
  massEEM->SetParameter(1,0.12);
  massEEM->SetParameter(2,0.02);  
  massDiGammaEEMHist_ ->Fit("massEEM",fitOption_.c_str(),"",0,1);
  pi0MassEEM_  = massEEM->GetParameter(1);
  pi0WidthEEM_ = massEEM->GetParameter(2);
}

// ---------------------------------------------------------------------------------------
// ------------------ Function to do double cluster resolution studies -------------------
void doDoubleClusterResolutionPlots(SetOfIntPairs myBCpairs, bool isAfterPi0Selection)
{
  // loop on pairs of basicClusters
  for(SetOfIntPairs::const_iterator pairItr = myBCpairs.begin(); pairItr != myBCpairs.end(); ++pairItr)
    {
      int bClusterA = pairItr->first;
      int bClusterB = pairItr->second;

      // find seed of either clusters
      int seedA=0; float ampliSeedA=-1e9;
      for(int cryInBC=0; cryInBC<treeVars_.nXtalsInCluster[bClusterA]; cryInBC++){
	if(treeVars_.xtalInBCAmplitudeADC[bClusterA][cryInBC]>ampliSeedA){
	  ampliSeedA = treeVars_.xtalInBCAmplitudeADC[bClusterA][cryInBC];	  seedA      = cryInBC;	}// end if
      }// loop on crystals of bClusterA
      int seedB=0; float ampliSeedB=-1e9;
      for(int cryInBC=0; cryInBC<treeVars_.nXtalsInCluster[bClusterB]; cryInBC++){
	if(treeVars_.xtalInBCAmplitudeADC[bClusterB][cryInBC]>ampliSeedB){
	  ampliSeedB = treeVars_.xtalInBCAmplitudeADC[bClusterB][cryInBC];	  seedB      = cryInBC;	}// end if
      }// loop on crystals of bClusterB
      
      // calculate time difference of cluster seeds
      float DtSeeds = treeVars_.xtalInBCTime[bClusterA][seedA] - treeVars_.xtalInBCTime[bClusterB][seedB];

      // calculate Aeff for the cluster seeds pair
      float aEffSeeds = ampliSeedA * ampliSeedB / sqrt(ampliSeedA*ampliSeedA + ampliSeedB*ampliSeedB);
      
      float eTA = treeVars_.clusterTransverseEnergy[bClusterA];
      float eTB = treeVars_.clusterTransverseEnergy[bClusterB];
      
      // re-build the pi0 candidate
      math::PtEtaPhiMLorentzVectorD gammaA (eTA, treeVars_.clusterEta[bClusterA], treeVars_.clusterPhi[bClusterA], 0);
      math::PtEtaPhiMLorentzVectorD gammaB (eTB, treeVars_.clusterEta[bClusterB], treeVars_.clusterPhi[bClusterB], 0);
      math::PtEtaPhiMLorentzVectorD pi0Candidate = gammaA + gammaB;
      
      bool isEB=false;
      if( fabs(pi0Candidate.Eta()) < BarrelLimit) isEB=true;
      //std::cout << "eta: " << pi0Candidate.Eta() << " isEB " << isEB << std::endl;//gfcomm

      // Make the time check between two clusters in the peaks
      //std::pair<float,float> timeAndUncertClusterA = timeAndUncertSingleCluster(bClusterA);
      ClusterTime timeAndUncertClusterA = timeAndUncertSingleCluster(bClusterA);
      //  Protect against clusters having zero crystals above amplitude threshold
      //  (which will cause the timeErr, time, etc. to be -999999)
      if(timeAndUncertClusterA.timeErr <= 0) // if something went wrong combining the times, bail out
	continue;

      //std::pair<float,float> timeAndUncertClusterB = timeAndUncertSingleCluster(bClusterB);
      ClusterTime timeAndUncertClusterB = timeAndUncertSingleCluster(bClusterB);
      if(timeAndUncertClusterB.timeErr <= 0) // if something went wrong combining the times, bail out
	continue;
      
      float Dt      = timeAndUncertClusterB.time-timeAndUncertClusterA.time;
      float errorDt = sqrt( pow(timeAndUncertClusterA.timeErr,2) + pow(timeAndUncertClusterB.timeErr,2));
      float Pt      = pi0Candidate.Et();
      
      //std::cout << "--A time: " << timeAndUncertClusterA.time
      //  	<< " timeErr: " << timeAndUncertClusterA.timeErr
      //  	<< " timeNcry: "<< timeAndUncertClusterA.numCry
      //  	<< " timechi2: "<< timeAndUncertClusterA.chi2
      //  	<< "\n--B time: " << timeAndUncertClusterB.time
      //  	<< " timeErr: " << timeAndUncertClusterB.timeErr
      //  	<< " timeNcry: "<< timeAndUncertClusterB.numCry
      //  	<< " timechi2: "<< timeAndUncertClusterB.chi2
      //  	<< std::endl;

      if(timeAndUncertClusterA.numCry < 1 || timeAndUncertClusterB.numCry < 1)  // need at least one cry above amp threshold
        continue;

      float chi2NormA = timeAndUncertClusterA.chi2/(timeAndUncertClusterA.numCry-1);
      float chi2NormB = timeAndUncertClusterB.chi2/(timeAndUncertClusterB.numCry-1);
      if(isAfterPi0Selection)
	{
	  //////////////////////////////////////////////////////////
	  // from here on I have pi0 candidates
	  // bClusterA and bClusterB are the two clusters making this candidate
	  diPhotonPeakOccupancyAny_  -> Fill(pi0Candidate.Eta(), pi0Candidate.Phi());
	  
	  dtDoubleClusterHistPi0Peak_            ->Fill(Dt);
	  if(isEB) dtDoubleClusterHistPi0PeakEB_ ->Fill(Dt);
	  else     dtDoubleClusterHistPi0PeakEE_ ->Fill(Dt);

	  dtPullDoubleClusterHistPi0Peak_            ->Fill(Dt/errorDt);
          if(chi2NormA < maxChi2NDF_ && chi2NormB < maxChi2NDF_ && chi2NormA > 0 && chi2NormB > 0)
            dtPullChi2CutDoubleClusterHistPi0Peak_            ->Fill(Dt/errorDt);

	  if(isEB)
          {
            dtPullDoubleClusterHistPi0PeakEB_ ->Fill(Dt/errorDt);
            if(chi2NormA < maxChi2NDF_ && chi2NormB < maxChi2NDF_ && chi2NormA > 0 && chi2NormB > 0)
            dtPullChi2CutDoubleClusterHistPi0PeakEB_ ->Fill(Dt/errorDt);
          }
	  else
          {
            dtPullDoubleClusterHistPi0PeakEE_ ->Fill(Dt/errorDt);
            if(chi2NormA < maxChi2NDF_ && chi2NormB < maxChi2NDF_ && chi2NormA > 0 && chi2NormB > 0)
              dtPullChi2CutDoubleClusterHistPi0PeakEE_ ->Fill(Dt/errorDt);
          }

	  dtVsPtDoubleClusterHistPi0Peak_            ->Fill(Pt,Dt);
	  if(isEB) dtVsPtDoubleClusterHistPi0PeakEB_ ->Fill(Pt,Dt);
	  else     dtVsPtDoubleClusterHistPi0PeakEE_ ->Fill(Pt,Dt);

	  float sigma=sigmaNoiseEE; 
	  if (isEB) sigma=sigmaNoiseEB;
	  dtSeedsVsAeffDoubleClusterHistPi0Peak_ -> Fill(aEffSeeds/sigma, DtSeeds);
	  dtVsAeffDoubleClusterHistPi0Peak_      -> Fill(aEffSeeds/sigma, Dt);
	  if(isEB) 	  dtSeedsVsAeffDoubleClusterHistPi0PeakEB_ -> Fill(aEffSeeds/sigma, DtSeeds);
	  else 	          dtSeedsVsAeffDoubleClusterHistPi0PeakEE_ -> Fill(aEffSeeds/sigma, DtSeeds);
	  if(isEB) 	  dtVsAeffDoubleClusterHistPi0PeakEB_      -> Fill(aEffSeeds/sigma, Dt);
	  else 	          dtVsAeffDoubleClusterHistPi0PeakEE_      -> Fill(aEffSeeds/sigma, Dt);

	} // isAfterPi0Selection
      else
	{// these are any cluster, not pi0 selected
	  dtDoubleClusterHistAny_     ->Fill(Dt);
	  dtPullDoubleClusterHistAny_ ->Fill(Dt/errorDt);
	  dtVsPtDoubleClusterHistAny_ ->Fill(Pt,Dt);
	  
	  float sigma=sigmaNoiseEE; 
	  if (isEB) sigma=sigmaNoiseEB;
	  dtSeedsVsAeffDoubleClusterHist_ -> Fill(aEffSeeds/sigma, DtSeeds);
	  dtVsAeffDoubleClusterHist_      -> Fill(aEffSeeds/sigma, Dt);
	  if(isEB) 	  dtSeedsVsAeffDoubleClusterHistEB_ -> Fill(aEffSeeds/sigma, DtSeeds);
	  else 	          dtSeedsVsAeffDoubleClusterHistEE_ -> Fill(aEffSeeds/sigma, DtSeeds);
	  if(isEB) 	  dtVsAeffDoubleClusterHistEB_      -> Fill(aEffSeeds/sigma, Dt);
	  else 	          dtVsAeffDoubleClusterHistEE_      -> Fill(aEffSeeds/sigma, Dt);

          if(chi2NormA < maxChi2NDF_ && chi2NormB < maxChi2NDF_ && chi2NormA > 0 && chi2NormB > 0)
            dtPullChi2CutDoubleClusterHistAny_ ->Fill(Dt/errorDt);

	}// !isAfterPi0Selection
      
    }//loop on cluster pairs
}// end doDoubleClusterResolutionPlots


// ---------------------------------------------------------------------------------------
// ------------------ Function to do slicing and fitting of plots ------------------------
void doFinalPlots()
{
  for (int sliceX=0; sliceX<AeffNBins_; sliceX++)  {//looping on the X axis, at constant Aeff
    for (int binY=0; binY<numDtBins_; binY++)  {// looping in Delta t bins
      dtSliceVSAeffAny_[sliceX]   ->SetBinContent( (binY+1), (dtVSAeffHistAny_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSliceVSAeffEB_[sliceX]    ->SetBinContent( (binY+1), (dtVSAeffHistEB_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSliceVSAeffEE_[sliceX]    ->SetBinContent( (binY+1), (dtVSAeffHistEE_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSliceVSAeffAnyPeak_[sliceX] ->SetBinContent( (binY+1), (dtVSAeffHistAnyPeak_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSliceVSAeffEBPeak_[sliceX] ->SetBinContent( (binY+1), (dtVSAeffHistEBPeak_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSliceVSAeffEEPeak_[sliceX] ->SetBinContent( (binY+1), (dtVSAeffHistEEPeak_->GetBinContent((sliceX+1),(binY+1))) ); 

      // Dt slices for seeds/clusters in Aeff bins
      dtSeedsSlicesVsAeffDoubleClusterHist_[sliceX]       ->SetBinContent( (binY+1), (dtSeedsVsAeffDoubleClusterHist_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSeedsSlicesVsAeffDoubleClusterHistEE_[sliceX]     ->SetBinContent( (binY+1), (dtSeedsVsAeffDoubleClusterHistEE_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSeedsSlicesVsAeffDoubleClusterHistEB_[sliceX]     ->SetBinContent( (binY+1), (dtSeedsVsAeffDoubleClusterHistEB_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSlicesVsAeffDoubleClusterHist_[sliceX]            ->SetBinContent( (binY+1), (dtVsAeffDoubleClusterHist_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSlicesVsAeffDoubleClusterHistEE_[sliceX]          ->SetBinContent( (binY+1), (dtVsAeffDoubleClusterHistEE_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSlicesVsAeffDoubleClusterHistEB_[sliceX]          ->SetBinContent( (binY+1), (dtVsAeffDoubleClusterHistEB_->GetBinContent((sliceX+1),(binY+1))) ); 

      // Dt slices for seeds/clusters in Aeff bins under pi0 peak
      dtSeedsSlicesVsAeffDoubleClusterHistPi0Peak_[sliceX]       ->SetBinContent( (binY+1), (dtSeedsVsAeffDoubleClusterHistPi0Peak_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEE_[sliceX]     ->SetBinContent( (binY+1), (dtSeedsVsAeffDoubleClusterHistPi0PeakEE_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEB_[sliceX]     ->SetBinContent( (binY+1), (dtSeedsVsAeffDoubleClusterHistPi0PeakEB_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSlicesVsAeffDoubleClusterHistPi0Peak_[sliceX]            ->SetBinContent( (binY+1), (dtVsAeffDoubleClusterHistPi0Peak_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSlicesVsAeffDoubleClusterHistPi0PeakEE_[sliceX]          ->SetBinContent( (binY+1), (dtVsAeffDoubleClusterHistPi0PeakEE_->GetBinContent((sliceX+1),(binY+1))) ); 
      dtSlicesVsAeffDoubleClusterHistPi0PeakEB_[sliceX]          ->SetBinContent( (binY+1), (dtVsAeffDoubleClusterHistPi0PeakEB_->GetBinContent((sliceX+1),(binY+1))) ); 

    }// end loop on Ybins 

    // do slices RMS and fitting for  Any 
    if( dtSliceVSAeffAny_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      // extract RMS and sigma for each Aeff=const slice
      float RMS       = dtSliceVSAeffAny_[sliceX] -> GetRMS();
      float RMSErr    = dtSliceVSAeffAny_[sliceX] -> GetRMSError();
      dtRMSVSAeffAny_ -> SetBinContent(sliceX+1, RMS);
      dtRMSVSAeffAny_ -> SetBinError(sliceX+1, RMSErr);
      
      TF1 *gauss = new TF1("dtFit","gaus",-DtMax_,DtMax_); // require min number entries
      gauss                    ->SetParLimits(1,-5,5); // limit on gaussian central 
      gauss                    ->SetParameter(1,0);    // initialize mean on central value
      //TFitResultPtr Fit(TF1 *f1 ,Option_t *option ,Option_t *goption, Double_t xxmin, Double_t xxmax)
      if (!limitFit_) dtSliceVSAeffAny_[sliceX]->Fit("dtFit",fitOption_.c_str());
      else            dtSliceVSAeffAny_[sliceX]->Fit("dtFit",fitOption_.c_str(),"",-3*RMS,+3*RMS);
      float sigma     = gauss -> GetParameter(2);
      float sigmaErr  = gauss -> GetParError(2);
      if(limitFit_)   dtSliceVSAeffAny_[sliceX]->Fit("dtFit",fitOption_.c_str(),"",-2*sigma,+2*sigma);
      sigma     = gauss -> GetParameter(2);
      sigmaErr  = gauss -> GetParError(2);
      // setting histogram and arrays with sigmas and their errros + abcsissa bins
      dtSigmaAeffAny_ -> SetBinContent(sliceX+1, sigma);    sigmaAeffAny_[sliceX]=sigma;
      dtSigmaAeffAny_ -> SetBinError(sliceX+1, sigmaErr);   sigmaAeffErrAny_[sliceX]=sigmaErr;
      //       AeffBinCentersAny_[sliceX]    = dtSigmaAeffAny_->GetXaxis()->GetBinCenter(sliceX+1);
      //       AeffBinCentersErrAny_[sliceX] = dtSigmaAeffAny_->GetXaxis()->GetBinWidth(sliceX+1)/2;
      AeffBinCentersAny_[sliceX]    = AeffSliceAny_[sliceX]->GetMean();
      AeffBinCentersErrAny_[sliceX] = AeffSliceAny_[sliceX]->GetMeanError();
    }// slices for Any

    // do slices RMS and fitting for EB
    if( dtSliceVSAeffEB_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      // extract RMS and sigma for each Aeff=const slice
      float RMS       = dtSliceVSAeffEB_[sliceX] -> GetRMS();
      float RMSErr    = dtSliceVSAeffEB_[sliceX] -> GetRMSError();
      dtRMSVSAeffEB_ -> SetBinContent(sliceX+1, RMS);
      dtRMSVSAeffEB_ -> SetBinError(sliceX+1, RMSErr);
      
      TF1 *gauss = new TF1("dtFitEB","gaus",-DtMax_,DtMax_); // require min number entries
      gauss                    ->SetParLimits(1,-5,5); // limit on gaussian central 
      gauss                    ->SetParameter(1,0);    // initialize mean on central value
      if (!limitFit_) dtSliceVSAeffEB_[sliceX]->Fit("dtFitEB",fitOption_.c_str());
      else            dtSliceVSAeffEB_[sliceX]->Fit("dtFitEB",fitOption_.c_str(),"",-3*RMS,+3*RMS);
      float sigma     = gauss -> GetParameter(2);
      float sigmaErr  = gauss -> GetParError(2);
      if (limitFit_)   dtSliceVSAeffEB_[sliceX]->Fit("dtFitEB",fitOption_.c_str(),"",-2*sigma,+2*sigma);
      sigma     = gauss -> GetParameter(2);
      sigmaErr  = gauss -> GetParError(2);
      // settnig histogram and arrays with sigmas and their errros + abcsissa bins
      dtSigmaAeffEB_ -> SetBinContent(sliceX+1, sigma);    sigmaAeffEB_[sliceX]=sigma;
      dtSigmaAeffEB_ -> SetBinError(sliceX+1, sigmaErr);   sigmaAeffErrEB_[sliceX]=sigmaErr;
      //       AeffBinCentersEB_[sliceX]    = dtSigmaAeffEB_->GetXaxis()->GetBinCenter(sliceX+1);
      //       AeffBinCentersErrEB_[sliceX] = dtSigmaAeffEB_->GetXaxis()->GetBinWidth(sliceX+1)/2;
      AeffBinCentersEB_[sliceX]    = AeffSliceEB_[sliceX]->GetMean();
      AeffBinCentersErrEB_[sliceX] = AeffSliceEB_[sliceX]->GetMeanError();
    }// slices for EB

    // do slices RMS and fitting for EE
    if( dtSliceVSAeffEE_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      // extract RMS and sigma for each Aeff=const slice
      float RMS       = dtSliceVSAeffEE_[sliceX] -> GetRMS();
      float RMSErr    = dtSliceVSAeffEE_[sliceX] -> GetRMSError();
      dtRMSVSAeffEE_ -> SetBinContent(sliceX+1, RMS);
      dtRMSVSAeffEE_ -> SetBinError(sliceX+1, RMSErr);
      
      TF1 *gauss = new TF1("dtFitEE","gaus",-DtMax_,DtMax_); // require min number entries
      gauss                    ->SetParLimits(1,-5,5); // limit on gaussian central 
      gauss                    ->SetParameter(1,0);    // initialize on central value
      if (!limitFit_) dtSliceVSAeffEE_[sliceX]->Fit("dtFitEE",fitOption_.c_str());
      else            dtSliceVSAeffEE_[sliceX]->Fit("dtFitEE",fitOption_.c_str(),"",-3*RMS,+3*RMS);
      float sigma     = gauss -> GetParameter(2);
      float sigmaErr  = gauss -> GetParError(2);
      if(limitFit_) dtSliceVSAeffEE_[sliceX]->Fit("dtFitEE",fitOption_.c_str(),"",-2*sigma,+2*sigma);
      sigma     = gauss -> GetParameter(2);
      sigmaErr  = gauss -> GetParError(2);
      // setting histogram and arrays with sigmas and their errros + abcsissa bins
      dtSigmaAeffEE_ -> SetBinContent(sliceX+1, sigma);  sigmaAeffEE_[sliceX]=sigma;
      dtSigmaAeffEE_ -> SetBinError(sliceX+1, sigmaErr); sigmaAeffErrEE_[sliceX]=sigmaErr;
      //       AeffBinCentersEE_[sliceX]    = dtSigmaAeffEE_->GetXaxis()->GetBinCenter(sliceX+1);
      //       AeffBinCentersErrEE_[sliceX] = dtSigmaAeffEE_->GetXaxis()->GetBinWidth(sliceX+1)/2;
      AeffBinCentersEE_[sliceX]    = AeffSliceEE_[sliceX]->GetMean();
      AeffBinCentersErrEE_[sliceX] = AeffSliceEE_[sliceX]->GetMeanError();
    }// slices for EE

    
    // **** cluster fits ****    

    // get sigma from fitting Cluster Seeds Any
    if( dtSeedsSlicesVsAeffDoubleClusterHist_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitClus","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSeedsSlicesVsAeffDoubleClusterHist_[sliceX]->Fit("dtFitClus",fitOption_.c_str());
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtSeedsSigmaAeffAny_ -> SetBinContent(sliceX+1, sigma); dtSeedsSigmaAeffAny_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster Seeds Any
    
    // get sigma from fitting Cluster Seeds EE
    if( dtSeedsSlicesVsAeffDoubleClusterHistEE_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitClusEE","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSeedsSlicesVsAeffDoubleClusterHistEE_[sliceX]->Fit("dtFitClusEE",fitOption_.c_str());
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtSeedsSigmaAeffEE_ -> SetBinContent(sliceX+1, sigma); dtSeedsSigmaAeffEE_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster Seeds EE

    // get sigma from fitting Cluster Seeds EB
    if( dtSeedsSlicesVsAeffDoubleClusterHistEB_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitClusEB","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSeedsSlicesVsAeffDoubleClusterHistEB_[sliceX]->Fit("dtFitClusEB",fitOption_.c_str());
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtSeedsSigmaAeffEB_ -> SetBinContent(sliceX+1, sigma); dtSeedsSigmaAeffEB_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster Seeds EB
    


    // get sigma from fitting Cluster (combined crystals) Any
    if( dtSlicesVsAeffDoubleClusterHist_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitClusComb","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSlicesVsAeffDoubleClusterHist_[sliceX]->Fit("dtFitClusComb",fitOption_.c_str());
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtClustersSigmaAeffAny_ -> SetBinContent(sliceX+1, sigma); dtClustersSigmaAeffAny_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster (combined crystals) Any
    
    // get sigma from fitting Cluster (combined crystals) EE
    if( dtSlicesVsAeffDoubleClusterHistEE_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitClusCombEE","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSlicesVsAeffDoubleClusterHistEE_[sliceX]->Fit("dtFitClusCombEE",fitOption_.c_str());
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtClustersSigmaAeffEE_ -> SetBinContent(sliceX+1, sigma); dtClustersSigmaAeffEE_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster (combined crystals) EE

    // get sigma from fitting Cluster (combined crystals) EB
    if( dtSlicesVsAeffDoubleClusterHistEB_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitClusCombEB","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSlicesVsAeffDoubleClusterHistEB_[sliceX]->Fit("dtFitClusCombEB",fitOption_.c_str());
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtClustersSigmaAeffEB_ -> SetBinContent(sliceX+1, sigma); dtClustersSigmaAeffEB_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster (combined crystals) EB
    

    // **** Peak fit for clusters ****
    
    // do slices RMS and fitting for any peak
    if( dtSliceVSAeffAnyPeak_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      // extract RMS and sigma for each Aeff=const slice
      float RMS       = dtSliceVSAeffAnyPeak_[sliceX] -> GetRMS();
      float RMSErr    = dtSliceVSAeffAnyPeak_[sliceX] -> GetRMSError();
      dtRMSVSAeffAnyPeak_ -> SetBinContent(sliceX+1, RMS);
      dtRMSVSAeffAnyPeak_ -> SetBinError(sliceX+1, RMSErr);
      
      TF1 *gauss = new TF1("dtFitPeak","gaus",-DtMax_,DtMax_); // require min number entries
      gauss                    ->SetParLimits(1,-5,5); // limit on gaussian central 
      gauss                    ->SetParameter(1,0);    // initialize on central value
      dtSliceVSAeffAnyPeak_[sliceX]->Fit("dtFitPeak",fitOption_.c_str());
      float sigma     = gauss -> GetParameter(2);
      float sigmaErr  = gauss -> GetParError(2);
      dtSigmaAeffAnyPeak_ -> SetBinContent(sliceX+1, sigma);
      dtSigmaAeffAnyPeak_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for any peak
    
    // do slices RMS and fitting for EB
    if( dtSliceVSAeffEBPeak_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      // extract RMS and sigma for each Aeff=const slice
      float RMS       = dtSliceVSAeffEBPeak_[sliceX] -> GetRMS();
      float RMSErr    = dtSliceVSAeffEBPeak_[sliceX] -> GetRMSError();
      dtRMSVSAeffEBPeak_ -> SetBinContent(sliceX+1, RMS);
      dtRMSVSAeffEBPeak_ -> SetBinError(sliceX+1, RMSErr);
      
      TF1 *gauss = new TF1("dtFitPeakEB","gaus",-DtMax_,DtMax_); // require min number entries
      gauss                    ->SetParLimits(1,-5,5); // limit on gaussian central 
      gauss                    ->SetParameter(1,0);    // initialize on central value
      dtSliceVSAeffEBPeak_[sliceX]->Fit("dtFitPeakEB",fitOption_.c_str());
      float sigma     = gauss -> GetParameter(2);
      float sigmaErr  = gauss -> GetParError(2);
      dtSigmaAeffEBPeak_ -> SetBinContent(sliceX+1, sigma);
      dtSigmaAeffEBPeak_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for EB
    
    // do slices RMS and fitting for EE
    if( dtSliceVSAeffEEPeak_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      // extract RMS and sigma for each Aeff=const slice
      float RMS       = dtSliceVSAeffEEPeak_[sliceX] -> GetRMS();
      float RMSErr    = dtSliceVSAeffEEPeak_[sliceX] -> GetRMSError();
      dtRMSVSAeffEEPeak_ -> SetBinContent(sliceX+1, RMS);
      dtRMSVSAeffEEPeak_ -> SetBinError(sliceX+1, RMSErr);
      
      TF1 *gauss = new TF1("dtFitPeakEE","gaus",-DtMax_,DtMax_); // require min number entries
      gauss                    ->SetParLimits(1,-5,5); // limit on gaussian central 
      gauss                    ->SetParameter(1,0);    // initialize on central value
      dtSliceVSAeffEEPeak_[sliceX]->Fit("dtFitPeakEE",fitOption_.c_str());
      float sigma     = gauss -> GetParameter(2);
      float sigmaErr  = gauss -> GetParError(2);
      dtSigmaAeffEEPeak_ -> SetBinContent(sliceX+1, sigma);
      dtSigmaAeffEEPeak_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for EE


    // get sigma from fitting Cluster Seeds Under the peak Any
    if( dtSeedsSlicesVsAeffDoubleClusterHistPi0Peak_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitSeed","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSeedsSlicesVsAeffDoubleClusterHistPi0Peak_[sliceX]->Fit("dtFitSeed",fitOption_.c_str());
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtSeedsSigmaAeffPi0Peak_ -> SetBinContent(sliceX+1, sigma); dtSeedsSigmaAeffPi0Peak_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster Seeds Any
    
    // get sigma from fitting Cluster Seeds Under the peak EE
    if( dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEE_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitSeedEE","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEE_[sliceX]->Fit("dtFitSeedEE");
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtSeedsSigmaAeffPi0PeakEE_ -> SetBinContent(sliceX+1, sigma); dtSeedsSigmaAeffPi0PeakEE_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster Seeds EE

    // get sigma from fitting Cluster Seeds Under the peak EB
    if( dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEB_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitSeedEB","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSeedsSlicesVsAeffDoubleClusterHistPi0PeakEB_[sliceX]->Fit("dtFitSeedEB");
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtSeedsSigmaAeffPi0PeakEB_ -> SetBinContent(sliceX+1, sigma); dtSeedsSigmaAeffPi0PeakEB_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster Seeds EB
    


    // get sigma from fitting Cluster (combined crystals) Under the peak Any
    if( dtSlicesVsAeffDoubleClusterHistPi0Peak_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitSeedComb","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSlicesVsAeffDoubleClusterHistPi0Peak_[sliceX]->Fit("dtFitSeedComb");
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtClustersSigmaAeffPi0Peak_ -> SetBinContent(sliceX+1, sigma); dtClustersSigmaAeffPi0Peak_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster (combined crystals) Any
    
    // get sigma from fitting Cluster (combined crystals) Under the peak EE
    if( dtSlicesVsAeffDoubleClusterHistPi0PeakEE_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitSeedCombEE","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSlicesVsAeffDoubleClusterHistPi0PeakEE_[sliceX]->Fit("dtFitSeedCombEE");
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtClustersSigmaAeffPi0PeakEE_ -> SetBinContent(sliceX+1, sigma); dtClustersSigmaAeffPi0PeakEE_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster (combined crystals) EE

    // get sigma from fitting Cluster (combined crystals) Under the peak EB
    if( dtSlicesVsAeffDoubleClusterHistPi0PeakEB_[sliceX] -> Integral()  >= minEntriesForFit_ ){
      TF1 *gauss = new TF1("dtFitSeedCombEB","gaus",-DtMax_,DtMax_); // require min number entries
      gauss ->SetParLimits(1,-5,5);   gauss ->SetParameter(1,0);
      dtSlicesVsAeffDoubleClusterHistPi0PeakEB_[sliceX]->Fit("dtFitSeedCombEB");
      float sigma = gauss -> GetParameter(2);                 float sigmaErr  = gauss -> GetParError(2);
      dtClustersSigmaAeffPi0PeakEB_ -> SetBinContent(sliceX+1, sigma); dtClustersSigmaAeffPi0PeakEB_ -> SetBinError(sliceX+1, sigmaErr);
    }// slices for Cluster (combined crystals) EB
    





  }// end loop on Xslices
  
}// end doFinalPlots


// ---------------------------------------------------------------------------------------
// ------------------ Function to make unique set from set of pairs ----------------------
std::set<int> makeUniqueList1D(SetOfIntPairs myPairs)
{
  std::set<int> returnSet;
  // Loop
  for(SetOfIntPairs::const_iterator pairItr = myPairs.begin(); pairItr != myPairs.end();
      ++pairItr)
  {
    returnSet.insert(pairItr->first);
    returnSet.insert(pairItr->second);
  }
  
  return returnSet;
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
  TChain * chain = new TChain ("EcalTimePi0Analysis") ;
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

  // FIXME
  // fit to mass to be made robust
  // re masses to a-priori values for now 
  pi0MassEB_   = 0.111;
  pi0WidthEB_  = 0.013; 
  pi0MassEE_   = 0.126;
  pi0WidthEE_  = 0.030;

  int eventCounter = 0;
  /////////////////////////////////////////////////////
  // Main loop over entries
  for (int entry = 0 ; (entry < nEntries && eventCounter < numEvents_); ++entry)
  {
    chain->GetEntry (entry) ;
    
    // do analysis if the run is in the desired range  
    if( treeVars_.runId<minRun_  || maxRun_<treeVars_.runId) continue;
    
    // do analysis if the LS is in the desired range  
    if( treeVars_.lumiSection<minLS_  || maxLS_<treeVars_.lumiSection) continue;
    
    // if evet being actually processed, increment counter of analyzed events
    eventCounter++;

    speak_=false;
    if (entry<10 || entry%10000==0) speak_=true;

    if (speak_)  std::cout << "\n\n------> reading entry " << entry << " <------\n" ; 
    if (speak_)  std::cout << "  found " << treeVars_.nSuperClusters << " superclusters" << std::endl ;
    if (speak_)  std::cout << "  found " << treeVars_.nClusters << " basic clusters" << std::endl ;
    if (speak_)  std::cout << "  found " << treeVars_.nXtals << " crystals\n" ;    

    // Plot the control hists
    doControlHists();

    // Make pairs of all BasicClusters
    SetOfIntPairs allBCPairs;
    for (int bCluster=0; bCluster < treeVars_.nClusters; bCluster++)
    {
      for (int bClusterA=bCluster+1; bClusterA < treeVars_.nClusters; bClusterA++)
      {
        allBCPairs.insert(std::make_pair<int,int>(bCluster,bClusterA));
      }
    }
    // Do singleCluster plots -- all BC pairs (no pi-zero selection)
    std::set<int> allMyBasicClusterIndicies = makeUniqueList1D(allBCPairs);
    doSingleClusterResolutionPlots(allMyBasicClusterIndicies,false);
    // Do doubleCluster plots -- all BC pairs (no pi-zero selection)
    doDoubleClusterResolutionPlots(allBCPairs,false);

    // ---------------- Select Pi-zeros
    SetOfIntPairs myPi0BasicClusterPairs = selectPi0Candidates();

    // Do double cluster plots
    doDoubleClusterResolutionPlots(myPi0BasicClusterPairs,true);

    // Make unique list of BasicClusters from the list of pairs
    std::set<int> myPi0BasicClusters = makeUniqueList1D(myPi0BasicClusterPairs);
    // Do the singleCluster again on the pi0 BasicClusters
    doSingleClusterResolutionPlots(myPi0BasicClusters,true);


  }   // end of loop over entries
  // now you have di-mass plots filled => get masses


  // Fit the invariant mass spectra
  fitMassSpectra();
  // FIXME
  // fit to mass to be made robust
  // re-set masses to a-priori values for now 
  pi0MassEB_  = 0.111;
  pi0WidthEB_ = 0.013; 
  pi0MassEE_  = 0.126;
  pi0WidthEE_  = 0.030;

  // Do plots that need histograms to be filled with events
  doFinalPlots();

  // now save the plots
  writeHists();
  saving_->Close();

  delete chain ;
  //delete BCNumPerEventHist_;
  //delete BCEnergyHist_;
  //delete BCEtHist_;
  //delete BCNumCrysHist_;
  //delete BCOccupancyHistEB_;
  //delete xtalEnergyHist_; 
  //delete xtalTimeHist_;
  //delete xtalStatusHist_;
  //delete xtalOccupancyHistEB_;
  
  return 0 ;
}
