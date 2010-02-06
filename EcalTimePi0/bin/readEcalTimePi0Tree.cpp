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
#define numAeffBins  250

struct clusterTime {
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
char buffer_ [50];
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
float eTGammaMinEB_   = 0.2;
float s4s9GammaMinEB_ = 0.85;
float eTPi0MinEB_     = 0.65;
float eTGammaMinEE_   = 0.250;
float s4s9GammaMinEE_ = 0.85;
float eTPi0MinEE_     = 0.800;

float minAmpliOverSigma_   = 7;    // dimensionless

// parameters for histograms and ranges
int AeffMax_     =250;
int numDtBins_   =50;
int DtMax_       =10;

// Consts
const float sigmaNoiseEB        = 1.06; // ADC 
const float sigmaNoiseEE        = 2.10; // ADC
const float timingResParamN     = 35.1; // ns
const float timingResParamConst = 0.020; //ns

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
TH2F*   dtVSAeffHistEB_;
TH2F*   dtVSAeffHistEE_;
TProfile* dtVSAeffProfAny_;
TProfile* dtVSAeffProfEB_;
TProfile* dtVSAeffProfEE_;

TH1F*    dtRMSVSAeffAny_;
TH1F*    dtSigmaAeffAny_;

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
TH2F*     dtVSAeffHistEBPeak_;
TH2F*     dtVSAeffHistEEPeak_;
TProfile* dtVSAeffProfAnyPeak_;
TProfile* dtVSAeffProfEBPeak_;
TProfile* dtVSAeffProfEEPeak_;

// double cluster resolution
TH1F* dtDoubleClusterHistAny_;
TH1F* dtPoolDoubleClusterHistAny_;
TH2F* dtVsPtDoubleClusterHistAny_;

TH1F* dtDoubleClusterHistPi0Peak_;
TH1F* dtDoubleClusterHistPi0PeakEE_;
TH1F* dtDoubleClusterHistPi0PeakEB_;
TH1F* dtPoolDoubleClusterHistPi0Peak_;
TH1F* dtPoolDoubleClusterHistPi0PeakEE_;
TH1F* dtPoolDoubleClusterHistPi0PeakEB_;
TH2F* dtVsPtDoubleClusterHistPi0Peak_;
TH2F* dtVsPtDoubleClusterHistPi0PeakEE_;
TH2F* dtVsPtDoubleClusterHistPi0PeakEB_;



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
      std::cout << " --eTGammaMinEB: min eT for EB gammas" << std::endl;
      std::cout << " --s4s9GammaMinEB: min EB shower shape" << std::endl;
      std::cout << " --eTPi0MinEB min eT for EB pi0 candidate" << std::endl;
      std::cout << " --eTGammaMinEE: min eT for EE gammas" << std::endl;
      std::cout << " --s4s9GammaMinEE: min EE shower shape" << std::endl;
      std::cout << " --eTPi0MinEE min eT for EE pi0 candidate" << std::endl;
      std::cout <<  "--minAOverSigma min ampli considered for time" << std::endl;
      std::cout << " --i <list of strings> list of input files" << std::endl ;     
      exit(1);      }


    else if (argv[v] == stringNumEvents) { // set number of events
      std::cout << "events number" << std::endl;
      numEvents_=atoi(argv[v+1]);
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
void initializeHists()
{
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

  dtVSAeffHistAny_ = new TH2F("#Delta(t) VS A_{eff}/#sigma_{N}","#Delta(t) VS A_{eff}/#sigma_{N}",numAeffBins,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);
  dtVSAeffHistEB_  = new TH2F("EB:  #Delta(t)  VS  A_{eff}/#sigma_{N}","EB:  #Delta(t)  VS  A_{eff}/#sigma_{N}",numAeffBins,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);
  dtVSAeffHistEE_  = new TH2F("EE:  #Delta(t)  VS  A_{eff}/#sigma_{N}","EE:  #Delta(t)  VS  A_{eff}/#sigma_{N}",numAeffBins,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);
  dtVSAeffProfAny_ = new TProfile("#Delta(t)  VS  A_{eff}/#sigma_{N} prof","#Delta(t) VS A_{eff}/#sigma_{N} prof",numAeffBins,0.,AeffMax_,-DtMax_,DtMax_);
  for (int v=0; v<numAeffBins; v++){// build histograms for RMS and sigma of DeltaT
    float binLeft=(v*AeffMax_/numAeffBins); float binRight=((v+1)*AeffMax_/numAeffBins);
    sprintf (buffer_, "#Delta t bin %d, [%4.1f,%4.1f]", v+1, binLeft, binRight);
    dtSliceVSAeffAny_[v] = new TH1F(buffer_,buffer_,numDtBins_,-DtMax_,DtMax_);  }
  dtRMSVSAeffAny_  = new TH1F("RMS(#Delta(t)) VS   A_{eff}", "RMS(#Delta(t)) VS   A_{eff}",numAeffBins,0.,AeffMax_);  
  dtSigmaAeffAny_  = new TH1F("#sigma(#Delta(t)) VS   A_{eff}", "#sigma(#Delta(t)) VS   A_{eff}",numAeffBins,0.,AeffMax_);  
  dtVSAeffProfEB_  = new TProfile("EB:  #Delta(t)   VS  A_{eff}/#sigma_{N} prof","EB:  #Delta(t)  VS  A_{eff}/#sigma_{N} prof",numAeffBins,0.,AeffMax_,-DtMax_,DtMax_);
  dtVSAeffProfEE_  = new TProfile("EE:  #Delta(t)   VS  A_{eff}/#sigma_{N} prof","EE:  #Delta(t)  VS  A_{eff}/#sigma_{N} prof",numAeffBins,0.,AeffMax_,-DtMax_,DtMax_);

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

  dtVSAeffHistAnyPeak_ = new TH2F("Peak: #Delta(t) VS A_{eff}/#sigma_{N}","Peak: #Delta(t) VS A_{eff}/#sigma_{N}",numAeffBins,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);
  dtVSAeffHistEBPeak_  = new TH2F("EBPeak: #Delta(t) VS A_{eff}/#sigma_{N}","EBPeak #Delta(t) VS A_{eff}/#sigma_{N}",numAeffBins,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);
  dtVSAeffHistEEPeak_  = new TH2F("EEPeak: E#Delta(t) VS A_{eff}/#sigma_{N}","EEPeak: #Delta(t) VS A_{eff}/#sigma_{N}",numAeffBins,0.,AeffMax_,numDtBins_,-DtMax_,DtMax_);
  dtVSAeffProfAnyPeak_ = new TProfile("Peak: #Delta(t) VS A_{eff}/#sigma_{N} prof","Peak: #Delta(t) VS A_{eff}/#sigma_{N} prof",numAeffBins,0.,AeffMax_,-DtMax_,DtMax_);
  dtVSAeffProfEBPeak_  = new TProfile("EBPeak: #Delta(t) VS A_{eff}/#sigma_{N} prof","EBPeak #Delta(t) VS A_{eff}/#sigma_{N} prof",numAeffBins,0.,AeffMax_,-DtMax_,DtMax_);
  dtVSAeffProfEEPeak_  = new TProfile("EEPeak: #Delta(t) VS A_{eff}/#sigma_{N} prof","EEPeak: #Delta(t) VS A_{eff}/#sigma_{N} prof",numAeffBins,0.,AeffMax_,-DtMax_,DtMax_);

  // should these DeltaT be vary in [-DtMax_, DtMax_] ? Once fixed/understood
  // Initialize histograms -- double cluster resolution
  dtDoubleClusterHistAny_     = new TH1F("DeltaTDoubleClusterAny","#Delta(t) between two clusters EB/EE",100,-25,25);
  dtPoolDoubleClusterHistAny_ = new TH1F("DeltaTPoolDoubleClusterAny","#Delta(t)/#sigma(t) between two clusters EB/EE (pool)",100,-5,5);
  dtVsPtDoubleClusterHistAny_ = new TH2F("DeltaTVSPtDoubleClusterAny","#Delta(t)  between two clusters EB/EE VS P_{t}(di-photon) ",50,0,10,50,-25,25);

  dtDoubleClusterHistPi0Peak_  = new TH1F("DeltaTDoubleClusterPi0Peak","#Delta(t) between two clusters in #pi^{0} mass peak EB/EE",100,-25,25);
  dtDoubleClusterHistPi0PeakEE_= new TH1F("DeltaTDoubleClusterPi0PeakEE","#Delta(t) between two EE clusters in #pi^{0} mass peak",100,-25,25);
  dtDoubleClusterHistPi0PeakEB_= new TH1F("DeltaTDoubleClusterPi0PeakEB","#Delta(t) between two EB clusters in #pi^{0} mass peak",100,-25,25);
  dtPoolDoubleClusterHistPi0Peak_  = new TH1F("DeltaTPoolDoubleClusterPeak","#Delta(t)/#sigma(t) between two clusters EB/EE (pool)",100,-5,5);
  dtPoolDoubleClusterHistPi0PeakEE_= new TH1F("DeltaTPoolDoubleClusterPeakEE","#Delta(t)/#sigma(t) between two clusters EE (pool)",100,-5,5);
  dtPoolDoubleClusterHistPi0PeakEB_= new TH1F("DeltaTPoolDoubleClusterPeakEB","#Delta(t)/#sigma(t) between two clusters EB (pool)",100,-5,5);
  dtVsPtDoubleClusterHistPi0Peak_  = new TH2F("DeltaTVSPtDoubleClusterPeak","#Delta(t) between two clusters EB/EE VS P_{t}(#pi_{0}) ",50,0,10,50,-25,25);
  dtVsPtDoubleClusterHistPi0PeakEE_= new TH2F("DeltaTVSPtDoubleClusterPeakEE","#Delta(t) between two clusters EE VS P_{t}(#pi_{0}) ",50,0,10,50,-25,25);
  dtVsPtDoubleClusterHistPi0PeakEB_= new TH2F("DeltaTVSPtDoubleClusterPeakEB","#Delta(t) between two clusters EB VS P_{t}(#pi_{0}) ",50,0,10,50,-25,25);
}


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

  dtRMSVSAeffAny_-> Write();
  dtSigmaAeffAny_-> Write();

  // write out 1-d control plots for DeltaT RMS and sigma 
  TDirectory *singleClusResolutionSlices = singleClusResolution->mkdir("dtslices-any");
  singleClusResolutionSlices->cd();
  for (int v=0; v<numAeffBins; v++){    dtSliceVSAeffAny_[v] -> Write();  }


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
  
  // write out double cluster resolution plots
  TDirectory *doubleClusResolution = saving_->mkdir("double-resolution");
  doubleClusResolution->cd();

  dtDoubleClusterHistAny_    ->Write();
  dtPoolDoubleClusterHistAny_->Write();
  dtVsPtDoubleClusterHistAny_->Write();

  dtDoubleClusterHistPi0Peak_->Write();
  dtDoubleClusterHistPi0PeakEE_->Write();
  dtDoubleClusterHistPi0PeakEB_->Write();
  dtPoolDoubleClusterHistPi0Peak_  ->Write();
  dtPoolDoubleClusterHistPi0PeakEE_->Write();
  dtPoolDoubleClusterHistPi0PeakEB_->Write();
  dtVsPtDoubleClusterHistPi0Peak_  ->Write();
  dtVsPtDoubleClusterHistPi0PeakEE_->Write();
  dtVsPtDoubleClusterHistPi0PeakEB_->Write();

}

// ---------------------------------------------------------------------------------------
// ------------------ Function to compute time and error for a cluster -------------------
//std::pair<float,float> timeAndUncertSingleCluster(int bClusterIndex)
clusterTime timeAndUncertSingleCluster(int bClusterIndex)
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
    float ampliOfThis = treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] / sigmaNoiseOfThis; 
    if( ampliOfThis < minAmpliOverSigma_) continue;

    numCrystals++;
    float timeOfThis  = treeVars_.xtalInBCTime[bClusterIndex][thisCry];
    float sigmaOfThis = sqrt(pow(timingResParamN/ampliOfThis,2)+pow(timingResParamConst,2));

    //std::cout << "GFdeb eampli: " << treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] //gfdebug
    //          << " ampliOfThis: " << ampliOfThis
    //          << " timeOfThis: " << timeOfThis
    //          << " sigmaOfThis: " << sigmaOfThis
    //          << std::endl;//gfdebug

    weightTsum+=(timeOfThis/pow(sigmaOfThis,2));
    weightSum+=1/pow(sigmaOfThis,2);
  }
  float bestTime = weightTsum/weightSum;

  float chi2 = -1;
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
	
	float ampliOfThis = treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] / sigmaNoiseOfThis; 
	if( ampliOfThis < minAmpliOverSigma_) continue;
	
	float timeOfThis  = treeVars_.xtalInBCTime[bClusterIndex][thisCry];
	float sigmaOfThis = sqrt(pow(timingResParamN/ampliOfThis,2)+pow(timingResParamConst,2));
	
	chi2 += pow( (timeOfThis-bestTime)/sigmaOfThis, 2);
	
      }// end loop on cry
    chi2 /=(numCrystals-1);
  }//end if


  clusterTime theResult; //initialize
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
    // loop on the cry components of a basic cluster
    for(int thisCry=0; thisCry<treeVars_.nXtalsInCluster[bCluster]; thisCry++)
    {
      bool  thisIsInEB=false;
      float sigmaNoiseOfThis=0;
      if(treeVars_.xtalInBCIEta[bCluster][thisCry]      !=-999999)   {sigmaNoiseOfThis=sigmaNoiseEB; thisIsInEB=true;}
      else if (treeVars_.xtalInBCIy[bCluster][thisCry]  !=-999999)   {sigmaNoiseOfThis=sigmaNoiseEE; thisIsInEB=false;}
      else {std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}

      float ampliOfThis = treeVars_.xtalInBCAmplitudeADC[bCluster][thisCry] / sigmaNoiseOfThis; 
      if( ampliOfThis < minAmpliOverSigma_) continue;


      for(int thatCry=thisCry+1; thatCry<treeVars_.nXtalsInCluster[bCluster]; thatCry++)
      {
        float sigmaNoiseOfThat=0;
        if(treeVars_.xtalInBCIEta[bCluster][thatCry]      !=-999999)   sigmaNoiseOfThat=sigmaNoiseEB;
        else if (treeVars_.xtalInBCIy[bCluster][thatCry]  !=-999999)   sigmaNoiseOfThat=sigmaNoiseEE;
        else {std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}

        float ampliOfThat = treeVars_.xtalInBCAmplitudeADC[bCluster][thatCry] / sigmaNoiseOfThat; 
        if( ampliOfThat < minAmpliOverSigma_) continue;

        float Aeff = ampliOfThis * ampliOfThat / sqrt( pow(ampliOfThis,2) + pow(ampliOfThat,2) );
        float dt  = treeVars_.xtalInBCTime[bCluster][thisCry] - treeVars_.xtalInBCTime[bCluster][thatCry]; 

        // for debug
        //std::cout << "ampliOfThis: " << ampliOfThis << "\tampliOfThat: " << ampliOfThat
        //          << "\n timeOfThis: " << treeVars_.xtalInBCTime[bCluster][thisCry] << "\ttimeOfThat: " << treeVars_.xtalInBCTime[bCluster][thatCry]
        //          << "\n Aeff: " << Aeff << "\tdt: " << dt 
        //          << std::endl;

        if(!isAfterPi0Selection)
        {
          dtVSAeffHistAny_  -> Fill(Aeff, dt); 
          dtVSAeffProfAny_  -> Fill(Aeff, dt); 
          if (thisIsInEB) {
            if      (Aeff < 6)   dtUpToQuarterGeVEB_->Fill(dt);
            else if (Aeff < 12)  dtUpToHalfGeVEB_   ->Fill(dt);
            else if (Aeff < 24)  dtUpToOneGeVEB_    ->Fill(dt);
            else if (Aeff < 48)  dtUpToTwoGeVEB_    ->Fill(dt);
            else                 dtUpOverTwoGeVEB_  ->Fill(dt);

            dtVSAeffHistEB_ -> Fill(Aeff, dt); 
            dtVSAeffProfEB_ -> Fill(Aeff, dt); 
          }
          else      {
            if      (Aeff < 6)   dtUpToThreeQuarterGeVEE_->Fill(dt);
            else if (Aeff < 12)  dtUpToOneAndHalfGeVEE_  ->Fill(dt);
            else if (Aeff < 24)  dtUpToThreeGeVEE_       ->Fill(dt);
            else if (Aeff < 48)  dtUpToSixGeVEE_         ->Fill(dt);
            else                 dtUpOverSixGeVEE_       ->Fill(dt);

            dtVSAeffHistEE_ -> Fill(Aeff, dt); 
            dtVSAeffProfEE_ -> Fill(Aeff, dt); }
        }
        else // clusters matching the pi0 mass
        {
          dtVSAeffHistAnyPeak_  -> Fill(Aeff, dt); 
          dtVSAeffProfAnyPeak_  -> Fill(Aeff, dt); 
          if (thisIsInEB) {
            if      (Aeff < 6)   dtUpToQuarterGeVEBPeak_->Fill(dt);
            else if (Aeff < 12)  dtUpToHalfGeVEBPeak_   ->Fill(dt);
            else if (Aeff < 24)  dtUpToOneGeVEBPeak_    ->Fill(dt);
            else if (Aeff < 48)  dtUpToTwoGeVEBPeak_    ->Fill(dt);
            else                 dtUpOverTwoGeVEBPeak_  ->Fill(dt);

            dtVSAeffHistEBPeak_ -> Fill(Aeff, dt); 
            dtVSAeffProfEBPeak_ -> Fill(Aeff, dt); 
          }
          else      {
            if      (Aeff < 6)   dtUpToThreeQuarterGeVEEPeak_->Fill(dt);
            else if (Aeff < 12)  dtUpToOneAndHalfGeVEEPeak_  ->Fill(dt);
            else if (Aeff < 24)  dtUpToThreeGeVEEPeak_       ->Fill(dt);
            else if (Aeff < 48)  dtUpToSixGeVEEPeak_         ->Fill(dt);
            else                 dtUpOverSixGeVEEPeak_       ->Fill(dt);

            dtVSAeffHistEEPeak_ -> Fill(Aeff, dt); 
            dtVSAeffProfEEPeak_ -> Fill(Aeff, dt); }
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
  massDiGammaEBHist_ ->Fit("massEB","","",0,1);
  pi0MassEB_  = massEB->GetParameter(1);
  pi0WidthEB_ = massEB->GetParameter(2);

  // get mass and width for EE 
  TF1 *massEE = new TF1("massEE","gaus(0)+pol3(3)",0,1); 
  massEE->SetParameter(0,500);
  massEE->SetParameter(1,0.12);
  massEE->SetParameter(2,0.02);  
  massDiGammaEEHist_ ->Fit("massEE","","",0,1);
  pi0MassEE_  = massEE->GetParameter(1);
  pi0WidthEE_ = massEE->GetParameter(2);

  // get mass and width for EP 
  TF1 *massEEP = new TF1("massEEP","gaus(0)+pol3(3)",0,1); 
  massEEP->SetParameter(0,500);
  massEEP->SetParameter(1,0.12);
  massEEP->SetParameter(2,0.02);  
  massDiGammaEEPHist_ ->Fit("massEEP","","",0,1);
  pi0MassEEP_  = massEEP->GetParameter(1);
  pi0WidthEEP_ = massEEP->GetParameter(2);
  
  // get mass and width for EM 
  TF1 *massEEM = new TF1("massEEM","gaus(0)+pol3(3)",0,1); 
  massEEM->SetParameter(0,500);
  massEEM->SetParameter(1,0.12);
  massEEM->SetParameter(2,0.02);  
  massDiGammaEEMHist_ ->Fit("massEEM","","",0,1);
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
      clusterTime timeAndUncertClusterA = timeAndUncertSingleCluster(bClusterA);
      if(timeAndUncertClusterA.timeErr <= 0) // if something went wrong combining the times, bail out
	continue;
      //std::pair<float,float> timeAndUncertClusterB = timeAndUncertSingleCluster(bClusterB);
      clusterTime timeAndUncertClusterB = timeAndUncertSingleCluster(bClusterB);
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

      if(isAfterPi0Selection)
	{
	  
	  //////////////////////////////////////////////////////////
	  // from here on I have pi0 candidates
	  // bClusterA and bClusterB are the two clusters making this candidate
	  diPhotonPeakOccupancyAny_  -> Fill(pi0Candidate.Eta(), pi0Candidate.Phi());
	  
	  dtDoubleClusterHistPi0Peak_            ->Fill(Dt);
	  if(isEB) dtDoubleClusterHistPi0PeakEB_ ->Fill(Dt);
	  else     dtDoubleClusterHistPi0PeakEE_ ->Fill(Dt);

	  dtPoolDoubleClusterHistPi0Peak_            ->Fill(Dt/errorDt);
	  if(isEB) dtPoolDoubleClusterHistPi0PeakEB_ ->Fill(Dt/errorDt);
	  else     dtPoolDoubleClusterHistPi0PeakEE_ ->Fill(Dt/errorDt);

	  dtVsPtDoubleClusterHistPi0Peak_            ->Fill(Pt,Dt);
	  if(isEB) dtVsPtDoubleClusterHistPi0PeakEB_ ->Fill(Pt,Dt);
	  else     dtVsPtDoubleClusterHistPi0PeakEE_ ->Fill(Pt,Dt);

	} // isAfterPi0Selection
      else
	{
	  dtDoubleClusterHistAny_     ->Fill(Dt);
	  dtPoolDoubleClusterHistAny_ ->Fill(Dt/errorDt);
	  dtVsPtDoubleClusterHistAny_ ->Fill(Pt,Dt);
	}// !isAfterPi0Selection
      
    }//loop on pairs
}// end doDoubleClusterResolutionPlots


// ---------------------------------------------------------------------------------------
// ------------------ Function to do double cluster resolution studies -------------------
void doFinalPlots()
{
  for (int sliceX=0; sliceX<numAeffBins; sliceX++)  {//looping on the X axis, at constant Aeff
    for (int binY=0; binY<numDtBins_; binY++)  {// looping in Delta t bins
      dtSliceVSAeffAny_[sliceX] ->SetBinContent( (binY+1), (dtVSAeffHistAny_->GetBinContent((sliceX+1),(binY+1))) );    
    }// end loop on Ybins 
  
    // std::cout << "integral is: " << dtSliceVSAeffAny_[sliceX] -> Integral() << std::endl;
    if( dtSliceVSAeffAny_[sliceX] -> Integral()  < 10 ) continue;
    
    // extract RMS and sigma for each Aeff=const slice
    float RMS       = dtSliceVSAeffAny_[sliceX] -> GetRMS();
    float RMSErr    = dtSliceVSAeffAny_[sliceX] -> GetRMSError();
    dtRMSVSAeffAny_ -> SetBinContent(sliceX+1, RMS);
    dtRMSVSAeffAny_ -> SetBinError(sliceX+1, RMSErr);
    
    TF1 *gauss = new TF1("dtFit","gaus",-DtMax_,DtMax_); // require min number entries
    gauss                    ->SetParLimits(1,-5,5); // limit on gaussian central 
    dtSliceVSAeffAny_[sliceX]->Fit("dtFit");
    float sigma     = gauss -> GetParameter(2);
    float sigmaErr  = gauss -> GetParError(2);
    dtSigmaAeffAny_ -> SetBinContent(sliceX+1, sigma);
    dtSigmaAeffAny_ -> SetBinError(sliceX+1, sigmaErr);
  }// end loop on Xslices
  
  // gf: add here sigma/RMS VS Aeff also for clusters matching the pi0 mass

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
	
  setBranchAddresses (chain, treeVars_);

  // Initialize output root file
  saving_ = new TFile(outputRootName_.c_str (),"recreate");

  // Initialize the histograms
  initializeHists();

  // FIXME
  // fit to mass to be made robust
  // re masses to a-priori values for now 
  pi0MassEB_  = 0.111;
  pi0WidthEB_ = 0.013; 
  pi0MassEE_  = 0.126;
  pi0WidthEE_  = 0.030;

  /////////////////////////////////////////////////////
  // Main loop over entries
  for (int entry = 0 ; (entry < nEntries && entry < numEvents_); ++entry)
  {
    chain->GetEntry (entry) ;

    speak_=false;
    if (entry<10 || entry%1000==0) speak_=true;

    if (speak_)  std::cout << "------> reading entry " << entry << " <------\n" ; 
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
