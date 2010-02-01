#include <iostream>
#include <string>

#include <map>
#include <vector>
#include <algorithm> 
#include <functional>

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



#define BarrelLimit 1.479
#define EndcapLimit 3.0
#define ADCtoGeVEB 0.039
#define ADCtoGeVEE 0.063

// initial authors P. Govoni et al
// authors: S. Cooper and G. Franzoni (UMN)

//! main program
int main (int argc, char** argv)
{

  // default output file
  std::string outputRootName = "outputHistos.root";
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
  std::string stringNumEvents        = "--n";

  std::vector<std::string> listOfFiles;
  int   numEvents    =-1;
  float	eTGammaMinEB   = 0.2;
  float s4s9GammaMinEB = 0.85;
  float eTPi0MinEB     = 0.65;
  float	eTGammaMinEE   = 0.250;
  float s4s9GammaMinEE = 0.85;
  float eTPi0MinEE     = 0.800;

  float sigmaNoiseEB      = 1.06; // ADC 
  float sigmaNoiseEE      = 2.10; // ADC
  float minAmpliOverSigma = 7;    // dimensionless
  
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
	std::cout << " --i <list of strings> list of input files" << std::endl ;     
	exit(1);      }

      
      else if (argv[v] == stringNumEvents) { // set number of events
	std::cout << "events number" << std::endl;
	numEvents=atoi(argv[v+1]);
	v++;
      }

      
      else if (argv[v] == stringETGammaMinEB) { // choose et cut for EB single cluster
	eTGammaMinEB = atof(argv[v+1]);
	v++;
      }
      
      else if (argv[v] == strings4s9GammaMinEB) { // choose cut for EB shower shape
	s4s9GammaMinEB = atof(argv[v+1]);
	v++;
      }
      
      else if (argv[v] == stringeTPi0MinEB) { // choose et cut for EB pi0 candidate
	eTPi0MinEB = atof(argv[v+1]);
	v++;
      }

      
      else if (argv[v] == stringETGammaMinEE) { // choose et cut for EE single cluster
	eTGammaMinEE = atof(argv[v+1]);
	v++;
      }
      
      else if (argv[v] == strings4s9GammaMinEE) { // choose cut for EE shower shape
	s4s9GammaMinEE = atof(argv[v+1]);
	v++;
      }
      
      else if (argv[v] == stringeTPi0MinEE) { // choose et cut for EE pi0 candidate
	eTPi0MinEE = atof(argv[v+1]);
	v++;
      }

      else if (argv[v] == stringOutFileName) { // set output file
	outputRootName = argv[v+1];
	v++;
      }

      // handle here the case of multiple arguments for input files
      else if (argv[v] == stringInputFileName){// && v<(argc-1) ) {

	for (int u=v+1; u<argc; u++) {
	  
	  if ( 0==std::string(argv[u]).find( stringGenericOption ) ){
	    if ( 0==listOfFiles.size())  {std::cout << "no input files listed" << std::cout;}
	    //else  {std::cout << "no more files listed, found: " << argv[u] << std::cout;}
	    break;
	  }

	  else {  listOfFiles.push_back(argv[u]);
	    v++;
	  }

	}// loop on arguments following --i

	continue;

      }//end 'if input files'

      
      else
	{std::cout << "input format unrecognized" << std::endl; exit(1);}

    }// loop over arguments input to the program


  
  if (listOfFiles.size()==0){
    std::cout << "\tno input file found" << std::endl;
    return(1);
  }
  else{
    std::cout << "\tfound " << listOfFiles.size() << " input files: " << std::endl;
    for(std::vector<std::string>::const_iterator  file_itr=listOfFiles.begin(); file_itr!=listOfFiles.end(); file_itr++){
      std::cout << "\t" << (*file_itr) << std::endl;
    }
  }
  




  // Tree construction
  TChain * chain = new TChain ("EcalTimePi0Analysis") ;
  std::vector<std::string>::const_iterator file_itr;
  for(file_itr=listOfFiles.begin(); file_itr!=listOfFiles.end(); file_itr++){
    chain->Add( (*file_itr).c_str() );
  }
  int nEntries = chain->GetEntries () ;
  if (numEvents==-1) numEvents = nEntries;
  std::cout << "\n\tFOUND "         <<  listOfFiles.size() << " input files" << std::endl ;    
  std::cout << "\n\tFOUND "         <<  nEntries << " events" << std::endl ;    
  std::cout << "\tWILL run on: "    <<  numEvents << " events" << std::endl;
  std::cout << "\tOutput file: "    <<  outputRootName << std::endl;
  std::cout << "\teTGammaMinEB: "   <<  eTGammaMinEB << std::endl;
  std::cout << "\ts4s9GammaMinEB: " <<  s4s9GammaMinEB << std::endl;
  std::cout << "\teTPi0MinEB: "     <<  eTPi0MinEB << std::endl;
  std::cout << "\teTGammaMinEE: "   <<  eTGammaMinEE << std::endl;
  std::cout << "\ts4s9GammaMinEE: " <<  s4s9GammaMinEE << std::endl;
  std::cout << "\teTPi0MinEE: "     <<  eTPi0MinEE << std::endl;
	
  EcalTimePi0TreeContent treeVars ; 
  setBranchAddresses (chain, treeVars) ;

  // Initialize output root file
  TFile saving (outputRootName.c_str (),"recreate") ;
  saving.cd () ;

  // Initialize histograms -- xtals
  TH1F* xtalEnergyHist_ = new TH1F("XtalEnergy","Crystal energy;GeV",110,-1,10);
  TH1F* xtalTimeHist_ = new TH1F("XtalTime","Time of all crystals;ns",150,-75,75);
  TH1F* xtalIEtaHist_ = new TH1F("xtalIEta","i#eta of crystal",171,-85,86);
  TH1F* xtalIPhiHist_ = new TH1F("xtalIPhi","i#phi of crystal",361,1,361);
  TH1F* xtalIXHist_ = new TH1F("xtalIX","ix of crystal",101,1,101);
  TH1F* xtalIYHist_ = new TH1F("xtalIY","iy of crystal",101,1,101);
  // TH1F* xtalStatusHist_ = new TH1F("XtalStatus","Crystal status flag",16,0,15);
  // TH2F* xtalOccupancyHistEB_ = new TH2F("XtalOccupancyEB","Crystal occupancy;i#phi;i#eta",360,1.,361.,172,-86,86);

  // Initialize histograms -- BasicClusters
  TH1F* BCNumPerEventHist_ = new TH1F("BCNumPerEvent","Number of BC per event",100,0,100);
  TH1F* BCNumCrysHist_ = new TH1F("BCNumCrys","Number of crystals per BC",10,0,10);
  TH1F* BCEnergyHist_ = new TH1F("BCEnergy","Energy of BCs;GeV",100,0,25);
  TH1F* BCEtHist_ = new TH1F("BCEt","E_{T} of BCs;GeV",100,0,25);
  TH2F* BCOccupancyEBHist_  = new TH2F("BCOccupancyEB","BC occupancy;i#eta;i#phi",171,-85,86,361,1.,361.);
  TH2F* BCOccupancyEEHist_  = new TH2F("BCOccupancyEE","BC occupancy;ix;iy",101,1.,101.,101,1,101);
  TH2F* BCOccupancyEEPHist_  = new TH2F("BCOccupancyEEP","BC occupancy;ix;iy",101,1.,101.,101,1,101);
  TH2F* BCOccupancyEEMHist_  = new TH2F("BCOccupancyEEM","BC occupancy;ix;iy",101,1.,101.,101,1,101);
  TH2F* BCOccupancyHistAny_ = new TH2F("BCOccupancyAny","BC occupancy;#eta;#phi",50,-3.5,3.5,50,-1*TMath::Pi(),TMath::Pi());
  TH1F* BCEtaHist_ = new TH1F("Cluster #eta","#eta of cluster",171,-3.5,3.5);
  TH1F* BCPhiHist_ = new TH1F("Cluster #phi","#phi of cluster",50,-1*TMath::Pi(),TMath::Pi());
  TH1F* BCClusterShapeEEHist_ = new TH1F("EE cluster shape","e2x2 / e3x3",65,-0.1,1.2);
  TH1F* BCClusterShapeEEPHist_ = new TH1F("EEP cluster shape","e2x2 / e3x3",65,-0.1,1.2);
  TH1F* BCClusterShapeEEMHist_ = new TH1F("EEM cluster shape","e2x2 / e3x3",65,-0.1,1.2);
  TH1F* BCClusterShapeEBHist_  = new TH1F("EB cluster shape","e2x2 / e3x3",65,-0.1,1.2);

  // Initialize histograms -- diphotons control plots
  TH1F* massDiGammaHist_      = new TH1F("massDiGamma","m(#gamma#gamma)",50,0,0.500);
  TH1F* massDiGammaEBHist_    = new TH1F("massDiGamma EB","m(#gamma#gamma) EB",50,0,0.500);
  TH1F* massDiGammaEEHist_    = new TH1F("massDiGamma EE","m(#gamma#gamma) EE",50,0,0.500);
  TH1F* massDiGammaEEPHist_   = new TH1F("massDiGamma EEP","m(#gamma#gamma) EEP",50,0,0.500);
  TH1F* massDiGammaEEMHist_   = new TH1F("massDiGamma EEM","m(#gamma#gamma) EEM",50,0,0.500);
  TH2F* diPhotonOccupancyAny_ = new TH2F("di-photon occupancy","di-photon occupancy;#eta;#phi",50,-3.5,3.5,50,-1*TMath::Pi(),TMath::Pi());

  // Initialize histograms -- single cluster resolution

  TH1F*     dtUpToQuarterGeVEB   = new TH1F("EB #delta(t),   A_{eff} up to 1/4 GeV", "EB #delta(t),   ~5<A_{eff}/#sigma_{N}<6", 400, -20, 20); 
  TH1F*     dtUpToHalfGeVEB      = new TH1F("EB #delta(t),   A_{eff} up to half GeV", "EB #delta(t),   6<A_{eff}/#sigma_{N}<12", 400, -20, 20); 
  TH1F*     dtUpToOneGeVEB       = new TH1F("EB #delta(t),   A_{eff} up to one GeV", "EB #delta(t),   12<A_{eff}/#sigma_{N}<24", 200, -10, 10); 
  TH1F*     dtUpToTwoGeVEB       = new TH1F("EB #delta(t),   A_{eff} up to two GeV", "EB #delta(t),   24<A_{eff}/#sigma_{N}<48", 200, -10, 10); 
  TH1F*     dtUpOverTwoGeVEB     = new TH1F("EB #delta(t),   A_{eff} over two GeV", "EB #delta(t),   A_{eff}/#sigma_{N}>48", 200, -10, 10); 
          				    
  TH1F*     dtUpToThreeQuarterGeVEE = new TH1F("EE #delta(t),   A_{eff} up to 3/4 GeV", "EE #delta(t),   A_{eff}/#sigma_{N}<6", 200, -10, 10); 
  TH1F*     dtUpToOneAndHalfGeVEE   = new TH1F("EE #delta(t),   A_{eff} up to one&1/2 GeV", "EE #delta(t),   6<A_{eff}/#sigma_{N}<12", 200, -10, 10); 
  TH1F*     dtUpToThreeGeVEE        = new TH1F("EE #delta(t),   A_{eff} up to three GeV", "EE #delta(t),   12<A_{eff}/#sigma_{N}<24", 200, -10, 10); 
  TH1F*     dtUpToSixGeVEE          = new TH1F("EE #delta(t),   A_{eff} up to six GeV", "EE #delta(t),   24<A_{eff}/#sigma_{N}<48", 200, -10, 10); 
  TH1F*     dtUpOverSixGeVEE        = new TH1F("EE #delta(t),   A_{eff} over six GeV", "EE #delta(t),   A_{eff}/#sigma_{N}>48", 200, -10, 10); 

  TH2F*     dtVSAeffHistAny_ = new TH2F("#delta(t) VS A_{eff}/#sigma_{N}","#delta(t) VS A_{eff}/#sigma_{N}",250,0.,250,100,0,10);
  TH2F*     dtVSAeffHistEB_  = new TH2F("EB:  #delta(t)  VS  A_{eff}/#sigma_{N}","EB:  #delta(t)  VS  A_{eff}/#sigma_{N}",250,0.,250,100,0,10);
  TH2F*     dtVSAeffHistEE_  = new TH2F("EE:  #delta(t)  VS  A_{eff}/#sigma_{N}","EE:  #delta(t)  VS  A_{eff}/#sigma_{N}",250,0.,250,100,0,10);
  TProfile* dtVSAeffProfAny_ = new TProfile("#delta(t)  VS  A_{eff}/#sigma_{N} prof","#delta(t) VS A_{eff}/#sigma_{N} prof",125,0.,250,0,10);
  TProfile* dtVSAeffProfEB_  = new TProfile("EB:  #delta(t)   VS  A_{eff}/#sigma_{N} prof","EB:  #delta(t)  VS  A_{eff}/#sigma_{N} prof",125,0.,250,0,10);
  TProfile* dtVSAeffProfEE_  = new TProfile("EE:  #delta(t)   VS  A_{eff}/#sigma_{N} prof","EE:  #delta(t)  VS  A_{eff}/#sigma_{N} prof",125,0.,250,0,10);

  // Initialize histograms -- selection on pi0 candidates 
  TH2F* diPhotonPeakOccupancyAny_     = new TH2F("#pi_{0} occupancy","di-photon peak;#eta;#phi",50,-3.5,3.5,50,-1*TMath::Pi(),TMath::Pi());
  TH2F* diPhotonSidesOccupancyAny_    = new TH2F("di-photon side-bands","di-photon side-bands;#eta;#phi",50,-3.5,3.5,50,-1*TMath::Pi(),TMath::Pi());


  // Initialize histograms -- single cluster resolution in pi0 peak

  TH1F*     dtUpToQuarterGeVEBPeak   = new TH1F("EBPeak: #delta(t) A_{eff} up to 1/4 GeV", "EBPeak: #delta(t) ~5<A_{eff}/#sigma_{N}<6", 400, -20, 20); 
  TH1F*     dtUpToHalfGeVEBPeak      = new TH1F("EBPeak: #delta(t) A_{eff} up to half GeV", "EBPeak: #delta(t) 6<A_{eff}/#sigma_{N}<12", 400, -20, 20); 
  TH1F*     dtUpToOneGeVEBPeak       = new TH1F("EBPeak: #delta(t) A_{eff} up to one GeV", "EBPeak: #delta(t) 12<A_{eff}/#sigma_{N}<24", 200, -10, 10); 
  TH1F*     dtUpToTwoGeVEBPeak       = new TH1F("EBPeak: #delta(t) A_{eff} up to two GeV", "EBPeak: #delta(t) 24<A_{eff}/#sigma_{N}<48", 200, -10, 10); 
  TH1F*     dtUpOverTwoGeVEBPeak     = new TH1F("EBPeak: #delta(t) A_{eff} over two GeV", "EBPeak: #delta(t) A_{eff}/#sigma_{N}>48", 200, -10, 10); 
          				    
  TH1F*     dtUpToThreeQuarterGeVEEPeak = new TH1F("EEPeak: #delta(t) A_{eff} up to 3/4 GeV", "EEPeak: #delta(t) A_{eff}/#sigma_{N}<6", 200, -10, 10); 
  TH1F*     dtUpToOneAndHalfGeVEEPeak    = new TH1F("EEPeak: #delta(t) A_{eff} up to one&1/2 GeV", "EEPeak: #delta(t) 6<A_{eff}/#sigma_{N}<12", 200, -10, 10); 
  TH1F*     dtUpToThreeGeVEEPeak        = new TH1F("EEPeak: #delta(t) A_{eff} up to three GeV", "EEPeak: #delta(t) 12<A_{eff}/#sigma_{N}<24", 200, -10, 10); 
  TH1F*     dtUpToSixGeVEEPeak          = new TH1F("EEPeak: #delta(t) A_{eff} up to six GeV", "EEPeak: #delta(t) 24<A_{eff}/#sigma_{N}<48", 200, -10, 10); 
  TH1F*     dtUpOverSixGeVEEPeak        = new TH1F("EEPeak: #delta(t) A_{eff} over six GeV", "EEPeak: #delta(t) A_{eff}/#sigma_{N}>48", 200, -10, 10); 

  TH2F*     dtVSAeffHistAnyPeak_ = new TH2F("Peak: #delta(t) VS A_{eff}/#sigma_{N}","Peak: #delta(t) VS A_{eff}/#sigma_{N}",250,0.,250,100,0,10);
  TH2F*     dtVSAeffHistEBPeak_  = new TH2F("EBPeak: #delta(t) VS A_{eff}/#sigma_{N}","EBPeak #delta(t) VS A_{eff}/#sigma_{N}",250,0.,250,100,0,10);
  TH2F*     dtVSAeffHistEEPeak_  = new TH2F("EEPeak: E#delta(t) VS A_{eff}/#sigma_{N}","EEPeak: #delta(t) VS A_{eff}/#sigma_{N}",250,0.,250,100,0,10);
  TProfile* dtVSAeffProfAnyPeak_ = new TProfile("Peak: #delta(t) VS A_{eff}/#sigma_{N} prof","Peak: #delta(t) VS A_{eff}/#sigma_{N} prof",125,0.,250,0,10);
  TProfile* dtVSAeffProfEBPeak_  = new TProfile("EBPeak: #delta(t) VS A_{eff}/#sigma_{N} prof","EBPeak #delta(t) VS A_{eff}/#sigma_{N} prof",125,0.,250,0,10);
  TProfile* dtVSAeffProfEEPeak_  = new TProfile("EEPeak: #delta(t) VS A_{eff}/#sigma_{N} prof","EEPeak: #delta(t) VS A_{eff}/#sigma_{N} prof",125,0.,250,0,10);


  /////////////////////////////////////////////////////
  // FIRST loop over entries
  for (int entry = 0 ; (entry < nEntries && entry < numEvents); ++entry)
    {
      chain->GetEntry (entry) ;

      bool speak=false;
      if (entry<10 || entry%1000==0) speak=true;

      if (speak) std::cout << "------> reading entry " << entry << " <------\n" ; 


      // loop on calorimetric quantities

      if (speak)  std::cout << "  found " << treeVars.nSuperClusters << " superclusters" << std::endl ;
      if (speak)  std::cout << "  found " << treeVars.nClusters << " basic clusters" << std::endl ;
      BCNumPerEventHist_->Fill(treeVars.nClusters);


        /////////////////////////////////////////////////////////////////
	// (FIRST) loop on basic cluster - to fill basic cluster histograms
	for (int bCluster=0; bCluster < treeVars.nClusters; bCluster++)
	  {

	    float eBC=0; // calculate energy of BC for validation
	    for (int cryInBC=0; cryInBC < treeVars.nXtalsInCluster[bCluster]; cryInBC++){
	      eBC+= treeVars.xtalInBCEnergy[bCluster][cryInBC];}

            BCEnergyHist_->Fill(treeVars.clusterEnergy[bCluster]);
            BCEtHist_->Fill(treeVars.clusterTransverseEnergy[bCluster]);
            BCNumCrysHist_->Fill(treeVars.nXtalsInCluster[bCluster]);

	    // basic cluster occupancy in physics coordinates
	    BCOccupancyHistAny_ -> Fill(treeVars.clusterEta[bCluster],treeVars.clusterPhi[bCluster]);
	    BCEtaHist_ -> Fill(treeVars.clusterEta[bCluster]);
	    BCPhiHist_ -> Fill(treeVars.clusterPhi[bCluster]);

	    //  basic cluster occupancy in detector coordinates, using first cry of BC as a representative
            if(treeVars.xtalInBCIEta[bCluster][0] != -999999)                                        // ieta=-999999 tags EE
              BCOccupancyEBHist_->Fill(treeVars.xtalInBCIEta[bCluster][0],treeVars.xtalInBCIPhi[bCluster][0]);

            else if (treeVars.xtalInBCIx[bCluster][0] != -999999 && treeVars.clusterEta[bCluster]>0){ // ix=-999999 tags EB
              BCOccupancyEEHist_ ->Fill(treeVars.xtalInBCIx[bCluster][0],treeVars.xtalInBCIy[bCluster][0]);
              BCOccupancyEEPHist_->Fill(treeVars.xtalInBCIx[bCluster][0],treeVars.xtalInBCIy[bCluster][0]);}

            else if (treeVars.xtalInBCIx[bCluster][0] != -999999 && treeVars.clusterEta[bCluster]<0){ // ix=-999999 tags EB
              BCOccupancyEEHist_ ->Fill(treeVars.xtalInBCIx[bCluster][0],treeVars.xtalInBCIy[bCluster][0]);
              BCOccupancyEEMHist_->Fill(treeVars.xtalInBCIx[bCluster][0],treeVars.xtalInBCIy[bCluster][0]);}

	    if (speak)  std::cout << "\tbCluster: num "               << bCluster 
				  << "\t eBC: "                      << treeVars.clusterEnergy[bCluster]
				  << "\t eBC_predicted: "            << eBC
	      //<< "\n\t et: "                     << treeVars.clusterTransverseEnergy[bCluster]
	      //<< "\t predicted et: "             << treeVars.clusterEnergy[bCluster]*sin(2*atan(exp(-1* treeVars.clusterEta[bCluster] )) )
				  << " eta: "                        << treeVars.clusterEta[bCluster]
				  << "\n\t num crystals: "           << treeVars.nXtalsInCluster[bCluster]
				  << "\n\t\tfirst crystal:  \tieta " << treeVars.xtalInBCIEta[bCluster][0] 
				  << "\teta "                        << treeVars.xtalInBCEta[bCluster][0] 
				  << " \t energy "                   << treeVars.xtalInBCEnergy[bCluster][0] 
				  << " \t ADC "                      << treeVars.xtalInBCAmplitudeADC[bCluster][0] 
				  << " \t time "                     << treeVars.xtalInBCTime[bCluster][0] 
				  << std::endl;

	    for(int thisCry=0; thisCry<treeVars.nXtalsInCluster[bCluster]; thisCry++)
	      {
		if (treeVars.xtalInBCIEta[bCluster][thisCry]!=-999999)  xtalIEtaHist_ -> Fill (treeVars.xtalInBCIEta[bCluster][thisCry]);
		if (treeVars.xtalInBCIPhi[bCluster][thisCry]!=-999999)  xtalIPhiHist_ -> Fill (treeVars.xtalInBCIPhi[bCluster][thisCry]);
		if (treeVars.xtalInBCIx[bCluster][thisCry]  !=-999999)  xtalIXHist_   -> Fill (treeVars.xtalInBCIx[bCluster][thisCry]);
		if (treeVars.xtalInBCIy[bCluster][thisCry]  !=-999999)  xtalIYHist_   -> Fill (treeVars.xtalInBCIy[bCluster][thisCry]);
	      }
	    
	  }// end  (FIRST) to fill basic cluster histograms
	
	if (speak) std::cout << "  found " << treeVars.nXtals << " crystals\n" ;    


	
      float eTA, eTB ;
      float e22A, e33A,    e22B, e33B;
      float eTGammaMinA,   eTGammaMinB;
      float s4s9GammaMinA, s4s9GammaMinB;
      bool  AisEB,         BisEB;
      float eTPi0Min;
      
      ////////////////////////////////////////////////////////////////////      
      // (FIRST) loop on basic cluster - to build pi0 candidates and get the mass
      for (int bCluster=0; bCluster < treeVars.nClusters; bCluster++)
	for (int bClusterA=0; bClusterA < treeVars.nClusters; bClusterA++)
	  {
	    eTA = treeVars.clusterTransverseEnergy[bClusterA];
	    
	    e22A = treeVars.clusterE2x2[bClusterA];
	    e33A = treeVars.clusterE3x3[bClusterA];
	    
	    // discriminate between EE and EB and set thresholds accordingly
	    if ( fabs(treeVars.clusterEta[bClusterA]) < BarrelLimit) {
	      AisEB         = true;
	      eTGammaMinA   = eTGammaMinEB;
	      s4s9GammaMinA = s4s9GammaMinEB;
	    }
	    else{
	      AisEB         = false;
	      eTGammaMinA   = eTGammaMinEE;
	      s4s9GammaMinA = s4s9GammaMinEB;
	    }
	    

	  if(treeVars.clusterEta[bClusterA]<-1.4)     {
	    BCClusterShapeEEHist_  -> Fill(e22A/e33A);
	    BCClusterShapeEEMHist_ -> Fill(e22A/e33A);}
	  else if(treeVars.clusterEta[bClusterA]>1.4) {
	    BCClusterShapeEEHist_  -> Fill(e22A/e33A);
	    BCClusterShapeEEPHist_ -> Fill(e22A/e33A);}
	  else	                                      
	    BCClusterShapeEBHist_  -> Fill(e22A/e33A);


	  // first selecton cut: photon candidate Et
	  if( eTA < eTGammaMinA ) continue;
	  
	  // second selection cut: cluster shape
	  if ( e22A/e33A < s4s9GammaMinA ) continue;
	  
	  for (int bClusterB=(bClusterA+1); bClusterB < treeVars.nClusters; bClusterB++)
	    {

	      eTB = treeVars.clusterTransverseEnergy[bClusterB];
	      
	      e22B = treeVars.clusterE2x2[bClusterB];
	      e33B = treeVars.clusterE3x3[bClusterB];
	      
	      // discriminate between EE and EB and set thresholds accordingly
	      if ( fabs(treeVars.clusterEta[bClusterB]) < BarrelLimit) {
		BisEB         = true;
		eTGammaMinB   = eTGammaMinEB;
		s4s9GammaMinB = s4s9GammaMinEB;
	      }
	      else{
		BisEB         = false;
		eTGammaMinB   = eTGammaMinEE;
		s4s9GammaMinB = s4s9GammaMinEB;
	      }
	      

	      // first selecton cut: photon candidate Et
	      if( eTB < eTGammaMinB ) continue;

	      // second selection cut: cluster shape
	      if ( e22B/e33B < s4s9GammaMinB ) continue;
	      
	      math::PtEtaPhiMLorentzVectorD gammaA (eTA, treeVars.clusterEta[bClusterA], treeVars.clusterPhi[bClusterA], 0);
	      math::PtEtaPhiMLorentzVectorD gammaB (eTB, treeVars.clusterEta[bClusterB], treeVars.clusterPhi[bClusterB], 0);
	      
	      math::PtEtaPhiMLorentzVectorD pi0Candidate = gammaA + gammaB;
	      
	      // std::cout << "gammaA: " << gammaA << " " << gammaA.M() << "\t\t gammaB: " << gammaB << " " << gammaB.M() << std::endl;
	      // std::cout << "pi0Candidate: " << pi0Candidate << " " << pi0Candidate.M() << std::endl;

	      if ( fabs(pi0Candidate.Eta()) < BarrelLimit) {
		eTPi0Min      = eTPi0MinEB;}
	      else{
		eTPi0Min      = eTPi0MinEE;}
	      


	      // third selection cut: pi0 candidate Et
	      if(pi0Candidate.Et() < eTPi0Min ) continue;
	      
	      massDiGammaHist_ -> Fill(pi0Candidate.M());
	      
	      if(treeVars.clusterEta[bClusterA]<-1.4){
		massDiGammaEEHist_  -> Fill(pi0Candidate.M());
		massDiGammaEEMHist_ -> Fill(pi0Candidate.M());}
	      else if(treeVars.clusterEta[bClusterA]>1.4) {
		massDiGammaEEHist_  -> Fill(pi0Candidate.M());
		massDiGammaEEPHist_ -> Fill(pi0Candidate.M());}
	      else	     
		massDiGammaEBHist_  -> Fill(pi0Candidate.M());
	      
	      // occupancy of all candidates (this is FIRST loop)
	      diPhotonOccupancyAny_ -> Fill(pi0Candidate.Eta(), pi0Candidate.Phi());
	      
	      
	    }//loop on candidateB
	  
	  }//loop on candidateA - (FIRST) to build pi0 candidates and get the mass
      
    }   // end of first loop over entries
        // now you have di-mass plots filled => get masses
  
  

  // get mass and width for EB 
  TF1 *massEB = new TF1("massEB","gaus(0)+pol3(3)",0,1); 
  //massEB->SetParameter(0,500);
  massEB->SetParameter(1,0.12);
  massEB->SetParameter(2,0.015);  
  massDiGammaEBHist_ ->Fit("massEB","","",0,1);
  float pi0MassEB  = massEB->GetParameter(1);
  float pi0WidthEB = massEB->GetParameter(2);

  // get mass and width for EE 
  TF1 *massEE = new TF1("massEE","gaus(0)+pol3(3)",0,1); 
  massEE->SetParameter(0,500);
  massEE->SetParameter(1,0.12);
  massEE->SetParameter(2,0.02);  
  massDiGammaEEHist_ ->Fit("massEE","","",0,1);
  float pi0MassEE  = massEE->GetParameter(1);
  float pi0WidthEE = massEE->GetParameter(2);

  // get mass and width for EP 
  TF1 *massEEP = new TF1("massEEP","gaus(0)+pol3(3)",0,1); 
  massEEP->SetParameter(0,500);
  massEEP->SetParameter(1,0.12);
  massEEP->SetParameter(2,0.02);  
  massDiGammaEEPHist_ ->Fit("massEEP","","",0,1);
  float pi0MassEEP  = massEEP->GetParameter(1);
  float pi0WidthEEP = massEEP->GetParameter(2);
  
  // get mass and width for EM 
  TF1 *massEEM = new TF1("massEEM","gaus(0)+pol3(3)",0,1); 
  massEEM->SetParameter(0,500);
  massEEM->SetParameter(1,0.12);
  massEEM->SetParameter(2,0.02);  
  massDiGammaEEMHist_ ->Fit("massEEM","","",0,1);
  float pi0MassEEM  = massEEM->GetParameter(1);
  float pi0WidthEEM = massEEM->GetParameter(2);
  
  // FIX
  // fit to mass to be made robust
  // set masses to a-priori values for now 
  pi0MassEB  = 0.111;
  pi0WidthEB = 0.013; 
  pi0MassEE  = 0.126;
  pi0WidthEE  = 0.030;


  /////////////////////////////////////////////////////
  // SECOND loop over entries
  for (int entry = 0 ; (entry < nEntries && entry < numEvents); ++entry)
    {// SECOND loop over entries
      chain->GetEntry (entry) ;

	////////////////////////////////////////////////////////////////////      
	// (SECOND) loop on basic cluster - to perform single-cluster resolutuon studies
	// perform single-cluster resolution study
	for (int bCluster=0; bCluster < treeVars.nClusters; bCluster++)
	  {// loop on bc

	    // loop on the cry components of a basic cluster
	    // put single cluster time resolution studies  here!
	    for(int thisCry=0; thisCry<treeVars.nXtalsInCluster[bCluster]; thisCry++)
	      {
		bool  thisIsInEB=false;
		float sigmaNoiseOfThis=0;
		if(treeVars.xtalInBCIEta[bCluster][thisCry]      !=-999999)   {sigmaNoiseOfThis=sigmaNoiseEB; thisIsInEB=true;}
		else if (treeVars.xtalInBCIy[bCluster][thisCry]  !=-999999)   {sigmaNoiseOfThis=sigmaNoiseEE; thisIsInEB=false;}
		else {std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}

		float ampliOfThis = treeVars.xtalInBCAmplitudeADC[bCluster][thisCry] / sigmaNoiseOfThis; 
		if( ampliOfThis < minAmpliOverSigma) continue;


		for(int thatCry=thisCry+1; thatCry<treeVars.nXtalsInCluster[bCluster]; thatCry++)
		  {
		    float sigmaNoiseOfThat=0;
		    if(treeVars.xtalInBCIEta[bCluster][thatCry]      !=-999999)   sigmaNoiseOfThat=sigmaNoiseEB;
		    else if (treeVars.xtalInBCIy[bCluster][thatCry]  !=-999999)   sigmaNoiseOfThat=sigmaNoiseEE;
		    else {std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}
		    
		    float ampliOfThat = treeVars.xtalInBCAmplitudeADC[bCluster][thatCry] / sigmaNoiseOfThat; 
		    if( ampliOfThat < minAmpliOverSigma) continue;
		    
		    float Aeff = ampliOfThis * ampliOfThat / sqrt( pow(ampliOfThis,2) + pow(ampliOfThat,2) );
		    float dt  = treeVars.xtalInBCTime[bCluster][thisCry] - treeVars.xtalInBCTime[bCluster][thatCry]; 

		    // for debug
		    //std::cout << "ampliOfThis: " << ampliOfThis << "\tampliOfThat: " << ampliOfThat
		    //          << "\n timeOfThis: " << treeVars.xtalInBCTime[bCluster][thisCry] << "\ttimeOfThat: " << treeVars.xtalInBCTime[bCluster][thatCry]
		    //          << "\n Aeff: " << Aeff << "\tdt: " << dt 
		    //          << std::endl;
		    
		    dtVSAeffHistAny_  -> Fill(Aeff, dt); 
		    dtVSAeffProfAny_  -> Fill(Aeff, dt); 
		    if (thisIsInEB) {
		      if      (Aeff < 6)   dtUpToQuarterGeVEB->Fill(dt);
		      else if (Aeff < 12)  dtUpToHalfGeVEB   ->Fill(dt);
		      else if (Aeff < 24)  dtUpToOneGeVEB    ->Fill(dt);
		      else if (Aeff < 48)  dtUpToTwoGeVEB    ->Fill(dt);
		      else                 dtUpOverTwoGeVEB  ->Fill(dt);
		      
		      dtVSAeffHistEB_ -> Fill(Aeff, fabs(dt)); 
		      dtVSAeffProfEB_ -> Fill(Aeff, fabs(dt)); 
		    }
		    else      {
		      if      (Aeff < 6)   dtUpToThreeQuarterGeVEE->Fill(dt);
		      else if (Aeff < 12)  dtUpToOneAndHalfGeVEE  ->Fill(dt);
		      else if (Aeff < 24)  dtUpToThreeGeVEE       ->Fill(dt);
		      else if (Aeff < 48)  dtUpToSixGeVEE         ->Fill(dt);
		      else                 dtUpOverSixGeVEE       ->Fill(dt);

		      dtVSAeffHistEE_ -> Fill(Aeff, fabs(dt)); 
		      dtVSAeffProfEE_ -> Fill(Aeff, fabs(dt)); }
		    
		  }// loop on thatCry
		
	      }// loop on thisCry
	  }// end loop on bc
	
	

	
      float eTA, eTB ;
      float e22A, e33A,    e22B, e33B;
      float eTGammaMinA,   eTGammaMinB;
      float s4s9GammaMinA, s4s9GammaMinB;
      bool  AisEB,         BisEB;
      float eTPi0Min;
      std::vector<int> clustersFromPi0;

      ////////////////////////////////////////////////////////////////////      
      // (SECOND) loop on basic cluster - to perform di-cluster studies
      // loop on cluster A
      for (int bClusterA=0; bClusterA < treeVars.nClusters; bClusterA++)
	{
	  eTA = treeVars.clusterTransverseEnergy[bClusterA];

	  e22A = treeVars.clusterE2x2[bClusterA];
	  e33A = treeVars.clusterE3x3[bClusterA];

	  // discriminate between EE and EB and set thresholds accordingly
	  if ( fabs(treeVars.clusterEta[bClusterA]) < BarrelLimit) {
	    AisEB         = true;
	    eTGammaMinA   = eTGammaMinEB;
	    s4s9GammaMinA = s4s9GammaMinEB;
	  }
	  else{
	    AisEB         = false;
	    eTGammaMinA   = eTGammaMinEE;
	    s4s9GammaMinA = s4s9GammaMinEB;
	  }

	  // first selecton cut: photon candidate Et
	  if( eTA < eTGammaMinA ) continue;
	  
	  // second selection cut: cluster shape
	  if ( e22A/e33A < s4s9GammaMinA ) continue;
	  
	  for (int bClusterB=(bClusterA+1); bClusterB < treeVars.nClusters; bClusterB++)
	    {

	      eTB = treeVars.clusterTransverseEnergy[bClusterB];
	      
	      e22B = treeVars.clusterE2x2[bClusterB];
	      e33B = treeVars.clusterE3x3[bClusterB];
	      
	      // discriminate between EE and EB and set thresholds accordingly
	      if ( fabs(treeVars.clusterEta[bClusterB]) < BarrelLimit) {
		BisEB         = true;
		eTGammaMinB   = eTGammaMinEB;
		s4s9GammaMinB = s4s9GammaMinEB;
	      }
	      else{
		BisEB         = false;
		eTGammaMinB   = eTGammaMinEE;
		s4s9GammaMinB = s4s9GammaMinEB;
	      }
	      

	      // first selecton cut: photon candidate Et
	      if( eTB < eTGammaMinB ) continue;

	      // second selection cut: cluster shape
	      if ( e22B/e33B < s4s9GammaMinB ) continue;

	      // now build the pi0 candidate
	      math::PtEtaPhiMLorentzVectorD gammaA (eTA, treeVars.clusterEta[bClusterA], treeVars.clusterPhi[bClusterA], 0);
	      math::PtEtaPhiMLorentzVectorD gammaB (eTB, treeVars.clusterEta[bClusterB], treeVars.clusterPhi[bClusterB], 0);
	      math::PtEtaPhiMLorentzVectorD pi0Candidate = gammaA + gammaB;
	      
	      bool candidateIsEB=false;
	      if ( fabs(pi0Candidate.Eta()) < BarrelLimit) {
		eTPi0Min      = eTPi0MinEB;	
		candidateIsEB = true;	      }
	      else{		
		eTPi0Min      = eTPi0MinEE;
		candidateIsEB = false;	      }

	      // third selection cut: pi0 candidate Et
	      if(pi0Candidate.Et() < eTPi0Min ) continue;



	      /////////////////////////////////////////////////////////////
	      // here I have di-gamma pairs that pass cuts
	      // now select pi0's based on the mass
	      /////////////////////////////////////////////////////////////
		
	      float nSigma =1;
	      // reject sidebands in EB
	      if( candidateIsEB &&
		  (pi0Candidate.M() < pi0MassEB-nSigma*pi0WidthEB ||
		   pi0MassEB+nSigma*pi0WidthEB > pi0Candidate.M())
		  ) {
		diPhotonSidesOccupancyAny_ -> Fill(pi0Candidate.Eta(), pi0Candidate.Phi());
		continue;}
	      
	      // reject sidebands in EE
	      if( (!candidateIsEB) &&
		  (pi0Candidate.M() < pi0MassEE-nSigma*pi0WidthEE ||
		   pi0MassEE+nSigma*pi0WidthEE > pi0Candidate.M())
		  ) {
		diPhotonSidesOccupancyAny_ -> Fill(pi0Candidate.Eta(), pi0Candidate.Phi());
		continue;}
	      

	      /////////////////////////////////////////////////////////////
	      // from here on I have pi0 candidates
	      // bClusterA and bClusterB are the two clusters making this candidate
	      diPhotonPeakOccupancyAny_  -> Fill(pi0Candidate.Eta(), pi0Candidate.Phi());
	

	      // add cluster to the list of those coming from pi0 if not already added
	      if( std::find( clustersFromPi0.begin(), clustersFromPi0.end(), bClusterA)  == clustersFromPi0.end()  )
		 clustersFromPi0.push_back(bClusterA);
	      if( std::find( clustersFromPi0.begin(), clustersFromPi0.end(), bClusterB)  == clustersFromPi0.end()  )
		 clustersFromPi0.push_back(bClusterB);

	    }//loop on clusterB
	}//loop on clusterA






      

      // now do single cluster resulution study on the clcusters which passed pi0 peak selection
      std::vector<int>::iterator bcFromPi0Itr;
      for(bcFromPi0Itr  = clustersFromPi0.begin(); 
	  bcFromPi0Itr != clustersFromPi0.end();
	  bcFromPi0Itr++)
	{
	  int bCluster = (*bcFromPi0Itr);
	  
	  // loop on the cry components of a basic cluster
	  for(int thisCry=0; thisCry<treeVars.nXtalsInCluster[bCluster]; thisCry++)
	    {
	      bool  thisIsInEBPeak=false;
	      float sigmaNoiseOfThis=0;
	      if(treeVars.xtalInBCIEta[bCluster][thisCry]      !=-999999)   {sigmaNoiseOfThis=sigmaNoiseEB; thisIsInEBPeak=true;}
	      else if (treeVars.xtalInBCIy[bCluster][thisCry]  !=-999999)   {sigmaNoiseOfThis=sigmaNoiseEE; thisIsInEBPeak=false;}
	      else {std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}
	      
	      float ampliOfThis = treeVars.xtalInBCAmplitudeADC[bCluster][thisCry] / sigmaNoiseOfThis; 
	      if( ampliOfThis < minAmpliOverSigma) continue;
	      
	      
	      for(int thatCry=thisCry+1; thatCry<treeVars.nXtalsInCluster[bCluster]; thatCry++)
		{
		  float sigmaNoiseOfThat=0;
		  if(treeVars.xtalInBCIEta[bCluster][thatCry]      !=-999999)   sigmaNoiseOfThat=sigmaNoiseEB;
		  else if (treeVars.xtalInBCIy[bCluster][thatCry]  !=-999999)   sigmaNoiseOfThat=sigmaNoiseEE;
		  else {std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}
		  
		  float ampliOfThat = treeVars.xtalInBCAmplitudeADC[bCluster][thatCry] / sigmaNoiseOfThat; 
		  if( ampliOfThat < minAmpliOverSigma) continue;
		  
		  float Aeff = ampliOfThis * ampliOfThat / sqrt( pow(ampliOfThis,2) + pow(ampliOfThat,2) );
		  float dt  = treeVars.xtalInBCTime[bCluster][thisCry] - treeVars.xtalInBCTime[bCluster][thatCry]; 
		  
		  dtVSAeffHistAnyPeak_  -> Fill(Aeff, dt); 
		  dtVSAeffProfAnyPeak_  -> Fill(Aeff, dt); 
		  if (thisIsInEBPeak) {
		    if      (Aeff < 6)   dtUpToQuarterGeVEBPeak->Fill(dt);
		    else if (Aeff < 12)  dtUpToHalfGeVEBPeak   ->Fill(dt);
		    else if (Aeff < 24)  dtUpToOneGeVEBPeak    ->Fill(dt);
		    else if (Aeff < 48)  dtUpToTwoGeVEBPeak    ->Fill(dt);
		    else                 dtUpOverTwoGeVEBPeak  ->Fill(dt);
		    
		    dtVSAeffHistEBPeak_ -> Fill(Aeff, fabs(dt)); 
		    dtVSAeffProfEBPeak_ -> Fill(Aeff, fabs(dt)); 
		  }
		  else      {
		    if      (Aeff < 6)   dtUpToThreeQuarterGeVEEPeak->Fill(dt);
		    else if (Aeff < 12)  dtUpToOneAndHalfGeVEEPeak  ->Fill(dt);
		    else if (Aeff < 24)  dtUpToThreeGeVEEPeak       ->Fill(dt);
		    else if (Aeff < 48)  dtUpToSixGeVEEPeak         ->Fill(dt);
		    else                 dtUpOverSixGeVEEPeak       ->Fill(dt);
		    
		    dtVSAeffHistEEPeak_ -> Fill(Aeff, fabs(dt)); 
		    dtVSAeffProfEEPeak_ -> Fill(Aeff, fabs(dt)); }
		  
		}// loop on thatCry
	      
	    }// loop on thisCry

	}// end loop on bc selected from pi0


    }   // end of SECOND loop over entries



  // now save the plots

  // write out control histograms
  TDirectory *controlPlots = saving.mkdir("control");
  controlPlots->cd();

  BCNumPerEventHist_->Write();
  BCEnergyHist_->Write();
  BCEtHist_->Write();
  BCNumCrysHist_->Write();
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
  massDiGammaEBHist_-> Write(); 
  massDiGammaEEHist_-> Write(); 
  massDiGammaEEPHist_-> Write(); 
  massDiGammaEEMHist_-> Write(); 

  diPhotonOccupancyAny_      -> Write();
  diPhotonPeakOccupancyAny_  -> Write();
  diPhotonSidesOccupancyAny_ -> Write();



  // write out single cluster resolution plots 
  TDirectory *singleClusResolution = saving.mkdir("single-resolution");
  singleClusResolution->cd();
  
  dtVSAeffHistAny_-> Write(); 
  dtVSAeffHistEB_ -> Write(); 
  dtVSAeffHistEE_ -> Write(); 
  dtVSAeffProfAny_-> Write(); 
  dtVSAeffProfEB_ -> Write(); 
  dtVSAeffProfEE_ -> Write(); 
  
  dtUpToQuarterGeVEB-> Write(); 
  dtUpToHalfGeVEB   -> Write(); 
  dtUpToOneGeVEB    -> Write(); 
  dtUpToTwoGeVEB    -> Write(); 
  dtUpOverTwoGeVEB  -> Write(); 
				    
  dtUpToThreeQuarterGeVEE-> Write(); 
  dtUpToOneAndHalfGeVEE  -> Write(); 
  dtUpToThreeGeVEE       -> Write(); 
  dtUpToSixGeVEE         -> Write(); 
  dtUpOverSixGeVEE       -> Write(); 


  // write out single cluster resolution plots after pi0 selection
  TDirectory *singleClusResolutionPi0Clusters = saving.mkdir("single-resolution-pi0cluster");
  singleClusResolutionPi0Clusters->cd();

  dtUpToQuarterGeVEBPeak -> Write(); 
  dtUpToHalfGeVEBPeak    -> Write(); 
  dtUpToOneGeVEBPeak     -> Write(); 
  dtUpToTwoGeVEBPeak     -> Write(); 
  dtUpOverTwoGeVEBPeak  -> Write(); 
				    
  dtUpToThreeQuarterGeVEEPeak -> Write(); 
  dtUpToOneAndHalfGeVEEPeak   -> Write(); 
  dtUpToThreeGeVEEPeak        -> Write(); 
  dtUpToSixGeVEEPeak          -> Write(); 
  dtUpOverSixGeVEEPeak        -> Write(); 
  
  dtVSAeffHistAnyPeak_ -> Write(); 
  dtVSAeffHistEBPeak_  -> Write(); 
  dtVSAeffHistEEPeak_  -> Write(); 
  dtVSAeffProfAnyPeak_ -> Write(); 
  dtVSAeffProfEBPeak_  -> Write(); 
  dtVSAeffProfEEPeak_  -> Write(); 




  saving.Close () ;

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
