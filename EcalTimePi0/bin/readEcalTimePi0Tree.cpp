#include <iostream>
#include <string>

#include <map>
#include <vector>
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
#include "TFile.h"


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
  std::string stringETGammaMin       = "--eTGammaMin";
  std::string strings4s9GammaMin     = "--s4s9GammaMin";
  std::string stringeTPi0Min         = "--stringeTPi0Min";
  std::string stringNumEvents        = "--n";

  std::vector<std::string> listOfFiles;
  int   numEvents    =-1;
  float	eTGammaMin   =0.2;
  float s4s9GammaMin = 0.85;
  float eTPi0Min     = 0.65;

  //gf: support development
  //std::cout << "\nargc:       " << argc << std::endl;
  //for (int v=0; v<argc; v++ ){      std::cout << "argument: " << v << " argv: " << argv[v] << std::endl;    }

  
  // if no arguments are passed, suggest help
  if (argc < 2){
    std::cerr << "\n\tERROR: specify arguments, e.g. --help\n" << std::endl ;
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
	std::cout << " --eTGammaMin: min eT for gammas" << std::endl;
	std::cout << " --s4s9GammaMin: min shower shape" << std::endl;
	std::cout << " --eTPi0Min min eT for pi0 candidate" << std::endl;
	std::cout << " --i <list of strings> list of input files" << std::endl ;     
	exit(1);      }

      
      else if (argv[v] == stringNumEvents) { // set number of events
	std::cout << "events number" << std::endl;
	numEvents=atoi(argv[v+1]);
	v++;
      }

      
      else if (argv[v] == stringETGammaMin) { // choose et cut for single cluster
	eTGammaMin = atof(argv[v+1]);
	v++;
      }
      
      else if (argv[v] == strings4s9GammaMin) { // choose cut for shower shape
	s4s9GammaMin = atof(argv[v+1]);
	v++;
      }
      
      else if (argv[v] == stringeTPi0Min) { // choose et cut for pi0 candidate
	eTPi0Min = atof(argv[v+1]);
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
  std::cout << "\n\tFOUND " << nEntries << " events" << std::endl ;    
  std::cout << "\tWILL run on: " <<  numEvents << " events" << std::endl;
  std::cout << "\tOutput file: " <<  outputRootName << std::endl;
  std::cout << "\teTGammaMin: "  <<  eTGammaMin << std::endl;
  std::cout << "\ts4s9GammaMin: "  <<  s4s9GammaMin << std::endl;

  EcalTimePi0TreeContent treeVars ; 
  setBranchAddresses (chain, treeVars) ;

  // Initialize output root file
  TFile saving (outputRootName.c_str (),"recreate") ;
  saving.cd () ;
  // Initialize histograms -- BasicClusters
  TH1F* BCNumPerEventHist_ = new TH1F("BCNumPerEvent","Number of BC per event",50,0,49);
  TH1F* BCEnergyHist_ = new TH1F("BCEnergy","Energy of BCs;GeV",100,0,24);
  TH1F* BCEtHist_ = new TH1F("BCEt","E_{T} of BCs;GeV",100,0,24);
  TH1F* BCNumCrysHist_ = new TH1F("BCNumCrys","Number of crystals per BC",10,0,10);
  TH2F* BCOccupancyHistEB_ = new TH2F("BCOccupancyEB","BC occupancy;i#phi;i#eta",360,1.,361.,172,-86,86);
  // Initialize histograms -- xtals
  TH1F* xtalEnergyHist_ = new TH1F("XtalEnergy","Crystal energy;GeV",110,-1,10);
  TH1F* xtalTimeHist_ = new TH1F("XtalTime","CrystalTime;ns",60,-25,35);
  // TH1F* xtalStatusHist_ = new TH1F("XtalStatus","Crystal status flag",16,0,15);
  // TH2F* xtalOccupancyHistEB_ = new TH2F("XtalOccupancyEB","Crystal occupancy;i#phi;i#eta",360,1.,361.,172,-86,86);
  // Initialize histograms -- diphotons
  TH1F* massDiGammaHist_ = new TH1F("massDiGamma","m(#gamma#gamma)",50,0,0.500);
  



  //loop over entries
  for (int entry = 0 ; (entry < nEntries && entry < numEvents); ++entry)
    {
      chain->GetEntry (entry) ;

      bool speak=false;
      if (entry<10 || entry%100==0) speak=true;

      if (speak) std::cout << "------> reading entry " << entry << " <------\n" ; 


      // loop on calorimetric quantities

      if (speak)  std::cout << "  found " << treeVars.nSuperClusters << " superclusters" << std::endl ;
      if (speak)  std::cout << "  found " << treeVars.nClusters << " basic clusters" << std::endl ;
      BCNumPerEventHist_->Fill(treeVars.nClusters);


        /////////////////////////////////////////////////////
	//loop on basic cluster

	for (int bCluster=0; bCluster < treeVars.nClusters; bCluster++)
	  {

	    float eBC=0; // calculate energy of BC for validation
	    for (int cryInBC=0; cryInBC < treeVars.nXtalsInCluster[bCluster]; cryInBC++){
	      eBC+= treeVars.xtalInBCEnergy[bCluster][cryInBC];}

            BCEnergyHist_->Fill(treeVars.clusterEnergy[bCluster]);
            BCEtHist_->Fill(treeVars.clusterTransverseEnergy[bCluster]);
            BCNumCrysHist_->Fill(treeVars.nXtalsInCluster[bCluster]);
            if(treeVars.xtalInBCIEta[bCluster][0] != 0)
              BCOccupancyHistEB_->Fill(treeVars.xtalInBCIPhi[bCluster][0],treeVars.xtalInBCIEta[bCluster][0]);  // first cry for now

	    if (speak)  std::cout << "\tbCluster: num"               << bCluster 
				  << "\t eBC: "                      << treeVars.clusterEnergy[bCluster]
				  << "\t eBC_predicted: "            << eBC
				  << "\n\t et: "                     << treeVars.clusterTransverseEnergy[bCluster]
				  << "\t predicted et: "             << treeVars.clusterEnergy[bCluster]*sin(2*atan(exp(-1* treeVars.clusterEta[bCluster] )) )
				  << " eta: "                        << treeVars.clusterEta[bCluster]
				  << "\n\t num crystals: "           << treeVars.nXtalsInCluster[bCluster]
				  << "\n\t\tfirst crystal:  \tieta " << treeVars.xtalInBCIEta[bCluster][0] 
				  << "\teta "                        << treeVars.xtalInBCEta[bCluster][0] 
				  << " \t energy "                   << treeVars.xtalInBCEnergy[bCluster][0] 
				  << " \t ADC "                      << treeVars.xtalInBCAmplitudeADC[bCluster][0] 
				  << " \t time "                     << treeVars.xtalInBCTime[bCluster][0] 
				  << std::endl;
	  }

	if (speak) std::cout << "  found " << treeVars.nXtals << " crystals\n" ;    
      //PG loop over crystals
      for (int XTLindex = 0 ; XTLindex < treeVars.nXtals ; ++XTLindex)
        {
          //TODO: For this to work, we need another way to tell EB from EE,
          //      as hashedIndex itself is not enough!
          //if(EBDetId::validHashIndex(treeVars.xtalHashedIndex[XTLindex]))
          //{
          //  EBDetId ebDet = EBDetId::unhashIndex(treeVars.xtalHashedIndex[XTLindex]);   
          //  xtalOccupancyHistEB_->Fill(ebDet.iphi(),ebDet.ieta());
          //}
          //else if(EEDetId::validHashIndex(treeVars.xtalHashedIndex[XTLindex]))
          //  EEDetId eeDet = EEDetId::unhashIndex(treeVars.xtalHashedIndex[XTLindex]);
          //else
          //  std::cout << "Crystal with hash: " << 

          xtalEnergyHist_->Fill(treeVars.xtalEnergy[XTLindex]);
          xtalTimeHist_->Fill(treeVars.xtalTime[XTLindex]);
          //TODO
          //xtalStatusHist_->Fill

        } //PG loop over crystals






      for (int bClusterA=0; bClusterA < treeVars.nClusters; bClusterA++)
	{

	  float eTA = treeVars.clusterTransverseEnergy[bClusterA];
	  if( eTA < eTGammaMin ) continue;
	  
	  float e22A = treeVars.clusterE2x2[bClusterA];
	  float e33A = treeVars.clusterE3x3[bClusterA];
	  if ( e22A/e33A < s4s9GammaMin ) continue;
        
	  for (int bClusterB=(bClusterA+1); bClusterB < treeVars.nClusters; bClusterB++)
	    {

	      float eTB = treeVars.clusterTransverseEnergy[bClusterB];
	      if( eTB < eTGammaMin) continue;

	      float e22B = treeVars.clusterE2x2[bClusterB];
	      float e33B = treeVars.clusterE3x3[bClusterB];
	      if ( e22B/e33B < s4s9GammaMin ) continue;
	      
	      math::PtEtaPhiMLorentzVectorD gammaA (eTA, treeVars.clusterEta[bClusterA], treeVars.clusterPhi[bClusterA], 0);
	      math::PtEtaPhiMLorentzVectorD gammaB (eTB, treeVars.clusterEta[bClusterB], treeVars.clusterPhi[bClusterB], 0);
	      
	      math::PtEtaPhiMLorentzVectorD pi0Candidate = gammaA + gammaB;

	      std::cout << gammaA << " " << gammaA.M() << std::endl;
	      std::cout << gammaB << " " << gammaB.M() << std::endl;
	      std::cout << pi0Candidate << " " << pi0Candidate.M() << std::endl;
	      
	      if(pi0Candidate.Et() < eTPi0Min) continue;

	      massDiGammaHist_ -> Fill(pi0Candidate.M());

	    }//loop on candidateB
	}//loop on candidateA




    } //PG loop over entries


  BCNumPerEventHist_->Write();
  BCEnergyHist_->Write();
  BCEtHist_->Write();
  BCNumCrysHist_->Write();
  BCOccupancyHistEB_->Write();
  xtalEnergyHist_->Write(); 
  xtalTimeHist_->Write();
  // xtalStatusHist_->Write();
  // xtalOccupancyHistEB_->Write();
  massDiGammaHist_-> Write(); 

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
