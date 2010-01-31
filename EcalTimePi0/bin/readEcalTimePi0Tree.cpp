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
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"

#define BarrelLimit 1.479
#define EndcapLimit 3.0

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

  float sigmaNoiseEB     = 1.06; // ADC 
  float sigmaNoiseEE     = 2.10; // ADC
  float minApliOverSigma = 7;    // dimensionless
  
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
  std::cout << "\n\tFOUND "         << nEntries << " events" << std::endl ;    
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

  // Initialize histograms -- diphotons
  TH1F* massDiGammaHist_ = new TH1F("massDiGamma","m(#gamma#gamma)",50,0,0.500);
  TH1F* massDiGammaEBHist_ = new TH1F("massDiGamma EB","m(#gamma#gamma) EB",50,0,0.500);
  TH1F* massDiGammaEEHist_ = new TH1F("massDiGamma EE","m(#gamma#gamma) EE",50,0,0.500);
  TH1F* massDiGammaEEPHist_ = new TH1F("massDiGamma EEP","m(#gamma#gamma) EEP",50,0,0.500);
  TH1F* massDiGammaEEMHist_ = new TH1F("massDiGamma EEM","m(#gamma#gamma) EEM",50,0,0.500);

  // Initialize histograms -- single cluster resolution
  TH2F* dtVSAeffHistAny_ = new TH2F("#delta(t) VS A_{eff}","#delta(t) VS A_{eff}",250,0.,250,100,0,10);
  TH2F* dtVSAeffHistEB_  = new TH2F("EB: #delta(t) VS A_{eff}","EB #delta(t) VS A_{eff}",250,0.,250,100,0,10);
  TH2F* dtVSAeffHistEE_  = new TH2F("EE: E#delta(t) VS A_{eff}","EE: #delta(t) VS A_{eff}",250,0.,250,100,0,10);
  TProfile* dtVSAeffProfAny_ = new TProfile("#delta(t) VS A_{eff} prof","#delta(t) VS A_{eff} prof",250,0.,250,0,10);
  TProfile* dtVSAeffProfEB_  = new TProfile("EB: #delta(t) VS A_{eff} prof","EB #delta(t) VS A_{eff} prof",250,0.,250,0,10);
  TProfile* dtVSAeffProfEE_  = new TProfile("EE: E#delta(t) VS A_{eff} prof","EE: #delta(t) VS A_{eff} prof",250,0.,250,0,10);



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
		eTPi0Min = eTPi0MinEB;	      }
	      else{		eTPi0Min = eTPi0MinEE;
	      }


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
  



  /////////////////////////////////////////////////////
  // SECOND loop over entries
  for (int entry = 0 ; (entry < nEntries && entry < numEvents); ++entry)
    {
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
		if( ampliOfThis < minApliOverSigma) continue;


		for(int thatCry=thisCry+1; thatCry<treeVars.nXtalsInCluster[bCluster]; thatCry++)
		  {
		    float sigmaNoiseOfThat=0;
		    if(treeVars.xtalInBCIEta[bCluster][thatCry]      !=-999999)   sigmaNoiseOfThat=sigmaNoiseEB;
		    else if (treeVars.xtalInBCIy[bCluster][thatCry]  !=-999999)   sigmaNoiseOfThat=sigmaNoiseEE;
		    else {std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}
		    
		    float ampliOfThat = treeVars.xtalInBCAmplitudeADC[bCluster][thatCry] / sigmaNoiseOfThat; 
		    if( ampliOfThat < minApliOverSigma) continue;
		    
		    float Aeff = ampliOfThis * ampliOfThat / sqrt( pow(ampliOfThis,2) + pow(ampliOfThat,2) );
		    float dt  = treeVars.xtalInBCTime[bCluster][thisCry] - treeVars.xtalInBCTime[bCluster][thatCry]; 
		    dt        = fabs(dt);

		    // for debug
		    std::cout << "ampliOfThis: " << ampliOfThis << "\tampliOfThat: " << ampliOfThat
			      << "\n timeOfThis: " << treeVars.xtalInBCTime[bCluster][thisCry] << "\ttimeOfThat: " << treeVars.xtalInBCTime[bCluster][thatCry]
			      << "\n Aeff: " << Aeff << "\tdt: " << dt 
			      << std::endl;
		    

		    dtVSAeffHistAny_  -> Fill(Aeff, dt); 
		    dtVSAeffProfAny_  -> Fill(Aeff, dt); 
		    if (thisIsInEB) {
		      dtVSAeffHistEB_ -> Fill(Aeff, dt); 
		      dtVSAeffProfEB_ -> Fill(Aeff, dt); }
		    else      {
		      dtVSAeffHistEE_ -> Fill(Aeff, dt); 
		      dtVSAeffProfEE_ -> Fill(Aeff, dt); }
		    
		  }// loop on thatCry
		
	      }// loop on thisCry
	  }// end loop on bc
	
	

	
      float eTA, eTB ;
      float e22A, e33A,    e22B, e33B;
      float eTGammaMinA,   eTGammaMinB;
      float s4s9GammaMinA, s4s9GammaMinB;
      bool  AisEB,         BisEB;
      float eTPi0Min;
      
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
	      

	      if ( fabs(pi0Candidate.Eta()) < BarrelLimit) {
		eTPi0Min = eTPi0MinEB;	      }
	      else{		eTPi0Min = eTPi0MinEE;
	      }


	      // third selection cut: pi0 candidate Et
	      if(pi0Candidate.Et() < eTPi0Min ) continue;

	      // here I have pi0's that pass the definition   


	      
	    }//loop on candidateB
	}//loop on candidateA
      

    }   // end of first loop over entries
        // now you have di-mass plots filled => get masses































  // write out histograms
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
  massDiGammaHist_-> Write(); 
  massDiGammaEBHist_-> Write(); 
  massDiGammaEEHist_-> Write(); 
  massDiGammaEEPHist_-> Write(); 
  massDiGammaEEMHist_-> Write(); 
  
  dtVSAeffHistAny_-> Write(); 
  dtVSAeffHistEB_ -> Write(); 
  dtVSAeffHistEE_ -> Write(); 
  dtVSAeffProfAny_-> Write(); 
  dtVSAeffProfEB_ -> Write(); 
  dtVSAeffProfEE_ -> Write(); 


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
