#include <iostream> 
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMinuit.h>

#include "MathFunctions.h"
#include "AnaInput.h"
#include "PhotonAna.h" 
#include "timeVsAmpliCorrector.h"

using namespace std; 

int main( int argc, const char* argv[] ) { 

  string datacardfile = ( argc > 1 ) ? argv[1] : "DataCard.txt";
  AnaInput        *Input = new AnaInput( datacardfile );
  PhotonAna   *photonAna = new PhotonAna( datacardfile );
  //vector<int> module;
  //Input->GetParameters( "Module", &module );

  Input->GetForest("TheData", "EcalTimeAnalysis");

  photonAna->ReadTree() ;
  photonAna->ScalarPlotList() ;
  return 0;

}

