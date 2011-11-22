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

using namespace std; 

int main() { 

  AnaInput        *Input = new AnaInput();
  PhotonAna   *photonAna = new PhotonAna();
  //vector<int> module;
  //Input->GetParameters( "Module", &module );

  Input->GetForest("TheData", "EcalTimeAnalysis");

  photonAna->ReadTree() ;
  photonAna->ScalarPlotList() ;
  return 0;

}

