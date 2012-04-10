#ifndef PhotonAna_H
#define PhotonAna_H

#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TLorentzVector.h"

#include "AnaInput.h"
#include "Histograms.h"
#include "MathFunctions.h"
#include "timeVsAmpliCorrector.h"

#define MAXSC 50
#define MAXC 200
#define MAXXTALINC 25 
#define MAXOBJ 10
#define MAXVTX 40

class PhotonAna : public TObject {

public:

   PhotonAna( string datacardfile = "DataCard.txt");     
   ~PhotonAna();     
   
   void ReadTree();
   
   vector<double> BinningFit( TH2D* h2, string hName, int xbinMin, int xbinMax, bool debugPlots = false );
   void BinningFitScan( TH2D* h2, vector<double>& xV, vector<double>& yV, vector<double>& yWidthV, vector<double>& yErrV, 
                        bool debugPlots = false,  int rbin = 1, int startBin = 1, int finalBin = -1) ;

   // drawOpt = 0 : default , 1:SetLogx , 2:SetLogy 
   void ScalarPlotter( TH2D* h2, TString hname, double yMin, double yMax, int rbin = 1, int drawOpt = 0, bool debugPlots = false );
   void ScalarPlotList() ;

   double DeltaR( float eta1, float phi1, float eta2, float phi2 ) ;

private:

   AnaInput*       Input;

   string hfolder ;
   string plotType ;
   int ProcessEvents ;

   hJetTime hJets ;
   string debugStr ;

   vector<TLorentzVector> jetV ;
   vector<TLorentzVector> eleV ;
   vector<float> eIdV ;

   vector<double> photonCuts ;
   vector<double> electronCuts ;
   vector<double> jetCuts ;
   vector<double> BasicClusterCuts ;
   vector<double> XtalCuts ;
   int selectBackground ;
   int split ;
   int doTimeCorrection ;

   timeCorrector theTimeCorrector_;


   //ClassDef(PhotonAna, 1);
};

//#if !defined(__CINT__)
//    ClassImp(PhotonAna);
#endif

