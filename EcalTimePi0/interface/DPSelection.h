#ifndef DPSelection_H
#define DPSelection_H

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

#include "ECALTime/EcalTimePi0/interface/AnaInput.h"


#define MAXSC 50
#define MAXC 200
#define MAXXTALINC 25 
// defined variables to follow are by default set up in ECALTime/EcalTimePi0/interface/EcalTimePhyTreeContent.h
#ifndef MAXOBJ
#define MAXOBJ 10
#endif
#ifndef MAXPHO
#define MAXPHO 40
#endif
#ifndef MAXVTX
#define MAXVTX 40
#endif


typedef pair<int, TLorentzVector> objID ;


//class DPSelection : public TObject {
class DPSelection {

public:

   DPSelection( string datacardfile = "DataCard.txt");     
   ~DPSelection();     
   //friend class Sync ;   

   void Init( TTree* tr ) ;

   bool HLTFilter();
   bool PhotonFilter( bool doIso = true );
   bool JetMETFilter( bool useJetID = false );
   bool VertexFilter();
   bool ElectronFilter();
   bool MuonFilter();

   bool GammaJetsBackground() ; 

   bool SignalSelection( bool isTightPhoton = true ) ; 
   bool GammaJetsControlSample( bool isTightPhoton = true ) ; 
   bool QCDControlSample() ; 

   void ResetCuts( string cutName, vector<int>& cutId, vector<double>& newValue ) ;
   void ResetCuts( string cutName, int cutId, double newValue ) ;
   void ResetCuts( string cutName = "All" ) ; // set cuts to default values from datacard 

   void GetCollection( string collName, vector<objID>& coll ) ;
   void ResetCollection( string cutName = "All" ) ; // clean the storage containers 

private:

   AnaInput*       Input;

   vector<objID> phoV ;
   vector<objID> jetV ;
   vector<objID> eleV ;
   vector<objID> muV ;

   vector<double> photonCuts ;
   vector<double> photonIso ;
   vector<double> vtxCuts ;
   vector<double> jetCuts ;
   vector<double> electronCuts ;
   vector<double> muonCuts ;

   unsigned int eventId ;
   float vtxX[MAXVTX],    vtxY[MAXVTX],  vtxZ[MAXVTX],   vtxChi2[MAXVTX], vtxNdof[MAXVTX];
   float jetPx[MAXOBJ],   jetPy[MAXOBJ], jetPz[MAXOBJ],  jetE[MAXOBJ] ;
   float jetNDau[MAXOBJ], jetCM[MAXOBJ], jetCEF[MAXOBJ], jetNHF[MAXOBJ], jetNEF[MAXOBJ];
   float phoPx[MAXPHO],      phoPy[MAXPHO],      phoPz[MAXPHO],     phoE[MAXPHO], phoTime[MAXPHO] ;
   float phoEcalIso[MAXPHO], phoHcalIso[MAXPHO], phoTrkIso[MAXPHO], phoHovE[MAXPHO], phoSmin[MAXPHO], phoSmaj[MAXPHO] ;
   float elePx[MAXOBJ], elePy[MAXOBJ], elePz[MAXOBJ], eleE[MAXOBJ] ;
   float eleEcalIso[MAXOBJ], eleHcalIso[MAXOBJ], eleTrkIso[MAXOBJ] ;
   int   eleNLostHits[MAXOBJ] ;
   float muPx[MAXOBJ], muPy[MAXOBJ], muPz[MAXOBJ], muE[MAXOBJ] ;
   float muEcalIso[MAXOBJ], muHcalIso[MAXOBJ], muTrkIso[MAXOBJ] ;
   float metPx, metPy, metE ;
   int   nJets, nPhotons, nElectrons, nVertices, nMuons, trgCut, nSuperClusters, nClusters ;
   /*
   float superClusterPhi[MAXSC], superClusterEta[MAXSC] ;
   float CPIdx[MAXC], clusterTime[MAXC], clusterEnergy[MAXC];
   float xtalInBCEnergy[MAXC][MAXXTALINC], xtalInBCTime[MAXC][MAXXTALINC], xtalInBCTimeErr[MAXC][MAXXTALINC];
   float xtalInBCEta[MAXC][MAXXTALINC],    xtalInBCPhi[MAXC][MAXXTALINC],  xtalADC[MAXC][MAXXTALINC] ;
   float xtalChi2[MAXC][MAXXTALINC], xtalOutTimeChi2[MAXC][MAXXTALINC];
   float clusterPhi[MAXC], clusterEta[MAXC] ;
   int   nXtalInBC[MAXC], clusterXtals[MAXC] ;
   */
   //ClassDef(DPSelection, 1);
};

//#if !defined(__CINT__)
//    ClassImp(DPSelection);
#endif

