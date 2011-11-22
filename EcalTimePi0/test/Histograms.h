#ifndef Histograms_H
#define Histograms_H

#include <TObject.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

using namespace std;

class hObj {

  public:

  hObj(){}
  hObj( string fsuffix ){}
  virtual ~hObj(){}
  virtual void Draw( string hfolder, string plotType ) = 0 ;
  //virtual void Fill1D( TString hName, double weight, double scale = 1. )  = 0;
  //virtual void Fill2D( TString hName, double weight, double scale = 1. )  = 0;
  //virtual vector<TH2D*> Output2D() = 0;

};

class hJetTime : public hObj {

  public:

  hJetTime(){

   nJets    = new TH1D("nJets",    " N of jets ", 11,  0, 11);
   nPhotons = new TH1D("nPhotons", " N of photons ", 11,  0, 11);
   nElectrons = new TH1D("nElectrons", " N of electrons ", 11,  0, 11);
   ele0Pt   = new TH1D("ele0Pt",   " leading electron Pt spectrum", 50,   0., 1000.);

   pho0Pt   = new TH1D("pho0Pt",   " leading photon Pt spectrum ", 50,   0., 1000.);
   pho0T    = new TH1D("pho0T",    " leading photon Time spectrum ", 160, -4., 4.);
   pho0xT   = new TH1D("pho0xT",   " leading photon Time spectrum from xtals ", 160, -4., 4.);
   pho0lT   = new TH1D("pho0lT",   " leading photon Time spectrum from leading xtals ", 160, -4., 4.);

   jet0Pt   = new TH1D("jet0Pt",   " leading jet Pt spectrum ", 50,   0., 1000.);
   jet0T    = new TH1D("jet0T",    " leading jet Time spectrum from BC", 160, -4., 4.);
   jet0xT   = new TH1D("jet0xT",   " leading jet Time spectrum from xtals", 160, -4., 4.);
   jet0lT   = new TH1D("jet0lT",   " leading jet Time spectrum from leading xtals", 160, -4., 4.);

   jet1Pt   = new TH1D("jet1Pt",   " 2nd jet Pt spectrum ", 50,   0., 1000.);
   jet1T    = new TH1D("jet1T",    " 2nd jet Time spectrum from BC", 160, -4., 4.);
   jet1xT   = new TH1D("jet1xT",   " 2nd jet Time spectrum from xtals", 160, -4., 4.);
   jet1lT   = new TH1D("jet1lT",   " 2nd jet Time spectrum from leading xtals", 160, -4., 4.);

   phoXtalTErr  = new TH1D("phoXtalTErr",  " xtal Time Error from photon", 200,   -1., 9.);
   phoXtalTime  = new TH1D("phoXtalTime",  " xtal Time       from photon", 80, -4., 4.);
   jetXtalTErr  = new TH1D("jetXtalTErr",  " xtal Time Error from jets", 200,   -1., 9.);
   jetXtalTime  = new TH1D("jetXtalTime",  " xtal Time       from jets", 80, -4., 4.);
   dT_JetPho    = new TH1D("dT_JetPho",    " T_pho - T_jet", 80, -4., 4. );

   jetXtalPos   = new TH2D("jetXtalPos",   " xtal Position   from jets", 125, -3.125, 3.125, 126, -3.15, 3.15 );

   jet0Pt_T     = new TH2D("jet0Pt_T",  " leading jet Pt vs Time", 40,   0., 800., 80,  -4., 4.) ;
   jet0Eta_T    = new TH2D("jet0Eta_T", " leading jet Eta vs Time", 41, -5.125, 5.125, 80,  -4., 4.) ;
   jet_emF_T    = new TH2D("jet_emF_T", " jet emF vs Time from leading 2 jets", 46, -0.1, 2.2, 80,  -4., 4.) ;
   jet_Pt_nXtal = new TH2D("jet_Pt_nXtal", "jet Pt vs N of xtals in a jet", 40, 0., 800., 51,  -0.5, 50.5) ;
   xE_cE        = new TH2D("xE_cE",     " xtal E vs cluster E from a jet ", 40, 0., 800, 40,  0., 800.) ;
   jet0Pt_ADC   = new TH2D("jet0Pt_ADC",  " leading jet Pt vs xtal ADC sum",  40,     0.,  800.,  50,  300., 6300.) ;
   jet0Eta_ADC  = new TH2D("jet0Eta_ADC", " leading jet Eta vs xtal ADC sum", 41, -5.125, 5.125,  50,  300., 6300.) ;
   jet0T_ADC    = new TH2D("jet0T_ADC", "xtal ADC sum vs leading jet xtal Time", 50,  300., 6300., 80, -4, 4. ) ;
   j0_xtalEBT_ADC = new TH2D("j0_xtalEBT_ADC", "EB xtal ADC vs xtal Time from leading jet", 100,  0., 2500., 80, -4, 4. ) ;
   j0_xtalEET_ADC = new TH2D("j0_xtalEET_ADC", "EE xtal ADC vs xtal Time from leading jet", 100,  0., 2500., 80, -4, 4. ) ;

   jPt_dR_gj    = new TH2D("jPt_dR_gj", " jet Pt vs min_dR( photon, jets )", 40, 0., 800., 51, -0.05, 5.05);
   ePt_dR_ge    = new TH2D("ePt_dR_ge", " electron Pt vs min_dR( photon, electrons )", 40, 0., 800., 51, -0.05, 5.05);

   c1 = new TCanvas("c1","", 800, 600);
   c2 = new TCanvas("c2","", 800, 600);
   c3 = new TCanvas("c3","", 800, 600);
   c4 = new TCanvas("c4","", 800, 600);

   leg2 = new TLegend(.8, .4, .97, .55 );
   leg3 = new TLegend(.15, .7, .35, .9 );

  }

  virtual ~hJetTime(){

    delete nJets ;
    delete nPhotons ;
    delete nElectrons ;
    delete ele0Pt ;

    delete pho0Pt ;
    delete pho0T ;
    delete pho0xT ;
    delete pho0lT ;

    delete jet0Pt ;
    delete jet0T ;
    delete jet0xT ;
    delete jet0lT ;

    delete jet1Pt ;
    delete jet1T ;
    delete jet1xT ;
    delete jet1lT ;

    delete phoXtalTErr ;
    delete phoXtalTime ;
    delete jetXtalTErr ;
    delete jetXtalTime ;
    delete dT_JetPho ;

    delete jetXtalPos ;
    delete jet0Pt_T ;    
    delete jet0Eta_T ;    
    delete jet_emF_T ;    
    delete jet_Pt_nXtal ;    
    delete xE_cE ;
    delete jet0Pt_ADC ;    
    delete jet0Eta_ADC ;    
    delete jet0T_ADC ;    
    delete j0_xtalEBT_ADC ;    
    delete j0_xtalEET_ADC ;    

    delete jPt_dR_gj ;
    delete ePt_dR_ge ;

    delete c1 ;
    delete c2 ;
    delete c3 ;
    delete c4 ;
   
    delete leg2;
    delete leg3;
  }

  void Draw( string hfolder, string plotType );

  TH1D* nJets ;
  TH1D* nPhotons ;
  TH1D* nElectrons ;
  TH1D* jet0Pt ;
  TH1D* jet1Pt ;
  TH1D* pho0Pt ;
  TH1D* ele0Pt ;

  TH1D* pho0T ;
  TH1D* jet0T ;
  TH1D* jet1T ;
  TH1D* pho0xT ;
  TH1D* jet0xT ;
  TH1D* jet1xT ;
  TH1D* pho0lT ;
  TH1D* jet0lT ;
  TH1D* jet1lT ;

  TH1D* phoXtalTErr ;
  TH1D* phoXtalTime ;
  TH1D* jetXtalTErr ;
  TH1D* jetXtalTime ;
  TH1D* dT_JetPho ;

  TH2D* jetXtalPos ;
  TH2D* jet0Pt_T ;
  TH2D* jet0Eta_T ;
  TH2D* jet_emF_T ;
  TH2D* jet_Pt_nXtal ;
  TH2D* xE_cE ;
  TH2D* jet0Pt_ADC ;
  TH2D* jet0Eta_ADC ;
  TH2D* jet0T_ADC ;
  TH2D* j0_xtalEBT_ADC ;
  TH2D* j0_xtalEET_ADC ;

  TH2D* jPt_dR_gj ;
  TH2D* ePt_dR_ge ;

  TCanvas* c1 ;
  TCanvas* c2 ;
  TCanvas* c3 ;
  TCanvas* c4 ;

  TLegend* leg2 ;
  TLegend* leg3 ;

  TString plotname1 ;
  TString plotname2 ;
  TString plotname3 ;
  TString plotname4 ;

  //ClassDef(hJetTime, 1);
  private:

};
#endif
