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
   hMET     = new TH1D("hMET",     " MET spectrum", 50,   0., 500.);

   pho0Pt   = new TH1D("pho0Pt",   " leading photon Pt spectrum ", 50,   0., 1000.);
   pho0Eta  = new TH1D("pho0Eta",  " leading photon Eta ditribution ",  125, -3.125, 3.125);
   pho0xTB  = new TH1D("pho0xTB",  " leading Barrel photon Time spectrum from leading xtals in BC ", 250, -4., 21.);
   pho0lTB  = new TH1D("pho0lTB",  " leading photon Time spectrum from leading xtals in EB", 250, -4., 21.);
   pho0xTE  = new TH1D("pho0xTE",  " leading Endcap photon Time spectrum from leading xtals in BC", 250, -4., 21.);
   pho0lTE  = new TH1D("pho0lTE",  " leading photon Time spectrum from leading xtals in EE", 250, -4., 21.);

   jet0Pt   = new TH1D("jet0Pt",   " leading jet Pt spectrum ", 50,   0., 1000.);
   jet0Eta  = new TH1D("jet0Eta",  " leading jet Eta ditribution ",  125, -3.125, 3.125);
   jet0xTB  = new TH1D("jet0xTB",   " leading Barrel jet Time spectrum from leading xtals in BC", 250, -4., 21.);
   jet0lTB  = new TH1D("jet0lTB",   " leading jet Time spectrum from leading xtals in EB", 250, -4., 21.);
   jet0xTE  = new TH1D("jet0xTE",   " leading Endcap jet Time spectrum from leading xtals in BC", 250, -4., 21.);
   jet0lTE  = new TH1D("jet0lTE",   " leading jet Time spectrum from leading xtals in EE", 250, -4., 21.);

   elelTB  = new TH1D("elelTB",   " electron Time spectrum from leading xtals in EB", 250, -4., 21.);
   elelTE  = new TH1D("elelTE",   " electron Time spectrum from leading xtals in EE", 250, -4., 21.);

   jet1Pt   = new TH1D("jet1Pt",   " 2nd jet Pt spectrum ", 50,   0., 1000.);
   jet1xT   = new TH1D("jet1xT",   " 2nd jet Time spectrum from leading xtals in BC", 250, -4., 21.);
   jet1lT   = new TH1D("jet1lT",   " 2nd jet Time spectrum from leading xtals", 250, -4., 21.);

   phoXtalTErr  = new TH1D("phoXtalTErr",  " xtal Time Error from photon", 200,   -1., 9.);
   phoXtalTime  = new TH1D("phoXtalTime",  " xtal Time       from photon", 80, -4., 4.);
   jetXtalTErr  = new TH1D("jetXtalTErr",  " xtal Time Error from jets", 200,   -1., 9.);
   jetXtalTime  = new TH1D("jetXtalTime",  " xtal Time       from jets", 80, -4., 4.);
   dT_JetPho    = new TH1D("dT_JetPho",    " T_photon - T_jet", 125, -4., 21. );

   xtalEB_Chi2_T  = new TH2D("xtalEB_Chi2_T", "EB xtal Chi2 vs xtal Time", 90, 0, 90, 250, -4., 21. );
   xtalEE_Chi2_T  = new TH2D("xtalEE_Chi2_T", "EE xtal Chi2 vs xtal Time", 90, 0, 90, 250, -4., 21. );
   xtal_Chi2_E  = new TH2D("xtal_Chi2_E", "xtal Chi2 vs xtal Energy", 90, 0, 90, 160, 0., 800. );
   xtal_T_E     = new TH2D("xtal_T_E", "xtal Time vs xtal Energy",  250, -4., 21., 160, 0., 800. );
   jetXtalPos   = new TH2D("jetXtalPos",   " xtal Position   from jets", 125, -3.125, 3.125, 126, -3.15, 3.15 );

   jet0Pt_T     = new TH2D("jet0Pt_T",  " leading jet Pt vs Time", 40,   0., 800., 80,  -4., 4.) ;
   jet0Eta_T    = new TH2D("jet0Eta_T", " leading jet Eta vs Time", 41, -5.125, 5.125, 80,  -4., 4.) ;
   jet_emF_T    = new TH2D("jet_emF_T", " jet emF vs Time from leading 2 jets", 46, -0.1, 2.2, 80,  -4., 4.) ;
   jet_Pt_nXtal = new TH2D("jet_Pt_nXtal", "jet Pt vs N of xtals in a jet", 40, 0., 800., 51,  -0.5, 50.5) ;
   xE_cE        = new TH2D("xE_cE",     " xtal E vs cluster E from a jet ", 40, 0., 800, 40,  0., 800.) ;

   g0_xtalEta_ADC  = new TH2D("g0_xtalEta_ADC",  "xtal Eta vs xtal ADC from leading photon", 41, -5.125, 5.125, 100,  0., 1000.) ;
   g0_xtalEB_ADC_T = new TH2D("g0_xtalEB_ADC_T", "EB xtal ADC vs xtal Time from leading photon", 100,  0., 2500., 80, -4, 4. ) ;
   g0_xtalEB_ADC_E = new TH2D("g0_xtalEB_ADC_E", "EB xtal ADC vs xtal Energy from leading photon", 100,  0., 2500., 100, 0, 200. ) ;

   j0_xtalEta_ADC  = new TH2D("j0_xtalEta_ADC",  "xtal Eta vs xtal ADC", 41, -5.125, 5.125, 100,  0., 1000.) ;
   j0_xtalEB_ADC_T = new TH2D("j0_xtalEB_ADC_T", "EB xtal ADC vs xtal Time from leading jet", 100,  0., 2500., 80, -4, 4. ) ;
   j0_xtalEE_ADC_T = new TH2D("j0_xtalEE_ADC_T", "EE xtal ADC vs xtal Time from leading jet", 100,  0., 2500., 80, -4, 4. ) ;
   j0_xtalEB_ADC_E = new TH2D("j0_xtalEB_ADC_E", "EB xtal ADC vs xtal Energy from leading jet", 100,  0., 2500., 100, 0, 200. ) ;
   j0_xtalEE_ADC_E = new TH2D("j0_xtalEE_ADC_E", "EE xtal ADC vs xtal Energy from leading jet", 100,  0., 2500., 100, 0, 200. ) ;

   jPt_dR_gj    = new TH2D("jPt_dR_gj", " jet Pt vs min_dR( photon, jets )", 40, 0., 800., 51, -0.05, 5.05);
   ePt_dR_ge    = new TH2D("ePt_dR_ge", " electron Pt vs min_dR( photon, electrons )", 40, 0., 800., 51, -0.05, 5.05);

   c1 = new TCanvas("c1","", 800, 600);
   c2 = new TCanvas("c2","", 800, 600);
   c3 = new TCanvas("c3","", 800, 600);
   c4 = new TCanvas("c4","", 800, 600);
   c5 = new TCanvas("c5","", 800, 600);

   leg2 = new TLegend(.8, .4, .97, .55 );
   leg3 = new TLegend(.42, .75, .6, .9 );

  }

  virtual ~hJetTime(){

    delete nJets ;
    delete nPhotons ;
    delete nElectrons ;
    delete ele0Pt ;
    delete hMET ;

    delete pho0Pt ;
    delete pho0Eta ;
    delete pho0xTB ;
    delete pho0lTB ;
    delete pho0xTE ;
    delete pho0lTE ;

    delete jet0Pt ;
    delete jet0Eta ;
    delete jet0xTB ;
    delete jet0lTB ;
    delete jet0xTE ;
    delete jet0lTE ;

    delete jet1Pt ;
    delete jet1xT ;
    delete jet1lT ;

    delete elelTB ;
    delete elelTE ;

    delete phoXtalTErr ;
    delete phoXtalTime ;
    delete jetXtalTErr ;
    delete jetXtalTime ;
    delete dT_JetPho ;

    delete xtalEB_Chi2_T;
    delete xtalEE_Chi2_T;
    delete xtal_Chi2_E;
    delete xtal_T_E;
    delete jetXtalPos ;
    delete jet0Pt_T ;    
    delete jet0Eta_T ;    
    delete jet_emF_T ;    
    delete jet_Pt_nXtal ;    
    delete xE_cE ;

    delete g0_xtalEta_ADC ;    
    delete g0_xtalEB_ADC_T ;    
    delete g0_xtalEB_ADC_E ;    
    delete j0_xtalEta_ADC ;    
    delete j0_xtalEB_ADC_T ;    
    delete j0_xtalEE_ADC_T ;    
    delete j0_xtalEB_ADC_E ;    
    delete j0_xtalEE_ADC_E ;    

    delete jPt_dR_gj ;
    delete ePt_dR_ge ;

    delete c1 ;
    delete c2 ;
    delete c3 ;
    delete c4 ;
    delete c5 ;
   
    delete leg2;
    delete leg3;
  }

  void Draw( string hfolder, string plotType );
  void FitNDraw( string hfolder, string plotType );

  TH1D* nJets ;
  TH1D* nPhotons ;
  TH1D* nElectrons ;
  TH1D* jet0Pt ;
  TH1D* jet1Pt ;
  TH1D* pho0Pt ;
  TH1D* ele0Pt ;
  TH1D* hMET ;
  TH1D* jet0Eta ;
  TH1D* pho0Eta ;

  TH1D* pho0xTB ;
  TH1D* pho0lTB ;
  TH1D* pho0xTE ;
  TH1D* pho0lTE ;
  TH1D* jet0xTB ;
  TH1D* jet0lTB ;
  TH1D* jet0xTE ;
  TH1D* jet0lTE ;
  TH1D* jet1xT ;
  TH1D* jet1lT ;

  TH1D* elelTB;
  TH1D* elelTE;

  TH1D* phoXtalTErr ;
  TH1D* phoXtalTime ;
  TH1D* jetXtalTErr ;
  TH1D* jetXtalTime ;
  TH1D* dT_JetPho ;

  TH2D* xtalEB_Chi2_T ;
  TH2D* xtalEE_Chi2_T ;
  TH2D* xtal_Chi2_E ;
  TH2D* xtal_T_E ;
  TH2D* jetXtalPos ;
  TH2D* jet0Pt_T ;
  TH2D* jet0Eta_T ;
  TH2D* jet_emF_T ;
  TH2D* jet_Pt_nXtal ;
  TH2D* xE_cE ;
  TH2D* g0_xtalEta_ADC ;
  TH2D* g0_xtalEB_ADC_T ;
  TH2D* g0_xtalEB_ADC_E ;
  TH2D* j0_xtalEta_ADC ;
  TH2D* j0_xtalEB_ADC_T ;
  TH2D* j0_xtalEE_ADC_T ;
  TH2D* j0_xtalEB_ADC_E ;
  TH2D* j0_xtalEE_ADC_E ;

  TH2D* jPt_dR_gj ;
  TH2D* ePt_dR_ge ;

  TCanvas* c1 ;
  TCanvas* c2 ;
  TCanvas* c3 ;
  TCanvas* c4 ;
  TCanvas* c5 ;

  TLegend* leg2 ;
  TLegend* leg3 ;

  TString plotname1 ;
  TString plotname2 ;
  TString plotname3 ;
  TString plotname4 ;
  TString plotname5 ;

  //ClassDef(hJetTime, 1);
  private:

};
#endif
