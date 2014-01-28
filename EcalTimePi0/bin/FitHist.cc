//** My simple Macro to Fit Z mass,time, and Plots, hists ***/
//Run AS 
//>> root -l FitHist.cc
// OR 
//>>> ./compile 
//>>>./HFit
// Designed by 10Sr @2013 // norbe072@umn.edu
#include "TObject.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Riostream.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TGaxis.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include  "TMath.h"
#include  "TPaveText.h"
#include  "TVirtualFitter.h"
#include  "TMatrixT.h"
#include  "TMatrixD.h"
#include  "TROOT.h"
#include  <iostream>
#include  <string>
#include  <cmath>


//** Define Fit Constants **//
#define FitLowRange (-2)
#define FitHighRange (2)
#define LRange (-3)
#define HRange (3)
#define Title "Seed Time[ns]"
#define YTitle "numb. of seeds/0.05ns"
#define XTitle "t_{seed}[ns]"
#define FITTYPE "GAUS"
#define DET "EB"
#define HIST "EBEB/seed time"
#define CNAME "EB_EB_Seed_Time_DoubleElectron_Run2012A.png" 


using namespace std;
/*
const signed FitLowRange  = -3; 
const signed FitHighRange = 3;
const signed LRange = -3; 
const signed HRange = 3; 
*/
/*** 2D Gaus Fit ***/
double gauss2D(double *x, double *par) {
	double z1 = 0;
	if(par[2] != 0) z1 = double((x[0]-par[1])/par[2]);
	double z2 = 0;
	if(par[4] != 0 ) z2 = double((x[1]-par[3])/par[4]);          
  
    return par[0]*exp(-0.5*(z1*z1+z2*z2));
} 

/**private Defined  gaus fxn **/
double mygaus(double *t, double *par){
	double x = t[0];
	double arg0 = par[0];
	double arg1 = 0;
	if(par[2] != 0) arg1 =double( ( x - par[1])/par[2] );
	double  f   = arg0*TMath::Exp(-0.5*arg1*arg1);
	
    return f;
}



//## Style of Hist ###/////////
void HistS( TH1F* i_hist ){

   i_hist->SetTitle(Title);   
   i_hist->SetMarkerStyle(20);
   i_hist->SetMarkerSize(1.0);
   i_hist->SetMarkerColor(1);
   i_hist->SetLineStyle(1);
   i_hist->SetLineWidth(3);
   i_hist->SetLineColor(1);
   i_hist->SetStats(1);
   i_hist->SetTitleSize(0.08, "x");   
   i_hist->SetTitleOffset(1.0, "x");    
   i_hist->SetTitleSize(0.06, "y"); 
   i_hist->SetTitleOffset(0.95, "y");    
   i_hist->SetYTitle(YTitle); 
   i_hist->SetXTitle(XTitle); 
   i_hist->GetXaxis()->SetRangeUser(LRange, HRange);   
   

}

////### Fit Function ##### ////
void FfxnS( TF1* fxnfit){
   fxnfit->SetNpx(500);
   fxnfit->SetLineWidth(4);
   fxnfit->SetLineStyle(5);
   fxnfit->SetLineColor(kBlue);

}

////### Canvas Style ###////
void CanS ( TCanvas* ct ){
   ct->SetGridx();
   ct->SetGridy();
   ct->GetFrame()->SetFillColor(21);
   ct->GetFrame()->SetBorderMode(-1);
   ct->GetFrame()->SetBorderSize(5);
/* c1->Divide(2,1);  */
}

///### Legend Style ### ///

void LegS ( TLegend* lg, TF1* ffxn, TH1F* hf ){
//  draw the legend
   lg->SetTextFont(72);
   lg->SetTextSize(0.04);
   lg->AddEntry(hf,DET,"lpe");
   lg->AddEntry(ffxn,FITTYPE,"l");

}


//***############## main fitting Fxn ################ *****//
void Fit_Hist( TH1F* ihist, TF1* fitfxn, TCanvas*c1, TLegend *leg ){
            
   HistS( ihist );
   FfxnS( fitfxn );
   /// Set parms as parms of Fit Fxn///
   fitfxn->SetParameters(500, ihist->GetMean(), ihist->GetRMS() );
   fitfxn->SetParNames("CONST", "#mu(ns)", "#sigma(ns)");
   ihist->Fit("fitFcn", "LL"); /**Fit with improved LL**/
   //ihist->Fit("gaus", "LL"); /**Fit with improved LL**/
   std::cout << "Printing Fit Parameters for EBEB ......   " << std::endl;
   printf("Integral of function in EBEB = %g\n", fitfxn->Integral( FitLowRange, FitHighRange));

   //*** retrive fit results***//
   int npar = fitfxn->GetNpar();
   TVirtualFitter *fit = TVirtualFitter::GetFitter();
   fit->PrintResults(2,0.);
   TMatrixD *CovMatrix = new TMatrixD ( npar, npar, fit->GetCovarianceMatrix() );
   CovMatrix->Print();

   // Draw Plot with style ///  	
   CanS( c1 );	
   c1->cd();
   ihist->Draw();
   fitfxn->Draw("sames");
   c1->SetLogy(0);
   LegS(leg, fitfxn, ihist );
   leg->Draw();

   cout <<"Saving Canvas..." << endl;
   TString plotname = "EB_EB_Seed_Time_DoubleElectron_Run2012A.pdf" ;

   c1->Print( plotname );
   c1->SaveAs(CNAME);
}       


void FitHist( )
{


   /** Plot Options***/	
   //gROOT->Reset();
   // gROOT->Clear();
   gROOT->SetStyle("Plain") ;
   gROOT->SetBatch(kFALSE);
   gStyle->SetOptTitle(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1);
   gStyle->SetStatX(.89);
   gStyle->SetStatY(.89) ;
   gStyle->SetStatBorderSize(0);
   //gStyle->SetOptStat(1111111)
   gStyle->SetCanvasColor(kWhite);   // background is no longer mouse-dropping white
   gStyle->SetPalette(1);        // blue to red false color palette. Use 9 for b/w
   gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
   gStyle->SetPadBorderMode(0);
   gStyle->SetPaintTextFormat("5.2f");  // What precision to put numbers if plotted with "TEXT"
  // For publishing:
   gStyle->SetLineWidth(2);
   gStyle->SetTextSize(1.1);
   gStyle->SetLabelSize(0.06,"xy");
   gStyle->SetTitleSize(0.08,"xy");
   gStyle->SetTitleOffset(1.2,"x");
   gStyle->SetTitleOffset(1.0,"y");
   gStyle->SetPadTopMargin(0.1);
   gStyle->SetPadRightMargin(0.1);
   gStyle->SetPadBottomMargin(0.16);
   gStyle->SetPadLeftMargin(0.12);
   TGaxis::SetMaxDigits(2); // Set Axis to be of the form 0.11 10^N

   // Begin Fitting
   TFile* ifile  = new TFile("TimePerf-plots.root","READ");
   TF1 *fitf  = new TF1("fitFcn", mygaus, FitLowRange, FitHighRange, 3 );
   TH1F*h_sEB = (TH1F*)ifile->Get(HIST);
   if(h_sEB == 0){ std::cout  <<"!! Histogram Does not exist!!" << std::endl; throw 1;}
   TLegend *lgEB = new TLegend(0.15,0.72,0.3,0.85);
   TCanvas *C = new TCanvas("C",DET,200,10,800,900);

   cout <<" Calling Fitting Fxn" << endl;
   Fit_Hist(h_sEB, fitf, C, lgEB );
}

#ifndef __CINT__
int main() {
   
     FitHist();
    return 0; 
}
#endif
