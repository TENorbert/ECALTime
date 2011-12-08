#include "Histograms.h"
#include "MathFunctions.h"

void hJetTime::Draw( string hfolder, string plotType ) {

      // plot1 , nJets, nPhotons, neElectrons
      c1->SetFillColor(10);
      c1->SetFillColor(10);

      c1->cd();
      nJets->Draw() ;
      c1->Update();
      plotname1 = hfolder + "nJets" + "."+plotType ;
      c1->Print( plotname1 );

      c1->cd();
      nPhotons->Draw() ;
      c1->Update();
      plotname1 = hfolder + "nPhotons" + "."+plotType ;
      c1->Print( plotname1 );

      c1->cd();
      nElectrons->Draw() ;
      c1->Update();
      plotname1 = hfolder + "nElectrons" + "."+plotType ;
      c1->Print( plotname1 );

      // plot2, jet0,jet1, photon0 pt spectrum 
      c2->SetFillColor(10);
      c2->SetFillColor(10);
      c2->SetLogy();

      c2->cd();
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(2);
      jet0Pt->SetLineColor(2) ;
      jet0Pt->Draw() ;
      c2->Update();

      gStyle->SetStatY(0.72);
      gStyle->SetStatTextColor(4);
      jet1Pt->SetLineColor(4) ;
      jet1Pt->DrawCopy("sames") ;
      c2->Update();

      leg2->Clear();
      leg2->AddEntry(jet0Pt, "leading Jet ",  "L");
      leg2->AddEntry(jet1Pt, "2nd jet",       "L");
      leg2->Draw("same");
      c2->Update();

      plotname2 = hfolder + "JetPtSpectrum" + "."+plotType ;
      c2->Print( plotname2 );

      c2->cd();
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(1);
      pho0Pt->Draw() ;
      c2->Update();
      plotname2 = hfolder + "PhotonPtSpectrum" + "."+plotType ;
      c2->Print( plotname2 );

      c2->cd();
      ele0Pt->Draw() ;
      c2->Update();
      plotname2 = hfolder + "ElectronPtSpectrum" + "."+plotType ;
      c2->Print( plotname2 );

      c2->cd();
      hMET->Draw() ;
      c2->Update();
      plotname2 = hfolder + "METSpectrum" + "."+plotType ;
      c2->Print( plotname2 );

      // plot2, jet0, photon0 Eta distribution
      c2->SetFillColor(10);
      c2->SetFillColor(10);
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(1);

      c2->cd();
      pho0Eta->Draw() ;
      c2->Update();
      plotname2 = hfolder + "PhotonEta" + "."+plotType ;
      c2->Print( plotname2 );

      c2->cd();
      jet0Eta->Draw() ;
      c2->Update();
      plotname2 = hfolder + "Jet0Eta" + "."+plotType ;
      c2->Print( plotname2 );

      // plot3  jet and photon time spectrum
      gStyle->SetOptStat("eroum");
      c3->SetFillColor(10);
      c3->SetFillColor(10);
      c3->SetLogy(); 

      c3->cd();
      gStyle->SetStatY(0.9);
      phoXtalTErr->Draw() ;
      c3->Update();
      plotname3 = hfolder + "phoXtalTimeErr" + "."+plotType ;
      c3->Print( plotname3 );

      c3->cd();
      jetXtalTErr->Draw() ;
      c3->Update();
      plotname3 = hfolder + "jetXtalTimeErr" + "."+plotType ;
      c3->Print( plotname3 );

      c3->SetLogy(0); // set linear y-axis
      c3->cd();
      phoXtalTime->Draw() ;
      c3->Update();
      plotname3 = hfolder + "phoXtalTime" + "."+plotType ;
      c3->Print( plotname3 );

      c3->cd();
      jetXtalTime->Draw() ;
      c3->Update();
      plotname3 = hfolder + "jetXtalTime" + "."+plotType ;
      c3->Print( plotname3 );

      // plot 4 , 2D: jet0 Pt vs BC_Time
      gStyle->SetOptStat("eroum");
      c4->SetFillColor(10);
      c4->SetFillColor(10);

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      jet0Pt_T->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet0_Pt_Time" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      jet0Eta_T->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet0_Eta_Time" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      jet_emF_T->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet_emF_Time" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      jet_Pt_nXtal->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet_Pt_nXtal" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      g0_xtalEta_ADC->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "photon_xtalEta_ADC" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      g0_xtalEB_ADC_T->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "photon_xtalEB_ADC_T" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      g0_xtalEB_ADC_E->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "photon_xtalEB_ADC_E" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      j0_xtalEta_ADC->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "j0_xtalEta_ADC" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      j0_xtalEB_ADC_T->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet0_xtalEB_ADC_T" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      j0_xtalEE_ADC_T->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet0_xtalEE_ADC_T" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      j0_xtalEB_ADC_E->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet0_xtalEB_ADC_E" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      j0_xtalEE_ADC_E->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet0_xtalEE_ADC_E" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      xE_cE->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "xtalE_bcE" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      c4->SetLogz(1); // set log z-axis
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(10);
      gStyle->SetStatX(0.9);
      jetXtalPos->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "xtal_eta_phi" + "."+plotType ;
      c4->Print( plotname4 );
      c4->SetLogz(0); // set linear y-axis

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      jPt_dR_gj->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jPt_dR_gj" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      ePt_dR_ge->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "ePt_dR_ge" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      xtalEB_Chi2_T->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "xtalEB_Chi2_T" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      xtalEE_Chi2_T->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "xtalEE_Chi2_T" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      xtal_Chi2_E->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "xtal_Chi2_E" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      xtal_T_E->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "xtal_T_E" + "."+plotType ;
      c4->Print( plotname4 );

}

void hJetTime::FitNDraw( string hfolder, string plotType ) {

      TF1* func ;
      double p1_, p2_, yMax, yMin ;

      gStyle->SetStatX(0.95);

      c5->SetFillColor(10);
      c5->SetFillColor(10);
      c5->SetLogy();

      p1_ = jet0xTB->GetMean(1) ;
      p2_ = jet0xTB->GetRMS(1) ;
      yMin = p1_ - (2.*p2_) ;
      yMax = p1_ + (2.*p2_) ;
      func = new TF1("func",MathFunctions::fitGS, yMin, yMax, 3 );
      func->SetParLimits(1, p1_-p2_ , p1_+p2_ );
      func->SetParLimits(2, p2_*0.5 , p2_*1.5 );

      c5->cd();
      gStyle->SetOptStat("eroum");
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(2);
      jet0lTB->SetLineColor(2) ;
      jet0lTB->Draw() ;
      ///c5->Update();

      gStyle->SetOptFit(111);
      jet0xTB->Fit(func, "RQ0", "", yMin, yMax);

      gStyle->SetStatY(0.6);
      gStyle->SetStatTextColor(4);
      jet0xTB->SetLineColor(4) ;
      //jet0xTB->DrawCopy("sames") ;
      jet0xTB->Draw("sames") ;
      c5->Update();

      func->SetLineColor(4);
      func->DrawCopy("same") ;
      c5->Update();

      leg3->Clear();
      leg3->AddEntry(jet0lTB, "SeedXtal",   "L");
      leg3->AddEntry(jet0xTB, "Weighted Mean", "L");
      leg3->Draw("same");
      c5->Update();

      plotname5 = hfolder + "BarrelJet0Time" + "."+plotType ;
      c5->Print( plotname5 );

      c5->cd();
      gStyle->SetOptStat("eroum");
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(2);
      jet0lTE->SetLineColor(2) ;
      jet0lTE->Draw() ;
      ///c5->Update();

      p1_ = jet0xTE->GetMean(1) ;
      p2_ = jet0xTE->GetRMS(1) ;
      yMin = p1_ - (2.*p2_) ;
      yMax = p1_ + (2.*p2_) ;
      func->SetParLimits(1, p1_-p2_ , p1_+p2_ );
      func->SetParLimits(2, p2_*0.5 , p2_*1.5 );

      gStyle->SetOptFit(111);
      jet0xTE->Fit(func, "RQ0", "", yMin, yMax);

      gStyle->SetStatY(0.6);
      gStyle->SetStatTextColor(4);
      jet0xTE->SetLineColor(4) ;
      jet0xTE->Draw("sames") ;
      c5->Update();

      func->SetLineColor(4);
      func->DrawCopy("same") ;
      c5->Update();

      leg3->Clear();
      leg3->AddEntry(jet0lTE, "SeedXtal",   "L");
      leg3->AddEntry(jet0xTE, "Weighted Mean", "L");
      leg3->Draw("same");
      c5->Update();

      plotname5 = hfolder + "EndcapJet0Time" + "."+plotType ;
      c5->Print( plotname5 );

      c5->cd();
      gStyle->SetOptStat("eroum");
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(2);
      jet1lT->SetLineColor(2) ;
      jet1lT->Draw() ;
      //c5->Update();

      p1_ = jet1xT->GetMean(1) ;
      p2_ = jet1xT->GetRMS(1) ;
      yMin = p1_ - (2.*p2_) ;
      yMax = p1_ + (2.*p2_) ;
      func->SetParLimits(1, p1_-p2_ , p1_+p2_ );
      func->SetParLimits(2, p2_*0.5 , p2_*1.5 );

      gStyle->SetOptFit(111);
      jet1xT->Fit(func, "RQ0", "", yMin, yMax);

      gStyle->SetStatY(0.6);
      gStyle->SetStatTextColor(4);
      jet1xT->SetLineColor(4) ;
      //jet1xT->DrawCopy("sames") ;
      jet1xT->Draw("sames") ;
      c5->Update();

      func->SetLineColor(4);
      func->DrawCopy("same") ;
      c5->Update();

      leg3->Clear();
      leg3->AddEntry(jet1lT, "SeedXtal",   "L");
      leg3->AddEntry(jet1xT, "Weighted Mean", "L");
      leg3->Draw("same");
      c5->Update();

      plotname5 = hfolder + "jet1Time" + "."+plotType ;
      c5->Print( plotname5 );

      // photon/xtal in Barrel
      c5->cd();
      p1_ = pho0lTB->GetMean(1) ;
      p2_ = pho0lTB->GetRMS(1) ;
      yMin = p1_ - (2.*p2_) ;
      yMax = p1_ + (2.*p2_) ;
      func->SetParLimits(1, p1_-p2_ , p1_+p2_ );
      func->SetParLimits(2, p2_*0.2 , p2_*1.5 );

      gStyle->SetOptStat("eou");
      gStyle->SetOptFit(111);
      pho0lTB->Fit(func, "RQ0", "", yMin, yMax);
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(2);
      pho0lTB->SetLineColor(2) ;
      pho0lTB->Draw() ;
      c5->Update();

      gStyle->SetOptStat("eroum");
      gStyle->SetStatY(0.55);
      gStyle->SetStatTextColor(4);
      pho0xTB->SetLineColor(4) ;
      pho0xTB->DrawCopy("sames") ;
      c5->Update();

      func->SetLineColor(2);
      func->DrawCopy("same") ;
      c5->Update();

      leg3->Clear();
      leg3->AddEntry(pho0lTB, "SeedXtal",   "L");
      leg3->AddEntry(pho0xTB, "Weighted Mean", "L");
      leg3->Draw("same");
      c5->Update();

      plotname5 = hfolder + "BarrelPhotonTime" + "."+plotType ;
      c5->Print( plotname5 );

      // photon/xtal in Endcap
      p1_ = pho0lTE->GetMean(1) ;
      p2_ = pho0lTE->GetRMS(1) ;
      yMin = p1_ - (2.*p2_) ;
      yMax = p1_ + (2.*p2_) ;
      func->SetParLimits(1, p1_-p2_ , p1_+p2_ );
      func->SetParLimits(2, p2_*0.5 , p2_*1.5 );

      gStyle->SetOptStat("eou");
      gStyle->SetOptFit(111);
      pho0lTE->Fit(func, "RQ0", "", yMin, yMax);
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(2);
      pho0lTE->SetLineColor(2) ;
      pho0lTE->Draw() ;
      c5->Update();

      gStyle->SetOptStat("eroum");
      gStyle->SetStatY(0.55);
      gStyle->SetStatTextColor(4);
      pho0xTE->SetLineColor(4) ;
      pho0xTE->DrawCopy("sames") ;
      c5->Update();

      func->SetLineColor(2);
      func->DrawCopy("same") ;
      c5->Update();

      leg3->Clear();
      leg3->AddEntry(pho0lTE, "SeedXtal",   "L");
      leg3->AddEntry(pho0xTE, "Weighted Mean", "L");
      leg3->Draw("same");
      c5->Update();

      plotname5 = hfolder + "EndcapPhotonTime" + "."+plotType ;
      c5->Print( plotname5 );

      c5->cd();
      gStyle->SetOptStat("eroum");
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(2);
      elelTB->SetLineColor(2) ;
      elelTB->Draw() ;

      p1_ = elelTB->GetMean(1) ;
      p2_ = elelTB->GetRMS(1) ;
      yMin = p1_ - (2.*p2_) ;
      yMax = p1_ + (2.*p2_) ;
      func->SetParLimits(1, p1_-p2_ , p1_+p2_ );
      func->SetParLimits(2, p2_*0.5 , p2_*1.5 );

      gStyle->SetOptFit(111);
      elelTB->Fit(func, "RQ0", "", yMin, yMax);
      func->SetLineColor(2);
      func->DrawCopy("same") ;
      c5->Update();

      plotname5 = hfolder + "BarrelEleTime" + "."+plotType ;
      c5->Print( plotname5 );

      c5->cd();
      gStyle->SetOptStat("eroum");
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(2);
      elelTE->SetLineColor(2) ;
      elelTE->Draw() ;

      p1_ = pho0lTE->GetMean(1) ;
      p2_ = pho0lTE->GetRMS(1) ;
      yMin = p1_ - (2.*p2_) ;
      yMax = p1_ + (2.*p2_) ;
      func->SetParLimits(1, p1_-p2_ , p1_+p2_ );
      func->SetParLimits(2, p2_*0.5 , p2_*1.5 );

      gStyle->SetOptFit(111);
      elelTE->Fit(func, "RQ0", "", yMin, yMax);
      func->SetLineColor(2);
      func->DrawCopy("same") ;
      c5->Update();

      plotname5 = hfolder + "EndcapEleTime" + "."+plotType ;
      c5->Print( plotname5 );

      // dT( T_photon - T_jet ) 
      c5->SetLogy();
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(1);
      c5->cd();

      p1_ = dT_JetPho->GetMean(1) ;
      p2_ = dT_JetPho->GetRMS(1) ;
      yMin = p1_ - (1.5*p2_) ;
      yMax = p1_ + (1.5*p2_) ;
      //func = new TF1("func",MathFunctions::fitGS, yMin, yMax, 3 );
      //cout<<" ***  m = "<<p1_<<" s = "<< p2_ <<" h2 m:"<< h2->GetMean(2) << endl ;
      //TF1 *func = new TF1("func",MathFunctions::fitGS, yMin, yMax, 3 );
      func->SetParLimits(1, p1_-p2_ , p1_+p2_ );
      func->SetParLimits(2, p2_*0.3 , p2_*1.5 );

      dT_JetPho->Fit(func, "R", "", yMin, yMax);

      dT_JetPho->Draw() ;
      c3->Update();
      func->SetLineColor(2) ;
      func->Draw("same");

      plotname5 = hfolder + "dT_JetPho" + "."+plotType ;
      c5->Print( plotname5 );
 
      delete func ;

}
