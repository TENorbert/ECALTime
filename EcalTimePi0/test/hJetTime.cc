#include "Histograms.h"

void hJetTime::Draw( string hfolder, string plotType ) {

      // plot1 , nJets, nPhotons, neElectrons
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      //c1->Divide(2,2);

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


      // plot3  jet and photon time spectrum
      gStyle->SetOptStat("eroum");
      c3->SetFillColor(10);
      c3->SetFillColor(10);
      c3->SetLogy();

      c3->cd();
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(1);
      jet0lT->SetLineColor(1) ;
      jet0lT->Draw() ;
      c3->Update();

      gStyle->SetStatY(0.7);
      gStyle->SetStatTextColor(2);
      jet0T->SetLineColor(2) ;
      jet0T->DrawCopy("sames") ;
      c3->Update();

      gStyle->SetStatY(0.5);
      gStyle->SetStatTextColor(4);
      jet0xT->SetLineColor(4) ;
      jet0xT->DrawCopy("sames") ;
      c3->Update();

      leg3->Clear();
      leg3->AddEntry(jet0lT, "SeedXtal",   "L");
      leg3->AddEntry(jet0T,  "Average",    "L");
      leg3->AddEntry(jet0xT, "Weighted Mean", "L");
      leg3->Draw("same");
      c3->Update();

      plotname3 = hfolder + "jet0Time" + "."+plotType ;
      c3->Print( plotname3 );

      c3->cd();
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(1);
      jet1lT->SetLineColor(1) ;
      jet1lT->Draw() ;
      c3->Update();

      gStyle->SetStatY(0.7);
      gStyle->SetStatTextColor(2);
      jet1T->SetLineColor(2) ;
      jet1T->DrawCopy("sames") ;
      c3->Update();

      gStyle->SetStatY(0.5);
      gStyle->SetStatTextColor(4);
      jet1xT->SetLineColor(4) ;
      jet1xT->DrawCopy("sames") ;
      c3->Update();

      leg3->Clear();
      leg3->AddEntry(jet1lT, "SeedXtal",   "L");
      leg3->AddEntry(jet1T,  "Average",    "L");
      leg3->AddEntry(jet1xT, "Weighted Mean", "L");
      leg3->Draw("same");
      c3->Update();

      plotname3 = hfolder + "jet1Time" + "."+plotType ;
      c3->Print( plotname3 );

      c3->cd();
      gStyle->SetStatY(0.9);
      gStyle->SetStatTextColor(1);
      pho0lT->SetLineColor(1) ;
      pho0lT->Draw() ;
      c3->Update();

      gStyle->SetStatY(0.7);
      gStyle->SetStatTextColor(2);
      pho0T->SetLineColor(2) ;
      pho0T->DrawCopy("sames") ;
      c3->Update();

      gStyle->SetStatY(0.5);
      gStyle->SetStatTextColor(4);
      pho0xT->SetLineColor(4) ;
      pho0xT->DrawCopy("sames") ;
      c3->Update();

      leg3->Clear();
      leg3->AddEntry(pho0lT, "SeedXtal",   "L");
      leg3->AddEntry(pho0T,  "Average",    "L");
      leg3->AddEntry(pho0xT, "Weighted Mean", "L");
      leg3->Draw("same");
      c3->Update();

      plotname3 = hfolder + "photonTime" + "."+plotType ;
      c3->Print( plotname3 );

      c3->SetLogy(0); // set linear y-axis

      c3->cd();
      gStyle->SetStatY(0.9);
      phoXtalTErr->Draw() ;
      c3->Update();
      plotname3 = hfolder + "phoXtalTimeErr" + "."+plotType ;
      c3->Print( plotname3 );

      c3->cd();
      phoXtalTime->Draw() ;
      c3->Update();
      plotname3 = hfolder + "phoXtalTime" + "."+plotType ;
      c3->Print( plotname3 );

      c3->cd();
      jetXtalTErr->Draw() ;
      c3->Update();
      plotname3 = hfolder + "jetXtalTimeErr" + "."+plotType ;
      c3->Print( plotname3 );

      c3->cd();
      jetXtalTime->Draw() ;
      c3->Update();
      plotname3 = hfolder + "jetXtalTime" + "."+plotType ;
      c3->Print( plotname3 );

      c3->cd();
      dT_JetPho->Draw() ;
      c3->Update();
      plotname3 = hfolder + "dT_JetPho" + "."+plotType ;
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
      jet0Pt_ADC->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet0Pt_ADC" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      jet0Eta_ADC->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet0Eta_ADC" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      jet0T_ADC->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet0T_ADC" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      j0_xtalEBT_ADC->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet0_xtalEBT_ADC" + "."+plotType ;
      c4->Print( plotname4 );

      c4->cd();
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.9);
      j0_xtalEET_ADC->Draw("COLZ") ;
      c4->Update();
      plotname4 = hfolder + "jet0_xtalEET_ADC" + "."+plotType ;
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
}
