#include "PhotonAna.h"

PhotonAna::PhotonAna(){

  Input = new AnaInput();
  
  Input->GetParameters("PlotType",      &plotType ) ; 
  Input->GetParameters("Path",          &hfolder ) ; 
  Input->GetParameters("ProcessEvents", &ProcessEvents ) ; 

  Input->GetParameters( "PhotonCuts",   &photonCuts );
  Input->GetParameters( "ElectronCuts", &electronCuts );
  Input->GetParameters( "JetCuts",      &jetCuts );

  Input->GetParameters( "BasicClusterCuts",  &BasicClusterCuts );
  Input->GetParameters( "XtalCuts",           &XtalCuts );
  //Input->GetParameters("Debug", &debugStr ) ; 
  //if ( debugStr == "True" ) debug = true ;
}

PhotonAna::~PhotonAna(){

  delete Input;
  cout<<" done ! "<<endl ;

}

// analysis template
void PhotonAna::test() { 

  double lumi = 0 ;
  Input->GetParameters( "Lumi", &lumi ) ;
  cout <<" lumi = "<< lumi <<endl;

}
 
void PhotonAna::ReadTree() { 

   //vector<string> rFiles ;
   //Input->GetParameters( "TheData", &rFiles );
   //TTree* tr = Input->GetTree( rFiles[0],"EcalTimeAnalysis" );
   TTree* tr = Input->TreeMap( "data+" );

   float jetPx[MAXOBJ], jetPy[MAXOBJ], jetPz[MAXOBJ], jetE[MAXOBJ] ;
   float metPx[MAXOBJ], metPy[MAXOBJ], metE[MAXOBJ] ;
   float phoPx[MAXOBJ], phoPy[MAXOBJ], phoPz[MAXOBJ], phoE[MAXOBJ] ;
   float elePx[MAXOBJ], elePy[MAXOBJ], elePz[MAXOBJ], eleE[MAXOBJ] ;
   float CPIdx[MAXC], clusterTime[MAXC], clusterEnergy[MAXC]; 
   float xtalInBCEnergy[MAXC][MAXXTALINC], xtalInBCTime[MAXC][MAXXTALINC], xtalInBCTimeErr[MAXC][MAXXTALINC];
   float xtalInBCEta[MAXC][MAXXTALINC],    xtalInBCPhi[MAXC][MAXXTALINC],  xtalADC[MAXC][MAXXTALINC] ;
   float xtalChi2[MAXC][MAXXTALINC];
   int   nXtalInBC[MAXC], clusterXtals[MAXC] ;
   float clusterPhi[MAXC], clusterEta[MAXC] ;
   float vtxX[MAXVTX], vtxY[MAXVTX], vtxZ[MAXVTX], vtxChi2[MAXVTX], vtxNTracks[MAXVTX];
   int   nJets, nPhotons, nElectrons, eventId, nClusters ;

   tr->SetBranchAddress("eventId",    &eventId);
   tr->SetBranchAddress("nJets",      &nJets);
   tr->SetBranchAddress("nPhotons",   &nPhotons);
   tr->SetBranchAddress("nElectrons", &nElectrons);
   tr->SetBranchAddress("nClusters", &nClusters);
   tr->SetBranchAddress("metPx",       metPx );
   tr->SetBranchAddress("metPy",       metPy );
   tr->SetBranchAddress("met",         metE );
   tr->SetBranchAddress("jetPx",       jetPx );
   tr->SetBranchAddress("jetPy",       jetPy );
   tr->SetBranchAddress("jetPz",       jetPz );
   tr->SetBranchAddress("jetE",        jetE );
   tr->SetBranchAddress("phoPx",       phoPx );
   tr->SetBranchAddress("phoPy",       phoPy );
   tr->SetBranchAddress("phoPz",       phoPz );
   tr->SetBranchAddress("phoE",        phoE );
   tr->SetBranchAddress("elePx",       elePx );
   tr->SetBranchAddress("elePy",       elePy );
   tr->SetBranchAddress("elePz",       elePz );
   tr->SetBranchAddress("eleE",        eleE );

   tr->SetBranchAddress("CPIdx",        CPIdx );
   tr->SetBranchAddress("clusterTime",  clusterTime );
   tr->SetBranchAddress("clusterEnergy",  clusterEnergy );
   tr->SetBranchAddress("clusterXtals",   clusterXtals );
   tr->SetBranchAddress("clusterEta",   clusterEta );
   tr->SetBranchAddress("clusterPhi",   clusterPhi );
   tr->SetBranchAddress("clusterXtalsAbove3Sigma",  nXtalInBC );
  
   tr->SetBranchAddress("xtalInBCEnergy",  xtalInBCEnergy );
   tr->SetBranchAddress("xtalInBCTime",    xtalInBCTime );
   tr->SetBranchAddress("xtalInBCTimeErr", xtalInBCTimeErr );
   tr->SetBranchAddress("xtalInBCEta",     xtalInBCEta );
   tr->SetBranchAddress("xtalInBCPhi",     xtalInBCPhi );
   tr->SetBranchAddress("xtalInBCAmplitudeADC", xtalADC );
   tr->SetBranchAddress("xtalInBCChi2",    xtalChi2 );

   tr->SetBranchAddress("vtxX",       vtxX );
   tr->SetBranchAddress("vtxY",       vtxY );
   tr->SetBranchAddress("vtxZ",       vtxZ );
   tr->SetBranchAddress("vtxChi2",    vtxChi2 );
   tr->SetBranchAddress("vtxNTracks", vtxNTracks );

   int totalN = tr->GetEntries();
   cout<<" total entries = "<< totalN <<" Process "<< ProcessEvents <<endl;
    
   for ( int i=0; i< tr->GetEntries(); i++ ) {
       if ( ProcessEvents > 0 && i > ( ProcessEvents - 1 ) ) break;
       tr->GetEntry( i );
     
       //cout<<" Event "<< eventId <<endl ;
       hJets.nJets->Fill( nJets );
       hJets.nPhotons->Fill( nPhotons );
       hJets.nElectrons->Fill( nElectrons );

       // 1. met selection
       TLorentzVector metp4( metPx[0], metPy[0], 0., metE[0] ) ;
       if ( metp4.Et()  <  jetCuts[3] ) continue ;

       // 2. jet selection
       int nu_Jets = 0 ;
       jetV.clear() ;
       for ( int j=0 ; j< nJets; j++ ) {
           TLorentzVector jp4( jetPx[j], jetPy[j], jetPz[j], jetE[j] ) ;
           if ( jp4.Pt() < jetCuts[0] || fabs(jp4.Eta()) > jetCuts[1] ) continue ;
           nu_Jets++ ;
           jetV.push_back( jp4 );
       }
       if ( nu_Jets < jetCuts[2] ) continue ;

       // 2. get the leading electron info, just check
       eleV.clear() ;
       for ( int i=0 ; i< nElectrons; i++ ) {
           TLorentzVector eP4( elePx[i], elePy[i], elePz[i], eleE[i] ) ;
           if ( eP4.Pt() < electronCuts[0] || fabs( eP4.Eta()) > electronCuts[1] ) continue ;
           // check the isolation -- using dR
           double dR = 999. ;
           for ( int j=0 ; j< nJets; j++ ) {
               TLorentzVector jp4( jetPx[j], jetPy[j], jetPz[j], jetE[j] ) ;
               if ( jp4.Pt() < jetCuts[0] || fabs(jp4.Eta()) > jetCuts[1] ) continue ;
                dR  = (  eP4.DeltaR( jp4 ) < dR ) ?  eP4.DeltaR( jp4 ) : dR;
           }
           if ( dR < electronCuts[3] ) continue ;
           hJets.ele0Pt->Fill( eP4.Pt() );
           eleV.push_back( eP4.Pt() ) ;
       }

       // 3. photon cuts
       int nu_Photons = 0 ;
       int gid = -1 ;
       for ( int i=0 ; i< nPhotons; i++ ) {
           TLorentzVector phoP4( phoPx[i], phoPy[i], phoPz[i], phoE[i] ) ;
           if ( phoP4.Pt() < photonCuts[0] || fabs(phoP4.Eta()) > photonCuts[1] ) continue ;
           // check the isolation -- using dR_gj
           double dR_gj = 999. ;
           double jpt = 0 ;
           for ( int j=0 ; j< nJets; j++ ) {
               TLorentzVector jp4( jetPx[j], jetPy[j], jetPz[j], jetE[j] ) ;
               if ( jp4.Pt() < jetCuts[0] || fabs(jp4.Eta()) > jetCuts[1] ) continue ;
               if ( phoP4.DeltaR( jp4 ) < dR_gj ) {
                  dR_gj  = phoP4.DeltaR( jp4 ) ;
                  jpt = jp4.Pt() ;
               }
           }
           if ( dR_gj < photonCuts[3] ) continue ;

           // dR_ge 
           double dR_ge = 999. ;
           double ept = 0 ;
           for (size_t k =0 ; k < eleV.size() ; k++) {
               if ( phoP4.DeltaR( eleV[k] ) < dR_ge )  {
                  dR_ge = phoP4.DeltaR( eleV[k] ) ;
                  ept = eleV[k].Pt() ;
               }
           }
           if ( dR_ge < photonCuts[3] ) continue ;
           hJets.jPt_dR_gj->Fill( jpt, dR_gj ) ;
           hJets.ePt_dR_ge->Fill( ept, dR_ge ) ;

           gid = ( phoP4.Pt() > photonCuts[0] ) ? i : gid ;
           nu_Photons++ ;
       }
       if ( nu_Photons < photonCuts[4] || gid < 0 ) continue ;
       //cout<<" EventID : "<< eventId ;
       //cout<<" Jet size = "<< jetV.size() <<" Photon ID = "<< gid <<endl;

       TLorentzVector gp4( phoPx[gid], phoPy[gid], phoPz[gid], phoE[gid] ) ;
       hJets.pho0Pt->Fill( gp4.Pt() );
       //** Photon Time info
       //cout<<" photon 0  pt : "<< gp4.Pt() <<endl;
       float gTime  = 0. ;
       float gxTime = 0. ;
       float gleadEnergy = 0 ;
       float gleadTime   = 0 ;
       float ngX = 0. ;
       float ngC = 0. ;
       for ( int x =0; x< nClusters; x++ ) {
           if ( CPIdx[x] < 22. || CPIdx[x] > 23. ) continue ;
	   float gidx = 22.1 + gid ;
           if ( clusterEnergy[x] < BasicClusterCuts[0] )
	   if ( CPIdx[x] != gidx || fabs(clusterTime[x]) > BasicClusterCuts[1]  ) continue ;
           
	   //cout<<"  cluster "<< CPIdx[x] <<" time : "<< clusterTime[x] ;
	   //cout<<" E: "<< clusterEnergy[x]<<" N: "<< clusterXtals[x]  << endl;

	   float seedEnergy = 0 ;
	   float seedTime   = 0 ;
	   float seedTErr   = 999999 ;
           for (int y=0; y< nXtalInBC[x]; y++) {
               if ( xtalChi2[x][y] > XtalCuts[2] ) continue ;
	       if ( xtalInBCEnergy[x][y] < XtalCuts[0] ) continue ;
	       if ( xtalInBCTimeErr[x][y] < 0.2 || xtalInBCTimeErr[x][y] > XtalCuts[1] ) continue ;
	       hJets.phoXtalTime->Fill( xtalInBCTime[x][y] ); 
	       hJets.phoXtalTErr->Fill( xtalInBCTimeErr[x][y] ); 
	       if ( xtalInBCEnergy[x][y] > gleadEnergy ) {                         
                  gleadEnergy = xtalInBCEnergy[x][y] ;
                  gleadTime   = xtalInBCTime[x][y] ;
               }
               if ( xtalInBCEnergy[x][y] > seedEnergy ) {
                  seedEnergy = clusterEnergy[x] ; 
                  seedTime   = xtalInBCTime[x][y] ;
                  seedTErr   = xtalInBCTimeErr[x][y] ;
               }
	       gxTime += xtalInBCTime[x][y] / pow( xtalInBCTimeErr[x][y], 2) ;
               ngX += 1./ pow( xtalInBCTimeErr[x][y],2) ;
           }
           gTime +=  seedTime / pow( seedTErr, 2) ;
           ngC   +=       1. / pow( seedTErr, 2) ;
       }
       if ( ngX > 0. && ngC > 0. ) {
          gTime  = gTime  / ngC ;
          gxTime = gxTime / ngX ;
          hJets.pho0T->Fill( gTime );
	  hJets.pho0xT->Fill( gxTime ); 
	  hJets.pho0lT->Fill( gleadTime ); 
       }

        
       //** 2. get the associated bc for jet
       for ( int j=0 ; j< (int)(jetV.size()); j++) {
           if ( j > 1 ) break ; 
           TLorentzVector jp4( jetV[j].Px(), jetV[j].Py(), jetV[j].Pz(), jetV[j].E() ) ;
           //cout<<" jet ("<< jetV[j].Px() <<","<<jetV[j].Py()<<","<<jetV[j].Pz()<<","<< jetV[j].E()<<")" <<endl;
	   if ( j == 0 ) hJets.jet0Pt->Fill( jp4.Pt() );
	   if ( j == 1 ) hJets.jet1Pt->Fill( jp4.Pt() );

	   float emF = 0 ;
	   float nC  = 0 ;
	   float nX  = 0 ;
	   float xE  = 0 ;
	   float cE  = 0 ;
	   float jleadEnergy = 0 ;
	   float jleadTime   = 0 ;
	   float jTime  = 0 ;
	   float jxTime = 0 ;
	   float xADC   = 0 ;  
	   int   nXtal  = 0 ; 

           // loop the clusters
	   for (int x =0; x< nClusters; x++) {

	       if ( CPIdx[x] < 100. || CPIdx[x] > 101. ) continue ;

	       float jidx =  100. + (j+1.)*0.1 ;
	       if ( clusterEnergy[x] < BasicClusterCuts[0] ) continue;
	       if ( CPIdx[x] != jidx  || fabs(clusterTime[x]) > BasicClusterCuts[1] ) continue ;

	       double dRcj = DeltaR( clusterEta[x], clusterPhi[x], jetV[j].Eta(), jetV[j].Phi() ) ;
	       if ( dRcj > BasicClusterCuts[2] ) continue ;
	       //cout<<" BC"<<x<<" CPIdx = "<< CPIdx[x] << " dRcj = "<< dRcj <<endl ;
	       float seedEnergy = -999999 ;
	       float seedTime   = -999999 ;
	       float seedTErr   =  999999 ;
	       //cout<<"  cluster "<<CPIdx[x] <<" time : "<< clusterTime[x] ;
	       //cout<<" E: "<< clusterEnergy[x]<<" N: "<< clusterXtals[x]  << endl;
	       float xADC_BC = 0. ;
	       // loop the xtals
	       for (int y=0; y< nXtalInBC[x]; y++) {
                   if ( xtalInBCEnergy[x][y] < XtalCuts[0] ) continue ;
		   if ( xtalInBCTimeErr[x][y] < 0.2 || xtalInBCTimeErr[x][y] > XtalCuts[1] ) continue ;
		   if ( xtalChi2[x][y] > XtalCuts[2] ) continue ;
		   double dRxj = DeltaR( xtalInBCEta[x][y], xtalInBCPhi[x][y] , jetV[j].Eta(), jetV[j].Phi()  ) ;
		   if ( dRxj > XtalCuts[3] ) continue ;
		   if ( j == 0 ) hJets.jetXtalTime->Fill( xtalInBCTime[x][y] ); 
		   if ( j == 0 ) hJets.jetXtalTErr->Fill( xtalInBCTimeErr[x][y] ); 
		   if ( j == 0 ) hJets.jetXtalPos->Fill( xtalInBCEta[x][y], xtalInBCPhi[x][y] );
		   if ( xtalInBCEnergy[x][y] > jleadEnergy ) {                         
		      jleadEnergy = xtalInBCEnergy[x][y] ;
		      jleadTime   = xtalInBCTime[x][y] ;
                   }
		   if ( xtalInBCEnergy[x][y] > seedEnergy ) {                         
                      seedEnergy = xtalInBCEnergy[x][y] ;
		      seedTime   = xtalInBCTime[x][y] ;
		      seedTErr   = xtalInBCTimeErr[x][y] ;
                   }
		   jxTime += xtalInBCTime[x][y] / pow( xtalInBCTimeErr[x][y], 2) ;
		   nX     += 1. / pow( xtalInBCTimeErr[x][y],2) ;
		   //double sinTheta   = fabs( sin( 2. *atan( exp(-1*xtalInBCEta[x][y]  ) ) ) );
		   emF += xtalInBCEnergy[x][y] ;
		   xE  += xtalInBCEnergy[x][y] ;
		   xADC_BC = ( xtalADC[x][y] > xADC_BC ) ? xtalADC[x][y] : xADC_BC ;
		   nXtal++ ;
		   if ( j == 0 && fabs(xtalInBCEta[x][y])  < 1.479 ) hJets.j0_xtalEBT_ADC->Fill( xtalADC[x][y], xtalInBCTime[x][y] );
                   if ( j == 0 && fabs(xtalInBCEta[x][y]) >= 1.479 ) hJets.j0_xtalEET_ADC->Fill( xtalADC[x][y], xtalInBCTime[x][y] );
               }
	       jTime += seedTime/ pow( seedTErr, 2 ) ;
	       nC    += 1./ pow( seedTErr, 2) ;
	       xADC  += xADC_BC  ;
	       if ( nXtal > 0 ) cE += clusterEnergy[x] ;
           }
	   if ( nC > 0. && nXtal > 0 ) {
	      jTime  = jTime / nC ;
	      jxTime = jxTime / nX ;
	      emF    = emF / jp4.E() ;
	      hJets.jet_emF_T->Fill( emF , jTime ) ;
	      hJets.xE_cE->Fill( xE , cE ) ;
	      hJets.jet_Pt_nXtal->Fill( jp4.Pt() , nXtal ) ;
	      if ( j == 0 ) { 
                 hJets.jet0T->Fill( jTime ); 
		 hJets.jet0xT->Fill( jxTime ); 
		 hJets.jet0lT->Fill( jleadTime ); 
		 hJets.jet0Pt_T->Fill( jp4.Pt() , jTime ) ;
		 hJets.jet0Eta_T->Fill( jp4.Eta() , jTime ) ;
		 hJets.jet0Pt_ADC->Fill( jp4.Pt() , xADC ) ;
		 hJets.jet0Eta_ADC->Fill( jp4.Eta() , xADC ) ;
		 hJets.jet0T_ADC->Fill( xADC, jxTime ) ;
                 if ( ngX > 0. ) hJets.dT_JetPho->Fill( gleadTime - jxTime ) ;  
		 //cout<<"  jxTime: "<< jxTime <<" nX: "<< nX <<" nC: "<<nC<<" jTime:"<< jTime << endl ;
	      }
              if ( j == 1 ) {
                 hJets.jet1T->Fill( jTime );
		 hJets.jet1xT->Fill( jxTime );
		 hJets.jet1lT->Fill( jleadTime ); 
              }
	   }
       }
       

   }

   hJets.Draw( hfolder, plotType ) ;
}

void PhotonAna::ScalarPlotList() {

   bool debugPlots = false ;
   ScalarPlotter( hJets.j0_xtalEBT_ADC, "j0_xtalEBT_ADC", -2., 2., 5 ) ;
   ScalarPlotter( hJets.j0_xtalEET_ADC, "j0_xtalEET_ADC", -2., 2., 5 ) ;
   ScalarPlotter( hJets.jet0T_ADC,    "jet_T_ADC",    -2., 2., 2 ) ;
   ScalarPlotter( hJets.jet0Eta_ADC,  "jet_Eta_ADC",   300., 6300., 2 ) ;
   ScalarPlotter( hJets.jet0Pt_ADC,   "jet_Pt_ADC",    300., 6300., 2, debugPlots ) ;
   ScalarPlotter( hJets.jet_Pt_nXtal, "jet_Pt_nXtal",  0., 50., 2 ) ;
   ScalarPlotter( hJets.jet0Eta_T,    "jet_Eta_Time",  -2,   2 ) ;
   ScalarPlotter( hJets.jet0Pt_T,     "jet_Pt_Time",   -2,   2 ) ;

}

void PhotonAna::ScalarPlotter( TH2D* h2, TString hname, double yMin, double yMax, int rbin, bool debugPlots ) {

    TCanvas* c1 = new TCanvas("c1","", 800, 600);
    c1->SetFillColor(10);
    c1->SetFillColor(10);
    gPad->SetGridx();
    gPad->SetGridy();

    vector<double> xV ;
    vector<double> yV ;
    vector<double> yErrV ;
    BinningFitScan( h2, xV, yV, yErrV, debugPlots, rbin ) ; 

    const int sz = static_cast<const int>( xV.size() );
    float  x[sz] , xErr[sz], y[sz], yErr[sz] ;
    float yErrMax = 0 ;
    for ( int i=0; i< sz; i++) {
        x[i] = xV[i] ;
        y[i] = yV[i] ;
        yErr[i] = yErrV[i] ;
        xErr[i] = 0.  ;
        yErrMax = ( yErr[i] > yErrMax ) ? yErr[i] : yErrMax ;
    }

    c1->cd();
    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptFit(111);

    TGraphErrors* g1 = new TGraphErrors( sz, x, y,  xErr, yErr );
    g1->SetMaximum( yMax );
    g1->SetMinimum( yMin );
    g1->SetMarkerColor(4);
    g1->SetMarkerStyle(20);
    g1->SetLineWidth(2);
    g1->SetLineColor(4);
    g1->SetTitle( hname );
    g1->Draw("AP");
 
    c1->Update();
    TString plotname = hfolder + hname +  "_pj."+plotType ;
    c1->Print( plotname );

    TGraph* g2 = new TGraph( sz, x, yErr );
    g2->SetMaximum( yErrMax*2 );
    g2->SetMinimum( 0. );
    g2->SetMarkerColor(4);
    g2->SetMarkerStyle(20);
    g2->SetLineWidth(2);
    g2->SetLineColor(4);
    TString plotName = hname + "Error Distribution" ;
    g2->SetTitle( hname );
    g2->Draw("ACP");
 
    c1->Update();
    plotname = hfolder + hname +  "_pjerr."+plotType ;
    c1->Print( plotname );

    delete g1 ;
    delete g2 ;
    delete c1 ;
}


void PhotonAna::BinningFitScan( TH2D* h2, vector<double>& xV, vector<double>& yV, vector<double>& yErrV, bool debugPlots, int rbin, int startBin, int finalBin ) {

    int      sz = h2->GetNbinsX();
    for ( int i=0; i< sz; i++) {
        int b1 = startBin + (i*rbin) ;
        int b2 = b1 + rbin - 1 ;
        if ( b1 == finalBin ) break;
        vector<double> results = BinningFit( h2, "h2Py", b1, b2, debugPlots ) ;

        if ( results[2] == 0. ) continue ;
        xV.push_back( h2->GetBinCenter(b1) ) ;
        yV.push_back( results[1] ) ;
        yErrV.push_back( results[2] ) ;
    }

}


vector<double> PhotonAna::BinningFit( TH2D* h2, string hName, int xbinMin, int xbinMax, bool debugPlots ) { 

   TH1D* h1_pjy = h2->ProjectionY( hName.c_str() , xbinMin, xbinMax ) ;

   double p1_ = h1_pjy->GetMean(1) ;
   double p2_ = h1_pjy->GetRMS(1) ;
   double yMin = p1_ - (3*p2_) ;
   double yMax = p1_ + (3*p2_) ;
   //cout<<" ***  m = "<<p1_<<" s = "<< p2_ <<" h2 m:"<< h2->GetMean(2) << endl ;
   TF1 *func = new TF1("func",MathFunctions::fitGS, yMin, yMax, 3 );
   func->SetParLimits(1, p1_-p2_ , p1_+p2_ );
   func->SetParLimits(2, p2_*0.8 , p2_*1.2 );

   h1_pjy->Fit(func, "R", "", yMin, yMax);
 
   double amp   = func->GetParameter(0) ;
   double mean  = func->GetParameter(1) ;
   double sigma = func->GetParameter(2) ;

   vector<double> result ;
   result.push_back( amp );
   result.push_back( mean );
   result.push_back( sigma );

   TString plotname0 ;
   if ( debugPlots ) {
      cout <<" debugging slice&fit"<< endl; 
      gSystem->cd( hfolder.c_str() ) ;
      gSystem->mkdir( "debug" );
      gSystem->cd( "../" );
     
      TCanvas* c0 = new TCanvas("c0","", 800, 600);
      c0->SetFillColor(10);
      c0->SetFillColor(10);
      c0->cd();
      h1_pjy->Draw() ;
      c0->Update();
      ostringstream sfxStr ;
      sfxStr << "_py" ;
      sfxStr << xbinMin ;

      plotname0= hfolder + "debug/" + hName +  sfxStr.str() +"."+plotType ;
      c0->Print( plotname0 );
      delete c0 ;
   }

   delete h1_pjy ;
   delete func ;
   return result ;

}

double PhotonAna::DeltaR( float eta1, float phi1, float eta2, float phi2 ) {

      double df = fabs( phi1 - phi2 );
      if ( df > 3.1416 ) df =  6.2832 - df  ;
      double dh = fabs( eta1 - eta2 );

      double dR = sqrt( (df*df) + (dh*dh) ) ;

      return dR ;
}

