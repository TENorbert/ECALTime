#include "PhotonAna.h"

PhotonAna::PhotonAna( string datacardfile ){

  Input = new AnaInput( datacardfile );
  
  Input->GetParameters("PlotType",      &plotType ) ; 
  Input->GetParameters("Path",          &hfolder ) ; 
  Input->GetParameters("ProcessEvents", &ProcessEvents ) ; 

  Input->GetParameters( "PhotonCuts",   &photonCuts );
  Input->GetParameters( "ElectronCuts", &electronCuts );
  Input->GetParameters( "JetCuts",      &jetCuts );

  Input->GetParameters( "SelectBackground",  &selectBackground );
  Input->GetParameters( "SplitEvent",        &split );

  Input->GetParameters( "BasicClusterCuts",  &BasicClusterCuts );
  Input->GetParameters( "XtalCuts",          &XtalCuts );
  Input->GetParameters( "DoTimeCorrection",  &doTimeCorrection );
  //Input->GetParameters("Debug", &debugStr ) ; 
  //if ( debugStr == "True" ) debug = true ;
  theTimeCorrector_.initEB("EB");
  theTimeCorrector_.initEE("EE");

}

PhotonAna::~PhotonAna(){

  delete Input;
  cout<<" done ! "<<endl ;

}

// analysis template
void PhotonAna::ReadTree() { 

   //vector<string> rFiles ;
   //Input->GetParameters( "TheData", &rFiles );
   //TTree* tr = Input->GetTree( rFiles[0],"EcalTimeAnalysis" );
   TTree* tr = Input->TreeMap( "data+" );

   float jetPx[MAXOBJ], jetPy[MAXOBJ], jetPz[MAXOBJ], jetE[MAXOBJ] ;
   float phoPx[MAXOBJ], phoPy[MAXOBJ], phoPz[MAXOBJ], phoE[MAXOBJ] ;
   float elePx[MAXOBJ], elePy[MAXOBJ], elePz[MAXOBJ], eleE[MAXOBJ] ;
   float CPIdx[MAXC], clusterTime[MAXC], clusterEnergy[MAXC]; 
   float xtalInBCEnergy[MAXC][MAXXTALINC], xtalInBCTime[MAXC][MAXXTALINC], xtalInBCTimeErr[MAXC][MAXXTALINC];
   float xtalInBCEta[MAXC][MAXXTALINC],    xtalInBCPhi[MAXC][MAXXTALINC],  xtalADC[MAXC][MAXXTALINC] ;
   float xtalChi2[MAXC][MAXXTALINC];
   float clusterPhi[MAXC], clusterEta[MAXC] ;
   float vtxX[MAXVTX], vtxY[MAXVTX], vtxZ[MAXVTX], vtxChi2[MAXVTX], vtxNTracks[MAXVTX];
   int   nXtalInBC[MAXC], clusterXtals[MAXC] ;
   float metPx, metPy, metE ;
   int   nJets, nPhotons, nElectrons, eventId, nClusters ;

   tr->SetBranchAddress("eventId",    &eventId);
   tr->SetBranchAddress("nJets",      &nJets);
   tr->SetBranchAddress("nPhotons",   &nPhotons);
   tr->SetBranchAddress("nElectrons", &nElectrons);
   tr->SetBranchAddress("nClusters",  &nClusters);
   tr->SetBranchAddress("metPx",      &metPx );
   tr->SetBranchAddress("metPy",      &metPy );
   tr->SetBranchAddress("met",        &metE );
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

   tr->SetBranchAddress("CPIdx",         CPIdx );
   tr->SetBranchAddress("clusterTime",   clusterTime );
   tr->SetBranchAddress("clusterEnergy", clusterEnergy );
   tr->SetBranchAddress("clusterXtals",  clusterXtals );
   tr->SetBranchAddress("clusterEta",    clusterEta );
   tr->SetBranchAddress("clusterPhi",    clusterPhi );
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

       if ( i%2 == split ) continue ;    
 
       //cout<<" Event "<< eventId <<endl ;
       // 1. met selection
       TLorentzVector metp4( metPx, metPy, 0., metE ) ;
       if ( jetCuts[4] > 0 &&  metp4.Et() < jetCuts[4] ) continue ;
       if ( jetCuts[4] < 0 &&  metp4.Et() > fabs( jetCuts[4] ) ) continue ;

       // 2. jet selection
       int nu_Jets = 0 ;
       jetV.clear() ;
       for ( int j=0 ; j< nJets; j++ ) {
           TLorentzVector jp4( jetPx[j], jetPy[j], jetPz[j], jetE[j] ) ;
           if ( jp4.Pt() < jetCuts[0] || fabs(jp4.Eta()) > jetCuts[1] ) continue ;
           nu_Jets++ ;
           jetV.push_back( jp4 );
       }
       if ( nu_Jets < jetCuts[2] || nu_Jets > jetCuts[3] ) continue ;

       // 2. get the leading electron info, just check
       eleV.clear() ;
       eIdV.clear() ;
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
           eIdV.push_back( 11.1 + (0.1*i) );
       }

       // 3. photon cuts
       int nu_Photons = 0 ;
       int gid = -1 ;
       double maxPhotonPt = photonCuts[0] ;
       for ( int i=0 ; i< nPhotons; i++ ) {
           bool SignalLike = false ;
           TLorentzVector phoP4( phoPx[i], phoPy[i], phoPz[i], phoE[i] ) ;
           if ( phoP4.Pt() < photonCuts[0] || fabs(phoP4.Eta()) > photonCuts[1] ) continue ;
           if ( phoP4.Pt() > 100. ) SignalLike = true ;
           if ( SignalLike && selectBackground == 1 ) continue ;
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

           gid = ( phoP4.Pt() > maxPhotonPt ) ? i : gid ;
           nu_Photons++ ;
       }
       if ( nu_Photons < photonCuts[4] || gid < 0 ) continue ;
       //cout<<" EventID : "<< eventId ;
       //cout<<" Jet size = "<< jetV.size() <<" Photon ID = "<< gid <<endl;

       // event profile
       TLorentzVector gp4( phoPx[gid], phoPy[gid], phoPz[gid], phoE[gid] ) ;
       hJets.pho0Pt->Fill( gp4.Pt() );
       hJets.pho0Eta->Fill( gp4.Eta() );
       hJets.hMET->Fill( metE );
       hJets.nJets->Fill( nu_Jets );
       hJets.nPhotons->Fill( nu_Photons );
       hJets.nElectrons->Fill( eleV.size() );

       //** Photon Time info
       float gleadEnergy = 0 ;
       float gleadTime   = 0 ;
       float gleadEta    = 999 ;
       float gTime  = 0. ;
       float ngC = 0. ;
       for ( int x =0; x< nClusters; x++ ) {
           if ( CPIdx[x] < 22. || CPIdx[x] > 23. ) continue ;
	   float gidx = 22.1 + gid ;
           if ( clusterEnergy[x] < BasicClusterCuts[0] )
	   if ( CPIdx[x] != gidx || fabs(clusterTime[x]) > BasicClusterCuts[1]  ) continue ;
           
	   //cout<<"  cluster "<< CPIdx[x] <<" time : "<< clusterTime[x] ;
	   //cout<<" E: "<< clusterEnergy[x]<<" N: "<< clusterXtals[x]  << endl;

	   float seedEnergy = -999999 ;
	   float seedTime   = -999999 ;
	   float seedTErr   =  999999 ;
           for (int y=0; y< nXtalInBC[x]; y++) {

               double xtalTime  =  xtalInBCTime[x][y] 
                         + (xtalInBCTime[x][y]*theTimeCorrector_.getCorrection( (float) xtalInBCEnergy[x][y], xtalInBCEta[x][y]) );
               if ( doTimeCorrection == 0 )  xtalTime = xtalInBCTime[x][y] ;

	       hJets.phoXtalTime->Fill( xtalTime); 
	       hJets.phoXtalTErr->Fill( xtalInBCTimeErr[x][y] ); 
               if ( fabs(xtalInBCEta[x][y])  < 1.479 ) hJets.xtalEB_Chi2_T->Fill( xtalChi2[x][y], xtalTime );
               if ( fabs(xtalInBCEta[x][y]) >= 1.479 ) hJets.xtalEE_Chi2_T->Fill( xtalChi2[x][y], xtalTime );
               hJets.xtal_Chi2_E->Fill( xtalChi2[x][y], xtalInBCEnergy[x][y]);
               hJets.xtal_T_E->Fill( xtalTime, xtalInBCEnergy[x][y]);

	       if ( xtalInBCEnergy[x][y] < XtalCuts[0] ) continue ;
	       if ( xtalInBCTimeErr[x][y] < 0.2 || xtalInBCTimeErr[x][y] > XtalCuts[1] ) continue ;
               if ( fabs(xtalInBCEta[x][y])  < 1.479 && xtalChi2[x][y] > XtalCuts[2] ) continue ;
               if ( fabs(xtalInBCEta[x][y]) >= 1.479 && xtalChi2[x][y] > XtalCuts[3] ) continue ;

	       if ( xtalInBCEnergy[x][y] > gleadEnergy ) {                         
                  gleadEnergy = xtalInBCEnergy[x][y] ;
                  gleadTime   = xtalTime ;
                  gleadEta    = fabs( xtalInBCEta[x][y] );
               }
               if ( xtalInBCEnergy[x][y] > seedEnergy ) {
                  seedEnergy = xtalInBCEnergy[x][y] ; 
                  seedTime   = xtalTime ;
                  seedTErr   = xtalInBCTimeErr[x][y] ;
               }
               hJets.g0_xtalEta_ADC->Fill( xtalInBCEta[x][y], xtalADC[x][y] );
	       if ( fabs(xtalInBCEta[x][y])  < 1.479 ) hJets.g0_xtalEB_ADC_T->Fill( xtalADC[x][y], xtalTime );
	       if ( fabs(xtalInBCEta[x][y])  < 1.479 ) hJets.g0_xtalEB_ADC_E->Fill( xtalADC[x][y], xtalInBCEnergy[x][y] );
           }
           if ( seedEnergy != -999999 ) {
              gTime +=  seedTime / pow( seedTErr, 2) ;
              ngC   +=        1. / pow( seedTErr, 2) ;
           }
       } // end of photon loop
       if ( ngC > 0. ) {
          gTime  = gTime  / ngC ;
	  if ( gp4.Eta()  < 1.479 ) hJets.pho0xTB->Fill( gTime ); 
	  if ( gleadEta   < 1.479 ) hJets.pho0lTB->Fill( gleadTime ); 
	  if ( gp4.Eta() >= 1.479 ) hJets.pho0xTE->Fill( gTime ); 
	  if ( gleadEta  >= 1.479 ) hJets.pho0lTE->Fill( gleadTime ); 
       }

        
       //** 2. get the associated bc for jet
       int k = 0 ;
       float jTime  = 0 ;
       for ( int j=0 ; j< (int)(jetV.size()); j++) {
           if ( k > 2 ) break ; 
           TLorentzVector jp4( jetV[j].Px(), jetV[j].Py(), jetV[j].Pz(), jetV[j].E() ) ;
           //cout<<" jet ("<< jetV[j].Px() <<","<<jetV[j].Py()<<","<<jetV[j].Pz()<<","<< jetV[j].E()<<")" <<endl;
	   if ( j == 0 ) hJets.jet0Pt->Fill( jp4.Pt() );
	   if ( j == 0 ) hJets.jet0Eta->Fill( jp4.Eta() );
	   if ( j == 1 ) hJets.jet1Pt->Fill( jp4.Pt() );

           // loop the clusters
	   float xE  = 0 ;
	   float cE  = 0 ;
	   float jleadEnergy = 0 ;
	   float jleadTime   = 0 ;
           float jleadEta    = 999 ;
	   float jxTime  = 0 ;
	   int   nXtal   = 0 ; 
	   float nC  = 0 ;
	   for (int x =0; x< nClusters; x++) {

	       if ( CPIdx[x] < 100. || CPIdx[x] > 101. ) continue ;

	       float jidx =  100. + (j+1.)*0.1 ;
	       if ( clusterEnergy[x] < BasicClusterCuts[0] ) continue;
	       if ( CPIdx[x] != jidx  || fabs(clusterTime[x]) > BasicClusterCuts[1] ) continue ;

	       double dRcj = DeltaR( clusterEta[x], clusterPhi[x], jetV[j].Eta(), jetV[j].Phi() ) ;
	       if ( dRcj > BasicClusterCuts[2] ) continue ;

	       //cout<<" BC"<<x<<" CPIdx = "<< CPIdx[x] << " dRcj = "<< dRcj <<endl ;
	       float seedEnergy = -999999 ;
	       float seedTime   =  0 ;
	       float seedTErr   =  999999 ;
	       for (int y=0; y< nXtalInBC[x]; y++) {

                   double xtalTime  =  xtalInBCTime[x][y] 
                       + (xtalInBCTime[x][y]*theTimeCorrector_.getCorrection( (float) xtalInBCEnergy[x][y], xtalInBCEta[x][y]) );
                   if ( doTimeCorrection == 0 )  xtalTime = xtalInBCTime[x][y] ;

                   if ( fabs(xtalInBCEta[x][y])  < 1.479 ) hJets.xtalEB_Chi2_T->Fill( xtalChi2[x][y], xtalTime );
                   if ( fabs(xtalInBCEta[x][y]) >= 1.479 ) hJets.xtalEE_Chi2_T->Fill( xtalChi2[x][y], xtalTime );
                   hJets.xtal_Chi2_E->Fill( xtalChi2[x][y], xtalInBCEnergy[x][y]);
                   hJets.xtal_T_E->Fill( xtalTime, xtalInBCEnergy[x][y]);

		   if ( j == 0 ) hJets.jetXtalTime->Fill( xtalTime ); 
		   if ( j == 0 ) hJets.jetXtalTErr->Fill( xtalInBCTimeErr[x][y] ); 
		   if ( j == 0 ) hJets.jetXtalPos->Fill( xtalInBCEta[x][y], xtalInBCPhi[x][y] );

                   if ( xtalInBCEnergy[x][y] < XtalCuts[0] ) continue ;
		   if ( xtalInBCTimeErr[x][y] < 0.2 || xtalInBCTimeErr[x][y] > XtalCuts[1] ) continue ;
                   if ( fabs(xtalInBCEta[x][y])  < 1.479 && xtalChi2[x][y] > XtalCuts[2] ) continue ;
                   if ( fabs(xtalInBCEta[x][y]) >= 1.479 && xtalChi2[x][y] > XtalCuts[3] ) continue ;
		   double dRxj = DeltaR( xtalInBCEta[x][y], xtalInBCPhi[x][y] , jetV[j].Eta(), jetV[j].Phi()  ) ;
		   if ( dRxj > XtalCuts[4] ) continue ;

		   xE  += xtalInBCEnergy[x][y] ;

		   if ( xtalInBCEnergy[x][y] > jleadEnergy ) {                         
		      jleadEnergy = xtalInBCEnergy[x][y] ;
		      jleadTime   = xtalTime ;
                      jleadEta    = xtalInBCEta[x][y] ;
                   }
		   if ( xtalInBCEnergy[x][y] > seedEnergy ) {                         
                      seedEnergy= xtalInBCEnergy[x][y] ;
                      seedTime  =   xtalTime ;
		      seedTErr   = xtalInBCTimeErr[x][y] ;
                   }
		   //double sinTheta   = fabs( sin( 2. *atan( exp(-1*xtalInBCEta[x][y]  ) ) ) );
		   nXtal++ ;
                   hJets.j0_xtalEta_ADC->Fill( xtalInBCEta[x][y], xtalADC[x][y] );
		   if ( j == 0 && fabs(xtalInBCEta[x][y])  < 1.479 ) hJets.j0_xtalEB_ADC_T->Fill( xtalADC[x][y], xtalTime );
                   if ( j == 0 && fabs(xtalInBCEta[x][y])  < 1.479 ) hJets.j0_xtalEB_ADC_E->Fill( xtalADC[x][y], xtalInBCEnergy[x][y] );
                   if ( j == 0 && fabs(xtalInBCEta[x][y]) >= 1.479 ) hJets.j0_xtalEE_ADC_T->Fill( xtalADC[x][y], xtalTime );
                   if ( j == 0 && fabs(xtalInBCEta[x][y]) >= 1.479 ) hJets.j0_xtalEE_ADC_E->Fill( xtalADC[x][y], xtalInBCEnergy[x][y] );
               }
               if ( seedEnergy != -999999 ) {
      	          jxTime += seedTime/ pow( seedTErr, 2 ) ;
	          nC    += 1./ pow( seedTErr, 2) ;
               }
	       if ( xE > 0 ) cE += clusterEnergy[x] ;
           }
	   float emF    = xE / jp4.E() ;
	   if ( nC > 0. && nXtal > 0 ) { 
              k++ ;
	      jxTime  = jxTime / nC ;
	      hJets.jet_emF_T->Fill( emF , jxTime ) ;
	      hJets.xE_cE->Fill( xE , cE ) ;
	      hJets.jet_Pt_nXtal->Fill( jp4.Pt() , nXtal ) ;
	      if ( j == 0 ) { 
	         if ( jp4.Eta()  < 1.479 ) hJets.jet0xTB->Fill( jxTime ); 
		 if ( jleadEta   < 1.479 ) hJets.jet0lTB->Fill( jleadTime ); 
	         if ( jp4.Eta() >= 1.479 ) hJets.jet0xTE->Fill( jxTime ); 
		 if ( jleadEta  >= 1.479 ) hJets.jet0lTE->Fill( jleadTime ); 
		 hJets.jet0Pt_T->Fill( jp4.Pt() , jxTime ) ;
		 hJets.jet0Eta_T->Fill( jp4.Eta() , jxTime ) ;
	      }
              if ( j == 1 ) {
	         hJets.jet1xT->Fill( jxTime );
    	         hJets.jet1lT->Fill( jleadTime ); 
              }
              if ( k == 1 ) jTime = jxTime ;
	   } 
       } // end of jet loop
       if ( k > 0 && ngC > 0. ) hJets.dT_JetPho->Fill( gleadTime - jTime ) ;  
       
       //** 3. get the associated bc for electron - a cross check 
       for ( int j=0 ; j< (int)(eleV.size()); j++) {
           TLorentzVector ep4( eleV[j].Px(), eleV[j].Py(), eleV[j].Pz(), eleV[j].E() ) ;

           // loop the clusters
	   float eleadEnergy = 0 ;
	   float eleadTime   = 0 ;
           float eleadEta    = 999 ;
	   int   nXtal       = 0 ; 
	   for (int x =0; x< nClusters; x++) {

	       if ( CPIdx[x] < 11. || CPIdx[x] > 12. ) continue ;

	       //float jidx =  11. + (j+1.)*0.1 ;
	       if ( clusterEnergy[x] < BasicClusterCuts[0] ) continue;
	       if ( CPIdx[x] != eIdV[j]  || fabs(clusterTime[x]) > BasicClusterCuts[1] ) continue ;
	       //cout<<" BC"<<x<<" CPIdx = "<< CPIdx[x]  <<endl ;

	       for (int y=0; y< nXtalInBC[x]; y++) {

                   double xtalTime  =  xtalInBCTime[x][y] 
                       + (xtalInBCTime[x][y]*theTimeCorrector_.getCorrection( (float) xtalInBCEnergy[x][y], xtalInBCEta[x][y]) );
                   if ( doTimeCorrection == 0 )  xtalTime = xtalInBCTime[x][y] ;

                   if ( xtalInBCEnergy[x][y] < XtalCuts[0] ) continue ;
		   if ( xtalInBCTimeErr[x][y] < 0.2 || xtalInBCTimeErr[x][y] > XtalCuts[1] ) continue ;
                   if ( fabs(xtalInBCEta[x][y])  < 1.479 && xtalChi2[x][y] > XtalCuts[2] ) continue ;
                   if ( fabs(xtalInBCEta[x][y]) >= 1.479 && xtalChi2[x][y] > XtalCuts[3] ) continue ;
		   //double dRxj = DeltaR( xtalInBCEta[x][y], xtalInBCPhi[x][y] , eleV[j].Eta(), eleV[j].Phi()  ) ;
		   //if ( dRxj > XtalCuts[4] ) continue ;
                   nXtal++ ;
		   if ( xtalInBCEnergy[x][y] > eleadEnergy ) {                         
		      eleadEnergy = xtalInBCEnergy[x][y] ;
		      eleadTime   = xtalTime ;
                      eleadEta    = xtalInBCEta[x][y] ;
                   }
               }
           }
	   if ( nXtal > 0 ) { 
	       if ( eleadEta   < 1.479 ) hJets.elelTB->Fill( eleadTime ); 
	       if ( eleadEta  >= 1.479 ) hJets.elelTE->Fill( eleadTime ); 
	   } 
           //if ( nXtal == 0 ) cout<<" sucks! bad electron :"<< ep4.Pt() <<endl ;
       } // end of electron loop

   } // end of event looping

   hJets.Draw( hfolder, plotType ) ;
   hJets.FitNDraw( hfolder, plotType ) ;
}

void PhotonAna::ScalarPlotList() {

   bool debugPlots = false ;
   ScalarPlotter( hJets.g0_xtalEB_ADC_T, "g0_xtalEB_ADC_T", -1., 1., 5, 1, debugPlots ) ;
   ScalarPlotter( hJets.j0_xtalEB_ADC_T, "j0_xtalEB_ADC_T", -1., 1., 5, 1 ) ;
   ScalarPlotter( hJets.j0_xtalEE_ADC_T, "j0_xtalEE_ADC_T", -1., 1., 5, 1 ) ;
   ScalarPlotter( hJets.j0_xtalEta_ADC,  "j0_xtalEta_ADC",   0., 1000., 1  ) ;
   ScalarPlotter( hJets.jet0Eta_T,    "jet_Eta_Time",  -1,   1 ) ;
   ScalarPlotter( hJets.jet0Pt_T,     "jet_Pt_Time",   -1,   1, 1, 1 ) ;

}

void PhotonAna::ScalarPlotter( TH2D* h2, TString hname, double yMin, double yMax, int rbin, int drawOpt, bool debugPlots  ) {

    TCanvas* c1 = new TCanvas("c1","", 800, 600);
    c1->SetFillColor(10);
    c1->SetFillColor(10);
    if ( drawOpt == 1) c1->SetLogx();
    if ( drawOpt == 2) c1->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();

    vector<double> xV ;
    vector<double> yV ;
    vector<double> yWidthV ;
    vector<double> yErrV ;
    BinningFitScan( h2, xV, yV, yWidthV, yErrV, debugPlots, rbin ) ; 

    const int sz = static_cast<const int>( xV.size() );
    float  x[sz] , xErr[sz], y[sz], yWidth[sz], yErr[sz] ;
    float yErrMax = 0 ;
    for ( int i=0; i< sz; i++) {
        x[i] = xV[i] ;
        y[i] = yV[i] ;
        xErr[i] = 0.  ;
        yWidth[i] = yWidthV[i] ;
        yErr[i] = yErrV[i] ;
        yErrMax = ( yErr[i] > yErrMax ) ? yErr[i] : yErrMax ;
    }

    c1->cd();
    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptFit(111);

    //TGraphErrors* g1 = new TGraphErrors( sz, x, y,  xErr, yWidth );
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


void PhotonAna::BinningFitScan( TH2D* h2, vector<double>& xV, vector<double>& yV, vector<double>& yWidthV, 
                                vector<double>& yErrV ,bool debugPlots, int rbin, int startBin, int finalBin ) {

    int      sz = h2->GetNbinsX();
    for ( int i=0; i< sz; i++) {
        int b1 = startBin + (i*rbin) ;
        int b2 = b1 + rbin - 1 ;
        int bcen = (b1 + b2) / 2 ; 
        if ( b1 == finalBin || b1 > sz ) break;
        vector<double> results = BinningFit( h2, "h2Py", b1, b2, debugPlots ) ;

        if ( results[2] == 0. ) continue ;
        xV.push_back( h2->GetBinCenter(bcen) ) ;
        yV.push_back( results[1] ) ;
        yWidthV.push_back( results[2] ) ;
        yErrV.push_back( results[3] ) ;
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

   h1_pjy->Fit(func, "RQ", "", yMin, yMax);
 
   double amp   = func->GetParameter(0) ;
   double mean  = func->GetParameter(1) ;
   double sigma = func->GetParameter(2) ;
   double n = h1_pjy->Integral() ;

   vector<double> result ;
   result.push_back( amp );
   result.push_back( mean );
   result.push_back( sigma  );
   result.push_back( sigma / sqrt(n) );

   TString plotname0 ;
   if ( debugPlots ) {
      //cout <<" debugging slice&fit"<< endl; 
      gSystem->cd( hfolder.c_str() ) ;
      gSystem->mkdir( "debug" );
      gSystem->cd( "../" );
     
      gStyle->SetOptStat("ieroum") ;
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

