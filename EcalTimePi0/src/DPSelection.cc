#include "ECALTime/EcalTimePi0/interface/DPSelection.h"

DPSelection::DPSelection( string datacardfile ){

  // SC's getParameters method
  // If you don't like to use the Datacard.txt , markout this section and use CMSSW method
  Input = new AnaInput( datacardfile );
  
  Input->GetParameters( "VertexCuts",   &vtxCuts );
  Input->GetParameters( "PhotonCuts",   &photonCuts );
  Input->GetParameters( "PhotonIso",    &photonIso );
  Input->GetParameters( "ElectronCuts", &electronCuts );
  Input->GetParameters( "JetCuts",      &jetCuts );
  Input->GetParameters( "MuonCuts",     &muonCuts );

  /*
  // CMSSW getParameter method
  vtxCuts              = iConfig.getParameter<std::vector<double> >("vtxCuts");
  jetCuts              = iConfig.getParameter<std::vector<double> >("jetCuts");
  photonCuts           = iConfig.getParameter<std::vector<double> >("photonCuts");
  photonIso            = iConfig.getParameter<std::vector<double> >("photonIso");
  electronCuts         = iConfig.getParameter<std::vector<double> >("electronCuts");
  */

}

DPSelection::~DPSelection(){

  delete Input;
  cout<<" done with selection ! "<<endl ;

}

void DPSelection::Init( TTree* tr ) { 

   tr->SetBranchAddress("eventId",    &eventId);
   tr->SetBranchAddress("trgCut",     &trgCut);
   tr->SetBranchAddress("nJets",      &nJets);
   tr->SetBranchAddress("nPhotons",   &nPhotons);
   tr->SetBranchAddress("nElectrons", &nElectrons);
   tr->SetBranchAddress("nVertices",  &nVertices );
   tr->SetBranchAddress("nSuperClusters", &nSuperClusters);
   tr->SetBranchAddress("nClusters",  &nClusters);

   tr->SetBranchAddress("metPx",      &metPx );
   tr->SetBranchAddress("metPy",      &metPy );
   tr->SetBranchAddress("met",        &metE );

   tr->SetBranchAddress("jetPx",       jetPx );
   tr->SetBranchAddress("jetPy",       jetPy );
   tr->SetBranchAddress("jetPz",       jetPz );
   tr->SetBranchAddress("jetE",        jetE );
   tr->SetBranchAddress("jetNDau",     jetNDau );
   tr->SetBranchAddress("jetCM",       jetCM );
   tr->SetBranchAddress("jetCEF",      jetCEF );
   tr->SetBranchAddress("jetNHF",      jetNHF );
   tr->SetBranchAddress("jetNEF",      jetNEF );

   tr->SetBranchAddress("phoPx",       phoPx );
   tr->SetBranchAddress("phoPy",       phoPy );
   tr->SetBranchAddress("phoPz",       phoPz );
   tr->SetBranchAddress("phoE",        phoE );
   tr->SetBranchAddress("phoEcalIso",  phoEcalIso );
   tr->SetBranchAddress("phoHcalIso",  phoHcalIso );
   tr->SetBranchAddress("phoTrkIso",   phoTrkIso );
   tr->SetBranchAddress("phoHovE",     phoHovE );
   tr->SetBranchAddress("phoSmin",     phoSmin );
   tr->SetBranchAddress("phoSmaj",     phoSmaj );
   tr->SetBranchAddress("phoTime",     phoTime );

   tr->SetBranchAddress("elePx",        elePx );
   tr->SetBranchAddress("elePy",        elePy );
   tr->SetBranchAddress("elePz",        elePz );
   tr->SetBranchAddress("eleE",         eleE );
   tr->SetBranchAddress("eleEcalIso",   eleEcalIso );
   tr->SetBranchAddress("eleHcalIso",   eleHcalIso );
   tr->SetBranchAddress("eleTrkIso",    eleTrkIso );
   tr->SetBranchAddress("eleNLostHits", eleNLostHits );

   tr->SetBranchAddress("muPx",        muPx );
   tr->SetBranchAddress("muPy",        muPy );
   tr->SetBranchAddress("muPz",        muPz );
   tr->SetBranchAddress("muE",         muE );
   tr->SetBranchAddress("muEcalIso",   muEcalIso );
   tr->SetBranchAddress("muHcalIso",   muHcalIso );
   tr->SetBranchAddress("muTrkIso",    muTrkIso );

   tr->SetBranchAddress("vtxX",       vtxX );
   tr->SetBranchAddress("vtxY",       vtxY );
   tr->SetBranchAddress("vtxZ",       vtxZ );
   tr->SetBranchAddress("vtxChi2",    vtxChi2 );
   tr->SetBranchAddress("vtxNdof",    vtxNdof );
}

// analysis template
bool DPSelection::HLTFilter( ) { 
    
     bool pass = true  ; // temporarily turn off the HLT selection
     if ( trgCut ==  75 ) pass = true ; 
     if ( trgCut ==  90 ) pass = true ; 
     if ( trgCut ==  165 ) pass = true ; 
     return pass ;
}


bool DPSelection::PhotonFilter( bool doIso ) { 

       bool pass =  true ;

       // 0. photon cuts
       phoV.clear() ;
       for ( int j=0 ; j< nPhotons; j++ ) {
           TLorentzVector phoP4( phoPx[j], phoPy[j], phoPz[j], phoE[j] ) ;
           if (        phoP4.Pt() < photonCuts[0] )          continue ;
           if ( fabs(phoP4.Eta()) > photonCuts[1] )          continue ;

           if ( ( phoHcalIso[j]/phoP4.Pt() + phoHovE[j] ) * phoE[j] >= photonCuts[2] )      continue ;
           if ( phoSmin[j] <= photonCuts[5] || phoSmin[j] >= photonCuts[6] )                continue ;
 
           if ( phoTime[j] < photonCuts[7] ) continue ; 
           // Isolation
           if ( doIso ) {
              if ( phoTrkIso[j] / phoP4.Pt()  >= photonIso[0] )                                continue ;
              if ( phoEcalIso[j] >= photonIso[1] || phoEcalIso[j] / phoE[j] >= photonIso[2] )  continue ;
              if ( phoHcalIso[j] >= photonIso[3] || phoHcalIso[j] / phoE[j] >= photonIso[4] )  continue ;
           }

           // check the isolation -- using dR_gj
           double dR_gj = 999. ;
           for ( size_t k=0 ; k< jetV.size() ; k++ ) {
               if ( phoP4.DeltaR( jetV[k].second ) < dR_gj )  dR_gj  = phoP4.DeltaR( jetV[k].second ) ;
           }
           if ( dR_gj < photonCuts[3] ) continue ;
           
 
           phoV.push_back( make_pair( j , phoP4) );
       }
       if ( (int)phoV.size() < photonCuts[4] ) pass = false ;

       return pass ;
}

bool DPSelection::VertexFilter() { 

     bool pass =  true ;
     // 1. vertex cuts
     int nVtx = 0 ;
     for ( int j=0 ; j< nVertices; j++ ) {

         double vtxRho = sqrt( (vtxX[j]*vtxX[j]) + (vtxY[j]*vtxY[j]) ); 
	 if ( nVertices < 1 )                continue ;
	 if ( vtxNdof[j]     < vtxCuts[0] )  continue ;
	 if ( fabs(vtxZ[j]) >= vtxCuts[1] )  continue ;
	 if ( vtxRho        >= vtxCuts[2] )  continue ;
         nVtx++ ;
     }
     if ( nVtx < 1 ) pass = false ;
     return pass ;
}

bool DPSelection::JetMETFilter( bool useJetID ) { 

     bool pass =  true ;
     // 1. jet selection
     jetV.clear() ;
     for ( int j=0 ; j< nJets; j++ ) {
         TLorentzVector jp4( jetPx[j], jetPy[j], jetPz[j], jetE[j] ) ;
	 if ( jp4.Pt()        < jetCuts[0] ) continue ;
	 if ( fabs(jp4.Eta()) > jetCuts[1] ) continue ;
           
         // Jet ID cuts
         if ( useJetID ) {
            if ( jetNDau[j] <    2 )  continue ;
	    if ( jetCEM[j] >= 0.99 )  continue ;
	    if ( jetNHF[j] >= 0.99 )  continue ;
	    if ( jetNEF[j] >= 0.99 )  continue ;
	    if ( fabs( jp4.eta() ) < 2.4 && jetCM[j]  <=0 ) continue ;
         } 

	 double dR_gj = 999. ;
	 for ( size_t k=0; k< phoV.size() ; k++) {
	     if ( phoV[k].second.DeltaR( jp4 ) < dR_gj )  dR_gj  = phoV[k].second.DeltaR( jp4 ) ;
	 }
	 if ( dR_gj < photonCuts[3] ) continue ;
	 jetV.push_back( make_pair( j, jp4 ) );

     }
     int nu_Jets = (int)jetV.size() ;
     if ( nu_Jets < jetCuts[2] || nu_Jets > jetCuts[3] )      pass = false ;

     // 2. met selection
     TLorentzVector metp4( metPx, metPy, 0., metE ) ;
     if ( jetCuts[4] >= 0 &&  metp4.Et() < jetCuts[4] )         pass = false ;
     if ( jetCuts[4]  < 0 &&  metp4.Et() > fabs( jetCuts[4] ) ) pass = false ;
 
     return pass ;
}

bool DPSelection::ElectronFilter() { 

     bool pass =  true ;
     eleV.clear() ;
     for ( int j=0 ; j< nElectrons; j++ ) {
         TLorentzVector eP4( elePx[j], elePy[j], elePz[j], eleE[j] ) ;
	 if ( eP4.Pt() < electronCuts[0] )          continue ;
	 if ( fabs( eP4.Eta()) > electronCuts[1] )  continue ;

         double relIso   = (eleEcalIso[j] + eleHcalIso[j] + eleTrkIso[j] ) / eP4.Pt() ;

         if ( relIso > electronCuts[2] )            continue ;
         //if ( eleNLostHits[j] >= electronCuts[4]  ) continue ;

         double dR_ej = 999. ;
         for ( size_t k=0 ; k< jetV.size() ; k++ ) {
             if ( eP4.DeltaR( jetV[k].second ) < dR_ej )  dR_ej  = eP4.DeltaR( jetV[k].second ) ;
         }
         if ( dR_ej < electronCuts[3] ) continue ;

         eleV.push_back( make_pair( j, eP4.Pt()) ) ;
     }
     if ( eleV.size() < 1 ) pass = false ;
     return pass = false ;
}

bool DPSelection::MuonFilter() { 

     bool pass =  true ;
     muV.clear() ;
     for ( int j=0 ; j< nMuons; j++ ) {
         TLorentzVector mP4( muPx[j], muPy[j], muPz[j], muE[j] ) ;
	 if ( mP4.Pt() < muonCuts[0] )          continue ;
	 if ( fabs( mP4.Eta()) > muonCuts[1] )  continue ;

         double relIso   = (muEcalIso[j] + muHcalIso[j] + muTrkIso[j] ) / mP4.Pt() ;
         if ( relIso > muonCuts[2] )            continue ;

         double dR_mj = 999. ;
         for ( size_t k=0 ; k< jetV.size() ; k++ ) {
             if ( mP4.DeltaR( jetV[k].second ) < dR_mj )  dR_mj  = mP4.DeltaR( jetV[k].second ) ;
         }
         if ( dR_mj < muonCuts[3] ) continue ;

         muV.push_back( make_pair( j, mP4.Pt() ) ) ;
     }
     if ( muV.size() < 1 ) pass = false ;
     return pass = false ;
}

bool DPSelection::GammaJetsBackground( ) {

    bool pass = true ;

    if ( jetV.size() < 1 ) return false ;
    if ( phoV.size() < 1 ) return false ;

    double dR   = phoV[0].second.DeltaR( jetV[0].second ) ;
    double ratio1 = jetV[0].second.Pt() / phoV[0].second.Pt() ;

    if ( dR <= (2*3.141593/3) )           pass = false ;
    if ( ratio1 >= 1.3 || ratio1 <= 0.7 ) pass = false ;

    if ( jetV.size() > 1 ) {
       double ratio2 = jetV[1].second.Pt() / phoV[0].second.Pt() ;
       if ( ratio2 >= 0.1 ) pass = false ;
    }

    return pass ;

}

bool DPSelection::GammaJetsControlSample( bool isTightPhoton ) {

       bool passHLT  = HLTFilter();
       bool passVtx  = VertexFilter();

       // reset cuts to tight photon selection
       if ( isTightPhoton ) {
          ResetCuts( "PhotonCuts", 0, 100. ) ;  // pt
	  ResetCuts( "PhotonCuts", 1, 1.4 ) ;   // eta
	  ResetCuts( "PhotonCuts", 7, -2. ) ;   // seed time
	  ResetCuts( "PhotonIso",  0, 0.1 ) ;   // Trk Iso
	  ResetCuts( "PhotonIso",  1, 2.4 ) ;   // Ecal Et 
	  ResetCuts( "PhotonIso",  2, 0.05 ) ;  // Ecal Ratio
	  ResetCuts( "PhotonIso",  3, 2.4 ) ;   // Hcal Et
	  ResetCuts( "PhotonIso",  4, 0.05 ) ;  // Hcal Ratio
       }
       bool passPho = PhotonFilter( true );

       ResetCuts( "JetCuts",  2, 1 ) ;  // Max Num of Jets
       ResetCuts( "JetCuts",  3, 3 ) ;  // Max Num of Jets
       bool passJet = JetMETFilter();

       bool passGJets = GammaJetsBackground() ;

       ResetCuts() ;  // reset cuts from Datacard

       bool isGJet = ( passHLT && passVtx  && passPho && passJet && passGJets ) ? true : false ;
       return isGJet ;
}

bool DPSelection::SignalSelection( bool isTightPhoton ) {

       bool passHLT  = HLTFilter();
       bool passVtx  = VertexFilter();

       // reset cuts to tight photon selection
       if ( isTightPhoton ) {
          ResetCuts( "PhotonCuts", 0, 100. ) ;  // pt
	  ResetCuts( "PhotonCuts", 1, 1.4 ) ;   // eta
	  ResetCuts( "PhotonCuts", 7, -2. ) ;   // seed time
	  ResetCuts( "PhotonIso",  0, 0.1 ) ;   // Trk Iso
	  ResetCuts( "PhotonIso",  1, 2.4 ) ;   // Ecal Et 
	  ResetCuts( "PhotonIso",  2, 0.05 ) ;  // Ecal Ratio
	  ResetCuts( "PhotonIso",  3, 2.4 ) ;   // Hcal Et
	  ResetCuts( "PhotonIso",  4, 0.05 ) ;  // Hcal Ratio
       }
       bool passPho = PhotonFilter( true );

       ResetCuts( "JetCuts",  2, 3 ) ;  // Min Num of Jets
       bool passJet = JetMETFilter();

       bool passGJets = GammaJetsBackground() ;

       ResetCuts() ;  // reset cuts from Datacard

       bool isSignal = ( passHLT && passVtx  && passPho && passJet && !passGJets ) ? true : false ;
       return isSignal ;

}

bool DPSelection::QCDControlSample() {

       bool passHLT  = HLTFilter();
       bool passVtx  = VertexFilter();
       // reset cuts to tight photon selection
       ResetCuts( "PhotonCuts", 0, 100. ) ;  // pt
       ResetCuts( "PhotonCuts", 1, 1.4 ) ;   // eta
       ResetCuts( "PhotonCuts", 7, -2. ) ;   // seed time
       ResetCuts( "PhotonIso",  0, 0.1 ) ;   // Trk Iso
       ResetCuts( "PhotonIso",  1, 2.4 ) ;   // Ecal Et 
       ResetCuts( "PhotonIso",  2, 0.05 ) ;  // Ecal Ratio
       ResetCuts( "PhotonIso",  3, 2.4 ) ;   // Hcal Et
       ResetCuts( "PhotonIso",  4, 0.05 ) ;  // Hcal Ratio
       bool passTight = PhotonFilter( true );
 
       ResetCuts( "PhotonCuts", 0, 70. ) ;  // pt
       ResetCuts( "PhotonCuts", 1, 9.9 ) ;  // eta
       ResetCuts( "PhotonIso",  0, 0.2 ) ;  // Trk Iso
       ResetCuts( "PhotonIso",  1, 4.5 ) ;  // Ecal Et 
       ResetCuts( "PhotonIso",  2, 0.1 ) ;  // Ecal Ratio
       ResetCuts( "PhotonIso",  3, 4.0 ) ;  // Hcal Et
       ResetCuts( "PhotonIso",  4, 0.1 ) ;  // Hcal Ratio
       bool passPho = PhotonFilter( true );

       ResetCuts( "JetCuts",  2, 3 ) ;  // Min Num of Jets
       bool passJet = JetMETFilter();

       //bool passGJets = GammaJetsBackground() ;

       ResetCuts() ;  // reset cuts from Datacard

       bool isQCD = ( passHLT && passVtx  && passPho && passJet && !passTight ) ? true : false ;

       return isQCD ;
}

void DPSelection::ResetCuts( string cutName, vector<int>& cutId, vector<double>& newValue ) {

     for ( size_t i=0; i< cutId.size() ; i++ ) {
         if ( cutName == "PhotonCuts" )   photonCuts[ cutId[i] ]   = newValue[i] ;
         if ( cutName == "PhotonIso" )    photonIso[ cutId[i] ]    = newValue[i] ;
         if ( cutName == "VertexCuts" )   vtxCuts[ cutId[i] ]      = newValue[i] ;
         if ( cutName == "ElectronCuts" ) electronCuts[ cutId[i] ] = newValue[i] ;
         if ( cutName == "JetCuts" )      jetCuts[ cutId[i] ]      = newValue[i] ;
         if ( cutName == "MuonCuts" )     muonCuts[ cutId[i] ]     = newValue[i] ;
     }

}

void DPSelection::ResetCuts( string cutName, int cutId, double newValue ) {

         if ( cutName == "PhotonCuts" )   photonCuts[ cutId ]   = newValue ;
         if ( cutName == "PhotonIso" )    photonIso[ cutId ]    = newValue ;
         if ( cutName == "VertexCuts" )   vtxCuts[ cutId ]      = newValue ;
         if ( cutName == "ElectronCuts" ) electronCuts[ cutId ] = newValue ;
         if ( cutName == "JetCuts" )      jetCuts[ cutId ]      = newValue ;
         if ( cutName == "MuonCuts" )     muonCuts[ cutId ]     = newValue ;

}

void DPSelection::ResetCuts( string cutName ) {

    if ( cutName == "VertexCuts"   || cutName == "All" ) Input->GetParameters( "VertexCuts",   &vtxCuts );
    if ( cutName == "PhotonCuts"   || cutName == "All" ) Input->GetParameters( "PhotonCuts",   &photonCuts );
    if ( cutName == "PhotonIso"    || cutName == "All" ) Input->GetParameters( "PhotonIso",    &photonIso );
    if ( cutName == "ElectronCuts" || cutName == "All" ) Input->GetParameters( "ElectronCuts", &electronCuts );
    if ( cutName == "JetCuts"      || cutName == "All" ) Input->GetParameters( "JetCuts",      &jetCuts );
    if ( cutName == "MuonCuts"     || cutName == "All" ) Input->GetParameters( "MuonCuts",     &muonCuts );

}

void DPSelection::ResetCollection( string cutName ) {

    if ( cutName == "Photon"   || cutName == "All" ) phoV.clear() ;
    if ( cutName == "Electron" || cutName == "All" ) eleV.clear() ;
    if ( cutName == "Jet"      || cutName == "All" ) jetV.clear() ;
    if ( cutName == "Muon"     || cutName == "All" ) muV.clear() ;

}

void DPSelection::GetCollection( string collName, vector<objID>& coll ) {

   if ( collName == "Photon" ) {
      for ( size_t i=0; i< phoV.size() ; i++ ) coll.push_back(  phoV[i] )  ;
   }
   else if ( collName == "Jet" ) {
      for ( size_t i=0; i< jetV.size() ; i++ ) coll.push_back( jetV[i] )  ;
   }
   else if ( collName == "Muon" ) {
      for ( size_t i=0; i< muV.size() ; i++ ) coll.push_back( muV[i] )  ;
   }
   else if ( collName == "Electron" ) {
      for ( size_t i=0; i< eleV.size() ; i++ ) coll.push_back( eleV[i] )  ;
   }
   else {
      cout <<" no collection matched ! " <<endl ;
   }

}

