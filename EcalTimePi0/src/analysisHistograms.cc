#include "ECALTime/EcalTimePi0/interface/analysisHistograms.h"

void HistSet::setTree(EcalTimePhyTreeContent * treeVars){
  treeVars_ = treeVars;
}


void HistSet::book(TFileDirectory subDir, const std::string& post) {
  
  treeVars_ = 0;

  nVertices_=subDir.make<TH1F>("num vertices","num vertices; num vertices",41,-0.5,40.5);
  mass_               =(TH1F*) subDir.make<TH1F>("mass","mass; m(ele,ele) [GeV]",80,50,130);
  dZvertices_         =(TH1F*) subDir.make<TH1F>("dZvertices","dZvertices; #DeltaZ(ele_{1},ele_{2}) [cm]",250,0,25);
  redAmplitudePair_   =(TH1F*) subDir.make<TH1F>("redAmplitudePair","redAmplitudePair; A_{1}A_{2}/#sqrt(A_{1}^{2}+A_{2}^{2})  [ADC]",500,0,5000);
  Zvertices_          =(TH1F*) subDir.make<TH1F>("Zvertices","Zvertices; z vertex [cm]",250,-25,25);

  // Initialize histograms -- xtals
  chi2_                =(TH1F*) subDir.make<TH1F>("cluster chi2 ","cluster chi2 ; #chi^{2}",100,0,10);
  chi2ndf_             =(TH1F*) subDir.make<TH1F>("cluster chi2 per d.o.f.","cluster chi2 ; #chi^{2}/N_{dof}",100,0,10);

  seedTime_            =(TH1F*) subDir.make<TH1F>("seed time","seed time; t_{seed} [ns]; num. seeds/0.05ns",binsTDistro_*incrementForDelayedGamma_,-rangeTDistro_*incrementForDelayedGamma_,rangeTDistro_*incrementForDelayedGamma_);
  seedTimeTOF_            =(TH1F*) subDir.make<TH1F>("seed time TOF-corr","seed time TOF-corr; t_{seed} [ns]; num. seeds/0.05ns",binsTDistro_*incrementForDelayedGamma_,-rangeTDistro_*incrementForDelayedGamma_,rangeTDistro_*incrementForDelayedGamma_);
  secondTime_          =(TH1F*) subDir.make<TH1F>("second time","second time; t_{second} [ns]; num. secs/0.05ns",binsTDistro_*incrementForDelayedGamma_,-rangeTDistro_*incrementForDelayedGamma_,rangeTDistro_*incrementForDelayedGamma_);
  secondTimeTOF_          =(TH1F*) subDir.make<TH1F>("second time TOF-corr","second time TOF-corr; t_{second} [ns]; num. secs/0.05ns",binsTDistro_*incrementForDelayedGamma_,-rangeTDistro_*incrementForDelayedGamma_,rangeTDistro_*incrementForDelayedGamma_);
  clusterTime_         =(TH1F*) subDir.make<TH1F>("cluster time","cluster time; t_{cluster} [ns]; num. clusters/0.05ns",binsTDistro_*incrementForDelayedGamma_,-rangeTDistro_*incrementForDelayedGamma_,rangeTDistro_*incrementForDelayedGamma_);
  clusterTimeTOF_         =(TH1F*) subDir.make<TH1F>("cluster time TOF-corr","cluster time TOF-corr; t_{cluster} [ns]; num. clusters/0.05ns",binsTDistro_*incrementForDelayedGamma_,-rangeTDistro_*incrementForDelayedGamma_,rangeTDistro_*incrementForDelayedGamma_);

  TOFsingle_           = (TH1F*) subDir.make<TH1F>("TOF single","TOF single correction; #Delta TOF [ns]; num. clusters/0.05ns",binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);
  TOFcorrections_      = (TH1F*) subDir.make<TH1F>("TOF difference","TOF difference; #Delta TOF [ns]; num. clusters/0.05ns",binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);
  TOFcorrectionsVSdeltaEta_=(TH2F*) subDir.make<TH2F>("TOF corrections VS #Delta#eta","TOF corrections VS #Delta#eta; #Delta#eta; #Delta TOF [ns]; ",30,0,3.,binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);
  clusTimeDiffHistTOFVSdeltaEtaRightVertex_=(TH2F*) subDir.make<TH2F>("TOF-corr cluster time difference VS #Delta#eta RightVertex","TOF-corr cluster time difference VS #Delta#eta RightVertex; |#Delta#eta|; (t_{clus1} - t_{clus2}) TOF-corrected [ns]; ",30,0,3.,binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);
  clusTimeDiffHistTOFVSdeltaEtaWrongVertex_=(TH2F*) subDir.make<TH2F>("TOF-corr cluster time difference VS #Delta#eta WrongVertex","TOF-corr cluster time difference VS #Delta#eta WrongVertex; |#Delta#eta|;  (t_{clus1} - t_{clus2}) TOF-corrected [ns]; ",30,0,3.,binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);

  zPairVsTPair_        =(TH2F*) subDir.make<TH2F>("z_{pair} VS t_{coll}","z_{pair} VS t_{coll}; t_{coll} [ns]; z_{pair}  [cms]",75,-1.5,1.5,50,-25,25);

  seed2secSingleClus_  =(TH1F*) subDir.make<TH1F>("single cluster: t_{seed}-t_{second}","single cluster: t_{seed}-t_{second}; t_{seed} - t_{second} [ns]; num. clusters/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);

  tCollVSdeltaEtaRightVertex_=(TH2F*) subDir.make<TH2F>("t_{coll} VS #Delta#eta RightVertex","t_{coll} VS #Delta#eta RightVertex; |#Delta#eta|; t_{coll} [ns]; ",30,0,3.,binsTDistro_,-rangeTDistro_/4.,rangeTDistro_/4.);
  tCollVStimeDiffHistTOF_=(TH2F*) subDir.make<TH2F>("TOF-corrected: (t_{clus1} + t_{clus2})/2  VS  (t_{clus1} - t_{clus2})/2","TOF-corrected: (t_{clus1} + t_{clus2})/2  VS  (t_{clus1} - t_{clus2})/2 [ns]; (t_{clus1} - t_{clus2})/2 ; (t_{clus1} + t_{clus2})/2 VS [ns]; ",binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.,binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);

  seedTimeVSchi2_  =(TH2F*) subDir.make<TH2F>("t_{seed} VS #chi^{2}/N_{dof}","t_{seed} VS #chi^{2}/N_{dof}; t_{seed} [ns];  #chi^{2}/N_{dof}; ",binsTDistro_,-rangeTDistro_/4.,rangeTDistro_/4.,20,0,10);

  seedTimeVSseedMsecond_ =(TH2F*) subDir.make<TH2F>("t_{seed} VS (t_{seed}-t_{second})","t_{seed} VS  (t_{seed}-t_{second}); t_{seed} [ns];  (t_{seed}-t_{second}) [ns]; ",binsTDistro_,-rangeTDistro_/4.,rangeTDistro_/4.,binsTDistro_/2.,-rangeTDistro_/8.,rangeTDistro_/8.);

  seedTimeDiffHist_    =(TH1F*) subDir.make<TH1F>("time difference of seeds","seeds time difference; t_{seed1} - t_{seed2} [ns]; num. seed pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  seedTimeDiffHistTOF_ =(TH1F*) subDir.make<TH1F>("TOF-corr time difference of seeds","TOF-corr seed time difference; (t_{seed1} - t_{seed2}) TOF-corrected   [ns]; num. seed pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);

  secondTimeDiffHist_  =(TH1F*) subDir.make<TH1F>("time difference of seconds","second time difference; t_{second1} - t_{second2} [ns]; num. seed pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);//GF new
  secondTimeDiffHistTOF_ =(TH1F*) subDir.make<TH1F>("TOF-corr time difference of seconds","TOF-corr seconds time difference;  (t_{second1} - t_{second2}) TOF-corrected [ns]; num. sec pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);


 seed2secondTimeDiffHist_ =(TH1F*) subDir.make<TH1F>("time difference of seed1 to second2","seed2second time difference; t_{seed1} - t_{second2} [ns]; num. seed-second pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  seed2secondTimeDiffHistTOF_ =(TH1F*) subDir.make<TH1F>("TOF-corr time difference of seed1 to second2","TOF-corr seed2second time difference;  (t_{seed1} - t_{second2}) TOF-corrected [ns]; num. seed-second pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);


  clusTimeDiffHist_    =(TH1F*) subDir.make<TH1F>("cluster time difference","cluster time difference;  t_{clus1} - t_{clus2} [ns]; num. cluster pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  clusTimeDiffHistTOF_ =(TH1F*) subDir.make<TH1F>("TOF-corr cluster time difference","TOF-corr cluster time difference; (t_{clus1} - t_{clus2}) TOF-corrected [ns]; num. cluster pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);

  clusTimeDiffHistTOFwrongVertex_=(TH1F*) subDir.make<TH1F>("TOF-corr cluster time difference wrong vertex","TOF-corr cluster time difference wronge vertex; (t_{clus1} - t_{clus2}) TOF-corrected [ns]; num. cluster pairs/0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);

  tColl_     =(TH1F*) subDir.make<TH1F>("t_{coll}","t_{coll} [ns];  (t_{clus1} + t_{clus2})/2 [ns]",binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);
  tCollseeds_=(TH1F*) subDir.make<TH1F>("t_{coll-seeds}","t_{coll-seeds} [ns];  (t_{seed1} + t_{seed2})/2 [ns]",binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);
  tCollsecs_ =(TH1F*) subDir.make<TH1F>("t_{coll-seconds}","t_{coll-seconds} [ns];  (t_{second1} + t_{second2})/2 [ns]",binsTDistro_,-rangeTDistro_/2.,rangeTDistro_/2.);


  numCryBC1            =(TH1F*) subDir.make<TH1F>("num cry in bc1","num cry in bc1; num cry",25,0,25);
  numCryBC2            =(TH1F*) subDir.make<TH1F>("num cry in bc2","num cry in bc2; num cry",25,0,25);
  timeVsEta_           =(TH2F*) subDir.make<TH2F>("timeVsEta","timeVsEta;; #eta t [ns]",50,-2.5,2.5,150,-1.5,1.5);
  timeVsEtaLead_       =(TH2F*) subDir.make<TH2F>("timeVsEtaLead","timeVsEtaLead;#eta_{lead}; t [ns]",50,-2.5,2.5,150,-1.5,1.5);
  timeVsEtaSub_        =(TH2F*) subDir.make<TH2F>("timeVsEtaSub","timeVsEtaSub; #eta_{sublead}; t [ns]",50,-2.5,2.5,150,-1.5,1.5);
  outliersVsEtaPhi_    =(TH2F*) subDir.make<TH2F>("outliersVsEtaPhi","outliersVsEtaPhi; #eta; #phi",50,-2.5,2.5,72,-3.14,3.14);
  seedAmpli_           =(TH1F*) subDir.make<TH1F>("E(seed)  ","E(seed) ; E [GeV]",250,0,250);
  secondAmpli_         =(TH1F*) subDir.make<TH1F>("E(second)  ","E(second) ; E [GeV]",250,0,250);
  diffSeedOther_       =(TH1F*) subDir.make<TH1F>("t_{seed}-t_{others}","t_{seed}-t_{others}; t_{seed}-t_{others} [ns]; num./0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  diffSeedOtherOverErr_ =(TH1F*) subDir.make<TH1F>("(t_{seed}-t_{others})/#sigma","(t_{seed}-t_{others})/#sigma; (t_{seed}-t_{others})/#sigma; num./0.05ns",binsTDistro_*5./3,-rangeTDistro_*5./3,rangeTDistro_*5./3);

  diffSeedSecond_      =(TH1F*) subDir.make<TH1F>("t_{seed}-t_{second}","t_{seed}-t_{second}; t_{seed}-t_{second} [ns]; num./0.05ns",binsTDistro_,-rangeTDistro_,rangeTDistro_);
  diffSeedSecondOverErr_ =(TH1F*) subDir.make<TH1F>("(t_{seed}-t_{second})/#sigma","(t_{seed}-t_{second})/#sigma; (t_{seed}-t_{second})/#sigma; num./0.05ns",binsTDistro_*5./3,-rangeTDistro_*5./3,rangeTDistro_*5./3);
  seedVSSecond_        =(TH2F*) subDir.make<TH2F>("t_{seed} VS t_{second}","t_{seed} VS t_{second}; t_{seed} [ns]; t_{second} [ns]",75,-1.5,1.5,75,-1.5,1.5);

}
  
int HistSet::fill(int sc1, int sc2, int bc1, int bc2 ){
  // if launched with no arguments for cuts, 
  // call the generalized 'fil' setting type=0 and cut=0,
  // which means cuts are ignored 
  int type=0; float cut=0; 
  std::pair<float,float> dummyPhase = std::make_pair(0.,0.);
  return fill(sc1, sc2, bc1, bc2, type, cut, dummyPhase);
}

int HistSet::fill(int sc1, int sc2, int bc1, int bc2, std::pair<float,float>& globalPhases, timeCorrector& theCorrector ){
  // if launched with no arguments for cuts, 
  // call the generalized 'fil' setting type=0 and cut=0,
  // which means cuts are ignored 
  int type=0; float cut=0; 
  return fill(sc1, sc2, bc1, bc2, type, cut, globalPhases, theCorrector);
}

int HistSet::fill(int sc1, int sc2, int bc1, int bc2, std::pair<float,float>& globalPhases ){
  // if launched with no arguments for cuts, 
  // call the generalized 'fil' setting type=0 and cut=0,
  // which means cuts are ignored 
  int type=0; float cut=0; 
  return fill(sc1, sc2, bc1, bc2, type, cut, globalPhases);
}

int HistSet::fill(int sc1, int sc2, int bc1, int bc2, int type, float cut) {
  std::pair<float,float> dummyPhase = std::make_pair(0.,0.);
  return fill(sc1, sc2, bc1, bc2, type,  cut, dummyPhase );
}


int HistSet::fill(int sc1, int sc2, int bc1, int bc2, int type, float cut, std::pair<float,float>& globalPhases ){
  timeCorrector theCorrector;
  return fill(sc1, sc2, bc1, bc2, type, cut, globalPhases, theCorrector );
}

int HistSet::fill(int sc1, int sc2, int bc1, int bc2, int type, float cut, std::pair<float,float>& globalPhases, timeCorrector& theCorrector ){
  
  if(!treeVars_)   {
    std::cout << " treeVars_ is null when calling HistSet::fill! Bailing out."<< std::endl;
    assert(-1);
  }
  else             {if (0) {std::cout << " treeVars_ is not null:  " << treeVars_ << std::endl;}   }
  
  float et1 = treeVars_->superClusterRawEnergy[sc1]/cosh( treeVars_->superClusterEta[sc1] );
  math::PtEtaPhiELorentzVectorD  el1(et1  ,
				     treeVars_->superClusterEta[sc1],
				     treeVars_->superClusterPhi[sc1],
				     treeVars_->superClusterRawEnergy[sc1] );  
  float et2 = treeVars_->superClusterRawEnergy[sc2]/cosh( treeVars_->superClusterEta[sc2] );
  math::PtEtaPhiELorentzVectorD  el2(et2 ,
				     treeVars_->superClusterEta[sc2],
				     treeVars_->superClusterPhi[sc2],
				     treeVars_->superClusterRawEnergy[sc2] );
  math::PtEtaPhiELorentzVectorD diEle = el1;
  diEle += el2;
  
  // ////////////////////////
  mass_      ->Fill(diEle.M()); 
  float dvertex = pow(treeVars_->superClusterVertexZ[sc1]-treeVars_->superClusterVertexZ[sc2],2);

  dvertex       = sqrt(dvertex);
  dZvertices_->Fill(dvertex); 
  Zvertices_->Fill( (treeVars_->superClusterVertexZ[sc1]+treeVars_->superClusterVertexZ[sc2])/2 ); 
  nVertices_->Fill(treeVars_->nVertices); 


  float phase(0.);
  if( fabs(treeVars_->superClusterEta[sc1])<1.45 ) phase=globalPhases.first;
  else                                             phase=globalPhases.second;
  ClusterTime bcTime1 = timeAndUncertSingleCluster(bc1,phase,theCorrector,(*treeVars_));

  if( fabs(treeVars_->superClusterEta[sc2])<1.45 ) phase=globalPhases.first;
  else                                             phase=globalPhases.second;
  ClusterTime bcTime2 = timeAndUncertSingleCluster(bc2,phase,theCorrector,(*treeVars_));
  

  ///////////////////////////////////// cuts //////////////////////////////////////////////////////////////
  if (type==0){ // no cuts are being used
    ;
  }
  
  else if (type==1){ // cut straight on the agreement between .seedtime and .secondtime within either cluster
    if(bcTime1.second<0 || bcTime2.second<0) return 0;  // check that there's crystals beyond seed // ALL doubles 
    if(  fabs( bcTime1.seedtime - bcTime1.secondtime ) > cut 
	 ||
	 fabs( bcTime2.seedtime - bcTime2.secondtime ) > cut ) 
      return 0;   }
  
  else if (type==2){ // cut on the agreement between .seedtime and .otherstime within either cluster
    if(bcTime1.second<0 || bcTime2.second<0) return 0;  // check that there's crystals beyond seed // ALL doubles 
    if(  fabs( bcTime1.seedtime - bcTime1.otherstime ) > cut 
	 ||
	 fabs( bcTime2.seedtime - bcTime2.otherstime ) > cut ) 
      return 0;   }

  
  else if (type==3){ // cut on chi2 of either cluster
    if(   bcTime1.second<0 || bcTime2.second<0)           return 0;  // check that there's crystals beyond seed // ALL doubles 
    if(  (bcTime1.numCry<2)  	 ||  (bcTime2.numCry<2) ) return 0;
    
    if(  fabs( bcTime1.chi2/bcTime1.numCry ) > cut 
	 ||
	 fabs( bcTime2.chi2/bcTime2.numCry ) > cut ) 
      return 0;
  }

  else{
    std::cout << "+++ You're using HistSet with a cut type: " << type << " which is not implemented; bailing out" << std::endl;
    assert(-1);
  }
  ///////////////////////////////////// cuts are done////////////////////////////////////////////////////////  

  //std::cout << "analysisHisograms is launching fillSingle; bcTime1.seedtime: " << bcTime1.seedtime << std::endl;
  // FIRST, tale care of all single-cluster histograms; do for each of the two clusters  
  HistSet::fillSingle(sc1, bc1,  bcTime1, type, cut, globalPhases, theCorrector);
  HistSet::fillSingle(sc2, bc2,  bcTime2, type, cut, globalPhases, theCorrector);

  float A1 = treeVars_->xtalInBCAmplitudeADC[bc1][bcTime1.seed]; 
  float A2 = treeVars_->xtalInBCAmplitudeADC[bc2][bcTime2.seed]; 
  redAmplitudePair_ -> Fill( A1*A2 / sqrt( A1*A1 + A2*A2 ));

  TOFcorrections_           -> Fill(extraTravelTime(sc2,(*treeVars_)) - extraTravelTime(sc1,(*treeVars_)) ); // double
  TOFcorrectionsVSdeltaEta_ -> Fill(  fabs(treeVars_->superClusterEta[sc2]-treeVars_->superClusterEta[sc1]) , extraTravelTime(sc2,(*treeVars_)) - extraTravelTime(sc1,(*treeVars_))  ); // double

  int vtxOfThisEle=-99;
  // look for the vertex which electrons  are attached to: 
  for(int u=0; u<treeVars_->nVertices; u++){
    // matching done with 1mm tolerance
    if( fabs(treeVars_->superClusterVertexZ[sc2]-treeVars_->vtxZ[u]) < 0.1) {       vtxOfThisEle=u;     }
  }
  //std::cout << "\n\tdebugging vertices" << std::endl;
  if(vtxOfThisEle >-99){

    for(int u=0; u<treeVars_->nVertices; u++){
      if(u==vtxOfThisEle) {
	clusTimeDiffHistTOFVSdeltaEtaRightVertex_ -> Fill(  fabs(treeVars_->superClusterEta[sc2]-treeVars_->superClusterEta[sc1]) , 
							    (bcTime1.time-extraTravelTime(sc1,u,(*treeVars_)))  - (bcTime2.time-extraTravelTime(sc2,u,(*treeVars_)))
							    //(bcTime1.time-extraTravelTime(sc1,(*treeVars_)))  - (bcTime2.time-extraTravelTime(sc2,(*treeVars_)))
							    );

	tCollVSdeltaEtaRightVertex_-> Fill(  fabs(treeVars_->superClusterEta[sc2]-treeVars_->superClusterEta[sc1]) , 
					     ( (bcTime1.time-extraTravelTime(sc1,u,(*treeVars_))) + (bcTime2.time-extraTravelTime(sc2,u,(*treeVars_))) ) /2.
					     );
	tCollVStimeDiffHistTOF_    -> Fill(  ( (bcTime1.time-extraTravelTime(sc1,u,(*treeVars_))) - (bcTime2.time-extraTravelTime(sc2,u,(*treeVars_))) ) /2. ,
					     ( (bcTime1.time-extraTravelTime(sc1,u,(*treeVars_))) + (bcTime2.time-extraTravelTime(sc2,u,(*treeVars_))) ) /2.
					     );
	tColl_                     -> Fill( ( (bcTime1.time-extraTravelTime(sc1,u,(*treeVars_))) + (bcTime2.time-extraTravelTime(sc2,u,(*treeVars_))))  /2. );
	tCollseeds_                -> Fill( ( (bcTime1.seedtime-extraTravelTime(sc1,u,(*treeVars_))) + (bcTime2.seedtime-extraTravelTime(sc2,u,(*treeVars_))))  /2. );
	tCollsecs_                 -> Fill( ( (bcTime1.secondtime-extraTravelTime(sc1,u,(*treeVars_))) + (bcTime2.secondtime-extraTravelTime(sc2,u,(*treeVars_))))  /2. );

	zPairVsTPair_              -> Fill( ( (bcTime1.time-extraTravelTime(sc1,u,(*treeVars_))) + (bcTime2.time-extraTravelTime(sc2,u,(*treeVars_))) ) /2. ,
					    ( treeVars_->superClusterVertexZ[sc1] + treeVars_->superClusterVertexZ[sc2] )/2. ); 

      }// if correct vertex
      else   {
	clusTimeDiffHistTOFVSdeltaEtaWrongVertex_ -> Fill(  fabs(treeVars_->superClusterEta[sc2]-treeVars_->superClusterEta[sc1]) , 
							    (bcTime1.time-extraTravelTime(sc1,u,(*treeVars_)))  - (bcTime2.time-extraTravelTime(sc2,u,(*treeVars_)))
							    //(bcTime1.time-extraTravelTime(sc1,(*treeVars_)))  - (bcTime2.time-extraTravelTime(sc2,(*treeVars_)))
							    ); //double
	clusTimeDiffHistTOFwrongVertex_           -> Fill( (bcTime1.time-extraTravelTime(sc1,u,(*treeVars_)))  - (bcTime2.time-extraTravelTime(sc2,u,(*treeVars_)) )); //double
      }// if wrong vertex
    }//loop on vertices

  } // if vertex matching succeeded
  //  else std::cout << "vertex was not found which matches electrons track... " << std::endl;
  
  // take care of the seeds
  seedTimeDiffHist_    -> Fill( bcTime1.seedtime - bcTime2.seedtime ); // double
  seedTimeDiffHistTOF_ -> Fill( (bcTime1.seedtime-extraTravelTime(sc1,(*treeVars_)))  - (bcTime2.seedtime-extraTravelTime(sc2,(*treeVars_)))  ); // double
  

  if(bcTime1.second>-1 && bcTime2.second>-1){  // check that there's crystals beyond seed // ALL doubles 
    secondTime_                     -> Fill( bcTime1.secondtime);  secondTime_->Fill(bcTime2.secondtime); 
    secondTimeDiffHist_             -> Fill( bcTime1.secondtime - bcTime2.secondtime );
    secondTimeDiffHistTOF_          -> Fill( (bcTime1.secondtime-extraTravelTime(sc1,(*treeVars_)))  - (bcTime2.secondtime-extraTravelTime(sc2,(*treeVars_)))  );
    seed2secondTimeDiffHist_        -> Fill( bcTime1.seedtime - bcTime2.secondtime );
    seed2secondTimeDiffHist_        -> Fill( bcTime2.seedtime - bcTime1.secondtime );
    seed2secondTimeDiffHistTOF_     -> Fill( (bcTime1.seedtime-extraTravelTime(sc1,(*treeVars_)))  - (bcTime2.secondtime-extraTravelTime(sc2,(*treeVars_)))  );
    seed2secondTimeDiffHistTOF_     -> Fill( (bcTime2.seedtime-extraTravelTime(sc2,(*treeVars_)))  - (bcTime1.secondtime-extraTravelTime(sc1,(*treeVars_)))  );
  }
  
  
  clusTimeDiffHist_    -> Fill(bcTime1.time - bcTime2.time ); // double
  clusTimeDiffHistTOF_ -> Fill( (bcTime1.time - extraTravelTime(sc1,(*treeVars_)) ) - (bcTime2.time -extraTravelTime(sc2,(*treeVars_))) );

  timeVsEtaLead_ -> Fill( treeVars_->superClusterEta[sc1] , (bcTime1.time - extraTravelTime(sc1,(*treeVars_)) ) ); // double 
  timeVsEtaSub_  -> Fill( treeVars_->superClusterEta[sc2] , (bcTime2.time - extraTravelTime(sc2,(*treeVars_)) ) ); // double

  return 1;
}
// end HistSet::fill


int HistSet::fillSingle(int sc1, int bc1, ClusterTime bcTime1){  
  // if launched with no arguments for cuts, 
  // call the generalized 'fil' setting type=0 and cut=0,
  // which means cuts are ignored 
  int type  =0;
  float cut =0;
  std::pair<float,float> dummyPhase = std::make_pair(0.,0.);
  return fillSingle(sc1, bc1, bcTime1,  type,  cut , dummyPhase);
}

int HistSet::fillSingle(int sc1, int bc1, ClusterTime bcTime1, std::pair<float,float>& globalPhases){  
  // if launched with no arguments for cuts, 
  // call the generalized 'fil' setting type=0 and cut=0,
  // which means cuts are ignored 
  int type  =0;
  float cut =0;
  timeCorrector theCorrector;
  return fillSingle(sc1, bc1, bcTime1,  type,  cut , globalPhases, theCorrector);
}


int HistSet::fillSingle(int sc1, int bc1, ClusterTime bcTime1, std::pair<float,float>& globalPhases, timeCorrector& theCorrector){  
  // if launched with no arguments for cuts, 
  // call the generalized 'fil' setting type=0 and cut=0,
  // which means cuts are ignored 
  int type  =0;
  float cut =0;
  return fillSingle(sc1, bc1, bcTime1,  type,  cut , globalPhases, theCorrector);
}


int HistSet::fillSingle(int sc1, int bc1, ClusterTime bcTime1, int type, float cut){
  std::pair<float,float> dummyPhase = std::make_pair(0.,0.);
  return fillSingle(sc1, bc1, bcTime1, type, cut , dummyPhase );
}


int HistSet::fillSingle(int sc1, int bc1, ClusterTime bcTime1, int type, float cut , std::pair<float,float>& globalPhases){  
  timeCorrector theCorrector;
  return fillSingle(sc1, bc1, bcTime1, type, cut , globalPhases, theCorrector);
}


int HistSet::fillSingle(int sc1, int bc1, ClusterTime bcTime1, int type, float cut , std::pair<float,float>& globalPhases, timeCorrector& theCorrector){  
  

  if(!bcTime1.isvalid) return -1;
  
  
  // check that treeVars_ actually holds data
  if(!treeVars_)   {
    std::cout << " treeVars_ is null when calling HistSet::fill! Bailing out."<< std::endl;
    assert(-1);
  }
  else             {if (0) {std::cout << " treeVars_ is not null:  " << treeVars_ << std::endl;}   }
  
  
  ///////////////////////////////////// cuts //////////////////////////////////////////////////////////////
  if (type==0){ // no cuts are being used
    ;
  }
  
  else if (type==1){ // cut straight on the agreement between .seedtime and .secondtime within either cluster
    if(  fabs( bcTime1.seedtime - bcTime1.secondtime ) > cut )        return 0;   }
  
  
  else if (type==2){ // cut on the agreement between .seedtime and .otherstime within either cluster
    if(  fabs( bcTime1.seedtime - bcTime1.otherstime ) > cut )        return 0;   }
  
  
  else if (type==3){ // cut on chi2 of either cluster
    if(  (bcTime1.numCry<2) )       return 0;
    if(  fabs( bcTime1.chi2/bcTime1.numCry ) > cut )       return 0;
  }
  
  else{
    std::cout << "+++ You're using HistSet with a cut type: " << type << " which is not implemented; bailing out" << std::endl;
    assert(-1);
  }
  ///////////////////////////////////// cuts are done////////////////////////////////////////////////////////  


  TOFsingle_ -> Fill( extraTravelTime(sc1,(*treeVars_)) ); //single

  if(bcTime1.numCry>1){
    chi2_   ->Fill(bcTime1.chi2);      chi2ndf_->Fill(bcTime1.chi2/bcTime1.numCry);   } //single

  // take care of the seeds
  seedTime_               -> Fill(bcTime1.seedtime); //single
  seedTimeTOF_            -> Fill(bcTime1.seedtime-extraTravelTime(sc1,(*treeVars_))); //single


  // take care of the second-highest amplitude crystal
  if(bcTime1.second>-1) {
    secondAmpli_        -> Fill(treeVars_->xtalInBCEnergy[bc1][bcTime1.second]);  // check that there's crystals beyond seed // single
    seed2secSingleClus_ -> Fill( bcTime1.seedtime - bcTime1.secondtime ); // single
    secondTime_         -> Fill( bcTime1.secondtime);
    secondTimeTOF_      -> Fill( bcTime1.secondtime-extraTravelTime(sc1,(*treeVars_)));
  }
  else                  secondAmpli_->Fill(0);  
  
  clusterTime_         -> Fill(bcTime1.time);
  clusterTimeTOF_      -> Fill(bcTime1.time-extraTravelTime(sc1,(*treeVars_)));
  seedAmpli_->Fill(treeVars_->xtalInBCEnergy[bc1][bcTime1.seed]); // single

 
  if(bcTime1.otherstime>-999) // check that there's crystals beyond seed
    {
      diffSeedOther_ -> Fill(bcTime1.seedtime-bcTime1.otherstime); // single
      diffSeedOtherOverErr_->Fill( (bcTime1.seedtime-bcTime1.otherstime) / sqrt( pow(treeVars_->xtalInBCTimeErr[bc1][bcTime1.seed],2) -0.6*0.6+timingResParamConstEB*timingResParamConstEB + pow(bcTime1.otherstimeErr,2)) ); // single
    }
  if( bcTime1.second>-1 
      && treeVars_->xtalInBCTimeErr[bc1][bcTime1.second]>-10 
      && treeVars_->xtalInBCTimeErr[bc1][bcTime1.seed]>-1  )
    {
      diffSeedSecond_ -> Fill(bcTime1.seedtime-treeVars_->xtalInBCTime[bc1][bcTime1.second]); // single
      seedVSSecond_ -> Fill(treeVars_->xtalInBCTime[bc1][bcTime1.second],bcTime1.seedtime);  // single
      diffSeedSecondOverErr_    -> Fill( (bcTime1.seedtime-treeVars_->xtalInBCTime[bc1][bcTime1.second]) 
					 / sqrt( pow(treeVars_->xtalInBCTimeErr[bc1][bcTime1.seed],2) 
						 +  pow(treeVars_->xtalInBCTimeErr[bc1][bcTime1.second],2)
						 - 2* 0.6*0.6 + 2*timingResParamConstEB*timingResParamConstEB 
						 )   // single
					 ); 
  }
  
  
  numCryBC1->Fill(bcTime1.numCry); // single
  
  timeVsEta_ -> Fill( treeVars_->superClusterEta[sc1] , (bcTime1.time - extraTravelTime(sc1,(*treeVars_)) ) ); // single
  
  // catch location of time outliers
  if ( fabs(    (bcTime1.seedtime - extraTravelTime(sc1,(*treeVars_)))    )>1.5 )  outliersVsEtaPhi_ -> Fill( treeVars_->superClusterEta[sc1] , treeVars_->superClusterPhi[sc1]); // single
  
  return 1;
}
