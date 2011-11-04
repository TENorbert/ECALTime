#include "ECALTime/EcalTimePi0/interface/EcalObjectTime.h"
#include <iostream>
#include <math.h> 
#include <assert.h>


// ---------------------------------------------------------------------------------------
// --------------- Function to compute time and error for a physics object ---------------

ClusterTime timeAndUncertyPhoton(int bClusterIndex, EcalTimeTreeContent treeVars_)
{
  return timeAndUncertSingleCluster( bClusterIndex, treeVars_);
}

ClusterTime timeAndUncertyJet(int bClusterIndex, EcalTimeTreeContent treeVars_)
{
  return timeAndUncertSingleCluster( bClusterIndex, treeVars_);
}



// ---------------------------------------------------------------------------------------
// ------------------ Function to compute time and error for a cluster -------------------

ClusterTime timeAndUncertSingleCluster(int bClusterIndex, EcalTimeTreeContent treeVars_)
{
  ClusterTime theResult; //initialize
  theResult.isvalid = false;
  theResult.numCry = -999999;   theResult.time   = -999999;
  theResult.timeErr= -999999;   theResult.chi2   = -999999;
  theResult.seed   = -999999;   theResult.seedtime=-999999;
  theResult.second = -999999;   theResult.otherstime=-999999;
  theResult.otherstimeErr=-999999;

  float weightTsum  = 0;
  float weightSum   = 0;
  float weightTOthersum  = 0;
  float weightOtherSum   = 0;
  int   numCrystals = 0;
  float timingResParamN    =0;
  float timingResParamConst=0;
  

  
  bool  thisIsInEB=false;
  float sigmaNoiseOfThis=0;
  if(treeVars_.xtalInBCIEta[bClusterIndex][0]!=-999999)    sigmaNoiseOfThis   =sigmaNoiseEB;
  else                                                     sigmaNoiseOfThis   =sigmaNoiseEE;

  int seed(-1); float tmpEne=-9999; // cluster seed
  for (int cry=0; cry<treeVars_.nXtalsInCluster[bClusterIndex]; cry++){
    if(treeVars_.xtalInBCEnergy[bClusterIndex][cry]>tmpEne 
       && treeVars_.xtalInBCAmplitudeADC[bClusterIndex][cry]/sigmaNoiseOfThis>minAmpliOverSigma_
       ){
      tmpEne=treeVars_.xtalInBCEnergy[bClusterIndex][cry];
      seed=cry;
    } 	}
  int second(-1); tmpEne=-9999;   // second most energetic crystal
  for (int cry=0; cry<treeVars_.nXtalsInCluster[bClusterIndex]; cry++){
    if(treeVars_.xtalInBCEnergy[bClusterIndex][cry]>tmpEne 
       && treeVars_.xtalInBCAmplitudeADC[bClusterIndex][cry]/sigmaNoiseOfThis>minAmpliOverSigma_
       && cry!=seed )
      {
	tmpEne=treeVars_.xtalInBCEnergy[bClusterIndex][cry];
	second=cry;
      } 	}
  if(second==-1 && 0) std::cout << "second not found" << std::endl;

  
  if(0) std::cout << "\n++ BC statrs (eta: " << treeVars_.clusterEta[bClusterIndex] << ") : "  << std::endl;

  // loop on the cry components of a basic cluster; get timeBest and uncertainty 
  for(int thisCry=0; thisCry<treeVars_.nXtalsInCluster[bClusterIndex]; thisCry++)
    {
    if(treeVars_.xtalInBCIEta[bClusterIndex][thisCry]!=-999999)       {
      sigmaNoiseOfThis   =sigmaNoiseEB;
      timingResParamN    =timingResParamNEB;
      timingResParamConst=timingResParamConstEB;
      thisIsInEB=true;    }
    else if(treeVars_.xtalInBCIy[bClusterIndex][thisCry]!=-999999)    {
      sigmaNoiseOfThis=sigmaNoiseEE;
      timingResParamN    =timingResParamNEE;
      timingResParamConst=timingResParamConstEE;
      thisIsInEB=false;    }
    else    {  std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}

    // remove hits beyond gain switch
    if ( treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] > 3950 )  continue;

    float ampliOverSigOfThis = treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] / sigmaNoiseOfThis; 
    // minimum amplitude and spike rejection (could be updated)
    if( ampliOverSigOfThis < minAmpliOverSigma_) continue;
    if( treeVars_.xtalInBCSwissCross[bClusterIndex][thisCry] > 0.95) continue;

    numCrystals++;
    float timeOfThis  = treeVars_.xtalInBCTime[bClusterIndex][thisCry];
    //    old estimated: fully parameterized
    //    float sigmaOfThis = sqrt(pow(timingResParamN/ampliOverSigOfThis,2)+pow(timingResParamConst,2));
    //    new estimate: error from ratio + constant term  
    // float sigmaOfThis = pow( treeVars_.xtalInBCTimeErr[bClusterIndex][thisCry], 2 ) + pow( timingResParamConst, 2);
    // supposedly a large time constant is already included in the cxtalInBCTimeErr of 600 ps; rough estimate is that it's actually 300 (EB)

    // remove 0.6 constant term, put in timingResParamConstEB 
    float sigmaOfThis = pow( treeVars_.xtalInBCTimeErr[bClusterIndex][thisCry], 2) - 0.6*0.6 + timingResParamConstEB*timingResParamConstEB ;
    sigmaOfThis       = sqrt(sigmaOfThis);

    if(0) std::cout << "t: " << timeOfThis << " a/s: " << ampliOverSigOfThis << " sig: " << sigmaOfThis << "\t\t";
    weightTsum+=(timeOfThis/pow(sigmaOfThis,2));
    weightSum+=1/pow(sigmaOfThis,2);
    if(thisCry!=seed) {
      weightTOthersum+=(timeOfThis/pow(sigmaOfThis,2));
      weightOtherSum+=1/pow(sigmaOfThis,2);    }
    }
  float bestTime(-999999);
  if   (weightSum>0) bestTime=weightTsum/weightSum;
  else theResult.isvalid = false;
  float bestOtherTime(-999999);
  if   (weightOtherSum>0) bestOtherTime= weightTOthersum/weightOtherSum;
  // else std::cout << "bestOtherTime not made" << std::endl;
    

  float chi2 = -999999;
  // loop on the cry components to get chi2
  // do this only if you have at least 2 crystals over threshold and not spiky
  if(numCrystals>1){
    chi2=0;
    for(int thisCry=0; thisCry<treeVars_.nXtalsInCluster[bClusterIndex]; thisCry++)
      {
  	//bool  thisIsInEB=false;
  	float sigmaNoiseOfThis=0;
  	if(treeVars_.xtalInBCIEta[bClusterIndex][thisCry]!=-999999)       {
  	  sigmaNoiseOfThis=sigmaNoiseEB;
  	  //thisIsInEB=true;
  	}
  	else if(treeVars_.xtalInBCIy[bClusterIndex][thisCry]!=-999999)    {
  	  sigmaNoiseOfThis=sigmaNoiseEE;
  	  //thisIsInEB=false;    
  	}
  	else    {  std::cout << "crystal neither in eb nor in ee?? PROBLEM." << std::endl;}
  	
  	float ampliOverSigOfThis = treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] / sigmaNoiseOfThis; 
  	if( ampliOverSigOfThis < minAmpliOverSigma_) continue;
  	
	// remove hits beyond gain switch
	if ( treeVars_.xtalInBCAmplitudeADC[bClusterIndex][thisCry] > 3950 )  continue;
	if( treeVars_.xtalInBCSwissCross[bClusterIndex][thisCry] > 0.95) continue;

  	float timeOfThis  = treeVars_.xtalInBCTime[bClusterIndex][thisCry];
  	float sigmaOfThis = sqrt(pow(timingResParamN/ampliOverSigOfThis,2)+pow(timingResParamConst,2));
  	
  	chi2 += pow( (timeOfThis-bestTime)/sigmaOfThis, 2);
  	
      }// end loop on cry
  }//end if

  
  if(weightSum <= 0) {
    if(0) std::cout << "bestTime n.a. " << std::endl;
    theResult.isvalid = false;
    return theResult;}
  else{
    if(0) std::cout << "bestTime = " << bestTime << " error: " << sqrt(1/weightSum) << " chi2: " << chi2 << std::endl;//gfdebug
    theResult.isvalid    = true;
    theResult.numCry     = numCrystals;
    theResult.seed       = seed;
    theResult.second     = second;
    theResult.seedtime   = treeVars_.xtalInBCTime[bClusterIndex][seed];
    if(second>-1) {
      theResult.secondtime = treeVars_.xtalInBCTime[bClusterIndex][second];}
    theResult.time       = bestTime;
    theResult.timeErr    = sqrt(1/weightSum);
    theResult.otherstime = bestOtherTime;
    theResult.otherstimeErr=sqrt(1/weightOtherSum);
    theResult.chi2       = chi2;
    return theResult;
  }

}// end timeAndUncertSingleCluster



// ---------------------------------------------------------------------------------------
// ----------------  Functios to handle time of flight of supercluster -------------------

float travelDistance(int sc_num, EcalTimeTreeContent treeVars_) {
  return
    sqrt (	  pow( (treeVars_.superClusterVertexX[sc_num]-treeVars_.superClusterX[sc_num]), 2) +
		  pow( (treeVars_.superClusterVertexY[sc_num]-treeVars_.superClusterY[sc_num]), 2) +   
		  pow( (treeVars_.superClusterVertexZ[sc_num]-treeVars_.superClusterZ[sc_num]), 2)
		  );
}


float extraTravelTime(int sc_num, EcalTimeTreeContent & treeVars_) { // extra travel time with respect to collision at IP, in ns
  
  float travelled = sqrt (	  pow( (treeVars_.superClusterX[sc_num]-treeVars_.superClusterVertexX[sc_num]), 2) +
				  pow( (treeVars_.superClusterY[sc_num]-treeVars_.superClusterVertexY[sc_num]), 2) +   
				  pow( (treeVars_.superClusterZ[sc_num]-treeVars_.superClusterVertexZ[sc_num]), 2)
				  );
  float nominal = sqrt (	  pow( (treeVars_.superClusterX[sc_num]), 2) +
				  pow( (treeVars_.superClusterY[sc_num]), 2) +   
				  pow( (treeVars_.superClusterZ[sc_num]), 2)
				  );

  //std::cout << "extraTravelTime [ns]: " <<  (travelled-nominal)/100./lightSpeed*1e9 << std::endl;
  return  (travelled-nominal)/100./lightSpeed*1e9;

}


float extraTravelTime(int sc_num, int vtx_num, EcalTimeTreeContent & treeVars_) { // extra travel time with respect to an arbitrary vertex
  if(vtx_num<0 || vtx_num>=treeVars_.nVertices){
    std::cout<< "Usnig invalid vtx_num "<<vtx_num<<" within extraTravelTime(int sc_num, int vtx_num, EcalTimeTreeContent & treeVars_). Stopping the program." << std::endl;
    assert(0);
  }
  
  float travelled = sqrt (	  pow( (treeVars_.superClusterX[sc_num]-treeVars_.vtxX[vtx_num]), 2) +
				  pow( (treeVars_.superClusterY[sc_num]-treeVars_.vtxY[vtx_num]), 2) +   
				  pow( (treeVars_.superClusterZ[sc_num]-treeVars_.vtxZ[vtx_num]), 2)
				  );
  float nominal = sqrt (	  pow( (treeVars_.superClusterX[sc_num]), 2) +
				  pow( (treeVars_.superClusterY[sc_num]), 2) +   
				  pow( (treeVars_.superClusterZ[sc_num]), 2)
				  );

  return  (travelled-nominal)/100./lightSpeed*1e9;

}



/*
	// computing extra time of flight due to bending of electrons
	float r_curv1 = et1 / 0.3 / 3.8 ; // in meters
	float distShowerToVertex =        // IP - shower distance in transverse plane (~ ECAL Radius)
	pow( (treeVars_.superClusterVertexX[sc1]-treeVars_.superClusterX[sc1]), 2) +
	pow( (treeVars_.superClusterVertexY[sc1]-treeVars_.superClusterY[sc1]), 2);  
	// pow( (treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterZ[sc1]), 2);
	distShowerToVertex = sqrt(distShowerToVertex) / 100;  // in meters
	//std::cout << "distShowerToVertex: " << distShowerToVertex << std::endl;
	float alpha = asin( 0.5 * distShowerToVertex / r_curv1 ) ;
	float road  = r_curv1 * 2 * alpha;
	road        = sqrt( road*road + pow(  (treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterZ[sc1])  ,2)/100./100. );
	// actual time of flight of electron subtracted of straight line 
	float roadMoreThanPhoton = road - 
	sqrt( pow( (treeVars_.superClusterVertexX[sc1]-treeVars_.superClusterX[sc1]), 2) +
	pow( (treeVars_.superClusterVertexY[sc1]-treeVars_.superClusterY[sc1]), 2) +
	pow( (treeVars_.superClusterVertexZ[sc1]-treeVars_.superClusterZ[sc1]), 2)
	)/100.;
	std::cout << "eta: " << treeVars_.superClusterEta[sc1] << " pt: " << treeVars_.superClusterRawEnergy[sc1] 
	<< " r_curv1: " << r_curv1 << " transv-distShowerToVertex: " << distShowerToVertex 
	<< " path: " << road << " roadMoreThanPhoton: " << roadMoreThanPhoton << std::endl;
*/
