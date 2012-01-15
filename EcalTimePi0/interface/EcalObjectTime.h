//#include "CalibCalorimetry/EcalTiming/interface/EcalTimeTreeContent.h"
#include "ECALTime/EcalTimePi0/interface/EcalTimePhyTreeContent.h"

#define lightSpeed 299792458

// Consts
const float sigmaNoiseEB          = 1.06;  // ADC ; using total single-sample noise
const float sigmaNoiseEE          = 2.10;  // ADC ; using total single-sample noise
const float timingResParamNEB     = 28.51; // ns ; plots approved http://indico.cern.ch/conferenceDisplay.py?confId=92739
const float timingResParamConstEB = 0.1;   // ns ; 
const float timingResParamNEE     = 31.84; // ns ; Fig. 2 from CFT-09-006
const float timingResParamConstEE = 0.3;   // ns ; rough, probably conservative estimate
const float minAmpliOverSigma_    = 30;    // dimensionless, and holds consistently for EB and EE
const float maxChi2NDF_           = 20;  //TODO: gf configurable

// container to hold time-relaed observables for clusters or objects
struct ClusterTime {
  bool  isvalid;
  int   numCry;
  int   seed;
  int   second;
  float seedtime;
  float secondtime;
  float time;
  float timeErr;
  float otherstime;
  float otherstimeErr;
  float chi2;
} ;

ClusterTime timeAndUncertSingleCluster (int bClusterIndex, EcalTimePhyTreeContent treeVars_);
ClusterTime timeAndUncertyPhoton       (int bClusterIndex, EcalTimePhyTreeContent treeVars_);
ClusterTime timeAndUncertyJet          (int bClusterIndex, EcalTimePhyTreeContent treeVars_);


float travelDistance(int sc_num, EcalTimePhyTreeContent treeVars_);
float extraTravelTime(int sc_num, EcalTimePhyTreeContent & treeVars_);
float extraTravelTime(int sc_num, int vtx_num, EcalTimePhyTreeContent & treeVars_);
