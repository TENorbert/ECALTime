#include "ECALTime/EcalTimePi0/interface/timeVsAmpliCorrector.h"

timeCorrector::timeCorrector(){

  // the  constructor sets the time corrections to 0 everywhere => the timeVsAmpliCorrector can be used in 'off mode' 
  // provided initEB and initEE are not invoked

  theCorrectionsEB.push_back(0.0);   theBinsEB.push_back(-10.0);
  theCorrectionsEB.push_back(0.0);   theBinsEB.push_back(21610.);

  theCorrectionsEE.push_back(0.0);   theBinsEE.push_back(-10.0);
  theCorrectionsEE.push_back(0.0);   theBinsEE.push_back(10626.2);

}

////////////////////////////////////////////////////////////////////////////////////////////
// implementation of the class methods

float timeCorrector::getCorrectionEB(float A){
  // pick a representative pseudorapidity to invoke the desired correction
  return getCorrection(A, 1.);
}


float timeCorrector::getCorrectionEE(float A){
  // pick a representative pseudorapidity to invoke the desired correction
  return getCorrection(A, 2.);
}


float timeCorrector::getCorrection(float A, float eta){

  std::vector<float> theCorrections;
  std::vector<float> theBins;

  // chose between EB and EE
  if( fabs(eta) < 1.45){
     theCorrections=theCorrectionsEB; 
     theBins=theBinsEB;
     if (theCorrections.size()==0) {std::cout << "you want corrections for EB, but the size of the correction vector is 0. Strange! " << std::endl;}
   }
  else if ( 1.45< fabs(eta) && fabs(eta) < 3.0){
     theCorrections=theCorrectionsEE; 
     theBins=theBinsEE;
     if (theCorrections.size()==0) {std::cout << "you want corrections for EE, but the size of the correction vector is 0. Strange! " << std::endl;}
  }
  else{
  std::cout << "wrong eta value: " << eta << std::endl; 
  return -999;
  }

  double theCorrection=0;

  int myBin = -1;
  for (int bin=0; bin<(int)theBins.size(); bin++ ){
    if(A > theBins.at(bin)) {
      myBin = bin;     }
    else break;
  }
  
  if (myBin == -1)
    {
      // if A below the covered range, return time correction for the lowest A bin; sign will be flipped at the end, since correction = -1*effect 
      theCorrection = theCorrections.at(0);
    }    
  else if  ( myBin == ((int)(theBins.size()-1))   ) 
    {
      // if A above the covered range, return time correction for the highest A bin sign; will be flipped at the end, since correction = -1*effect 
      theCorrection = theCorrections.at( myBin );      
    }    
  else if  ( -1 < myBin   &&   myBin <  ((int)theBins.size()-1) )
    {
      // interpolate linearly between two assingned points
      theCorrection  = ( theCorrections.at(myBin+1) - theCorrections.at(myBin) );
      theCorrection *= ( ((float)A) -  theBins.at(myBin) ) / ( theBins.at(myBin+1) - theBins.at(myBin) );
      theCorrection += theCorrections.at(myBin);
    }
  else
    {
      std::cout << "Assigning time correction impossible. Bailing out." << std::endl;
      assert(-1);
    }
  
  // debug
  // std::cout << "\nA: " << A << " eta: " << eta<< " nmyBin is: " << myBin << " and my corr is: " << theCorrection << "\n"<< std::endl; //GF 

  // flip the sing since: correction = -1*effect
  // return correction in ns
  return -1 * theCorrection;

}
  


void timeCorrector::initEE(std::string version){
   std::cout << "\ninitializing corrections for EE" << std::endl;
   theCorrectionsEE.clear();          theBinsEE.clear();
 
   // this is binned measured bias in nanosecond; correction needs be the opposited. 
   if(version==std::string("EElow")){

     // using Tambe's result: http://tambe.web.cern.ch/tambe/TimingResolutionPlots/PhysPhoton/CorrFactorsFromLowEEAndEB/

     //theCorrectionsEE.push_back(0.25);   theBinsEE.push_back(9.47062);
     //theCorrectionsEE.push_back(0.25);   theBinsEE.push_back(10.6262);
     //theCorrectionsEE.push_back(0.25);   theBinsEE.push_back(11.9228);
     //theCorrectionsEE.push_back(0.25);   theBinsEE.push_back(13.3776);
     //theCorrectionsEE.push_back(0.25);   theBinsEE.push_back(15.0099);
     //theCorrectionsEE.push_back(0.25);   theBinsEE.push_back(16.8414);
     //theCorrectionsEE.push_back(0.25);   theBinsEE.push_back(18.8964);
     //theCorrectionsEE.push_back(0.25);   theBinsEE.push_back(21.2021);
     //theCorrectionsEE.push_back(0.25);   theBinsEE.push_back(23.7891);
     theCorrectionsEE.push_back(-0.155735);   theBinsEE.push_back(26.6918);
     theCorrectionsEE.push_back(-0.109962);   theBinsEE.push_back(29.9487);
     theCorrectionsEE.push_back(-0.0715073);   theBinsEE.push_back(33.603);
     theCorrectionsEE.push_back(-0.060626);   theBinsEE.push_back(37.7032);
     theCorrectionsEE.push_back(-0.0463754);   theBinsEE.push_back(42.3037);
     theCorrectionsEE.push_back(-0.0414362);   theBinsEE.push_back(47.4656);
     theCorrectionsEE.push_back(-0.0339144);   theBinsEE.push_back(53.2572);
     theCorrectionsEE.push_back(-0.0138342);   theBinsEE.push_back(59.7556);
     theCorrectionsEE.push_back(-0.0199016);   theBinsEE.push_back(67.0469);
     theCorrectionsEE.push_back(-0.0199928);   theBinsEE.push_back(75.2278);
     theCorrectionsEE.push_back(-0.0193039);   theBinsEE.push_back(84.407);
     theCorrectionsEE.push_back(-0.0150271);   theBinsEE.push_back(94.7062);
     theCorrectionsEE.push_back(-0.00883894);   theBinsEE.push_back(106.262);
     theCorrectionsEE.push_back(-0.0160218);   theBinsEE.push_back(119.228);
     theCorrectionsEE.push_back(-0.0179918);   theBinsEE.push_back(133.776);
     theCorrectionsEE.push_back(-0.0203784);   theBinsEE.push_back(150.099);
     theCorrectionsEE.push_back(-0.0238918);   theBinsEE.push_back(168.414);
     theCorrectionsEE.push_back(-0.0172227);   theBinsEE.push_back(188.964);
     theCorrectionsEE.push_back(-0.0153935);   theBinsEE.push_back(212.021);
     theCorrectionsEE.push_back(-0.0166607);   theBinsEE.push_back(237.891);
     theCorrectionsEE.push_back(-0.0168141);   theBinsEE.push_back(266.918);
     theCorrectionsEE.push_back(-0.00898462);   theBinsEE.push_back(299.487);
     theCorrectionsEE.push_back(-0.0128605);   theBinsEE.push_back(336.03);
     theCorrectionsEE.push_back(-0.0110128);   theBinsEE.push_back(377.032);
     theCorrectionsEE.push_back(-0.0148466);   theBinsEE.push_back(423.037);
     theCorrectionsEE.push_back(-0.0124556);   theBinsEE.push_back(474.656);
     theCorrectionsEE.push_back(-0.0180545);   theBinsEE.push_back(532.572);
     theCorrectionsEE.push_back(-0.00312731);   theBinsEE.push_back(597.556);
     theCorrectionsEE.push_back(-0.00553613);   theBinsEE.push_back(670.469);
     theCorrectionsEE.push_back(-0.00063198);   theBinsEE.push_back(752.278);
     theCorrectionsEE.push_back(-0.017085);   theBinsEE.push_back(844.07);
     theCorrectionsEE.push_back(-0.0231402);   theBinsEE.push_back(947.062);
     theCorrectionsEE.push_back(-0.00411875);   theBinsEE.push_back(1062.62);
     theCorrectionsEE.push_back(-0.00740504);   theBinsEE.push_back(1192.28);
     theCorrectionsEE.push_back(-0.0123253);   theBinsEE.push_back(1337.76);
     theCorrectionsEE.push_back(-0.00596664);   theBinsEE.push_back(1500.99);
     theCorrectionsEE.push_back(0.0134503);   theBinsEE.push_back(1684.14);
     theCorrectionsEE.push_back(0.0338138);   theBinsEE.push_back(1889.64);
     theCorrectionsEE.push_back(0.0485774);   theBinsEE.push_back(2120.21);
     theCorrectionsEE.push_back(0.0493507);   theBinsEE.push_back(2378.91);
     theCorrectionsEE.push_back(0.0326397);   theBinsEE.push_back(2669.18);
     theCorrectionsEE.push_back(0.0410606);   theBinsEE.push_back(2994.87);
     theCorrectionsEE.push_back(0.019095);   theBinsEE.push_back(3360.3);
     theCorrectionsEE.push_back(0.0183005);   theBinsEE.push_back(3770.32);
     theCorrectionsEE.push_back(0.070756);   theBinsEE.push_back(4230.37);
     theCorrectionsEE.push_back(0.0951582);   theBinsEE.push_back(4746.56);
     theCorrectionsEE.push_back(0.136976);   theBinsEE.push_back(5325.72);
     theCorrectionsEE.push_back(0.113812);   theBinsEE.push_back(5975.56);
     theCorrectionsEE.push_back(0.156085);   theBinsEE.push_back(6704.69);
     theCorrectionsEE.push_back(0.229406);   theBinsEE.push_back(7522.78);
     theCorrectionsEE.push_back(0.677131);   theBinsEE.push_back(8440.7);
     theCorrectionsEE.push_back(0.522906);   theBinsEE.push_back(9470.62);
     theCorrectionsEE.push_back(0.251158);   theBinsEE.push_back(10626.2);

     //theCorrectionsEE.push_back(-1.71373);   theBinsEE.push_back(11922.8);
     //theCorrectionsEE.push_back(-4.75);      theBinsEE.push_back(13377.6);
     //theCorrectionsEE.push_back(-4.75);      theBinsEE.push_back(15009.9);
     //theCorrectionsEE.push_back(0);   theBinsEE.push_back(0);

   } else if (version==std::string("EE")){

     // obsolete, Mon Dec  5 12:12:29 CET 2011
     //theCorrectionsEE.push_back(0);   theBinsEE.push_back(18.9123);
     //theCorrectionsEE.push_back(0);   theBinsEE.push_back(21.218);
     //theCorrectionsEE.push_back(0);   theBinsEE.push_back(23.805);
     //theCorrectionsEE.push_back(-0.167688);   theBinsEE.push_back(26.7077);
     //theCorrectionsEE.push_back(-0.136405);   theBinsEE.push_back(29.9646);
     //theCorrectionsEE.push_back(-0.10191);   theBinsEE.push_back(33.6189);
     //theCorrectionsEE.push_back(-0.0818226);   theBinsEE.push_back(37.7191);
     //theCorrectionsEE.push_back(-0.0596632);   theBinsEE.push_back(42.3196);
     //theCorrectionsEE.push_back(-0.0452422);   theBinsEE.push_back(47.4814);
     //theCorrectionsEE.push_back(-0.0204676);   theBinsEE.push_back(53.2731);
     //theCorrectionsEE.push_back(-0.00708579);   theBinsEE.push_back(59.7715);
     //theCorrectionsEE.push_back(-0.00502538);   theBinsEE.push_back(67.0628);
     //theCorrectionsEE.push_back(-0.00289768);   theBinsEE.push_back(75.2437);
     //theCorrectionsEE.push_back(-0.00188651);   theBinsEE.push_back(84.4229);
     //theCorrectionsEE.push_back(0.00172187);   theBinsEE.push_back(94.7221);
     //theCorrectionsEE.push_back(0.00422572);   theBinsEE.push_back(106.278);
     //theCorrectionsEE.push_back(0.00352601);   theBinsEE.push_back(119.244);
     //theCorrectionsEE.push_back(0.000234735);   theBinsEE.push_back(133.792);
     //theCorrectionsEE.push_back(-0.00497225);   theBinsEE.push_back(150.115);
     //theCorrectionsEE.push_back(-0.00255359);   theBinsEE.push_back(168.43);
     //theCorrectionsEE.push_back(0.000276509);   theBinsEE.push_back(188.98);
     //theCorrectionsEE.push_back(-0.000676717);   theBinsEE.push_back(212.037);
     //theCorrectionsEE.push_back(0.00375508);   theBinsEE.push_back(237.907);
     //theCorrectionsEE.push_back(0.00800281);   theBinsEE.push_back(266.934);
     //theCorrectionsEE.push_back(0.00939014);   theBinsEE.push_back(299.503);
     //theCorrectionsEE.push_back(0.00933337);   theBinsEE.push_back(336.046);
     //theCorrectionsEE.push_back(0.0127813);   theBinsEE.push_back(377.048);
     //theCorrectionsEE.push_back(0.0149195);   theBinsEE.push_back(423.053);
     //theCorrectionsEE.push_back(0.0169709);   theBinsEE.push_back(474.671);
     //theCorrectionsEE.push_back(0.0198988);   theBinsEE.push_back(532.588);
     //theCorrectionsEE.push_back(0.0265556);   theBinsEE.push_back(597.572);
     //theCorrectionsEE.push_back(0.0266586);   theBinsEE.push_back(670.485);
     //theCorrectionsEE.push_back(0.0308715);   theBinsEE.push_back(752.294);
     //theCorrectionsEE.push_back(0.0147255);   theBinsEE.push_back(844.086);
     //theCorrectionsEE.push_back(0.0155513);   theBinsEE.push_back(947.078);
     //theCorrectionsEE.push_back(0.0234095);   theBinsEE.push_back(1062.64);
     //theCorrectionsEE.push_back(0.0319073);   theBinsEE.push_back(1192.3);
     //theCorrectionsEE.push_back(0.0440146);   theBinsEE.push_back(1337.78);
     //theCorrectionsEE.push_back(0.0401477);   theBinsEE.push_back(1501.01);
     //theCorrectionsEE.push_back(0.0538419);   theBinsEE.push_back(1684.16);
     //theCorrectionsEE.push_back(0.06414);   theBinsEE.push_back(1889.65);
     //theCorrectionsEE.push_back(0.0769302);   theBinsEE.push_back(2120.22);
     //theCorrectionsEE.push_back(0.0733546);   theBinsEE.push_back(2378.93);
     //theCorrectionsEE.push_back(0.0753102);   theBinsEE.push_back(2669.2);
     //theCorrectionsEE.push_back(0.0786608);   theBinsEE.push_back(2994.89);
     //theCorrectionsEE.push_back(0.0960978);   theBinsEE.push_back(3360.32);
     //theCorrectionsEE.push_back(0.111197);   theBinsEE.push_back(3770.34);
     //theCorrectionsEE.push_back(0.108132);   theBinsEE.push_back(4230.39);
     //theCorrectionsEE.push_back(0.12744);   theBinsEE.push_back(4746.57);
     //theCorrectionsEE.push_back(0.177114);   theBinsEE.push_back(5325.74);
     //theCorrectionsEE.push_back(0.168588);   theBinsEE.push_back(5975.58);
     //theCorrectionsEE.push_back(0.178248);   theBinsEE.push_back(6704.7);
     //theCorrectionsEE.push_back(0.286695);   theBinsEE.push_back(7522.8);
     //theCorrectionsEE.push_back(0.649502);   theBinsEE.push_back(8440.72);
     //theCorrectionsEE.push_back(0.561537);   theBinsEE.push_back(9470.64);
     //theCorrectionsEE.push_back(0.315861);   theBinsEE.push_back(10626.2);
     //theCorrectionsEE.push_back(-0.260055);   theBinsEE.push_back(11922.8);
     //theCorrectionsEE.push_back(-0.485963);   theBinsEE.push_back(13377.6);
     //theCorrectionsEE.push_back(-0.526979);   theBinsEE.push_back(15009.9);
     //theCorrectionsEE.push_back(-0.320099);   theBinsEE.push_back(16841.4);
     std::cout << "you've selected a correction type (" << version << ") for EE which is OBSOLETE and won't be used." << std::endl;
   }
   else {
   std::cout << "you've selected a correction type (" << version << ") for EE which does not exist" << std::endl;
   }

  if(theCorrectionsEE.size() != theBinsEE.size()){
  std::cout << "theCorrectionsEE and theBinsEE don't have the same size; bailing out" << std::endl;
  assert(0);
  } else{
  std::cout << "number of bins for EE corrections: " << theBinsEE.size() << "\n"<< std::endl;
  }

}// initEE




void timeCorrector::initEB(std::string version){
  std::cout << "\ninitializing corrections for EB" << std::endl;

   theCorrectionsEB.clear();          theBinsEB.clear();

   // this is binned measured bias in nanosecond; correction needs be the opposited. 
   if(version==std::string("EBmod4")){

   theCorrectionsEB.push_back(0.0427097);   theBinsEB.push_back(30.5506);
   theCorrectionsEB.push_back(0.0324853);   theBinsEB.push_back(34.2752);
   theCorrectionsEB.push_back(0.0325199);   theBinsEB.push_back(38.4542);
   theCorrectionsEB.push_back(0.0332106);   theBinsEB.push_back(43.1432);
   theCorrectionsEB.push_back(0.034634);   theBinsEB.push_back(48.4044);
   theCorrectionsEB.push_back(0.0335669);   theBinsEB.push_back(54.3075);
   theCorrectionsEB.push_back(0.0295928);   theBinsEB.push_back(60.9309);
   theCorrectionsEB.push_back(0.0277353);   theBinsEB.push_back(68.3624);
   theCorrectionsEB.push_back(0.0221848);   theBinsEB.push_back(76.7008);
   theCorrectionsEB.push_back(0.0165266);   theBinsEB.push_back(86.0566);
   theCorrectionsEB.push_back(0.0122511);   theBinsEB.push_back(96.5539);
   theCorrectionsEB.push_back(0.0102246);   theBinsEB.push_back(108.332);
   theCorrectionsEB.push_back(0.0119881);   theBinsEB.push_back(121.548);
   theCorrectionsEB.push_back(0.0106686);   theBinsEB.push_back(136.375);
   theCorrectionsEB.push_back(0.00547731);   theBinsEB.push_back(153.013);
   theCorrectionsEB.push_back(0.00245125);   theBinsEB.push_back(171.68);
   theCorrectionsEB.push_back(0.00516569);   theBinsEB.push_back(192.625);
   theCorrectionsEB.push_back(0.00210641);   theBinsEB.push_back(216.126);
   theCorrectionsEB.push_back(0.000294294);   theBinsEB.push_back(242.494);
   theCorrectionsEB.push_back(0.00208257);   theBinsEB.push_back(272.079);
   theCorrectionsEB.push_back(0.00211399);   theBinsEB.push_back(305.275);
   theCorrectionsEB.push_back(7.8848e-05);   theBinsEB.push_back(342.521);
   theCorrectionsEB.push_back(-0.00245371);   theBinsEB.push_back(384.312);
   theCorrectionsEB.push_back(-0.00114789);   theBinsEB.push_back(431.202);
   theCorrectionsEB.push_back(-0.00419166);   theBinsEB.push_back(483.813);
   theCorrectionsEB.push_back(-0.0050314);   theBinsEB.push_back(542.844);
   theCorrectionsEB.push_back(0.00228042);   theBinsEB.push_back(609.078);
   theCorrectionsEB.push_back(-0.00554273);   theBinsEB.push_back(683.394);
   theCorrectionsEB.push_back(0.00669728);   theBinsEB.push_back(766.777);
   theCorrectionsEB.push_back(-0.00950755);   theBinsEB.push_back(860.335);
   theCorrectionsEB.push_back(-0.0101945);   theBinsEB.push_back(965.308);
   theCorrectionsEB.push_back(-0.0111277);   theBinsEB.push_back(1083.09);
   theCorrectionsEB.push_back(-0.0138652);   theBinsEB.push_back(1215.24);
   theCorrectionsEB.push_back(-0.0129775);   theBinsEB.push_back(1363.52);
   theCorrectionsEB.push_back(-0.0243837);   theBinsEB.push_back(1529.9);
   theCorrectionsEB.push_back(-0.0340678);   theBinsEB.push_back(1716.57);
   theCorrectionsEB.push_back(-0.0222574);   theBinsEB.push_back(1926.02);
   theCorrectionsEB.push_back(-0.0188277);   theBinsEB.push_back(2161.02);
   theCorrectionsEB.push_back(-0.0401259);   theBinsEB.push_back(2424.71);
   theCorrectionsEB.push_back(-0.0512177);   theBinsEB.push_back(2720.56);
   theCorrectionsEB.push_back(-0.0397173);   theBinsEB.push_back(3052.52);
   theCorrectionsEB.push_back(-0.0307929);   theBinsEB.push_back(3424.98);
   theCorrectionsEB.push_back(0.0123128);   theBinsEB.push_back(3842.89);
   theCorrectionsEB.push_back(-0.121747);   theBinsEB.push_back(4311.79);
   theCorrectionsEB.push_back(-0.139355);   theBinsEB.push_back(4837.9);
   theCorrectionsEB.push_back(-0.199441);   theBinsEB.push_back(5428.21);
   theCorrectionsEB.push_back(-0.382732);   theBinsEB.push_back(6090.55);
   theCorrectionsEB.push_back(-0.535643);   theBinsEB.push_back(6833.7);
   theCorrectionsEB.push_back(-0.703777);   theBinsEB.push_back(7667.54);
   theCorrectionsEB.push_back(-0.520158);   theBinsEB.push_back(8603.12);
   theCorrectionsEB.push_back(-0.729732);   theBinsEB.push_back(9652.85);
   theCorrectionsEB.push_back(-0.757868);   theBinsEB.push_back(10830.7);
   theCorrectionsEB.push_back(-0.732588);   theBinsEB.push_back(12152.2);
   theCorrectionsEB.push_back(-1.03576);   theBinsEB.push_back(13635);
   theCorrectionsEB.push_back(4.49757);   theBinsEB.push_back(15298.7);
   theCorrectionsEB.push_back(-0.978876);   theBinsEB.push_back(17165.5);
   //theCorrectionsEB.push_back(0);   theBinsEB.push_back(19259.9);
   //theCorrectionsEB.push_back(0);   theBinsEB.push_back(21610);
   //theCorrectionsEB.push_back(0);   theBinsEB.push_back(24246.8);
   //theCorrectionsEB.push_back(0);   theBinsEB.push_back(27205.4);
   //theCorrectionsEB.push_back(0);   theBinsEB.push_back(0.00104358);

   } else if (version==std::string("EB") ){

     // using Tambe's result: http://tambe.web.cern.ch/tambe/TimingResolutionPlots/PhysPhoton/CorrFactorsFromLowEEAndEB/
     theCorrectionsEB.push_back(0.129119);   theBinsEB.push_back(15.2987);
     theCorrectionsEB.push_back(0.0666236);   theBinsEB.push_back(17.1654);
     theCorrectionsEB.push_back(0.0750206);   theBinsEB.push_back(19.2599);
     theCorrectionsEB.push_back(0.0585525);   theBinsEB.push_back(21.61);
     theCorrectionsEB.push_back(0.0524942);   theBinsEB.push_back(24.2468);
     theCorrectionsEB.push_back(0.0410405);   theBinsEB.push_back(27.2054);
     theCorrectionsEB.push_back(0.0422359);   theBinsEB.push_back(30.5249);
     theCorrectionsEB.push_back(0.0323442);   theBinsEB.push_back(34.2495);
     theCorrectionsEB.push_back(0.0261736);   theBinsEB.push_back(38.4286);
     theCorrectionsEB.push_back(0.0184196);   theBinsEB.push_back(43.1176);
     theCorrectionsEB.push_back(0.0184786);   theBinsEB.push_back(48.3787);
     theCorrectionsEB.push_back(0.0145926);   theBinsEB.push_back(54.2818);
     theCorrectionsEB.push_back(0.00614763);   theBinsEB.push_back(60.9052);
     theCorrectionsEB.push_back(0.00712762);   theBinsEB.push_back(68.3368);
     theCorrectionsEB.push_back(0.00556952);   theBinsEB.push_back(76.6751);
     theCorrectionsEB.push_back(0.000358535);   theBinsEB.push_back(86.0309);
     theCorrectionsEB.push_back(-0.00035637);   theBinsEB.push_back(96.5283);
     theCorrectionsEB.push_back(-0.00591454);   theBinsEB.push_back(108.307);
     theCorrectionsEB.push_back(-0.00270605);   theBinsEB.push_back(121.522);
     theCorrectionsEB.push_back(-0.00176438);   theBinsEB.push_back(136.35);
     theCorrectionsEB.push_back(-0.00849706);   theBinsEB.push_back(152.987);
     theCorrectionsEB.push_back(-0.0116634);   theBinsEB.push_back(171.654);
     theCorrectionsEB.push_back(-0.00833955);   theBinsEB.push_back(192.599);
     theCorrectionsEB.push_back(-0.0153301);   theBinsEB.push_back(216.1);
     theCorrectionsEB.push_back(-0.0151621);   theBinsEB.push_back(242.468);
     theCorrectionsEB.push_back(-0.00896533);   theBinsEB.push_back(272.054);
     theCorrectionsEB.push_back(-0.01008);   theBinsEB.push_back(305.249);
     theCorrectionsEB.push_back(-0.0188685);   theBinsEB.push_back(342.495);
     theCorrectionsEB.push_back(-0.0171632);   theBinsEB.push_back(384.286);
     theCorrectionsEB.push_back(-0.0159843);   theBinsEB.push_back(431.176);
     theCorrectionsEB.push_back(-0.0185847);   theBinsEB.push_back(483.787);
     theCorrectionsEB.push_back(-0.0209603);   theBinsEB.push_back(542.818);
     theCorrectionsEB.push_back(-0.0130973);   theBinsEB.push_back(609.052);
     theCorrectionsEB.push_back(-0.0154015);   theBinsEB.push_back(683.368);
     theCorrectionsEB.push_back(-0.00529718);   theBinsEB.push_back(766.751);
     theCorrectionsEB.push_back(-0.0174984);   theBinsEB.push_back(860.309);
     theCorrectionsEB.push_back(-0.0284475);   theBinsEB.push_back(965.283);
     theCorrectionsEB.push_back(-0.0264779);   theBinsEB.push_back(1083.07);
     theCorrectionsEB.push_back(-0.0234928);   theBinsEB.push_back(1215.22);
     theCorrectionsEB.push_back(-0.0262998);   theBinsEB.push_back(1363.5);
     theCorrectionsEB.push_back(-0.0412084);   theBinsEB.push_back(1529.87);
     theCorrectionsEB.push_back(-0.043147);   theBinsEB.push_back(1716.54);
     theCorrectionsEB.push_back(-0.034127);   theBinsEB.push_back(1925.99);
     theCorrectionsEB.push_back(-0.0282815);   theBinsEB.push_back(2161);
     theCorrectionsEB.push_back(-0.0649723);   theBinsEB.push_back(2424.68);
     theCorrectionsEB.push_back(-0.0586241);   theBinsEB.push_back(2720.54);
     theCorrectionsEB.push_back(-0.0689535);   theBinsEB.push_back(3052.49);
     theCorrectionsEB.push_back(-0.0644636);   theBinsEB.push_back(3424.95);
     theCorrectionsEB.push_back(-0.0302987);   theBinsEB.push_back(3842.86);
     theCorrectionsEB.push_back(-0.231401);   theBinsEB.push_back(4311.76);
     theCorrectionsEB.push_back(-0.271354);   theBinsEB.push_back(4837.87);
     theCorrectionsEB.push_back(-0.334614);   theBinsEB.push_back(5428.18);
     theCorrectionsEB.push_back(-0.500895);   theBinsEB.push_back(6090.52);
     theCorrectionsEB.push_back(-0.671807);   theBinsEB.push_back(6833.68);
     theCorrectionsEB.push_back(-0.71232);   theBinsEB.push_back(7667.51);
     theCorrectionsEB.push_back(-0.69643);   theBinsEB.push_back(8603.09);
     theCorrectionsEB.push_back(-0.692654);   theBinsEB.push_back(9652.83);
     theCorrectionsEB.push_back(-0.810051);   theBinsEB.push_back(10830.7);
     theCorrectionsEB.push_back(-0.821316);   theBinsEB.push_back(12152.2);
     theCorrectionsEB.push_back(-0.786343);   theBinsEB.push_back(13635);
     // the points below are outside of the plot: 
     // http://tambe.web.cern.ch/tambe/TimingResolutionPlots/PhysPhoton/CorrFactorsFromLowEEAndEB/Mean_Time_(+250ps)_Vs_Amplitude_Crys%20In%20EB.png
     // check before confirming
     theCorrectionsEB.push_back(-1.28217);   theBinsEB.push_back(15298.7);
     theCorrectionsEB.push_back(-2.25185);   theBinsEB.push_back(17165.4);
     theCorrectionsEB.push_back(-3.48219);   theBinsEB.push_back(19259.9);
     theCorrectionsEB.push_back(-4.75);   theBinsEB.push_back(21610);
     // theCorrectionsEB.push_back(5.25);   theBinsEB.push_back(24246.8); // this looks unphysical
     
   }
   else {
   std::cout << "you've selected a correction type  (" << version << ")  for EB which does not exist" << std::endl;
   }  


  if(theCorrectionsEB.size() != theBinsEB.size()){
  std::cout << "theCorrectionsEB and theBinsEB don't have the same size; bailing out" << std::endl;
  assert(0);
  } else{
  std::cout << "number of bins for EB corrections: " << theBinsEB.size() << "\n" << std::endl;
  }

}//initEB
