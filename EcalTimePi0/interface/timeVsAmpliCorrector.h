#include <assert.h>
#include <stdio.h>

#include <iostream>
#include<string>
#include<vector>


#include <math.h>

////////////////////////////////////////////////////////////////////////////////////////////
// holds A-dependent corrections and provides interface
class timeCorrector {
  private:
    std::vector <float> theCorrectionsEB; 
    std::vector <float> theBinsEB;   

    std::vector <float> theCorrectionsEE; 
    std::vector <float> theBinsEE; 


  public:
    float getCorrection(float,float);
    float getCorrectionEB(float);
    float getCorrectionEE(float);
    void initEB(std::string);
    void initEE(std::string);
  };



////////////////////////////////////////////////////////////////////////////////////////////
// implementation of the class methods

float timeCorrector::getCorrectionEB(float A){
  return timeCorrector::getCorrection(A, 1.);
}


float timeCorrector::getCorrectionEE(float A){
  return timeCorrector::getCorrection(A, 2.);
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
      // if A below the covered range, return time correction for the lowest A bin
      theCorrection = theCorrections.at(0);
    }    
  else if  ( myBin == ((int)(theBins.size()-1))   ) 
    {
      // if A above the covered range, return time correction for the highest A bin
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

  // return correction in ns
  return theCorrection;

}
  


void timeCorrector::initEE(std::string version){
   std::cout << "\ninitializing corrections for EE" << std::endl;
   theCorrectionsEE.clear();          theBinsEE.clear();
 
   // these corrections are in nanosecond
   if(version==std::string("EElow")){

    //theCorrectionsEE.push_back(0);   theBinsEE.push_back(18.9123);
    //theCorrectionsEE.push_back(0);   theBinsEE.push_back(21.218);
    //theCorrectionsEE.push_back(0);   theBinsEE.push_back(23.805);
    theCorrectionsEE.push_back(-0.133039);   theBinsEE.push_back(26.7077);
    theCorrectionsEE.push_back(-0.103766);   theBinsEE.push_back(29.9646);
    theCorrectionsEE.push_back(-0.0728432);   theBinsEE.push_back(33.6189);
    theCorrectionsEE.push_back(-0.0578023);   theBinsEE.push_back(37.7191);
    theCorrectionsEE.push_back(-0.0422981);   theBinsEE.push_back(42.3196);
    theCorrectionsEE.push_back(-0.0345444);   theBinsEE.push_back(47.4814);
    theCorrectionsEE.push_back(-0.0129457);   theBinsEE.push_back(53.2731);
    theCorrectionsEE.push_back(-0.00379589);   theBinsEE.push_back(59.7715);
    theCorrectionsEE.push_back(-0.00313104);   theBinsEE.push_back(67.0628);
    theCorrectionsEE.push_back(-0.000504792);   theBinsEE.push_back(75.2437);
    theCorrectionsEE.push_back(0.000467353);   theBinsEE.push_back(84.4229);
    theCorrectionsEE.push_back(0.00421467);   theBinsEE.push_back(94.7221);
    theCorrectionsEE.push_back(0.00638916);   theBinsEE.push_back(106.278);
    theCorrectionsEE.push_back(0.00640062);   theBinsEE.push_back(119.244);
    theCorrectionsEE.push_back(0.00341919);   theBinsEE.push_back(133.792);
    theCorrectionsEE.push_back(-0.00301614);   theBinsEE.push_back(150.115);
    theCorrectionsEE.push_back(-0.00250331);   theBinsEE.push_back(168.43);
    theCorrectionsEE.push_back(-0.000503837);   theBinsEE.push_back(188.98);
    theCorrectionsEE.push_back(-0.00158694);   theBinsEE.push_back(212.037);
    theCorrectionsEE.push_back(0.00106072);   theBinsEE.push_back(237.907);
    theCorrectionsEE.push_back(0.00538018);   theBinsEE.push_back(266.934);
    theCorrectionsEE.push_back(0.00640622);   theBinsEE.push_back(299.503);
    theCorrectionsEE.push_back(0.00550587);   theBinsEE.push_back(336.046);
    theCorrectionsEE.push_back(0.00805342);   theBinsEE.push_back(377.048);
    theCorrectionsEE.push_back(0.00965012);   theBinsEE.push_back(423.053);
    theCorrectionsEE.push_back(0.0111103);   theBinsEE.push_back(474.671);
    theCorrectionsEE.push_back(0.0138681);   theBinsEE.push_back(532.588);
    theCorrectionsEE.push_back(0.020215);   theBinsEE.push_back(597.572);
    theCorrectionsEE.push_back(0.0194182);   theBinsEE.push_back(670.485);
    theCorrectionsEE.push_back(0.0249657);   theBinsEE.push_back(752.294);
    theCorrectionsEE.push_back(0.00798197);   theBinsEE.push_back(844.086);
    theCorrectionsEE.push_back(0.00655069);   theBinsEE.push_back(947.078);
    theCorrectionsEE.push_back(0.0131119);   theBinsEE.push_back(1062.64);
    theCorrectionsEE.push_back(0.0210328);   theBinsEE.push_back(1192.3);
    theCorrectionsEE.push_back(0.0336762);   theBinsEE.push_back(1337.78);
    theCorrectionsEE.push_back(0.0304386);   theBinsEE.push_back(1501.01);
    theCorrectionsEE.push_back(0.0417921);   theBinsEE.push_back(1684.16);
    theCorrectionsEE.push_back(0.0525334);   theBinsEE.push_back(1889.65);
    theCorrectionsEE.push_back(0.0647914);   theBinsEE.push_back(2120.22);
    theCorrectionsEE.push_back(0.0605786);   theBinsEE.push_back(2378.93);
    theCorrectionsEE.push_back(0.0594094);   theBinsEE.push_back(2669.2);
    theCorrectionsEE.push_back(0.0678715);   theBinsEE.push_back(2994.89);
    theCorrectionsEE.push_back(0.0787326);   theBinsEE.push_back(3360.32);
    theCorrectionsEE.push_back(0.0905882);   theBinsEE.push_back(3770.34);
    theCorrectionsEE.push_back(0.0856276);   theBinsEE.push_back(4230.39);
    theCorrectionsEE.push_back(0.0967914);   theBinsEE.push_back(4746.57);
    theCorrectionsEE.push_back(0.163573);   theBinsEE.push_back(5325.74);
    theCorrectionsEE.push_back(0.144924);   theBinsEE.push_back(5975.58);
    theCorrectionsEE.push_back(0.147008);   theBinsEE.push_back(6704.7);
    theCorrectionsEE.push_back(0.268835);   theBinsEE.push_back(7522.8);
    theCorrectionsEE.push_back(0.605681);   theBinsEE.push_back(8440.72);
    theCorrectionsEE.push_back(0.52275);   theBinsEE.push_back(9470.64);
    theCorrectionsEE.push_back(0.404734);   theBinsEE.push_back(10626.2);
    theCorrectionsEE.push_back(-0.155081);   theBinsEE.push_back(11922.8); // GF which of the following is reliable?
    theCorrectionsEE.push_back(0.138734);   theBinsEE.push_back(13377.6);
    theCorrectionsEE.push_back(-0.232964);   theBinsEE.push_back(15009.9);
    theCorrectionsEE.push_back(-0.0243007);   theBinsEE.push_back(16841.4);
    //theCorrectionsEE.push_back(0);   theBinsEE.push_back(0.000581012);

   } else if (version==std::string("EE")){

    //theCorrectionsEE.push_back(0);   theBinsEE.push_back(30.5506);
    //theCorrectionsEE.push_back(0);   theBinsEE.push_back(34.2752);
    //theCorrectionsEE.push_back(0);   theBinsEE.push_back(38.4542);
    theCorrectionsEE.push_back(-0.167688);   theBinsEE.push_back(43.1432);
    theCorrectionsEE.push_back(-0.136405);   theBinsEE.push_back(48.4044);
    theCorrectionsEE.push_back(-0.10191);   theBinsEE.push_back(54.3075);
    theCorrectionsEE.push_back(-0.0818226);   theBinsEE.push_back(60.9309);
    theCorrectionsEE.push_back(-0.0596632);   theBinsEE.push_back(68.3624);
    theCorrectionsEE.push_back(-0.0452422);   theBinsEE.push_back(76.7008);
    theCorrectionsEE.push_back(-0.0204676);   theBinsEE.push_back(86.0566);
    theCorrectionsEE.push_back(-0.00708579);   theBinsEE.push_back(96.5539);
    theCorrectionsEE.push_back(-0.00502538);   theBinsEE.push_back(108.332);
    theCorrectionsEE.push_back(-0.00289768);   theBinsEE.push_back(121.548);
    theCorrectionsEE.push_back(-0.00188651);   theBinsEE.push_back(136.375);
    theCorrectionsEE.push_back(0.00172187);   theBinsEE.push_back(153.013);
    theCorrectionsEE.push_back(0.00422572);   theBinsEE.push_back(171.68);
    theCorrectionsEE.push_back(0.00352601);   theBinsEE.push_back(192.625);
    theCorrectionsEE.push_back(0.000234735);   theBinsEE.push_back(216.126);
    theCorrectionsEE.push_back(-0.00497225);   theBinsEE.push_back(242.494);
    theCorrectionsEE.push_back(-0.00255359);   theBinsEE.push_back(272.079);
    theCorrectionsEE.push_back(0.000276509);   theBinsEE.push_back(305.275);
    theCorrectionsEE.push_back(-0.000676717);   theBinsEE.push_back(342.521);
    theCorrectionsEE.push_back(0.00375508);   theBinsEE.push_back(384.312);
    theCorrectionsEE.push_back(0.00800281);   theBinsEE.push_back(431.202);
    theCorrectionsEE.push_back(0.00939014);   theBinsEE.push_back(483.813);
    theCorrectionsEE.push_back(0.00933337);   theBinsEE.push_back(542.844);
    theCorrectionsEE.push_back(0.0127813);   theBinsEE.push_back(609.078);
    theCorrectionsEE.push_back(0.0149195);   theBinsEE.push_back(683.394);
    theCorrectionsEE.push_back(0.0169709);   theBinsEE.push_back(766.777);
    theCorrectionsEE.push_back(0.0198988);   theBinsEE.push_back(860.335);
    theCorrectionsEE.push_back(0.0265556);   theBinsEE.push_back(965.308);
    theCorrectionsEE.push_back(0.0266586);   theBinsEE.push_back(1083.09);
    theCorrectionsEE.push_back(0.0308715);   theBinsEE.push_back(1215.24);
    theCorrectionsEE.push_back(0.0147255);   theBinsEE.push_back(1363.52);
    theCorrectionsEE.push_back(0.0155513);   theBinsEE.push_back(1529.9);
    theCorrectionsEE.push_back(0.0234095);   theBinsEE.push_back(1716.57);
    theCorrectionsEE.push_back(0.0319073);   theBinsEE.push_back(1926.02);
    theCorrectionsEE.push_back(0.0440146);   theBinsEE.push_back(2161.02);
    theCorrectionsEE.push_back(0.0401477);   theBinsEE.push_back(2424.71);
    theCorrectionsEE.push_back(0.0538419);   theBinsEE.push_back(2720.56);
    theCorrectionsEE.push_back(0.06414);   theBinsEE.push_back(3052.52);
    theCorrectionsEE.push_back(0.0769302);   theBinsEE.push_back(3424.98);
    theCorrectionsEE.push_back(0.0733546);   theBinsEE.push_back(3842.89);
    theCorrectionsEE.push_back(0.0753102);   theBinsEE.push_back(4311.79);
    theCorrectionsEE.push_back(0.0786608);   theBinsEE.push_back(4837.9);
    theCorrectionsEE.push_back(0.0960978);   theBinsEE.push_back(5428.21);
    theCorrectionsEE.push_back(0.111197);   theBinsEE.push_back(6090.55);
    theCorrectionsEE.push_back(0.108132);   theBinsEE.push_back(6833.7);
    theCorrectionsEE.push_back(0.12744);   theBinsEE.push_back(7667.54);
    theCorrectionsEE.push_back(0.177114);   theBinsEE.push_back(8603.12);
    theCorrectionsEE.push_back(0.168588);   theBinsEE.push_back(9652.85);
    theCorrectionsEE.push_back(0.178248);   theBinsEE.push_back(10830.7);
    theCorrectionsEE.push_back(0.286695);   theBinsEE.push_back(12152.2);
    theCorrectionsEE.push_back(0.649502);   theBinsEE.push_back(13635);
    theCorrectionsEE.push_back(0.561537);   theBinsEE.push_back(15298.7);
    theCorrectionsEE.push_back(0.315861);   theBinsEE.push_back(17165.5);
    theCorrectionsEE.push_back(-0.260055);   theBinsEE.push_back(19259.9);
    theCorrectionsEE.push_back(-0.485963);   theBinsEE.push_back(21610);
    theCorrectionsEE.push_back(-0.526979);   theBinsEE.push_back(24246.8);
    theCorrectionsEE.push_back(-0.320099);   theBinsEE.push_back(27205.4);
    //theCorrectionsEE.push_back(0);   theBinsEE.push_back(0.000463988);
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

   // these corrections are in nanosecond
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

      theCorrectionsEB.push_back(0.0491541);   theBinsEB.push_back(30.5506);
      theCorrectionsEB.push_back(0.0366708);   theBinsEB.push_back(34.2752);
      theCorrectionsEB.push_back(0.0329733);   theBinsEB.push_back(38.4542);
      theCorrectionsEB.push_back(0.0294081);   theBinsEB.push_back(43.1432);
      theCorrectionsEB.push_back(0.0259356);   theBinsEB.push_back(48.4044);
      theCorrectionsEB.push_back(0.0231976);   theBinsEB.push_back(54.3075);
      theCorrectionsEB.push_back(0.020153);   theBinsEB.push_back(60.9309);
      theCorrectionsEB.push_back(0.0189855);   theBinsEB.push_back(68.3624);
      theCorrectionsEB.push_back(0.018099);   theBinsEB.push_back(76.7008);
      theCorrectionsEB.push_back(0.0150657);   theBinsEB.push_back(86.0566);
      theCorrectionsEB.push_back(0.013667);   theBinsEB.push_back(96.5539);
      theCorrectionsEB.push_back(0.0127761);   theBinsEB.push_back(108.332);
      theCorrectionsEB.push_back(0.0142797);   theBinsEB.push_back(121.548);
      theCorrectionsEB.push_back(0.0134406);   theBinsEB.push_back(136.375);
      theCorrectionsEB.push_back(0.00871709);   theBinsEB.push_back(153.013);
      theCorrectionsEB.push_back(0.00774489);   theBinsEB.push_back(171.68);
      theCorrectionsEB.push_back(0.00916078);   theBinsEB.push_back(192.625);
      theCorrectionsEB.push_back(0.0071625);   theBinsEB.push_back(216.126);
      theCorrectionsEB.push_back(0.00582299);   theBinsEB.push_back(242.494);
      theCorrectionsEB.push_back(0.00974502);   theBinsEB.push_back(272.079);
      theCorrectionsEB.push_back(0.00911389);   theBinsEB.push_back(305.275);
      theCorrectionsEB.push_back(0.00685906);   theBinsEB.push_back(342.521);
      theCorrectionsEB.push_back(0.00427364);   theBinsEB.push_back(384.312);
      theCorrectionsEB.push_back(0.00683086);   theBinsEB.push_back(431.202);
      theCorrectionsEB.push_back(0.00458362);   theBinsEB.push_back(483.813);
      theCorrectionsEB.push_back(0.00442922);   theBinsEB.push_back(542.844);
      theCorrectionsEB.push_back(0.0111124);   theBinsEB.push_back(609.078);
      theCorrectionsEB.push_back(0.00586062);   theBinsEB.push_back(683.394);
      theCorrectionsEB.push_back(0.0156288);   theBinsEB.push_back(766.777);
      theCorrectionsEB.push_back(0.00229034);   theBinsEB.push_back(860.335);
      theCorrectionsEB.push_back(0.000420746);   theBinsEB.push_back(965.308);
      theCorrectionsEB.push_back(0.00168491);   theBinsEB.push_back(1083.09);
      theCorrectionsEB.push_back(-0.00289166);   theBinsEB.push_back(1215.24);
      theCorrectionsEB.push_back(0.000502757);   theBinsEB.push_back(1363.52);
      theCorrectionsEB.push_back(-0.0113229);   theBinsEB.push_back(1529.9);
      theCorrectionsEB.push_back(-0.0228308);   theBinsEB.push_back(1716.57);
      theCorrectionsEB.push_back(-0.0108323);   theBinsEB.push_back(1926.02);
      theCorrectionsEB.push_back(-0.0109995);   theBinsEB.push_back(2161.02);
      theCorrectionsEB.push_back(-0.034372);   theBinsEB.push_back(2424.71);
      theCorrectionsEB.push_back(-0.0414099);   theBinsEB.push_back(2720.56);
      theCorrectionsEB.push_back(-0.0361992);   theBinsEB.push_back(3052.52);
      theCorrectionsEB.push_back(-0.0293669);   theBinsEB.push_back(3424.98);
      theCorrectionsEB.push_back(0.00775366);   theBinsEB.push_back(3842.89);
      theCorrectionsEB.push_back(-0.150378);   theBinsEB.push_back(4311.79);
      theCorrectionsEB.push_back(-0.179345);   theBinsEB.push_back(4837.9);
      theCorrectionsEB.push_back(-0.221932);   theBinsEB.push_back(5428.21);
      theCorrectionsEB.push_back(-0.411671);   theBinsEB.push_back(6090.55);
      theCorrectionsEB.push_back(-0.565842);   theBinsEB.push_back(6833.7);
      theCorrectionsEB.push_back(-0.710823);   theBinsEB.push_back(7667.54);
      theCorrectionsEB.push_back(-0.560967);   theBinsEB.push_back(8603.12);
      theCorrectionsEB.push_back(-0.724621);   theBinsEB.push_back(9652.85);
      theCorrectionsEB.push_back(-0.786503);   theBinsEB.push_back(10830.7);
      theCorrectionsEB.push_back(-0.823316);   theBinsEB.push_back(12152.2);
      theCorrectionsEB.push_back(-0.868813);   theBinsEB.push_back(13635);
      theCorrectionsEB.push_back(-1.36518);   theBinsEB.push_back(15298.7);
      theCorrectionsEB.push_back(-1.06207);   theBinsEB.push_back(17165.5);
      //theCorrectionsEB.push_back(0);   theBinsEB.push_back(19259.9);
      //theCorrectionsEB.push_back(0);   theBinsEB.push_back(21610);
      //theCorrectionsEB.push_back(0);   theBinsEB.push_back(24246.8);
      //theCorrectionsEB.push_back(0);   theBinsEB.push_back(27205.4);
      //theCorrectionsEB.push_back(0);   theBinsEB.push_back(0.000589421);

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
