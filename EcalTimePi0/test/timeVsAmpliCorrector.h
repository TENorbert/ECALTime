#ifndef timeCorrector_H
#define timeCorrector_H

#include <assert.h>
#include <stdio.h>

#include <iostream>
#include <string>
#include <vector>

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

#endif

