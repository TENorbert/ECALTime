#!/bin//bash
g++ FitHist.cc -o HFit `root-config --cflags --libs` 
#g++ -c `root−config --cflags --glibs` FitHist.cc
#g++ -o HFit `root-config --glibs` FitHist.o
