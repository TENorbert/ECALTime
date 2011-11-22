#! /bin/sh

input=$1
output=$2

if [ -z $1 ]; then
   echo "Please specify the source file you want to compile"
   exit
fi

if [ -z $2 ]; then
   output=a.out
fi

cflags=`sh $ROOTSYS/bin/root-config --cflags`
libs=`sh $ROOTSYS/bin/root-config --libs`
glibs=`sh $ROOTSYS/bin/root-config --glibs`
#cxxflags="-g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2"
cxxflags="-g -fPIC -Wall -O -ansi -D_GNU_SOURCE -g -O2"
cxx="-m64"

g++ MathFunctions.cc -c -o MathFunctions.o $libs $cflags $cxx $cxxflags
g++ AnaInput.cc -c  -o AnaInput.o $libs $cflags $cxx $cxxflags
g++ hJetTime.cc -c  -o hJetTime.o $libs $cflags $cxx $cxxflags
g++ PhotonAna.cc -c  -o PhotonAna.o $libs $cflags $cxx $cxxflags

g++ $input -o $output MathFunctions.o AnaInput.o hJetTime.o PhotonAna.o $libs $cflags $cxx $cxxflags $glibs -lMinuit
