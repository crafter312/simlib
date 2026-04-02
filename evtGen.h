#ifndef _evtGen
#define _evtGen
#include <fstream>
#include <string>
#include "random.h"
#include "frame.h"
#include <iostream>
#include "constants.h"


using namespace std;

class evtGen
{
 public:
  evtGen(int,int,string,CRandom*,double,CFrame**,double,double,double);
  ~evtGen();
  void getEvent(bool uniform, double ET1=-1.);
  int Nx,Ny;
  int NN;
  double ET;
  double gamma;
  CRandom * ran;

  double Mass_p;
  double MassCore;

  double Epp;
  double cosThetak;
  CFrame ** real;
  double ETT;
 
  private:
  double * array;
};
#endif
