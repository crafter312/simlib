#include <string>
#include <fstream>
#include <iostream>
#include "frame.h"

using namespace std;


class moscow
{
  private:
  int ntot;

  static const float constant;
  ifstream  ifFile;
 public:
  double A1;
  double A2;
  double A3;
  moscow(string);
  bool getEvent(float,CFrame*,CFrame*,CFrame*);
  bool getEventP(CFrame*,CFrame*,CFrame*);
  bool getEvent8C(float,CFrame*,CFrame*,CFrame*);
};
