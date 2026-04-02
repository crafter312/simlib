#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;

class range{
 public:

  double *x;
  double *y;
  double *z;
  double *zz;
  int N;

  range(string filename);
  ~range();

  double getRange(double);
  double getLateralStraggle();
  double getLongitudinalStraggle();

  double straggle;
  double longStraggle;
};
