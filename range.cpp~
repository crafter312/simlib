#include "range.h"

range::range(string filename)
{
  filename = "/home/Argon30/simlib/loss/"+filename+".range";

  ifstream file(filename.c_str());
  if (!file)
    {
      cout << filename << " not found" << endl;
      abort();
    }

  string stuff;

  getline(file,stuff);
  cout << stuff << endl;

  getline(file,stuff);

  file >> N;

  x = new double [N];
  y = new double [N];
  z = new double [N];
  zz = new double [N];

  for (int i=0;i<N;i++)
    {
      file >> x[i] >> y[i] >> z[i] >> zz[i];
    }  
  file.close();
}

//******************************************************************

range::~range()
{
  delete [] x;
  delete [] y;
  delete [] z;
  delete [] zz;
}
//*************************************************************
double range::getRange(double E)
{
  if (E < x[0] ) 
    {
      cout << "energy = " << E << " below tabulation in range " << x[0] << endl;
      abort();
    }
  if (E > x[N-1])
    {
      cout << "energy above tabulation in range " << E <<  endl;
      abort();
    }

  int i = 1;
  for (;;)
    {
      if (E <= x[i]) break;
      i++;
    }


   straggle= (z[i]-z[i-1])/(x[i]-x[i-1])*(E-x[i-1]) + z[i-1];

   longStraggle = (zz[i]-zz[i-1])/(x[i]-x[i-1])*(E-x[i-1]) + zz[i-1];

  return     (y[i]-y[i-1])/(x[i]-x[i-1])*(E-x[i-1]) + y[i-1];
}
//***********************
double range::getLateralStraggle()
{
  return straggle;
}
//*************************
double range::getLongitudinalStraggle()
{
  return longStraggle;
}
