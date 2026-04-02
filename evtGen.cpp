#include "evtGen.h"

evtGen::evtGen(int Nx0,int Ny0,string filename, CRandom * ran0, double ET0, CFrame ** real0, double Mass_p0, double MassCore0, double gamma0)
{
  gamma = gamma0;
  Mass_p = Mass_p0;
  MassCore = MassCore0;
  real = real0;
  Nx = Nx0;
  Ny = Ny0;
  NN = Nx*Ny;
  ran = ran0;
  int n = 0;
  ET= ET0;
  array = new double [Nx*Ny];


  ifstream fileIn(filename.c_str());

  string name;
  getline(fileIn,name);
  cout << name << endl;

  getline(fileIn,name);
  cout << name << endl;

  getline(fileIn,name);
  cout << name << endl;

  getline(fileIn,name);
  cout << name << endl;


  for (int j=0;j<Ny;j++)
    {
      for (int i=0;i<Nx;i++)
	{
	  fileIn >> array[n];
          n++;
	}       
    }


  for (int i=1;i<NN;i++)
    {
      array[i] += array[i-1];
    }

  for (int i=0;i<NN;i++)
    {
      array[i] /= array[NN-1];
    }
}

//*****************************************

evtGen::~evtGen()
{
  delete array;
}
//****************************************
void evtGen::getEvent(bool uniform, double ET1)
{

  int iix,iiy;

  if (uniform)
    {
    iix = 50*ran->Rndm();
    iiy = 50*ran->Rndm();
    }
  else 
    {
  double x = ran->Rndm();
  int ii = 0;
  for (;;)
    {
      if (x  < array[ii]) break;
      ii++;
    }

  iix = ii/Nx;
  iiy = ii%Nx;
    }

  if (ET1 <= 0)
    {
  //cout << ET << " " << gamma << endl;
  ETT = 0.;
  for (;;)
    {
     ETT = ran->BreitWigner(ET,gamma);
     if (ETT > 0. && ETT < 15.) break;
    }
  //cout << "ETT = " << ETT << endl;
    }
  else ETT = ET1;

  Epp = ((double)iix+ran->Rndm())/(double)Nx *ETT;
  //if (Epp/ETT > .99)
  //  {
  //    cout << "here" << endl;
  //  }
  cosThetak = ((double)iiy+ran->Rndm())/(double)Ny*2.-1.;
  double thetak = acos(cosThetak);


  //core-diproton breakup
  double mass1 = MassCore/m0;
  double mass2 = (2.*Mass_p + Epp)/m0;
  double mu = mass1*mass2/(mass1+mass2);
  double Vrel = sqrt(2.*(ETT-Epp)/mu)*vfact;
  double v1  = mass2/(mass1+mass2)*Vrel;
  double v2  = Vrel - v1;


   double gamma1 = 1./sqrt(1.-pow(v1/c,2));
   float pc = gamma1*v1*mass1*m0/c;

   for (;;)
     {
       double E1 = sqrt(pow(pc,2)+pow(mass1*m0,2))- mass1*m0;
       double E2 = sqrt(pow(pc,2)+pow(mass2*m0,2))- mass2*m0;
       double y = E1 + E2  - (ETT-Epp);
       double dE1 = pc/sqrt(pow(pc,2)+pow(mass1*m0,2));
       double dE2 = pc/sqrt(pow(pc,2)+pow(mass2*m0,2));
       double dy = dE1 + dE2;
       double dpc = - y/dy;
       if (fabs(dpc) < .00001) break;
       pc += dpc;
     }
   v1 = pc/sqrt(pow(pc,2)+pow(mass1*m0,2))*c;
   v2 = pc/sqrt(pow(pc,2)+pow(mass2*m0,2))*c;



   real[2]->v[0] = 0.;
   real[2]->v[1] = 0.;
   real[2]->v[2] = -v1;

   CFrame di(2.+Epp/m0);
   di.v[0] = 0.;
   di.v[1] = 0.;
   di.v[2] = v2;


   double v2p[3] = {0.,0.,v2};

   mass1 = Mass_p/m0;
   mass2 = Mass_p/m0;
   mu = mass1*mass2/(mass1+mass2);
   Vrel = sqrt(2.*Epp/mu)*vfact;
   v1  = mass2/(mass1+mass2)*Vrel;
   v2  = Vrel - v1;


   gamma1 = 1./sqrt(1.-pow(v1/c,2));
   pc = gamma1*v1*mass1*m0/c;

   for (;;)
     {
       double E1 = sqrt(pow(pc,2)+pow(mass1*m0,2))- mass1*m0;
       double E2 = sqrt(pow(pc,2)+pow(mass2*m0,2))- mass2*m0;
       double y = E1 + E2  - Epp;
       double dE1 = pc/sqrt(pow(pc,2)+pow(mass1*m0,2));
       double dE2 = pc/sqrt(pow(pc,2)+pow(mass2*m0,2));
       double dy = dE1 + dE2;
       double dpc = - y/dy;
       if (fabs(dpc) < .00001) break;
       pc += dpc;
     }
   v1 = pc/sqrt(pow(pc,2)+pow(mass1*m0,2))*c;
   v2 = pc/sqrt(pow(pc,2)+pow(mass2*m0,2))*c;


   double phi = ran->Rndm()*2.*acos(-1.);
   
   real[0]->v[2] = v1*cosThetak;
   real[0]->v[0] = v1*cos(phi)*sin(thetak);
   real[0]->v[1] = v1*sin(phi)*sin(thetak);
  
   for (int i=0;i<3;i++) real[1]->v[i] = -real[0]->v[i];
   real[1]->transformVelocity(v2p);
   real[0]->transformVelocity(v2p);


   //check energy conservation
   // double EE1 = real[0]->getEnergy() + real[1]->getEnergy() + real[2]->getEnergy();

   //cout << " one " << EE1 << " " << ETT << endl;


   
   
   // randomize orientation of decay configuration
   //first in zx plane 
   double dtheta = acos(1.-2.*ran->Rndm());

   for (int k=0;k<3;k++)
     {

       double v2 = real[k]->v[2];
       double v0 = real[k]->v[0];

       real[k]->v[2] = v2*cos(dtheta) - v0*sin(dtheta);

       real[k]->v[0] = v0*cos(dtheta) + v2*sin(dtheta);
     }

   
   //then in the x y plane
   double dphi = acos(-1.)*2.*ran->Rndm();

   for (int k=0;k<3;k++)
     {
       double v0 = real[k]->v[0];
       double v1 = real[k]->v[1];

       real[k]->v[0] = v0*cos(dphi) + v1*sin(dphi);

       real[k]->v[1] = v1*cos(dphi) - v0*sin(dphi);
     }
   


   /*
   
   //check energy conservation
   double EE = real[0]->getEnergy() + real[1]->getEnergy() + real[2]->getEnergy();
  double ppc[3];
  for (int i=0;i<3;i++) ppc[i] = real[0]->pc[i] + real[1]->pc[i] + real[2]->pc[i];

  cout << EE << " " << ETT << " " << ppc[0] << " " << ppc[1] << " " << ppc[3] << endl;
   */
}
