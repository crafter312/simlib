#include "frag.h"
#include <cmath>


float const CFrag::pi=acos(-1.);
CRandom CFrag::ran;

/**
 * Constructor
 \param Z0 is the atomic number of fragment
 \param mass0 is mass of fragment in amu
 \param filename specifies the name of the file with energy loss information
 \param CsI_res0 fractional energy resolution of this fragment in CsI
 \param thickness is target thickness in mg/cm2 
 \param scaleSmallAngle - scale the predicted width for small angle scattering
 \param  
*/

CFrag::CFrag(float Z0,float mass0, string filename, float CsI_res0,
	     float thickness, float scaleSmallAngle0, bool useRealP0):fragment(Z0,mass0)
{
  ifstream ifile;
  ifile.open("/home/Oxygen11/simlib/teles.dat");

  extra = false;
  loss = new CLoss(filename,mass);
  CsI_res = CsI_res0;
  //be9 target  
  float thick = thickness/1000./9.*6.02e23; // atoms/ cem2
  multScat = new CMultScat((int)Z,4,thick);

  float fromCenter = 27.2;//35.;
  float radius = 63.; //50.;
  Array = new CArray(radius,fromCenter,6.42,7.6);
  shadowArray = new CArray(radius+2.,fromCenter,7.4,8.2);
  real = new CFrame(mass);
  recon = new CFrame(mass);
  scaleSmallAngle = scaleSmallAngle0;
  useRealP = useRealP0;

  const double dz = -10.44;
  //  Pixels = new pixels (dz);
  //Pixels->prepareSim();


  //defining variables for the teles
  float xcenter, ycenter, zcenter;
  float xhoriz, yhoriz, zhoriz;
  float xdiag, ydiag, zdiag;
  int itele;

  float activex = 6.4; //in cm
  float activey = 6.4; //in cm

  float target_offset = -0.0180025; //in m
  target_offset = -0.0225; //to be consistent with data

  float * rcenter[14];
  float * rback[14];
  float * rdiag[14];
  float * rnormal[14];
  float * rfront[14];

  for (int i=0; i<14; i++)
    {
      rcenter[i] = new float[3];
      rback[i] = new float[3];
      rdiag[i] = new float[3];
      rnormal[i] = new float[3];
      rfront[i] = new float[3];
      
      for (int j=0; j<3; j++)
	{
	  rcenter[i][j] = 0.;
	  rback[i][j] = 0.;
	  rdiag[i][j] = 0.;
	  rnormal[i][j] = 0.;
	  rfront[i][j] = 0.;
	}
    }

  while (ifile >> itele >> xcenter >> ycenter >> zcenter >> xhoriz >> yhoriz >> zhoriz >> xdiag >> ydiag >> zdiag)
    {
      zcenter += target_offset;
      
      findVectors(rcenter[itele], rback[itele], rdiag[itele], rnormal[itele], rfront[itele], xcenter, ycenter, zcenter, xhoriz, yhoriz, zhoriz, xdiag, ydiag, zdiag);
      
      Teles[itele] = new CTele(rcenter[itele], rfront[itele], rback[itele], activex, activey);
      // Teles[itele]->getCsICenters(itele);
    }


  //proton thresholds
  ifstream pfile("/home/Oxygen11/simlib/pthres.dat");
  int one;
  float two;
  for (int i=0;i<56;i++)
    {
      pfile >> one >> two;
      p_thres[i] = two;
      if (one != i) cout << " bad " << endl;
    }

  pfile.close();


}


//*********************************************************
/**
*Destructor
*/
CFrag::~CFrag()
{
  delete Array;
  delete real;
  delete recon;
  //delete Pixels;
}

//**********************************************************
  /**
   * logical function determines if a particle is detected
   */
int CFrag::hit()
{
  is_hit = Array->hit(real->theta,real->phi,(float)0.,(float)0.) ;

  if (is_hit) 
    {

      recon->theta = Array->thetaRecon;
      recon->phi = Array->phiRecon;

      recon->energy = real->energy + sqrt(real->energy)*CsI_res*
       ran.Gaus(0.,1.);




     recon->getVelocity();

    }

  return is_hit;
}
//*******************************************************
int CFrag::hitShadow(float xtarget , float ytarget)
{
  is_hit = shadowArray->hit(real->theta,real->phi,xtarget,ytarget) ;

  return is_hit;
}


//**********************************************************
  /**
   * logical function determines if a particle is detected
\param xTarget is in-plane location of interection in target form nomimal center
\param yTarget is out-of-plane location (cm)
  */
int CFrag::hit(float xTarget, float yTarget)
{
  is_hit = Array->hit(real->theta,real->phi,xTarget, yTarget) ;

  if (is_hit)
    {
      
      recon->theta = Array->thetaRecon;
      recon->phi = Array->phiRecon;
     recon->energy = real->energy + sqrt(real->energy)*CsI_res*
       ran.Gaus(0.,1.);
     
     /*
      recon->theta = real->theta;
      recon->phi = real->phi;
      
      recon->energy = real->energy;
     */
     recon->getVelocity();

    }

  return is_hit;
}
//******************************************************************
  /**
   * Add a velocity vector to the fragments velocity vector.
   * Used to transform between reference frames
   */
void CFrag::AddVelocity(double *Vplf)
{
  real->transformVelocity(Vplf);
}
//****************************************************************
  /** 
   * returns the energy after the fragment has exited the target
\param thick is the thickness of target that the particle has to traverse (mg/cm2)
  */
float CFrag::Eloss(float thick)
{
  if (real->energy <=0. || real->energy != real->energy) return 0.;
  real->energy = loss->getEout(real->energy,thick);
  return real->energy;
}
//*******************************************************************
  /**
   * corrects energy of a detected particle for the energy loss
   * in the target.
\param thick is the thickness of target material though which the particle passed (mg/cm2)
  */
float CFrag::Egain(float thick)
{
  if (thick > 0.)
    recon->energy = loss->getEin(recon->energy,thick/cos(recon->theta));


  recon->getVelocity();

   return recon->energy;
}
//***********************************************
//include multiple scattering
  /**
   * Monte Carlo choice of small angle scattering due to passage through the target
\param fractionalThick is the fractional thick of the target through which the particle passed
  */
void CFrag::MultiScat(float fractionalThick)
{

  float thetaRMS = multScat->thetaRMS(real->energy,fractionalThick);
  float sigma = thetaRMS/sqrt(2.)*scaleSmallAngle;
  //cout << "thetaRMS= " << thetaRMS << endl;
  float deltaTheta = sqrt(2.)*sigma*sqrt(-log(ran.Rndm()));
  //cout << "deltaTheta= " << deltaTheta << endl;
  float deltaPhi = 2.*pi*ran.Rndm();
  //cout << "delta Phi= " << deltaPhi << endl;

  float x = sin(deltaTheta)*cos(deltaPhi);
  float y = sin(deltaTheta)*sin(deltaPhi);
  float z = cos(deltaTheta);



  //rotate in z-x plane by theta
  float xx = x*cos(real->theta) + z*sin(real->theta);
  float yy = y;
  float zz = z*cos(real->theta) - x*sin(real->theta);


  //rotate in x-y plane
  float xxx = xx*cos(real->phi) - yy*sin(real->phi);
  float yyy = yy*cos(real->phi) + xx*sin(real->phi);
  float zzz = zz;


  float thetaNew = acos(zzz);
  float phiNew = atan2(yyy,xxx);


  real->theta = thetaNew;
  real->phi = phiNew;
}
//*********************
  /**
   * accounts for multiscattering and energy loss in the target
   \param dthick is thickness of target though the particle must pass (mg/cm2)
\param thickness is total target thickness (mg/cm2)
   */
bool CFrag::targetInteraction(float dthick, float thickness)
{

  bool stopped = 0;
  if (dthick == 0.)
    {
    return stopped;
    }
  float thick = dthick/cos(real->theta);
  Eloss(thick);
  if (real->energy <= 0. || real->energy != real->energy) 
    { 
      stopped = 1;

      return stopped;
    }


  MultiScat(thick/thickness);

  return stopped;
}
//**********************************************************
  /**
   * logical function determines if a particle is detected
\param xTarget is in-plane location of interection in target form nomimal center
\param yTarget is out-of-plane location (cm)
  */
int CFrag::hit2(float xTarget, float yTarget)
{
  //cout << "my wife" << endl;
  for(int i=0; i<14; i++)
    {
      is_hit = Teles[i]->hit(real->theta,real->phi,xTarget,yTarget);
      if(is_hit)
	{
	  hitTele = i;
	  break;
	}
    }

  if (!is_hit) return is_hit;



 if (fabs(real->A- 1.) < 0.5)
    {
      if (!protonHole(hitTele,Teles[hitTele]->ixStrip,Teles[hitTele]->iyStrip)) 
       is_hit = 0;
     }



 if (real->A > 8.5 && real->A < 11.5)
    {
      if (!carbonHole(hitTele,Teles[hitTele]->ixStrip,Teles[hitTele]->iyStrip)) 
       is_hit = 0;
     }


 if (real->A > 11.5)
    {
      if (!nitrogenHole(hitTele,Teles[hitTele]->ixStrip,Teles[hitTele]->iyStrip)) 
       is_hit = 0;
     }

   if (fabs(real->A-4.) < 0.5)
     {
       if (!alphaHole(hitTele,Teles[hitTele]->ixStrip,Teles[hitTele]->iyStrip)) is_hit = 0;
     }

  if(real->A >10.)
    {
      //    if (itele !=6  && itele !=7) 
      //is_hit =0;
      //only have else if using inner 4 CsI not inner 8
//       else
// 	{
// 	  if(itele ==6)
// 	    {
// 	      if(Pixels->ICsI == 0 || Pixels->ICsI ==3) is_hit=0;
// 	    }
// 	  else
// 	    {
// 	      if(Pixels->ICsI == 1 || Pixels->ICsI ==2) is_hit=0;
// 	    }
// 	}

    }
  if (is_hit)
    {
      
      recon->theta = Teles[hitTele]->thetaRecon;
      recon->phi = Teles[hitTele]->phiRecon;

      recon->energy = real->energy + sqrt(real->energy)*CsI_res*
	ran.Gaus(0.,1.);
      
      if (useRealP)
	{
         recon->theta = real->theta;
         recon->phi = real->phi;
         recon->energy = real->energy;
	}

      recon->getVelocity();

    }

  // if(fabs(real->A-1.)< 0.5 && recon->energy < p_thres[itele]) is_hit = 0;

  return is_hit;




}
//**************************************************
bool CFrag::alphaHole(int itele, int ifront, int iback)
{


  if (itele == 0) //det0 done
    {
      float x = 1.;

      if (iback == 15) x *= .82;
      if (iback == 16) x *= .88;
      if (ifront == 14) x *= .95;
      if (ifront == 15) x *= .68;
      if (ifront == 16) x *= .92;
      if (ifront == 31) x *= .94;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 1)//det1 done
    {
      float x = 1.;
      if (iback == 15) x *= .92;
      if (iback == 16) x *= .69;
      if (ifront == 14) x *= .88;
      if (ifront == 15) x *= .75;
      if (ifront == 16) x *= .95;
      if (ifront == 31) x *= .88;
      if (ran.Rndm() > x) return false;
    }    

  else if (itele == 2)//det2 done
    {
      float x = 1.;
      if (iback == 15) x *= .67;
      if (iback == 16) x *= .95;
      if (ifront == 4) return false;
      if (ifront == 7) return false;
      if (ifront == 15) x*= .70;
      if (ifront == 16) x*= .88;
      if (ifront == 18) return false;
      if (ifront == 19) return false;
      if (ifront == 31) x*=.96;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 3) //det3 done
    {
      float x = 1.;
      if (iback == 15) x*= .84;
      if (iback == 16) x*= .78;
      if (ifront == 15) x*= .82;
      if (ifront == 16) x*= .79;
      if (ifront == 31) x*= .98;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 4) //det4  done
    {
      float x = 1.;
      if (iback == 15) x*= .92;
      if (iback == 16) x*= .67;
      if (ifront == 15) x*=.74;
      if (ifront == 16) x*=.88;
      if (ifront == 30) x*= .65;
      if (ifront == 31) x*= .51;
      if (ran.Rndm() > x) return false;

    }
  else if (itele == 5) //det5 done
    {
      float x = 1.;
      if (iback == 15) x*=.89;
      if (iback == 16) x*=.67;
      if (iback == 20) return false;
      if (iback == 25) return false;
      if (iback == 27) return false;
      if (ifront == 15) x*=.83;
      if (ifront == 16) x*=.77;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 6) //det 6 done 
    {
      float x = 1.;
      if (iback == 15) x*=.64;
      if (iback == 16) x*=.87;
      if (ifront == 0) x*=.97;
      if (ifront == 15) x*=.82;
      if (ifront == 15) x*=.66;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 7) //det 7 done
    {
      float x = 1.;
      if (iback == 15) x*=.90;
      if (iback == 16) x*=.64;
      if (iback == 21) return false;
      if (ifront == 15) x*=.72;
      if (ifront == 16) x*=.75;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 8) //det 8 done
    {
      float x = 1.;
      if (iback == 15) x*=.82;
      if (iback == 16) x*= .90;
      if (ifront == 4) return false;
      if (ifront == 6) return false;
      if (ifront == 7) return false;

      if (ifront == 15) x*= .7;
      if (ifront == 16) x*= .88;
      if (ifront == 31) x*= .98;

      if (ran.Rndm() > x) return false;
    }
  else if (itele == 9) //det 9 done
    {

      float x = 1.;
      if (iback == 15) x *= .59;
      if (iback == 16) x*= .92;
      if (ifront == 15)x*=.97;
      if (ifront == 16)x*=.60;
      if (ifront == 31)x*=.95;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 10) //det 10 done 
    {
      float x = 1.;
      if (iback == 15) x*= .81;
      if (iback == 16) x*= .77;
      if (iback == 31) x*= .98;
      if (ifront == 0) x*=.94;
      if (ifront == 2) x*=.094;
      if (ifront == 15) x*=.92;
      if (ifront == 16) x*=.83;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 11) //det 11 done
    {
      float x = 1.;
      if (iback == 15) x*= .94;
      if (iback == 16) x*= .67;
      if (ifront == 15) x*=.96;
      if (ifront == 16) x*= .66;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 12) //det 12  done
    {
      float x = 1.;
      if (iback == 16) x*=.70;
      if (ifront == 0) x*=.97;
      if (ifront == 16) x*=.70;
      if (ifront == 16) x*=.98;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 13) //det 13 
    {
      float x = 1.;
      if (iback == 15) x*= .94;
      if (iback == 16) x*=.67;
      if (iback == 19) return false;
      if (iback == 31) return false;

      if (ifront == 0) x*=.96;
      if (ifront == 16) x*=.68;
      if (ifront == 31) return false;
      if (ran.Rndm() > x) return false;
    }

  return true;
}
//**************************************************
bool CFrag::protonHole(int itele, int ifront, int iback)
{


  if (itele == 0) //det0 done
    {
      float x = 1.;

      if (iback == 14) x *= .98;
      if (iback == 15) x *= .77;
      if (iback == 16) x *= .82;
      if (iback == 17) x *= .99;
      if (iback == 31) x *= .96;
      if (ifront == 13) x *= .98;
      if (ifront == 14) x *= .83;
      if (ifront == 15) x *= .68;
      if (ifront == 16) x *= .95;
      if (ifront == 31) x *= .95;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 1)//det1 done
    {
      float x = 1.;
      if (iback == 0) x *= .98;
      if (iback == 15) x *= .86;
      if (iback == 16) x *= .65;
      if (iback == 17) x *= .98;
      if (iback == 31) x *= .98;
      if (ifront == 13) x *= .98;
      if (ifront == 14) x *= .76;
      if (ifront == 15) x *= .75;
      if (ifront == 15) x *= .96;
      if (ifront == 31) x *= .92;
      if (ran.Rndm() > x) return false;
    }    

  else if (itele == 2)//det2 done
    {
      float x = 1.;
      if (iback == 0) x *= .97;
      if (iback == 14) x *= .93;
      if (iback == 15) x *= .70;
      if (iback == 16) x *= .92;
      if (iback == 31) x*= .98;
      if (ifront == 4) return false;
      if (ifront == 7) return false;
      if (ifront == 14) x*= .97;
      if (ifront == 15) x*= .67;
      if (ifront == 16) x*= .86;
      if (ifront == 18) return false;
      if (ifront == 19) return false;
      if (ifront == 31) x*=.97;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 3) //det3 done
    {
      float x = 1.;
      if (iback == 14) x *= .99;
      if (iback == 15) x*= .81;
      if (iback == 16) x*= .74;
      if (iback == 17) x*= .99;
      if (iback == 31) x*= .98;
      if (ifront == 0) x*= .95;
      if (ifront == 14) x*= .98;
      if (ifront == 15) x*= .74;
      if (ifront == 16) x*= .75;
      if (ifront == 31) x*= .98;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 4) //det4 done
    {
      float x = 1.;
      if (iback == 15) x *= .88;
      if (iback == 16) x*= .64;
      if (iback == 17) x*= .76;
      if (ifront == 14) x*=.97;
      if (ifront == 15) x*=.69;
      if (ifront == 16) x*=.84;
      if (ran.Rndm() > x) return false;

    }
  else if (itele == 5) //det5 done
    {
      float x = 1.;
      if (iback == 13) x*= .98;
      if (iback == 14) x*= .75;
      if (iback == 15) x*=.77;
      if (iback == 16) x*=.95;
      if (iback == 20) return false;
      if (iback == 25) return false;
      if (iback == 27) return false;
      if (iback == 31) x*=.95;
      if (ifront == 0) x*=.95;
      if (ifront == 14) x*=.98;
      if (ifront == 15) x*=.76;
      if (ifront == 16) x*=.74;
      if (ifront == 17) x*=.98;
      if (ifront == 31) x*=.94;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 6) //det 6 done
    {
      float x = 1.;
      if (iback == 14) x*= .96;
      if (iback == 15) x*=.66;
      if (iback == 16) x*=.82;
      if (iback == 17) x*=.98;
      if (iback == 31) x*=.99;
      if (ifront == 0) x*=.97;
      if (ifront == 14) x*=.99;
      if (ifront == 15) x*=.77;
      if (ifront == 16) x*=.66;
      if (ifront == 17) x*=.99;
      if (ifront == 31) x*=.99;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 7) //det 7
    {

      float x = 1.;
      if (iback == 15) x*=.87;
      if (iback == 16) x*=.60;
      if (iback == 17) x*=.94;
      if (iback == 21) return false;
      if (ifront == 14) x*=.98;
      if (ifront == 15) x*=.68;
      if (ifront == 16) x*=.73;
      if (ifront == 17) x*=.99;

      if (ran.Rndm() > x) return false;
    }
  else if (itele == 8) //det 8 done
    {
      float x = 1.;
      if (iback == 4) x *= .67;
      if (iback == 15) x*=.95;
      if (iback == 16) x*=.77;
      if (iback == 17) x*=.79;
      if (iback == 26) x*=.68;
      if (ifront == 4) return false;
      if (ifront == 6) return false;
      if (ifront == 7) return false;
      if (ifront == 15) x*= .66;
      if (ifront == 16) x*= .84;
      if (ifront == 31) x*= .96;

      if (ran.Rndm() > x) return false;
    }
  else if (itele == 9) //det 9 done
    {

      float x = 1.;
      if (iback == 0) x *= .97;
      if (iback == 13) x *= .99;
      if (iback == 14) x *= .90;
      if (iback == 15) x *= .62;
      if (iback == 16) x*= .91;
      if (iback == 22) return false;
      if (iback == 31) x*=.96;
      if (ifront == 0)x*=.95;
      if (ifront == 15)x*=.89;
      if (ifront == 16)x*=.60;
      if (ifront == 17)x*=.97;
      if (ifront == 31)x*=.93;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 10) //det 10 done
    {
      float x = 1.;
      if (iback == 14) x *= .98;
      if (iback == 15) x*= .75;
      if (iback == 16) x*= .75;
      if (iback == 17) x*= .99;
      if (ifront == 0) x*=.96;
      if (ifront == 2) x*=.12;
      if (ifront == 14) x*=.98;
      if (ifront == 15) x*=.86;
      if (ifront == 16) x*=.68;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 11) //det 11 done
    {
      float x = 1.;
      if (iback == 15) x*= .9;
      if (iback == 16) x*=.63;
      if (iback == 17) x*=.95;
      if (ifront == 15) x*=.88;
      if (ifront == 16) x*= .64;
      if (ifront == 17) x*=.96;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 12) //det 12  done
    {
      float x = 1.;
      if (iback == 14) x*= .97;
      if (iback == 15) x*= .74;
      if (iback == 16) x*=.86;
      if (ifront == 14) x*=.97;
      if (ifront == 5) x*=.72;
      if (ifront == 16) x*=.89;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 13) //det 13 done
    {
      float x = 1.;
      if (ifront == 0) x*=.9;
      if (ifront == 16) x*=.78;
      if (ifront == 17) x*=.96;
      if (ifront == 31) return false;
      if (ran.Rndm() > x) return false;
    }

  return true;
}



//**************************************************
bool CFrag::carbonHole(int itele, int ifront, int iback)
{


  if (itele == 0) //det0 done
    {
      float x = 1.;

      if (iback == 14) x *= .90;
      if (iback == 15) x *= .90;
      if (iback == 20) x *= .70;
      if (iback == 23) x *= .54;
      if (ifront == 31) x *= .82;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 1)//det1 done
    {
      float x = 1.;
      if (iback == 16) x *= .79;
      if (ifront == 0) x *= .71;
      if (ifront == 15) x *= .74;
      if (ifront == 31) x *= .81;
      if (ran.Rndm() > x) return false;
    }    

  else if (itele == 2)//det2 done
    {
      float x = 1.;
      if (iback == 15) x *= .88;
      if (ifront == 4) return false;
      if (ifront == 6) return false;
      if (ifront == 18) return false;
      if (ifront == 19) return false;
      if (ifront == 31) x*=.84;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 3) //det3 done
    {
      float x = 1.;
      if (iback == 0) x *= .96;
      if (iback == 15) x*= .92;
      if (iback == 16) x*= .87;
      if (iback == 17) x*= .98;
      if (ifront == 0) x*= .95;
      if (ifront == 14) x*= .4;
      if (ifront == 23) x*= .34;
      if (ifront == 31) x*= .94;
      if (ran.Rndm() > x) return false;

      if (extra)
	{
	  if (ifront == 31 ) return false;
          if (ifront == 16 ) return false;
          if (iback == 16) return false;
	}


    }
  else if (itele == 4) //det4 done
    {
      float x = 1.;
      if (iback == 15) x *= .85;
      if (ifront == 15) x*=.88;
      if (ifront == 23) x*=.82;
      if (ifront == 30) x*=.71;
      if (ifront == 31) return false;;
      if (ran.Rndm() > x) return false;

    }
  else if (itele == 5) //det5 done
    {
      float x = 1.;
      if (iback == 16) x*= .85;
      if (iback == 27) return false;
      if (ifront == 0) x*=.83;
      if (ifront == 16) x*=.86;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 6) //det 6 done
    {
      float x = 1.;
      if (iback == 15) x*=.9;
      if (iback == 16) x*=.88;
      if (ifront == 15) x*=.91;
      if (ifront == 16) x*=.81;
      if (ifront == 31) x*=.9;
      if (ran.Rndm() > x) return false;
      if (extra)
	{
	  if (iback == 31) return false;
	  if (iback == 15) return false; 
	}

    }
  else if (itele == 7) //det 7
    {

      float x = 1.;
      if (ifront == 0) x*=.9;
      if (ifront == 15) x*=.81;
      if (ifront == 16) x*=.90;
      if (ifront == 31) x*=.96;

      if (ran.Rndm() > x) return false;

      if (extra)
	{
	  if (iback == 0) return false;
	  if (iback == 16) return false;
	  if (ifront == 15) return false;
	  if (ifront == 31) return false;
	}
    }
  else if (itele == 8) //det 8 done
    {
      float x = 1.;
      if (iback == 4) x *= .65;
      if (iback == 16) x*=.86;
      if (iback == 26) x*=.69;
      if (ifront == 4) return false;
      if (ifront == 5) x*= .74;
      if (ifront == 6) return false;
      if (ifront == 7) return false;
      if (ifront == 15) x*= .57;
      if (ifront == 31) x*= .85;

      if (ran.Rndm() > x) return false;
    }
  else if (itele == 9) //det 9 done
    {

      float x = 1.;
      if (ifront == 16)x*=.74;
      if (ifront == 31)x*=.62;
      if (ran.Rndm() > x) return false;
      if (extra)
	{
	  if (iback == 31) return false;
	}
    }
  else if (itele == 10) //det 10 done
    {
      float x = 1.;
      if (iback == 15) x*= .98;
      if (iback == 16) x*= .8;
      if (iback == 31) x*= .96;
      if (ifront == 0) x*=.94;
      if (ifront == 2) x*=.09;
      if (ifront == 16) x*=.85;
      if (ifront == 23) x*=.84;
      if (ran.Rndm() > x) return false;
      if (extra)
	{
	  if (iback == 0) return false;
          if (iback == 16) return false;
          if (ifront == 0) return false;
	}
    }
  else if (itele == 11) //det 11 done
    {
      float x = 1.;
      if (iback == 16) x*=.79;
      if (ifront == 0) x*=.92;
      if (ifront == 16) x*= .83;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 12) //det 12  done
    {
      float x = 1.;
      if (iback == 15) x*= .92;
      if (iback == 16) x*=.90;
      if (ifront == 1) x*=.66;
      if (ifront == 2) x*=.78;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 13) //det 13 done
    {
      float x = 1.;
      if (iback == 16) x*= .7;
      if (iback == 19) x*= .5;
      if (iback == 31) return false;;
      if (iback == 30) x*= .58;
      if (ifront == 0) x*=.77;
      if (ifront == 7) x*=.78;
      if (ifront == 19) x*=.65;
      if (ifront == 31) return false;
      if (ran.Rndm() > x) return false;
    }

  return true;
}


//**************************************************
bool CFrag::nitrogenHole(int itele, int ifront, int iback)
{


  if (itele == 0) //det0 done
    {
      float x = 1.;

      if (iback == 14) x *= .90;
      if (iback == 15) x *= .90;
      if (iback == 20) return false;
      if (iback == 23) return false;
      if (ifront == 31) x *= .82;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 1)//det1 done
    {
      float x = 1.;
      if (iback == 16) x *= .79;
      if (ifront == 0) x *= .71;
      if (ifront == 15) x *= .74;
      if (ifront == 31) x *= .81;
      if (ran.Rndm() > x) return false;
    }    

  else if (itele == 2)//det2 done
    {
      float x = 1.;
      if (iback == 15) x *= .88;
      if (ifront == 4) return false;
      if (ifront == 6) return false;
      if (ifront == 18) return false;
      if (ifront == 19) return false;
      if (ifront == 31) x*=.84;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 3) //det3 done
    {
      float x = 1.;
      if (iback == 0) x *= .96;
      if (iback == 15) x*= .79;
      if (iback == 31) x*= .90;
      if (ifront == 0) x*= .95;
      if (ifront == 14) return false;
      if (ifront == 23) x*= .08;
      if (ifront == 31) x*= .94;
      if (ran.Rndm() > x) return false;

      if (extra)
	{
	  if (ifront == 31 ) return false;
          if (ifront == 16 ) return false;
          if (iback == 16) return false;
	}


    }
  else if (itele == 4) //det4 done
    {
      float x = 1.;
      if (iback == 15) x *= .85;
      if (ifront == 15) x*=.88;
      if (ifront == 23) x*=.82;
      if (ifront == 30) x*=.71;
      if (ifront == 31) return false;;
      if (ran.Rndm() > x) return false;

    }
  else if (itele == 5) //det5 done
    {
      float x = 1.;
      if (iback == 16) x*= .85;
      if (iback == 27) return false;
      if (iback == 20) return false;
      if (ifront == 0) x*=.71;
      if (ifront == 16) x*=.88;
      if (ifront == 31) x*=.67;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 6) //det 6 done
    {
      float x = 1.;
      if (iback == 15) x*=.9;
      if (iback == 16) x*=.88;
      if (ifront == 15) x*=.89;
      if (ifront == 16) x*=.82;
      if (ifront == 31) x*=.85;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 7) //det 7
    {

      float x = 1.;
      if (iback == 16) x *= .78;
      if (iback == 21) return false;
      if (ifront == 0) x*=.91;
      if (ifront == 15) x*=.83;
      if (ifront == 16) x*=.92;
      if (ifront == 31) x*=.96;

      if (ran.Rndm() > x) return false;

    }
  else if (itele == 8) //det 8 done
    {
      float x = 1.;
      if (iback == 4) return false;

      if (iback == 16) x*=.86;
      if (iback == 26) return false;

      if (ran.Rndm() > x) return false;
    }
  else if (itele == 9) //det 9 done
    {

      float x = 1.;
      if (iback == 31) x*=.89;
      if (iback == 15) x*= .79;
      if (ifront == 16)x*=.74;
      if (ran.Rndm() > x) return false;

    }
  else if (itele == 10) //det 10 done
    {
      float x = 1.;
      if (iback == 15) x*= .79;
      if (iback == 16) x*= .8;
      if (iback == 31) x*= .96;
      if (ifront == 2) return false;
      if (ifront == 16) x*=.85;
      if (ifront == 23) x*=.84;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 11) //det 11 done
    {
      float x = 1.;
      if (iback == 16) x*=.72;
      if (ifront == 1) x*=.34;
      if (ifront == 16) x*= .83;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 12) //det 12  done
    {
      float x = 1.;
      if (iback == 15) x*= .92;
      if (iback == 16) x*=.90;
      if (ifront == 1) x*=.66;
      if (ifront == 2) x*=.78;
      if (ran.Rndm() > x) return false;
    }
  else if (itele == 13) //det 13 done
    {
      float x = 1.;
      if (iback == 16) x*= .7;
      if (iback == 19) return false;
      if (iback == 31) return false;
      if (iback == 30) return false;
      if (ran.Rndm() > x) return false;
    }

  return true;
}


void CFrag::findVectors(float *rcenter, float *rback, float *rdiag, float *rnormal, float *rfront, float xcenter0, float ycenter0, float zcenter0, float xhoriz0, float yhoriz0, float zhoriz0, float xdiag0, float ydiag0, float zdiag0)
{
  float m_to_cm = 100.; //to convert from meters to centimeters

  rcenter[0] = m_to_cm*xcenter0;
  rcenter[1] = m_to_cm*ycenter0;
  rcenter[2] = m_to_cm*zcenter0;

  float rfront_mag = sqrt(pow(xhoriz0,2) + pow(yhoriz0,2) + pow(zhoriz0,2));

  rfront[0] = xhoriz0/rfront_mag; //go left to right (beam right point - beam left point)
  rfront[1] = yhoriz0/rfront_mag; 
  rfront[2] = zhoriz0/rfront_mag;

  rdiag[0] = -xdiag0; // beam left point  - beam right point (this is wrong, this is a later comment)
  rdiag[1] = -ydiag0; 
  rdiag[2] = -zdiag0;

  rnormal[0] = rfront[1]*rdiag[2] - rfront[2]*rdiag[1];
  rnormal[1] = rfront[2]*rdiag[0] - rfront[0]*rdiag[2];
  rnormal[2] = rfront[0]*rdiag[1] - rfront[1]*rdiag[0];

  float rnormal_mag = sqrt(pow(rnormal[0],2) + pow(rnormal[1],2) + pow(rnormal[2],2));

  rnormal[0] = rnormal[0]/rnormal_mag;
  rnormal[1] = rnormal[1]/rnormal_mag;
  rnormal[2] = rnormal[2]/rnormal_mag;

  //rback = rfront x rnormal

  rback[0] = rfront[1]*rnormal[2] - rfront[2]*rnormal[1];
  rback[1] = rfront[2]*rnormal[0] - rfront[0]*rnormal[2];
  rback[2] = rfront[0]*rnormal[1] - rfront[1]*rnormal[0];

  float rback_mag = sqrt(pow(rback[0],2) + pow(rback[1],2) + pow(rback[2],2));

  rback[0] = rback[0]/rback_mag;
  rback[1] = rback[1]/rback_mag;  
  rback[2] = rback[2]/rback_mag;
}
