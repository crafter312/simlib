#include <math.h>



#include "Profile.h"

const double Profile::hbarc=197.326938;
//const double Profile::amu=931.494088;
const double Profile::amu=931.5;
const double Profile::pi=acos(-1.);

double Profile::Sommerfeld(double E, double mu, int z1, int z2)
{
  return z1 * z2 * sqrt(mu / E) / 6.3;
}

double Profile::P_l(double E, double r, double l, double mu, int z1, int z2)
{

  double eta = Sommerfeld(E, mu, z1, z2);
  double rho = r * sqrt(2. * mu * amu * E) / hbarc;
  //cout << rho << endl;
  gsl_sf_result F, G, Fp, Gp;
  double fe, ge;
  int err = gsl_sf_coulomb_wave_FG_e(eta, rho, l, 0, &F, &Fp, &G, &Gp, &fe, &ge);
  if(err != 0)
    {
      return 0;
    }
  return rho/(pow(F.val,2) + pow(G.val,2));
}

double Profile::Gamma(double E, double Et, double r, double l, double mu, int z1, int z2, double rwidth2)
{
  return 2 * P_l(Et - E, r, l, mu, z1, z2) * rwidth2;
}

double Profile::Shift(double E, double r, double l, double mu, int z1, int z2)
{
  if (E > 0.) 
    {
     double eta = Sommerfeld(E, mu, z1, z2);
     double rho = r * sqrt(2.*mu*amu*E)/hbarc;
     gsl_sf_result F, G, Fp, Gp;
     double fe, ge;


      if(gsl_sf_coulomb_wave_FG_e(eta, rho, l, 0, &F, &Fp, &G, &Gp, &fe, &ge) != 0)
       return 0;

     return rho * (F.val * Fp.val + G.val * Gp.val) / (pow(F.val,2)+pow(G.val,2));
    }
  else 
    {
     double eta = Sommerfeld(-E, mu, z1, z2);
     double rho = r * sqrt(-2.*mu*amu*E)/hbarc;
     

     double U0 = gsl_sf_hyperg_U((double)l+1.+eta,2.+2.*(double)l,2*rho);
     double U1 = gsl_sf_hyperg_U((double)l+2.+eta,3.+2.*(double)l,2*rho);
     //cout << "eta " << eta << " " << U0 << " " << U1 << " " << E << " " << mu << " " << z1 << " " << z2 << endl;
 
     return (((double)l+1.-rho)*U0-2.*((double)l+1.+eta)*rho*U1)/U0;
    }

}


//profile of a sequential decay
Profile::Profile(double Et0, double Er, int z1, int z2, int z3, double mu123, double ac123, int l1, double mu23, double ac23, int l2, double rwidth2_23, double B23)
{

  Et = Et0;
  gsl_set_error_handler_off();
  dE = .001;
  

  if (B23 == 999.) 
    {
     B23= Shift(Er,ac23,l2,mu23,z2,z3);

    }
     cout << "B23= " << B23 <<endl;

  N = (int)floor(Et/dE);
  profile = new double[N];

  //find GammaTot
  double ymax = 0.;
  double Emax = 0;
  double sumY = 0.;
  for (int i=0;i<10*N;i++)
    {
      double E = ((double)i+.5)*dE;
      double G2 = Gamma(0.,E,ac23,l2,mu23,z2,z3,rwidth2_23);
      double D2 = -rwidth2_23*(Shift(E,ac23,l2,mu23,z2,z3)- B23);

      double Y = G2/(pow(E-Er-D2,2)+pow(G2/2.,2));
      if (Y > ymax)
	{
         ymax = Y;
	 Emax = E;
	}
      sumY += Y*dE;
      //cout << "E= " << E << " Y= " << Y << " G2= "<< G2 << " D2= " << D2 <<   endl;
      //cout << E << " " << Y << endl;
    }

  double c = 1./sumY; // normalization const for rho
  //cout << "c= " << c << endl;
  
  double gamma1Tot = 0.;
  double S1Tot = 0.;
  double B1Tot = 0.;
  double delta1Tot = 0.;
  double yy_old = 0.;
  double E1,E2;
  for(int i = 0; i < N; i++)
    {

      double E = double(i) * dE;

      double P = P_l(E, ac123, l1, mu123, z1, z2 + z3);

      double g = Gamma(E, Et, ac23, l2, mu23, z2, z3, rwidth2_23);
      double delta = -rwidth2_23*(Shift(Et-E,ac23,l2,mu23,z2,z3)- B23);
      double yy = g / (pow(Et - E - Er - delta,2) + pow(g,2)/4);
      if ((yy-ymax/2.)*(yy_old-ymax/2.) < 0. && yy > yy_old)
	{
          double slope = (yy-yy_old)/dE;
          E1 = E - dE +(ymax/2.-yy_old)/slope;
	}
      else if ((yy-ymax/2.)*(yy_old-ymax/2.) < 0. && yy < yy_old)
	{
	  double slope = (yy-yy_old)/dE;
          E2 = E - dE + (ymax/2.-yy_old)/slope;
	}
      yy_old = yy;


      double y = 2.*c*P * yy;

      /*
      if (fabs(Et- 2.074) < 0.001)
	cout << Et << " " << E << " " << Et - E << " " << P << " " << 
	  c * g / (pow(Et - E - Er - delta,2) + pow(g,2)/4) << " " << y  << endl;
      */
      gamma1Tot += y*dE;
      
      double S1 = Shift(E,ac123, l1, mu123, z1, z2 + z3);
      double B1 = Shift(2.074-Et+E,ac123, l1, mu123, z1, z2 + z3);
      delta1Tot += -4.*(S1-B1)*c*g*dE/ (pow(Et - E - Er - delta,2) 
                          + pow(g,2)/4);

      /*
  if (fabs(Et- 2.074) < 0.001)
    cout << E << " " << 2.074-Et+E << " " << S1 << " " << B1 << " " << delta1Tot << endl; 
      */
      S1Tot += S1*c*g*dE/ (pow(Et - E - Er - delta,2) + pow(g,2)/4);
      B1Tot += B1*c*g*dE/ (pow(Et - E - Er - delta,2) + pow(g,2)/4);

      /*
           if (i%10==0)
	     cout << E << " " << Et-E << " " << P << " " 
                  << g / (pow(Et - E - Er -  delta,2) + pow(g,2)/4) 
	          << " " << y << " " << delta << endl;
      */
      profile[i] = ((i == 0)?0:profile[i-1]) + y*dE;

    }
  //cout << "Et= " << Et << " gamma1Tot= " << S1Tot << endl;
  //cout << Et << "  " << gamma1Tot << " " << delta1Tot << " " << S1Tot <<  " " << B1Tot << endl;
  //cout << Et << " " << gamma1Tot/(pow(Et-2.016-delta1Tot,2)+pow(gamma1Tot,2)) << endl;

  cout << "Emax = " << Emax << endl;
  cout << E1 << " " << E2 << " " << "FWHM = " << E2-E1 << endl;


  double scale = 1 / profile[N-1];

  for(int i = 0; i < N; i++)
    {
      profile[i] *= scale;
    }
}

Profile::~Profile()
{
  delete[] profile;
}

double Profile::rand(double xran)
{

  int n = -1;
  for(int i = 0; i < N; i++)
    {
      if(xran <= profile[i])
	{
	  n = ((i==0)?0:i-1);
	  break;
	}
    }
  if(n == -1) return -1;

  double deltaP = xran - profile[n];
  double dEdP = dE / (profile[n+1] - profile[n]);

  return n*dE + dEdP * deltaP;
}

//profile of a single decay
Profile::Profile(double Er, int z1, int z2, double mu, double ac, int l,
		 double rwidth2, double B0, double range)
{
  B = B0;

  if (B == 999.) 
    {
     B= Shift(Er,ac,l,mu,z1,z2);

    }
  cout << "Er= " << Er << " B= " << B <<endl;
  gsl_set_error_handler_off();
  dE = .03;
  if (Er/dE < 400) dE = Er/400.;
  if (range== 0)N = (int)floor(Er*2.155/dE);
  else N = (int)floor(range/dE);
  profile = new double[N];

  cout << "N = " << N << endl;

  Earray = new double[N];
  yarray = new double[N];
  Parray = new double[N];

  double P = P_l(Er, ac, l, mu, z1, z2);


  double Emax = 0;

  double g = 2.*rwidth2*P;
  cout << " Gamma at Er = " << g << " P = " << P << endl;

  double ymax = 0.;
  for(int i = 1; i <= N; i++)
    {

      double E = (double(i)) * dE;

      double P = P_l(E, ac, l, mu, z1, z2);
      Parray[i-1] = P;
      double g = 2.*rwidth2*P;
      double delta = -rwidth2*(Shift(E,ac,l,mu,z1,z2)- B);
      double y =  g / (pow(E - Er - delta,2) + pow(g,2)/4);
      
      //if (i%100 == 0)cout << E << " " << y << " " << delta << " " << g << " " << P << endl;

  

      if (i == 1) 
	{
	  cout << g << " " << E << " " << Er << " " << delta << " " << y << " " << B <<  endl;
	}

      Earray[i-1] = E;
      yarray[i-1] = y;
      profile[i-1] = y;
      if (y > ymax)
	{
           ymax = y;     
           Emax = E;
	}
      
    }

  double yy_old = 0.;
  double E1,E2;
  for (int i=1;i<=N;i++)
    {
      double E = Earray[i];
      double yy = profile[i];
      if ((yy-ymax/2.)*(yy_old-ymax/2.) < 0. && yy > yy_old)
	{
          double slope = (yy-yy_old)/dE;
          E1 = E - dE +(ymax/2.-yy_old)/slope;
	}
      else if ((yy-ymax/2.)*(yy_old-ymax/2.) < 0. && yy < yy_old)
	{
	  double slope = (yy-yy_old)/dE;
          E2 = E - dE + (ymax/2.-yy_old)/slope;
	}
      yy_old = yy;

      if(i != 1)profile[i-1] += profile[i-2];
    }

  cout << "Emax = " << Emax << endl;
  cout << E1 << " " << E2 << " " << "FWHM = " << E2-E1 << endl;
  //abort();
  double scale = 1. / profile[N-1];


  for(int i = 1; i <= N; i++)
    {
      profile[i-1] *= scale;
    }

}



//profile of a single decay with two branches
Profile::Profile(double Er, double DE0, int z1, int z2, double mu, double ac, int l1, int l2, double rwidth2_1, double rwidth2_2, double B_1, double B_2)
{

  DE = DE0;
  if (B_1 == 999.) 
    {
     B_1= Shift(Er,ac,l1,mu,z1,z2);
     cout << "B_1= " << B_1 <<endl;
    }
  if (B_2 == 999.)
    {
     B_2= Shift(Er-DE,ac,l2,mu,z1,z2);
     cout << "B_2= " << B_2 <<endl;
    }
  gsl_set_error_handler_off();
  dE = .03;
  
  N = (int)floor(Er*2.55/dE);

  cout << "N = " << N << endl;

  Earray = new double[N];
  y1array = new double[N];
  y2array = new double[N];
  profile = new double [N];
  //P1array = new double[N];
  //P2array = new double[N];

  double P1 = P_l(Er, ac, l1, mu, z1, z2);
  double P2 = P_l(Er-DE,ac,l2,mu,z1,z2);

  double Emax = 0;

  double g1 = 2.*rwidth2_1*P1;
  double g2 = 2.*rwidth2_2*P2;
  cout << " Gamma1 at Er = " << g1 << " P1 = " << P1 << endl;
  cout << " Gamma2 at Er = " << g2 << " P2 = " << P2 << endl;

  double ymax = 0.;
  for(int i = 1; i <= N; i++)
    {

      double E = (double(i)) * dE;

      double P1 = P_l(E, ac, l1, mu, z1, z2);
      double P2 = 0.;
      if (E-DE > 0) P2 = P_l(E-DE,ac,l2,mu,z1,z2);
      //P1array[i-1] = P1;
      //P2array[i-1] = P2;

      double g1 = 2.*rwidth2_1*P1;
      double g2 = 2.*rwidth2_2*P2;
      double delta1 = -rwidth2_1*(Shift(E,ac,l1,mu,z1,z2)- B_1);
      double delta2 = 0;
      if (E-DE > 0.) delta2 = -rwidth2_2*(Shift(E-DE,ac,l2,mu,z1,z2)- B_2);
      else  delta2 = -rwidth2_2*(Shift(.001,ac,l2,mu,z1,z2)- B_2);
      double y1 =  g1 / (pow(E - Er - delta1 - delta2,2) + pow((g1+g2),2)/4);
      double y2 =  g2 / (pow(E - Er - delta1 - delta2,2) + pow((g1+g2),2)/4);
  

      Earray[i-1] = E;
      y1array[i-1] = y1;
      y2array[i-1] = y2;
      profile[i-1] = y1 + y2;
      if (y1 > ymax) ymax = y1;     

      
    }

  double yy_old = 0.;
  double E1,E2;
  for (int i=1;i<=N;i++)
    {
      double E = Earray[i];
      double yy = profile[i];
      if ((yy-ymax/2.)*(yy_old-ymax/2.) < 0. && yy > yy_old)
	{
          double slope = (yy-yy_old)/dE;
          E1 = E - dE +(ymax/2.-yy_old)/slope;
	}
      else if ((yy-ymax/2.)*(yy_old-ymax/2.) < 0. && yy < yy_old)
	{
	  double slope = (yy-yy_old)/dE;
          E2 = E - dE + (ymax/2.-yy_old)/slope;
	}
      yy_old = yy;

      if(i != 1)profile[i-1] += profile[i-2];
    }

  cout << E1 << " " << E2 << " " << "FWHM = " << E2-E1 << endl;
  //abort();
  double scale = 1. / profile[N-1];


  for(int i = 1; i <= N; i++)
    {
      profile[i-1] *= scale;
    }

}

double Profile::rand_2branches(double xran1, double xran2)
{

  int n = -1;
  for(int i = 0; i < N; i++)
    {
      if(xran1 <= profile[i])
	{
	  n = ((i==0)?0:i-1);
	  break;
	}
    }
  if(n == -1) return -1;

  double deltaP = xran1 - profile[n];
  double dEdP = dE / (profile[n+1] - profile[n]);

  double dEE = dEdP * deltaP;
  double E = n*dE + dEE;

  // branching ratio
  double branch1 = (dE-dEE)/DE*y1array[n] + dEE/dE*y1array[n+1];
  double branch2 = (dE-dEE)/DE*y2array[n] + dEE/dE*y2array[n+1];
  double prob1 = branch1/(branch1+branch2);
  if (xran2< prob1) return E;
  else return (E-DE);
 

}



double Profile::getWigner(double  mu, double ac)
{
  return 41.863/mu/pow(ac,2);
}
