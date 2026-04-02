#if !defined(MYLIB_CONSTANTS_H)
#define MYLIB_CONSTANTS_H 1

#include <cmath>

const double m0 = 931.478;  // atomic mass unit
const double c = 30.;   //speed of light
const double vfact = c/sqrt(m0);  //velocity (cm/ns) = vfact *sqrt(2.*E(MeV)/A(amu))
const double pi = acos(-1.);
//masses from AME2016 compilation
const double mass_p = 7.28897;
const double mass_d = 13.13572;
const double mass_t = 14.9498;
const double mass_3He = 14.93121;
const double mass_alpha = 2.42491;
const double mass_6He = 17.5920;
const double mass_8He = 31.6096;
const double mass_5Li = 11.678886;
const double mass_6Li = 14.0868;
const double mass_7Li = 14.9071;
const double mass_8Li = 20.9458;
const double mass_9Li = 24.9549;
const double mass_6Be = 18.375033;
const double mass_7Be = 15.768999;
const double mass_8Be = 4.9416;
const double mass_9Be = 11.3484;
const double mass_10Be = 12.6074;
const double mass_11Be = 20.1771;
const double mass_8B = 22.9215;
const double mass_9B = 12.416488;
const double mass_10B = 12.0506;
const double mass_11B = 8.6677;
const double mass_9C = 28.910972;
const double mass_10C = 15.698672;
const double mass_11C = 10.649396;
const double mass_12C = 0.;
const double mass_13C = 3.12500888;
const double mass_14C = 3.019892;
const double mass_11N = 24.303559;
const double mass_12N = 17.338068;
const double mass_13N = 5.345481;
const double mass_14N = 2.863416;
const double mass_15N = 0.101438;
const double mass_13O = 23.115432;
const double mass_14O = 8.007781;
const double mass_15O = 2.855605;
const double mass_16O = -4.737001;
const double mass_17O = -0.808763;
const double mass_14F = 31.964402;
const double mass_15F = 16.566751;
const double mass_17F = 1.951702;
const double mass_18F = .873113;
const double mass_17Ne = 16.500447;
const double mass_18Ne = 5.317614;

#endif
