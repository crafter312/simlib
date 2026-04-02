#ifndef _frag
#define _frag

#include "array.h"
#include "loss.h"
#include <string>
#include "mScat.h"
#include "random.h"
#include "frame.h"
#include "fragment.h"
#include "pixels.h"
#include "tele.h"
using namespace std;

/**
 *!\brief information about a single fragment and its interection with the detector
 *
 */
class CFrag : public fragment
{
 public:
  static CRandom ran;
  static float const pi; //!< 3.14159....

  float CsI_res;   

  CFrag(float,float,string,float,float,float=1.,bool=false);
  ~CFrag();
  CArray * Array;
  CArray * shadowArray;
  int hit();
  int hit2(float,float);
  int hit(float,float);
  int hitShadow(float,float);
  int is_hit;
  int hitTele;

  void SetVelocity(float*,float*);
  void AddVelocity(double*);
  void AddVelocityRel(float*);
  float Eloss(float); // energy loss in target
  float Egain(float); // corrects for energy loss in target
  void MultiScat(float);
  void findVectors(float*,float*,float*,float*,float*,float,float,float,float,float,float,float,float,float);

  bool targetInteraction(float,float);

    CLoss *loss;
    CMultScat *multScat;

    CFrame *real;      //<!real particles energy, velocity, etc
    CFrame *recon;    //<!reconstructed properties

    pixels * Pixels;
    CTele * Teles[14];
    float scaleSmallAngle;

    float p_thres[56];

    bool alphaHole(int itele, int ifront, int iback);
    bool protonHole(int itele, int ifront, int iback);
    bool carbonHole(int itele, int ifront, int iback);
    bool nitrogenHole(int itele, int ifront, int iback);
    bool useRealP;   //<! use the real momentum, instead of detected value
    bool extra;
};

#endif
