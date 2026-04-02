#ifndef histo_sort_
#define histo_sort_
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "TH1I.h"
#include "TH2I.h"
#include "TH3I.h"
#include "TFile.h"
#include "TDirectory.h"

using namespace std;

class histo_sort
{
 protected:

  
  TFile * file; //!< output root file

  //correlations
  //de
  TDirectoryFile * dirdee_sort; //!< directory for dee maps
  

  //Summary
  TDirectoryFile * dirSum; //!< directory for summary spectra

  //Front
  TDirectoryFile * dirFront; //!< directory for front strip info
  TDirectory * dirFrontRaw; //!< directory for raw front strip info
  TDirectory * dirFrontV;   //!< directory for front voltages
  TDirectory * dirFrontLV; //!< directory for low gain front voltages
  TDirectory * dirFrontCal; //!< directory for calibrated front strip info
  TDirectory * dirFrontLGR; //!< directory for raw front string low gain
  TDirectory * dirFrontLGC; //!< directory for calibrated front string low gain
  TDirectory * dirFrontTime; //!< directory for times
  TDirectory * dirFrontTLG; //!< directory for low gains times
  TDirectory * dirFrontHiLo;
  TDirectory * dirTotCal;

  //Back
  TDirectoryFile * dirBack; //!< directory for back strip info
  TDirectory * dirBackRaw; //!< directory for raw back strip info
  TDirectory * dirBackV;   //!directory for back voltages
  TDirectory * dirBackLV; //!directory for back low gain voltages
  TDirectory * dirBackCal; //!< directory for calibrated back strip info
  TDirectory * dirBackLGR; //!< directory for raw back string low gain
  TDirectory * dirBackLGC; //!< directory for calibrated back string low gain
  TDirectory * dirBackLGCC; //!< directory for calibrated back string low gain
  TDirectory * dirBackTime; //!< directory for times
  TDirectory * dirBackTLG; //!<directory for low gain times
  TDirectory * dirBackHiLo;

  //CsI
  TDirectoryFile *dirCsI; //!< directory for the CsI info
  TDirectory * dirCsIRaw; //!< CsI energies
  TDirectory * dirCsITime; //!< CsI time
  TDirectory * dirCsISum; //!< CsI Sum
  TDirectory * dirCsIBlocker; //!< Blocker CsI
  TDirectory * dirCsIGate; //!< CsI gated on Z
  TDirectory * dirCsIMult; //!< CsI Multiplicity




 public:
  histo_sort();                  //!< constructor
  ~histo_sort(){};
  void write(); //!< write the root spectra to file
  
  
  int Ntele;
  int Nstrip;
  int Nceasar;
  int NCsI;
  int Nring;
  
  
  
  TH1I * Txfp;
  TH1I * TRF;
  
  TH2I ** EFCSum;
  TH2I ** EFCLSum;
  TH2I ** EBCSum;
  TH2I ** EFSum;
  TH2I ** EBSum;
  
  TH2I ** EFTSum;
  TH2I ** EBTSum;

  TH2I ** EFHiLo;
  TH2I ** EBHiLo;
  
  TH1I * SiFTime;
  TH1I * SiBTime;
  
  //Front spectra
  TH1I *** EfrontR; //high gain raw
  TH1I *** EfrontV; //high gain voltages
  TH1I *** EfrontLV; //low gain voltages
  TH1I *** EfrontC; //high gain calib
  TH1I *** EfrontLR; //low gain raw
  TH1I *** EfrontLC; // low gain calib
  TH1I *** TfrontR; // time raw

  //back spectra
  TH1I *** EbackR;
  TH1I *** EbackV;
  TH1I *** EbackLV; 
  TH1I *** EbackC;
  TH1I *** EbackLR;
  TH1I *** EbackLC;
  TH1I *** TbackR;
  TH1I *** TbackLG;

  TH1I ** ECsI;
  TH1I ** TCsI;
  TH2I * ECsISum;
  TH2I * ECsICSum;
  TH1I ** CsIOver;
  TH1I ** ECsI_Zgate;
  TH1I * CsIMult;
  TH1I * Etot;



    
  TH2I ** CaldEE;


  //Blocker CsI
  TH1I* Blocker_E1;
  TH1I* Blocker_E2;
  TH1I* Blocker_E3;
  TH1I* Blocker_E4;
  TH1I* Blocker_ESum;
  TH2I* Blocker_ETOF;


  TH1I * T_RFSCIN;
  TH1I * T_RFCYC;
  TH1I * T_A1900;
    
  TH2I **  dee_sort;
  TH2I * HL_24;
    
};
#endif
