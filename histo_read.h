#ifndef histo_read_
#define histo_read_
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "TH1F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TH3I.h"
#include "TFile.h"
#include "TDirectory.h"

using namespace std;

class histo_read
{
protected:

  TFile * file_read; //!< output root file

  //correlations
  TDirectoryFile * dirCorr; //!< directory for correlations
  TDirectory * dir11O; //!< directory for 11O correlations
  TDirectory * dir12O; //!< directory for 12O correlations
  TDirectory * dir13O; //!< directory for 13O correlations
  TDirectory * dir11N; //!< directory for 11N correlations
  TDirectory * dir12N; //!< directory for 12N correlations
  TDirectory * dir6Li; //!< directory for 6Li correlations
  TDirectory * dir7Li; //!< direcotry for 7Li correlations 
  TDirectory * dir7Be; //!< direcotry for 7Be correlations 
  TDirectory * dir12C;

  TDirectory * dir11O_CsI;
  TDirectory * dir11O_HitMaps;
  TDirectory * dir11O_stripcuts;
  TDirectory * dir12O_stripcuts;
  TDirectory * dir13O_stripcuts;

  TDirectory * dir12O_stripcuts_trans;
  TDirectory * dir12O_smallCuts;


  TDirectoryFile *dirdEE; //!< directory for the dEE maps

  TDirectoryFile *dirProjections;
  TDirectoryFile *dirProjectionsRaw;
  



  TDirectoryFile *dirCsI;
  TDirectory *dirCsIGate;

  TDirectoryFile *dirSi;
  TDirectory *dirELoss;

  //Front & Back
  TDirectoryFile * dirFB; //!< directory for front/back mixed spectra


    
public:
    histo_read();                  //!< constructor
    ~histo_read(){};
    void write(); //!< write the root spectra to file


    int Ntele;
    int Nstrip;
    int Nceasar;
    int NCsI;
    int Nring;

    TH2I ** dEE;
    TH2I * HitMap;

    TH1I *** dEE_Projections;
    TH1I *** dEE_ProjectionsRaw;

    TH1I * dEE_CsIReactions_25;
    TH1I * dEE_CsIReactions_26;
    TH1I * dEE_CsIReactions_28;
    TH1I * dEE_CsIReactions_31;

    //CsI spectra
    TH1I ** Light;

    //ELoss spectra
    TH1I ** ELoss;
      
    //Erel for 11O by CsI spectra
    TH1I ** Erel_11O_2p9C_CsI;

    TH2I ** FvB_Tele;
    
    TH1I* TEC_13O_added;

    //FB
    TH2I ** FBMult;
    TH1I ** FBDiff;
    TH1I ** FBDiffLG;
    TH2I ** FB;

   
    
    //correlations
    //11O
    TH1I* Erel_11O_2p9C;
    TH1I* Erel_11O_2p9C_maskp0;
    TH1I* Erel_11O_2p9C_maskp1;
    TH1I* Erel_11O_2p9C_CsI25;
    TH1I* Erel_11O_2p9C_CsI26;
    TH1I* Erel_11O_2p9C_CsI28;
    TH1I* Erel_11O_2p9C_CsI31;
    TH1I* Erel_11O_2p9C_nonCentral;
    TH1I* Erel_11O_2p9C_Trans;
    TH2I* Erel_11O_2p9C_icsi;
    TH1F* Erel_11O_2p9C_Fake;
    TH1F* Erel_11O_2p9C_Fake_Channeling;
    TH1F* Erel_11O_2p9C_Fake_CsIreaction;
    TH1F* Erel_11O_2p9C_Fake_CsIreaction_C11;
    TH1F* Erel_11O_2p9C_Fake_CsIreaction_Trans;
    TH1F* Erel_11O_2p9C_Fake_CsIreaction_C11_Trans;
    TH1I* vel_11O_2p9C;
    TH2I * JacobiY_11O_2p9C;
    TH2I * JacobiT_11O_2p9C;
    TH2I * JacobiY_11O_2p9C_Fake_CsIreaction;
    TH2I * JacobiT_11O_2p9C_Fake_CsIreaction;
    TH2I * JacobiY_11O_2p9C_Fake_CsIreaction_C11;
    TH2I * JacobiT_11O_2p9C_Fake_CsIreaction_C11;

    TH1I* Erel_11O_2p9C_stripcuts;
    TH1I* Erel_11O_2p9C_Trans_stripcuts;
    TH1F* Erel_11O_2p9C_Fake_CsIreaction_stripcuts;
    TH1F* Erel_11O_2p9C_Fake_CsIreaction_C11_stripcuts;
    TH1F* Erel_11O_2p9C_Fake_CsIreaction_Trans_stripcuts;
    TH1F* Erel_11O_2p9C_Fake_CsIreaction_C11_Trans_stripcuts;
    TH2I * JacobiY_11O_2p9C_stripcuts;
    TH2I * JacobiT_11O_2p9C_stripcuts;
    TH2I * JacobiY_11O_2p9C_Fake_CsIreaction_stripcuts;
    TH2I * JacobiT_11O_2p9C_Fake_CsIreaction_stripcuts;
    TH2I * JacobiY_11O_2p9C_Fake_CsIreaction_C11_stripcuts;
    TH2I * JacobiT_11O_2p9C_Fake_CsIreaction_C11_stripcuts;

    TH2I * Erel_11O_total_vs_mask;

    TH2I * C10HitMap;

    TH2I * Erel_EnergyBleed;
    TH1I * EnergyBleed_Projections;

    TH1I* Erel_11O_2p9C_minusMiddleStrips;


    //13O
    TH1I* Erel_13O_2p11C;
    TH1I* Erel_13O_2p11C_maskp0;
    TH1I* Erel_13O_2p11C_maskp1;
    TH1I* Erel_13O_2p11C_Trans;
    TH1I* vel_13O_2p11C;
    TH2I * JacobiY_13O_2p11C;
    TH2I * JacobiT_13O_2p11C;

    TH1I* Erel_13O_2p11C_stripcuts;
    TH1I* Erel_13O_2p11C_Trans_stripcuts;
    TH2I * JacobiY_13O_2p11C_stripcuts;
    TH2I * JacobiT_13O_2p11C_stripcuts;

    TH2I * Erel_13O_total_vs_mask;


    //12O
    TH1I* Erel_12O_2p10C;
    TH1I* Erel_12O_2p10C_Trans;
    TH1I* vel_12O_2p10C;
    TH2I * JacobiY_12O_2p10C;
    TH2I * JacobiT_12O_2p10C;
    TH2I * JacobiY_12O_2p10C_star1;
    TH2I * JacobiT_12O_2p10C_star1;    
    TH2I * JacobiY_12O_2p10C_star2;
    TH2I * JacobiT_12O_2p10C_star2;  
    TH1F * Erel_12O_2p10C_Fake_CsIreaction;
    TH1F * Erel_12O_2p10C_Fake_CsIreaction_Trans; 
    TH1I* Erel_12O_4p2a;

    TH1I* Erel_12O_2p10C_stripcuts;
    TH1I* Erel_12O_2p10C_Trans_stripcuts;
    TH2I * JacobiY_12O_2p10C_stripcuts;
    TH2I * JacobiT_12O_2p10C_stripcuts;
    TH2I * JacobiY_12O_2p10C_star1_stripcuts;
    TH2I * JacobiT_12O_2p10C_star1_stripcuts;    
    TH2I * JacobiY_12O_2p10C_star2_stripcuts;
    TH2I * JacobiT_12O_2p10C_star2_stripcuts;  

    TH2I * JacobiY_12O_2p10C_Trans_stripcuts;
    TH2I * JacobiT_12O_2p10C_Trans_stripcuts;
    TH2I * JacobiY_12O_2p10C_Trans_star1_stripcuts;
    TH2I * JacobiT_12O_2p10C_Trans_star1_stripcuts;    
    TH2I * JacobiY_12O_2p10C_Trans_star2_stripcuts;
    TH2I * JacobiT_12O_2p10C_Trans_star2_stripcuts;  

    TH2I * JacobiY_Cut1;
    TH2I * JacobiY_Cut2;
    TH2I * JacobiY_Cut3;
    TH2I * JacobiY_Cut4;
    TH2I * JacobiY_Cut5;
    TH2I * JacobiY_Cut6;


    TH1F * Erel_12O_2p10C_Fake_CsIreaction_stripcuts;
    TH1F * Erel_12O_2p10C_Fake_CsIreaction_Trans_stripcuts; 
    TH2I * JacobiY_12O_2p10C_Fake_stripcuts;
    TH2I * JacobiT_12O_2p10C_Fake_stripcuts;
    TH2I * JacobiY_12O_2p10C_star1_Fake_stripcuts;
    TH2I * JacobiT_12O_2p10C_star1_Fake_stripcuts;    
    TH2I * JacobiY_12O_2p10C_star2_Fake_stripcuts;
    TH2I * JacobiT_12O_2p10C_star2_Fake_stripcuts;  

    TH1I * JacobiY_12O_2p10C_1D_stripcuts;
    TH1I * JacobiY_12O_2p10C_star1_1D_stripcuts;
    TH1I * JacobiY_12O_2p10C_star2_1D_stripcuts;

    //11N
    TH1I * Erel_11N_3p2a;
    TH1I * vel_11N_3p2a;
    TH1I* Ex_9B_11N_3p2a;
    TH1I * Erel_11N_3p2a_IAS;
    TH1I * Ex_8Be_11N_3p2a;
    TH2I * B9_vs_Total;

    //12N
    TH1I* Erel_12N_2p10B;
    TH1I* Erel_12N_2p10B_Trans;
    TH1I* vel_12N_2p10B;
    TH1I* angle_12N_2p10B;
    TH2I * JacobiY_12N_2p10B;
    TH2I * JacobiT_12N_2p10B;
    TH1I* Erel_12N_1p11C;
    TH1I* Erel_12N_1p11C_Trans;
    TH1I* Erel_12N_1p11C_AngCut;
    TH1I* Erel_12N_1p11C_Trans_AngCut;
    TH1I* vel_12N_1p11C;
    TH1I* angle_12N_1p11C;

     //Blocker CsI
    TH2D* Blocker_ETOF;

    //7Li
    TH1I * Erel_7Li;

    TH1I * Erel_Back_Diff;
    TH1I * Erel_For_Diff;
    TH1I * Erel_Trans_Diff;
    TH1I * Erel_Back_Same;
    TH1I * Erel_For_Same;
    TH1I * Erel_Trans_Same;
    TH1I * Costheta;
    TH1I * vel_7Li;

    //6Li

    TH1I * Erel_6Li;
    TH1I * Ex_6Li_da;
    
    //Be7
    TH1I * Ex7Be_alpha3He;
    TH1I * vel_7Be_alpha3He;
    TH1I * vel_7Be_alpha3He_left;
    TH1I * vel_7Be_alpha3He_right;
    TH1I * TEC_7Be_alpha3He;
    TH1I * Ex7Be_alpha3He_gamma;
    TH1I * Ex7Be_alpha3He_gammaback1;
    TH1I * Ex7Be_alpha3He_gammaback2;
    TH1I * Ex7Be_alpha3He_gammaback3;
    TH2I * Caesarmult;

    //C12
    TH1I * Ex12C_da6Li;
};
#endif
