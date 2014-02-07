/*
author: Andreas Ritter
created: 04/12/2007
last edited: 05/05/2007
compiler: g++ 4.0
basis machine: Ubuntu Linux 6.06
*/

#include "MExtract.h"

using namespace std;

int main(int argc, char *argv[])
{
  cout << "MExtract::main: argc = " << argc << endl;
  if (argc < 6)
  {
    cout << "MExtract::main: ERROR: Not enough parameters specified!" << endl;
    cout << "USAGE: optextract(char[] FitsFileName_In, char[] DatabaseFileName_In, char[] FitsFileName_Out, double Gain, double ReadOutNoise[, TELLURIC=int[0 - none, 1 - Piskunov, 2 - LinFit]][, MAX_ITER_SF=int][, MAX_ITER_SKY=int][, MAX_ITER_SIG=int][, SWATH_WIDTH=int][, ERR_IN=char[]][, ERR_OUT_2D=char[]][, ERR_OUT_EC=char[]][, SKY_OUT_EC=char[]][, SKY_OUT_2D=char[]][, SKY_ERR_OUT_EC=char[]][, PROFILE_OUT=char[]][, IM_REC_OUT=char[]][, REC_FIT_OUT=char[]][, MASK_OUT=char[]][, SPFIT_OUT_EC=char[]][, EC_FROM_PROFILE_OUT=char[]][, ERR_FROM_PROFILE_OUT=char[]])" << endl;
    exit(EXIT_FAILURE);
  }

  Array<CString, 1> CS_A1_Args(10);
  CS_A1_Args = CString("\0");
  void **PP_Args;
  PP_Args = (void**)malloc(sizeof(void*) * 10);
  CString CS(" ");
  CString CS_comp(" ");
  CString *P_CS = new CString(" ");
  CString *P_CS_ErrIn = new CString(" ");
  CString *P_CS_ErrOut = new CString(" ");
  CString *P_CS_ErrOutEc = new CString(" ");
  CString *P_CS_SkyOut = new CString(" ");
  CString *P_CS_SkyArrOut = new CString(" ");
  CString *P_CS_SkyErrOut = new CString(" ");
  CString *P_CS_ImOut = new CString(" ");
  CString *P_CS_RecFitOut = new CString(" ");
  CString *P_CS_ProfileOut = new CString(" ");
  CString *P_CS_MaskOut = new CString(" ");
  CString *P_CS_SPFitOut = new CString(" ");
  CString *P_CS_EcFromProfileOut = new CString(" ");
  CString *P_CS_ErrFromProfileOut = new CString(" ");

  int I_SwathWidth = 0;
  int I_MaxIterSF = 8;
  int I_MaxIterSky = 12;
  int I_MaxIterSig = 2;
  char *P_CharArr_In = (char*)argv[1];
  char *P_CharArr_DB = (char*)argv[2];
  char *P_CharArr_Out = (char*)argv[3];
  double D_Gain = (double)(atof((char*)argv[4]));
  cout << "MExtract::main: D_Gain set to " << D_Gain << endl;
  double D_ReadOutNoise = (double)(atof((char*)argv[5]));
  cout << "MExtract::main: D_ReadOutNoise set to " << D_ReadOutNoise << endl;
  int I_Telluric=0;

  /// read optional parameters
  for (int i = 6; i <= argc; i++){
    CS.Set((char*)argv[i]);
    cout << "MExtract: Reading Parameter " << CS << endl;
    CS_comp.Set("TELLURIC");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS);
        P_CS = CS.SubString(CS_comp.GetLength()+1);
        I_Telluric = (int)(atoi(P_CS->Get()));
        cout << "MExtract::main: I_Telluric set to " << I_Telluric << endl;
      }
    }

    CS_comp.Set("SWATH_WIDTH");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS);
        P_CS = CS.SubString(CS_comp.GetLength()+1);
        I_SwathWidth = (int)(atoi(P_CS->Get()));
        cout << "MExtract::main: I_SwathWidth set to " << I_SwathWidth << endl;
      }
    }

    CS_comp.Set("MAX_ITER_SF");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS);
        P_CS = CS.SubString(CS_comp.GetLength()+1);
        I_MaxIterSF = (int)(atoi(P_CS->Get()));
        cout << "MExtract::main: I_MaxIterSF set to " << I_MaxIterSF << endl;
      }
    }

    CS_comp.Set("MAX_ITER_SKY");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS);
        P_CS = CS.SubString(CS_comp.GetLength()+1);
        I_MaxIterSky = (int)(atoi(P_CS->Get()));
        cout << "MExtract::main: I_MaxIterSky set to " << I_MaxIterSky << endl;
      }
    }

    CS_comp.Set("MAX_ITER_SIG");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS);
        P_CS = CS.SubString(CS_comp.GetLength()+1);
        I_MaxIterSig = (int)(atoi(P_CS->Get()));
        cout << "MExtract::main: I_MaxIterSig set to " << I_MaxIterSig << endl;
      }
    }

    /// 2D
    CS_comp.Set("ERR_IN");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_ErrIn);
        P_CS_ErrIn = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: ERR_IN set to " << *P_CS_ErrIn << endl;
      }
    }

    /// 2D
    CS_comp.Set("ERR_OUT_2D");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_ErrOut);
        P_CS_ErrOut = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: ERR_OUT_2D set to " << *P_CS_ErrOut << endl;
      }
    }

    /// 1D
    CS_comp.Set("ERR_OUT_EC");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_ErrOutEc);
        P_CS_ErrOutEc = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: ERR_OUT_EC set to " << *P_CS_ErrOut << endl;
      }
    }

    /// 1D
    CS_comp.Set("SKY_OUT_EC");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_SkyOut);
        P_CS_SkyOut = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: SKY_OUT_EC set to " << *P_CS_SkyOut << endl;
      }
    }

    /// 2D
    CS_comp.Set("SKY_OUT_2D");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_SkyArrOut);
        P_CS_SkyArrOut = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: SKY_OUT_2D set to " << *P_CS_SkyArrOut << endl;
      }
    }

    /// 1D
    CS_comp.Set("SKY_ERR_OUT_EC");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_SkyErrOut);
        P_CS_SkyErrOut = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: SKY_ERR_OUT_EC set to " << *P_CS_SkyErrOut << endl;
      }
    }

    /// 2D
    CS_comp.Set("IM_REC_OUT");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_ImOut);
        P_CS_ImOut = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: IM_REC_OUT set to " << *P_CS_ImOut << endl;
      }
    }

    /// 2D
    CS_comp.Set("REC_FIT_OUT");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_RecFitOut);
        P_CS_RecFitOut = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: REC_FIT_OUT set to " << *P_CS_RecFitOut << endl;
      }
    }

    /// 2D
    CS_comp.Set("PROFILE_OUT");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_ProfileOut);
        P_CS_ProfileOut = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: PROFILE_OUT set to " << *P_CS_ProfileOut << endl;
      }
    }

    /// 2D
    CS_comp.Set("MASK_OUT");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_MaskOut);
        P_CS_MaskOut = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: MASK_OUT set to " << *P_CS_MaskOut << endl;
      }
    }

    /// 1D
    CS_comp.Set("SPFIT_OUT_EC");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_SPFitOut);
        P_CS_SPFitOut = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: SPFIT_OUT_EC set to " << *P_CS_SPFitOut << endl;
      }
    }

    /// 1D
    CS_comp.Set("EC_FROM_PROFILE_OUT");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_EcFromProfileOut);
        P_CS_EcFromProfileOut = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: EC_FROM_PROFILE_OUT set to " << *P_CS_EcFromProfileOut << endl;
      }
    }

    /// 1D
    CS_comp.Set("ERR_FROM_PROFILE_OUT");
    if (CS.GetLength() > CS_comp.GetLength()){
      delete(P_CS);
      P_CS = CS.SubString(0,CS.CharPos('=')-1);
      cout << "MExtract::main: *P_CS set to " << *P_CS << endl;
      if (P_CS->EqualValue(CS_comp)){
        delete(P_CS_ErrFromProfileOut);
        P_CS_ErrFromProfileOut = CS.SubString(CS_comp.GetLength()+1);
        cout << "MExtract::main: ERR_FROM_PROFILE_OUT set to " << *P_CS_ErrFromProfileOut << endl;
      }
    }
  }

//  return false;

  if (I_Telluric > 0){
    CS_A1_Args(1).Set("TELLURIC");
    PP_Args[1] = &I_Telluric;
  }
  time_t seconds;
//  if (argc == 8)
//  {
//    I_SwathWidth = (int)(atoi((char*)argv[7]));
//    cout << "MExtract::main: I_SwathWidth set to " << I_SwathWidth << endl;
  if (I_SwathWidth > 0.){
    CS_A1_Args(0).Set("SWATH_WIDTH");
    PP_Args[0] = &I_SwathWidth;
  }
  else
  {
    CS_A1_Args(0) = CString("");
    PP_Args[0] = &I_SwathWidth;
  }
  CString CS_FitsFileName_In;
  CS_FitsFileName_In.Set(P_CharArr_In);
  CString CS_FitsFileName_Out;
  CS_FitsFileName_Out.Set(P_CharArr_Out);
  CString CS_DatabaseFileName_In;
  CS_DatabaseFileName_In.Set(P_CharArr_DB);

  CFits F_Image;
  cout << "MExtract::main: Starting F_Image.SetFileName(" << CS_FitsFileName_In.Get() << ")" << endl;
  if (!F_Image.SetFileName(CS_FitsFileName_In))
  {
    cout << "MExtract::main: ERROR: F_Image.SetFileName(" << CS_FitsFileName_In.Get() << ") returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// Set ReadOutNoise
  cout << "MExtract::main: Starting F_Image.Set_ReadOutNoise(" << D_ReadOutNoise << ")" << endl;
  if (!F_Image.Set_ReadOutNoise( D_ReadOutNoise ))
  {
    cout << "MExtract::main: ERROR: F_Image.Set_ReadOutNoise(" << D_ReadOutNoise << ") returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// Set Gain
  cout << "MExtract::main: Starting F_Image.Set_Gain(" << D_Gain << ")" << endl;
  if (!F_Image.Set_Gain( D_Gain ))
  {
    cout << "MExtract::main: ERROR: F_Image.Set_Gain(" << D_Gain << ") returned FALSE" << endl;
    exit(EXIT_FAILURE);
  }

  /// Set Oversample to 10
  cout << "MExtract::main: Starting F_Image.Set_OverSample(10)" << endl;
  if (!F_Image.Set_OverSample( 10 ))
  {
    cout << "MExtract::main: ERROR: F_Image.Set_OverSample() returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// Set I_MaxIterSF
  cout << "MExtract::main: Starting F_Image.Set_MaxIterSF(" << I_MaxIterSF << ")" << endl;
  if (!F_Image.Set_MaxIterSF( I_MaxIterSF ))
  {
    cout << "MExtract::main: ERROR: F_Image.Set_MaxIterSF(" << I_MaxIterSF << ") returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// Set I_MaxIterSky
  cout << "MExtract::main: Starting F_Image.Set_MaxIterSky(" << I_MaxIterSky << ")" << endl;
  if (!F_Image.Set_MaxIterSky( I_MaxIterSky ))
  {
    cout << "MExtract::main: ERROR: F_Image.Set_MaxIterSky(" << I_MaxIterSky << ") returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// Set I_MaxIterSig
  cout << "MExtract::main: Starting F_Image.Set_MaxIterSig(" << I_MaxIterSig << ")" << endl;
  if (!F_Image.Set_MaxIterSig( I_MaxIterSig ))
  {
    cout << "MExtract::main: ERROR: F_Image.Set_MaxIterSig(" << I_MaxIterSig << ") returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// Read FitsFile
  cout << "MExtract::main: Starting F_Image.ReadArray()" << endl;
  if (!F_Image.ReadArray())
  {
    cout << "MExtract::main: ERROR: F_Image.ReadArray() returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// Set DatabaseFileName_In
  cout << "MExtract::main: Starting F_Image.SetDatabaseFileName(" << CS_DatabaseFileName_In << ")" << endl;
  if (!F_Image.SetDatabaseFileName(CS_DatabaseFileName_In))
  {
    cout << "MExtract::main: ERROR: F_Image.SetDatabaseFileName(" << CS_DatabaseFileName_In << ") returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// Read DatabaseFileName_In
  cout << "MExtract::main: Starting F_Image.ReadDatabaseEntry()" << endl;
  if (!F_Image.ReadDatabaseEntry())
  {
    cout << "MExtract::main: ERROR: F_Image.ReadDatabaseEntry() returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// Calculate Trace Functions
  cout << "MExtract::main: Starting F_Image.CalcTraceFunctions()" << endl;
  if (!F_Image.CalcTraceFunctions())
  {
    cout << "MExtract::main: ERROR: F_Image.CalcTraceFunctions() returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  cout << "MExtract::main: P_CS_ErrIn = " << *P_CS_ErrIn << ")" << endl;
  if (P_CS_ErrIn->GetLength() > 1){
    /// Set ErrFileName_In
    cout << "MExtract::main: Starting F_Image.SetErrFileName(" << *P_CS_ErrIn << ")" << endl;
    if (!F_Image.SetErrFileName(*P_CS_ErrIn))
    {
      cout << "MExtract::main: ERROR: F_Image.SetErrFileName(" << *P_CS_ErrIn << ") returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }

    /// Read Error image
    cout << "MExtract::main: Starting F_Image.ReadErrArray()" << endl;
    if (!F_Image.ReadErrArray())
    {
      cout << "MExtract::main: ERROR: F_Image.ReadErrArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }


  CString *P_CS_Temp = CS_FitsFileName_In.DToA(23.3453,2);
  cout << "MExtract: P_CS_Temp set to " << *P_CS_Temp << endl;
  /// Write aperture header information
  F_Image.WriteApHead(CString("aphead_")+CS_FitsFileName_In+CString(".head"));

  /// Calculate Profile Image
  seconds = time(NULL);
  cout << "MExtract::main: Starting F_Image.MkProfIm(): time = " << seconds << endl;

  CS_A1_Args(2).Set("OLD");

  if (!F_Image.MkProfIm(CS_A1_Args, PP_Args))
  {
    cout << "MExtract::main: ERROR: F_Image.MkProfIm() returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }
  seconds = time(NULL);
  cout << "MExtract::main: MkProfIm returned true at " << seconds << endl;

  CFits F_OutImage;
  /// Set CS_FitsFileName_In
  cout << "MExtract::main: Starting F_OutImage.SetFileName(" << CS_FitsFileName_In << ")" << endl;
  if (!F_OutImage.SetFileName(CS_FitsFileName_In))
  {
    cout << "MExtract::main: ERROR: F_OutImage.SetFileName(" << CS_FitsFileName_In << ") returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  ///Read FitsFile
  cout << "MExtract::main: Starting F_OutImage.ReadArray()" << endl;
  if (!F_OutImage.ReadArray())
  {
    cout << "MExtract::main: ERROR: F_OutImage.ReadArray() returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// Set CS_FitsFileName_Out
  cout << "MExtract::main: Starting F_OutImage.SetFileName(" << CS_FitsFileName_Out << ")" << endl;
  if (!F_OutImage.SetFileName(CS_FitsFileName_Out))
  {
    cout << "MExtract::main: ERROR: F_OutImage.SetFileName(" << CS_FitsFileName_Out << ") returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// change size of F_OutImage to (NApertures x NRows)
  if (!F_OutImage.SetNCols(F_Image.GetNRows()))
  {
    cout << "MExtract::main: ERROR: F_OutImage.SetNCols(" << F_Image.GetNRows() << ") returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  if (!F_OutImage.SetNRows(F_Image.Get_NApertures()))
  {
    cout << "MExtract::main: ERROR: F_OutImage.SetNRows(" << F_Image.Get_NApertures() << ") returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  F_OutImage.GetPixArray() = F_Image.GetSpec();
  F_OutImage.WriteArray();

  /// Write EcFromProfileOut 1D
  bool B_WithSky = false;
  if (I_Telluric > 0)
    B_WithSky = true;
  if (P_CS_EcFromProfileOut->GetLength() > 1)
  {
    cout << "MExtract::main: Starting to write EcFromProfileOut" << endl;
    if (!F_Image.ExtractSpecFromProfile(B_WithSky))
    {
      cout << "MExtract::main: ERROR: F_Image.ExtractSpecFromProfile() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    if (!F_OutImage.SetFileName(*P_CS_EcFromProfileOut))
    {
      cout << "MExtract::main: ERROR: F_OutImage.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    F_OutImage.GetPixArray() = F_Image.GetLastExtracted();//.transpose(secondDim, firstDim);
    if (!F_OutImage.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_OutImage.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// Write SP_Fit 1D
  if (P_CS_SPFitOut->GetLength() > 1)
  {
    if (!F_OutImage.SetFileName(*P_CS_SPFitOut))
    {
      cout << "MExtract::main: ERROR: F_OutImage.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    cout << "MExtract::main: Starting to write SPFit" << endl;
    F_OutImage.GetPixArray() = F_Image.GetSpecFit();//.transpose(secondDim, firstDim);
    if (!F_OutImage.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_OutImage.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// Write ImOut 2D
  if (P_CS_ImOut->GetLength() > 1){
    if (!F_Image.SetFileName(*P_CS_ImOut))
    {
      cout << "MExtract::main: ERROR: F_Image.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    cout << "MExtract::main: Starting to write ImOut" << endl;
    F_Image.GetPixArray() = F_Image.GetRecArray();
    if (!F_Image.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_Image.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// Write RecFitOut 2D
  if (P_CS_RecFitOut->GetLength() > 1){
    if (!F_Image.SetFileName(*P_CS_RecFitOut))
    {
      cout << "MExtract::main: ERROR: F_Image.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    cout << "MExtract::main: Starting to write RecFitOut" << endl;
    F_Image.GetPixArray() = F_Image.GetRecFitArray();
    if (!F_Image.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_Image.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// Write ProfileOut 2D
  if (P_CS_ProfileOut->GetLength() > 1){
    if (!F_Image.SetFileName(*P_CS_ProfileOut))
    {
      cout << "MExtract::main: ERROR: F_Image.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    cout << "MExtract::main: Starting to write ProfileOut" << endl;
    F_Image.GetPixArray() = F_Image.GetProfArray();
    if (!F_Image.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_Image.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// Write MaskOut 2D
  if (P_CS_MaskOut->GetLength() > 1){
    if (!F_Image.SetFileName(*P_CS_MaskOut))
    {
      cout << "MExtract::main: ERROR: F_Image.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    cout << "MExtract::main: Starting to write MaskOut" << endl;
    Array<int, 2> I_A2_MaskArray(F_Image.GetNRows(), F_Image.GetNCols());
    I_A2_MaskArray = F_Image.GetMaskArray();
    Array<double, 2> D_A2_MaskArray(F_Image.GetNRows(), F_Image.GetNCols());
    D_A2_MaskArray = 1. * I_A2_MaskArray;
    F_Image.GetPixArray() = D_A2_MaskArray;
    if (!F_Image.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_Image.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// Write SkyArrOut 2D
  if (P_CS_SkyArrOut->GetLength() > 1){
    if (!F_Image.SetFileName(*P_CS_SkyArrOut))
    {
      cout << "MExtract::main: ERROR: F_Image.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    cout << "MExtract::main: Starting to write SkyArrOut" << endl;
    F_Image.GetPixArray() = F_Image.GetRecSkyArray();
    cout << "MExtract::main: F_Image.GetRecSkyArray().rows() = " << F_Image.GetRecSkyArray().rows() << endl;
    cout << "MExtract::main: F_Image.GetRecSkyArray().cols() = " << F_Image.GetRecSkyArray().cols() << endl;
    if (!F_Image.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_Image.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// Write ErrOut 2D
  if (P_CS_ErrOut->GetLength() > 1){
    cout << "MExtract::main: Writing F_Image.GetErrArray()" << endl;
    if (!F_Image.SetFileName(*P_CS_ErrOut))
    {
      cout << "MExtract::main: ERROR: F_Image.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    cout << "MExtract::main: Starting to write ErrOut" << endl;
    cout << "MExtract::main: F_Image.GetErrArray().rows() = " << F_Image.GetErrArray().rows() << endl;
    cout << "MExtract::main: F_Image.GetErrArray().cols() = " << F_Image.GetErrArray().cols() << endl;
    F_Image.GetPixArray() = F_Image.GetErrArray();
    if (!F_Image.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_Image.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// Write aperture header information
  F_Image.WriteApHead(CString("aphead_")+CS_FitsFileName_In+CString(".head"));

  /// Create ErrOutEc
  if (P_CS_ErrFromProfileOut->GetLength() > 1){// || P_CS_ErrOutEc->GetLength() > 1){
    cout << "MExtract::main: Starting F_Image.ExtractErrors()" << endl;
    if (!F_Image.ExtractErrors())
    {
      cout << "MExtract::main: ERROR: F_Image.ExtractErrors() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// output extracted spectrum 1D
  cout << "MExtract::main: Starting to write EcOut" << endl;
  F_OutImage.GetPixArray() = F_Image.GetSpec();//.transpose(secondDim, firstDim);

  //cout << "MExctract: F_Image.GetSpec = " << F_Image.GetSpec() << endl;

  // Write Profile Image
/*  if (!F_OutImage.SetFileName(CS_FitsFileName_Out))
  {
    cout << "MExtract::main: ERROR: F_OutImage.SetFileName(" << CS_FitsFileName_Out << ") returned FALSE!" << endl;
    exit(EXIT_FAILURE);
}*/
  cout << "MExtract::main: Starting F_OutImage.WriteArray()" << endl;
  if (!F_OutImage.WriteArray())
  {
    cout << "MExtract::main: ERROR: F_OutImage.WriteArray() returned FALSE!" << endl;
    exit(EXIT_FAILURE);
  }

  /// Write ErrOutEc 1D
  cout << "MExtract::main: Writing F_Image.GetErrorsEc()" << endl;
  if (P_CS_ErrOutEc->GetLength() > 1)
  {
    if (!F_OutImage.SetFileName(*P_CS_ErrOutEc))
    {
      cout << "MExtract::main: ERROR: F_OutImage.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    cout << "MExtract::main: Starting to write ErrOutEc" << endl;
    F_OutImage.GetPixArray() = F_Image.GetErrorsEc();//.transpose(secondDim, firstDim);
    if (!F_OutImage.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_OutImage.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// Write ErrFromProfile 1D
  if (P_CS_ErrFromProfileOut->GetLength() > 1)
  {
    if (!F_OutImage.SetFileName(*P_CS_ErrFromProfileOut))
    {
      cout << "MExtract::main: ERROR: F_OutImage.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    cout << "MExtract::main: Starting to write ErrFromProfile" << endl;
    F_OutImage.GetPixArray() = F_Image.GetLastExtracted();//.transpose(secondDim, firstDim);
    if (!F_OutImage.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_OutImage.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// Write SkyOut 1D
  if (I_Telluric > 0 && P_CS_SkyOut->GetLength() > 1)
  {
    if (!F_OutImage.SetFileName(*P_CS_SkyOut))
    {
      cout << "MExtract::main: ERROR: F_OutImage.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    cout << "MExtract::main: Starting to write SkyOut" << endl;
    F_OutImage.GetPixArray() = F_Image.GetSky();//.transpose(secondDim, firstDim);
    if (!F_OutImage.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_OutImage.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /// Write SkyErrOut 1D
  if (I_Telluric > 0 && P_CS_SkyErrOut->GetLength() > 1)
  {
    if (!F_OutImage.SetFileName(*P_CS_SkyErrOut))
    {
      cout << "MExtract::main: ERROR: F_OutImage.SetFileName() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
    cout << "MExtract::main: Starting to write SkyErrOut" << endl;
    F_OutImage.GetPixArray() = F_Image.GetSkyError();//.transpose(secondDim, firstDim);
    if (!F_OutImage.WriteArray())
    {
      cout << "MExtract::main: ERROR: F_OutImage.WriteArray() returned FALSE!" << endl;
      exit(EXIT_FAILURE);
    }
  }
/*  fitsfile *P_FitsFileIn;
  fitsfile *P_FitsFileOut;
  int bitpixA, bitpixB, hdunumA, hdunumB, nhdusA, nhdusB, hdutypeA, hdutypeB;
  int anynulA, anynulB, extendA, extendB, simpleA, simpleB;
  int naxisA, naxisB;
  long pcountA, pcountB, gcountA, gcountB;
  long naxesA[2], naxesB[2];
  long fpixelA, fpixelB, nelementsA, nelementsB;
  int      countA, countB, Status;
  double *p_ArrayA, *p_ArrayB;
  float nullvalA, nullvalB;
  char strbufA[256], strbufB[256];

  Status=0;
  fits_open_file(&P_FitsFileIn, CS_FitsFileName_In.Get(), READONLY, &Status);
  fits_read_imghdr(P_FitsFileIn, 2, &simpleA , &bitpixA, &naxisA, naxesA,
                   &pcountA, &gcountA, &extendA, &Status);
  if (Status !=0)
  {
    printf("CFits::ReadArray: Error %d opening file %s\n", Status, CS_FitsFileName_In.Get());
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }
  cout << "CFits::ReadArray: FitsFileName <" << CS_FitsFileName_In.Get() << "> opened" << endl;
  cout << "CFits::ReadArray: FitsFileName contains <" << naxesA[1]  << "> rows!!!!!!! and <" << naxesA[0] << "> columns!!!!!!!!! naxes = <" << naxesA << ">" << endl;

  fits_open_file(&P_FitsFileOut, CS_FitsFileName_Out.Get(), READWRITE, &Status);
//  fits_read_imghdr(P_FitsFileOut, 2, &simpleB , &bitpixB, &naxisB, naxesB,
//                   &pcountB, &gcountB, &extendB, &Status);
  if (Status !=0)
  {
    printf("CFits::ReadArray: Error %d opening file %s\n", Status, CS_FitsFileName_Out.Get());
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }
  cout << "CFits::ReadArray: FitsFileName <" << CS_FitsFileName_Out.Get() << "> opened" << endl;
//  cout << "CFits::ReadArray: FitsFileName contains <" << naxesB[1]  << "> rows!!!!!!! and <" << naxesB[0] << "> columns!!!!!!!!! naxes = <" << naxesB << ">" << endl;
  nhdusB = fits_get_num_hdus(P_FitsFileOut, &hdunumB, &Status);
  if (Status !=0)
  {
    cout << "CFits::ReadArray: Error " << Status << " reading hdunum" << endl;
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }
  cout << "MExctract::main: nhdusB = " << nhdusB << ", hdunumB = " << hdunumB << endl;

  fits_get_hdu_type(P_FitsFileOut, &hdutypeB, &Status);
  if (Status !=0)
  {
    cout << "CFits::ReadArray: Error " << Status << " reading hdutype" << endl;
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }
  cout << "MExctract::main: hdutypeB = " << hdutypeB << endl;

  fits_movabs_hdu(P_FitsFileOut, 1, &hdutypeB, &Status);
  if (Status !=0)
  {
    cout << "CFits::ReadArray: Error " << Status << " moving to hdunum 1 " << endl;
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }
  cout << "MExctract::main: moved to hdunum 1" << endl;

  char card[FLEN_CARD];
  int  nkeysA, nkeysB, ii;
  fits_get_hdrspace(P_FitsFileIn, &nkeysA, NULL, &Status);
  cout << "MExctract::main: nkeysA = " << nkeysA << endl;
  fits_get_hdrspace(P_FitsFileOut, &nkeysB, NULL, &Status);
  cout << "MExctract::main: nkeysB = " << nkeysB << endl;

  for (ii = 1; ii <= nkeysA; ii++)  {
    fits_read_record(P_FitsFileIn, ii, card, &Status);
    cout << card << endl;
    fits_write_record(P_FitsFileOut, card, &Status);
    if (Status !=0)
    {
      cout << "CFits::ReadArray: Error " << Status << " writing header" << endl;
      char* P_ErrMsg = new char[255];
      ffgerr(Status, P_ErrMsg);
      cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
      delete[] P_ErrMsg;
      return false;
    }
  }


/*  fits_copy_hdu(P_FitsFileIn, P_FitsFileOut, &Status);
  if (Status !=0)
  {
    cout << "CFits::ReadArray: Error " << Status << " copying header" << endl;
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }
  cout << "MExctract::main: Header copied" << endl;
  /*

  nhdusA = fits_get_num_hdus(P_FitsFileIn, &hdunumA, &Status);
  if (Status !=0)
  {
    cout << "CFits::ReadArray: Error " << Status << " reading hdunum" << endl;
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }
  cout << "MExctract::main: nhdusA = " << nhdusA << ", hdunumA = " << hdunumA << endl;

  fits_get_hdu_type(P_FitsFileIn, &hdutypeA, &Status);
  if (Status !=0)
  {
    cout << "CFits::ReadArray: Error " << Status << " reading hdutype" << endl;
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }
  cout << "MExctract::main: hdutypeA = " << hdutypeA << endl;

  fits_movabs_hdu(P_FitsFileIn, 1, &hdutypeA, &Status);
  if (Status !=0)
  {
    cout << "CFits::ReadArray: Error " << Status << " moving to hdunum 1 " << endl;
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }
  cout << "MExctract::main: moved to hdunum 1" << endl;

  fits_copy_file(P_FitsFileIn, P_FitsFileOut, 0, 1, 1, &Status);
  if (Status !=0)
  {
    cout << "CFits::ReadArray: Error " << Status << " copying header" << endl;
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }
  cout << "MExctract::main: Header copied" << endl;
  */
/*  fits_close_file(P_FitsFileIn, &Status);
  if (Status !=0)
  {
    cout << "CFits::ReadArray: Error " << Status << " closing file " << CS_FitsFileName_In.Get() << endl;
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }

  fits_close_file(P_FitsFileOut, &Status);
  if (Status !=0)
  {
    cout << "CFits::ReadArray: Error " << Status << " closing file " << CS_FitsFileName_Out.Get() << endl;
    char* P_ErrMsg = new char[255];
    ffgerr(Status, P_ErrMsg);
    cout << "CFits::ReadArray: <" << P_ErrMsg << "> => Returning FALSE" << endl;
    delete[] P_ErrMsg;
    return false;
  }
  */
  delete(P_CS);
  delete(P_CS_ErrIn);
  delete(P_CS_ErrOut);
  delete(P_CS_SkyOut);
  delete(P_CS_SkyErrOut);
  delete(P_CS_ImOut);
  delete(P_CS_ProfileOut);
  delete(P_CS_ErrOutEc);
  delete(P_CS_SkyArrOut);
  delete(P_CS_MaskOut);
  delete(P_CS_SPFitOut);
  delete(P_CS_EcFromProfileOut);
  delete(P_CS_ErrFromProfileOut);
  return EXIT_SUCCESS;
}
