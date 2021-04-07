//
//
// 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include "Riostream.h"
#include <map>
#include <string>
//#include <boost/algorithm/string.hpp>
//#include <boost/algorithm/string/trim.hpp>
#include <vector>
#include <math.h>
//#include <TCint.h>
//#include <TGenericClassInfo.h> 
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TTree.h>
#include "TBranch.h"
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h> 
#include <TF1.h>  
#include <TF2.h> 
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TDSet.h"
#include "TChain.h"
#include <time.h> 
#include <TSystemDirectory.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TMinuit.h>
#include "Math/WrappedMultiTF1.h"
#include "TRandom.h" 
#include "TRandom3.h" 
#include  <TStopwatch.h>
#include "TH1F.h"
#include "TH2F.h"			// unused?
#include "TStyle.h"
#include "TCanvas.h"
#include <TGraphAsymmErrors.h>
#include <TFrame.h>
#include <TFitResult.h>
#include <TFitter.h>
#include "Fit/Fitter.h"
#include <TMatrixDSym.h>
#include <TBinomialEfficiencyFitter.h>
#include <TKDTreeBinning.h>
#include <TH2Poly.h>
#include <Math/Vector4Dfwd.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
//#include <RooMath/GenVector/PtEtaPhiM4D.h>
//#if !defined(__CINT__) || defined(__MAKECINT__)
#include <RooRandom.h>
#include <RooFit.h>
#include <RooMinuit.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooPlot.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooArgSet.h>
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsCategory.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooProduct.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooPolynomial.h>
#include <RooChebychev.h>
#include <RooWorkspace.h>
#include <RooExponential.h>
#include <RooErrorVar.h>
#include <RooFitResult.h>
#include <RooRangeBinning.h>
#include <RooBinning.h>
#include <RooNumGenConfig.h>
#include <RooBernstein.h>
#include <RooPolynomial.h>
#include <RooExtendPdf.h>
#include "RooFoamGenerator.h"
#include "RooAcceptReject.h"
#include "GBRMath.h"
#include "RooDoubleCBFast.h"
#include "RooBernsteinSideband.h"
#include "RooMCStudy.h"
#include "TFoam.h"
#include "TRatioPlot.h"

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>
// System stuff
#include <fstream> 
#include <sys/time.h>
#include <sys/times.h>

timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 

using namespace std; 
using namespace ROOT;
using namespace RooFit;


void FitSBModel();
void CreateInputHistoFile();
void replaceAll(std::string& str, const std::string& from, const std::string& to) ;
void replaceChar(char * txt, const char * txt1, const char * txt2) ;
//int RunEra=2017;
TTree* makeGenSBTTree(RooDataSet* geneDataNew, RooDataSet* geneMassNew,RooBernsteinSideband* BernSideBand,  TFile*OutFileNtupla); 
double FitMassSpectrumRoofit(RooDataSet* RooDataMass, TCanvas* c2, TH1D* masHist, TH1D* pdfHist, TH1D*sigHist, TH1D* bkgHist, int MaxDegreeBckg);
RooGaussian* _constrainVar(RooRealVar *var,RooWorkspace *w=0);
float*  _getFittedVar(const char* varName,RooWorkspace *w=0);
std::string grep(std::string& filename, std::string keyword);
bool wrongTagged = false;
bool SetMinuit2  = false;
bool Folded      = false;
bool integral    = false;
bool GenMiniMC   = false;
bool MCStudy     = false;
bool AnaMiniMC   = false;
bool ReFit       = false;
bool FirstMC     = false;

  char GooFitDir[300]               = "/scratch/users/dini/GooFit.2.1/BUILDNOTPOS-128-2017-dev3new-root60806/examples/";
  char AnaMiniMCDir[300]            = "";
  char RecoDir[100]		    = "/eos/cms/store/user/fiorendi/p5prime/[RunEra]/skims/newphi";
  char InputRecoB0TreeName[10]	    = "ntuple";
  char OutputRecoB0TreeName[10]     = "ntuple";
  char OutputGeneB0TreeName[10]     = "ntuple";
  char InputFileNameRecoB0[300]     = "[RunEra]Data_All_finalSelection.root";
  char ListParName[400] 	    =  "ListParValues-[RunEra]-Q2Bin-2-Bins-.txt";
  char ListPloName[400] 	    =  "ListParValues-[RunEra]-Q2Bin-2-Bins-.plo";
  char ListParNorm[410] 	    =  "ListParValues-[RunEra]-Q2Bin-2-Bins-.txt_norm";
  char ListPloNorm[410] 	    =  "ListParValues-[RunEra]-Q2Bin-2-Bins-.plo_norm";
  char FitStraName[400] 	    =  "namelist-[RunEra]-Q2Bin-2-Bins-.stra";
  char OutFileName[400] 	    =  "";
  char OutFileNameInputHisto[300]   =  "";
  char OutFileNameMiniMCHisto[300]  =  "";
  char OutSaveFileName[400]	    =  "";
  char LogFileName[400]	            =  "";
  char PDFNameRecoHisto[350]	    =  "B0-RecoHist-[RunEra]-Q2Bin-1.pdf";
  char PDFNameGeneHisto[350]	    =  "B0-GeneHist-[RunEra]-Q2Bin-1.pdf";
  char PNGNameMassHist[350]	    =  "B0-MassHist-[RunEra]-Q2Bin-1.png";
  char PNGNameMassQ2Hist[350]	    =  "B0-MassHist-[RunEra]-Q2Bin-1.png";
  char PNGNameCoeffPulls[350]	    =  "Coeff-Pulls-[RunEra]-Q2Bin-1.png";
  char PNGNameModErrPulls[350]	    =  "Coeff-ModErrPulls-[RunEra]-Q2Bin-1.png";
  char PNGNameCoeffNotGen[350]      =  "Coeff-NotGen-[RunEra]-Q2Bin-1.png";
  char PNGNameChi2PVal[350]         =  "Coeff-Chi2Pval-[RunEra]-Q2Bin-1.png";
  char PNGNameStudyMCPulls[350]	    =  "Coeff-StudyMCPulls-[RunEra]-Q2Bin-1.png";
  char PNGNameStudyMCParam[350]	    =  "Coeff-StudyMCParam-[RunEra]-Q2Bin-1.png";
  char PNGNameStudyMCError[350]     =  "Coeff-StudyMCError-[RunEra]-Q2Bin-1.png";
  char ProjectTXT[300]		    =  "";
  char SigmaMethodTXT[100]	    =  "";
  char TaggedVarTXT[100]	    =  "";
  char FoldedTXT[100]		    =  "";
  char MiniMCTXT[100]		    =  "";
  char AnaMiniTXT[100]		    =  "";
  char *RunEra;
//  char fitMassFileName[100]	    =  "results_fits_[RunEra]_newCB.root";
  char fitMassFileName[100]         =  "results_fits_[RunEra]_fM_newbdt.root";// RunEra  will be set after...
  char fitMassFileNameJpsi[100]     =  "results_fits_[RunEra]_fM_Jpsi_newbdt.root";// RunEra  will be set after...
  char FMTNSigma1L[10]		    ="";
  char FMTNSigma2L[10]		    ="";
  char FMTNSigma1R[10]		    ="";
  char FMTNSigma2R[10]		    ="";
//  char fitMassFileName[100]	    =  "results_fits_2017.root";
//  char fitMassFileName[100]	    =  "rf607_fitresult.root";

  
  TFile*OutFile = 0;

 

  char PDFNameMass[350]          = "B0-Mass-[RunEra]-global.pdf";
  char PNGNameFitSB3D[350]       = "B0-FitSB3D-[RunEra]-global.png";
  char PNGNamePlotSB3D[350]      = "B0-PlotSB3D-[RunEra]-global.png";
  char PNGNameReFitSB3D[350]       = "B0-ReFitSB3D-[RunEra]-global.png";
  char testo[300]     = "" ;
  float MarkerSizeSet = 0.35;
  int   PlotLineWidth = 1.;
  int   NCPU = 64;

//============================
// maxDegree START
// now defined in NAMELIST 
//============================
  int maxDegree1 =0;
  int maxDegree2 =0;
  int maxDegree3 =0;
// Mass Spectrum Bernstein   
  int maxDegree  =0;
  
  int fixParam   =1;
//============================
// now defined in NAMELIST 
// maxDegree END
//============================


//
// Il Bin in Q^2 !!!!
//   
  double Q2Min = 0.; 
  double Q2Max = 0.; 
  int    Q2Bin = -1;
  double Q2Fake= -1.;
//=================  
// Number of Normalization Integrals
//=================
  int    NormInteg = 11;
//=================
//=================
  int SETNumBinsX=100;
  int SETNumBinsY=100;
  int SETNumBinsZ=100;
//=================
//=================
//   
//=================
//=================
// NFact!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  int   NFact = 1; 
//=================
//     MCToys
//=================

  int   NFactGen = 1;
  int   NIniGen = 1;
  int   NSampleMC = 1000;
//=== Initialize counter ===  
  int   ivolte  = 0;
  int   SigmaProbSign=1;
// NFact!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//=================
//=================

  int FitPrintLevel=1;
  bool boolHesse = true;
  bool AutoFixPar= false;
  bool FirstFit  = false;
  int  iFitLoop  = 3;
  double xCoeffNorm =0.0;
  int    xCoeffIndex=-1;
//=================
//=================
//=================
//=================
  double ParMin =  -100000.;
  double ParMax =   100000.;
  double RndMin =  0.;
  double RndMax =  1.;
//=================
//=================
//=================
//=================

  double  tagged_mass_rangeMin=5.0;
  double  tagged_mass_rangeMax=5.6;

  double XMinSign = 4.9;
  double XMaxSign = 5.6;
  double B0Mass   = 5.27962;
  double B0Sigma  = 0.030;
  double JPsiMass = 3.096916;
  double PsiPMass = 3.686109;
  float  HistMassL1 = 4.935;
  float  HistMassL2 = 5.65;
  double SetMinRatio=0;
  double SetMaxRatio=3;
  double SetMinProj=0;

  double XMinSBL  = 0.;
  double XMaxSBL  = 0.;
  double XMinSBR  = 0.;
  double XMaxSBR  = 0.;
//
  double NSigma1L = 0.;
  double NSigma2L = 0.;
  double NSigma1R = 0.;
  double NSigma2R = 0.;
//
//   double NSigmaSBL = -2.;
//   double NSigmaSBR = 0;
//   double BiasSB   = 6;
  double XMinCosThetaL	       = -1.;
  double XMaxCosThetaL	       =  1.;
  double XMinCosThetaK	       = -1.;
  double XMaxCosThetaK	       =  1.;
//  double XMinPhi	       =0.;
  double XMinPhi	       =-TMath::Pi();
  double XMaxPhi	       = TMath::Pi();

  double XMinCosThetaLUnfolded = -1.;
  double XMaxCosThetaLUnfolded =  1.;
  double XMinCosThetaKUnfolded = -1.;
  double XMaxCosThetaKUnfolded =  1.;
  double XMinPhiUnfolded       =-TMath::Pi();
  double XMaxPhiUnfolded       = TMath::Pi();

  int	 xCosLHBin =  25;
  int	 xCosKHBin =  25;
  int	 xPhiHBin  =  25;
  double BinWCosThetaL= 0;
  double BinWCosThetaK= 0;
  double BinWPhi      = 0;
  
  
  double xMinQMuMu = 1.;
  double xMaxQMuMu = 19.;
  double NSigma  = 3.;
//  double XMinSignW = XMinSign;
//  double XMaxSignW = XMaxSign;
//  double XMinSignW = B0Mass - NSigma*B0Sigma;
//  double XMaxSignW = B0Mass + NSigma*B0Sigma;
  double NBckgInt3Sigma =0.;
  double NSignInt2Sigma =0.;
  double NBckgInt2Sigma =0.;
  double XLeftSet =0.;
  double XRightSet =0.;
  double XStepSign = 0.0025;
  double XStepMinuit = 0.00001;
  float xMassHBin = (XMaxSign -XMinSign)/XStepSign;
  float xQ2HBin   = (xMaxQMuMu -xMinQMuMu)/0.1;
  double XHScale = 10;
  
//   double yieldSignal = 1.41895e+06;
//   double yieldBckg   = 4.32096e+05;
  double yieldSignal = 1.37238e+06;
  double yieldBckg   = 6.98786e+05;
//   double ParMin = -1000;
//   double ParMax =  1000;
//   double RndMin = -0.1;
//   double RndMax =  0.1;
//  double c_const       = 0.0299792458;

  

  float xMassHBin2   =  xMassHBin /5; // plot only!
  
  TH1D* HxMass         = new TH1D( "HxMass"     , "B^{0} Mass"		 ,	      xMassHBin2, XMinSign, XMaxSign);
  TH1D* HxMassQ2       = new TH1D( "HxMassQ2"   , "B^{0} Mass"		 ,	      xMassHBin2, XMinSign, XMaxSign);
  TH1D* HxMassQ2SB     = new TH1D( "HxMassQ2SB" , "B^{0} Mass"		 ,	      xMassHBin2, XMinSign, XMaxSign);
  TH1D* pdfHxMass      = new TH1D( "pdfHxMass"  , "B^{0} Mass Fit"	 ,  XHScale * xMassHBin , XMinSign, XMaxSign);
  TH1D* sigHxMass      = new TH1D( "sigHxMass"  , "B^{0} Mass Fit"	 ,  XHScale * xMassHBin , XMinSign, XMaxSign);
  TH1D* bkgHxMass      = new TH1D( "bkgHxMass"  , "B^{0} Mass Fit"	 ,  XHScale * xMassHBin , XMinSign, XMaxSign);
  TH1D* pdfHxMassQ2    = new TH1D( "pdfHxMassQ2", "B^{0} Mass Fit Q2 Bin",  XHScale * xMassHBin , XMinSign, XMaxSign);
  TH1D* sigHxMassQ2    = new TH1D( "sigHxMassQ2", "B^{0} Mass Fit Q2 Bin",  XHScale * xMassHBin , XMinSign, XMaxSign);
  TH1D* bkgHxMassQ2    = new TH1D( "bkgHxMassQ2", "B^{0} Mass Fit Q2 Bin",  XHScale * xMassHBin , XMinSign, XMaxSign);
//
  TH2D* HxMassVsCosL   = new TH2D( "HxMassVsCosL","B^{0} Mass%CosL"    ,(int)xMassHBin2, XMinSign,  XMaxSign, NFact*xCosLHBin, XMinCosThetaL, XMaxCosThetaL);
  TH2D* HxMassVsCosK   = new TH2D( "HxMassVsCosK","B^{0} Mass%CosK"    ,(int)xMassHBin2, XMinSign,  XMaxSign, NFact*xCosKHBin, XMinCosThetaK, XMaxCosThetaK);
  TH2D* HxMassVsPhi    = new TH2D( "HxMassVsPhi", "B^{0} Mass%Phi"     ,(int)xMassHBin2, XMinSign,  XMaxSign, NFact*xPhiHBin , XMinPhi, XMaxPhi);
//TH1D* pdfHist        = new TH1D( "pdfHist", "B^{0} Mass Fit",     xMassHBin2, XMinSign,  XMaxSign);
//TH1D* sigHist        = new TH1D( "sigHist", "B^{0} Mass Fit",     xMassHBin2, XMinSign,  XMaxSign);
//TH1D* bkgHist        = new TH1D( "bkgHist", "B^{0} Mass Fit",     xMassHBin2, XMinSign,  XMaxSign);

// RooRealVar x("x", "x",  XMinCosThetaL,XMaxCosThetaL);
// RooRealVar y("y", "y",  XMinCosThetaK,XMaxCosThetaK);
// RooRealVar z("z", "z",  XMinPhi,XMaxPhi);
RooRealVar* ctL = new RooRealVar("ctL", "ctL",  XMinCosThetaK,XMaxCosThetaK);
RooRealVar* ctK = new RooRealVar("ctK", "ctK",  XMinCosThetaL,XMaxCosThetaL);
RooRealVar* phi = new RooRealVar("phi", "phi",  XMinPhi,XMaxPhi);
//
//
RooRealVar   BWidthX("BWidthX"  ,"BWidthX"  ,0., 2.*XMaxCosThetaL);
RooRealVar   BWidthY("BWidthY"  ,"BWidthY"  ,0., 2.*XMaxCosThetaK);
RooRealVar   BWidthZ("BWidthZ"  ,"BWidthZ"  ,0., 2.*XMaxPhi	 );

RooRealVar *tagged_mass=0;
RooRealVar *mumuMass=0;
RooRealVar *mumuMassE=0;
RooDataSet *geneMassNew=0;
RooDataSet *geneDataNew=0;
RooBernsteinSideband* BernSideBand=0;
RooBernsteinSideband* BernSideFits=0;
RooBernsteinSideband* Bern        =0;
RooAbsPdf  *bkg_mass_sb           =0;
RooAbsPdf  *bkg_exp               =0;
//RooBernstein        * Bern        =0;
int NumMassNewGen =0;
//
TTree* RecoB0TreeOut  =0;
TTree* GenMiniMCB0TreeOut =0;
TFile*OutFileInputHisto;
double EffiFunc3D(Double_t *var, Double_t *par);
double EffiFunc2D(Double_t *var, Double_t *par);
std::map<std::string, std::string>  ReadNamelist(int argc, char** argv);

TCanvas *csignstudy=0;

TRatioPlot* RatioDataModel3DX = 0;
TRatioPlot* RatioDataModel3DY = 0; 
TRatioPlot* RatioDataModel3DZ = 0; 


//==========================================
// Adaptive Binning...
//==========================================
TKDTreeBinning* RecoAdaptBinsC = 0;
TKDTreeID* TKDTreeIDC =0;
int   xAdaptNumBinC = 1;
int   MinContAdaptBin = 5;  
int   NDim = 3;
//
//==========================================
//==========================================
//=========    MAIN    =====================
//==========================================
//==========================================

int main (int argc, char** argv) {
//gSystem->Load("libRooDoubleCBFast.so");
// if (argc>1 ){
//     Q2Bin = (int) (*argv)[1];
// }
//     cout<<Q2Bin<<endl;
//     cout<<argv[0]<<endl;
//     exit(1);

if (argc<=1 ){
    cout<<"Q2Bin not set"<<endl;
    cout<<"Usage: "<<argv[0]<< " [QBin2] {where QBin2=0,1,2,3,4,5,6,7,8} [RunEra] {where RunEra=2016,2017,2018} [g] {option: generate 10 x datastatistic}\n"<<endl;
    cout<<"example: "<<argv[0]<< " 3 2017   {plot projections data&model}\n"<<endl;
    cout<<"example: "<<argv[0]<< " 3 2017 f {plot projections data&model &refit}\n"<<endl;
    cout<<"example: "<<argv[0]<< " 3 2017 g {plot projections data&model, generate MiniMC from fit parameterisation}\n"<<endl;
    exit(1);
}   


 
switch ( *argv[1] ) {

  case '0' : 
   Q2Min = 1.; 
   Q2Max = 2.;
   Q2Bin = 0;
   Q2Fake= 1.5;
    break;
  case '1' : 
   Q2Min = 2.; 
   Q2Max = 4.3; 
   Q2Bin = 1;
   Q2Fake= 3.;
    break;
  case '2' : 
   Q2Min = 4.3; 
   Q2Max = 6.; 
   Q2Bin = 2;
   Q2Fake= 5.;
    break;
  case '3' : 
   Q2Min = 6.;  
   Q2Max = 8.68; 
   Q2Bin = 3;
   Q2Fake= 7.;
    break;
  case '4' : 
   Q2Min = 8.68; 
   Q2Max = 10.09; 
   Q2Bin = 4;
   Q2Fake= 9.;
   sprintf(fitMassFileName,fitMassFileNameJpsi);
    break;
  case '5' : 
   Q2Min = 10.09; 
   Q2Max = 12.86; 
   Q2Bin = 5;
   Q2Fake= 11.;
    break;
  case '6' : 
   Q2Min = 12.86; 
   Q2Max = 14.18; 
   Q2Bin = 6;
   Q2Fake= 13.;
    break;
  case '7' : 
   Q2Min = 14.18; 
   Q2Max = 16.; 
   Q2Bin = 7;
   Q2Fake= 15.;
    break;
  case '8' : 
   Q2Min = 16; 
   Q2Max = 19.; 
   Q2Bin = 8;
   Q2Fake= 17.;
    break;

  default : 
    // Process for all other cases.
    cout<<"Q2Bin not set correctly!!!"<<endl;
    cout<<"Usage: "<<argv[0]<< " QBin2 [where QBin2=0,1,2,3,4,5,6,7,8]\n"<<endl;
    exit(1);

}
//======================================================================
//======================================================================
//===========		InputFileNameMCGene              ===============
//======================================================================
//======================================================================
//  sprintf(InputFileNameMCGene,"testGene-2017-Q2Bin-%d.root",Q2Bin);
//======================================================================
//======================================================================
//======================================================================
    std::cout<<"====================================="<<endl;
    
  if (argc>2 && ((strcmp(argv[2],"2016") == 0)||(strcmp(argv[2],"2017") == 0)||(strcmp(argv[2],"2018") == 0)) ){
     RunEra = (char *) malloc(strlen(argv[2])+1);
     strcpy(RunEra,argv[2]);
     std::cout<<"========================================================================="<<endl;
     std::cout<<Form("Setting the [RunEra]: %s",argv[2])<<std::endl;
     std::cout<<"========================================================================="<<endl;
  }else{
    std::cout<<"Error!!! Please set the [RunEra] [2016,2017,2018]...exit."<<std::endl;
    exit(1);
  }
    
  if (argc>3 && (strcmp(argv[3],"g") == 0) ){
//    sprintf(MiniMCTXT,"MiniMC-NFactGen%d",NFactGen);
    GenMiniMC = true;
    std::cout<<"========================================================================="<<endl;
    std::cout<<"Setting the option: GenMiniMC events                                     "<<std::endl;
    std::cout<<"========================================================================="<<endl;
  }
  if (argc>3 && (strcmp(argv[3],"f") == 0) ){
    ReFit = true;
    std::cout<<"========================================================================="<<endl;
    std::cout<<"Setting the option: ReFit events                                         "<<std::endl;
    std::cout<<"========================================================================="<<endl;
  }
  if (argc>3 && (strcmp(argv[3],"m") == 0) ){
    MCStudy = true;
    std::cout<<"========================================================================="<<endl;
    std::cout<<"Setting the option: Roofit MC Study                                      "<<std::endl;
    std::cout<<"========================================================================="<<endl;
  }
//GooFitDir
  sprintf(GooFitDir,"%s/testSB3DB0-%s-newglobal/",GooFitDir,RunEra);
//
  if (argc>3 && (strcmp(argv[3],"a") == 0) ){
    if(argc>4 && (strcmp(argv[4],"auto") == 0) ){
//     sprintf(MiniMCTXT,"MiniMC-NFactGen%d",NFactGen);
     sprintf(AnaMiniMCDir,"%s/AutoMiniMC",GooFitDir);
    }else{
//     sprintf(MiniMCTXT,"MiniMC-NFactGen%d",NFactGen);
     sprintf(AnaMiniMCDir,"%s/MiniMC",GooFitDir);
    } 
    AnaMiniMC = true;
    std::cout<<"========================================================================="<<endl;
    std::cout<<"Setting the option: AnaMiniMC event                                      "<<std::endl;
    std::cout<<"========================================================================="<<endl;
  }
  
  
  
  char NameList[300];;
  sprintf(NameList,"namelist-SB3DB0-%s-Q2Bin-%d.lis",RunEra, Q2Bin);
//  
//   if(Folded){
//     sprintf(FoldedTXT,"-PhiFolded");
//     XMinPhi = 0.;
//     ParMin = 0.;
//     RndMin = 0.;
//     std::cout<<"===================================================================="<<endl;
//     std::cout<<"======== SETTING: Phi Ang.Variable FOLDED             =============="<<std::endl;
//     std::cout<<"===================================================================="<<endl;
//   };
  
  
  char*argn[]={NameList};
  
  std::map<std::string, std::string> mappa = ReadNamelist(1,argn );
//
  maxDegree1	    =	 atoi (mappa["maxDegree1"].c_str() ) ;
  maxDegree2	    =	 atoi (mappa["maxDegree2"].c_str() ) ;
  maxDegree3	    =	 atoi (mappa["maxDegree3"].c_str() ) ;
  xCosLHBin	    =	 atof (mappa["xCosLHBin" ].c_str() ) ;
  xCosKHBin	    =	 atof (mappa["xCosKHBin" ].c_str() ) ;
  xPhiHBin	    =	 atof (mappa["xPhiHBin"  ].c_str() ) ;
  NSigma1L	    =	 atof (mappa["NSigma1L"  ].c_str() ) ;
  NSigma2L	    =	 atof (mappa["NSigma2L"  ].c_str() ) ;
  NSigma1R	    =	 atof (mappa["NSigma1R"  ].c_str() ) ;
  NSigma2R	    =	 atof (mappa["NSigma2R"  ].c_str() ) ;
  maxDegree	    =	 atoi (mappa["maxDegree" ].c_str() ) ;
  fixParam	    =	 atoi (mappa["fixParam"  ].c_str() ) ; 
  map<string,string>::iterator  it= mappa.find("MinContAdaptBin");
  if(it != mappa.end()) {
   MinContAdaptBin   =    atoi (mappa["MinContAdaptBin"  ].c_str() ) ;
  }
  map<string,string>::iterator  it2= mappa.find("AutoFixPar");
  if(it2 != mappa.end()) {
   AutoFixPar   =    atol (mappa["AutoFixPar"  ].c_str() ) ;
  }
  if(AutoFixPar){
   std::cout<<"Warning: setting search for Param to fix = "<<AutoFixPar<<std::endl;
  }
  map<string,string>::iterator  it3= mappa.find("NFactGen");
  if(it3 != mappa.end()) {
   NFactGen   =    atoi (mappa["NFactGen"  ].c_str() ) ;
   std::cout<<"Warning: setting  NFactGen from namelist= "<<NFactGen<<std::endl;     
   if(NFactGen<6&&AnaMiniMC==true){
    sprintf(AnaMiniMCDir,"%s/DataLike",AnaMiniMCDir);
    std::cout<<Form("Warning: setting  DataLike")<<std::endl;     
    std::cout<<Form("DataLike ==> AnaMiniMCDir=%s",AnaMiniMCDir)<<std::endl;     
   }  

  }
  map<string,string>::iterator  it4= mappa.find("SigmaProbSign");
  if(it4 != mappa.end()) {
   SigmaProbSign   =    atoi (mappa["SigmaProbSign"  ].c_str() ) ;
   std::cout<<"Warning: setting  SigmaProbSign from namelist= "<<SigmaProbSign<<std::endl;
  }
  map<string,string>::iterator  it5= mappa.find("NSampleMC");
  if(it5 != mappa.end()) {
   NSampleMC   =    atoi (mappa["NSampleMC"  ].c_str() ) ;
   std::cout<<"Warning: setting  NSampleMC from namelist= "<<NSampleMC<<std::endl;
  }
  map<string,string>::iterator  it6= mappa.find("NIniGen");
  if(it6 != mappa.end()) {
   NIniGen   =    atoi (mappa["NIniGen"  ].c_str() ) ;
   std::cout<<"Warning: setting  NIniGen from namelist= "<<NIniGen<<std::endl;
  }
  map<string,string>::iterator  it7= mappa.find("tagged_mass_rangeMin");
  if(it7 != mappa.end()) {
   tagged_mass_rangeMin   =    atof (mappa["tagged_mass_rangeMin"  ].c_str() ) ;
   std::cout<<"Warning: setting  tagged_mass_rangeMin from namelist= "<<tagged_mass_rangeMin<<std::endl;
   if(tagged_mass_rangeMin<XMinSign){
    std::cout<<Form("Error: setting  tagged_mass_rangeMin=%f < XMinSign=%f",tagged_mass_rangeMin,XMinSign)<<std::endl;
    exit(0);
   }
  }
  map<string,string>::iterator  it8= mappa.find("tagged_mass_rangeMax");
  if(it8 != mappa.end()) {
   tagged_mass_rangeMax   =    atof (mappa["tagged_mass_rangeMax"  ].c_str() ) ;
   std::cout<<"Warning: setting  tagged_mass_rangeMax from namelist= "<<tagged_mass_rangeMax<<std::endl;
   if(tagged_mass_rangeMax>XMaxSign){
    std::cout<<Form("Error: setting  tagged_mass_rangeMax=%f < XMaxSign=%f",tagged_mass_rangeMax,XMaxSign)<<std::endl;
    exit(0);
   }
  }

  std::cout<<" Num Param Bernstein polynomial CosL:	"<<maxDegree1<<std::endl;
  std::cout<<" Num Param Bernstein polynomial CosK:	"<<maxDegree2<<std::endl;
  std::cout<<" Num Param Bernstein polynomial Phi :	"<<maxDegree3<<std::endl;
  std::cout<<" Binning choice for CosL            :	"<<xCosLHBin<<std::endl;
  std::cout<<" Binning choice for CosK            :	"<<xCosKHBin<<std::endl;
  std::cout<<" Binning choice for Phi             :	"<<xPhiHBin <<std::endl;
//  
  std::cout<<" Min CosL XMinCosThetaL             :	"<<XMinCosThetaL<<std::endl;
  std::cout<<" Max CosL XMaxCosThetaL             :	"<<XMaxCosThetaL<<std::endl;
  std::cout<<" Min CosK XMinCosThetaK             :	"<<XMinCosThetaK<<std::endl;
  std::cout<<" Min CosK XMaxCosThetaK             :	"<<XMaxCosThetaK<<std::endl;
  std::cout<<" Min Phi  XMinPhi		          :	"<<XMinPhi<<std::endl;
  std::cout<<" Min Phi  XMaxPhi		          :	"<<XMaxPhi<<std::endl;
//
  if(SigmaProbSign==0){
   sprintf(SigmaMethodTXT,"-SigmaGauss");
   std::cout<<" NSigma1L [sigma gauss model]:  "<<NSigma1L<<std::endl;
   std::cout<<" NSigma2L [sigma gauss model]:  "<<NSigma2L<<std::endl;
   std::cout<<" NSigma1R [sigma gauss model]:  "<<NSigma1R<<std::endl;
   std::cout<<" NSigma2R [sigma gauss model]:  "<<NSigma2R<<std::endl;
  }else if(SigmaProbSign==-1){
   sprintf(SigmaMethodTXT,"-SigmaBare");
   std::cout<<" NSigma1L [bare min limit left ]:"<<NSigma1L<<std::endl;
   std::cout<<" NSigma2L [bare max limit left ]:"<<NSigma2L<<std::endl;
   std::cout<<" NSigma1R [bare min limit right]:"<<NSigma1R<<std::endl;
   std::cout<<" NSigma2R [bare max limit right]:"<<NSigma2R<<std::endl;
  }else if (SigmaProbSign==1){
   sprintf(SigmaMethodTXT,"-SigmaProb");
   std::cout<<" NSigma1L [Limit in GeV Left ]:  "<<NSigma1L<<std::endl;
   std::cout<<" NSigma2L [n. gauss stand. dev. sign]:  "<<NSigma2L<<std::endl;
   std::cout<<" NSigma1R [n. gauss stand. dev. sign]:  "<<NSigma1R<<std::endl;
   std::cout<<" NSigma2R [Limit in GeV Right]:  "<<NSigma2R<<std::endl;
  }else{
   std::cout<<Form(" SigmaProbSign: 	INVALID OPTION: %d !!! Exit...",SigmaProbSign)<<std::endl;
   exit(1);
  }
  std::cout<<" Num Param Bernstein polynomial Mass:	"<<maxDegree<<std::endl;
  std::cout<<" Parameter to Fix for normalization :     "<<fixParam<<std::endl;
  std::cout<<" Num of events inside adaptive bins :     "<<MinContAdaptBin<<std::endl;
  
//
  if (NSigma1L==0.0 ||
      NSigma2L==0.0 ||
      NSigma1R==0.0 ||
      NSigma2R==0.0){
     std::cout<<"====================================================================================="<<std::endl;
     std::cout<<"====================================================================================="<<std::endl;
     std::cout<<"====================================================================================="<<std::endl;
     std::cout<<"Error reading sideband limits: at least one NSigma[]==0 found !!!!"<<std::endl;
     std::cout<<"====> EXIT from Main!!!"<<std::endl;
     std::cout<<"====================================================================================="<<std::endl;
     std::cout<<"====================================================================================="<<std::endl;
     std::cout<<"====================================================================================="<<std::endl;
     exit(0);
  }    


//
  char TXTNSigma1L[10]="";
  char TXTNSigma2L[10]="";
  char TXTNSigma1R[10]="";
  char TXTNSigma2R[10]="";

  std::stringstream sss;
  std::string ss;
  std::string dot1=".";
  std::string dot2="dot";

  sss<<NSigma1L;
  ss=sss.str();
  cout<<ss<<endl;
  strcpy(FMTNSigma1L,ss.c_str());
  replaceAll( ss,  dot1, dot2);
  strcpy(TXTNSigma1L,ss.c_str());
  sss.str("");
  sss.clear();
  sss<<NSigma1R;
  ss=sss.str();
  cout<<ss<<endl;
  strcpy(FMTNSigma1R,ss.c_str());
  replaceAll( ss,  dot1, dot2);
  strcpy(TXTNSigma1R,ss.c_str());
  sss.str("");
  sss.clear();
  sss<<NSigma2L;
  ss=sss.str();
  cout<<ss<<endl;
  strcpy(FMTNSigma2L,ss.c_str());
  replaceAll( ss,  dot1, dot2);
  strcpy(TXTNSigma2L,ss.c_str());
  sss.str("");
  sss.clear();
  sss<<NSigma2R;
  ss=sss.str();
  cout<<ss<<endl;
  strcpy(FMTNSigma2R,ss.c_str());
  replaceAll( ss,  dot1, dot2);
  strcpy(TXTNSigma2R,ss.c_str());
//   
  sprintf(MiniMCTXT,"MiniMC-NFactGen%d",NFactGen);
//
  
  std::cout<<"--------------------------------------------\n"<<endl;
  std::cout<<" Setting selection for Q^2 bin: "<<*argv[1]<<" ==> "<<Q2Min<<"<Q^2<"<<Q2Max<<std::endl;
  std::cout<<"--------------------------------------------\n"<<endl;
  if(AnaMiniMC){
   sprintf(OutFileName,"testSidebandFit-AnaMiniMC-%s%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s.root.%s",RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
   sprintf(LogFileName,"MiniMCToys-Fit-%s%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s.log.%s",RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
    
//   sprintf(LogFileName,"MiniMCToys-Fit-%s-%d.log.%s%d_",RunEra,Q2Bin,MiniMCTXT); 
  }else if(MCStudy){
   sprintf(OutFileName,"testSidebandFit-MCStudy-%s%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s.root.%s",RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
  }else if(GenMiniMC){
   sprintf(OutFileName,"testSidebandFit-GenStudy-%s%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s.root.%s",RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
  }else{
   sprintf(OutFileName,"testGoofitSB3DB0-%s%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s.root",RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT); 
  }
  replaceChar(RecoDir,"[RunEra]",RunEra);
//  sprintf(RecoDir,"~/p5prime/%s/skims/newphi",RunEra);
//  sprintf(InputFileNameRecoB0,"%sData_All_finalSelection.root",RunEra);
  replaceChar(InputFileNameRecoB0,"[RunEra]",RunEra);
//  sprintf(fitMassFileName,"results_fits_%s_fixDataPdf.root",RunEra);
//   std::stringstream sssname,sssera,ssstxt;
//   sssname<<fitMassFileName;
//   std::string ssname=sssname.str();
//   sssera<<RunEra;
//   ssstxt<<"[RunEra]";
//   replaceAll(ssname,ssstxt.str(),sssera.str());
//   strcpy(fitMassFileName,ssname.c_str());
  replaceChar(fitMassFileName,"[RunEra]",RunEra);
// 
  sprintf(ListParName,"ListParValues-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s.txt"   ,RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT); 
  sprintf(ListPloName,"ListParValues-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s.plo"   ,RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT); 
  sprintf(PNGNameFitSB3D,"B0-FitSB3D-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s-Adapt-%d.png",RunEra,ProjectTXT,  Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT,MinContAdaptBin); 
  sprintf(PNGNamePlotSB3D,"B0-PlotSB3D-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s-Adapt-%d.png",RunEra,ProjectTXT,  Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT,MinContAdaptBin); 
  sprintf(PNGNameReFitSB3D,"B0-ReFitSB3D-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s-Adapt-%d.png",RunEra,ProjectTXT,  Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT,MinContAdaptBin); 
  sprintf(PNGNameMassQ2Hist,"B0-MassQ2-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s.png" ,RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
  sprintf(PNGNameMassHist,"B0-MassTot-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s.png"  ,RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
  sprintf(PNGNameCoeffPulls ,"Coeff-Pulls-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s.png"  ,RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
  sprintf(PNGNameModErrPulls ,"Coeff-ModErrPulls-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s.png"  ,RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
  sprintf( PNGNameStudyMCPulls,"Coeff-StudyMCPulls-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s.png"  ,RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
  sprintf( PNGNameStudyMCParam,"Coeff-StudyMCParam-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s.png"  ,RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
  sprintf( PNGNameStudyMCError,"Coeff-StudyMCError-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s.png"  ,RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
  sprintf(PNGNameCoeffNotGen,"Coeff-notGen-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s.png"  ,RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
  sprintf(PNGNameChi2PVal   ,"Coeff-Chi2PVal-%s-%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s%s.png"  ,RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT,MiniMCTXT); 
//  sprintf(fitMassFileName,"save-w-2017-bin-%d.root",Q2Bin);
  sprintf(OutSaveFileName,"savesb-%s%s-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-WSBL-%s-%s-WSBR-%s-%s%s%s%s.root",RunEra,ProjectTXT, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3,TXTNSigma1L,TXTNSigma2L,TXTNSigma1R,TXTNSigma2R,TaggedVarTXT,SigmaMethodTXT,FoldedTXT); 
  sprintf(OutFileNameInputHisto,"testGoofitSB3DB0-%s-InputHisto-Q2Bin-%d-Bins-%d-%d-%d-masspectrum%s.root",RunEra, Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,FoldedTXT); 

  if(GenMiniMC){
    sprintf(OutFileNameMiniMCHisto,"testGoofitSB3DB0-%s-InputHisto-Q2Bin-%d-Bins-%d-%d-%d-masspectrum%s.root.%s",RunEra,Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,FoldedTXT,MiniMCTXT); 
  }


  
  BinWCosThetaL=(XMaxCosThetaL-XMinCosThetaL)/double(xCosLHBin);
  BinWCosThetaK=(XMaxCosThetaK-XMinCosThetaK)/double(xCosKHBin);
  BinWPhi      =(XMaxPhi-XMinPhi)/double(xPhiHBin);

//  TApplication app("app",&argc, argv);
 
  TStopwatch TimeWatch;
  TimeWatch.Start();

  FitSBModel(); 
//  app.Run() ;
  cout<<"esco..." <<endl;
  TimeWatch.Stop();
  TimeWatch.Print();
  
  return 0 ;
}


void FitSBModel(){


  gettimeofday(&startTime, NULL);
  startCPU = times(&startProc);

   gROOT->Reset();
//   gROOT->SetStyle("Plain");
//   gROOT->ForceStyle();
   //gStyle->SetOptStat("e");
//   gStyle->SetOptStat("MNReou");
   //gStyle->SetOptFit(0000111);
   gStyle->SetPadBorderMode(0);
   gStyle->SetOptStat(000000);
   gStyle->SetOptFit(000000);

   TCanvas* c1 = new TCanvas("c1","Mass",200,10,900,780);
   TCanvas* c2 = new TCanvas("c2","Fit Mass Spectrum",200,10,900,780);
//    TCanvas* c5 = new TCanvas("c5","                      ",200,10,900,780);
   TCanvas* c6 = new TCanvas("c6","Sideband",200,10,900,780);
   TCanvas* cxy = new TCanvas("cxy","Sideband",200,10,750,750);
   TCanvas* cyz = new TCanvas("cyz","Sideband",200,10,750,780);
   TCanvas* cxz = new TCanvas("czy","Sideband",200,10,750,780);
   csignstudy = new TCanvas("csignstudy","Mass Signal study",200,10,900,450);

//    TCanvas* c7 = new TCanvas("c7","Efficiencies" ,200,10,900,780);
//    TCanvas* c8 = new TCanvas("c8","Effi-Reco Closure Test" ,200,10,900,780);
//   c2->Divide(2,2);  
//    c5->Divide(2,2);  
   c6->Divide(2,2);  
   csignstudy->Divide(2,1);
//    c7->Divide(2,2);  
//    c8->Divide(2,2);  
//    TPad* pad2 = (TPad*)c2->GetPad(0);
//    pad2->SetLeftMargin(0.15); 
//    pad2->SetRightMargin(0.15); 
//    TPad* pad6 = (TPad*)c6->GetPad(0);
//    pad6->SetLeftMargin(0.15); 
//    pad6->SetRightMargin(0.15); 

   //TPad* pad1 = (TPad*)c1->GetPad(0);
//   TPad* pad2 = (TPad*)c2->GetPad(0);
   //pad1->SetLeftMargin(0.15); 
//   pad2->SetLeftMargin(0.15); 



//  gSystem->Exec(Form("mv %s %s.tmp",OutFileName,OutFileName));
//  OutFileInputHisto = TFile::Open(OutFileNameInputHisto,"READ");
  if (!TFile::Open(OutFileNameInputHisto,"READ"))
  {
    cout<<"File:"<<OutFileNameInputHisto<<" not found!!! create..."<<endl;
    CreateInputHistoFile();
    OutFileInputHisto = TFile::Open(OutFileNameInputHisto,"READ");
  }else{
   OutFileInputHisto = TFile::Open(OutFileNameInputHisto,"READ");
   cout<<"File:"<<OutFileNameInputHisto<<" FOUND !!!"<<endl;
  } 
  gSystem->Exec(Form("mv %s %s.tmp",OutFileName,OutFileName));
  OutFile = TFile::Open(OutFileName,"RECREATE");
  TTree *RecoB0TreeOut     = (TTree*)OutFileInputHisto->Get(OutputRecoB0TreeName);
   if(!RecoB0TreeOut ){
     cout<<"TTree Reco Data: "<< OutputRecoB0TreeName <<" not found!!! Suggestion: remove this file e try again..."<<endl;
     exit(1);
   }else{
     cout<<"TTree Reco Data: "<< OutputRecoB0TreeName <<" OK FOUND!!!"<<endl;
   }  
  HxMass    = (TH1D*)OutFileInputHisto->Get("HxMass");
  if(!HxMass ){
    cout<<"HxMass Histo: not found!!! Exit..."<<endl;
    exit(1);
  }else{
    cout<<"HxMass Histo: OK FOUND!!! Entries: "<<HxMass->GetEntries()<<endl;
  }
//    HxMassQ2    = (TH1D*)OutFileInputHisto->Get("HxMassQ2");
//   if(!HxMassQ2 ){
//     cout<<"HxMassQ2 Histo: not found!!! Exit..."<<endl;
//     exit(1);
//   }else{
//     cout<<"HxMassQ2 Histo: OK FOUND!!! Entries: "<<HxMassQ2->GetEntries()<<endl;
//   } 
 

  
//   HxReco    = (TH3D*)OutFileInputHisto->Get("HxReco");
//   if(!HxReco ){
//     cout<<"HxReco Histo: not found!!! Exit..."<<endl;
//     exit(1);
//   }else{
//     cout<<"HxReco Histo: OK FOUND!!! Entries: "<<HxReco->GetEntries()<<endl;
//     if(HxReco->GetNbinsX()!=xCosLHBin){cout<<"Error HxReco NBinsX = "<<HxReco->GetNbinsX()<<" != xCosLHBin = "<<xCosLHBin<<endl;exit(1);}
//     if(HxReco->GetNbinsY()!=xCosKHBin){cout<<"Error HxReco NBinsY = "<<HxReco->GetNbinsY()<<" != xCosKHBin = "<<xCosKHBin<<endl;exit(1);}
//     if(HxReco->GetNbinsZ()!=xPhiHBin ){cout<<"Error HxReco NBinsZ = "<<HxReco->GetNbinsZ()<<" != xPhiHBin  = "<<xPhiHBin <<endl;exit(1);}
//   }  
  cout<<"======================="<<endl;
  cout<<"xCosLHBin = "<<xCosLHBin<<endl;
  cout<<"xCosKHBin = "<<xCosKHBin<<endl;
  cout<<"xPhiHBin  = "<<xPhiHBin <<endl; 
  cout<<"======================="<<endl;
//
  TH3D* HxReco = new   TH3D( "HxReco"    , "B^{0} Reco correct tagged",  xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
 									 xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
									 xPhiHBin , XMinPhi, XMaxPhi );
//
  int	numParameters = (maxDegree1+1)*(maxDegree2+1)*(maxDegree3+1);

  int icc=0;
  char ParCheck[numParameters][30];
  for(int i = 0; i <= maxDegree1 ; ++i) {
    for(int j = 0; j <= maxDegree2 ; ++j) {
     for(int k = 0; k <= maxDegree3 ; ++k) {
            
      	    sprintf(ParCheck[icc], "cosL=%d cosK=%d phi=%d", i,j,k); 
//	    cout<<Form("ParCheck(%d)= %s",icc,ParCheck[icc])<<endl;
	    icc++;
     }
    }
  }
 
  double cos_theta_l	;
  double cos_theta_k	;
  double phi_kst_mumu	;
  double xtagged_mass	;
  double xmumuMass	;
  double xmumuMassE     ;
  double mmk1	        ;
  double mmk2           ;
  RecoB0TreeOut->SetBranchAddress("cos_theta_l"   ,&cos_theta_l );
  RecoB0TreeOut->SetBranchAddress("cos_theta_k"   ,&cos_theta_k );
  RecoB0TreeOut->SetBranchAddress("phi_kst_mumu"  ,&phi_kst_mumu);
  RecoB0TreeOut->SetBranchAddress("tagged_mass"   ,&xtagged_mass );
  RecoB0TreeOut->SetBranchAddress("mumuMass"      ,&xmumuMass );
  RecoB0TreeOut->SetBranchAddress("mumuMassE"     ,&xmumuMassE );
  RecoB0TreeOut->SetBranchAddress("mmk1"          ,&mmk1 );
  RecoB0TreeOut->SetBranchAddress("mmk2"          ,&mmk2 );
  int nentries = (int)RecoB0TreeOut->GetEntries();
//  int nentries = 0;
//  nentries = 3000.;
  cout<<"nentries: "<< nentries<<endl;
  tagged_mass = new RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", XMinSign, XMaxSign, "GeV");
  mumuMass    = new RooRealVar("mumuMass"    , "mumuMass" , 0, 6);
  mumuMassE   = new RooRealVar("mumuMassE"   , "mumuMassE", 0, 10000);
  RooDataSet *fulldata   = new RooDataSet("fulldata", "fulldataset",  RooArgSet(*tagged_mass,*mumuMass,*mumuMassE));
  for (Int_t i=0;i<nentries;i++) {
         RecoB0TreeOut->GetEntry(i);
//  	 if(mmk2>3.6&&mmk2<4.2&&mmk1>4.7&&mmk1<4.9&&Q2Bin==4)continue;
         if( (xtagged_mass>XMinSign&&xtagged_mass<XMaxSign) ){
	   HxMassQ2->Fill(xtagged_mass);
	   tagged_mass->setVal(xtagged_mass);
	   mumuMass->setVal(xmumuMass);
	   mumuMassE->setVal(xmumuMassE);
	   fulldata->add(RooArgSet(*tagged_mass,*mumuMass,*mumuMassE));
	 } 
  }
//  
  RooArgList *coefLis = new RooArgList();
  RooArgList *coefFit = new RooArgList();
//
  int NumParamFree = 0;
  std::vector<RooRealVar> parLis; 
  std::vector<RooRealVar> errLis; 
  std::vector<RooRealVar> parFit; 
  std::vector<RooRealVar> errFit; 
  TH1F * Hparams         = new TH1F ( "Hparams"     , "N. parameter fit not gen",numParameters, 0., float(numParameters));

  char varName[100];
  char errName[100];
//  double ParMaxTemp=0.0;
  double parIni=0.0;
  double errIni=0.0;
  std::string line;
  std::cout<<"Try to open list of parameters :"<< ListParName <<std::endl;
  std::fstream *parListInput = new std::fstream(ListParName,std::ifstream::in);
  if(parListInput->is_open()){
     std::cout<<"List of initial parameters :"<< ListParName <<" FOUND!!!"<<std::endl;
     for (int i=0;i<numParameters;++i){
      	    sprintf(varName, "p%02d_%s", i, RunEra);
      	    sprintf(errName, "e%02d_%s", i, RunEra);
// 	    *parListInput >> parIni;
            std::getline(*parListInput, line);
	    char* pEnd;
       	    char* pStop;
	    parIni =  strtod(line.c_str(), &pEnd);
	    pStop=pEnd+3;
    	    errIni =  strtod(pStop, NULL);
// 	                 if (errIni==0.00){
// 			   ParMaxTemp=0.0;
// 			 }else{
// 			   ParMaxTemp=ParMax;
// 			 }  
     			 parLis.emplace_back(varName,varName, parIni,ParMin,ParMax);
     			 errLis.emplace_back(errName,errName, errIni,ParMin,ParMax);     
			 if (errIni==0.00) parLis[i].setConstant(kTRUE);
//
     			 if(parIni==1 || parIni==0.00){
			  parFit.emplace_back(varName,varName, parIni,ParMin,ParMax);
     			  errFit.emplace_back(errName,errName, errIni,ParMin,ParMax); 
			 }else{ 
			  parFit.emplace_back(varName,varName, 0.001,ParMin,ParMax);
			 }     
			 if (errIni==0.00) parFit[i].setConstant(kTRUE);

			 std::cout<<Form("Coeff (%s) = %f+/-%f   => %s",varName,parLis[i].getVal(),errLis[i].getVal(), ParCheck[i])<<std::endl;
     			 if(fabs(parIni)>0.0 ) NumParamFree ++;
//     			if(fabs(parIni)>0.0009999 && parIni!=1.0000) NumParamFree ++;
     }
    parListInput->close();
  }else{
     std::cout<<"List of initial parameters "<< ListParName <<" not found!"<<std::endl;
     exit(1);
  }   
  for (int i=0;i<numParameters;++i){
   coefLis->add(parLis[i]);
   coefFit->add(parFit[i]);
  }	      
  BernSideBand    = new RooBernsteinSideband(Form("BernSideBand_bin%d_%s",Q2Bin,RunEra),Form("BernSideBand_bin%d_%s",Q2Bin,RunEra),*ctL,*ctK,*phi,*coefLis,maxDegree1,maxDegree2,maxDegree3);
  BernSideFits    = new RooBernsteinSideband("BernSideFits","BernSideFits",*ctL,*ctK,*phi,*coefFit,maxDegree1,maxDegree2,maxDegree3);
//
//   phi->setRange("redu",-3.14159, 3.14159 );
//   RooAbsReal* BernCheckNorm = BernSideBand->createIntegral(RooArgSet(*ctL,*ctK,*phi),RooArgSet(*ctL,*ctK,*phi));
//   RooAbsReal* BernCheckNormRedu = BernSideBand->createIntegral(RooArgSet(*ctL,*ctK,*phi),RooArgSet(*ctL,*ctK,*phi),"redu");
//   std::cout<<Form("Check Normalization for Angular Bkg Model (Bernstein) = %1.20f (redu=%1.20f)\n",BernCheckNorm->getVal(),BernCheckNormRedu->getVal())<<std::endl;
//
  if(GenMiniMC) FirstMC=true;
  B0Sigma = FitMassSpectrumRoofit(fulldata, c2, HxMassQ2,pdfHxMassQ2,sigHxMassQ2,bkgHxMassQ2, maxDegree);
//  
//  exit(1);
  std::cout<< "Setting B0Sigma = "<<B0Sigma<<" from the fit to the mass spectrum\n"<<std::endl; 
//  c2->Write();
  gSystem->Exec(Form("mv %s %s.tmp",PNGNameMassQ2Hist,PNGNameMassQ2Hist));
  c2->Print(PNGNameMassQ2Hist);
//   XMinSBL = B0Mass - NSigma1L*B0Sigma;
//   XMaxSBL = B0Mass - NSigma2L*B0Sigma;
//   XMinSBR = B0Mass + NSigma1R*B0Sigma;
//   XMaxSBR = B0Mass + NSigma2R*B0Sigma;
  std::cout<<" XMinSBL  			 :     "<<XMinSBL<<std::endl;
  std::cout<<" XMaxSBL  			 :     "<<XMaxSBL<<std::endl;
  std::cout<<" XMinSBR  			 :     "<<XMinSBR<<std::endl;
  std::cout<<" XMaxSBR  			 :     "<<XMaxSBR<<std::endl;
  std::cout<<Form("(tagged_mass>%2.40f&&tagged_mass<%2.40f)||(tagged_mass>%2.40f&&tagged_mass<%2.40f)",XMinSBL,XMaxSBL,XMinSBR,XMaxSBR) <<std::endl;
  
  std::vector<double> CorreAdaptX;
  std::vector<double> CorreAdaptY;
  std::vector<double> CorreAdaptZ;
  
  TH1D* HWMass         = new TH1D( "HWMass"     , "B^{0} Mass"		 ,	      xMassHBin2, XMinSign, XMaxSign);
  TCanvas* cw = new TCanvas("cw","Fit Mass Spectrum",200,10,900,780);
//  RooDataSet* dataSBMass  = (RooDataSet*)fulldata->reduce(RooArgSet(*tagged_mass,*mumuMass,*mumuMassE),"xtagged_mass>XMinSBL&&xtagged_mass<XMaxSBL)||xtagged_mass>XMinSBR&&xtagged_mass<XMaxSBR") ;
  RooDataSet *dataSBMass   = new RooDataSet("dataSBMass", "dataSBMass",  RooArgSet(*tagged_mass,*mumuMass,*mumuMassE));
  RooDataSet *dataSBAng    = new RooDataSet("dataSBAng", "dataSBAng",  RooArgSet(*ctL,*ctK,*phi));
  for (Int_t i=0;i<nentries;i++) {
  	  RecoB0TreeOut->GetEntry(i);
//  	 if(mmk2>3.6&&mmk2<4.2&&mmk1>4.7&&mmk1<4.9&&Q2Bin==4)continue;
	  if(cos_theta_l==-99.) continue;
          if( (xtagged_mass>XMinSBL&&xtagged_mass<XMaxSBL)|| 
	      (xtagged_mass>XMinSBR&&xtagged_mass<XMaxSBR)){
// 	  if((mmk2< 4)) continue;
// 	   if(cos_theta_l==-99.){
// 	    cout<<Form("99 tagged_mass=%f %f %f %f",xtagged_mass,cos_theta_l,cos_theta_k,phi_kst_mumu)<<endl;
// 	    continue;
// 	   } ;
//  	   if(cos_theta_k<0.6||cos_theta_k>0.8) continue;
//  	   if(xtagged_mass>5.2) continue;
	   HWMass->Fill(xtagged_mass);
	   HxMassQ2SB->Fill(xtagged_mass);
	   HxReco->Fill(cos_theta_l,cos_theta_k,phi_kst_mumu);
	   HxMassVsCosL->Fill(xtagged_mass,cos_theta_l);
	   HxMassVsCosK->Fill(xtagged_mass,cos_theta_k);
	   HxMassVsPhi ->Fill(xtagged_mass,phi_kst_mumu);
 	   dataSBMass->add(RooArgSet(*tagged_mass,*mumuMass,*mumuMassE));
	   ctL->setVal(cos_theta_l);
	   ctK->setVal(cos_theta_k);
	   phi->setVal(phi_kst_mumu);
 	   dataSBAng->add(RooArgSet(*ctL,*ctK,*phi));
//	   std::cout<<xL<<" "<<yK<<" "<<zP<<std::endl;
	   CorreAdaptX.push_back(cos_theta_l);
	   CorreAdaptY.push_back(cos_theta_k);
	   CorreAdaptZ.push_back(phi_kst_mumu);
//          }else{
// 	  if(cos_theta_l!=-99.) cout<<Form("tagged_mass=%f %f %f %f",xtagged_mass,cos_theta_l,cos_theta_k,phi_kst_mumu)<<endl;
	 }
  }   
  
  cw->cd();
  HWMass->Draw();
  cw->Print("hwmass.png");
  int GenEntries = NFactGen*dataSBAng->sumEntries();
  std::cout<<"Found SB entries = "<<HxReco->GetEntries()<<std::endl;
  if(HxReco->GetEntries()<10){
   std::cout<<"Error!! too few SB entries for a Fit: SB entries =  "<<HxReco->GetEntries()<<" EXIT!!!"<<std::endl;
   exit(0);
  }
//  double xBinw =  HxReco->GetXaxis()->GetBinWidth(1) ;
//  double yBinw =  HxReco->GetYaxis()->GetBinWidth(1) ;
//  double zBinw =  HxReco->GetZaxis()->GetBinWidth(1) ;
//   TH1D* HxRecoCosL = (TH1D*) HxReco->Project3D("x");
//   TH1D* HxRecoCosK = (TH1D*) HxReco->Project3D("y");
//   TH1D* HxRecoPhi  = (TH1D*) HxReco->Project3D("z");
// 
//  TH2D* HxRecoCosLK =(TH2D*) HxReco->Project3D("xy");
  OutFile->cd();
    
  
  TH3D* HSBFunc 	  = new TH3D( "HSBFunc"          , "HSBFunc",		NFact*xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
 										NFact*xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
										NFact*xPhiHBin , XMinPhi, XMaxPhi );

//   TH3D* HSideBandRecoTest = new TH3D( "HSideBandRecoTest", "HSideBandRecoTest", NFact*xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
//  										NFact*xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
// 										NFact*xPhiHBin , XMinPhi, XMaxPhi );
//   totalParams=0;

  RooDataSet *dataSBMass_Plot   = new RooDataSet("dataSBMass_Plot", "dataSBMass_Plot",  RooArgSet(*ctL,*ctK,*phi));
  if (ReFit){ 
    RooAbsReal* nll = BernSideFits->createNLL(*dataSBAng,RooFit::NumCPU(NCPU),RooFit::Verbose(kTRUE)); 
    RooMinuit Minuit(*nll) ;
//    P_3.setConstant(kTRUE);
//    P4p.setConstant(kTRUE);
//    P5p.setConstant(kTRUE);
//    P6p.setConstant(kTRUE);
//    P8p.setConstant(kTRUE);
//     parLis[4].setConstant(kTRUE); 
    Minuit.migrad() ;
    Minuit.hesse() ;
    for (int i=0;i<numParameters;++i){
     if(parFit[i].getVal()!=0.0){
      double coeffy = parFit[i].getVal();
      double errory = parFit[i].getError();
      std::cout<<Form("RESULTS ==>  p(%d)=%f+/-%f rate=%f => [%s]",i,coeffy,errory,coeffy/errory,ParCheck[i])<<std::endl;
     }
    }

   double SBEntries= HxReco->GetEntries();
   TCanvas* cSideRooRefit = new TCanvas("cSideRooRefit","Sideband Projections",200,10,900,780);
   cSideRooRefit->Divide(2,2);
   RooPlot* xframeFits=ctL->frame(RooFit::Bins(xCosLHBin));
   dataSBAng->plotOn(xframeFits, RooFit::MarkerStyle(kPlus));
   BernSideFits->plotOn( xframeFits,RooFit::LineColor(kRed),Normalization(SBEntries,RooAbsReal::NumEvent),RooFit::ProjWData(RooArgSet(*ctL),*dataSBAng));
   cSideRooRefit->cd(1);
   xframeFits->Draw();
   RooPlot* yframeFits=ctK->frame(RooFit::Bins(xCosKHBin));
   dataSBAng->plotOn(yframeFits, RooFit::MarkerStyle(kPlus));
   BernSideFits->plotOn( yframeFits,RooFit::LineColor(kRed),Normalization(SBEntries,RooAbsReal::NumEvent),RooFit::ProjWData(RooArgSet(*ctK),*dataSBAng));
   cSideRooRefit->cd(2);
   yframeFits->Draw();
   RooPlot* zframeFits=phi->frame(RooFit::Bins(xPhiHBin));
   dataSBAng->plotOn(zframeFits, RooFit::MarkerStyle(kPlus));
   BernSideFits->plotOn( zframeFits,RooFit::LineColor(kRed),Normalization(SBEntries,RooAbsReal::NumEvent),RooFit::ProjWData(RooArgSet(*phi),*dataSBAng));
   cSideRooRefit->cd(3);
   zframeFits->Draw();
   gSystem->Exec(Form("mv %s %s.tmp",PNGNameReFitSB3D,PNGNameReFitSB3D));
   cSideRooRefit->Print(PNGNameReFitSB3D);
   cSideRooRefit->Write();
  }   
   
  if (MCStudy){
//      BernSideBand->defaultGeneratorConfig()->methodND(kFALSE,kFALSE).setLabel("RooAcceptReject") ;
//      BernSideBand->defaultGeneratorConfig()->getConfigSection("RooAcceptReject").setRealValue("nTrial0D",1e9);
//      BernSideBand->defaultGeneratorConfig()->getConfigSection("RooAcceptReject").setRealValue("nTrial1D",1e9);
//      BernSideBand->defaultGeneratorConfig()->getConfigSection("RooAcceptReject").setRealValue("nTrial3D",1e11);
//      BernSideBand->defaultGeneratorConfig()->getConfigSection("RooAcceptReject").Print("V");

   BernSideBand->defaultGeneratorConfig()->methodND(kFALSE,kFALSE).setLabel("RooFoamGenerator") ;
   BernSideBand->defaultGeneratorConfig()->getConfigSection("RooFoamGenerator").setRealValue("chatLevel",1);
   BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").setRealValue("nSample",50000);
   BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").setRealValue("nCell3D",100000);
   BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").setRealValue("nCell1D",100000);
   BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").setRealValue("EvPerBin",100);
   TCanvas* cmc = new TCanvas("cmc","Pulls",900,900) ;
   TCanvas* cmp = new TCanvas("cmp","Parameters",900,900) ;
   TCanvas* cme = new TCanvas("cme","Errors",900,900) ;
   cmc->Divide(2,3) ;
   cmp->Divide(2,3) ;
   cme->Divide(2,3) ;
//
//   RooPolynomial* Poly = new RooPolynomial("Poly","Poly",*phi,*coefLis);
//   Bern = new RooBernstein("Bern","Bern",*phi,*coefLis);
//   Bern = new RooBernstein("Bern","Bern",*phi,*coefLis);
//   RooRealVar  NSample("NSample","number of events",GenEntries,0,2*GenEntries) ; 
//   RooExtendPdf BernExt("BernEx","extended Bernstein PDF",*Bern,NSample) ; 

   Bern    = new RooBernsteinSideband("Bern","BernSideBand MiniMC",*ctL,*ctK,*phi,*coefLis,maxDegree1,maxDegree2,maxDegree3); 
//   
   std::vector<RooPlot*> rootPullPlots;
   std::vector<RooPlot*> rootParmPlots;
   std::vector<RooPlot*> rootErroPlots;
//    RooMCStudy* mcstudy = new RooMCStudy(*BernSideBand,RooArgSet(*ctL,*ctK,*phi),RooFit::NumCPU(NCPU),Binned(kFALSE),Silence(),
//                      FitOptions(Save(kTRUE),PrintEvalErrors(1))) ;
//   RooMCStudy* mcstudy = new RooMCStudy(*Poly,RooArgSet(*phi),Binned(kFALSE),Extended(kFALSE),FitOptions(Save(kTRUE),PrintEvalErrors(1),NumCPU(NCPU))) ;
//    RooMCStudy* mcstudy = new RooMCStudy(*BernSideBand,RooArgSet(*ctL,*ctK,*phi),Binned(kFALSE),
//                       FitOptions(Save(kTRUE),PrintEvalErrors(1),NumCPU(NCPU))) ;

//ultima
//   RooMCStudy* mcstudy = new RooMCStudy(*Bern,RooArgSet(*phi),Binned(kFALSE),
//                      FitOptions(Save(kTRUE),PrintEvalErrors(1),NumCPU(NCPU))) ;
//
   RooMCStudy* mcstudy = new RooMCStudy(*Bern,RooArgSet(*ctL,*ctK,*phi),Binned(kFALSE),
                      FitOptions(Save(kTRUE),PrintEvalErrors(1),NumCPU(NCPU))) ;
		      
		      
   mcstudy->generateAndFit(NSampleMC,GenEntries);
   int ipp=0;
   for (int i=0;i<numParameters;++i){
    if(parLis[i].getVal()!=0.0&&parLis[i].getVal()!=1.0){
      rootPullPlots.emplace_back(mcstudy->plotPull(parLis[i],Bins(40),FitGauss(kTRUE))) ;
      rootParmPlots.emplace_back(mcstudy->plotParam(parLis[i],Bins(40))) ;
      rootErroPlots.emplace_back(mcstudy->plotError(parLis[i],Bins(40))) ;
      cmc->cd(ipp+1);gPad->SetLeftMargin(0.15) ; 
      rootPullPlots[ipp]->GetYaxis()->SetTitleOffset(1.4) ;
      rootPullPlots[ipp]->Draw();
      cmp->cd(ipp+1);gPad->SetLeftMargin(0.15) ; 
      rootParmPlots[ipp]->GetYaxis()->SetTitleOffset(1.4) ;
      rootParmPlots[ipp]->Draw();
      cme->cd(ipp+1);gPad->SetLeftMargin(0.15) ; 
      rootErroPlots[ipp]->GetYaxis()->SetTitleOffset(1.4) ;
      rootErroPlots[ipp]->Draw();
      ipp++; 
    }
   } 
   gSystem->Exec(Form("mv %s %s.tmp",PNGNameStudyMCPulls,PNGNameStudyMCPulls));  
   cmc->Print(PNGNameStudyMCPulls); 
   gSystem->Exec(Form("mv %s %s.tmp",PNGNameStudyMCParam,PNGNameStudyMCParam));  
   cmp->Print(PNGNameStudyMCParam); 
   gSystem->Exec(Form("mv %s %s.tmp",PNGNameStudyMCError,PNGNameStudyMCError));  
   cme->Print(PNGNameStudyMCError); 
  }
  if (GenMiniMC){
   std::vector<TH1D*> HParFit; 
   char varName[100];
   for (int i=0;i<numParameters;++i){
    if(parLis[i].getVal()!=0.000&&parLis[i].getVal()!=1.000){
     sprintf(varName, "p%d", i);
//     HParVec.emplace_back(varName,varName,100.,-5.,5.);
     HParFit.push_back( new TH1D (Form("HParFit_p%d",i),Form("p%d",i),100.,-5.,5.) );
    }
   }  
//     BernSideBand->defaultGeneratorConfig()->methodND(kFALSE,kFALSE).setLabel("RooAcceptReject") ;
//     BernSideBand->defaultGeneratorConfig()->getConfigSection("RooAcceptReject").setRealValue("nTrial3D",1e11);

   BernSideBand->defaultGeneratorConfig()->methodND(kFALSE,kFALSE).setLabel("RooFoamGenerator") ;
   BernSideBand->defaultGeneratorConfig()->getConfigSection("RooFoamGenerator").setRealValue("chatLevel",2);


//   RooFoamGenerator * NumGenerator = (RooFoamGenerator *)(BernSideBand->specialGeneratorConfig(kTRUE)); 
//   NumGenerator->Print("2");
//   BernSideBand->defaultGeneratorConfig()->getConfigSection("RooFoamGenerator").setRealValue("OptDrive",2);
//   BernSideBand->defaultGeneratorConfig()->methodND(kFALSE,kFALSE).setLabel("RooAcceptReject") ;
   //BernSideBand->defaultGeneratorConfig()->getConfigSection("RooAcceptReject").Print("V");
// BernSideBand->specialGeneratorConfig(kTRUE)->methodND(kFALSE,kFALSE).setLabel("RooAcceptReject") ;
// BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooAcceptReject").Print("V"); 
// BernSideBand->specialGeneratorConfig(kTRUE)->methodND(kFALSE,kFALSE).setLabel("RooFoamGenerator") ;
// BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").setRealValue("chatLevel",1);
//   BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").setRealValue("OptRej",0);
//   BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").setRealValue("MaxWtRej",1.6);
// BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").Print("V");
  BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").setRealValue("nSample",50000);
  BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").setRealValue("nCell3D",100000);
  BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").setRealValue("EvPerBin",25);
// //BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").setRealValue("nCellND",100000);
// BernSideBand->specialGeneratorConfig(kTRUE)->getConfigSection("RooFoamGenerator").Print("V");
  
   int SEED0 = 1234+NIniGen+NSampleMC;
   std::cout<<Form("==> SEED0  %d", SEED0) <<std::endl;
   RooRandom::randomGenerator()->SetSeed(SEED0);
//   RooAbsPdf::GenSpec* genSpec = BernSideBand->prepareMultiGen(RooArgSet(*ctL,*ctK,*phi),ProtoData(*dataSBAng)) ;
   RooAbsPdf::GenSpec* genSpec = BernSideBand->prepareMultiGen(RooArgSet(*ctL,*ctK,*phi),NumEvents(GenEntries)) ;
//   RooAbsPdf::GenSpec* genSpec = BernSideBand->prepareMultiGen(RooArgSet(*ctL,*ctK,*phi),NumEvents(NumMassNewGen)) ;
//   RooAbsPdf::GenSpec* genSpec = BernSideBand->prepareMultiGen(RooArgSet(*ctL,*ctK,*phi),NumEvents(1000),ProtoData(*dataSBMass)) ;
   for (int iLoop=NIniGen;iLoop<=NSampleMC;iLoop++){
    std::cout<<Form("==== Gen MiniMC Loop %d",iLoop) <<std::endl;
//    gRandom = new TRandom3(0);
    int SEED = iLoop*1000+iLoop;
    std::cout<<Form("==> SEED Loop %d %d",iLoop, SEED) <<std::endl;
    RooRandom::randomGenerator()->SetSeed(SEED);
    char OutNameloop[310]="";
    char PNGNameMassQ2Hist20[310]="";
    sprintf(OutNameloop,Form("%s_%d",OutFileNameMiniMCHisto,iLoop));
    TFile* OutFileMiniMCHisto = TFile::Open(OutNameloop,"RECREATE");
    TCanvas* c20 = new TCanvas("c20","Fit Mass Spectrum MC",200,10,900,780);
        
//    TH1D* HxMassQ2MC       = new TH1D( "HxMassQ2MC"	, "B^{0} Mass MC"		  ,	       xMassHBin2, XMinSign, XMaxSign);
    TH1D* HxMassQ2MC = (TH1D*)geneMassNew->createHistogram("B^{0} Mass MC",*tagged_mass,xMassHBin2);
    TH1D* pdfHxMassQ2MC    = new TH1D( "pdfHxMassQ2MC", "B^{0} Mass Fit Q2 Bin MC",  XHScale * xMassHBin , XMinSign, XMaxSign);
    TH1D* sigHxMassQ2MC    = new TH1D( "sigHxMassQ2MC", "B^{0} Mass Fit Q2 Bin MC",  XHScale * xMassHBin , XMinSign, XMaxSign);
    TH1D* bkgHxMassQ2MC    = new TH1D( "bkgHxMassQ2MC", "B^{0} Mass Fit Q2 Bin MC",  XHScale * xMassHBin , XMinSign, XMaxSign);
    FirstMC=false;
    double B0SigmaMC = FitMassSpectrumRoofit(geneMassNew, c20, HxMassQ2MC,pdfHxMassQ2MC,sigHxMassQ2MC,bkgHxMassQ2MC, maxDegree);
    std::cout<<Form("==== Gen MiniMC Loop %d",iLoop) <<std::endl;
    sprintf(PNGNameMassQ2Hist20,Form("%s.MiniMC-NFactGen%d_%d",PNGNameMassQ2Hist,NFactGen,iLoop));
    
    c20->Print(PNGNameMassQ2Hist20);
    XMinSBL = B0Mass - NSigma1L*B0SigmaMC;
    XMaxSBL = B0Mass - NSigma2L*B0SigmaMC;
    XMinSBR = B0Mass + NSigma1R*B0SigmaMC;
    XMaxSBR = B0Mass + NSigma2R*B0SigmaMC;
    int NumGeneNew = \
    geneMassNew->sumEntries();
    int NumGeneNewSB =geneMassNew->sumEntries(Form("(tagged_mass>%f&&tagged_mass<%40f)||(tagged_mass>%2.40f&&tagged_mass<%2.40f)",XMinSBL,XMaxSBL,XMinSBR,XMaxSBR));
    std::cout<<Form("Start to generate %d sideband events ...\n",NumGeneNew) <<std::endl;
    std::cout<<Form("Start to generate %d sideband events inside range:\n",NumGeneNewSB) <<std::endl;
    std::cout<<Form("(tagged_mass>%2.40f&&tagged_mass<%2.40f)||(tagged_mass>%2.40f&&tagged_mass<%2.40f)",XMinSBL,XMaxSBL,XMinSBR,XMaxSBR) <<std::endl;
// ctL->setBins(1000);
// ctK->setBins(1000);
// phi->setBins(1000);
//     const RooNumGenConfig * MCGen = BernSideBand->getGeneratorConfig();
//     const RooCategory ND = MCGen->methodND(kFALSE, kFALSE);
// //     RooArgSet MCGenConf1 = (MCGen->getConfigSection("RooAcceptReject"));
// //     MCGenConf1.setRealValue("nTrialND", 20000);
//     RooArgSet MCGenConf2 = (MCGen->getConfigSection("RooFoamGenerator"));
//     MCGenConf2.setRealValue("chatLevel", 2);
//     MCGenConf2.setRealValue("nCells", 10000);



//    MCGen->methodND(kFALSE, kFALSE);
//    MCGen->Print("V");

    //RooAbsPdf::defaultGeneratorConfig()->methodND(kFALSE, kFALSE).setLabel("RooAcceptReject");
    geneDataNew = BernSideBand->generate(*genSpec) ;
    for (int i=0;i<numParameters;++i){
     if(parLis[i].getVal()!=0.0){
      parFit[i].setVal(parLis[i].getVal());
     } 
    }
    
    RooAbsReal* nll = BernSideFits->createNLL(*geneDataNew,RooFit::NumCPU(NCPU),RooFit::Verbose(kTRUE)); 
    RooMinuit Minuit(*nll) ;
    Minuit.migrad() ;
    Minuit.minos() ;
    int ipp=0;
    for (int i=0;i<numParameters;++i){
     if(parFit[i].getVal()!=0.0&&parFit[i].getVal()!=1.0){
      double coeffy = parFit[i].getVal();
      double errory = parFit[i].getError();
      std::cout<<Form("RESULTS ==>  p(%d)=%f+/-%f rate=%f => [%s] iLoop=%d",i,coeffy,errory,coeffy/errory,ParCheck[i],iLoop)<<std::endl;
      HParFit[ipp]->Fill((parLis[i].getVal()-coeffy)/errory);
      ipp++;
     }
    }
//    geneDataNew = BernSideBand->generate(RooArgSet(*ctL,*ctK,*phi),*dataSBMass,NumGeneNew) ;
//    geneDataNew = BernSideBand->generate(RooArgSet(*ctL,*ctK,*phi),NumGeneNew, Verbose(kTRUE)) ;
//    geneDataNew = BernSideBand->generate(RooArgSet(*ctL,*ctK,*phi),NumGeneNew, Verbose(kTRUE),Extended(kTRUE)) ;
    std::cout<<"End generate MiniMC\n" <<std::endl;
    GenMiniMCB0TreeOut = makeGenSBTTree(geneDataNew,geneMassNew,BernSideBand,OutFileMiniMCHisto);
    OutFileMiniMCHisto->cd();
    GenMiniMCB0TreeOut->Write();
    HxMass->Write();
    HxMassQ2MC->Write();
    pdfHxMassQ2MC->Write();
    sigHxMassQ2MC->Write();
    bkgHxMassQ2MC->Write();
    c20->Write();
    std::cout<<Form("Write generated Ntupla in file: %s\n",OutFileNameMiniMCHisto) <<std::endl;
    std::cout<<"==================="<<std::endl;
    std::cout<<"==     END	==="<<std::endl;
    std::cout<<"==================="<<std::endl;
    delete geneDataNew;
    delete geneMassNew;
    delete GenMiniMCB0TreeOut;
    delete nll;
    if (iLoop<NSampleMC){
     FirstMC=true;
     B0Sigma = FitMassSpectrumRoofit(fulldata, c2, HxMassQ2,pdfHxMassQ2,sigHxMassQ2,bkgHxMassQ2, maxDegree);
    }
    OutFileMiniMCHisto->Close();
    OutFile->cd();
   } 
   int ipp=0;
   for (int i=0;i<numParameters;++i){
    if(parFit[i].getVal()!=0.0&&parFit[i].getVal()!=1.0){
     HParFit[ipp]->Write();
     HParFit[ipp]->Fit("gaus");
     ipp++;
    }
   }
  }else if(AnaMiniMC){
   TH1D* HxChi2NDOF_reduced	= new TH1D( "HxChi2NDOF_reduced"     , "Chi2/NDOF [reduced]",	     100, 0.,  4.);
   TH1D* HxChi2NDOF		= new TH1D( "HxChi2NDOF"    , "Chi2/NDOF",	     100, 0.,  4.);
   TH1D* HxPValMin		= new TH1D( "HxPValMin"     , "P-Value Min Estimation",       50, 0.,  1.);
   TH1D* HxPValMax		= new TH1D( "HxPValMax"     , "P-Value Max Estimation",       50, 0.,  1.);

   TH1D* HxChi2NDOF_reduced_diff= new TH1D( "HxChi2NDOF_reduced [model diff]"     , "Chi2/NDOF [reduced]",	     100, 0.,  4.);
   TH1D* HxChi2NDOF_diff	= new TH1D( "HxChi2NDOF[model diff]"    , "Chi2/NDOF",	     100, 0.,  4.);
   TH1D* HxPValMin_diff		= new TH1D( "HxPValMin [model diff]"     , "P-Value Min Estimation",       50, 0.,  1.);
   TH1D* HxPValMax_diff		= new TH1D( "HxPValMax [model diff]"     , "P-Value Max Estimation",       50, 0.,  1.);
//   TF1 *fit_pull = new TF1("fit_pull","TMath::Gaus([0]*x + [1])", 0.01, 0.1);
   std::vector<TH1D*> HParVec; 
   std::vector<TH1D*> HParModErrVec; 
   std::vector<TH1D*> HSignifVec; 
   std::vector<double> PullVec; 
   std::vector<double> PullModErrVec; 
   std::vector<double> SignifVec; 
   char varName[100];
   int iParFree=0;
   for (int i=0;i<numParameters;++i){
    if(parLis[i].getVal()!=0.000&&parLis[i].getVal()!=1.000){
     sprintf(varName, "p%d", i);
//     HParVec.emplace_back(varName,varName,100.,-5.,5.);
     HParVec.push_back( new TH1D (varName,varName,100.,-5.,5.) );
     sprintf(varName, "p%d_ModErr", i);
     HParModErrVec.push_back( new TH1D (varName,varName,100.,-5.,5.) );
     sprintf(varName, "p%d_Signif", i);
     HSignifVec.push_back( new TH1D (varName,varName,150.,0.,15.) );
     PullVec.push_back(-100000);
     PullModErrVec.push_back(-100000);
     SignifVec.push_back(-100000);
     iParFree++;
    }
   }  
   char ListParNameLoop[650]="";
   int iCountError=0;
   int iModelError=0;
   int iCountFit=0;
   std::string gline;
//   std::string spvalue;
   std::string gkeywordPMin("PD P-Value		Min	[Adapt] = ");
   std::string gkeywordPMax("PD P-Value		Max	[Adapt] = ");
   std::string gkeywordChi2NDOF("PowerDivergence/NDOF");
   std::string gkeywordChi2NDOFR("PowerDivergence/NDOF_Reduced");
   for (int iLoop=NIniGen;iLoop<=NSampleMC;iLoop++){
    std::string gfile(Form("%s/%s_%d",AnaMiniMCDir,LogFileName,iLoop));
    bool LoopError = false;
    bool ModError  = false;
    sprintf(ListParNameLoop,Form("%s/%s.%s_%d",AnaMiniMCDir,ListParName,MiniMCTXT,iLoop));
//    char varName[100];
    int NumParamFree = 0;
    double parRead=0.;
    double errRead=0.;
//    bool FoundParamNotZero = false;
    std::string line;
    std::cout<<"Try to open list of initial parameters :"<< ListParNameLoop<<std::endl;
    std::fstream *parListInputLoop = new std::fstream(ListParNameLoop,std::ifstream::in);
    if(parListInputLoop->is_open()){
       std::cout<<"List of initial parameters :"<< ListParNameLoop <<" FOUND!!!"<<std::endl;
       int ii=0;
       std::string gline1=grep(gfile,gkeywordPMin);
       std::string spvalueMin(gline1.begin()+26,gline1.begin()+38);
       double pvalueMin = atof(spvalueMin.c_str());
       std::string gline2=grep(gfile,gkeywordPMax);
       std::string spvalueMax(gline2.begin()+26,gline2.begin()+38);
       double pvalueMax = atof(spvalueMax.c_str());
       std::string gline3=grep(gfile,gkeywordChi2NDOF);
       std::string sChi2NDOF(gline3.begin()+31,gline3.begin()+43);
//       std::cout<<Form("gline3=%s  %s\n ",gline3.c_str(),sChi2NDOF.c_str());exit(1);
       double Chi2NDOF = atof(sChi2NDOF.c_str());
       std::string gline4=grep(gfile,gkeywordChi2NDOFR);
       std::string sChi2NDOFR(gline4.begin()+39,gline4.begin()+53);
       double Chi2NDOFR = atof(sChi2NDOFR.c_str());
//       std::cout<<Form("gline4=%s  %s\n ",gline4.c_str(),sChi2NDOFR.c_str());exit(1);
       std::cout<<Form("Chi2/P-Value iLoop=%d  pvalueMin=%f pvalueMax=%f Chi2/NDOF=%f Chi2/NDOF_R=%f\n"\
       ,iLoop,pvalueMin,pvalueMax,Chi2NDOF,Chi2NDOFR);
       for (int i=0;i<numParameters;++i){
    	 sprintf(varName, "p%d", i);
    	 std::getline(*parListInputLoop, line);
    	 char* pEnd;
    	 char* pStop;
    	 parRead =  strtod(line.c_str(), &pEnd);
	 pStop=pEnd+3;
    	 errRead =  strtod(pStop, NULL);
	 double rate = parRead/errRead;
//	 std::cout<<"LOG file Fit: "<<gfile<<"\n"<<std::endl;
         if(parLis[i].getVal()!=0.0000&&parLis[i].getVal()!=1.0000){
	  PullVec[ii] = (parLis[i].getVal()-parRead)/errRead;
	  SignifVec[ii] = rate;
	  if (parRead==0.0000||parRead>=1.000){
           std::cout<<Form("Error iLoop=%d  p(%d)=%f+/-%f [ref %f sign %f] [%s]!!!",
	   iLoop,i,parRead,errRead,parLis[i].getVal(),parLis[i].getVal()/errLis[i].getVal(),ParCheck[i])<<std::endl;
	   LoopError=true;
	  }else{ 
           std::cout<<Form("Pull=%f iLoop=%d  rate[p(%d)/err(%d)]=%f !!! ",PullVec[ii],iLoop,i,i,rate)<<std::endl;
	  } 
//	  HParVec[ii]->SetStats(true);
//  	  if(!LoopError){
// 	   HParVec[ii]->Fill(pull);
// 	   HParVec[ii]->Fit("gaus");
// 	  }else{
// 	   HParModErrVec[ii]->Fill(pull);
// 	  } 
          ii++;
	 }else{ 
	  if (parRead!=0.0000&&parRead!=1.0000){
	    Hparams->Fill(i);
	    std::cout<<Form("Error iLoop=%d  p(%d)=%f+/-%f rate[p(%d)/err(%d)]=%f [ref %f]  [%s]!!! p-value Min=%f",
	    iLoop,i,parRead,errRead,i,i,rate,parLis[i].getVal(),ParCheck[i],pvalueMin)<<std::endl;
	    LoopError=true;
	    if(parLis[i].getVal()==0.0000) ModError =true;
	  }  
	 } 
     	 if(fabs(parRead)>0.0 ) NumParamFree ++;
	 std::cout<<Form("%s = %f +/- %f",varName,parRead,errRead)<<std::endl;
       }
       iCountFit++;
       for (int ip=0;ip<iParFree;++ip){
        if(LoopError){
         HParModErrVec[ip]->Fill(PullVec[ip]);
        }else{  
         HParVec[ip]->Fill(PullVec[ip]);
         HSignifVec[ip]->Fill(SignifVec[ip]);
	} 
       } 
       if(LoopError) {
	HxPValMin_diff->Fill(pvalueMin);
	HxPValMax_diff->Fill(pvalueMax);
	HxChi2NDOF_diff ->Fill(Chi2NDOF);
	HxChi2NDOF_reduced_diff ->Fill(Chi2NDOFR);
       }else{	 
	HxPValMin ->Fill(pvalueMin);
	HxPValMax ->Fill(pvalueMax);
	HxChi2NDOF ->Fill(Chi2NDOF);
	HxChi2NDOF_reduced ->Fill(Chi2NDOFR);
       }	
       if(LoopError) iCountError++;
       if(ModError ) iModelError++;
    parListInputLoop->close();
    }else{
     std::cout<<"List :"<<Form("%s", ListParNameLoop) <<" Not Found!!!\n"<<std::endl;
    }
   }
   std::cout<<Form("=====> Total count for Errors=%d ModelError=%d [N. Fits=%d] [rateErr=%f]  [rateModelDiff=%f]",iCountError,iModelError,iCountFit,float(iCountError)/float(iCountFit),float(iModelError)/float(iCountFit))<<std::endl;
   gStyle->SetStatW(.2);
   gStyle->SetStatFont(63);
   gStyle->SetStatFontSize(11);
   int Den =2;
   int NumParamFree1=NumParamFree;
   if (NumParamFree1>8) Den =3;
   int Div = round(NumParamFree1/Den);
   if(NumParamFree1%Den) Den= Den+1;
   int XWin=800;
   if(NumParamFree>3) XWin=1600;
   
   TCanvas* cp = new TCanvas("cp","Fit Coeff. Pulls",200,200,XWin,800);
   TCanvas* dp = new TCanvas("dp","Fit Coeff. Pulls [!= Gen model]",200,200,XWin,800);
   TCanvas* xp = new TCanvas("xp","Fit Coeff. not Gen.",200,200,1600,800);
   cp->Divide(Div,Den);   
   dp->Divide(Div,Den);   
   int ii=0 ; 
   for (int i=0;i<numParameters;++i){
    if(parLis[i].getVal()!=0.000&&parLis[i].getVal()!=1.000){
     cp->cd(ii+1);
     HParVec[ii]->Fit("gaus");
     HParVec[ii]->Draw();
     cp->Update();
     dp->cd(ii+1);
     HParModErrVec[ii]->Fit("gaus");
     HParModErrVec[ii]->Draw();
     dp->Update();
     OutFile->cd();
     HParVec[ii]->Write();
     HParModErrVec[ii]->Write();
     HSignifVec[ii]->Write();
     ii++;
    }
   } 
   TCanvas* ch = new TCanvas("ch","Fit Chi2 & PVal",200,200,1600,800);
   ch->Divide(4,2);   
   ch->cd(1);
   gStyle->SetStatW(.5);
   gStyle->SetStatH(.08);
   gStyle->SetOptStat("MNRe");
   HxPValMin ->Draw();
   ch->cd(2);
   HxPValMax ->Draw();
   ch->cd(3);
   HxChi2NDOF ->Draw();
   ch->cd(4);
   HxChi2NDOF_reduced ->Draw();
   ch->cd(5);
   HxPValMin_diff->Draw();
   ch->cd(6);
   HxPValMax_diff->Draw();
   ch->cd(7);
   HxChi2NDOF_diff ->Draw();
   ch->cd(8);
   HxChi2NDOF_reduced_diff ->Draw();
    
   xp->cd();
   Hparams->Draw();
   OutFile->cd();
   HxPValMin ->Write();
   HxPValMax ->Write();
   HxPValMin_diff->Write();
   HxPValMax_diff->Write();
   HxChi2NDOF ->Write();
   HxChi2NDOF_reduced ->Write();
   HxChi2NDOF_diff ->Write();
   HxChi2NDOF_reduced_diff ->Write();
   cp->Write();
   dp->Write();
   xp->Write();
   ch->Write();
   gSystem->Exec(Form("mv %s %s.tmp",PNGNameCoeffPulls,PNGNameCoeffPulls));
   cp->Print(PNGNameCoeffPulls);
//   
   gSystem->Exec(Form("mv %s %s.tmp",PNGNameModErrPulls,PNGNameModErrPulls));
   dp->Print(PNGNameModErrPulls);
//   
   gSystem->Exec(Form("mv %s %s.tmp",PNGNameCoeffNotGen,PNGNameCoeffNotGen));
   xp->Print(PNGNameCoeffNotGen);
//   
   gSystem->Exec(Form("mv %s %s.tmp",PNGNameChi2PVal,PNGNameChi2PVal));
   ch->Print(PNGNameChi2PVal);
// gStyle->SetOptStat(000000);
// gStyle->SetOptFit(000000);
  }else{
   RooWorkspace wsb("wsb","workspace sideband");
// set to the fit range!!!
   RooRealVar *max_sbl=new RooRealVar(Form("max_sbl_bin%d_%s",Q2Bin,RunEra),Form("max_sbl_bin%d_%s",Q2Bin,RunEra),XMaxSBL);
   RooRealVar *min_sbr=new RooRealVar(Form("min_sbr_bin%d_%s",Q2Bin,RunEra),Form("min_sbr_bin%d_%s",Q2Bin,RunEra),XMinSBR);

   wsb.import(*BernSideBand);
   wsb.import(*bkg_mass_sb);
   wsb.import(*max_sbl);
   wsb.import(*min_sbr);
   wsb.writeToFile(OutSaveFileName);
   cout<<"save  workspace ==> wsb in "<<OutSaveFileName<<"\n"<<endl;
  }
//   tagged_mass->setRange(XMinSign,XMaxSign);


  std::vector<double> DataBinContent;
  std::vector<double> PdfBinContent;

//  int icount=0;
  double pdfVal = 0.;
  double totalPdf = 0.;
  double xBinwPlot = HSBFunc->GetXaxis()->GetBinWidth(1) ;
  double yBinwPlot = HSBFunc->GetYaxis()->GetBinWidth(1) ;
  double zBinwPlot = HSBFunc->GetZaxis()->GetBinWidth(1) ;
  for(int i = 0; i < xCosLHBin*NFact ; ++i) {
   double xi = XMinCosThetaL + xBinwPlot/2.+ i*xBinwPlot;
   for(int j = 0 ; j < xCosKHBin*NFact ; ++j) {
    double yj = XMinCosThetaK + yBinwPlot/2.+ j*yBinwPlot;
    for(int k = 0 ; k < xPhiHBin*NFact  ; ++k) {
    double zk = XMinPhi + zBinwPlot/2.+ k*zBinwPlot;
//    if (xi<=XMaxCosThetaL && yj<=XMaxCosThetaK && zk<=XMaxPhi){
//      x.setVal(xi);
//      y.setVal(yj);
//      z.setVal(zk);
     ctL->setVal(xi);
     ctK->setVal(yj);
     phi->setVal(zk);
     if(MCStudy) {
      pdfVal = Bern->getVal();
     }else{
      pdfVal = BernSideBand->getVal();
     } 
//     double pdfVal = BernSideBand->getVal();
//     totalPdf += pdfVal*xBinwPlot*yBinwPlot*zBinwPlot;
//     double pdfVal = BernSideBand->evaluateInt(xBinwPlot,yBinwPlot,zBinwPlot);
//     cout<<Form("pdfVal(%d) = %f",icount,pdfVal)<<endl;
//     icount++;
     dataSBMass_Plot->add(RooArgSet(*ctL,*ctK,*phi));
     DataBinContent.push_back(HxReco->GetBinContent(i,j,k));
     PdfBinContent.push_back(pdfVal);
     totalPdf += pdfVal;
    }
   }
  }
  double dataPlotEntries=dataSBMass_Plot->sumEntries();
//  icount=0;
  double Chi2 =0.;
  int iskip =0;
  int NDegreeofFreedomN =0;
//  int NDegreeofFreedom =-(NumParamFree+1);
  double SBEntries= HxReco->GetEntries();
//  cout <<"totalPdf = "<<totalPdf*Vol3D<<endl;
  cout <<"dataPlotEntries = "<<dataPlotEntries<<endl;
  cout <<"totalPdf (*vol bin) = "<<totalPdf*xBinwPlot*yBinwPlot*zBinwPlot<<endl;
  cout <<"totalPdf            = "<<totalPdf*xBinwPlot*yBinwPlot*zBinwPlot/(8*TMath::Pi())<<endl;
  cout <<"SBEntries= "<<SBEntries<<endl;
  cout <<"SBEntries= "<<dataSBAng->sumEntries()<<"  (Check)"<<endl;
  cout <<"nentgrid = "<<dataPlotEntries<<endl;
  double LLChi2 =0.;
  double ProbFuncTot =0;
  const RooArgSet* pdfObs = 0 ;
  if(MCStudy){
    pdfObs = Bern->getObservables(dataSBMass_Plot) ;
  }else{
    pdfObs = BernSideBand->getObservables(dataSBMass_Plot) ;
  }  
//  const RooArgSet* pdfObs = BernSideBand->getObservables(dataSBMass_Plot) ;
   for (int i = 0; i < dataPlotEntries; ++i) {
     pdfObs = dataSBMass_Plot->get(i);
     *ctL = pdfObs->getRealValue("ctL" );
     *ctK = pdfObs->getRealValue("ctK" );
     *phi = pdfObs->getRealValue("phi" );
     double pdfVal=PdfBinContent[i];
//     double pdfVal = BernSideBand->evaluateInt(xBinwPlot,yBinwPlot,zBinwPlot);
//     cout<<Form("pdfVal(%d) = %f dopo",icount,pdfVal)<<endl;
//     icount++;
     if(pdfVal<0.0) std::cout<<"Warning!!!: SB Model "<<pdfVal <<"0 in CosL="<<ctL->getVal()<<" CosK="<<ctK->getVal()<<" Phi="<<phi->getVal()<<std::endl;
     HSBFunc->Fill(ctL->getVal(),ctK->getVal(),phi->getVal(), SBEntries*pdfVal/totalPdf);
     if(pdfVal>0.) {
      double ProbFunc = SBEntries*pdfVal/totalPdf;
//      double Chi2Temp = ( DataBinContent[i] - ProbFunc)*(DataBinContent[i] - ProbFunc)/ProbFunc;
      Chi2 = Chi2+( DataBinContent[i] - ProbFunc)*(DataBinContent[i] - ProbFunc)/ProbFunc;
      if(DataBinContent[i]!=0){
        LLChi2 = LLChi2+DataBinContent[i]*(log(DataBinContent[i])-log(ProbFunc));
      } else {
        iskip++;
      }	
//      if(DataBinContent[i]!=0) LLChi2 = LLChi2+DataBinContent[i]*(log(DataBinContent[i]/ProbFunc))+ProbFunc-DataBinContent[i];
      ProbFuncTot += ProbFunc;
//       if(Chi2Temp>1.) {
//        icount++;
//        std::cout<<"======== "<<icount<<" ======================================================================\n"<<std::endl;
//        std::cout<<"==> Chi2Temp = "<<Chi2Temp<<" CosL="<<xCosL_x.getValue()<<" CosK"<<xCosK_y.getValue()<<" Phi="<<xPhiK_z.getValue()<<"\n"<<std::endl;
//        std::cout<<"==> ProbFunc ="<<ProbFunc<<" DataCont="<<DataBinContent[i]<<"\n"<<std::endl;
//       } 
      NDegreeofFreedomN++; 
     }
     //cout <<xGene_w.getValue() <<endl;
  }
  LLChi2=2*LLChi2;
  int NDegreeofFreedomR = NDegreeofFreedomN-NumParamFree;
  cout <<"Integrated Probability Function "<<ProbFuncTot <<endl;
  cout <<"In LLchi2 skipped cell = "<<iskip<<endl;
  double  PValueModelMin = TMath::Prob(Chi2, NDegreeofFreedomR);
  double  PValueModelMax = TMath::Prob(Chi2, NDegreeofFreedomN);
  double  PValueLLModelMin  = TMath::Prob(LLChi2, NDegreeofFreedomR);
  double  PValueLLModelMax  = TMath::Prob(LLChi2, NDegreeofFreedomN);
  cout <<"------------------------------------------------------- "<<endl;
  cout <<"Even Binning ("<<xCosLHBin*xCosKHBin*xPhiHBin <<" Bins) "<<endl;
  cout <<"------------------------------------------------------- "<<endl;
  cout <<"Chi2 Sideband			= "<<Chi2<<endl;
  cout <<"Chi2/NDOF			= "<<Chi2/NDegreeofFreedomN<<endl;
  cout <<"Chi2/NDOF_Reduced		= "<<Chi2/NDegreeofFreedomR<<endl;
  cout <<"P-Value		Min	= "<<PValueModelMin<<endl;
  cout <<"P-Value		Max	= "<<PValueModelMax<<endl;
  cout <<"LLChi2 Sideband		= "<<LLChi2<<endl;
  cout <<"LLChi2/NDOF			= "<<LLChi2/(NDegreeofFreedomN)<<endl;
  cout <<"LLChi2/NDOF_Reduced		= "<<LLChi2/(NDegreeofFreedomR)<<endl;
  cout <<"LL P-Value		Min	= "<<PValueLLModelMin<<endl;
  cout <<"LL P-Value		Max	= "<<PValueLLModelMax<<endl;
  cout <<"NDOF	(reduced)		= "<<NDegreeofFreedomR<<endl;
  cout <<"NDOF				= "<<NDegreeofFreedomN<<endl;
  cout <<"Num Free Param.		= "<<NumParamFree<<endl;
  cout <<"------------------------------------------------------- "<<endl;
  HSBFunc->Sumw2();

//======================================================================================
// Adaptive Binning GOF...
//======================================================================================  
TCanvas* ca = new TCanvas("ca","Adaptive Binning Histograms",200,200,800,800);
  ca->Divide(2,2);   
  
//
  int iCorreTagOrig = CorreAdaptX.size();
  xAdaptNumBinC = int(iCorreTagOrig/MinContAdaptBin);
  int iCorreTag=xAdaptNumBinC*MinContAdaptBin;
  std::cout<<"TKDTreeBinning Start   "<<std::endl;
  std::cout<<"TKDTreeBinning set iCorreTag	  = "<<iCorreTag<<std::endl;
  std::cout<<"TKDTreeBinning xAdaptNumBinC    = "<<xAdaptNumBinC<<std::endl;
  double *RecoAdaptC = new double[NDim*iCorreTag];
  for (int iC=0;iC<iCorreTag;iC++) { 
   RecoAdaptC[iC]	     =  CorreAdaptX[iC];
   RecoAdaptC[iC+  iCorreTag]=  CorreAdaptY[iC];
   RecoAdaptC[iC+2*iCorreTag]=  CorreAdaptZ[iC];
  } 
  RecoAdaptBinsC = new TKDTreeBinning(iCorreTag, NDim, RecoAdaptC, xAdaptNumBinC);
  int nbinsC =RecoAdaptBinsC->GetNBins();
  std::cout<<"TKDTreeBinning nbinsC    = "<<nbinsC<<std::endl;
  TH2Poly* h2polxyContC = new TH2Poly("h2polxyContC", "adapt. binning contents [CosL CosK]", RecoAdaptBinsC->GetDataMin(0), RecoAdaptBinsC->GetDataMax(0), RecoAdaptBinsC->GetDataMin(1), RecoAdaptBinsC->GetDataMax(1));
  TH2Poly* h2polxzContC = new TH2Poly("h2polxzContC", "adapt. binning contents [CosL Phi ]", RecoAdaptBinsC->GetDataMin(0), RecoAdaptBinsC->GetDataMax(0), RecoAdaptBinsC->GetDataMin(2), RecoAdaptBinsC->GetDataMax(2));
  TH2Poly* h2polxyDensC = new TH2Poly("h2polxyDensC", "adapt. binning density [CosL CosK]", RecoAdaptBinsC->GetDataMin(0), RecoAdaptBinsC->GetDataMax(0), RecoAdaptBinsC->GetDataMin(1), RecoAdaptBinsC->GetDataMax(1));
  TH2Poly* h2polxzDensC = new TH2Poly("h2polxzDensC", "adapt. binning density [CosL Phi ]", RecoAdaptBinsC->GetDataMin(0), RecoAdaptBinsC->GetDataMax(0), RecoAdaptBinsC->GetDataMin(2), RecoAdaptBinsC->GetDataMax(2));
  const double* binsMinEdgesC = RecoAdaptBinsC->GetBinsMinEdges();
  const double* binsMaxEdgesC = RecoAdaptBinsC->GetBinsMaxEdges();
  int edgeDim=0;       
  std::vector<double> DataAdaptBinContent;
  std::vector<double> Vol3DAdaptBin;
  const double* xyzvar;
  const double* xyzbinw;
  RooDataSet *dataAdapt   = new RooDataSet("dataAdapt", "dataAdapt",  RooArgSet(*ctL,*ctK,*phi,BWidthX,BWidthY,BWidthZ));
  std::vector<int > AdaptExcludedEvents;
  for (int i = iCorreTag; i < iCorreTagOrig; ++i) {
  
   double point[3] = {CorreAdaptX[i],CorreAdaptY[i],CorreAdaptZ[i]};
   
   AdaptExcludedEvents.push_back(RecoAdaptBinsC->FindBin(point));
   cout<<Form("CosL=%f CosK=%f Phi=%f point = %d",CorreAdaptX[i],CorreAdaptY[i],CorreAdaptZ[i],RecoAdaptBinsC->FindBin(point))<<endl;
   
  }
  double vol =0;
  totalPdf = 0;
  vector<double> pdfVals_Model3DAdapt;
  for (int i = 0; i < nbinsC; ++i) {
     edgeDim = i * NDim;
     h2polxyContC->AddBin(binsMinEdgesC[edgeDim], binsMinEdgesC[edgeDim + 1], binsMaxEdgesC[edgeDim], binsMaxEdgesC[edgeDim + 1]);
     h2polxzContC->AddBin(binsMinEdgesC[edgeDim], binsMinEdgesC[edgeDim + 2], binsMaxEdgesC[edgeDim], binsMaxEdgesC[edgeDim + 2]);
     h2polxyDensC->AddBin(binsMinEdgesC[edgeDim], binsMinEdgesC[edgeDim + 1], binsMaxEdgesC[edgeDim], binsMaxEdgesC[edgeDim + 1]);
     h2polxzDensC->AddBin(binsMinEdgesC[edgeDim], binsMinEdgesC[edgeDim + 2], binsMaxEdgesC[edgeDim], binsMaxEdgesC[edgeDim + 2]);
     xyzvar  = RecoAdaptBinsC->GetBinCenter(i);
     xyzbinw = RecoAdaptBinsC->GetBinWidth(i);
     ctL->setVal(xyzvar[0]);
     ctK->setVal(xyzvar[1]);
     phi->setVal(xyzvar[2]);
     BWidthX.setVal(xyzbinw[0]);
     BWidthY.setVal(xyzbinw[1]);
     BWidthZ.setVal(xyzbinw[2]);
     vol += xyzbinw[0]*xyzbinw[1]*xyzbinw[2];
//      cout<<"================================"<<endl;
//      cout<<"NBIN    = "<<i<<endl;
//      cout<<"BWidthX = "<<xyzbinw[0]<<endl;
//      cout<<"BWidthY = "<<xyzbinw[1]<<endl;
//      cout<<"BWidthZ = "<<xyzbinw[2]<<endl;
//      cout<<"BWidthX = "<<binsMaxEdgesC[edgeDim]-binsMinEdgesC[edgeDim]<<endl;
//      cout<<"BWidthY = "<<binsMaxEdgesC[edgeDim + 1]-binsMinEdgesC[edgeDim + 1]<<endl;
//      cout<<"BWidthZ = "<<binsMaxEdgesC[edgeDim + 2]-binsMinEdgesC[edgeDim + 2]<<endl;
//   double pdfVal = BernSideBand->getVal(RooArgSet(*ctL,*ctK,*phi))*xyzbinw[0]*xyzbinw[1]*xyzbinw[2];
     double pdfVal = BernSideBand->evaluateInt(xyzbinw[0],xyzbinw[1],xyzbinw[2]);
     dataAdapt->add(RooArgSet(*ctL,*ctK,*phi,BWidthX,BWidthY,BWidthZ));
     DataAdaptBinContent.push_back(RecoAdaptBinsC->GetBinContent(i));
     Vol3DAdaptBin.push_back(xyzbinw[0]*xyzbinw[1]*xyzbinw[2]);
     pdfVals_Model3DAdapt.push_back(pdfVal);
     totalPdf += pdfVal;
  }
  cout <<"Vol tot  [Adapt] = "<<vol<<endl;
  
//  
//

  cout <<"totalPdf [Adapt] = "<<totalPdf<<endl;
  double AdaptProbFuncTot =0;
//double   icount = 0;
  double lambda = 2./3.;
  double AdaptChi2 = 0;
  double AdaptLLChi2 = 0;
  double AdaptPowerDivergence = 0. ;
  int NDegreeofFreedomAdaptN = 0;  
//  int NDegreeofFreedomAdapt =-(NumParamFree);  
  for (int ii = 0; ii < nbinsC; ++ii){
       h2polxyContC->SetBinContent(ii+1, RecoAdaptBinsC->GetBinContent(ii));
       h2polxzContC->SetBinContent(ii+1, RecoAdaptBinsC->GetBinContent(ii));
       h2polxyDensC->SetBinContent(ii+1, RecoAdaptBinsC->GetBinDensity(ii));
       h2polxzDensC->SetBinContent(ii+1, RecoAdaptBinsC->GetBinDensity(ii));
       double pdfVal = pdfVals_Model3DAdapt[ii];
       for (int iii = 0; iii < iCorreTagOrig-iCorreTag; ++iii) {
        if(AdaptExcludedEvents[iii]==ii){
	 DataAdaptBinContent[ii]++;
        }
       }	 	
       if( DataAdaptBinContent[ii]!=MinContAdaptBin) std::cout<<Form("DataAdaptBinContent[%d]=%f",ii,DataAdaptBinContent[ii])<<std::endl;
       if(pdfVal!=0.) {
//        double ProbFunc = pdfVal/totalPdf;
//        double ProbFunc = iCorreTag*pdfVal*Vol3DAdaptBin[ii]/totalPdf;
        double ProbFunc = SBEntries*pdfVal/totalPdf;
	AdaptPowerDivergence = AdaptPowerDivergence+ (DataAdaptBinContent[ii])*(pow(DataAdaptBinContent[ii]/ProbFunc,lambda)-1);

//      double AdaptChi2Temp = ( DataAdaptBinContent[ii] - ProbFunc)*(DataAdaptBinContent[ii] - ProbFunc)/(ProbFunc);
//      AdaptLLChi2 = AdaptLLChi2+DataAdaptBinContent[ii]*(log(DataAdaptBinContent[ii])-log(ProbFunc));
      AdaptLLChi2 = AdaptLLChi2+DataAdaptBinContent[ii]*(log(DataAdaptBinContent[ii])-log(ProbFunc))+ProbFunc-DataAdaptBinContent[ii];
//      AdaptChi2 = AdaptChi2+ AdaptChi2Temp;
      AdaptChi2 = AdaptChi2+ ( DataAdaptBinContent[ii] - ProbFunc)*(DataAdaptBinContent[ii] - ProbFunc)/(ProbFunc);
      AdaptProbFuncTot = ProbFunc+AdaptProbFuncTot;
//        if(AdaptChi2Temp>1.) {
// 	icount++;
// 	std::cout<<"======== "<<icount<<" ======================================================================\n"<<std::endl;
// 	std::cout<<"==> AdaptChi2Temp = "<<AdaptChi2Temp<<" CosL="<<xCosL_x.getValue()<<" CosK"<<xCosK_y.getValue()<<" Phi="<<xPhiK_z.getValue()<<"\n"<<std::endl;
// 	std::cout<<"==> ProbFunc ="<<ProbFunc<<" DataAdaptCont="<<DataAdaptBinContent[ii]<<"\n"<<std::endl;
//        } 
      NDegreeofFreedomAdaptN++;
     }
   }
  AdaptPowerDivergence = 2.*AdaptPowerDivergence/(lambda+1)/lambda ;	
  cout <<"Integrated Probability Function [Adapt]"<<ProbFuncTot <<endl;
  int NDegreeofFreedomAdaptR = NDegreeofFreedomAdaptN-NumParamFree;
  AdaptLLChi2=2*AdaptLLChi2;
  PValueModelMin = TMath::Prob(AdaptChi2, NDegreeofFreedomAdaptR);
  PValueModelMax = TMath::Prob(AdaptChi2, NDegreeofFreedomAdaptN);
  PValueLLModelMin =  TMath::Prob(AdaptLLChi2, NDegreeofFreedomAdaptR);
  PValueLLModelMax =  TMath::Prob(AdaptLLChi2, NDegreeofFreedomAdaptN);
  double PValuePDModelMin =  TMath::Prob(AdaptPowerDivergence, NDegreeofFreedomAdaptR);
  double PValuePDModelMax =  TMath::Prob(AdaptPowerDivergence, NDegreeofFreedomAdaptN);
  cout <<"------------------------------------------------------- "<<endl;
  cout <<"Adaptive Binning ("<<xAdaptNumBinC <<" Bins)\n"          <<endl;
  cout <<"MinContAdaptBin = "<<MinContAdaptBin                     <<endl;
  cout <<"------------------------------------------------------- "<<endl;
  cout <<"Chi2 Sideband			[Adapt] = "<<AdaptChi2<<endl;
  cout <<"Chi2/NDOF			[Adapt] = "<<AdaptChi2/NDegreeofFreedomAdaptR<<endl;
  cout <<"Chi2/NDOF_Reduced		[Adapt] = "<<AdaptChi2/NDegreeofFreedomAdaptN<<endl;
  cout <<"P-Value			Min	[Adapt] = "<<PValueModelMin<<endl;
  cout <<"P-Value			Max	[Adapt] = "<<PValueModelMax<<endl;
  cout <<"LL Chi2 Sideband		[Adapt] = "<<AdaptLLChi2<<endl;
  cout <<"LL Chi2/NDOF			[Adapt] = "<<AdaptLLChi2/NDegreeofFreedomAdaptN<<endl;
  cout <<"LL Chi2/NDOF_Reduced		[Adapt] = "<<AdaptLLChi2/NDegreeofFreedomAdaptR<<endl;
  cout <<"LL P-Value		Min	[Adapt] = "<<PValueLLModelMin<<endl;
  cout <<"LL P-Value		Max	[Adapt] = "<<PValueLLModelMax<<endl;
  cout <<"PowerDivergence			[Adapt] = "<<AdaptPowerDivergence<<endl;
  cout <<"PowerDivergence/NDOF		[Adapt] = "<<AdaptPowerDivergence/NDegreeofFreedomAdaptN<<endl;
  cout <<"PowerDivergence/NDOF_Reduced	[Adapt] = "<<AdaptPowerDivergence/NDegreeofFreedomAdaptR<<endl;
  cout <<"PD P-Value		Min	[Adapt] = "<<PValuePDModelMin<<endl;
  cout <<"PD P-Value		Max	[Adapt] = "<<PValuePDModelMax<<endl;
  cout <<"NDOF	(reduced)		[Adapt] = "<<NDegreeofFreedomAdaptR<<endl;
  cout <<"NDOF				[Adapt] = "<<NDegreeofFreedomAdaptN<<endl;
  cout <<"Num Free Param.	   	   = "<<NumParamFree<<endl;
  cout <<"------------------------------------------------------- "<<endl;

 //  h2polxyContC->Draw("lego");
  h2polxyContC->Draw("COLZ L");
  ca->Update();   
  ca->cd(2);
//  h2polxzContC->Draw("lego");
  h2polxzContC->Draw("COLZ L");
  ca->Update();   
  ca->cd(3);
  h2polxyDensC->Draw("COLZ L");
  ca->Update();   
  ca->cd(4);
  h2polxzDensC->Draw("COLZ L");
  ca->Update();   
 //    h2polxyC->Draw("LEGO");
//    ca->Update();   
//    ca->cd(4);
//    h2polxzC->Draw("LEGO");
//    ca->Update();   
//
// Adaptive Binning...End
//  



//  cout<<std::scientific << std::setprecision(40)<<"Norm PDF SideBand = "<<modelPlot->normalize()<<";"<<endl;
//==== closure test plots
//   modelPlot->setData(dataPlot);
//   vector<vector<double> > pdfVals_Model3D_Reco = modelPlot->getCompProbsAtDataPoints();
//   pdfVal = 0.0;
//   for (int i = 0; i < dataPlot->getNumEvents(); ++i) {
//     dataReco->loadEvent(i);
//     pdfVal = pdfVals_Model3D_Reco[0][i];
//     if (pdfVal<0.0) std::cout<<"Warning!!!: Effi Model Reco Test"<<pdfVal <<"0 in CosL="<<xCosL_x.getValue()<<" CosK="<<xCosK_y.getValue()<<" Phi="<<xPhiK_z.getValue()<<std::endl;
// //     double xL = xCosL_x.getValue();
// //     double yK = xCosK_y.getValue();
// //     double zP = xPhiK_z.getValue();
// //    HSBFunc->Fill(xL,yK,zP, pdfVal);
//      HSideBandRecoTest->Fill(xCosL_x.getValue(),xCosK_y.getValue(),xPhiK_z.getValue(), pdfVal*xGene_w.getValue());
// //     HSideBandCosLFunc->Fill(xL, pdfVal);
// //     HSideBandCosKFunc->Fill(yK, pdfVal);
// //     HSideBandPhiFunc ->Fill(zP, pdfVal);
//   }
//  HSideBandRecoTest->Sumw2();

  TH1D* HSBFuncX = (TH1D*) HSBFunc->ProjectionX("HSBFuncX",1,HSBFunc->GetNbinsY(),1,HSBFunc->GetNbinsZ());HSBFuncX->SetTitle(Form("Cos#theta_{L} Projection [Q^{2} bin %d run %s]",Q2Bin,RunEra));
  TH1D* HSBFuncY = (TH1D*) HSBFunc->ProjectionY("HSBFuncY",1,HSBFunc->GetNbinsX(),1,HSBFunc->GetNbinsZ());HSBFuncY->SetTitle(Form("Cos#theta_{K} Projection [Q^{2} bin %d run %s]",Q2Bin,RunEra));
  TH1D* HSBFuncZ = (TH1D*) HSBFunc->ProjectionZ("HSBFuncZ",1,HSBFunc->GetNbinsX(),1,HSBFunc->GetNbinsY());HSBFuncZ->SetTitle(Form("#phi Projection [Q^{2} bin %d run %s]",Q2Bin,RunEra));
  TH2D* HSBFuncXY = (TH2D*) HSBFunc->Project3D("xy");HSBFuncXY->SetTitle(Form("Profilo 2D Modello (Cos#theta_{l},Cos#theta_{k})    [q^{2} bin %d Run II %s]",Q2Bin,RunEra));
  TH2D* HSBFuncZY = (TH2D*) HSBFunc->Project3D("zy");HSBFuncZY->SetTitle(Form("Profilo 2D Modello (#varphi,Cos#theta_{k})    [q^{2} bin %d Run II %s]",Q2Bin,RunEra));
  TH2D* HSBFuncZX = (TH2D*) HSBFunc->Project3D("zx");HSBFuncZX->SetTitle(Form("Profilo 2D Modello (#varphi,Cos#theta_{l})   [q^{2} bin %d Run II %s]",Q2Bin,RunEra));
//  gStyle->SetOptStat(111111);
//  gStyle -> SetOptFit(111111);
  
   TH3D *HSideBand3D  = (TH3D*)HxReco->Clone(); HSideBand3D->SetName("HSideBand3D");HSideBand3D->Sumw2();
  
  
  
   TH1D* HSideBand3DX = (TH1D*) HSideBand3D->ProjectionX("HSideBand3DX",1,HSideBand3D->GetNbinsY(),1,HSideBand3D->GetNbinsZ());HSideBand3DX->SetTitle(Form("Cos#theta_{L} Projection [Q^{2} bin %d run %s]",Q2Bin,RunEra));
   TH1D* HSideBand3DY = (TH1D*) HSideBand3D->ProjectionY("HSideBand3DY",1,HSideBand3D->GetNbinsX(),1,HSideBand3D->GetNbinsZ());HSideBand3DY->SetTitle(Form("Cos#theta_{K} Projection [Q^{2} bin %d run %s]",Q2Bin,RunEra));
   TH1D* HSideBand3DZ = (TH1D*) HSideBand3D->ProjectionZ("HSideBand3DZ",1,HSideBand3D->GetNbinsX(),1,HSideBand3D->GetNbinsY());HSideBand3DZ->SetTitle(Form("#phi Projection [Q^{2} bin %d run %s]",Q2Bin,RunEra));
   TH2D* HSideBand3DXY = (TH2D*) HSideBand3D->Project3D("xy");HSideBand3DXY->SetTitle(Form("2D Projection (Cos#theta_{l},Cos#theta_{k})    [q^{2} bin %d Run II %s]",Q2Bin,RunEra));
   TH2D* HSideBand3DZY = (TH2D*) HSideBand3D->Project3D("zy");HSideBand3DZY->SetTitle(Form("2D Projection (#varphi,Cos#theta_{k})    [q^{2} bin %d Run II %s]",Q2Bin,RunEra));
   TH2D* HSideBand3DZX = (TH2D*) HSideBand3D->Project3D("zx");HSideBand3DZX->SetTitle(Form("2D Projection (#varphi,Cos#theta_{l})    [q^{2} bin %d Run II %s]",Q2Bin,RunEra));




//////////////////////////////////
//
// SideBand Plots
//
//////////////////////////////////
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFontSize(0.05) ;
  c6->cd(1);
//  HSideBand3DX->Draw("E1");
//  HSideBandX->Draw("E1");
  TH1D* HSBFuncX_ratio = (TH1D*)HSBFuncX->Clone();
  HSBFuncX_ratio->Rebin(NFact);
  HSBFuncX_ratio->SetLineWidth(2.);
  HSBFuncX_ratio->SetLineColor(kRed);
//  HSBFuncX->Scale(HSideBandX->Integral()/HSBFuncX->Integral()*NFact);
//  HSBFuncX->Scale(HSideBand3DX->Integral()/HSBFuncX->Integral());
//  HSBFuncX->Draw("same,HIST C");
  HSideBand3DX->SetMinimum(SetMinProj);
  RatioDataModel3DX = new TRatioPlot(HSideBand3DX,HSBFuncX_ratio);
  RatioDataModel3DX->SetGraphDrawOpt("L");
  RatioDataModel3DX->SetSeparationMargin(0.0);
  RatioDataModel3DX->SetH1DrawOpt("E1");
  RatioDataModel3DX->SetH2DrawOpt("HIST C");
  RatioDataModel3DX->Draw();
  RatioDataModel3DX->GetUpperPad()->cd();;
  HSBFuncX->Scale(NFact);
  HSBFuncX->Draw("same HIST C");
  RatioDataModel3DX->GetLowerRefGraph()->SetMinimum(SetMinRatio);
  RatioDataModel3DX->GetLowerRefGraph()->SetMaximum(SetMaxRatio);
  c6->Update();
  c6->cd(2);
//  HSideBandY->Draw("E1");
//  HSideBand3DY->Draw("E1");
  TH1D* HSBFuncY_ratio = (TH1D*)HSBFuncY->Clone();
  HSBFuncY_ratio->Rebin(NFact);
  HSBFuncY_ratio->SetLineWidth(2.);
  HSBFuncY_ratio->SetLineColor(kRed);
//  HSBFuncY->Scale(HSideBand3DY->Integral()/HSBFuncY->Integral());
//  HSBFuncY->Scale(HSideBandY->Integral()/HSBFuncY->Integral()*NFact);
//  HSBFuncY->Draw("same,HIST C");
  HSideBand3DY->SetMinimum(SetMinProj);
  RatioDataModel3DY = new TRatioPlot(HSideBand3DY,HSBFuncY_ratio);
  RatioDataModel3DY->SetGraphDrawOpt("L");
  RatioDataModel3DY->SetSeparationMargin(0.0);
  RatioDataModel3DY->SetH1DrawOpt("E1");
  RatioDataModel3DY->SetH2DrawOpt("HIST C");
  RatioDataModel3DY->Draw();
  RatioDataModel3DY->GetUpperPad()->cd();;
  HSBFuncY->Scale(NFact);
  HSBFuncY->Draw("same HIST C");
  RatioDataModel3DY->GetLowerRefGraph()->SetMinimum(SetMinRatio);
  RatioDataModel3DY->GetLowerRefGraph()->SetMaximum(SetMaxRatio);
  TLegend* leg_SBFunc = new TLegend(0.40,0.67,0.90,0.90);
  leg_SBFunc->SetTextSize(0.025) ;
  leg_SBFunc->SetTextAlign(13);
  leg_SBFunc->SetBorderSize(0.);
  leg_SBFunc->SetFillStyle(0);
//  leg_SBFunc->SetEntrySeparation(0.0001);
//  leg_SBFunc->SetNColumns(1);
//   leg_SBFunc->AddEntry(HSBFuncY ,Form( "%5.3f<#Chi^{2}_{/NDOF}<%5.3f", AdaptChi2/NDegreeofFreedomAdaptN,AdaptChi2/NDegreeofFreedomAdaptR),"");
//   leg_SBFunc->AddEntry(HSBFuncY ,Form( "reduced NDOF=%d [NDOF=%d] ", NDegreeofFreedomAdaptR,NDegreeofFreedomAdaptN),"");
//   leg_SBFunc->AddEntry(HSBFuncY ,Form( "%6.5f<P-Value Min<%6.5f ", PValueModelMin,PValueModelMax),"");
//   leg_SBFunc->AddEntry(HSBFuncY ,Form( "Events x  adapt. bin=%d ",MinContAdaptBin),"");
//   leg_SBFunc->AddEntry(HSBFuncY ,Form( "%5.3f<#Chi^{2}_{/NDOF}<%5.3f", AdaptPowerDivergence/NDegreeofFreedomAdaptN,AdaptPowerDivergence/NDegreeofFreedomAdaptR),"");
//   leg_SBFunc->AddEntry(HSBFuncY ,Form( "reduced NDOF=%d [NDOF=%d] ", NDegreeofFreedomAdaptR,NDegreeofFreedomAdaptN),"");
//   leg_SBFunc->AddEntry(HSBFuncY ,Form( "%6.5f<P-Value<%6.5f ", PValuePDModelMin,PValuePDModelMax),"");
//   leg_SBFunc->AddEntry(HSBFuncY ,Form( "num of Events x [adaptive] bin=%d ",MinContAdaptBin),"");
  leg_SBFunc->AddEntry(RatioDataModel3DY ,Form( "%5.3f<#Chi^{2}_{/NDOF}<%5.3f", AdaptPowerDivergence/NDegreeofFreedomAdaptN,AdaptPowerDivergence/NDegreeofFreedomAdaptR),"");
  leg_SBFunc->AddEntry(RatioDataModel3DY ,Form( "reduced NDOF=%d [NDOF=%d] ", NDegreeofFreedomAdaptR,NDegreeofFreedomAdaptN),"");
  leg_SBFunc->AddEntry(RatioDataModel3DY ,Form( "%6.5f<P-Value<%6.5f ", PValuePDModelMin,PValuePDModelMax),"");
  leg_SBFunc->AddEntry(RatioDataModel3DY ,Form( "num of Events x [adaptive] bin=%d ",MinContAdaptBin),"");
  leg_SBFunc->Draw();
  c6->Update();
  c6->cd(3);
//  HSideBandZ->Draw("E1");
//  HSideBand3DZ->Draw("E1");
  TH1D* HSBFuncZ_ratio = (TH1D*)HSBFuncZ->Clone();
  HSBFuncZ_ratio->Rebin(NFact);
  HSBFuncZ_ratio->SetLineWidth(2.);
  HSBFuncZ_ratio->SetLineColor(kRed);
//  HSBFuncZ->Scale(HSideBandZ->Integral()/HSBFuncZ->Integral()*NFact);
//  HSBFuncZ->Scale(HSideBand3DZ->Integral()/HSBFuncZ->Integral());
//  HSBFuncZ->Draw("same,HIST C");
  HSideBand3DZ->SetMinimum(SetMinProj);
  RatioDataModel3DZ = new TRatioPlot(HSideBand3DZ,HSBFuncZ_ratio);
  RatioDataModel3DZ->SetGraphDrawOpt("L");
  RatioDataModel3DZ->SetSeparationMargin(0.0);
  RatioDataModel3DZ->SetH1DrawOpt("E1");
  RatioDataModel3DZ->SetH2DrawOpt("HIST C");
  RatioDataModel3DZ->Draw();
  RatioDataModel3DZ->GetUpperPad()->cd();;
  HSBFuncZ->Scale(NFact);
  HSBFuncZ->Draw("same HIST C");
  RatioDataModel3DZ->GetLowerRefGraph()->SetMinimum(SetMinRatio);
  RatioDataModel3DZ->GetLowerRefGraph()->SetMaximum(SetMaxRatio);
  RatioDataModel3DZ->GetUpperRefXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  RatioDataModel3DZ->GetLowerRefGraph()->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  c6->Update();

  c6->cd(4); 
  TH1D* pdfHxMassQ2Set = (TH1D*)pdfHxMassQ2->Clone();
  HxMassQ2   ->GetXaxis()->SetRangeUser(HistMassL1,HistMassL2);
  HxMassQ2SB ->GetXaxis()->SetRangeUser(XMinSBL,XMaxSBR);
  pdfHxMassQ2->GetXaxis()->SetRangeUser(XMinSBL,XMaxSBR);
  sigHxMassQ2->GetXaxis()->SetRangeUser(XMinSBL,XMaxSBR);
  bkgHxMassQ2->GetXaxis()->SetRangeUser(XMinSBL,XMaxSBR);
  pdfHxMassQ2Set->GetXaxis()->SetRangeUser(XLeftSet,XRightSet);
//   HxMassQ2   ->GetXaxis()->SetLimits(HistMassL1,HistMassL2);
//   HxMassQ2SB ->GetXaxis()->SetLimits(HistMassL1,HistMassL2);
//   pdfHxMassQ2->GetXaxis()->SetLimits(HistMassL1,HistMassL2);
//   sigHxMassQ2->GetXaxis()->SetLimits(HistMassL1,HistMassL2);
//   bkgHxMassQ2->GetXaxis()->SetLimits(HistMassL1,HistMassL2);
//  pdfHxMassQ2Set->GetXaxis()->SetLimits(HistMassL1,HistMassL2);
//  sigHxMassQ2->SetLineStyle(kDashed);
//  bkgHxMassQ2->SetLineStyle(kDashed);
  pdfHxMassQ2Set->SetFillColor(kBlue);
  pdfHxMassQ2Set->SetFillStyle(3013);
  pdfHxMassQ2Set->SetLineWidth(1.0);
//  
  TLegend* leg_signSB = new TLegend(0.25,0.75,0.90,0.90);
  leg_signSB->SetTextSize(0.025) ;
  leg_signSB->SetTextAlign(11);
  leg_signSB->SetBorderSize(0.);
  leg_signSB->SetFillStyle(0);
//  leg_signSB->AddEntry(HxMassQ2 ,Form( "in red:"),"");
  leg_signSB->AddEntry(HxMassQ2 ,Form( "#color[2]{sideband for q^{2} bin %d [%2.1f<q^{2}<%2.1f]}", Q2Bin, Q2Min,Q2Max),"");
  if(SigmaProbSign==0){
   leg_signSB->AddEntry(HxMassQ2 ,Form( "sb entries=%4.0f [-%s#sigma,-%s#sigma]&[%s#sigma,%s#sigma]",\
    SBEntries,FMTNSigma1L,FMTNSigma2L,FMTNSigma1R,FMTNSigma2R),"");
   leg_signSB->AddEntry(HxMassQ2 ,Form( "bckg entries=%4.0f [-2#sigma,2#sigma]",NBckgInt2Sigma),"");
  }else{
   leg_signSB->AddEntry(HxMassQ2 ,Form( "sb entries=%4.0f",SBEntries),"");
   leg_signSB->AddEntry(HxMassQ2 ,Form( "bckg entries=%4.0f [95.5%% of signal]",NBckgInt2Sigma),"");
  }  
//  leg_signSB->AddEntry(HxMassQ2 ,Form( "bckg entries=%4.0f [-2#sigma,2#sigma]",NBckgInt2Sigma),"");
//  leg_signSB->AddEntry(HxMassQ2 ,Form( "#chi^{2}_{/NDOF}=%5.2f NDOF=%d", AdaptChi2/NDegreeofFreedomAdapt,NDegreeofFreedomAdapt),"pel");
// //  leg_signSB->SetHeader(Form( "B_{0} #rightarrow K^{*0}#mu^{+}#mu^{-} Mass Spectrum for bin %d [%2.1f<q^{2}<%2.1f] ", Q2Bin, Q2Min,Q2Max));
//   leg_signSB->AddEntry(HxMassQ2 ,Form( "B_{0} #rightarrow K^{*0}#mu^{+}#mu^{-} Mass Spectrum for bin %d [%2.1f<q^{2}<%2.1f] ", Q2Bin, Q2Min,Q2Max),"pel");
//   leg_signSB->AddEntry(HxMassQ2SB ,Form( "B_{0} #rightarrow K^{*0}#mu^{+}#mu^{-} Mass Spectrum for bin %d [%2.1f<q^{2}<%2.1f] ", Q2Bin, Q2Min,Q2Max),"pel");
//  gStyle->SetTitleFontSize(0.09) ;
  HxMassQ2->SetMaximum(1.5 * HxMassQ2->GetMaximum());
  HxMassQ2->SetTitle(Form("B_{0} #rightarrow K^{*0}#mu^{+}#mu^{-} Mass Spectrum [q^{2} bin %d run %s]",Q2Bin,RunEra));
//  HxMassQ2->SetTitle(Form( "B_{0} #rightarrow K^{*0}#mu^{+}#mu^{-} Mass Spectrum and sideband for q^{2} bin %d [%2.1f<q^{2}<%2.1f] ", Q2Bin, Q2Min,Q2Max));
  HxMassQ2->SetLineColor(kBlue);
  HxMassQ2->DrawCopy("E1,9");
  HxMassQ2->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
//  HxMassQ2SB->SetTitleSize(20);
  HxMassQ2SB->SetLineColor(kOrange);
  HxMassQ2SB->SetTitle(Form( "B_{0} #rightarrow K^{*0}#mu^{+}#mu^{-} Mass Spectrum  and sideband for q^{2} bin %d [%2.1f<q^{2}<%2.1f] ", Q2Bin, Q2Min,Q2Max));
  HxMassQ2SB->SetFillColor(kRed);
  HxMassQ2SB->DrawCopy("SAME,B,9");
  pdfHxMassQ2->DrawCopy("same,HIST C,9");
  sigHxMassQ2->DrawCopy("same,HIST C,9");
  bkgHxMassQ2->DrawCopy("same,HIST C,9");
  pdfHxMassQ2Set->DrawCopy("same,HIST C,9");
  leg_signSB->Draw();
  c6->Update();
  
//////////////////////////////////
//////////////////////////////////

  gStyle->SetPalette(57);
  cxy->cd();
  gPad->SetTheta(40.);
  gPad->SetPhi(40.);
  //HSideBand3DXY->Rebin2D(2,2);
  HSideBand3DXY->SetFillColor(38);
  HSideBand3DXY->Smooth();
//  HSideBand3DXY->Draw("SURF2 0");
  HSideBand3DXY->GetXaxis()->SetLabelSize(0);
  HSideBand3DXY->GetYaxis()->SetLabelSize(0);
  HSideBand3DXY->GetZaxis()->SetLabelSize(0);
  HSideBand3DXY->GetXaxis()->SetTitle("Cos#theta_{l}");
  HSideBand3DXY->GetYaxis()->SetTitle("Cos#theta_{k}");
  HSideBand3DXY->Draw("LEGO2 0 fbbb");
//  auto cutg = new TCutG("cutg",5);
//  cxy->cd(2);
//  HSBFuncXY->Scale(HSideBand3DXY->Integral()/HSBFuncXY->Integral()*NFact);
  HSBFuncXY->SetLineColor(kRed);
  HSBFuncXY->GetXaxis()->SetTitle("Cos#theta_{l}");
  HSBFuncXY->GetYaxis()->SetTitle("Cos#theta_{k}");
  HSBFuncXY->Draw("SURF A same fbbb");
  
//
  cyz->cd();
  gPad->SetTheta(40.);
  gPad->SetPhi(40.);
  //HSideBand3DZY->Rebin2D(2,2);
  HSideBand3DZY->Smooth();
  HSideBand3DZY->GetXaxis()->SetLabelSize(0);
  HSideBand3DZY->GetYaxis()->SetLabelSize(0);
  HSideBand3DZY->GetZaxis()->SetLabelSize(0);
  HSideBand3DZY->GetXaxis()->SetTitle("#varphi");
  HSideBand3DZY->GetYaxis()->SetTitle("Cos#theta_{k}");
  HSideBand3DZY->Draw("LEGO2 0 fbbb");
//  cyz->cd(2);
//  HSBFuncZY->Scale(HSideBand3DZY->Integral()/HSBFuncZY->Integral()*NFact);
  HSBFuncZY->SetLineColor(kRed);
  HSBFuncZY->GetXaxis()->SetTitle("#varphi");
  HSBFuncZY->GetYaxis()->SetTitle("Cos#theta_{k}");
  HSBFuncZY->Draw("SURF A same fbbb");
//
  cxz->cd();
  gPad->SetTheta(40.);
  gPad->SetPhi(40.);
  //HSideBand3DZX->Rebin2D(2,2);
  HSideBand3DZX->Smooth();
  HSideBand3DZX->GetXaxis()->SetLabelSize(0);
  HSideBand3DZX->GetYaxis()->SetLabelSize(0);
  HSideBand3DZX->GetZaxis()->SetLabelSize(0);
  HSideBand3DZX->GetXaxis()->SetTitle("#varphi");
  HSideBand3DZX->GetYaxis()->SetTitle("Cos#theta_{l}");
  HSideBand3DZX->Draw("LEGO2 0 fbbb");
//  cxz->cd(2);
//  HSBFuncZX->Scale(HSideBand3DZX->Integral()/HSBFuncZX->Integral()*NFact);
  HSBFuncZX->SetLineColor(kRed);
  HSBFuncZX->GetXaxis()->SetTitle("#varphi");
  HSBFuncZX->GetYaxis()->SetTitle("Cos#theta_{l}");
  HSBFuncZX->Draw("SURF A same fbbb");
 
  
  
//
// Mass Spectrum
//

  OutFile->cd();

  c1->cd();
//     TLegend* leg_sign = new TLegend(0.30,0.70,0.90,0.90);
//     leg_sign->SetTextSize(0.025) ;
//     leg_sign->SetTextAlign(31);
//     leg_sign->SetBorderSize(0.);
//     leg_sign->SetFillStyle(0);
//     leg_sign->SetHeader(Form( "B_{0} #rightarrow K^{*0}#mu^{+}#mu^{-} Mass Spectrum for bin %d [%2.1f<Q^{2}<%2.1f] ", Q2Bin, Q2Min,Q2Max));
//     leg_sign->AddEntry(HxMass ,"","");
//     if(signalYield->getError()!=0){
//       leg_sign->AddEntry(&HxMass ,Form( "Yield_{Sign} =    %5.0f  #pm %5.0f",signalYield->getValue(),signalYield->getError()),"");
//     }else{
//       leg_sign->AddEntry(&HxMass ,Form( "Yield_{Sign} =    %5.0f Fixed",signalYield->getValue()),"");
//     }
//     if(bckgYield->getError()!=0){
//       leg_sign->AddEntry(&HxMass ,Form( "Yield_{Bckg} =    %5.0f  #pm  %5.0f",bckgYield->getValue(),bckgYield->getError()),"");
//     }else{
//       leg_sign->AddEntry(&HxMass ,Form( "Yield_{Bckg} =    %5.0f  Fixed",bckgYield->getValue()),"");
//     }
//     
//     if(mean.getError()!=0){
//      leg_sign->AddEntry(&HxMass ,Form( "M_{B^{0}} =   %5.5f  #pm %5.5f",mean.getValue(),mean.getError()),"");
//     }else{
//      leg_sign->AddEntry(&HxMass ,Form( "M_{B^{0}} =   %5.5f Fixed",mean.getValue()),"");
//      }
//     if(sigma1.getError()!=0){
//      leg_sign->AddEntry(&HxMass ,Form( "#sigma#scale[0.6]{1}_{B^{0}} =   %5.5f  #pm %5.5f",sigma1.getValue(),sigma1.getError()),"");
//     }else{
//      leg_sign->AddEntry(&HxMass ,Form( "#sigma#scale[0.6]{1}_{B^{0}} =   %5.5f Fixed",sigma1.getValue()),"");
//     }
//     if(sigma2.getError()!=0){
//      leg_sign->AddEntry(&HxMass ,Form( "#sigma#scale[0.6]{2}_{B^{0}} =   %5.5f  #pm %5.5f",sigma2.getValue(),sigma2.getError()),"");
//     }else{
//      leg_sign->AddEntry(&HxMass ,Form( "#sigma#scale[0.6]{2}_{B^{0}} =   %5.5f Fixed",sigma2.getValue()),"");
//    }
  gStyle->SetTitleBorderSize(0);
//  gStyle->SetTitleFontSize(0.08) ;
  HxMass->SetFillStyle(0);
  HxMass->SetTitle("B_{0} #rightarrow K^{*0}#mu^{+}#mu^{-} Mass Spectrum");
  HxMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  HxMass->SetMarkerStyle(8);
  HxMass->SetMarkerSize(MarkerSizeSet);
  HxMass->SetStats(kFALSE);
  HxMass->Draw("E1");
//   HxMassQ2->SetLineColor(kBlue);
//   HxMassQ2->Draw("same"); 

  HxMass->Write();
  HxMassQ2SB->Write();
  HxMassQ2->Write();
  HxMassVsCosL->Write();
  HxMassVsCosK->Write();
  HxMassVsPhi->Write();
  HxReco->Write();
  HSideBand3D->Write();
  HSBFunc->Write();
  HSBFuncX->Write();
  HSBFuncY->Write();
  HSBFuncZ->Write();
  HSBFuncXY->Write();
  HSBFuncZY->Write();
  HSBFuncZX->Write();
  HSideBand3DX->Write();
  HSideBand3DY->Write();
  HSideBand3DZ->Write();
  HSideBand3DXY->Write();
  HSideBand3DZY->Write();
  HSideBand3DZX->Write();
  pdfHxMassQ2->Write();
  sigHxMassQ2->Write();
  bkgHxMassQ2->Write();
  
//   TH3D *H3Div = (TH3D *)HxReco->Clone(); 
//   H3Div->SetName("H3RecoDivGen");
//   H3Div->Divide(HxGene);

  
  
   c1->Write();
   c2->Write();
   gSystem->Exec(Form("mv %s %s.tmp",PNGNameMassQ2Hist,PNGNameMassQ2Hist));
   c2->Print(PNGNameMassQ2Hist);
   c6->Write();
   gSystem->Exec(Form("mv %s %s.tmp",PNGNameFitSB3D,PNGNameFitSB3D));
   c6->Print(PNGNameFitSB3D);
   h2polxyContC->Write();
   h2polxzContC->Write();
   h2polxyDensC->Write();
   h2polxzDensC->Write();
   ca->Write();
   
   TCanvas* cSideRooPlots = new TCanvas("cSideRooPlots","Sideband Projections",200,10,900,780);
   cSideRooPlots->Divide(2,2);
   RooPlot* xframe=ctL->frame(RooFit::Bins(xCosLHBin));
   dataSBAng->plotOn(xframe, RooFit::MarkerStyle(kPlus));
   BernSideBand->plotOn( xframe,RooFit::LineColor(kRed),Normalization(SBEntries,RooAbsReal::NumEvent),RooFit::ProjWData(RooArgSet(*ctL),*dataSBAng));
   cSideRooPlots->cd(1);
   xframe->Draw();
   RooPlot* yframe=ctK->frame(RooFit::Bins(xCosKHBin));
   dataSBAng->plotOn(yframe, RooFit::MarkerStyle(kPlus));
   BernSideBand->plotOn( yframe,RooFit::LineColor(kRed),Normalization(SBEntries,RooAbsReal::NumEvent),RooFit::ProjWData(RooArgSet(*ctK),*dataSBAng));
   cSideRooPlots->cd(2);
   yframe->Draw();
   RooPlot* zframe=phi->frame(RooFit::Bins(xPhiHBin));
   dataSBAng->plotOn(zframe, RooFit::MarkerStyle(kPlus));
   BernSideBand->plotOn( zframe,RooFit::LineColor(kRed),Normalization(SBEntries,RooAbsReal::NumEvent),RooFit::ProjWData(RooArgSet(*phi),*dataSBAng));
   cSideRooPlots->cd(3);
   zframe->Draw();
   gSystem->Exec(Form("mv %s %s.tmp",PNGNamePlotSB3D,PNGNamePlotSB3D));
   cSideRooPlots->Print(PNGNamePlotSB3D);
   cSideRooPlots->Write();


//   c7->Write();
//   c8->Write();
//   gSystem->Exec(Form("mv %s %s.tmp",PDFNameMass,PDFNameMass));
//   c1->Print(PDFNameMass);
//   sprintf(testo,"mv %s %s.tmp",PDFNameFitEffi3D,PDFNameFitEffi3D);
//   gSystem->Exec(Form("mv %s %s.tmp",PDFNameFitEffi3D,PDFNameFitEffi3D));
//   c6->Print(PDFNameFitEffi3D);
//   gSystem->Exec(Form("mv %s %s.tmp",PNGNameFitEffi3D,PNGNameFitEffi3D));
//   c6->Print(PNGNameFitEffi3D);
//   gSystem->Exec(Form("mv %s %s.tmp",PDFNameFitClosure,PDFNameFitClosure));
//   c8->Print(PDFNameFitClosure);
//   gSystem->Exec(Form("mv %s %s.tmp",PNGNameFitClosure,PNGNameFitClosure));
//   c8->Print(PNGNameFitClosure);
  OutFile->Close();
  std::cout<<"==========================================" <<std::endl;
  std::cout<<"==========================================" <<std::endl;
//=================================================================================  
//=================================================================================  
// START write par+norm txt files  
//=================================================================================  
//=================================================================================  
//   coeffy =0.0;
//   std::fstream *parListNormOutput =  new std::fstream(ListParNorm,ios::out);
//   std::fstream *parPlotNormOutput =  new std::fstream(ListPloNorm,ios::out);
//   if(parListNormOutput->is_open() && parPlotNormOutput->is_open() ){
//    std::cout<<"Open: "<<ListParNorm<<std::endl ;
//    std::cout<<"Open: "<<ListPloNorm<<std::endl ;
//    for(int i=0;i<numParameters;++i) {
//     if(fabs(coeffPoly[i].getValue())>fabs(coeffPoly[i].getError())){
//      coeffy=  coeffPoly[i].getValue();
//     }else{
//      coeffy=  0.0;
//     } 
//     *parListNormOutput <<std::scientific << std::setprecision(20)<< coeffy<<std::endl;
//     *parPlotNormOutput <<std::scientific << std::setprecision(20)<< coeffPoly[i].getValue()<<std::endl;
//    }
//
//    xCosL_x.setNumBins(SETNumBinsX);
//    xCosK_y.setNumBins(SETNumBinsY);
//    xPhiK_z.setNumBins(SETNumBinsZ);
//    cout<<"CosL integration NBins ="<<xCosL_x<<endl;
//    cout<<"CosK integration NBins ="<<xCosK_y<<endl;
//    cout<<"Phi  integration NBins ="<<xPhiK_z<<endl;
// //  
// // cout<<"cudaFreeHost(Norm)  = "<<(cudaFreeHost(model))<<endl;
//    cout<<"cudaFreeHost(modelPlot)  = "<<(cudaFreeHost(modelPlot))<<endl;
// //   cout<<"cudaFreeHost(model)      = "<<(cudaFreeHost(model))<<endl;
//    totalParams=0;
//    parListNormOutput->close();
//    parPlotNormOutput->close();
//    std::cout<<"Close: "<<ListParNorm<<std::endl ;
//    std::cout<<"Close: "<<ListPloNorm<<std::endl ;
//   }else{
//    if(!parListNormOutput->is_open()) std::cout<<"Error: can not open "<<ListParNorm<<std::endl ;
//    if(!parPlotNormOutput->is_open()) std::cout<<"Error: can not open "<<ListPloNorm<<std::endl ;
//    exit(1);
//   }
//=================================================================================  
//=================================================================================  
// END write par+norm txt files  
//=================================================================================  
//=================================================================================  
  std::cout<<"===================================================================="<<endl;
  std::cout<<"======== REMIND:  	NFact="<<NFact<<"		=============="<<std::endl;
  std::cout<<"====================================================================\n\n"<<endl;
  if(xCoeffIndex>0){
   std::cout<<"===================================================================="<<endl;
   std::cout<<Form("WARNING: Normalization could be fixed better if p(%d)=1",xCoeffIndex)<<std::endl ;
   std::cout<<"===================================================================="<<endl;
  }

  stopCPU = times(&stopProc);
  gettimeofday(&stopTime, NULL);
  // Print total minimization time
  double myCPU = stopCPU - startCPU;
  double totalCPU = myCPU; 

  timersub(&stopTime, &startTime, &totalTime);
  std::cout << "Wallclock time  : " << totalTime.tv_sec + totalTime.tv_usec/1000000.0 << " seconds." << std::endl;
  std::cout << "CPU time: " << (myCPU / CLOCKS_PER_SEC) << std::endl; 
  std::cout << "Total CPU time: " << (totalCPU / CLOCKS_PER_SEC) << std::endl; 
  myCPU = stopProc.tms_utime - startProc.tms_utime;
  std::cout << "Processor time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;
  std::cout<<"==========================================" <<std::endl;
  std::cout<<"==========================================" <<std::endl;
}

//============================================================================================================================================================
//============================================================================================================================================================
//============================================================================================================================================================
//
//  CreateInputHistoFile();
//
//============================================================================================================================================================
//============================================================================================================================================================
//============================================================================================================================================================
void CreateInputHistoFile(){   

  if(Folded) printf("****************************** WARNING: Folded ******************************\n");
  
  TFile*OutFileNtupla = TFile::Open(OutFileNameInputHisto,"RECREATE");
  RecoB0TreeOut = new TTree(OutputRecoB0TreeName,OutputRecoB0TreeName) ;
  RecoB0TreeOut -> SetAutoSave(5000000000);
  TCanvas* c2 = new TCanvas("c2","Fit Mass Spectrum",200,10,900,780);
  TCanvas* c3 = new TCanvas("c3","Reco Histograms",200,10,900,780);
  c3->Divide(2,2);  

  int nfile=0;
  TChain* RecoB0Tree = new TChain();  
  nfile = RecoB0Tree->Add(Form("%s/%s/%s",RecoDir,InputFileNameRecoB0,InputRecoB0TreeName));
  if( nfile==0 ||  !RecoB0Tree->GetFile() ){
    cout<<"Error:  no Reco files found!!!\n"<<endl;
    exit(1);
  }else{
    printf("Try to open %s/%s/%s\n",RecoDir,InputFileNameRecoB0,InputRecoB0TreeName);
    cout<<"Opening "<<nfile <<" Reco files found!!!\n"<<endl;
  }  
  if(!RecoB0Tree ){
    cout<<"TTree Reco Data: "<< InputRecoB0TreeName <<" not found!!!\n"<<endl;
    exit(1);
  }else{
    cout<<"TTree Reco Data: "<< InputRecoB0TreeName <<" OK FOUND!!!\n"<<endl;
  }  
  


  printf("(Mass Window     : xB0Mass>%8f && xB0Mass<%8f \n",XMinSign,XMaxSign);
//  printf("(SB Mass Windows : xB0Mass>%8f && xB0Mass<%8f || xB0Mass>%8f && xB0Mass<%8f)\n",XMinSBL,XMaxSBL,XMinSBR,XMaxSBR);




//   TH3D* HxReco = new   TH3D( "HxReco"    , "B^{0} Reco correct tagged",  xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
//  									 xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
// 									 xPhiHBin , XMinPhi, XMaxPhi );




	
						
 

//======================================================================
//======================================================================
//======================================================================
//
//			      RECONSTRUCTED EVENTS
//
//======================================================================
//======================================================================
//======================================================================
  double  tagged_mass	 ;
  double  cos_theta_l	 ;
  double  cos_theta_k	 ;
  double  phi_kst_mumu   ;
  double  mumuMass	 ;
  double  mumuMassE	 ;
  double  mumuPt   	  ;
  double  mumuPhi   	  ;
  double  mumuEta   	  ;
  double  recQ2 	 ;
  double  tagB0          ;
  double  mmk1		 ;
  double  mmk2		 ;
  double  bMass		 ;
  double  bBarMass	 ;
  double  dR_mum_trkm;
  double  dR_mup_trkp;
  bool    passB0Psi_lmnr ;
  bool    passB0Psi_jpsi ;
  bool    passB0Psi_psip ;
  bool    passB0Psi      ;
  double piMass = 0.13957039;
  double kMass = 0.493677;
  double BpMass = 5.2791;
  double B0Mass = 5.27962;
  double KstarMass = 0.892;
  double MassPsi2S=3.696;
  double MassJPsi =3.0969;
  double mmpip;
  double mmpim;
  double mmpipkm;
  double mmpimkp;
  double mmpiKaon;
  double mmpiKaon2;
  double mmpi1;
  double mmpi2;
  double mmpipi;
  double kstTrkmPt;
  double kstTrkmPhi;
  double kstTrkmEta;
  double kstTrkpPt;
  double kstTrkpPhi;
  double kstTrkpEta;
  double pi1Pt=0;
  double pi2Pt=0;
  double pimPt=0;
  double pipPt=0;
  double piKaon =0;
  double pipk =0;
  double pimk =0;
  double mmkp =0;
  double mmkm =0;
  double mmka1 =0;
  double mmka2 =0;
  double kstKp =0;
  double kstKm =0;
  double kst1  =0;
  double kst2  =0;
//
  RecoB0Tree->SetBranchAddress("tagB0"         ,&tagB0);
  RecoB0Tree->SetBranchAddress("tagged_mass"   ,&tagged_mass);
  RecoB0Tree->SetBranchAddress("cos_theta_l"   ,&cos_theta_l);
  RecoB0Tree->SetBranchAddress("cos_theta_k"   ,&cos_theta_k);
  RecoB0Tree->SetBranchAddress("phi_kst_mumu"  ,&phi_kst_mumu);
  RecoB0Tree->SetBranchAddress("mumuMass"      ,&mumuMass);
  RecoB0Tree->SetBranchAddress("mumuMassE"     ,&mumuMassE);
  RecoB0Tree->SetBranchAddress("mumuPt"        ,&mumuPt );
  RecoB0Tree->SetBranchAddress("mumuPhi"       ,&mumuPhi );
  RecoB0Tree->SetBranchAddress("mumuEta"       ,&mumuEta );
  RecoB0Tree->SetBranchAddress("kstTrkmPt"     ,&kstTrkmPt);
  RecoB0Tree->SetBranchAddress("kstTrkmPhi"    ,&kstTrkmPhi);
  RecoB0Tree->SetBranchAddress("kstTrkmEta"    ,&kstTrkmEta);
  RecoB0Tree->SetBranchAddress("kstTrkpPt"     ,&kstTrkpPt);
  RecoB0Tree->SetBranchAddress("kstTrkpPhi"    ,&kstTrkpPhi);
  RecoB0Tree->SetBranchAddress("kstTrkpEta"    ,&kstTrkpEta);
  RecoB0Tree->SetBranchAddress("mmk1"          ,&mmk1);
  RecoB0Tree->SetBranchAddress("mmk2"          ,&mmk2);
  RecoB0Tree->SetBranchAddress("bMass"         ,&bMass);
  RecoB0Tree->SetBranchAddress("bBarMass"      ,&bBarMass);
  RecoB0Tree->SetBranchAddress("dR_mum_trkm"   ,&dR_mum_trkm);
  RecoB0Tree->SetBranchAddress("dR_mup_trkp"   ,&dR_mup_trkp);
  RecoB0Tree->SetBranchAddress("passB0Psi_lmnr",&passB0Psi_lmnr);
  RecoB0Tree->SetBranchAddress("passB0Psi_jpsi",&passB0Psi_jpsi);
  RecoB0Tree->SetBranchAddress("passB0Psi_psip",&passB0Psi_psip);
  
  RecoB0TreeOut->Branch("cos_theta_l"   ,&cos_theta_l    ,   "cos_theta_l/D"   );
  RecoB0TreeOut->Branch("cos_theta_k"   ,&cos_theta_k    ,   "cos_theta_k/D"   );
  RecoB0TreeOut->Branch("phi_kst_mumu"  ,&phi_kst_mumu   ,   "phi_kst_mumu/D"  );
  RecoB0TreeOut->Branch("tagged_mass"   ,&tagged_mass    ,   "tagged_mass/D"   );
  RecoB0TreeOut->Branch("mumuMass"      ,&mumuMass       ,   "mumuMass/D"      );
  RecoB0TreeOut->Branch("mumuMassE"     ,&mumuMassE      ,   "mumuMassE/D"     );
  RecoB0TreeOut->Branch("mumuPhi"	,&mumuPhi        ,   "mumuPhiD"        );
  RecoB0TreeOut->Branch("mumuEta"	,&mumuEta        ,   "mumuEta/D"       );
  RecoB0TreeOut->Branch("kstTrkmPt"	,&kstTrkmPt	 ,   "kstTrkmPt/D"     );
  RecoB0TreeOut->Branch("kstTrkmPhi"	,&kstTrkmPhi	 ,   "kstTrkmPhi/D"    );
  RecoB0TreeOut->Branch("kstTrkmEta"	,&kstTrkmEta	 ,   "kstTrkmEta/D"    );
  RecoB0TreeOut->Branch("kstTrkpPt"	,&kstTrkpPt	 ,   "kstTrkpPt/D"     );
  RecoB0TreeOut->Branch("kstTrkpPhi"	,&kstTrkpPhi	 ,   "kstTrkpPhi/D"    );
  RecoB0TreeOut->Branch("kstTrkpEta"	,&kstTrkpEta	 ,   "kstTrkpEta/D"    );
  RecoB0TreeOut->Branch("mmk1"          ,&mmk1           ,   "mmk1/D"          );
  RecoB0TreeOut->Branch("mmk2"          ,&mmk2           ,   "mmk2/D"          );
//
  bool XCut=true;
//
  double x0Cut=-0.4;
  double y0Cut= 0.3;
  double x1Cut= 0.6;
  double y1Cut=-0.1;
//  
  double x_0Cut=3;
  double y_0Cut=3.8;
  double x_1Cut=3.6;
  double y_1Cut=4.8;
  
  double CutX1=3.2;
  double CutX2=3.6;
  double CutY1=4.7;
  double CutY2=4.9;
//
  int nentries = (int)RecoB0Tree->GetEntries();
  
   for (Int_t i=0;i<nentries;i++) { 
    RecoB0Tree->GetEntry(i);
    recQ2         =  mumuMass*mumuMass  ;
//    double theBMass = tagged_mass;
//     if (!(recQ2> 8.68&&recQ2<10.09)&&
//         !(recQ2>12.86&&recQ2<14.18)){
    if(Q2Bin==4){
     passB0Psi =passB0Psi_jpsi;
    }else if(Q2Bin==6){  
     passB0Psi =passB0Psi_psip;
    }else{  
     passB0Psi =passB0Psi_lmnr;
    }  
	
      if(passB0Psi){

       if(tagged_mass>=XMinSign&&tagged_mass<=XMaxSign){
//	if(dR_mum_trkm>0.0001&&
//	   dR_mup_trkp>0.0001
//	  ){
//
//         xMass.setValue(tagged_mass);
         HxMass  ->Fill(tagged_mass);
//         dataMass->addEvent();
         if(recQ2>Q2Min&&recQ2<Q2Max){
     	  if(cos_theta_l>=XMinCosThetaL&&cos_theta_l<=XMaxCosThetaL&&cos_theta_k>=XMinCosThetaK&&cos_theta_k<=XMaxCosThetaK){
	   if(Q2Bin==4){
 	     ROOT::Math::PtEtaPhiMVector JpsiVec (mumuPt,mumuEta,mumuPhi,mumuMass);

 	     ROOT::Math::PtEtaPhiMVector PipVec (kstTrkpPt,kstTrkpEta,kstTrkpPhi,piMass);
 	     ROOT::Math::PtEtaPhiMVector PimVec (kstTrkmPt,kstTrkmEta,kstTrkmPhi,piMass);
 
 	     ROOT::Math::PtEtaPhiMVector KpVec (kstTrkpPt,kstTrkpEta,kstTrkpPhi,kMass);
 	     ROOT::Math::PtEtaPhiMVector KmVec (kstTrkmPt,kstTrkmEta,kstTrkmPhi,kMass);
	     mmpipi = (JpsiVec+PipVec+PimVec).mass();
	     mmpip  = (JpsiVec+PipVec).mass();
	     mmpim  = (JpsiVec+PimVec).mass();
	     mmkp   = (JpsiVec+KpVec).mass();
	     mmkm   = (JpsiVec+KmVec).mass();
	     kstKp  = (KpVec+PimVec).mass();
	     kstKm  = (KmVec+PipVec).mass();
	     pipPt  = (PipVec).pt();
	     pimPt  = (PimVec).pt();
	     mmpipkm = (JpsiVec+PipVec+KmVec).mass();
	     mmpimkp = (JpsiVec+PimVec+KpVec).mass();
	     pipk = (PipVec+KmVec).mass();
	     pimk = (PimVec+KpVec).mass();
 	     if (tagB0>0) {
	      mmpi1=mmpip;
	      mmpi2=mmpim;
	      kst1 = kstKp;
	      kst2 = kstKm;
	      mmpiKaon=mmpipkm;
	      mmpiKaon2=mmpimkp;
	      pi1Pt=pipPt;
	      pi2Pt=pimPt;
	      piKaon = pipk;
	      mmka1=mmkp;
	      mmka2=mmkm;
	     }else{
	      mmpi1=mmpim;
	      mmpi2=mmpip;
	      kst1 = kstKm;
	      kst2 = kstKp;
	      mmpiKaon=mmpimkp;
	      mmpiKaon2=mmpipkm;
	      pi1Pt=pimPt;
	      pi2Pt=pipPt;
	      piKaon = pimk;
	      mmka1=mmkm;
	      mmka2=mmkp;
	     }
//	     XCut =((mmpi2>3.2&&mmpi2<3.6)&&mmka1>4.7&&mmka1<4.9)&&fabs(BpMass-mmpiKaon)<0.139&&pi1Pt>pi2Pt&&(kst2-KstarMass)>0.09;
//	     XCut =fabs(BpMass-mmpiKaon)<0.139&&pi1Pt>pi2Pt&&(kst2-KstarMass)>0.09&&mmpiKaon2<5.15;
	     XCut=(((BpMass-mmpiKaon)-y0Cut)/(y1Cut-y0Cut))<(((kst2-KstarMass)-x0Cut)/(x1Cut-x0Cut))&&pi1Pt>pi2Pt&&(kst2-KstarMass)>0
	          &&(mmpi2>CutX1&&mmpi2<CutX2)&&(mmka1>CutY1&&mmka1<CutY2)&&((mmka1-y_0Cut)/(y_1Cut-y_0Cut))>((mmpi2-x_0Cut)/(x_1Cut-x_0Cut));
	   }  
	   if(!XCut){
     	       RecoB0TreeOut->Fill();
           }
//     	   }
//     	 }
//        }
       }
      } 
     }
//     }
    }
   }
  
   
   
//  char TXT[200];
//  sprintf(TXT,"Mass Reco   Entries = %7f",HxMassQ2->GetEntries());
  cout<<"***********************************"<<endl;
  cout<<"***** RECONSTRUCTED EVENTS ********\n"<<endl;
  cout<<"***********************************\n"<<endl;
  cout<<"RecoB0Tree   Entries      = "<<nentries<<endl;
//  cout<<TXT<<endl;
  cout<<"\n***********************************"<<endl;
  cout<<"***********************************"<<endl;
//
//  
  
  
 


  
//   TH1D* HxRecoX  = (TH1D*) HxReco->ProjectionX("HxRecoX",1,HxReco->GetNbinsY(),1,HxReco->GetNbinsZ());HxRecoX->SetTitle("HxReco Projection Cos#theta_{L}");
//   TH1D* HxRecoY  = (TH1D*) HxReco->ProjectionY("HxRecoY",1,HxReco->GetNbinsX(),1,HxReco->GetNbinsZ());HxRecoY->SetTitle("HxReco Projection Cos#theta_{K}");
//   TH1D* HxRecoZ  = (TH1D*) HxReco->ProjectionZ("HxRecoZ",1,HxReco->GetNbinsX(),1,HxReco->GetNbinsY());HxRecoZ->SetTitle("HxReco Projection #phi");
  
  OutFileNtupla->cd();
  
  c2->Write();
  gSystem->Exec(Form("mv %s %s.tmp",PNGNameMassHist,PNGNameMassHist));
  c2->Print(PNGNameMassHist);
//   HxReco->Write();
//   HxRecoX->Write();
//   HxRecoY->Write();
//   HxRecoZ->Write();
  HxMass->Write();
  pdfHxMass->Write();
  sigHxMass->Write();
  bkgHxMass->Write();
  HxMassQ2->Write();
  RecoB0TreeOut->Write();
  OutFileNtupla->Close();
  
//   sprintf(testo,"mv %s %s.tmp",PDFNameRecoHisto,PDFNameRecoHisto);
//   gSystem->Exec(testo);
//   c3->cd(1);HxRecoX->Draw();
//   c3->cd(2);HxRecoY->Draw();
//   c3->cd(3);HxRecoZ->Draw();
//   c3->Print(PDFNameRecoHisto);
  cout<<"**********************************************************************\n"<<endl;
  cout<<"save  HxReco  in "<<OutFileNameInputHisto<<"\n"<<endl;
  cout<<"**********************************************************************\n"<<endl;
}
//==========================================================================================
//
//       FitMassSpectrumRoofit
//
//==========================================================================================
double FitMassSpectrumRoofit(RooDataSet* data, TCanvas* c2, TH1D* masHist, TH1D* pdfHist, TH1D*sigHist, TH1D* bkgHist, int MaxDegreeBckg){
// //
//   if(MaxDegreeBckg<=0) {
//    cout<<"**********************************************************************\n"<<endl;
//    cout<<"Error!! MaxDegree <=0 in   FitMassSpectrumRoofit		       *\n"<<endl;
//    cout<<"**********************************************************************\n"<<endl;
//   }
// //

    TFile* fitMassFile = new TFile( fitMassFileName, "READ" );
    if ( !fitMassFile || !fitMassFile ->IsOpen() ) {
      cout<<Form("File not found: %s\n",fitMassFileName)<<endl;
      exit(1);
    }
     RooWorkspace* w = (RooWorkspace*)fitMassFile->Get("w");
    if ( !w || w->IsZombie() ) {
     cout<<Form("Workspace not found in file:%s\n",fitMassFileName)<<endl;
     exit(1);
    } else {
     cout<<Form("Workspace Found!!! In file : %s\n",fitMassFileName)<<endl;
    }
//
    if(!(w->loadSnapshot(Form("reference_fit_RT_%d",Q2Bin)))){
      cout<<Form("Snapshot %s Workspace not found!!!\n",Form("reference_fit_RT_%d",Q2Bin))<<endl;
      exit(1);
    }else{
      w->loadSnapshot(Form("reference_fit_RT_%d",Q2Bin));
//      w->cd(Form("reference_fit_RT_%d",Q2Bin));
      cout<<Form("Snapshot %s Workspace found...\n",Form("reference_fit_RT_%d",Q2Bin))<<endl;
      cout<<"=========================================================================="<<endl;
      cout<<"=========================================================================="<<endl;
      cout<<Form(" DUMP WORKSPACE IN reference_fit_RT_%d",Q2Bin)<<endl;
      cout<<"=========================================================================="<<endl;
      cout<<Form("Snapshot %s Workspace found...\n",Form("reference_fit_RT_%d",Q2Bin))<<endl;
      w->Print("V");
      cout<<"=========================================================================="<<endl;
      cout<<"=========================================================================="<<endl;
    };

    w->ls();
    RooRealVar* tagged_mass= w->var("tagged_mass");
    double tagged_mass_rangeValMin=tagged_mass->getMin();
    double tagged_mass_rangeValMax=tagged_mass->getMax();
    std::cout<<Form("From workspace read tagged_mass in range [%f-%f]\n",tagged_mass_rangeValMin,tagged_mass_rangeValMax) <<std::endl;
    if( (tagged_mass_rangeValMin!=tagged_mass_rangeMin) || (tagged_mass_rangeValMax!=tagged_mass_rangeMax) ){
      std::cout<<Form("Warning! Force setting tagged_mass in range [%f-%f]\n",tagged_mass_rangeMin,tagged_mass_rangeMax) <<std::endl;
     tagged_mass->setRange(tagged_mass_rangeMin,tagged_mass_rangeMax);
    } 
//    tagged_mass->setRange(XMinSign,XMaxSign);
//    if(Q2Bin==0)  tagged_mass_rangeMin = 4.9;
    tagged_mass->setRange("full",tagged_mass_rangeMin,tagged_mass_rangeMax);
   
//    tagged_mass->setRange(XMinFull,XMaxFull);
    
//     RooRealVar *nsig_ref    = w->var("Yield");
//     RooRealVar *nbkg_ref    = w->var("nbkg");
// //
//     if(!nsig_ref){
//       cout<<"Yield from ref. fit  not found!!!\n"<<endl;
//       exit(1);
//     }else{
//       cout<<Form("Yield from ref. fit found = %f\n",nsig_ref->getVal())<<endl;
//     }
//     if(!nbkg_ref){
//       cout<<"BackYield form ref not found!!!\n"<<endl;
//       exit(1);
//     }else{
//       cout<<Form("BackYield  found = %f\n",nbkg_ref->getVal())<<endl;
//     }
   
    
    RooRealVar *yield_fromMC_RT    = w->var(Form("nRT_%d",Q2Bin));
    RooRealVar *yield_fromMC_WT    = w->var(Form("nWT_%d",Q2Bin));
//
    if(!yield_fromMC_RT){
      cout<<"yield_fromMC_RT  not found!!!\n"<<endl;
      exit(1);
    }else{
      cout<<Form("yield_fromMC_RT  found = %f\n",yield_fromMC_RT->getVal())<<endl;
    }
    if(!yield_fromMC_WT){
      cout<<"yield_fromMC_WT  not found!!!\n"<<endl;
      exit(1);
    }else{
      cout<<Form("yield_fromMC_WT  found = %f\n",yield_fromMC_WT->getVal())<<endl;
    }
    RooRealVar * fraction = new RooRealVar("fraction","fraction",0.,1.);
    fraction->setVal(yield_fromMC_RT->getVal()/(yield_fromMC_RT->getVal()+yield_fromMC_WT->getVal()));
    fraction->setError(sqrt(
    pow(yield_fromMC_RT->getError()*yield_fromMC_WT->getVal()/(yield_fromMC_RT->getVal()+yield_fromMC_WT->getVal())/(yield_fromMC_RT->getVal()+yield_fromMC_WT->getVal()),2)+
    pow(yield_fromMC_WT->getError()*yield_fromMC_RT->getVal()/(yield_fromMC_RT->getVal()+yield_fromMC_WT->getVal())/(yield_fromMC_RT->getVal()+yield_fromMC_WT->getVal()),2)));
    cout<<Form("fraction mistagged = %f +/- %f\n",fraction->getVal(),fraction->getError())<<endl;
//    
    RooRealVar    * mean_rt      = 0;
    RooAbsPdf     *theRTgauss    = 0;
    RooRealVar	  * sigma_rt1	 = 0;
    RooRealVar	  * sigma_rt2	 = 0;
    RooGaussian   * c_sigma_rt1  = 0;
    RooGaussian   * c_sigma_rt2  = 0;
    RooGaussian   * c_mean_rt	 = 0;
    RooGaussian   * c_f1rt	 = 0;
    RooProdPdf	  * c_RTgauss	 = 0;
    RooRealVar	  * alpha_rt1	 = 0;
    RooRealVar	  * alpha_rt2	 = 0;
    RooRealVar	  * n_rt1	 = 0;
    RooRealVar	  * n_rt2	 = 0;
    RooRealVar	  * f1rt	 = 0;
    RooArgSet	  * c_vars	 = 0;
    RooGaussian   * c_alpha_rt1  = 0;
    RooGaussian   * c_alpha_rt2  = 0;
    RooGaussian   * c_n_rt1	 = 0;
    RooGaussian   * c_n_rt2	 = 0;
//    RooFormulaVar * deltaPeaks   = 0;
//    RooGaussian   * c_deltaPeaks = 0;
    RooArgList    * c_pdfs       = 0;
    RooArgList    * c_pdfs_rt    = 0;
    RooArgList    * c_pdfs_wt    = 0;
    RooRealVar	  *mean_wt	 = 0;
    RooRealVar	  *sigma_wt	 = 0;
    RooRealVar	  *alpha_wt1	 = 0;
    RooRealVar	  *alpha_wt2	 = 0;
    RooRealVar	  *n_wt1	 = 0;
    RooRealVar	  *n_wt2	 = 0;




// RT Double Gaussian    
    if(w->pdf(Form("doublegaus_RT%d",Q2Bin))) theRTgauss   = w->pdf(Form("doublegaus_RT%d",Q2Bin));
    if(theRTgauss){
       cout<<Form("FitMassSpectrumRoofit: doublegaus_RT%d",Q2Bin)<<endl;
       mean_rt     = w->var(Form("mean^{RT%d}",Q2Bin));
       sigma_rt1   = w->var(Form("#sigma_{1}^{RT%d}",Q2Bin));
       sigma_rt2   = w->var(Form("#sigma_{2}^{RT%d}",Q2Bin));
       f1rt	   = w->var(Form("f^{RT%d}",Q2Bin));
       c_sigma_rt1 = _constrainVar(sigma_rt1,w);
       c_sigma_rt2 = _constrainVar(sigma_rt2,w);
       c_mean_rt   = _constrainVar(mean_rt,w);
       c_f1rt	   = _constrainVar(f1rt,w);
 
       ////// creating constraints for the RT component
//       c_RTgauss = new RooProdPdf("c_RTgauss" , "c_RTgauss" , RooArgList(*theRTgauss,*c_sigma_rt1,*c_sigma_rt2,*c_mean_rt,*c_f1rt  ) );
       c_pdfs    = new RooArgList(*c_sigma_rt1,*c_sigma_rt2,*c_mean_rt,*c_f1rt);
       c_pdfs_rt = new RooArgList(*theRTgauss,*c_sigma_rt1,*c_sigma_rt2,*c_mean_rt,*c_f1rt);
       c_vars    = new RooArgSet( *sigma_rt1, *sigma_rt2, *mean_rt ,*f1rt);
    }else{
       cout<<Form("FitMassSpectrumRoofit: doublecb_RT%d",Q2Bin)<<endl;
       theRTgauss  = w->pdf(Form("doublecb_RT%d",Q2Bin));
       mean_rt     = w->var(Form("mean_{RT}^{%d}",Q2Bin));
       sigma_rt1   = w->var(Form("#sigma_{RT1}^{%d}",Q2Bin));
       alpha_rt1   = w->var(Form("#alpha_{RT1}^{%d}",Q2Bin));
       alpha_rt2   = w->var(Form("#alpha_{RT2}^{%d}",Q2Bin));
       n_rt1	   = w->var(Form("n_{RT1}^{%d}",Q2Bin));
       n_rt2	   = w->var(Form("n_{RT2}^{%d}",Q2Bin));
       f1rt	   = w->var(Form("f^{RT%d}",Q2Bin));
//      deltaPeaks  = RooFormulaVar("deltaPeaks%s"%ibin, "@0 - @1", RooArgList(mean_rt, mean_wt));  
       c_sigma_rt1 = _constrainVar(sigma_rt1,w);
       c_alpha_rt1 = _constrainVar(alpha_rt1,w);
       c_alpha_rt2 = _constrainVar(alpha_rt2,w);
//    
//  workaround
//    
//       if(Q2Bin==7){
//        n_rt2->setMax(60);
//       }  
       c_n_rt1     = _constrainVar(n_rt1,w);
       c_n_rt2     = _constrainVar(n_rt2,w);
       if (Q2Bin < 4){ 
           cout<<Form("FitMassSpectrumRoofit: fast doublecb_RT%d",Q2Bin)<<endl;
           c_pdfs    = new RooArgList(*c_sigma_rt1, *c_alpha_rt1, *c_alpha_rt2, *c_n_rt1, *c_n_rt2);
           c_pdfs_rt = new RooArgList(*theRTgauss,*c_sigma_rt1, *c_alpha_rt1, *c_alpha_rt2, *c_n_rt1, *c_n_rt2);
           c_vars    = new RooArgSet( *sigma_rt1,      *alpha_rt1,    *alpha_rt2,    *n_rt1,    *n_rt2);
       }else{
           cout<<Form("FitMassSpectrumRoofit: old doublecb_RT%d",Q2Bin)<<endl;
           sigma_rt2 = w->var(Form("#sigma_{RT2}^{%d}",Q2Bin));
           c_sigma_rt2   = _constrainVar(sigma_rt2, w);
           c_pdfs    = new RooArgList(*c_sigma_rt1, *c_sigma_rt2, *c_alpha_rt1, *c_alpha_rt2, *c_n_rt1, *c_n_rt2);
           c_pdfs_rt = new RooArgList(*theRTgauss,*c_sigma_rt1, *c_sigma_rt2, *c_alpha_rt1, *c_alpha_rt2, *c_n_rt1, *c_n_rt2);
           c_vars    = new RooArgSet(   *sigma_rt1,    *sigma_rt2,    *alpha_rt1,    *alpha_rt2,    *n_rt1,    *n_rt2);
       }
    } 
//
    c_RTgauss = new RooProdPdf("c_RTgauss" , "c_RTgauss" , *c_pdfs_rt );
//
    if(!(w->loadSnapshot(Form("reference_fit_WT_%d",Q2Bin)))){
      cout<<Form("Snapshot %s Workspace not found!!!\n",Form("reference_fit_WT_%d",Q2Bin))<<endl;
      exit(1);
    }else{
      w->loadSnapshot(Form("reference_fit_WT_%d",Q2Bin));
//      w->cd(Form("reference_fit_WT_%d",Q2Bin));
      
      cout<<"=========================================================================="<<endl;
      cout<<"=========================================================================="<<endl;
      cout<<Form(" DUMP WORKSPACE IN reference_fit_WT_%d",Q2Bin)<<endl;
      cout<<"=========================================================================="<<endl;
      cout<<Form("Snapshot %s Workspace found...\n",Form("reference_fit_WT_%d",Q2Bin))<<endl;
      w->Print("V");
      cout<<"=========================================================================="<<endl;
      cout<<"=========================================================================="<<endl;
    };
    
//    
    mean_wt	= !(w->var(Form("mean^{WT%d}",Q2Bin)))?(w->var(Form("mean_{WT}^{%d}",Q2Bin))):(w->var(Form("mean^{WT%d}",Q2Bin)));
    sigma_wt	= !(w->var(Form("#sigma_{CB}^{WT%d}",Q2Bin)))?(w->var(Form("#sigma_{WT1}^{%d}",Q2Bin))):(w->var(Form("#sigma_{CB}^{WT%d}",Q2Bin)));
    alpha_wt1	= !(w->var(Form("#alpha_{1}^{WT%d}",Q2Bin)))?(w->var(Form("#alpha_{WT1}^{%d}",Q2Bin))):(w->var(Form("#alpha_{1}^{WT%d}",Q2Bin)));
    alpha_wt2	= !(w->var(Form("#alpha_{2}^{WT%d}",Q2Bin)))?(w->var(Form("#alpha_{WT2}^{%d}",Q2Bin))):(w->var(Form("#alpha_{2}^{WT%d}",Q2Bin)));
    n_wt1	= !(w->var(Form("n_{1}^{WT%d}",Q2Bin)))?(w->var(Form("n_{WT1}^{%d}",Q2Bin))):(w->var(Form("n_{1}^{WT%d}",Q2Bin)));
    n_wt2	= !(w->var(Form("n_{2}^{WT%d}",Q2Bin)))?(w->var(Form("n_{WT2}^{%d}",Q2Bin))):(w->var(Form("n_{2}^{WT%d}",Q2Bin)));
//    
//  workaround
//    
//     if(Q2Bin==7){
//       n_wt2->setMax(150);
//     }  
      
//    
    RooAbsPdf *theWTgauss = w->pdf(Form("doublecb_%d",Q2Bin));
    if(!theWTgauss)  {
     cout<<"pdf theWTgauss not found!!!\n"<<endl;
     exit(1);
    } else{
     cout<<Form("pdf %s  found...\n",theWTgauss->GetName())<<endl;
    }
    RooGaussian* c_mean_wt     = _constrainVar(mean_wt, w);
    RooGaussian* c_sigma_wt    = _constrainVar(sigma_wt, w);
    RooGaussian* c_alpha_wt1   = _constrainVar(alpha_wt1, w);
    RooGaussian* c_alpha_wt2   = _constrainVar(alpha_wt2, w);
    RooGaussian* c_n_wt1       = _constrainVar(n_wt1, w);
    RooGaussian* c_n_wt2       = _constrainVar(n_wt2, w);
//   if(w->obj( "deltaPeaks")) {
//    if(w->obj( Form("deltaPeaks%d",Q2Bin))) {
//     deltaPeaks  = (RooFormulaVar *)w->obj(Form("deltaPeaks%d",Q2Bin));
//     c_deltaPeaks = new RooGaussian(Form("deltaPeaks%d",Q2Bin) , "c_deltaPeaks", *deltaPeaks, RooConst( deltaPeaks->getVal() ), RooConst( 0.0005 )); // value to be checked
// //    c_vars->add(*deltaPeaks);
// //    c_pdfs ->add(*c_deltaPeaks);
//    }else if(w->pdf(Form("doublecb_RT%d",Q2Bin))){
//     cout<<Form("FitMassSpectrumRoofit: deltaPeaks NOT FOUND!!!")<<endl;
// //    exit(0);
//    }	    

 //     RooGaussian* c_sigma_wt2   = _constrainVar(sigma_wt2, w);
//     RooGaussian* c_f3          = _constrainVar(f3wt, w);
    

    ////// creating constraints for the WT component
//    RooProdPdf* c_WTgauss  = new RooProdPdf("c_WTgauss" , "c_WTgauss",RooArgList(*theWTgauss,*c_alpha_wt1,*c_n_wt1,*c_sigma_wt,*c_mean_wt,*c_alpha_wt2,*c_n_wt2  ) );     
    RooRealVar  frt("F_{RT}"			  , "frt"   , fraction->getVal() , 0, 1);
    RooGaussian c_frt("c_frt"           	   , "c_frt" , frt,  RooFit::RooConst(fraction->getVal()) , RooFit::RooConst(fraction->getError()) );
//    RooAddPdf	signalFunction("sumgaus"	  , "rt+wt" , RooArgList(*theRTgauss,*theWTgauss), RooArgList(frt));
//     RooGaussian c_frt("c_frt"           	   , "c_frt" , frt,  (fraction) , (*fraction_s) );
//    RooGaussian c_frt("c_frt"           	   , "c_frt" , frt,  RooFit::RooConst(0.87628877977) , RooFit::RooConst(0.000523435458235) );
 
    c_pdfs_wt = new RooArgList(*theWTgauss);
    c_pdfs_wt->add(*c_sigma_wt);
    c_pdfs_wt->add(*c_mean_wt);
    c_pdfs_wt->add(*c_alpha_wt1);
    c_pdfs_wt->add(*c_alpha_wt2);
    c_pdfs_wt->add(*c_n_wt1);
    c_pdfs_wt->add(*c_n_wt2);


    RooProdPdf* c_WTgauss  = new RooProdPdf("c_WTgauss" , "c_WTgauss",*c_pdfs_wt);
    RooAddPdf	signalFunction("sumgaus"	  , "rt+wt" , RooArgList(*c_RTgauss,*c_WTgauss), RooArgList(frt));
    RooProdPdf  c_signalFunction("c_signalFunction", "c_signalFunction", RooArgList(signalFunction, c_frt))   ;  
    c_pdfs->add(*c_sigma_wt);
    c_pdfs->add(*c_mean_wt);
    c_pdfs->add(*c_alpha_wt1);
    c_pdfs->add(*c_alpha_wt2);
    c_pdfs->add(*c_n_wt1);
    c_pdfs->add(*c_n_wt2);
    c_pdfs->add(c_frt);
//    c_pdfs->add(signalFunction);

    c_vars->add(*sigma_wt);
    c_vars->add(*mean_wt);
    c_vars->add(*alpha_wt1);
    c_vars->add(*alpha_wt2);
    c_vars->add(*n_wt1);
    c_vars->add(*n_wt2);
    c_vars->add(frt);

////// now create background parametrization
    RooRealVar*  slope= new RooRealVar("slope"      , "slope"           ,    0.5,   -10, 10);
//    RooExponential bkg_exp("bkg_exp"    , "exponential"     ,  *slope,   *tagged_mass  );
//     RooRealVar     pol_c1("p1"          , "coeff x^0 term"  ,    0.5,   -10, 10);
//     RooRealVar     pol_c2("p2"          , "coeff x^1 term"  ,    0.5,   -10, 10);
//     RooRealVar     pol_c3("p3"          , "coeff x^2 term"  ,    0.5,   -10, 10);
//     RooRealVar     pol_c4("p4"          , "coeff x^3 term"  ,    0.5,   -10, 10);
// 
//     RooChebychev   bkg_exp("bkg_exp"    , "2nd order pol"   ,  *tagged_mass, RooArgList(pol_c1,pol_c2,pol_c3,pol_c4));
    RooRealVar*     pol_b0= new RooRealVar("pol_b0"          , "b0"  ,    0.5,  0., 1.);
    RooRealVar*     pol_b1= new RooRealVar("pol_b1"          , "b1"  ,    0.1,  0., 1.);
    RooRealVar*     pol_b2= new RooRealVar("pol_b2"          , "b2"  ,    0.1,  0., 1.);
    RooRealVar*     pol_b3= new RooRealVar("pol_b3"          , "b3"  ,    0.01 , 0., 1.);
    RooRealVar*     pol_b4= new RooRealVar("pol_b4"          , "b4"  ,    0.1 , 0., 1.);
    if(Q2Bin!=4){
     bkg_exp = new RooExponential("bkg_exp"    , "exponential"     ,  *slope,   *tagged_mass  );
    }else{
//     pol_b0->setConstant(kTRUE);
//     pol_b4->setConstant(kTRUE);
//     pol_b3->setConstant(kTRUE);
     bkg_exp = new RooBernstein("bkg_exp"    , "bernstein pol"  ,  *tagged_mass, RooArgList(*pol_b0,*pol_b1,*pol_b2,*pol_b3,*pol_b4));
    }
    
    int NCPU=1;
    if(NFactGen>1) NCPU=10;
    double yieldIni = NFactGen*1000;
    double backgIni = NFactGen*1000;
    double yieldMin = 0.;
    double backgMin = 0.;
    double yieldMax = NFactGen*1000000.;
    double backgMax = NFactGen*1000000.;
   if(Q2Bin==4) {
       yieldIni = yieldSignal;
       backgIni = yieldBckg;
       yieldMin = 100000.;
       backgMin = 10000.;
       yieldMax = 2000000.;
       backgMax = 1000000.;
       NCPU=60;
    }   
   
    RooRealVar     nsig("Yield"         , "signal frac"    ,   yieldIni,     yieldMin, yieldMax  );
    RooRealVar     nbkg("nbkg"          , "bkg fraction"   ,   backgIni,     backgMin, backgMax  );

//    RooRealVar     nsig("Yield"         , "signal frac"    ,    4000,     0,   1000000);
//    RooRealVar     nbkg("nbkg"          , "bkg fraction"   ,    1000,     0,   550000);
// 
//    RooProdPdf  c_signalFunction("c_signalFunction", "c_signalFunction", RooArgList(signalFunction, c_frt))   ;  
//    RooProdPdf c_signalFunction("c_signalFunction", "c_signalFunction", *c_pdfs);
    RooAddPdf fitFunction("fitfunction" , "fit function"  ,  RooArgList(c_signalFunction, *bkg_exp), RooArgList(nsig, nbkg));
//
//    RooAddPdf fitFunction("fitfunction" , "fit function"  ,  RooArgList(signalFunction, bkg_pol), RooArgList(nsig, nbkg));
//     tagged_mass->setRange("fullRedefined",XMinSBL,XMaxSBR);
    RooFitResult* r = fitFunction.fitTo(*data, 
    		       RooFit::Extended(kTRUE), 
    		       RooFit::NumCPU(NCPU),
    		       RooFit::Save(), 
    		       RooFit::Range("full"), 
    		       RooFit::Verbose(kFALSE),
    		       RooFit::Constrain(*c_vars)
    		      );
     		      
   r->Print();		      
//
//
// save a clone of bkg_exp 
//

    if(Q2Bin!=4){
     bkg_mass_sb = (RooExponential*)bkg_exp->clone(Form("bkg_mass_sb_bin%d_%s",Q2Bin,RunEra) );
    }else{ 
     bkg_mass_sb = (RooBernstein*)  bkg_exp->clone(Form("bkg_mass_sb_bin%d_%s",Q2Bin,RunEra) );
     bkg_mass_sb->setNormRange("full");
    } 
//   RooAbsBinning binning = (tagged_mass->getBinning("full")) ;
//     data->getRange(*tagged_mass,tagged_mass_rangeMin,tagged_mass_rangeMax);
//    std::cout<<Form("%f<[fit mass range]<%f",tagged_mass_rangeMin,tagged_mass_rangeMax)<<std::endl; 
//    tagged_mass_rangeMin=tagged_mass->getMin("full");
//    tagged_mass_rangeMax=tagged_mass->getMax("full");
//    std::cout<<Form("%f<[fit mass range]<%f",tagged_mass_rangeMin,tagged_mass_rangeMax)<<std::endl; 

//      cout<<"----"<<endl;
//       exit(1);
//      data->Print("V");
//      c2->cd();
//      RooPlot* frame = tagged_mass->frame( );
//      data->plotOn(frame, Binning(35), MarkerSize(.7));
//     fitFunction.plotOn(frame);
//     drawPdfComponents(fitFunction, frame, ROOT.kAzure, RooFit.NormRange("full"), RooFit.Range("full"), isData = True);
// 
//     fitFunction.paramOn(frame,  RooFit.Layout(0.62,0.86,0.88));
//     frame.Draw();
//     niceFrame(frame, '')
//     frame. addObject(_writeFitStatus(r))
// 
//     if not args.year=='test':  writeCMS(frame, args.year, [ q2binning[ibin], q2binning[ibin+1] ])
//     frame.Draw()
 //    c2->Print("test.pdf");
 
     double B0SigmaTemp=0.;
     c2->cd();
     TLegend* leg_sign = new TLegend(0.30,0.48,0.90,0.90);
     leg_sign->SetTextSize(0.025) ;
     leg_sign->SetTextAlign(31);
     leg_sign->SetBorderSize(0.);
     leg_sign->SetFillStyle(0);
     leg_sign->SetHeader("B^{0} mass spectrum  Fit Projection");
     if(nsig.getError()!=0){
       leg_sign->AddEntry(masHist ,Form( "Yield_{Sign} =     %5.0f  #pm %5.0f",nsig.getVal(),nsig.getError()),"");
     }else{
       leg_sign->AddEntry(masHist ,Form( "Yield_{Sign} =     %5.0f Fixed",nsig.getVal()),"");
     }
     if(nbkg.getError()!=0){
       leg_sign->AddEntry(masHist ,Form( "Yield_{Bckg} =     %5.0f  #pm  %5.0f",nbkg.getVal(),nbkg.getError()),"");
     }else{
       leg_sign->AddEntry(masHist ,Form( "Yield_{Bckg} =     %5.0f  Fixed",nbkg.getVal()),"");
     }
     if(mean_rt->getError()!=0){
      leg_sign->AddEntry(masHist ,Form( "M_{B^{0}}[RT] =   %5.5f  #pm %5.5f",mean_rt->getVal(),mean_rt->getError()),"");
     }else{
      leg_sign->AddEntry(masHist ,Form( "M_{B^{0}}[RT] =   %5.5f Fixed",mean_rt->getVal()),"");
      }
//     if(mean_wt==0){exit(0);};
     if(mean_wt->getError()!=0){
      leg_sign->AddEntry(masHist ,Form( "M_{B^{0}}[WT] =   %5.5f  #pm %5.5f",mean_wt->getVal(),mean_wt->getError()),"");
     }else{
      leg_sign->AddEntry(masHist ,Form( "M_{B^{0}}[WT] =   %5.5f Fixed",mean_wt->getVal()),"");
      }
     if(sigma_rt1->getError()!=XStepMinuit){
      leg_sign->AddEntry(masHist ,Form( "#sigma_{1}^{RT%d}#scale[0.6]{1}_{B^{0}} =   %5.5f  #pm %5.5f",Q2Bin,sigma_rt1->getVal(),sigma_rt1->getError()),"");
     }else{
      leg_sign->AddEntry(masHist ,Form( "#sigma_{1}^{RT%d}#scale[0.6]{1}_{B^{0}} =   %5.5f Fixed",Q2Bin,sigma_rt1->getVal()),"");
     }
     if(sigma_rt2!=0){
      if(sigma_rt2->getError()!=XStepMinuit){
      leg_sign->AddEntry(masHist ,Form( "#sigma_{2}^{RT%d}#scale[0.6]{2}_{B^{0}} =   %5.5f  #pm %5.5f",Q2Bin,sigma_rt2->getVal(),sigma_rt2->getError()),"");
//   	}else{
//    	 leg_sign->AddEntry(masHist ,Form( "#sigma#scale[0.6]{2}_{B^{0}} =   %5.5f Fixed",sigma2.getVal()),"");
      }
     } 
     if(sigma_wt->getError()!=XStepMinuit){
      leg_sign->AddEntry(masHist ,Form( "#sigma^{WT%d}#scale[0.6]{1}_{B^{0}} =   %5.5f  #pm %5.5f",Q2Bin,sigma_wt->getVal(),sigma_wt->getError()),"");
     }else{
      leg_sign->AddEntry(masHist ,Form( "#sigma^{WT%d}#scale[0.6]{1}_{B^{0}} =   %5.5f Fixed",Q2Bin,sigma_wt->getVal()),"");
     }
//      if(sigma_wt2){
//       if(sigma_wt2->getError()!=XStepMinuit){
//       leg_sign->AddEntry(masHist ,Form( "#sigma_{2}^{WT%d}#scale[0.6]{2}_{B^{0}} =   %5.5f  #pm %5.5f",Q2Bin,sigma_wt2->getVal(),sigma_wt2->getError()),"");
// //   	}else{
// //    	 leg_sign->AddEntry(masHist ,Form( "#sigma#scale[0.6]{2}_{B^{0}} =   %5.5f Fixed",sigma2.getVal()),"");
//       }
//     } 
     double min_CBGaus_rt=0;
     double max_CBGaus_rt=0;
//      float x1zoom=5.1;
//      float x2zoom=5.4;
     float x1zoom=XMinSign;
     float x2zoom=XMaxSign;
//      RooPlot *rframe = tagged_mass->frame(Title("Signal models RT"));
//      RooPlot *wframe = tagged_mass->frame(Title("Signal models WT"));
     RooAbsPdf *gaussRT_study = 0;
     RooGaussian *gaussRT_study1 = 0;
     RooGaussian *gaussRT_study2 = 0;
     TF1 * Func_theRTgauss    = 0;
     TF1 * Func_gaussRT_study = 0;
     TF1 * Clone_gaussRT_study = 0;
//     RooRealVar  f3("f3","f3",0.);
     if(w->pdf(Form("doublecb_RT%d",Q2Bin))){
      if(Q2Bin<4){
       min_CBGaus_rt =  mean_rt->getVal()-alpha_rt1->getVal()*sigma_rt1->getVal();
       max_CBGaus_rt =  mean_rt->getVal()+alpha_rt2->getVal()*sigma_rt1->getVal();
       gaussRT_study = new RooGaussian("gaussRT_study","gauss RT study"    ,*tagged_mass,*mean_rt,*sigma_rt1);
      }else{
//       sigma_rt2_pos = new RooRealVar("sigma_rt2_pos","|(sigma_rt2)|",fabs(sigma_rt2->getVal()));
       min_CBGaus_rt =  mean_rt->getVal()-fabs(alpha_rt1->getVal())*sigma_rt1->getVal();
       max_CBGaus_rt =  mean_rt->getVal()+fabs(alpha_rt2->getVal())*sigma_rt2->getVal();
       gaussRT_study1 = new RooGaussian("gaussRT_study1","gauss RT study1"    ,*tagged_mass,*mean_rt,*sigma_rt1);
       gaussRT_study2 = new RooGaussian("gaussRT_study2","gauss RT study2"    ,*tagged_mass,*mean_rt,*sigma_rt2);
       gaussRT_study = new RooAddPdf("gaussRT_study","gauss RT study"    ,RooArgList(*gaussRT_study1,*gaussRT_study2),RooArgList(*f1rt));
      } 
      std::cout<<Form("RT Double CB Gaus = %f<mass<%f",min_CBGaus_rt,max_CBGaus_rt)<<std::endl;
      csignstudy->cd(1);
      gPad->SetLeftMargin(0.15);
      Func_theRTgauss	 = theRTgauss	->asTF( RooArgList(*tagged_mass) );
      Func_gaussRT_study = gaussRT_study->asTF( RooArgList(*tagged_mass) );
      Func_theRTgauss->SetTitle("RT Model");
//      rframe->GetYaxis()->SetTitleOffset(1.4);
//      theRTgauss->plotOn(rframe,LineColor(kRed));
//      gaussRT_study->plotOn(rframe,LineColor(kBlue));
//      theRTgauss->plotOn(rframe, Range(min_CBGaus_rt,max_CBGaus_rt,kFALSE),FillColor(kRed),DrawOption("F"),FillStyle(3013),VLines());
      TLegend* leg_rt = new TLegend(0.70,0.70,0.90,0.90);
      leg_rt->SetTextSize(0.025) ;
      leg_rt->SetTextAlign(31);
      leg_rt->SetBorderSize(0.);
      leg_rt->SetFillStyle(0);
      leg_rt->AddEntry(Func_theRTgauss ,Form("#color[2]{Double CB Model}"),"");
      leg_rt->AddEntry(Func_gaussRT_study,Form("#color[4]{(2)Gaussian  Model}"),"");
//      rframe->addObject(leg_rt);
//      rframe->Draw();
      Func_theRTgauss->SetLineColor(kRed);
      Func_gaussRT_study->SetLineColor(kBlue);
      Func_theRTgauss->SetLineWidth(1.);
      Func_gaussRT_study->SetLineWidth(1.);
//     Func_gaussRT_study->SetLineStyle(kDashed);
      Func_theRTgauss->SetRange(x1zoom,x2zoom);
      Func_gaussRT_study->SetRange(x1zoom,x2zoom);
      Func_theRTgauss->Draw();
      Func_gaussRT_study->Draw("SAME");
      Clone_gaussRT_study =  (TF1*)Func_gaussRT_study->Clone();
      Clone_gaussRT_study->SetRange(min_CBGaus_rt,max_CBGaus_rt);
      Clone_gaussRT_study->SetFillColor(kBlue);
      Clone_gaussRT_study->SetFillStyle(3013);
      Clone_gaussRT_study->SetLineWidth(1.);
      Clone_gaussRT_study->Draw("SAME FC");
      leg_rt->Draw("SAME");
     } 
 //
     double min_CBGaus_wt =  mean_wt->getVal()-alpha_wt1->getVal()*sigma_wt->getVal();
     double max_CBGaus_wt =  mean_wt->getVal()+alpha_wt2->getVal()*sigma_wt->getVal();
     RooGaussian *gaussWT_study = new RooGaussian("gaussWT_study","gauss WT study"    ,*tagged_mass,*mean_wt,*sigma_wt);
     std::cout<<Form("WT Double CB Gaus = %f<mass<%f",min_CBGaus_wt,max_CBGaus_wt)<<std::endl;
     csignstudy->cd(2);
     gPad->SetLeftMargin(0.15);
//     wframe->GetYaxis()->SetTitleOffset(1.4);
//     theWTgauss->plotOn(wframe,LineColor(kRed));
     TF1 * Func_theWTgauss    = theWTgauss   ->asTF( RooArgList(*tagged_mass) );
     TF1 * Func_gaussWT_study = gaussWT_study->asTF( RooArgList(*tagged_mass) );
      Func_theWTgauss->SetTitle("WT Model");
//     RooAbsReal* IntegtheWTgauss    = theWTgauss->createIntegral(*tagged_mass,*tagged_mass,"full");
//     double scal = sigma_wt->getVal()*sqrt(2*TMath::Pi())*theWTgauss->getVal(RooArgList(*mean_wt));
//     cout<<IntegtheWTgauss->getVal()<<" "<<f->GetMaximum()<<endl;exit(0);
//     double scal = sigma_wt->getVal()*sqrt(2*TMath::Pi())*theWTgauss->getVal(RooArgList(*mean_wt));
//     double scal = IntegtheWTgauss->getVal()/gaussWT_study->getVal(RooArgList(*mean_wt));
//     gaussWT_study->plotOn(wframe,LineColor(kBlue),Normalization(1/IntegtheWTgauss->getVal()));
//     theWTgauss->plotOn(wframe, Range(min_CBGaus_wt,max_CBGaus_wt,kFALSE),FillColor(kRed),DrawOption("F"),FillStyle(3013),VLines());
     TLegend* leg_wt = new TLegend(0.70,0.70,0.90,0.90);
     leg_wt->SetTextSize(0.025) ;
     leg_wt->SetTextAlign(31);
     leg_wt->SetBorderSize(0.);
     leg_wt->SetFillStyle(0);
     leg_wt->AddEntry(Func_theWTgauss ,Form("#color[2]{Double CB Model}"),"");
     leg_wt->AddEntry(Func_gaussWT_study ,Form("#color[4]{Gaussian  Model}"),"");
//     leg_wt->AddEntry(theWTgauss ,Form("#color[2]{Double CB Model}"),"");
//     leg_wt->AddEntry(gaussWT_study ,Form("#color[4]{Gaussian  Model}"),"");
//     wframe->addObject(leg_wt);
//     wframe->Draw();
     Func_theWTgauss->SetLineColor(kRed);
     Func_gaussWT_study->SetLineColor(kBlue);
     Func_theWTgauss->SetLineWidth(1.0);
     Func_gaussWT_study->SetLineWidth(1.0);
//     Func_gaussWT_study->SetLineStyle(kDashed);
     Func_theWTgauss->SetRange(x1zoom,x2zoom);
     Func_gaussWT_study->SetRange(x1zoom,x2zoom);
     Func_theWTgauss->Draw();
     Func_gaussWT_study->Draw("SAME");
     TF1* Clone_gaussWT_study =  (TF1*)Func_gaussWT_study->Clone();
     Clone_gaussWT_study->SetRange(min_CBGaus_wt,max_CBGaus_wt);
     Clone_gaussWT_study->SetFillColor(kBlue);
     Clone_gaussWT_study->SetFillStyle(3013);
     Clone_gaussWT_study->SetLineWidth(1.0);
//     Clone_gaussWT_study->GetXaxis()->SetRangeUser(5.1,5.4);
     Clone_gaussWT_study->Draw("SAME FC");
     leg_wt->Draw("SAME");
     
     char PNGSignStudy[300]="";sprintf(PNGSignStudy,Form("signal-mass-study-Q2Bin-%d.png",Q2Bin));
     gSystem->Exec(Form("mv %s %s.tmp",PNGSignStudy,PNGSignStudy));
     csignstudy->Print(PNGSignStudy);
     
     double B0SigmaRT=0;
     double Sigma1RT =sigma_rt1->getVal();
     if(sigma_rt2!=0){
      double Sigma2RT =sigma_rt2->getVal();
      double WG1=f1rt->getVal();
      leg_sign->AddEntry(masHist ,Form( "f^{RT%d} =   %5.5f  #pm %5.5f",Q2Bin,f1rt->getVal(),f1rt->getError()),"");
      B0SigmaRT = sqrt(Sigma1RT*Sigma1RT*WG1+(1.-WG1)*Sigma2RT*Sigma2RT);
     }else{
      B0SigmaRT = Sigma1RT;
     }

      double B0sigma_wt =sigma_wt->getVal();
//      if(sigma_wt2){
//       double Sigma2WT =sigma_wt2->getVal();
//       double WG1=f3wt->getVal();
//       leg_sign->AddEntry(masHist ,Form( "f^{WT%d} =   %5.5f  #pm %5.5f",Q2Bin,f3wt->getVal(),f3wt->getError()),"");
//       B0sigma_wt = sqrt(Sigma1WT*Sigma1WT*WG1+(1.-WG1)*Sigma2WT*Sigma2WT);
//      }else{
//       B0sigma_wt = Sigma1WT;
//      }
     
     
     B0SigmaTemp = sqrt(B0SigmaRT*B0SigmaRT*fraction->getVal()+(1.-fraction->getVal())*B0sigma_wt*B0sigma_wt);
     

     std::cout<<Form("B0SigmaRT = %f B0sigma_wt = %f B0SigmaTot = %f",B0SigmaRT,B0sigma_wt,B0SigmaTemp) <<std::endl;

//     tagged_mass->setRange("SignLeft" ,mean_rt->getVal(),XMaxSign);
//     tagged_mass->setRange("SignRight",XMinSign,mean_rt->getVal());
//      tagged_mass->setRange(5.1,5.4);
//
     TCanvas* ccdf_signal = new TCanvas("ccdf_signal","cdf",200,10,1200,600);
//     ccdf_signal->Divide(2,2);
     TCanvas* ccdf_rt = new TCanvas("ccdf_rt","cdf RT",200,10,1200,600);
     ccdf_rt->Divide(2,2);
     
     TCanvas* ccdf_signal_zoom = new TCanvas("ccdf_signal_zoom","cdf",200,10,1200,600);
     ccdf_signal_zoom->Divide(2,2);
     TCanvas* ccdf_rt_zoom = new TCanvas("ccdf_rt_zoom","cdf RT",200,10,1200,600);
     ccdf_rt_zoom->Divide(2,2);
     TF1 *  Func_signal = c_signalFunction.asTF( RooArgList(*tagged_mass) );
     RooAbsPdf* signalCdf = (RooAbsPdf*)c_signalFunction.createCdf(*tagged_mass);
     TF1 *  Func_signalCdf= signalCdf->asTF( RooArgList(*tagged_mass) );
//     RooAbsPdf* theRTgaussCdf = (RooAbsPdf*)theRTgauss->createCdf(*tagged_mass);
//     TF1 *  Func_theRTgaussCdf= theRTgaussCdf->asTF( RooArgList(*tagged_mass) );
     std::cout<<"Estimate %% probability for Full/RT Signal:"<<std::endl;
     double integTails= 0;
//     double XLeftLim  = 0;
//     double XRightLim = 0;
//      double SigmaEstL = 0;
//      double SigmaEstR = 0;
     Func_signal->SetLineColor(kRed);
     Func_signal->SetLineWidth(1.0);
//      for(double iSigma=1;iSigma<5;iSigma++){
//      
//        integTails= TMath::Erfc(iSigma/sqrt(2));//  Erfc(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between x and infinity; t=x/sqrt(2)
//        XLeftLim  = Func_signalCdf->GetX(integTails/2.);
//        XRightLim = Func_signalCdf->GetX(1-integTails/2.);
//        SigmaEstL = (fabs(mean_rt->getVal()-XLeftLim))/iSigma;
//        SigmaEstR = (fabs(mean_rt->getVal()-XRightLim))/iSigma;
//        std::cout<<Form("Full signal %f%% [%iSigma gauss sigma] in the Range %f<mass<%f average sigmaL=%f sigmaR=%f",1-integTails,int(iSigma),XLeftLim,XRightLim,SigmaEstL,SigmaEstR)<<std::endl;
//        ccdf_signal->cd(iSigma);
//        Func_signal->SetRange(4.8,5.6);
//        Func_signal->SetMaximum(1.);
//        Func_signal->Draw();
//        TF1* Clone_signal =  (TF1*) Func_signal->Clone();
//        Clone_signal->SetRange(XLeftLim,XRightLim);
//        Clone_signal->SetFillColor(kBlue);
//        Clone_signal->SetFillStyle(3013);
//        Clone_signal->SetLineWidth(1.);
//        Clone_signal->Draw("SAME FC");
// //
//        ccdf_signal_zoom->cd(iSigma);
//        Func_signal->SetRange(4.8,5.6);
//        Func_signal->SetMaximum(0.2);
//        Func_signal->Draw();
//        Clone_signal =  (TF1*) Func_signal->Clone();
//        Clone_signal->SetRange(XLeftLim,XRightLim);
//        Clone_signal->Draw("SAME FC");
//        XLeftLim  = Func_theRTgaussCdf->GetX(integTails/2.);
//        XRightLim = Func_theRTgaussCdf->GetX(1-integTails/2.);
//        SigmaEstL = (fabs(mean_rt->getVal()-XLeftLim))/iSigma;
//        SigmaEstR = (fabs(mean_rt->getVal()-XRightLim))/iSigma;
//        std::cout<<Form("RT   signal %f%% [%iSigma gauss sigma] in the Range %f<mass<%f average sigmaL=%f sigmaR=%f",1-integTails,int(iSigma),XLeftLim,XRightLim,SigmaEstL,SigmaEstR)<<std::endl;
//
//        ccdf_rt->cd(iSigma);
//        Func_theRTgauss->SetRange(5.0,5.6);
//        Func_theRTgauss->SetMaximum(1.);
//        Func_theRTgauss->Draw();
//        TF1* Clone_theRTgauss =  (TF1*) Func_theRTgauss->Clone();
//        Clone_theRTgauss->SetRange(XLeftLim,XRightLim);
//        Clone_theRTgauss->SetFillColor(kBlue);
//        Clone_theRTgauss->SetFillStyle(3013);
//        Clone_theRTgauss->SetLineWidth(1.);
//        Clone_theRTgauss->Draw("SAME FC");
//        ccdf_rt->Update();
//        //
//        ccdf_rt_zoom->cd(iSigma);
//        Func_theRTgauss->SetRange(5.0,5.6);
//        Func_theRTgauss->SetMaximum(0.3);
//        Func_theRTgauss->Draw();
//        Clone_theRTgauss =  (TF1*) Func_theRTgauss->Clone();
//        Clone_theRTgauss->SetRange(XLeftLim,XRightLim);
//        Clone_theRTgauss->Draw("SAME FC");
//        ccdf_rt_zoom->Update();
//     }
//      ccdf_signal  ->Print(Form("cdf_signal_Q2Bin_%d.png",Q2Bin));
//      ccdf_rt->Print(Form("cdf_rt_Q2Bin_%d.png",Q2Bin));
//      ccdf_signal_zoom  ->Print(Form("cdf_signal_zoom_Q2Bin_%d.png",Q2Bin));
//      ccdf_rt_zoom->Print(Form("cdf_rt_zoom_Q2Bin_%d.png",Q2Bin));
//     tagged_mass->setRange("testRange",XLeftLim,XRightLim);
//     RooAbsReal*  testRTRange= theRTgauss->createIntegral(*tagged_mass,NormSet(*tagged_mass),Range("testRange"));
//     std::cout<<Form(" %f%% in the Range %f<mass<%f",testRTRange->getVal(),XLeftLim,XRightLim)<<std::endl;
//      TCanvas* ccdf = new TCanvas("ccdf","cdf",200,10,750,800);
//      ccdf->cd();
//      Func_signalCdf->Draw();
//      ccdf->Print("cdf.png");
     


//      tagged_mass->setRange("3sigmaintegral",mean_rt->getVal()-3*B0SigmaTemp,mean_rt->getVal()+3*B0SigmaTemp);
//      RooAbsReal* BckgInt3Sigma = bkg_exp->createIntegral(*tagged_mass,*tagged_mass,"3sigmaintegral");
//      RooAbsReal* SignInt3Sigma = c_signalFunction.createIntegral(*tagged_mass,*tagged_mass,"3sigmaintegral");
//      NBckgInt3Sigma = BckgInt3Sigma->getVal()*nbkg.getVal();
//      NSignInt3Sigma = SignInt3Sigma->getVal()*nsig.getVal();
//      std::cout<<Form("Bckg event +/- 3 sigma_w from mean RT = %f",NBckgInt3Sigma) <<std::endl;
//      std::cout<<Form("Sign event +/- 3 sigma_w from mean RT = %f",SignInt3Sigma) <<std::endl;
//      std::cout<<Form("Sign %% in +/- 3 sigma_w from mean RT = %f",SignInt3Sigma->getVal()) <<std::endl;
     integTails= TMath::Erfc(sqrt(2));
//      XLeftLim  = Func_theRTgaussCdf->GetX(integTails/2.);
//      XRightLim = Func_theRTgaussCdf->GetX(1-integTails/2.);
     XLeftSet  = Func_signalCdf->GetX(integTails/2.);
     XRightSet = Func_signalCdf->GetX(1-integTails/2.);
     tagged_mass->setRange("2sigmaRTintegral",XLeftSet,XRightSet);
     RooAbsReal* BckgInt2SigmaRT = bkg_exp->createIntegral(*tagged_mass,*tagged_mass,"2sigmaRTintegral");
     RooAbsReal* SignInt2SigmaRT = c_signalFunction.createIntegral(*tagged_mass,*tagged_mass,"2sigmaRTintegral");
     NBckgInt2Sigma = BckgInt2SigmaRT->getVal()*nbkg.getVal();
     NSignInt2Sigma = SignInt2SigmaRT->getVal()*nsig.getVal();
     std::cout<<Form("Bckg event %f<mass<%f from mean RT = %f",XLeftSet,XRightSet,NBckgInt2Sigma) <<std::endl;
     std::cout<<Form("Sign event %f<mass<%f from mean RT = %f",XLeftSet,XRightSet,NSignInt2Sigma) <<std::endl;
     std::cout<<Form("Sign %% in %f<mass<%f from mean RT = %f",XLeftSet,XRightSet,SignInt2SigmaRT->getVal()) <<std::endl;
//==============================================================
//===
//===   MODIFICA della scelta delle SideBand
//===
//==============================================================
//
     if(SigmaProbSign==0){
      std::cout<<Form("FitMassSpectrumRoofit::WARNING: Limits from sigma gauss model")<<std::endl;
      XMinSBL = B0Mass - NSigma1L*B0SigmaTemp;
      XMaxSBL = B0Mass - NSigma2L*B0SigmaTemp;
//      
      XMinSBR = B0Mass + NSigma1R*B0SigmaTemp;
      XMaxSBR = B0Mass + NSigma2R*B0SigmaTemp;
//
     }else if(SigmaProbSign==-1){
      std::cout<<Form("FitMassSpectrumRoofit::WARNING: setting bare limits")<<std::endl;
      XMinSBL = NSigma1L;
      XMaxSBL = NSigma2L;
//      
      XMinSBR = NSigma1R;
      XMaxSBR = NSigma2R;
     }else if(SigmaProbSign==1){
      std::cout<<Form("FitMassSpectrumRoofit::WARNING: Limits from gauss sign prob")<<std::endl;
      integTails= TMath::Erfc(NSigma2L/sqrt(2));
      XMaxSBL	= Func_signalCdf->GetX(integTails/2.);//  Erfc(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between x and infinity; t=x/sqrt(2)
//
      std::cout<<Form("Sideband Left  Internal Limit&Prob [%f - %f]",XMaxSBL, integTails/2.)<<std::endl;
      
//       integTails= TMath::Erfc(fabs(NSigma2L-NSigma1L)/sqrt(2));
//       XMinSBL	= XMaxSBL-(fabs(mean_rt->getVal()-Func_signalCdf->GetX(integTails/2.)));
      XMinSBL   = NSigma1L;
      integTails= TMath::Erfc(NSigma1R/sqrt(2));
      XMinSBR	= Func_signalCdf->GetX(1-integTails/2.);
      std::cout<<Form("Sideband Right Internal Limit&Prob [%f - %f]",XMinSBR, integTails/2.)<<std::endl;
//       integTails= TMath::Erfc(fabs(NSigma2R-NSigma1R)/sqrt(2));
//       XMaxSBR	= XMinSBR+(fabs(mean_rt->getVal()-Func_signalCdf->GetX(1-integTails/2.)));
      XMaxSBR   = NSigma2R;
     }else{
      std::cout<<Form("Sideband Definition=>>SigmaProbSign: 	INVALID OPTION: %d !!! Exit...",SigmaProbSign)<<std::endl;
      exit(1);
     } 
//     
     std::cout<<Form("Sideband Left  [%f-%f]",XMinSBL, XMaxSBL)<<std::endl;
     std::cout<<Form("Sideband Right [%f-%f]",XMinSBR, XMaxSBR)<<std::endl;

     ccdf_signal ->cd();
     Func_signal->SetRange(XMinSign,XMaxSign);
//     Func_signal->SetNormalized(true);
     TF1* Clone_signal  =  (TF1*) Func_signal->Clone();
     TF1* Clone_signalL =  (TF1*) Func_signal->Clone();
     TF1* Clone_signalR =  (TF1*) Func_signal->Clone();
//     Clone_signal->SetNormalized(true);
//     double Func_signal_Integ = Clone_signal->Integral(5.0,5.6);
//     std::cout<<"Func_signal_Integ ="<<Func_signal_Integ<<std::endl;
//     Clone_signal->GetHistogram()->Scale(1/Func_signal_Integ);
     Clone_signal->SetMaximum(0.1);
     Clone_signal->Draw();
//     Clone_signalL->SetNormalized(true);
//     double Func_signal_IntegL = Clone_signal->Integral(XMinSBL,XMaxSBL);
     Clone_signalL->SetRange(XMinSBL,XMaxSBL);
//     std::cout<<"Func_signal_IntegL ="<<Func_signal_IntegL<<std::endl;
//     Clone_signalL->GetHistogram()->Scale(Func_signal_IntegL/Func_signal_Integ);
     Clone_signalL->SetFillColor(kBlue);
     Clone_signalL->SetFillStyle(3013);
     Clone_signalL->SetLineWidth(1.);
     Clone_signalL->Draw("SAME FC");
//     Clone_signalR->SetNormalized(true);
//     double Func_signal_IntegR = Clone_signal->Integral(XMinSBR,XMaxSBR);
     Clone_signalR->SetRange(XMinSBR,XMaxSBR);
//     std::cout<<"Func_signal_IntegR ="<<Func_signal_IntegR<<std::endl;
//     Clone_signalR->GetHistogram()->Scale(Func_signal_IntegR/Func_signal_Integ);
     Clone_signalR->SetFillColor(kBlue);
     Clone_signalR->SetFillStyle(3013);
     Clone_signalR->SetLineWidth(1.);
     Clone_signalR->Draw("SAME FC");
     ccdf_signal ->Print(Form("cdf_signal_Q2Bin_%d.png",Q2Bin));
//     std::cout<<Form("Func_signal_IntegL%=%f Func_signal_IntegR%=%f",Func_signal_IntegL/Func_signal_Integ,Func_signal_IntegR/Func_signal_Integ)<<std::endl;
//exit(0);
//     RooRealVar* tagged_massS = new RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", 5.0, 5.6, "GeV");
//     RooRealVar* tagged_massF = new RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", 4.9, 5.6, "GeV");
//     tagged_mass->setRange(XMinSign,XMaxSign);
//     tagged_mass->setRange("full",XMinSign,XMaxSign);
     tagged_mass->setRange( "SBLeft"  ,XMinSBL,XMaxSBL);
     tagged_mass->setRange( "SBRight" ,XMinSBR,XMaxSBR);
//      tagged_massS->setRange("SBLeftS" ,XMinSBL,XMaxSBL);
//      tagged_massS->setRange("SBRightS",XMinSBR,XMaxSBR);
//      tagged_massF->setRange("SBLeftF" ,XMinSBL,XMaxSBL);
//      tagged_massF->setRange("SBRightF",XMinSBR,XMaxSBR);
     
     double tagged_mass_rangeFitMin=tagged_mass->getMin("full");
     double tagged_mass_rangeFitMax=tagged_mass->getMax("full");
     tagged_mass_rangeValMin=tagged_mass->getMin();
     tagged_mass_rangeValMax=tagged_mass->getMax();
     std::cout<<Form("%f<[tagged_mass fit mass range]<%f",tagged_mass_rangeFitMin,tagged_mass_rangeFitMax)<<std::endl; 
     std::cout<<Form("%f<[tagged_mass val mass range]<%f",tagged_mass_rangeValMin,tagged_mass_rangeValMax)<<std::endl; 
     
//     RooExponential* bkgEXP = new RooExponential(bkg_exp);
  
//      RooAbsReal* BckgAll = bkg_exp->createIntegral(*tagged_mass);
//      RooAbsReal* BckgFull = bkg_exp->createIntegral(*tagged_massF);
//      RooAbsReal* BckgFits = bkg_exp->createIntegral(*tagged_massS);
//      std::cout<<Form("Integrals All=%f Fits=%f and Full=%f",BckgAll->getVal(),BckgFits->getVal(),BckgFull->getVal())<<std::endl; 
     
     RooAbsReal* BckgSBL = bkg_exp->createIntegral(*tagged_mass,*tagged_mass,"SBLeft");
//     RooAbsReal* BckgSBL = bkg_exp->createIntegral(*tagged_mass,NormSet(*tagged_mass),Range("SBLeft"));
     RooAbsReal* BckgSBR = bkg_exp->createIntegral(*tagged_mass,*tagged_mass,"SBRight");
     double BckgEventsSBL = BckgSBL->getVal()*nbkg.getVal();
     double BckgEventsSBR = BckgSBR->getVal()*nbkg.getVal();
//     

     RooAbsReal* SignSBL = c_signalFunction.createIntegral(*tagged_mass,*tagged_mass,"SBLeft");
     RooAbsReal* SignSBR = c_signalFunction.createIntegral(*tagged_mass,*tagged_mass,"SBRight");
     double SignEventsSBL = SignSBL->getVal()*nsig.getVal();
     double SignEventsSBR = SignSBR->getVal()*nsig.getVal();
//     
     RooAbsReal* SignSBL_wt = c_WTgauss->createIntegral(*tagged_mass,*tagged_mass,"SBLeft");
     RooAbsReal* SignSBR_wt = c_WTgauss->createIntegral(*tagged_mass,*tagged_mass,"SBRight");
     double SignEventsSBL_wt = SignSBL_wt->getVal()*nsig.getVal()*(1-fraction->getVal());
     double SignEventsSBR_wt = SignSBR_wt->getVal()*nsig.getVal()*(1-fraction->getVal());
//     
     RooAbsReal* SignSBL_rt = c_RTgauss->createIntegral(*tagged_mass,*tagged_mass,"SBLeft");
     RooAbsReal* SignSBR_rt = c_RTgauss->createIntegral(*tagged_mass,*tagged_mass,"SBRight");
     double SignEventsSBL_rt = SignSBL_rt->getVal()*nsig.getVal()*fraction->getVal();
     double SignEventsSBR_rt = SignSBR_rt->getVal()*nsig.getVal()*fraction->getVal();
//     
     RooAbsReal* ModelSBL = fitFunction.createIntegral(*tagged_mass,*tagged_mass,"SBLeft");
     RooAbsReal* ModelSBR = fitFunction.createIntegral(*tagged_mass,*tagged_mass,"SBRight");
     double ModelEventsSBL = ModelSBL->getVal()*(nsig.getVal()+nbkg.getVal());
     double ModelEventsSBR = ModelSBR->getVal()*(nsig.getVal()+nbkg.getVal());
     double RealEventsSBL  = data->sumEntries(Form("tagged_mass>%f&&tagged_mass<%f",XMinSBL,XMaxSBL)) ;
     double RealEventsSBR  = data->sumEntries(Form("tagged_mass>%f&&tagged_mass<%f",XMinSBR,XMaxSBR)) ;
//
     RooAbsReal* SignalFull  = c_signalFunction.createIntegral(*tagged_mass,*tagged_mass,"full");
//     RooAbsReal* BckgFull    = bkg_exp->createIntegral(*tagged_mass,*tagged_mass,"full");
     RooAbsReal* ModelFull   = fitFunction.createIntegral(*tagged_mass,*tagged_mass,"full");

//     
     std::cout<<Form("Signal Norm Check = %f",SignalFull->getVal())  <<std::endl;
     std::cout<<Form("Estimated Bckg	events SB Left = %f",BckgEventsSBL)  <<std::endl;
     std::cout<<Form("Estimated Bckg	events SB Right= %f",BckgEventsSBR)  <<std::endl;
     std::cout<<Form("Estimated Bckg	events SB Total= %f",BckgEventsSBL+BckgEventsSBR)  <<std::endl;
     std::cout<<Form("========> real	events SB Left = %f",RealEventsSBL ) <<std::endl;
     std::cout<<Form("========> real	events SB Right= %f",RealEventsSBR ) <<std::endl;
     std::cout<<Form("========> real	events SB Total= %f",RealEventsSBL+RealEventsSBR ) <<std::endl;
     std::cout<<Form("========> real	all events inside [full] range Total= %f",data->sumEntries(Form("tagged_mass>%f&&tagged_mass<%f",tagged_mass_rangeMin,tagged_mass_rangeMax)) ) <<std::endl;
     std::cout<<Form("========> Fit	all events inside [full] range Total= %f",ModelFull->getVal()*(nsig.getVal()+nbkg.getVal())) <<std::endl;
     std::cout<<Form("========> real	all events Total= %f",data->sumEntries(Form("tagged_mass>%f&&tagged_mass<%f",XMinSign,XMaxSign)) ) <<std::endl;
     std::cout<<Form("========> Fit	all events Total= %f",(nsig.getVal()+nbkg.getVal())) <<std::endl;
//
     std::cout<<Form("Estimated Sign	events SB Left = %f",SignEventsSBL)  <<std::endl;
     std::cout<<Form("Estimated Sign	events SB Right= %f",SignEventsSBR)  <<std::endl;
     std::cout<<Form("Estimated Sign WT events SB Left = %f",SignEventsSBL_wt)  <<std::endl;
     std::cout<<Form("Estimated Sign WT	events SB Right= %f",SignEventsSBR_wt)  <<std::endl;
     std::cout<<Form("Estimated Sign RT events SB Left = %f",SignEventsSBL_rt)  <<std::endl;
     std::cout<<Form("Estimated Sign RT	events SB Right= %f",SignEventsSBR_rt)  <<std::endl;
     std::cout<<Form("Estimated Model	events SB Left = %f",ModelEventsSBL) <<std::endl;
     std::cout<<Form("Estimated Model	events SB Right= %f",ModelEventsSBR) <<std::endl;
     std::cout<<Form("events signal in SB [all]/ signal Tot  = %f",(SignEventsSBL+SignEventsSBR)/ nsig.getVal())<<std::endl;
     std::cout<<Form("events signal in SB [RT] / signal Tot  = %f",(SignEventsSBL_rt+SignEventsSBR_rt)/ nsig.getVal())<<std::endl;
     std::cout<<Form("events signal in SB [WT] / signal Tot  = %f",(SignEventsSBL_wt+SignEventsSBR_wt)/ nsig.getVal())<<std::endl;
     std::cout<<Form("events signal in SB [all]/ SB     Tot  = %f",(SignEventsSBL+SignEventsSBR)/ (RealEventsSBL+RealEventsSBR))<<std::endl;
     std::cout<<Form("events signal in SB [RT] / SB     Tot  = %f",(SignEventsSBL_rt+SignEventsSBR_rt)/ (RealEventsSBL+RealEventsSBR))<<std::endl;
     std::cout<<Form("events signal in SB [WT] / SB     Tot  = %f",(SignEventsSBL_wt+SignEventsSBR_wt)/ (RealEventsSBL+RealEventsSBR))<<std::endl;
     std::cout<<Form("events signal in SB Left / signal Tot  = %f",(SignEventsSBL)/ nsig.getVal())<<std::endl;
     std::cout<<Form("events signal in SB Right/ signal Tot  = %f",(SignEventsSBR)/ nsig.getVal())<<std::endl;
     std::cout<<Form("events signal in SB Left / SB     Tot  = %f",(SignEventsSBL)/ (RealEventsSBL+RealEventsSBR))<<std::endl;
     std::cout<<Form("events signal in SB Right/ SB     Tot  = %f",(SignEventsSBR)/ (RealEventsSBL+RealEventsSBR))<<std::endl;
     std::cout<<Form("events signal in SB Left / SB     Left = %f",(SignEventsSBL)/ (RealEventsSBL))<<std::endl;
     std::cout<<Form("events signal in SB Right/ SB     Right= %f",(SignEventsSBR)/ (RealEventsSBR))<<std::endl;
     std::cout<<"\n"<<std::endl;
     std::cout<<Form("%d &  %3.2f-%3.2f & %3.2f-%3.2f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f & %3.2f\\\\ \n",Q2Bin,XMinSBL,XMaxSBL,XMinSBR,XMaxSBR,\
     (SignEventsSBL)/nsig.getVal(),(SignEventsSBR)/nsig.getVal(), (SignEventsSBL)/ (RealEventsSBL),\
     (SignEventsSBR)/ (RealEventsSBR),\
     (SignEventsSBL_rt+SignEventsSBR_rt)/(RealEventsSBL+RealEventsSBR),\
     (SignEventsSBL_wt+SignEventsSBR_wt)/(RealEventsSBL+RealEventsSBR),\
     (SignEventsSBL+SignEventsSBR)/(RealEventsSBL+RealEventsSBR),\
     (RealEventsSBL+RealEventsSBR)/NBckgInt2Sigma) <<std::endl;
//     
     
     double xbinw = pdfHist->GetXaxis()->GetBinWidth(1);
     cout<<"Binw pdfHist ="<<xbinw<<endl;
//      for (int i = 1; pdfHist->GetNbinsX(); ++i) {
//         double xmass = xbinw/2.+(i-1)*xbinw;
// //      const RooArgSet * dataLoad = data->get (i);
// //      double xmass = dataLoad->getRealValue(tagged_mass->GetName());
// 	pdfHist->Fill(xmass, fitFunction.evaluate() );
// 	sigHist->Fill(xmass, c_signalFunction.evaluate());
// 	bkgHist->Fill(xmass, c_frt.evaluate());
//      }
     double NStepMass  = pdfHist->GetNbinsX();
     double NBINFactor = NStepMass/masHist->GetNbinsX();
     fitFunction.fillHistogram(pdfHist,*tagged_mass);
     c_signalFunction.fillHistogram(sigHist,*tagged_mass);
     bkg_exp->fillHistogram(bkgHist,*tagged_mass);
     tagged_mass->setRange(XMinSign,XMaxSign);
     RooAbsReal* BckgFullW   = bkg_exp->createIntegral(*tagged_mass,*tagged_mass);
     RooAbsReal* SignalFullW = c_signalFunction.createIntegral(*tagged_mass,*tagged_mass);
     RooAbsReal* ModelFullW  = fitFunction.createIntegral(*tagged_mass,*tagged_mass);
     
//     tagged_mass->setRange("full1",5.0,5.6);
     RooAbsReal* BckgFullS   = bkg_exp->createIntegral(*tagged_mass,*tagged_mass,"full");
     RooAbsReal* SignalFullS = c_signalFunction.createIntegral(*tagged_mass,*tagged_mass,"full");
     RooAbsReal* ModelFullS  = fitFunction.createIntegral(*tagged_mass,*tagged_mass,"full");
     
     double scaleModelW  = ModelFullW->getVal()/ModelFullS->getVal();
     double scaleSignalW = SignalFullW->getVal()/SignalFullS->getVal();
     double scaleBckgW   = BckgFullW->getVal()/BckgFullS->getVal();

     int MinBinRangeFit = pdfHist->GetXaxis()->FindBin(tagged_mass_rangeMin);
     int MaxBinRangeFit = pdfHist->GetXaxis()->FindBin(tagged_mass_rangeMax);
     std::cout<<Form("pdfHist (=model of mass spectrum) MinBinRangeFit = %d MaxBinRangeFit = %d NumBins = %d",MinBinRangeFit,MaxBinRangeFit,pdfHist->GetNbinsX()) <<std::endl;
//      pdfHist->Scale(NBINFactor*(nsig.getVal()+nbkg.getVal())*scaleModelW);
//      sigHist->Scale(NBINFactor*(nsig.getVal())*scaleSignalW);
//      bkgHist->Scale(NBINFactor*(nbkg.getVal())*scaleBckgW);
     pdfHist->Scale(NBINFactor*(nsig.getVal()+nbkg.getVal()));
     sigHist->Scale(NBINFactor*(nsig.getVal()));
     bkgHist->Scale(NBINFactor*(nbkg.getVal()));
//     cout<< Form("Integ Model (4.9-5.6) =%f  Integ Model (5.0-5.6) =%f scaleModelW=%f", ModelFullW->getVal(),ModelFullS->getVal(),scaleModelW)<<endl;
//     cout<< Form("Integ Signal(4.9-5.6) =%f  Integ Signal(5.0-5.6) =%f scaleSignalW=%f", SignalFullW->getVal(),SignalFullS->getVal(),scaleSignalW)<<endl;
//     cout<< Form("Integ Bckg  (4.9-5.6) =%f  Integ Bckg  (5.0-5.6) =%f scaleBckgW=%f", BckgFullW->getVal(),BckgFullS->getVal(),scaleBckgW)<<endl;
     cout<< Form("Integ Model (%f-%f) =%f  Integ Model (%f-%f) =%f scaleModelW=%f" ,XMinSign,XMaxSign,ModelFullW->getVal() ,tagged_mass_rangeMin,tagged_mass_rangeMax,ModelFullS->getVal() ,scaleModelW)<<endl;
     cout<< Form("Integ Signal(%f-%f) =%f  Integ Signal(%f-%f) =%f scaleSignalW=%f",XMinSign,XMaxSign,SignalFullW->getVal(),tagged_mass_rangeMin,tagged_mass_rangeMax,SignalFullS->getVal(),scaleSignalW)<<endl;
     cout<< Form("Integ Bckg  (%f-%f) =%f  Integ Bckg  (%f-%f) =%f scaleBckgW=%f"  ,XMinSign,XMaxSign,BckgFullW->getVal()  ,tagged_mass_rangeMin,tagged_mass_rangeMax,BckgFullS->getVal(),scaleBckgW)<<endl;
//exit(1);
   //   int MinBinSBL = bkgHist->GetXaxis()->FindBin(mean_rt->getVal()-3*B0SigmaTemp);
//      int MaxBinSBL = bkgHist->GetXaxis()->FindBin(mean_rt->getVal()+3*B0SigmaTemp);
//      double BckgIntSBL= nbkg.getVal()*bkgHist->Integral(MinBinSBL,MaxBinSBL)/bkgHist->Integral(MinBinRangeFit,MaxBinRangeFit);
//      std::cout<<Form("pdfHist (=model of mass spectrum) MinBinSBL = %d MaxBinSBL = %d Integ = %f",MinBinSBL,MaxBinSBL,BckgIntSBL) <<std::endl;
//  
     masHist->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
     masHist->SetMarkerStyle(8);
     masHist->SetMarkerSize(MarkerSizeSet);
     masHist->SetTitle("");
     masHist->Draw("E1");
//     masHist.Draw("p");
     pdfHist->GetXaxis()->SetRangeUser(tagged_mass_rangeMin,tagged_mass_rangeMax);
     sigHist->GetXaxis()->SetRangeUser(tagged_mass_rangeMin,tagged_mass_rangeMax);
     bkgHist->GetXaxis()->SetRangeUser(tagged_mass_rangeMin,tagged_mass_rangeMax);
     pdfHist->SetLineWidth(PlotLineWidth);
     pdfHist->SetFillColor(0);
     pdfHist->SetLineColor(kBlue);
     pdfHist->Draw("same,HIST C");
     sigHist->SetLineWidth(PlotLineWidth);
     sigHist->SetLineColor(kMagenta);
     sigHist->SetLineStyle(kDashed);
     sigHist->SetFillColor(0);
     sigHist->Draw("same,HIST C");
     bkgHist->SetLineWidth(PlotLineWidth);
     bkgHist->SetLineColor(kRed);
     bkgHist->SetLineStyle(kDashed);
     bkgHist->SetFillColor(0);
     bkgHist->Draw("same,HIST C");
     leg_sign->Draw("same");
     
     if(GenMiniMC&&FirstMC){
      ivolte++;
      
      tagged_mass->setRange("FullMass",XMinSign,XMaxSign);
      RooAbsReal* BckgFullMass = bkg_exp->createIntegral(*tagged_mass,*tagged_mass,"FullMass");
      RooAbsReal* SignFullMass = c_signalFunction.createIntegral(*tagged_mass,*tagged_mass,"FullMass");
      int NBckgIntFull = BckgFullMass->getVal()*nbkg.getVal();      
      int NSignIntFull = SignFullMass->getVal()*nsig.getVal();      
      int  NumMassNew = (NSignIntFull+NBckgIntFull)*NFactGen;
      NumMassNewGen = NumMassNew;
      std::cout<<Form("Start to generate %d events Mass...\n",NumMassNew) <<std::endl;
      std::cout<<Form("Bckg generated = %d\n",NBckgIntFull*NFactGen) <<std::endl;
      std::cout<<Form("Sign generated = %d\n",NSignIntFull*NFactGen) <<std::endl;
      RooRandom::randomGenerator()->SetSeed(ivolte*1000+ivolte);
      geneMassNew = fitFunction.generate(*tagged_mass,NumMassNew) ;
     }
  
  fitMassFile->Close();
  return B0SigmaTemp;

}
//=========================================================================================

RooGaussian* _constrainVar(RooRealVar *var,RooWorkspace *w){
    
//    float constr[2] = *_getFittedVar(var.GetName(), w);
    RooRealVar c_val(Form("c_val_%s", var->GetName()),Form("c_val_%s", var->GetName()),var->getVal());
    RooRealVar c_err(Form("c_err_%s", var->GetName()),Form("c_err_%s", var->GetName()),var->getError());
    RooGaussian* gauss_constr =
                            new RooGaussian(   Form("c_%s", var->GetName()) , 
                                Form("c_%s", var->GetName()) , 
                                *var         ,  
                                RooConst( var->getVal() ), 
                                RooConst( var->getError() ) 
                                ) ;
    std::cout<< Form("constraining var %s: %f with uncertainty %f - limits [%f , %f]",var->GetName(),c_val.getVal(),c_err.getVal(),var->getMin(),var->getMax())<<std::endl;  
    double checkMax = var->getVal()+ 7*var->getError();			      
    double checkMin = var->getVal()- 7*var->getError();			      
    std::cout<< Form("Warning in _constrainVar: limits for %s, from [%f,%f] ==> [%f,%f]\n",var->GetName(),var->getMin(),var->getMax(),checkMin,checkMax);
//    var->setMax(checkMax) ;	     
//    var->setMin(checkMin) ;	     
//    std::cout<< Form("Warning: redifine limits var %s: %f with uncertainty %f - limits [%f , %f]",var->GetName(),c_val.getVal(),c_err.getVal(),var->getMin(),var->getMax())<<std::endl;  
                            
    return gauss_constr;
}                           

//=========================================================================================
//
//=========================================================================================

TTree* makeGenSBTTree(RooDataSet* geneDataNew,RooDataSet* geneMassNew,RooBernsteinSideband* BernsteinSideband, TFile*OutFileNtupla ) 
{
   std::cout<<"Start to fill ntupla\n" <<std::endl;
  // Create ROOT TTree filled
  OutFileNtupla->cd();
  TTree* GeneB0TreeOut = new TTree(OutputGeneB0TreeName,OutputGeneB0TreeName) ;
  GeneB0TreeOut -> SetAutoSave(5000000000);
  double  tagged_mass	 ;
  double  cos_theta_l	 ;
  double  cos_theta_k	 ;
  double  phi_kst_mumu   ;
  double  mumuMass	 ;
  double  mumuMassE	 ;
  mumuMass=sqrt(Q2Fake);
  mumuMassE=0.02;
  GeneB0TreeOut->Branch("cos_theta_l"   ,&cos_theta_l    ,   "cos_theta_l/D"   );
  GeneB0TreeOut->Branch("cos_theta_k"   ,&cos_theta_k    ,   "cos_theta_k/D"   );
  GeneB0TreeOut->Branch("phi_kst_mumu"  ,&phi_kst_mumu   ,   "phi_kst_mumu/D"  );
  GeneB0TreeOut->Branch("tagged_mass"   ,&tagged_mass    ,   "tagged_mass/D"   );
  GeneB0TreeOut->Branch("mumuMass"      ,&mumuMass       ,   "mumuMass/D"      );
  GeneB0TreeOut->Branch("mumuMassE"     ,&mumuMassE      ,   "mumuMassE/D"     );
  int NGenSB = geneDataNew->sumEntries();
  int ii=0;
  int MaxrowData=geneDataNew->sumEntries();
  const RooArgSet* rowMass =  0;;       
  const RooArgSet* rowData =  0;; 
//  double Pi =TMath::Pi();
  for (int i=0 ; i<geneMassNew->sumEntries() ; i++) {
//  for (int i=0 ; i<geneDataNew->sumEntries() ; i++) {
   
   rowMass =  new  RooArgSet(*ctL,*ctK,*phi);       
   rowMass = geneMassNew->get(i);
   tagged_mass=rowMass->getRealValue("tagged_mass");
   cos_theta_l = -99.;
   cos_theta_k = -99.;
   phi_kst_mumu = -99.;
//   if( (tagged_mass>XMinSBL&&tagged_mass<XMaxSBL)|| 
//       (tagged_mass>XMinSBR&&tagged_mass<XMaxSBR)){
     if(ii<=MaxrowData){
      rowData =  new  RooArgSet(*ctL,*ctK,*phi);       
      rowData = geneDataNew->get(ii);
      if(rowData){
      cos_theta_l =rowData->getRealValue("ctL");
      cos_theta_k =rowData->getRealValue("ctK");
      phi_kst_mumu=rowData->getRealValue("phi");
//      if(cos_theta_l<-1.0||cos_theta_l>1.0) std::cout<<Form("%d cos_theta_l = %f",ii,cos_theta_l) <<std::endl;
//      if(cos_theta_k<-1.0||cos_theta_k>1.0) std::cout<<Form("%d cos_theta_l = %f",ii,cos_theta_k) <<std::endl;
//      if(phi_kst_mumu<-Pi||phi_kst_mumu>Pi) std::cout<<Form("%d phi_kst_mumu = %f",ii,phi_kst_mumu) <<std::endl;
      ii++;
      }
//     }
   }  
   GeneB0TreeOut->Fill() ;
  
  }
  
  if(ii!=NGenSB){std::cout<<Form("Fill %d, but generated sb events  %d",ii,NGenSB) <<std::endl;};
//  GeneB0TreeOut->Write();
//  OutFileNtupla->Close();
  std::cout<<"Closing ntupla\n" <<std::endl;
  return GeneB0TreeOut ;
}
//=========================================================================================
//
//=========================================================================================
//
// Namelist Routine
//
//=========================================================================================
//
//=========================================================================================
std::map<std::string, std::string> ReadNamelist(int argc, char** argv){
   if ( argc>=1 && (strcmp(argv[0],"namelist")>=0) ){
     std::cout<<"Defined namelist: "<<argv[0]<<std::endl;
   }else{
     std::cout<<"Namelist:"<<argv[0]<<"  should be named/renamed namelist*.list "<<argc<<std::endl;
     exit(1);
   }
   std::vector<std::string> split( char *str, char c = ' ');
   ifstream indata;
   std::map<std::string, std::string> mappa;
   std::string line;
   std::vector<std::string>vstring ;
//
    indata.open(argv[0]);
   if(!indata) { // file couldn't be opened
   	std::cout <<"Line: "<<__LINE__ <<" "<<argv[0]<< " Error: fileList can not be opened" << std::endl;
   	exit(1);
   }
   while(std::getline(indata, line)) {
	 line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
	 line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
	 line.erase(std::remove(line.begin(), line.end(), ' ' ), line.end());

 	 char *cstr = new char [line.size()+1];


 	 strcpy (cstr, line.c_str());
//	 cout <<"stringa->"<< cstr << endl;
	 vstring = split(cstr,'=');
	 mappa.insert( std::pair<string,string>(vstring[0],vstring[1]) );
    }
    std::cout<<"//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
    for (map<string,string>::iterator imap = mappa.begin();
    			       imap != mappa.end();
    			       ++imap)
    {
   	std::cout <<"mappa->"<< (*imap).first<<" = "<<(*imap).second << std::endl;
    }
    std::cout<<"//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
    indata.close();	
  return mappa ;
}
std::vector<std::string> split( char *str, char c = ' ')
{
    std::vector<std::string> result;

    while(1)
    {
         char *begin = str;

        while(*str != c && *str)
                str++;

        result.push_back(string(begin, str));

        if(0 == *str++)
                break;
    }

    return result;
}

std::string grep(std::string& filename, std::string keyword)
{
    std::string out("GREP: NOT FOUND!!!");
    std::string line;
    std::size_t found;
    ifstream in(filename.c_str());
//    std::cout<<"keyword = "<<keyword<<"\n"<<std::endl;
//    std::cout<<"filename= "<<filename<<"\n"<<std::endl;
    if( in.is_open())
    {
	  while( getline(in,line) )
	  {
	      found=line.find(keyword);
//	           std::cout<<"line "<<line<<"\n"<<std::endl;
	      if( found!=std::string::npos)
		   return line ;
	  }
    }else{
     std::cout<<"GREP: can't open filename= "<<filename<<"\n"<<std::endl;
    }
    return out;
}
//===============================================================================================================
void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    if(from.empty())
        return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

//===============================================================================================================
void replaceChar(char * txt,const  char * txt1,const  char * txt2) {

  std::stringstream sss,sss1,sss2;
  sss<<txt;
  sss1<<txt1;
  sss2<<txt2;
  std::string ss=sss.str();
  replaceAll( ss,  sss1.str(), sss2.str());
  strcpy(txt,ss.c_str());
  sss.str("");
  sss.clear();
  sss1.str("");
  sss1.clear();
  sss2.str("");
  sss2.clear();
  printf ("replaceChar output=>%s\n",txt);
}  

