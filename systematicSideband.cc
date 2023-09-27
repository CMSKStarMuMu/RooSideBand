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
#include <TMatrixD.h>
#include <TBinomialEfficiencyFitter.h>
#include <TKDTreeBinning.h>
#include <TH2Poly.h>
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
#include <RooArgList.h>
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
#include <RooSimultaneous.h>
#include <RooMultiVarGaussian.h>
#include "RooBernsteinSideband.h"
#include "RooMCStudy.h"
#include "TRatioPlot.h"
//#include "PdfSigAngMass.h"
//#include "ShapeSigAng.h"
#include "RooRealMPFE.h"

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
void ReadSimultaneousModel();
bool CheckLogFile(std::string filename, std::string keyword);
std::string LogSaveDirName =  "/gwpool/users/dini/FittoneNew3/UML-fit-unblinding-changes2-sbsyst/ ";
std::string OutSaveDirName =  "/gwpool/users/dini/FittoneNew3/UML-fit-unblinding-changes2-sbsyst/simFitResults4d/";
std::string OutSaveFileName="";
std::string LogSaveFileName="";
//char OutSaveFileName[400] =  "simFitResult_data_fullAngularMass_Swave_201620172018_MCStat_b4_fullStat_noMeanCon.root";
char *ReadOutSaveFileName;

char PNGNameReadSB3D[400] = "systematicSideband";
int    Q2Bin     =  -99;
int    xCosLHBin =  25;
int    xCosKHBin =  25;
int    xPhiHBin  =  25;
double XMinCosThetaL	     = -1.;
double XMaxCosThetaL	     =  1.;
double XMinCosThetaK	     = -1.;
double XMaxCosThetaK	     =  1.;
double XMinPhi  	     =-TMath::Pi();
double XMaxPhi  	     = TMath::Pi();
double XMinSign = 5.0;
double XMaxSign = 5.6;
double XStepSign = 0.025;
float  xMassHBin = (XMaxSign -XMinSign)/XStepSign;

int MaxSystematic = 100;

TFile* OutSaveFile = 0;
//=================================================
int main (int argc, char** argv) {
  TStopwatch TimeWatch;
  TimeWatch.Start();
  
  if (argc>1 ){
   switch ( *argv[1] ) {

    case '0' :
     Q2Bin = 0;
      break;
    case '1' :
     Q2Bin = 1;
      break;
    case '2' :
     Q2Bin = 2;
      break;
    case '3' :
     Q2Bin = 3;
      break;
    case '4' :
     Q2Bin = 4;
      break;
    case '5' :
     Q2Bin = 5;
      break;
    case '6' :
     Q2Bin = 6;
      break;
    case '7' :
     Q2Bin = 7;
      break;
    case '8' :
     Q2Bin = 8;
      break;

    default :
       ReadOutSaveFileName  = (char *) malloc(strlen(argv[1])+1);
       strcpy(ReadOutSaveFileName,argv[1]);
       OutSaveFileName = Form("%s",ReadOutSaveFileName);
       std::cout<<"========================================================================="<<endl;
       std::cout<<Form("Setting the File: %s",argv[1])<<std::endl;
       std::cout<<"========================================================================="<<endl;

   }
   if(Q2Bin!=-99){//default Bin
    OutSaveFileName=  Form("simFitResult_data_fullAngularMass_Swave_201620172018_b%d-XGBv8_unbl4.root",Q2Bin);
    LogSaveFileName=  Form("Bin%d-201620172018-GenSB",Q2Bin);
   } 

  }else{
     OutSaveFileName=  Form("simFitResult_data_fullAngularMass_Swave_201620172018_b%d-XGBv8_unbl4.root",Q2Bin);
     LogSaveFileName=  Form("Bin%d-201620172018-GenSB",Q2Bin);
     std::cout<<"========================================================================="<<endl;
     std::cout<<Form("Setting DEFAULT File: %s",OutSaveFileName.c_str())<<std::endl;
     std::cout<<Form("Systematic  for Bin : %d",Q2Bin)<<std::endl;
     std::cout<<"========================================================================="<<endl;
  }
  //sprintf(PNGNameReadSB3D,"%s-%s.png",PNGNameReadSB3D,OutSaveFileName);
  

  ReadSimultaneousModel(); 
//  app.Run() ;
  cout<<"End..." <<endl;
  TimeWatch.Stop();
  TimeWatch.Print();
  
  return 0 ;
}
void  ReadSimultaneousModel(){
    RooRealVar* ctL = new RooRealVar("ctL", "ctL",  XMinCosThetaK,XMaxCosThetaK);
    RooRealVar* ctK = new RooRealVar("ctK", "ctK",  XMinCosThetaL,XMaxCosThetaL);
    RooRealVar* phi = new RooRealVar("phi", "phi",  XMinPhi,XMaxPhi);
    RooArgList observables (*ctK, *ctL, *phi);
//    RooRealVar* tagged_mass = new RooRealVar("tagged_mass" , "#mu^{+}#mu^{-}K#pi mass", XMinSign,XMaxSign, "GeV/c^2");
    
    std::vector <TH1D*> HCoeff;
    RooWorkspace* w = 0;
    RooArgSet *params =0; RooRealVar* ivar=0;
    TFile* OutSaveSystematic = 0;
    std::cout<<Form("Opening ref file: %s/%s.0",OutSaveDirName.c_str(),OutSaveFileName.c_str())<<std::endl;
    OutSaveFile = new TFile( Form("%s/%s.0",OutSaveDirName.c_str(),OutSaveFileName.c_str()), "READ" );
    if ( !OutSaveFile || !OutSaveFile ->IsOpen() ) {
      cout<<Form("File not found: %s\n",OutSaveFileName.c_str())<<endl;
      exit(1);
    }else{
     OutSaveSystematic = new TFile( Form("systematic-%s",OutSaveFileName.c_str()), "RECREATE" );
     OutSaveFile->cd();
    } 
    w = (RooWorkspace*)OutSaveFile->Get("wsp_out");
    if ( !w || w->IsZombie() ) {
     cout<<Form("Workspace not found in file:%s\n",OutSaveFile->GetName())<<endl;
     exit(1);
    } else {
     cout<<Form("Workspace Found!!! In file : %s\n",OutSaveFile->GetName())<<endl;
     w->Print();
    }
//    simPdf = (RooSimultaneous*)w->pdf("simPdf");
//    params	= (RooArgSet *)simPdf->getParameters(observables);
   RooArgSet  params1	= (RooArgSet )w->allVars();
   params=&params1;
   
    params1.Print();
    if(!params){
     cout<<Form("Params not found in file:%s\n",OutSaveFile->GetName())<<endl;
     exit(1);
    }
    vector<double> ValRef;
    vector<double> ErrRef;
    auto iter  = params->createIterator();
    ivar =  (RooRealVar*)iter->Next();
    int ii=0;
     while (ivar) {
       double xmin=ivar->getVal()-7*ivar->getError();
       double xmax=ivar->getVal()+7*ivar->getError();
//        if(strncmp (ivar->GetName(),"Fl",2) == 0){
//         xmin=0.;
//         xmax=1;
//        }
//        if(strncmp (ivar->GetName(),"P",1) == 0){
//         xmin=-2.2;
//         xmax=2.2;
//        }
       HCoeff.emplace_back(new TH1D(Form("HPar_%s",ivar->GetName()),ivar->GetName(),200,xmin,xmax)); 
       ValRef.push_back(ivar->getVal());
       ErrRef.push_back(ivar->getError());
       cout<<Form("create hist %s\n",HCoeff[ii]->GetName());
       cout<<Form(" %s = %f\n",ivar->GetName(),ivar->getVal());
       ivar =  (RooRealVar*)iter->Next();
       ii++;
     }
    int icount =0; 
    for(int i=0;i<MaxSystematic;i++){
//       std::string tmp = Form("%s/%s%d.log",LogSaveDirName.c_str(),LogSaveFileName.c_str(),i);
       bool checkLOg = CheckLogFile(Form("%s/%s%d.log",LogSaveDirName.c_str(),LogSaveFileName.c_str(),i),"Not converge");
       if(checkLOg) continue;
       std::cout<<Form("Opening file: %s/%s.%d",OutSaveDirName.c_str(),OutSaveFileName.c_str(),i)<<std::endl;
       OutSaveFile = new TFile( Form("%s/%s.%d",OutSaveDirName.c_str(),OutSaveFileName.c_str(),i), "READ" );
       if ( !OutSaveFile || !OutSaveFile ->IsOpen() ) {
         cout<<Form("File not found: %s... skip!!!\n",OutSaveFile->GetName())<<endl;
         continue;
//         exit(1);
       }
       w = (RooWorkspace*)OutSaveFile->Get("wsp_out");
       if ( !w || w->IsZombie() ) {
        cout<<Form("Workspace not found in file:%s\n",OutSaveFile->GetName())<<endl;
        cout<<Form("Skip file:  file:%s\n",OutSaveFile->GetName())<<endl;
        continue;
//        exit(1);
       } else {
        cout<<Form("Workspace Found!!! In file : %s\n",OutSaveFile->GetName())<<endl;
       }
//        simPdf = (RooSimultaneous*)w->pdf("simPdf");
//        params      = (RooArgSet *)simPdf->getParameters(observables);
     
       RooArgSet  params1	= (RooArgSet )w->allVars();
       params=&params1;
   
       if(!params){
        cout<<Form("Params not found in file:%s\n",OutSaveFile->GetName())<<endl;
        exit(1);
       }

//        if(first){
//        auto iter  = params->createIterator();
//        ivar =  (RooRealVar*)iter->Next();
//        int ii=0;
//         while (ivar) {
// //	  HTemp = new TH1D(Form("HPar_%d",ii),ivar->GetName(),100,ivar->getMin(),ivar->getMax());
// 	  HCoeff.emplace_back(new TH1D(Form("HPar_%d",ii),ivar->GetName(),100,ivar->getMin(),ivar->getMax())); 
// 	  cout<<Form("create hist %s\n",HCoeff[ii]->GetName());
//           ivar =  (RooRealVar*)iter->Next();
// 	  ii++;
// 	}
//        }
       icount++;
       int ii=0;
       auto iter  = params->createIterator();
       ivar =  (RooRealVar*)iter->Next();
       while (ivar) {
	  
//	  HTemp = new TH1D(Form("HPar_%d",ii),"htaemp",100,-10,10);
// 	  HTemp = new TH1D(Form("HPar_%s",ivar->GetName()),ivar->GetName(),100,ivar->getMin(),ivar->getMax());
// 	  cout<<Form("create hist %s\n",Form("HPar_%d",ii));
// 	  HCoeff.push_back(HTemp); 
         cout <<Form("Loop %d  %s=%f +/- %f",i,ivar->GetName(),ivar ->getVal(),ivar ->getError())<<endl;
	 HCoeff[ii]->Fill(ivar ->getVal());
//	 HCoeff[ii]->Fill(ivar ->getVal());
//	 if(HCoeff[ii]) HCoeff[ii]->Fill(ivar ->getVal());
         ivar =  (RooRealVar*)iter->Next();
	 ii++;
       }
       OutSaveFile->Close();
     }
    OutSaveSystematic->cd();
    int i=0;
    int j=0;
    for (std::vector <TH1D*>::iterator it = HCoeff.begin(); it != HCoeff.end(); ++it) {
      if(strncmp (((*it)->GetTitle()),"P",1) == 0 || strncmp (((*it)->GetTitle()),"Fl",2) == 0){
       j++;
       std::string P  = ((*it)->GetTitle());
       if(j<5){
        cout<<Form("$%s_%s$ & $ %5.4f $ \\",P.substr(0, 1).c_str(),P.substr(1, 1).c_str(),(*it)->GetStdDev())<<endl;
       }else{
        cout<<Form("$%s'_%s$ & $ %5.4f $ \\",P.substr(0, 1).c_str(),P.substr(1, 1).c_str(),(*it)->GetStdDev())<<endl;
       }
       
      }
      // $P_1$  & $ -0.0182$ & $ -0.0496 $  & 0.0099
       
     (*it)->Write();
     cout<<Form("%s StdDev=%f [Mean=%f   ValRef=%f DiffVal=%f]",(*it)->GetTitle(),(*it)->GetStdDev(),(*it)->GetMean(),ValRef[i],(*it)->GetMean()-ValRef[i])<<endl;
     i++;
    }  
    cout<<Form("Sample used in systematic %d",icount)<<endl;
    OutSaveSystematic->Close();  
} 
//========================================================
// bool CheckLogFile(std::string LogFileName, std::string substring){
// 
//    
//   std::ifstream file_text(LogFileName);  // creating file_text object of ifstream type.
//   
//   std::string x_line;   
//   
//   bool ans=false;   // declaring ans variable which will get information if substring present or not.
//   
//   int line=1;   // it will count line and give an answer on which line we get a substring.
//   
//   if (file_text.is_open())       //is_open open the text file.
//   {
//     while ( getline (file_text,x_line) )    
//     {
//       if (x_line.find(substring, 0) != std::string::npos) {
//         {
//           
//           std::cout<<"CheckLogFile:: Substring found at line "<<line<<" in file: "<<LogFileName.c_str()<<std::endl;
//           ans=true;     // if substring present make ans variable true.
//     }
//     
//     line++;
//     }
//     file_text.close(); //to close text file.
//   }
// }
//   else 
//   std::cout << Form("CheckLogFile:: Unable to open file: %s \n",LogFileName.c_str())<<std::endl; 
//   
//   if(!ans)   // if subtring not present.
//   std::cout<<Form("CheckLogFile:: Substring \"%s\" not found in file: %s \n",substring.c_str(),LogFileName.c_str())<<std::endl;
//   return 0; 
// }
//===============================================================================================================
bool CheckLogFile(std::string filename, std::string keyword)
{
    bool out= false;
    std::string line;
    std::size_t found;
    ifstream in(filename.c_str());
//    std::cout<<"keyword = "<<keyword<<"\n"<<std::endl;
//    std::cout<<"filename= "<<filename<<"\n"<<std::endl;
    if( in.is_open())
    {
	  while( getline(in,line) )
	  {
	      found=line.find(keyword.c_str());
	      if( found!=std::string::npos){
		   out=true ;
                   std::cout<<Form("CheckLogFile: keyword=%s found in filename=%s",keyword.c_str(),filename.c_str())<<"\n"<<std::endl;
		   break;
	      }   
	  }
    }else{
     std::cout<<"CheckLogFile: can't open filename= "<<filename<<"\n"<<std::endl;
    }
    in.close();
    return out;
}
