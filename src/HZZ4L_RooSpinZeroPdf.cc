#include "Riostream.h" 
#include "HZZ4L_RooSpinZeroPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "TH2F.h"

using namespace TMath;

ClassImp(HZZ4L_RooSpinZeroPdf) 

  HZZ4L_RooSpinZeroPdf::HZZ4L_RooSpinZeroPdf(const char *name, const char *title, 
					     RooAbsReal& _kd,
					     RooAbsReal& _kdint,
					     RooAbsReal& _fai,
					     vector<TH2F*>& _histos, 
					     bool      _withAcc):
   RooAbsPdf(name,title), 
   kd("kd","kd",this,_kd),
   kdint("kdint","kdint",this,_kdint),
   fai("fai","fai",this,_fai),
   histos(_histos),
   withAcc(_withAcc)
 { 
  if (histos.size()!=3){
    coutE(InputArguments) << "HZZ4L_RooSpinZeroPdf::HZZ4L_RooSpinZeroPdf(" << GetName() 
			  << ") number of histograms must be 3" << endl ;
    assert(0);
  }
 } 


 HZZ4L_RooSpinZeroPdf::HZZ4L_RooSpinZeroPdf(const HZZ4L_RooSpinZeroPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   kd("kd",this,other.kd),
   kdint("kdint",this,other.kdint),
   fai("fai",this,other.fa3),
   histos(other.histos),
   withAcc(other.withAcc)
 { 
 } 



 Double_t HZZ4L_RooSpinZeroPdf::evaluate() const 
 { 

   double value = 0.;

   int binx =  histos[0]->GetXaxis()->FindBin(kd);
   int biny =  histos[0]->GetYaxis()->FindBin(kdint);
   
   Double_t T1 = histos[0]->GetBinContent(binx, biny);
   Double_t T2 = histos[1]->GetBinContent(binx, biny);
   Double_t T4 = histos[2]->GetBinContent(binx, biny);

   double mysgn = 1;

   if(fai < 0.) 
     {
       mysgn = -1.;
     }
   
   value = (1.-fabs(fai)) * T1 + fabs(fai) * T2 + mysgn*sqrt((1.-fabs(fai))*fabs(fai)) * T4; 
   
   if ( value <= 0.) return 1.0e-200;
   
   return value ; 
   
 } 

Int_t HZZ4L_RooSpinZeroPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if (matchArgs(allVars,analVars,RooArgSet(*kd.absArg(), *kdint.absArg()))) return 3 ;
  if (matchArgs(allVars,analVars,kd)) return 1 ;
  if (matchArgs(allVars,analVars,kdint)) return 2 ;

  return 0 ;

}

Double_t HZZ4L_RooSpinZeroPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{

  int nbinsx = histos[0]->GetXaxis()->GetNbins();
  double binwidthx = histos[0]->GetXaxis()->GetBinWidth(1);
  double xMin = histos[0]->GetXaxis()->GetBinLowEdge(1);
  double xMax = histos[0]->GetXaxis()->GetBinUpEdge(nbinsx);
  double dx = (xMax - xMin) / nbinsx; 

  
  int nbinsy = histos[0]->GetYaxis()->GetNbins();
  double binwidthy = histos[0]->GetYaxis()->GetBinWidth(1);
  double yMin = histos[0]->GetYaxis()->GetBinLowEdge(1);
  double yMax = histos[0]->GetYaxis()->GetBinUpEdge(nbinsy);
  double dy = (yMax - yMin) / nbinsy; 

  /*
  std::cout << "nbinsx = " << nbinsx << ", binwidthx = " << binwidthx << ", xMin = " << xMin << ", xMax = " << xMax << "\n";
  std::cout << "nbinsy = " << nbinsy << ", binwidthy = " << binwidthy << ", yMin = " << yMin << ", yMax = " << yMax << "\n";
  */
  
   switch(code)
     {

       // integrate out kd, depend on kdint
     case 1: 
       {

	 int biny = histos[0]->GetYaxis()->FindBin(kdint);
	 
	 double Int_T1 = histos[0]->Integral(1, nbinsx, biny, biny);
	 double Int_T2 = histos[1]->Integral(1, nbinsx, biny, biny);
	 double Int_T4 = histos[2]->Integral(1, nbinsx, biny, biny);
	 // something related to phase factor, this is by guess

	 double mysgn = 1.;
	 if(fai < 0.) 
	   {
	     mysgn = -1.;
	   }

	 double integral = (1.-fabs(fai)) * Int_T1 + fabs(fai) * Int_T2  + mysgn*sqrt((1.-fabs(fai))*fabs(fai)) * Int_T4; ;

	 integral = integral * dx * 4.; 
	  
	 return integral; 
	 
       }
       
       // integrate out  kdint, depend on kd
     case 2: 
       {
	 
	 int binx = histos[0]->GetXaxis()->FindBin(kd);

	 double Int_T1 = histos[0]->Integral(binx, binx, 1, nbinsy);
	 double Int_T2 = histos[1]->Integral(binx, binx, 1, nbinsy);
	 double Int_T4 = histos[2]->Integral(binx, binx, 1, nbinsy);

	 double mysgn = 1.;
	 if(fai < 0.) 
	   {
	     mysgn = -1.;
	   }

	 double integral = (1.-fabs(fai)) * Int_T1 + fabs(fai) * Int_T2 +  mysgn*sqrt((1.-fabs(fai))*fabs(fai)) * Int_T4; 
	 
	 // something related to phase factor, this is by guess
	 integral = integral * dy * 4.;

	 return integral;
       }
       
     case 3: 
       {
	 double Int_T1 = histos[0]->Integral();
	 double Int_T2 = histos[1]->Integral();
	 double Int_T4 = histos[2]->Integral();

	 double mysgn = 1.;
	 if(fai < 0.) 
	   {
	     mysgn = -1.;
	   }

	 double integral = (1.-fabs(fai)) * Int_T1 + fabs(fai) * Int_T2 + mysgn*sqrt((1.-fabs(fai))*fabs(fai)) * Int_T4; ;
	 int nbinsx = histos[0]->GetXaxis()->GetNbins();
	 int nbinsy = histos[0]->GetYaxis()->GetNbins();


	 integral = integral * dx * dy * 4.;

	 // std::cout << __LINE__ << " "<< integral << "\n";

	 return integral;
       }
       
     }
   
   assert(0) ;
   return 0 ;
}

