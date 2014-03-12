#include "Riostream.h" 
#include "../interface/HZZ4L_OffShellWidthPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "TH2F.h"

using namespace TMath;

ClassImp(HZZ4L_OffShellWidthPdf) 

  HZZ4L_OffShellWidthPdf::HZZ4L_OffShellWidthPdf(const char *name, const char *title, 
						 RooAbsReal& _mass,
						 RooAbsReal& _widthKD,
						 RooAbsReal& _Gamma,
						 RooAbsReal& _r,
						 TH2F *_histoBkg,
						 TH2F *_histoSig,
						 TH2F *_histoInterf):
   RooAbsPdf(name,title), 
   mass("mass","mass",this,_mass),
   widthKD("widthKD","widthKD",this,_widthKD),
   Gamma("Gamma","Gamma",this,_Gamma),
   r("r","r",this,_r),
   histoBkg(_histoBkg),
   histoSig(_histoSig),
   histoInterf(_histoInterf)
 { 
 } 


 HZZ4L_OffShellWidthPdf::HZZ4L_OffShellWidthPdf(const HZZ4L_OffShellWidthPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   mass("mass",this,other.mass),
   widthKD("widthKD",this,other.widthKD),
   Gamma("Gamma",this,other.Gamma),
   r("r",this,other.r),
   histoBkg(other.histoBkg),
   histoSig(other.histoSig),
   histoInterf(other.histoInterf)
 { 
 } 



 Double_t HZZ4L_OffShellWidthPdf::evaluate() const 
 { 

   double value = 0.;

   int binx =  histoBkg->GetXaxis()->FindBin(mass);
   int biny =  histoBkg->GetYaxis()->FindBin(widthKD);
   
   Double_t T1 = histoBkg->GetBinContent(binx, biny);
   Double_t T2 = histoSig->GetBinContent(binx, biny);
   Double_t T4 = histoInterf->GetBinContent(binx, biny);

   
   value = T1 + Gamma*r*T2 + sqrt(Gamma*r)*T4; 
   
   if ( value <= 0.) return 1.0e-200;
   
   return value ; 
   
 } 

Int_t HZZ4L_OffShellWidthPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if (matchArgs(allVars,analVars,RooArgSet(*mass.absArg(), *widthKD.absArg()))) return 3 ;
  if (matchArgs(allVars,analVars,mass)) return 1 ;
  if (matchArgs(allVars,analVars,widthKD)) return 2 ;

  return 0 ;

}

Double_t HZZ4L_OffShellWidthPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{

  int nbinsx = histoBkg->GetXaxis()->GetNbins();
  
  double xMin = histoBkg->GetXaxis()->GetBinLowEdge(1);
  double xMax = histoBkg->GetXaxis()->GetBinUpEdge(nbinsx);
  double dx = (xMax - xMin) / nbinsx; 

  
  int nbinsy = histoBkg->GetYaxis()->GetNbins();
  
  double yMin = histoBkg->GetYaxis()->GetBinLowEdge(1);
  double yMax = histoBkg->GetYaxis()->GetBinUpEdge(nbinsy);
  double dy = (yMax - yMin) / nbinsy; 
  
   switch(code)
     {

       // integrate out mass, depend on widthKD
     case 1: 
       {

	 int biny = histoBkg->GetYaxis()->FindBin(widthKD);
	 
	 double Int_T1 = histoBkg->Integral(1, nbinsx, biny, biny);
	 double Int_T2 = histoSig->Integral(1, nbinsx, biny, biny);
	 double Int_T4 = histoInterf->Integral(1, nbinsx, biny, biny);

	 double integral = Int_T1 + Gamma*r*Int_T2  + sqrt(Gamma*r)*Int_T4;

	 integral = integral * dx; 
	  
	 return integral; 
	 
       }
       
       // integrate out  widthKD, depend on mass
     case 2: 
       {
	 
	 int binx = histoBkg->GetXaxis()->FindBin(mass);

	 double Int_T1 = histoBkg->Integral(binx, binx, 1, nbinsy);
	 double Int_T2 = histoSig->Integral(binx, binx, 1, nbinsy);
	 double Int_T4 = histoInterf->Integral(binx, binx, 1, nbinsy);

	 double integral = Int_T1 + Gamma*r*Int_T2 +  sqrt(Gamma*r)*Int_T4; 
	 
	 integral = integral * dy;

	 return integral;
       }
       
     case 3: 
       {
	 double Int_T1 = histoBkg->Integral();
	 double Int_T2 = histoSig->Integral();
	 double Int_T4 = histoInterf->Integral();

	 double integral = Int_T1 + Gamma*r*Int_T2 + sqrt(Gamma*r)*Int_T4; 

	 integral = integral * dx * dy;

	 return integral;
       }
       
     }
   
   assert(0) ;
   return 0 ;
}

