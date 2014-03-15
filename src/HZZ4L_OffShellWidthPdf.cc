#include "Riostream.h" 
#include "../interface/HZZ4L_OffShellWidthPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "TIterator.h"

using namespace TMath;

ClassImp(HZZ4L_OffShellWidthPdf) 

  HZZ4L_OffShellWidthPdf::HZZ4L_OffShellWidthPdf(const char *name, const char *title, 
						 RooAbsReal& _mass,
						 RooAbsReal& _widthKD,
						 RooAbsReal& _Gamma,
						 RooAbsReal& _mu,
						 RooAbsReal& _kbkg,
						 const RooArgList& inHistList):
   RooAbsPdf(name,title), 
   mass("mass","mass",this,_mass),
   widthKD("widthKD","widthKD",this,_widthKD),
   Gamma("Gamma","Gamma",this,_Gamma),
   mu("mu","mu",this,_mu),
   kbkg("kbkg","kbkg",this,_kbkg),
   _histList("histList","ListOfHists",this)
 { 
   TIterator* histIter = inHistList.createIterator();
   RooAbsArg* func;
   while((func = (RooAbsArg*)histIter->Next()))
     {
       if (!dynamic_cast<RooAbsReal*>(func))
	 {
	   coutE(InputArguments) << "ERROR: :HZZ4L_OffShellWidthPdf(" << GetName() << ") Hist " << func->GetName() << " is not of type RooAbsReal" << endl;
	   assert(0);
	 }
       _histList.add(*func);
     }
   delete histIter;
   _histIter = _histList.createIterator();
 } 


 HZZ4L_OffShellWidthPdf::HZZ4L_OffShellWidthPdf(const HZZ4L_OffShellWidthPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   mass("mass",this,other.mass),
   widthKD("widthKD",this,other.widthKD),
   Gamma("Gamma",this,other.Gamma),
   mu("mu",this,other.mu),
   kbkg("kbkg",this,other.kbkg),
   _histList("histList",this,other._histList)
 { 
   _histIter = _histList.createIterator();
 } 



 Double_t HZZ4L_OffShellWidthPdf::evaluate() const 
 { 

   double value = 0.;
   
   Double_t T1 = dynamic_cast<const RooHistFunc*>(_histList.at(0))->getVal();
   Double_t T2 = dynamic_cast<const RooHistFunc*>(_histList.at(1))->getVal();
   Double_t T4 = dynamic_cast<const RooHistFunc*>(_histList.at(2))->getVal();

   double mysgn = 1;
   if (kbkg < 0.)
     {
       mysgn = -1.;
     }
   
   value = kbkg*T1 + Gamma*mu*T2 + mysgn*sqrt(Gamma*mu*fabs(kbkg))*T4; 
   
   if ( value <= 0.) return 1.0e-200;
   
   return value ; 
   
 } 

Int_t HZZ4L_OffShellWidthPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if (matchArgs(allVars,analVars,RooArgSet(*mass.absArg(), *widthKD.absArg()))) return 3 ;
  //if (matchArgs(allVars,analVars,mass)) return 1 ;
  //if (matchArgs(allVars,analVars,widthKD)) return 2 ;

  return 0 ;

}

Double_t HZZ4L_OffShellWidthPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{
  
   switch(code)
     {

       
       
     case 3: 
       {
	 
	 Double_t Int_T1 = dynamic_cast<const RooHistFunc*>(_histList.at(0))->analyticalIntegral(1000);
	 Double_t Int_T2 = dynamic_cast<const RooHistFunc*>(_histList.at(1))->analyticalIntegral(1000);
	 Double_t Int_T4 = dynamic_cast<const RooHistFunc*>(_histList.at(2))->analyticalIntegral(1000);

	 double mysgn = 1;
	 if (kbkg < 0.)
	   {
	     mysgn = -1.;
	   }

	 double integral = kbkg*Int_T1 + Gamma*mu*Int_T2 + mysgn*sqrt(Gamma*mu*fabs(kbkg))*Int_T4; 

	 integral = integral;

	 return integral;
       }
       
     }
   
   assert(0) ;
   return 0 ;
}

