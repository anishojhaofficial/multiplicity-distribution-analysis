/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 ROOT Foundation,  CERN/PH-SFT                   *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/MinimizerOptions.h"

Double_t fun(Double_t *x, Double_t *par) {
  
  Int_t MG = 22;
  Int_t n, k;
  Double_t S, S1, S2, S3, R1, R2, R3, R4;
  Double_t P0, P1, Q, Q1, Q2;
  Double_t Sum, SM;
  n = x[0];
 // Double_t par[6] = {4.91,     3.01,   4.34,  20.2, 0.2, 2.0};
                 // k_p, \overline m, n^h,    N,   alp, Om
 Double_t R;
 Double_t D2;

  //qg-cascade
  S1=par[0] + par[1]; //kp+m
  S2=par[0]/S1;      //kp/(kp+m)
  S3=par[1]/S1;     //m/(kp+m)

  // hadronization 
  R1=par[2]/par[3];  //nh/N 
  R2=1.-R1;          // 1-nh/N
  R3=R1/R2;          //(nh/N)/(1-nh/N) 
  D2=2.*par[3];  //2N
  
  P0 = TMath::Power(S2,par[0]); //[kp/(m+kp)]^kp
  R4 = TMath::Power(R2,D2);      //(1-nh/N)^2N

   if(D2<=n) R=0;
  else
    {R=1;
      
 for (int ln=0; ln<n; ln++)
   {R = R*(2.*par[3] - ln)*R3/(ln +1.);// C_(2N)^n*((nh/N)/(1-nh/N))^n/n!
   }

  R = R*R4;  // R*(1-nh/N)^2N
    }
 Sum = 0.;
 for (int m=1; m<=MG; m++)
   {P1 = (2. +par[4]*m)*par[3];
     if(P1 < n) {Sum = 0.;}
 else {
   S=1.;
for (int m1=1; m1 <= m; m1++)
  {
    S=S*(par[0]+m1-1.)*S3/m1; // kp*(kp+1)...(kp+m-1)(m/(m+kp)^m/m!
  }
    Q=1;
    
    for (int l=1; l<=n; l++)
      {
	Q=Q*((2.+par[4]*m)*par[3] -l+1.)*R3/l;
      }
    
    Q2 = TMath::Power(R2,(2.+par[4]*m)*par[3]);
    SM=S*Q*Q2;
    Sum=Sum+SM;
 }}
 Q1=(Sum+R)*P0*par[5];

  return Q1;
}

void test_34(){
  const Int_t npar = 6;

    TF1 *f1 = new TF1("f1",fun,2,36,6);
   f1->SetParameter(0,2.5);
   //   f1->SetParLimits(0,4.3,6.);
   f1->SetParameter(1,6.0);
   //   f1->SetParLimits(1,3.1,3.1);
   f1->SetParameter(2,5.5);
   //   f1->SetParLimits(2,4.,5.);
   f1->SetParameter(3,50.0);
   f1->SetParameter(4,0.1);
   f1->SetParameter(5,2.0);
   f1->SetParNames("k_p","m_mean","n^h","N","alpha","Omega");
   //   Double_t xvalues1[4] = {-7.0,-5.0,7.0,5.0};
   //   Double_t yvalues1[4] = {1000.0,1000.0,1000.0,1000.0};
   //   Double_t evalues1[4] = {30.0,30.0,30.0,30.0};
   //   Data at 22 GeV
   
   Double_t xvalues1[18] = {2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 32., 34., 36.};
   Double_t yvalues1[18] = {0.0447, 0.5733, 3.1675, 8.3797, 15.3130, 19.7927, 19.2450, 14.4335, 9.1819, 5.0623, 2.7161, 1.2348, 0.5173, 0.1977, 0.0831, 0.0400, 0.0132, 0.0043};
   Double_t evalues1[18] = {0.0455, 0.0759, 0.1424, 0.2133, 0.2897, 0.3349, 0.3312, 0.2843, 0.2231, 0.1627, 0.1189, 0.0803, 0.0519, 0.0326, 0.0213, 0.0155, 0.0084, 0.0055};
                            
                            
   for (int i=0; i<18; i++)
     
     {
       yvalues1[i] = yvalues1[i]/100.;
       evalues1[i] = evalues1[i]/100.;
	} 
   
   TGraphErrors *gr1 = new TGraphErrors(18, xvalues1, yvalues1, 0, evalues1);
     ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   //    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Fumili2");
								    
   gr1->Fit("f1");
   gr1->Draw("ACP");
   //   f1->Draw("SAME");
   
   std::cout << "  chi2 = " << f1->  GetChisquare() << std::endl;
   std::cout << "k_p =" << f1->GetParameter(0) << "+/-" << f1 ->GetParError(0) << std::endl;
   
   Double_t var = 0., SK = 0., SL = 0., mean = 0., f2;
   for(int l=0; l<18; l++)
     {
       mean = mean+(2.*l+2.)*yvalues1[l];
         SK = SK + yvalues1[l];
            }
   //  mean = mean/SK;
     for(int l1=0; l1<18; l1++)
      {
	SL  = SL +  (2.*l1+2.-mean)*(2.*l1+2.-mean)*yvalues1[l1];
	     }
     //  mean = mean/SK;
     //   SL = SL/SK;
   f2 = SL - mean;
   std::cout << "Sum of Pn = " << SK << " var = " << SL << endl;
   std::cout << "mean = " << mean << endl;
   std::cout << "f2 = " << f2 << endl; 
  }
     
