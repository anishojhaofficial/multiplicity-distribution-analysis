/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 ROOT Foundation,  CERN/PH-SFT                   *
 *                         
                                           *
**********************************************************************/
#include <iostream>
#include <cmath>
#include <vector>

#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "Math/MinimizerOptions.h"

Double_t fun(Double_t *x, Double_t *par) {
  
  Int_t MG = 20;
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

  if(D2<=n) {R=0;}
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
     if(P1 <= n) {Sum = 0.;}
 else {
   S=1;
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

void test_14()
{
  const Int_t npar = 6;
  
    TF1 *f1 = new TF1("f1",fun,2,26,6);
    /*    
    f1-> SetParameter(0, 2.);
    f1-> SetParameter(1,1.);
    f1->SetParLimits(1,0.08,5.);
    f1-> SetParameter(2,4.);
    f1-> SetParameter(3,10.);
    f1-> SetParameter(4,0.2);
    f1-> SetParameter(5,2.);
    */
    
   f1->SetParameter(0,14.);
   //  f1->FixParameter(0,14.);
   f1->SetParameter(1,0.08);
   //   f1->SetParLimits(0,0.08,20.);
   //   f1->FixParameter(1,2.6);
   f1->SetParameter(2,3.87);
   f1->SetParameter(3,17.7);
   //   f1->SetParLimits(3,10.,28.);
   f1->SetParameter(4,0.97);
   f1->SetParameter(5,2.0);

   f1->SetParNames("k_p","mean_m","n^h","N","alpha","Omega");

   //   Double_t xvalues1[4] = {-7.0,-5.0,7.0,5.0};
   //   Double_t yvalues1[4] = {1000.0,1000.0,1000.0,1000.0};
   //   Double_t evalues1[4] = {30.0,30.0,30.0,30.0};

   Double_t xvalues1[13] = {2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26.};   
   Double_t yvalues1[13] = {0.5783, 5.3986, 16.4745, 26.6517, 25.0602, 14.9795, 6.7216,
			    2.6178, 0.8703,  0.4682,  0.1463,  0.0262,  0.0067};
   Double_t evalues1[13] = {0.1800, 0.4449,  0.7230,  0.9146,  0.8596,  0.6175, 0.3885,
                            0.2336, 0.1306,  0.1159,  0.0606,  0.0208,  0.0097};
   for (int i=0; i<13; i++)
     
     {
       yvalues1[i] = yvalues1[i]/100.;
       evalues1[i] = evalues1[i]/100.;
	} 
	    
  
   TGraphErrors *gr1 = new TGraphErrors(13, xvalues1, yvalues1, 0, evalues1);
   //  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Fumili2");
   Double_t a;

   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   //   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Fumili2");
   gr1->Fit("f1");
   gr1->Draw("ACP");
   f1->Draw("SAME");

   std::cout << "  chi2 = " << f1->  GetChisquare() << std::endl;
   std::cout << "k_p =" << f1->GetParameter(0) << "+/-" << f1 ->GetParError(0) << std::endl;
      Double_t var = 0., mean = 0., P3 = 0., f2;
   
   for(int l=0; l<13; l++)
     {
       mean = mean+(2.*l+2.)*yvalues1[l];
       P3 = P3 + yvalues1[l];
     }
   for(int l1 = 0; l1<13; l1++)
     {var  = var + (2.*l1+2. - mean)*(2.*l1+2. - mean)*yvalues1[l1];
            }
   f2 = var - mean;
   std::cout << "norm = " << P3 << endl;
   std::cout << "mean = " << mean << endl;
   std::cout << "f2 = " << f2 << endl;
}
