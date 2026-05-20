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

void Pn_14()
{
  Int_t MG = 20;
  Int_t n, k, nn[13];
  Double_t S, S1, S2, S3, R1, R2, R3, R4;
  Double_t P0, Q, Q1, Q2;
  Double_t Sum, SM;

  Double_t par[6] = {80.5703,  0.08171171, 4.4675, 27.8775, 0.976031, 1.99769};
                 // k_p, \overline m, n^h,    N,   alp, Om
  for(int j = 0; j < 14; j++)
    {nn[j] = 2*j+2;
      n = nn[j];

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

  R=1.;
 for (int ln=0; ln<n; ln++)
   {R = R*(2.*par[3] - ln)*R3/(ln +1.);// C_(2N)^n*((nh/N)/(1-nh/N))^n/n!
   }

  R = R*R4;  // R*(1-nh/N)^2N
 
 Sum = 0.;
 for (int m=1; m<=MG; m++)
   {S=1;
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
   }
 Q1=(Sum+R)*P0*par[5];
 
 Double_t Pn_14 = Q1;
 std::cout << "\n " << " n = " << n  << " Pn(" << n << ") = " << Pn_14 << endl;
    }

    Double_t D = 1., D0, D1, D3, n0 = 0.;
  D3 = par[3]*par[4];
  D1 = TMath::Power(R2,D3);
  for(int lg=0; lg<MG; lg++)
    {D = D*(par[0]+lg)*D1*S3/(lg+1.);
    }
  D0 = P0*(1.+D)*R4*par[5];
  std::cout << "\n " << " n = " << n0  << " Pn(" << n0 << ") = " << D0 << endl; 

}
