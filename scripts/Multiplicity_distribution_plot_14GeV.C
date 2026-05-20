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

void fig_14()
{
  Int_t MG = 20;
  Int_t n, k,nn[14];
  Double_t   nm[15] = {0.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26., 28.};
  Double_t Pn14[15] = {0.00010862, 0.00604774, 0.0521548, 0.167115, 0.26675, 0.248167,
		       0.150883,0.0674846, 0.0256665, 0.00939889, 0.00338389,
		       0.00115272,0.000369063, 0.00011391, 0.0000344632};
  Double_t     nm1[13] =  {2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26.};
  Double_t  Pn14_e[13] = {0.5783, 5.3986, 16.4745, 26.6517, 25.0602,
			 14.9795, 6.7216, 2.6178, 0.8703,  0.4682,  0.1463,  0.0262,  0.0067 };
  Double_t dPn14_e[13] = {0.1800, 0.4449,  0.7230,  0.9146,  0.8596,  0.6175, 0.3885,
                          0.2336, 0.1306,  0.1159,  0.0606,  0.0208,  0.0097 };

    for (int i=0; i<13; i++)
      {
       Pn14_e[i] =  Pn14_e[i]/100.;
      dPn14_e[i] = dPn14_e[i]/100.;
     }

       TGraph *gg1 =new TGraph(15,nm,Pn14);
       gg1->SetTitle("Multiplicity Distribution in e^{+}e^{-} at 14 GeV;n;Pn");
       TGraphErrors *gg2 =new TGraphErrors(13,nm1,Pn14_e,0,dPn14_e);
       gg2->SetTitle("Multiplicity Distribution;n;Pn");
      TCanvas *c1=new TCanvas("c1","Multiplicity distribution in e+e- at 14 GeV",800,600);

      gg1->SetLineColor(4);
      gg2->SetMarkerStyle(20);
      gg1->SetLineWidth(3);
      gg1->Draw("AC");
      gg2->SetLineWidth(3);
      gg2->SetMarkerSize(1);
      gg2->SetLineColor(2);
      gg2->SetMarkerStyle(21);
      gg2->SetMarkerColor(2);
      gg2->Draw("P");

      auto *leg = new TLegend(0.1,0.7, 0.48,0.9);
      leg->AddEntry(gg1,"  -   GDM", "lp");
      leg->AddEntry(gg2, "  -  TASSO", "ep"); //lep for h/v errors
      leg->Draw();
           
      /** 
  for(int j=0; j<10;j++)
    
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
      **/
}
