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

void fig_43()
{
  Int_t MG0 = 22;
  Int_t n, k,nn[17];
  Double_t   nm[20] = {0., 2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26., 28., 30., 32., 34., 36., 38. };
  Double_t Pn43[20] = {0.00000023271, 0.000302526, 0.00368101, 0.0189397, 0.0555489, 0.108551, 0.155745, 0.175265, 0.162086, 0.127327, 0.0869237, 0.0523579, 0.0280965, 0.0135145, 0.00585035, 0.0022862, 0.00080862, 0.000259527, 0.0000757812, 0.0000201858};

  Double_t    nm1[19] =  {2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.,30.,32.,34.,36., 38.};
  Double_t Pn43_e[19] = {0.0328, 0.3446, 2.2958, 5.1983, 10.7780, 15.7378, 17.7662, 16.4074, 12.3672, 8.6055, 5.2095, 2.9188, 1.3995, 0.6002, 0.2374, 0.0591, 0.0264, 0.0047, 0.0111};
  Double_t dPn43_e[19] = {0.0357, 0.1505, 0.2682, 0.3303, 0.4604, 0.5594, 0.5933, 0.5626, 0.4756, 0.3927, 0.2972, 0.2171, 0.1416, 0.0934, 0.0592, 0.0265, 0.0182, 0.0091, 0.0217};

    for (int i=0; i<19; i++)
      {
       Pn43_e[i] =  Pn43_e[i]/100.;
      dPn43_e[i] = dPn43_e[i]/100.;
     }

       TGraph *gg1 =new TGraph(20,nm,Pn43);
       gg1->SetTitle("Multiplicity Distribution in e^{+}e^{-} at 43.6 GeV;n;Pn");
       TGraphErrors *gg2 =new TGraphErrors(19,nm1,Pn43_e,0,dPn43_e);
       gg2->SetTitle("Multiplicity Distribution;n;Pn");
      TCanvas *c1=new TCanvas("c1","Multiplicity distribution in e^{+}e^{-} at 43.6 GeV",800,600);

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
