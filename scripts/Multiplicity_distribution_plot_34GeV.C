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

void fig_34()
{
  Int_t MG0 = 22;
  Int_t n, k,nn[16];
  Double_t   nm[19] = {0., 2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26., 28., 30., 32., 34., 36.};
  Double_t Pn34[19] = {0.0000043323, 0.000453571, 0.00580373, 0.0301386, 0.0854587, 0.154566, 0.1968, 0.189286, 0.145291, 0.0931342, 0.0518287, 0.025886, 0.0119307, 0.00518803, 0.0021643, 0.000876527, 0.000347421, 0.000135495, 0.0000521772 };

  Double_t    nm1[18] =  {2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.,30.,32.,34.,36.};
  Double_t Pn34_e[18] = {0.0447, 0.5733, 3.1675, 8.3797, 15.3130, 19.7927, 19.2450, 14.4335, 9.1819, 5.0623, 2.7161, 1.2348, 0.5173, 0.1927, 0.0831, 0.0400, 0.0132, 0.0043 };
  Double_t dPn34_e[18] = {0.0455, 0.0759, 0.1424, 0.2133, 0.2897, 0.3349, 0.3312, 0.2843, 0.2231, 0.1627, 0.1189, 0.0803, 0.0519, 0.0326, 0.0213, 0.0155, 0.0084, 0.0055};

    for (int i=0; i<18; i++)
      {
       Pn34_e[i] =  Pn34_e[i]/100.;
      dPn34_e[i] = dPn34_e[i]/100.;
     }

       TGraph *gg1 =new TGraph(19,nm,Pn34);
       gg1->SetTitle("Multiplicity Distribution in e^{+}e^{-} at 34.8 GeV;n;Pn");
       TGraphErrors *gg2 =new TGraphErrors(18,nm1,Pn34_e,0,dPn34_e);
       gg2->SetTitle("Multiplicity Distribution;n;Pn");
      TCanvas *c1=new TCanvas("c1","Multiplicity distribution in e^{+}e^{-} at 34.8 GeV",800,600);

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
