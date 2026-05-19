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

void fig_22()
{
  Int_t MG = 20;
  Int_t n, k,nn[15];
  Double_t   nm[17] = {0., 2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26., 28., 30., 32.};
  Double_t Pn22[17] = {0.0000135418, 0.00165414, 0.0181607, 0.0768303, 0.169488, 0.229648, 0.21323,
		       0.147448, 0.0812501, 0.0377213, 0.0154244, 0.00574488, 0.00199684,
		       0.000658835, 0.000208786, 0.0000640609, 0.0000135418 };

  Double_t     nm1[14] =  {2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.};
  Double_t  Pn22_e[14] = {0.1631, 1.7797, 7.8243, 16.7981, 22.9196, 21.5560, 14.5702,
			 8.2160, 3.6614, 1.6538, 0.5892, 0.1637, 0.0697, 0.0355};
  Double_t dPn22_e[14] = {0.0895, 0.2557, 0.5185, 0.7497, 0.8749, 0.8332, 0.6494,
                         0.4705, 0.2927, 0.1931, 0.1048, 0.0513, 0.0312, 0.0253};

    for (int i=0; i<14; i++)
      {
       Pn22_e[i] =  Pn22_e[i]/100.;
      dPn22_e[i] = dPn22_e[i]/100.;
     }

       TGraph *gg1 =new TGraph(17,nm,Pn22);
       gg1->SetTitle("Multiplicity Distribution in e^{+}e^{-} at 22 GeV;n;Pn");
       TGraphErrors *gg2 =new TGraphErrors(14,nm1,Pn22_e,0,dPn22_e);
       gg2->SetTitle("Multiplicity Distribution;n;Pn");
      TCanvas *c1=new TCanvas("c1","Multiplicity distribution in e^{+}e^{-} at 22 GeV",800,600);

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
