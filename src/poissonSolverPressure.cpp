#include <iostream>
#include <math.h>
#include "input.h"
using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Function Declarations                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Main Function-->Poisson Solver for Pressure Finite Volume Solver         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void PoissonPressure(double* Phi, int row, int col,
             double delX,double delY,double* source,
             int totCell){
double lam = 1;

int itr = 0;
int stop = 0;
while (stop==0){
itr++;
for(int i=1; i<(row-1); i++){
 for(int j=1; j<(col-1); j++){
   int k = i*col+j;
   double PhiP = Phi[k];
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];
   double AP   = (-2*delY/delX)-(2*delX/delY);
   double AS   = (delX/delY);
   double AW   = (delY/delX);
   double AE   = (delY/delX);
   double AN   = (delX/delY);

   double  R     = source[k]- AP*PhiP-AE*PhiE-AW*PhiW-AN*PhiN-AS*PhiS;
   double delPhi = R/AP;
   Phi[k] = Phi[k]+lam*delPhi;
  }
 }
if(itr>MAXitrPr){stop=1;}
}
}
