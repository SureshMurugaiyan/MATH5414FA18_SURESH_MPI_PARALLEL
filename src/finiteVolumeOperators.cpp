#include <iostream>
#include <math.h>
using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Divergence of a Vector with variable coefficient- term in momentum eqn   !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void Div(double* Dn, double* Phi, double* U, double* V, int row, int col,double delX,double delY){

for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
      int k    = i*col+j;
   double PhiP = Phi[k];
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];

   double UP = U[k];
   double UE = U[k+1];
   double UW = U[k-1];
   double UN = U[k-col];
   double US = U[k+col];

   double VP = V[k];
   double VE = V[k+1];
   double VW = V[k-1];
   double VN = V[k-col];
   double VS = V[k+col];

   double Ee  = 0.5*(UE*PhiE+UP*PhiP);
   double Ew  = 0.5*(UW*PhiW+UP*PhiP);
   double Fn  = 0.5*(VN*PhiN+VP*PhiP);
   double Fs  = 0.5*(VS*PhiS+VP*PhiP);
   Dn[k]      = delX*(Fn-Fs)+delY*(Ee-Ew);
      }
   }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Divergence of a Vector with No-coefficient- in continuity & source term  !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

void Divergence(double* Dn, double* U, double* V,int row, int col, double delX, double delY){

for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
       int k = i*col+j;
   double UP = U[k];
   double UE = U[k+1];
   double UW = U[k-1];
   double UN = U[k-col];
   double US = U[k+col];

   double VP = V[k];
   double VE = V[k+1];
   double VW = V[k-1];
   double VN = V[k-col];
   double VS = V[k+col];

   double Ue = 0.5*(UE+UP);
   double Uw = 0.5*(UW+UP);
   double Un = 0.5*(UN+UP);
   double Us = 0.5*(US+UP);

   double Ve = 0.5*(VE+VP);
   double Vw = 0.5*(VW+VP);
   double Vn = 0.5*(VN+VP);
   double Vs = 0.5*(VS+VP);

  Dn[k] = (Ue-Uw)*delY+(Vn-Vs)*delX;
   }
 }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Laplacian of a Scalar                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void Laplacian(double* Ln, double *Phi, int row, int col, double delX, double delY){
for(int i = 1; i<(row-1); i++){
 for(int j =1; j<(col-1); j++){
   int k = i*col+j;
   double PhiP = Phi[k];
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];

   double Ee  = (PhiE-PhiP)/delX;
   double Ew  = (PhiP-PhiW)/delX;
   double Fn  = (PhiN-PhiP)/delY;
   double Fs  = (PhiP-PhiS)/delY;
   Ln[k]      = delX*(Fn-Fs)+delY*(Ee-Ew);
     }
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Gradient                                                                 !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void gradient(double* gradxPhi,double* gradyPhi,double* Phi,
                        int row, int col, double delX, double delY){
for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){

   int       k = i*col+j;
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];
   double PhiP = Phi[k];

   double Phie = 0.5*(PhiE + PhiP);
   double Phiw = 0.5*(PhiW + PhiP);
   double Phin = 0.5*(PhiN + PhiP);
   double Phis = 0.5*(PhiS + PhiP);

   gradxPhi[k] = (Phie-Phiw)/delX;
   gradyPhi[k] = (Phin-Phis)/delY;
    }
  }
}
