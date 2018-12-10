#include <iostream>
#include <math.h>
#include <stdio.h>
#include <input.h>
using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Function Declarations                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void InitializeField(double *phi, int row, int col);
void Div(double* Dn, double* Phi, double* U, double* V, int row, int col,double delX,double delY);
void Laplacian(double* Ln, double *Phi, int row, int col, double delX, double delY);
void timeStep(double* delt,double* ux,double* uy);
void storeOldValue(double *phinew, double *phiOld,int totCell);
void eulerPredictor(double* ux, double *uy,double delT,
                    double* Cnx, double *Cny,
                    double* Dnx,double* Dny,
                    int col, int row,double vol);
void adamPredictor(double* ux, double *uy,double delT,
                      double* Cnx, double *Cny,
                      double* Dnx,double* Dny,
                      double* CnxOld, double *CnyOld,
                      double* DnxOld,double* DnyOld,
                      int col, int row,double vol);
void Divergence(double* Dn, double* U, double* V,int row, int col, double delX, double delY);
void updateCn(double* Cn,double dt, int col,int row);
void PoissonPressure(double* Phi, int row, int col,
             double delX,double delY,double* source,
             int totCell);
void corrector(double* ux, double* uy,
               double* gradxP,double* gradyP,
               double dt, int col, int row);
void gradient(double* gradxPhi,double* gradyPhi,double* Phi,
                        int row, int col, double delX, double delY);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Main Function                                                            !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

void doNavier(double* ux,double* uy,double* p,double* uxOld,double* uyOld,
                                        double* pOld, int it, int ncG, int nxcG, int nycG){

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Variable Declarations                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Declaration of arrays for storing Derived Variables
// allocating space for Diffusion term andConvection term
double *Dnx;Dnx = (double*) malloc(ncG * sizeof(double));
double *Dny;Dny = (double*) malloc(ncG * sizeof(double));
double *Cnx;Cnx = (double*) malloc(ncG * sizeof(double));
double *Cny;Cny = (double*) malloc(ncG * sizeof(double));

//allocating space for Convection term in poisson eqn & gradient
double *Cn;        Cn = (double*) malloc(ncG * sizeof(double));
double *gradxP;gradxP = (double*) malloc(ncG * sizeof(double));
double *gradyP;gradyP = (double*) malloc(ncG * sizeof(double));

// Storing previous timestep diffusion term and convection term
double *DnxOld;DnxOld = (double*) malloc(ncG * sizeof(double));
double *DnyOld;DnyOld = (double*) malloc(ncG * sizeof(double));
double *CnxOld;CnxOld = (double*) malloc(ncG * sizeof(double));
double *CnyOld;CnyOld = (double*) malloc(ncG * sizeof(double));

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Initialize all the matrices                                              !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
InitializeField(Dnx,nycG,nxcG);
InitializeField(Dny,nycG,nxcG);
InitializeField(Cnx,nycG,nxcG);
InitializeField(Cny,nycG,nxcG);
InitializeField(Cn,nycG,nxcG);
InitializeField(gradxP,nycG,nxcG);
InitializeField(gradyP,nycG,nxcG);
InitializeField(DnxOld,nycG,nxcG);
InitializeField(DnyOld,nycG,nxcG);
InitializeField(CnxOld,nycG,nxcG);
InitializeField(CnyOld,nycG,nxcG);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculate TimeStep at each iteration based on max velocity at each step  !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
double dt = 0.0001;   // Initializing time step
timeStep(&dt,ux,uy); // Calculate the actual timeStep

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Things to be done for first iteration alone                              !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
if(it==1){
// Calculate Laplacian and Divergence term from Initial conditions
Laplacian(Dnx,ux,nycG,nxcG,dx,dy);      //  Dnx Diffusion term
Laplacian(Dny,uy,nycG,nxcG,dx,dy);      //  Dny Diffusion term
Div(Cnx,ux,ux,uy,nycG,nxcG,dx,dy);      //  Cnx Convection term
Div(Cny,uy,ux,uy,nycG,nxcG,dx,dy);      //  Cny Convection term
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Store old values                                                         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
storeOldValue(ux,uxOld,ncG);
storeOldValue(uy,uyOld,ncG);
storeOldValue(p,pOld,ncG);

storeOldValue(Dnx,DnxOld,ncG);
storeOldValue(Dny,DnyOld,ncG);
storeOldValue(Cnx,CnxOld,ncG);
storeOldValue(Cny,CnyOld,ncG);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of Laplacian and Divergence for Predictor step               !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
 Laplacian(Dnx,ux,nycG,nxcG,dx,dy);      //  Dnx Diffusion term
 Laplacian(Dny,uy,nycG,nxcG,dx,dy);      //  Dny Diffusion term
 Div(Cnx,ux,ux,uy,nycG,nxcG,dx,dy);      //  Cnx Convection term
 Div(Cny,uy,ux,uy,nycG,nxcG,dx,dy);      //  Cny Convection term
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Predictor                                                                !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
//eulerPredictor(ux,uy,dt,Cnx,Cny,Dnx,Dny,nxcG,nycG,dV);
adamPredictor(ux,uy,dt,Cnx,Cny,Dnx,Dny,CnxOld,CnyOld,DnxOld,DnyOld,nxcG,nycG,dV);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of Source term in Poisson Equation                           !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
Divergence(Cn,ux,uy,nycG,nxcG,dx,dy); //  source term in poisson equation
updateCn(Cn,dt,nxcG,nycG);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Solver Poisson Equation For Pressure                                     !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
PoissonPressure(p,nycG,nxcG,dx,dy,Cn,ncG);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of pressure gradient                                         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
gradient(gradxP,gradyP,p,nycG,nxcG,dx,dy);  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Corrector Step                                                           !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
corrector(ux,uy,gradxP,gradyP,dt,nxcG,nycG);

free(Dnx);
free(Dny);
free(Cnx);
free(Cny);
free(Cn);
free(gradxP);
free(gradyP);
free(DnxOld);
free(DnyOld);
free(CnxOld);
free(CnyOld);

}
